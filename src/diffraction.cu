#include "diffraction.h"
#include "config.h"
#include <thrust/sort.h>
#include <thrust/reduce.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <thrust/pair.h>
#include <thrust/extrema.h>
#include <thrust/partition.h>
#include <thrust/count.h>

__global__ void CUDA_scattering_at_k(float* real_part,float * imag_part, int * atomic_number, const float * sf_cache,
				     const float * HKL_list,const float * pos,const int k,const int natoms);

#define ELEMENTS 100
static float atomsf[ELEMENTS][9] = 
#include "atomsf.cdata"

static int atomsf_initialized = 0;

static void fill_ff_tables(){
  FILE * ff;
  char line[1024];
  char * p;
  int z;
  int i;
  ff = fopen("atomsf.lib","r");
  if(!ff){
/*    fprintf(stderr,"Info: atomsf.lib not found in current directory, using internal atomsf.lib\n");*/
    return;
  }
  
  /* Init atomsf array */
  for(i = 0;i<ELEMENTS;i++){
    atomsf[i][0] = -1;
  }
  /* skip the comments */
  for(fgets(line,1024,ff);strstr(line,"AD") == line;fgets(line,1024,ff));

  /* we're at the beginning of an atom record */
  while(line[0]){
    fgets(line,1024,ff);
    /* get Z */
    sscanf(line,"%d",&z);
    if(atomsf[z][0] == -1){
      p = line+23;
      sscanf(p,"%f",&atomsf[z][8]);
      fgets(line,1024,ff);
      p = line;
      sscanf(p,"%f",&atomsf[z][0]);
      p+= 16; 
      sscanf(p,"%f",&atomsf[z][1]);
      p+= 16; 
      sscanf(p,"%f",&atomsf[z][2]);
      p+= 16; 
      sscanf(p,"%f",&atomsf[z][3]);
      fgets(line,1024,ff);
      p = line;
      sscanf(p,"%f",&atomsf[z][4]);
      p+= 16; 
      sscanf(p,"%f",&atomsf[z][5]);
      p+= 16; 
      sscanf(p,"%f",&atomsf[z][6]);
      p+= 16; 
      sscanf(p,"%f",&atomsf[z][7]);
      /* get the last line of the atom record */
      fgets(line,1024,ff);
      /* get the first line of the next atom record */
      line[0] = 0;
      fgets(line,1024,ff);      
    }else{
      /* skip record */
      fgets(line,1024,ff);      
      fgets(line,1024,ff);      
      fgets(line,1024,ff);      
      line[0] = 0;
      fgets(line,1024,ff);            
    }
  }
}



/* d should be the size of the scattering vector |H| */
/* d is in m^-1 but we need it in A^-1 so divide  by 1^10*/
static float  scatt_factor(float d,int Z,float B){
  float res = 0;
  int i;  
  d *= 1e-10;
  /* the 0.25 is there because the 's' used by the aproxumation is 'd/2' */
  for(i = 0;i<4;i++){
    res+= atomsf[Z][i]*exp(-(atomsf[Z][i+4]+B)*d*d*0.25);
  }                
  res += atomsf[Z][8]*exp(-B*d*d/0.25);
  return res;    
}

static float ilumination_function(Experiment * exper,float * pos){
  float dist2;
  float sigma;
  /* If no fwhm is defined just return 1 everywhere */
  if(!exper->beam_fwhm){
    return 1;
  }
  /* calculate distance from the center of the beam */
  dist2 = (pos[0]-exper->beam_center_x)*(pos[0]-exper->beam_center_x)+(pos[1]-exper->beam_center_y)*(pos[1]-exper->beam_center_y);
  sigma = exper->beam_fwhm/2.355;
  return exp(-dist2/(2*sigma*sigma));
}

Diffraction_Pattern * cuda_compute_pattern_on_list(Molecule * mol, float * HKL_list, int HKL_list_size, float B,Experiment * exp,Options * opts){
#ifndef _USE_CUDA
  sp_error_fatal("Can't use cuda when not compiled for CUDA");
#else
  int i,j;
  float scattering_vector_length;
  float scattering_factor_cache[ELEMENTS];
  int is_element_in_molecule[ELEMENTS];
  Diffraction_Pattern * res = (Diffraction_Pattern *)malloc(sizeof(Diffraction_Pattern));
  int points_per_percent;
  float * atom_ilumination = (float *)malloc(sizeof(float)*mol->natoms);
  int threads_per_block = 64;
  int number_of_blocks = (mol->natoms+threads_per_block-1)/threads_per_block;
 
  if(!atomsf_initialized){
    fill_ff_tables();
    atomsf_initialized = 1;
  }

  res->F = (Complex *)malloc(sizeof(Complex)*HKL_list_size);
  res->ints = (float *)malloc(sizeof(float)*HKL_list_size);
  res->HKL_list = (float *)malloc(sizeof(float)*3*HKL_list_size);
  memcpy(res->HKL_list,HKL_list,sizeof(float)*3*HKL_list_size);
  res->HKL_list_size = HKL_list_size;
  float * d_real_part;
  cutilSafeCall(cudaMalloc((void **)&d_real_part,sizeof(float)*mol->natoms)); 
  float * d_imag_part;
  cutilSafeCall(cudaMalloc((void **)&d_imag_part,sizeof(float)*mol->natoms));
  int * d_atomic_number;
  cutilSafeCall(cudaMalloc((void **)&d_atomic_number,sizeof(int)*mol->natoms));
  cutilSafeCall(cudaMemcpy(d_atomic_number,mol->atomic_number,sizeof(int)*mol->natoms,cudaMemcpyHostToDevice));
  float * d_atomic_pos;
  cutilSafeCall(cudaMalloc((void **)&d_atomic_pos,sizeof(float)*mol->natoms*3));
  cutilSafeCall(cudaMemcpy(d_atomic_pos,mol->pos,sizeof(float)*mol->natoms*3,cudaMemcpyHostToDevice));
  float * d_sf_cache;
  cutilSafeCall(cudaMalloc((void **)&d_sf_cache,sizeof(float)*ELEMENTS));
  float * d_HKL_list;
  cutilSafeCall(cudaMalloc((void **)&d_HKL_list,sizeof(float)*3*HKL_list_size));
  cutilSafeCall(cudaMemcpy(d_HKL_list,HKL_list,sizeof(float)*3*HKL_list_size,cudaMemcpyHostToDevice));

  for(j = 0 ;j< ELEMENTS;j++){
    is_element_in_molecule[j] = 0;
  }
  for(j = 0 ;j< mol->natoms;j++){
    is_element_in_molecule[mol->atomic_number[j]] = 1;
    atom_ilumination[j] = ilumination_function(exp,&(mol->pos[j*3]));
  }

  points_per_percent = 1+(HKL_list_size)/100;
  for(i = 0;i<HKL_list_size;i++){
    if(i % points_per_percent == 0){
      fprintf(stderr,"%f percent done\n",(100.0*i)/(HKL_list_size));
    }
    sp_real(res->F[i]) = 0;
    sp_imag(res->F[i]) = 0;
    scattering_vector_length = sqrt(HKL_list[3*i]*HKL_list[3*i]+HKL_list[3*i+1]*HKL_list[3*i+1]+HKL_list[3*i+2]*HKL_list[3*i+2]);
    for(j = 0;j<ELEMENTS;j++){
      if(is_element_in_molecule[j]){
	scattering_factor_cache[j] = scatt_factor(scattering_vector_length,j,B);
      }
    }
    cutilSafeCall(cudaMemcpy(d_sf_cache,scattering_factor_cache,sizeof(float)*ELEMENTS,cudaMemcpyHostToDevice));
    CUDA_scattering_at_k<<<number_of_blocks, threads_per_block>>>(d_real_part,d_imag_part,d_atomic_number,d_sf_cache,d_HKL_list,d_atomic_pos,i,mol->natoms);
    thrust::device_ptr<float> begin =  thrust::device_pointer_cast(d_real_part);
    thrust::device_ptr<float> end =  thrust::device_pointer_cast(d_real_part+mol->natoms);
    sp_real(res->F[i]) = thrust::reduce(begin, end);
    begin =  thrust::device_pointer_cast(d_imag_part);
    end =  thrust::device_pointer_cast(d_imag_part+mol->natoms);
    sp_imag(res->F[i]) = thrust::reduce(begin, end);
    res->ints[i] = sp_cabs(res->F[i])*sp_cabs(res->F[i]);
  }
  return res;  
#endif 
}

__global__ void CUDA_scattering_at_k(float* real_part,float * imag_part, int * atomic_number, const float * sf_cache,
				     const float * HKL_list,const float * pos,const int k,const int natoms){
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  if(!atomic_number[i]){
    return;
  }
  if(i<natoms){
    float sf = sf_cache[atomic_number[i]];
    float tmp = 2*3.14159265F*(HKL_list[3*k]*-pos[i*3]+HKL_list[3*k+1]*-pos[i*3+1]+HKL_list[3*k+2]*-pos[i*3+2]);
    real_part[i] = sf*cos(tmp);
    imag_part[i] = sf*sin(tmp);
  }
}
