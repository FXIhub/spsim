/*
 * Copyright (c) 2006 Filipe Maia
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */


#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <complex.h>
#include <spimage.h>
#include <nfft3.h>
//#include <nfft/options.h>
//#include <nfft/window_defines.h>
#include "config.h"
#include "diffraction.h"
#include "mpi_comm.h"
#include "io.h"
#include "box.h"
#include "sse_mathfun.h"

#ifndef M_PI
#define M_PI     3.1415926535897932384626433832795029L
#endif
#ifndef PI
#define PI     3.1415926535897932384626433832795029L
#endif

#define ELEMENTS 100
static float atomsf[ELEMENTS][9] = 
#include "atomsf.cdata"

static float atomed[ELEMENTS][11];

static int atomsf_initialized = 0;
static int atomed_initialized = 0;

void multiply_pattern_on_list_with_scattering_factor(complex double * f,int Z,float * HKL_list, int HKL_list_size,float B);


static float sp_mod(float a, float b){
  while(a < 0){
    a += b;
  }
  while(a>=b){
    a -= b;
  }
  return a;
}

/*
  Use atomsf.lib to obtain details for all atoms in a problem.
*/

/*
static void write_ff_tables(){
  int i;
  FILE * fp = fopen("atomsf.data","w");
  fprintf(fp," {\n");
  i = 0;
  fprintf(fp,"{%e,%e,%e,%e,%e,%e,%e,%e,%e}\n",atomsf[i][0],
	  atomsf[i][1],atomsf[i][2],atomsf[i][3],
	  atomsf[i][4],atomsf[i][5],atomsf[i][6],
	  atomsf[i][7],atomsf[i][8]);
  for(i =1;i<ELEMENTS;i++){
    fprintf(fp,",{%e,%e,%e,%e,%e,%e,%e,%e,%e}\n",atomsf[i][0],
	    atomsf[i][1],atomsf[i][2],atomsf[i][3],
	    atomsf[i][4],atomsf[i][5],atomsf[i][6],
	    atomsf[i][7],atomsf[i][8]);
  };
  fprintf(fp,"};\n");   
  fclose(fp);
}
*/


/* The B is assumed to be the same for all the atoms */
static void fill_ed_tables(float B){
  if(B <= 0){
    B = 1;
  }
  for(int Z = 0;Z<ELEMENTS;Z++){
    //    double min_b = sp_min(sp_min(sp_min(atomsf[Z][4],atomsf[Z][5]),atomsf[Z][6]),atomsf[Z][7]);
    for(int i = 0;i<4;i++){
      atomed[Z][i] = 8*pow(M_PI,1.5)*(atomsf[Z][i])/pow((atomsf[Z][i+4]+B),1.5);
    }
    atomed[Z][4] = 8*pow(M_PI,1.5)*(atomsf[Z][8])/pow(B,1.5);
  }
}

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
  res += atomsf[Z][8]*exp(-B*d*d*0.25);
  return res;    
}


/* d should be the distance from the center of the atom */
/* d is in m but we need it in A so we multiply by 1^10*/
/* this formula results for the *spherical* fourier transform
   of the 4 gaussian aproximation of the scattering factors.
   Check the basic formulas paper for more details 

   
   The formula is:
   Sum i=1..4 8*Pi^(3/2)*ai/bi^(3/2)/exp(4*Pi^2*r^2/bi)

   And now the tricky part is what to do with c.
   I decided to add it to the gaussian with the smallest b (the widest gaussian)

   But this has all been factored in the fill_ed_table
*/
static float  electron_density(float d,int Z){
  float res = 0;
  int i;  
  d *= 1e10;
  for(i = 0;i<4;i++){
    res+= atomed[Z][i]*exp(-(atomed[Z][i+4])*d*d);
  }
  res+= atomed[Z][4]*exp(-1*d*d);
  return res;    
  }


static double ilumination_function(Experiment * exper,float * pos){
  double dist2;
  double sigma;
  /* If no fwhm is defined just return 1 everywhere */
  if(!exper->beam_fwhm){
    return 1;
  }
  /* calculate distance from the center of the beam */
  dist2 = (pos[0]-exper->beam_center_x)*(pos[0]-exper->beam_center_x)+(pos[1]-exper->beam_center_y)*(pos[1]-exper->beam_center_y);
  sigma = exper->beam_fwhm/2.355;
  printf("here\n");
  return exp(-dist2/(2*sigma*sigma));
}

float * get_HKL_list_for_detector(CCD * det, Experiment * exp,int * HKL_list_size){
  /* number of pixels */
  int nx, ny;
  /* pixel index */
  int x,y;
  /* physical location of pixel*/
  double px,py;
  /* reciprocal coordinates */
  double rx,ry;
  double real_to_reciprocal = 1.0/(det->distance*exp->wavelength);
  double ewald_radius = 1.0/exp->wavelength;
  double distance_to_ewald_sphere_center;

  float * HKL_list;
  int index = 0;
  
  nx = det->nx;
  ny = det->ny;

  HKL_list = malloc(sizeof(float)*nx*ny*3);
  for(x = 0;x<nx;x++){
    for(y = 0;y<ny;y++){


      /* 
	 Calculate the pixel coordinates in reciprocal space 	 
	 by dividing the physical position by detector_distance*wavelength.
	 
	 CCD center at (nx-1)/2,(ny-1)/2

	 Upper left corner of the detector with negative x and positive y
	 
	 Add detector center as it might not be the same as the beam
      */
      px = ((x-(nx-1.0)/2.0)/nx)*det->width+det->center_x;
      py = (((ny-1.0)/2.0-y)/ny)*det->height+det->center_y;

      rx = px*real_to_reciprocal;
      ry = py*real_to_reciprocal;
      
      /* Project pixel into Ewald sphere. */
      if(!det->spherical){
	distance_to_ewald_sphere_center = sqrt(rx*rx+ry*ry+ewald_radius*ewald_radius);
	HKL_list[index++] = rx * ewald_radius/distance_to_ewald_sphere_center;
	HKL_list[index++] = ry * ewald_radius/distance_to_ewald_sphere_center;
	HKL_list[index++] = ewald_radius-(ewald_radius * ewald_radius/distance_to_ewald_sphere_center);
      }else{
	HKL_list[index++] = rx;
	HKL_list[index++] = ry;
	HKL_list[index++] = 0;
      }
	 
    }
  }
  *HKL_list_size = nx*ny;
  printf("Last HKL %e %e %e\n",HKL_list[nx*ny*3-3],HKL_list[nx*ny*3-2],HKL_list[nx*ny*3-1]);
  return HKL_list;
}


SpRotation * apply_orientation_to_HKL_list(float ** HKL_list, int * HKL_list_size,Options * opts){
  SpRotation * rot;
  int d = 3;
  int i,j;
  sp_vector * v = sp_vector_alloc(3);
  real * tmp = v->data;
  if(opts->random_orientation){
    rot = sp_rot_uniform();
  }else{
    rot = sp_rot_euler(opts->euler_orientation[0],opts->euler_orientation[1],opts->euler_orientation[2]);
  }
  for(i = 0;i<*HKL_list_size;i++){
    v->data = &((*HKL_list)[i*d]);
    sp_vector * u = sp_matrix_vector_prod(rot,v);
    for(j = 0;j<d;j++){
      (*HKL_list)[i*d+j] = u->data[j];
    }
    sp_vector_free(u);		          
  }
  v->data = tmp;
  sp_vector_free(v);
  return rot;
}

float * get_HKL_list_for_3d_detector(CCD * det, Experiment * exp,int * HKL_list_size){
  /* number of pixels */
  int nx, ny,nz;
  /* pixel index */
  int x,y,z;
  /* physical location of pixel*/
  double px,py,pz;
  /* reciprocal coordinates */
  double rx,ry,rz;
  double real_to_reciprocal = 1.0/(det->distance*exp->wavelength);

  float * HKL_list;
  int index = 0;
  
  nx = det->nx;
  ny = det->ny;
  nz = det->nz;

  HKL_list = malloc(sizeof(float)*nx*ny*nz*3);
  for(x = 0;x<nx;x++){
    for(y = 0;y<ny;y++){
      for(z = 0;z<nz;z++){
      /* 
	 Calculate the pixel coordinates in reciprocal space 	 
	 by dividing the physical position by detector_distance*wavelength.
	 
	 CCD center at (nx)/2,(ny)/2

	 Upper left corner of the detector with negative x and positive y
      */
      px = ((x-(nx)/2.0)/nx)*det->width;
      py = ((y-(ny)/2.0)/ny)*det->height;
      pz = ((z-(nz)/2.0)/nz)*det->depth;

      rx = px*real_to_reciprocal;
      ry = py*real_to_reciprocal;
      rz = pz*real_to_reciprocal;
      
      /* Project pixel into Ewald sphere. */
      HKL_list[index++] = rx;
      HKL_list[index++] = ry;
      HKL_list[index++] = rz;      
      }
    }
  }
  *HKL_list_size = nx*ny*nz;
  return HKL_list;
}




void multiply_pattern_with_scattering_factor(complex double * f,int Z,int nx, int ny, int nz, double rs_pixel_x,double rs_pixel_y,double rs_pixel_z, float B){
  /* f is assumed to be C ordered*/
  int i = 0;
  for(int xi = -nx/2;xi<nx/2;xi++){
    float x = (float)xi/(nx)/rs_pixel_x;
    for(int yi = -ny/2;yi<ny/2;yi++){
      float y = (float)yi/(ny)/rs_pixel_y;
      for(int zi = -nz/2;zi<nz/2;zi++){
	float z = (float)zi/(nz)/rs_pixel_z;
	double distance = sqrt(x*x+y*y+z*z);
	float sf = scatt_factor(distance,Z,B);
	f[i] *= sf;
	i++;
      }
    }
  }
}


void multiply_pattern_on_list_with_scattering_factor(complex double * f,int Z,float * HKL_list, int HKL_list_size, float B){
  /* f is assumed to be C ordered*/
  int i = 0;
  for(i = 0;i<HKL_list_size;i++){
    double distance = sqrt(HKL_list[i*3]*HKL_list[i*3]+HKL_list[i*3+1]*HKL_list[i*3+1]+HKL_list[i*3+2]*HKL_list[i*3+2]);
    float sf = scatt_factor(distance,Z,B);
    f[i] *= sf;
  }
}

Diffraction_Pattern * compute_pattern_by_nfft(Molecule * mol, CCD * det, Experiment * exp, float B,float * HKL_list,Options * opts){
  double alpha_x = atan(det->width/(2.0 * det->distance));
  double alpha_y = atan(det->height/(2.0 * det->distance));
  double alpha_z = atan(det->depth/(2.0 * det->distance));
  /* Assuming spherical detector (which is a bit weird for a 3D detector)*/
  double smax_x = tan(alpha_x)/exp->wavelength;
  double smax_y = tan(alpha_y)/exp->wavelength;
  double smax_z = tan(alpha_z)/exp->wavelength;
  float rs_pixel_x = 1/(smax_x*2);  
  float rs_pixel_y = 1/(smax_y*2);  
  float rs_pixel_z = 1/(smax_z*2);  
  int is_element_in_molecule[ELEMENTS];
  fprintf(stderr,"Pixel size x-%e y-%e z-%e\n",rs_pixel_x,rs_pixel_y,rs_pixel_z);
  int nx = det->nx;
  int ny = det->ny;
  int nz = det->nz;
  int n_el_in_mol = 0;

  int mpi_skip = 1;
  int mpi_skip_flag = 0;

/* in meters defines the limit up to which we compute the electron density */
  Diffraction_Pattern * res = malloc(sizeof(Diffraction_Pattern));
  res->HKL_list_size = nx*ny*nz;
  res->F = malloc(sizeof(Complex)*res->HKL_list_size);
  res->ints = malloc(sizeof(float)*res->HKL_list_size);
  res->HKL_list = malloc(sizeof(float)*3*res->HKL_list_size);

  nfft_plan p;
  
#ifdef MPI  
  /* Skip a number of elements equivalent to the number of computers used */
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_skip);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_skip_flag);
#endif
  
  /*  if(!atomed_initialized){
    fill_ed_tables(B);
    atomed_initialized = 1;
    }*/
  for(int j = 0 ;j< ELEMENTS;j++){
    is_element_in_molecule[j] = 0;
  }
  for(int j = 0 ;j< mol->natoms;j++){
    if(is_element_in_molecule[mol->atomic_number[j]] == 0){
      n_el_in_mol++;
    }
    is_element_in_molecule[mol->atomic_number[j]]++;
  }

  for(int k = 0;k<nx*ny*nz;k++){
    sp_real(res->F[k]) = 0;
    sp_imag(res->F[k]) = 0;
  }
  for(int Z = 0 ;Z< ELEMENTS;Z++){
    if(!is_element_in_molecule[Z]){
      continue;
    }else if(is_element_in_molecule[Z]){
      mpi_skip_flag++;
      if(mpi_skip_flag < mpi_skip){
	continue;
      }
      mpi_skip_flag = 0;
      fprintf(stderr,"Calculating Z = %d\n",Z);
      /* more than 100 atoms of that kind in the molecule use nfft*/
      nfft_init_3d(&p,nx,ny,nz,is_element_in_molecule[Z]);

      int k = 0;
      for(int j = 0 ;j< mol->natoms;j++){
	if(mol->atomic_number[j] == Z){

	  p.f[k] = 1;
	  /* We have to multiply the position with the dimension of the box because the
	     fourier sample are taken between 0..1 (equivalent to 0..0.5,-0.5..0) */

	  /* For some unknown reason I have to take the negative of the position otherwise the patterns came out inverted.
	     I have absolutely no idea why! It can even be a bug in the NFFT library. I should check this out better
	  */
	  p.x[k*3] = -mol->pos[j*3]/(rs_pixel_x)/nx; 
	  p.x[k*3+1] = -mol->pos[j*3+1]/(rs_pixel_y)/ny; 
	  p.x[k*3+2] = -mol->pos[j*3+2]/(rs_pixel_z)/nz; 
	  k++;
	}
      }
      
      if(p.nfft_flags & PRE_ONE_PSI){
	nfft_precompute_one_psi(&p);
      }
      if(is_element_in_molecule[Z] < 100){
	ndft_adjoint(&p);  
      }else{
	nfft_adjoint(&p);  
      }
      
      if(!opts->delta_atoms){
	multiply_pattern_with_scattering_factor(p.f_hat,Z,nx,ny,nz,
						rs_pixel_x,rs_pixel_y,rs_pixel_z,B);
      }
      for(int k = 0;k<nx*ny*nz;k++){
	sp_real(res->F[k]) += creal(p.f_hat[k]);
	sp_imag(res->F[k]) += cimag(p.f_hat[k]);
      }   
      nfft_finalize(&p);
    }
  }
  for(int i = 0;i<nx*ny*nz;i++){
    res->ints[i] = sp_cabs(res->F[i])*sp_cabs(res->F[i]);
  }
  sum_patterns(res);
  return res;
}


Diffraction_Pattern * compute_pattern_on_list_by_nfft(Molecule * mol,float * HKL_list, int HKL_list_size,CCD * det, float B,Options * opts){
  int is_element_in_molecule[ELEMENTS];
/* in meters defines the limit up to which we compute the electron density */
  Diffraction_Pattern * res = malloc(sizeof(Diffraction_Pattern));
  res->HKL_list_size = HKL_list_size;
  res->F = malloc(sizeof(Complex)*HKL_list_size);
  res->ints = malloc(sizeof(float)*HKL_list_size);
  res->HKL_list = malloc(sizeof(float)*3*HKL_list_size);
  /* N is the cutoff frequency */
  double max_x[3] = {0,0,0};
  double max_v[3] = {0,0,0};
  int N[3] = {1000,1000,10};
  int n[3] = {1000,1000,10};
  int n_el_in_mol = 0;

  int mpi_skip = 1;
  int mpi_skip_flag = 0;
  nnfft_plan p;

#ifdef MPI  
  /* Skip a number of elements equivalent to the number of computers used */
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_skip);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_skip_flag);
#endif
  
  /*  if(!atomed_initialized){
    fill_ed_tables(B);
    atomed_initialized = 1;
    }*/
  for(int j = 0 ;j< ELEMENTS;j++){
    is_element_in_molecule[j] = 0;
  }
  for(int j = 0 ;j< mol->natoms;j++){
    for(int k = 0;k<3;k++){
      if(2*fabs(mol->pos[j*3+k]) > max_x[k]){
	max_x[k] = 2*fabs(mol->pos[j*3+k]);
      }
    }
    if(is_element_in_molecule[mol->atomic_number[j]] == 0){
      n_el_in_mol++;
    }
    is_element_in_molecule[mol->atomic_number[j]]++;
  }



  for(int j = 0;j<HKL_list_size;j++){
    sp_real(res->F[j]) = 0;
    sp_imag(res->F[j]) = 0;
    for(int k = 0;k<3;k++){
      if(2*fabs(HKL_list[j*3+k]) > max_v[k]){
	max_v[k] = 2*fabs(HKL_list[j*3+k]);
      }
    }
  }
  printf("max_x %e %e %e\n",max_x[0],max_x[1],max_x[2]);
  printf("max_v %e %e %e\n",max_v[0],max_v[1],max_v[2]);
  printf("max_v*max_x %e %e %e\n",max_x[0]*max_v[0],max_x[1]*max_v[1],max_x[2]*max_v[2]);
  
  for(int k = 0;k<3;k++){
    N[k] = ceil(max_x[k]*max_v[k]);
    max_v[k] = N[k]/max_x[k];
    n[k] = 2*N[k];
  }
/*  int tmp = N[0];
  N[0] = N[3];
  N[3] = tmp;
  tmp = n[0];
  n[0] = n[3];
  n[3] = tmp;*/
  for(int Z = 0 ;Z< ELEMENTS;Z++){
    if(!is_element_in_molecule[Z]){
      continue;
    }else if(is_element_in_molecule[Z]){
      mpi_skip_flag++;
      if(mpi_skip_flag < mpi_skip){
	continue;
      }
      mpi_skip_flag = 0;
      fprintf(stderr,"Calculating Z = %d\n",Z);
      /* more than 100 atoms of that kind in the molecule use nfft*/
      nnfft_init_guru(&p, 3, HKL_list_size, is_element_in_molecule[Z], N, n, 2,
                      PRE_PSI| PRE_PHI_HUT|
                      MALLOC_X| MALLOC_V| MALLOC_F_HAT| MALLOC_F);
/*      nnfft_init(&p,3,HKL_list_size,is_element_in_molecule[Z],N);*/

      int k = 0;
      for(int j = 0 ;j< mol->natoms;j++){
	if(mol->atomic_number[j] == Z){

	  p.f[k] = 1;    
	  /* We have to multiply the position with the dimension of the box because the
	     fourier sample are taken between 0..1 (equivalent to 0..0.5,-0.5..0) */
	  /* For some unknown reason I have to take the negative of the position otherwise the patterns came out inverted.
	     I have absolutely no idea why! It can even be a bug in the NFFT library. I should check this out better.
	     It might have to do with the order of the output!
	  */
	  p.x[k*3] = -mol->pos[j*3]/max_x[0];
	  p.x[k*3+1] = -mol->pos[j*3+1]/max_x[1];
	  p.x[k*3+2] = -mol->pos[j*3+2]/max_x[2];
	  k++;
	}
      }
      k = 0;
      for(int k = 0 ;k<HKL_list_size;k++){
	p.v[k*3] =  HKL_list[k*3]/max_v[0];
	p.v[k*3+1] =  HKL_list[k*3+1]/max_v[1];
	p.v[k*3+2] =  HKL_list[k*3+2]/max_v[2];
      }
      
      /** precompute psi, the entries of the matrix B */
      if(p.nnfft_flags & PRE_PSI)
	nnfft_precompute_psi(&p);
      
      if(p.nnfft_flags & PRE_FULL_PSI)
	nnfft_precompute_full_psi(&p);
      
      if(p.nnfft_flags & PRE_LIN_PSI)
	nnfft_precompute_lin_psi(&p);
      
      /** precompute phi_hut, the entries of the matrix D */
      if(p.nnfft_flags & PRE_PHI_HUT)
	nnfft_precompute_phi_hut(&p);
    
      if(is_element_in_molecule[Z] < 10){
	nndft_adjoint(&p);  
      }else{
	nnfft_adjoint(&p);  
      }
      if(!opts->delta_atoms){
	multiply_pattern_on_list_with_scattering_factor(p.f_hat,Z,HKL_list,HKL_list_size,B);
      }
      for(int k = 0;k<HKL_list_size;k++){
	sp_real(res->F[k]) += creal(p.f_hat[k]);
	sp_imag(res->F[k]) += cimag(p.f_hat[k]);
      }   
      nnfft_finalize(&p);
    }
  }
  for(int i = 0;i<HKL_list_size;i++){
    res->ints[i] = sp_cabs(res->F[i])*sp_cabs(res->F[i]);
  }
  sum_patterns(res);
  return res;
}


Diffraction_Pattern * compute_pattern_by_fft(Molecule * mol, CCD * det, Experiment * exp, float B){
  double alpha_x = atan(det->width/(2.0 * det->distance));
  double alpha_y = atan(det->height/(2.0 * det->distance));
  double alpha_z = atan(det->depth/(2.0 * det->distance));
  double smax_x = sin(alpha_x)/exp->wavelength;
  double smax_y = sin(alpha_y)/exp->wavelength;
  double smax_z = sin(alpha_z)/exp->wavelength;
  double rs_pixel_x = 1/(smax_x*2);  
  double rs_pixel_y = 1/(smax_y*2);  
  double rs_pixel_z = 1/(smax_z*2);  
  fprintf(stderr,"Pixel size x-%e y-%e z-%e\n",rs_pixel_x,rs_pixel_y,rs_pixel_z);
  int nx = det->nx;
  int ny = det->ny;
  int nz = det->nz;
  double max_atoms_radius = 5e-10; /* in meters defines the limit up to which we compute the electron density */
  int x_grid_radius = max_atoms_radius/rs_pixel_x;
  int y_grid_radius = max_atoms_radius/rs_pixel_y;
  int z_grid_radius = max_atoms_radius/rs_pixel_z;
  int total_el = 0;
  if(!atomed_initialized){
    fill_ed_tables(B);
    atomed_initialized = 1;
  }
  Image * rs = sp_image_alloc(nx,ny,nz);
  for(int j = 0 ;j< mol->natoms;j++){
    if(!mol->atomic_number[j]){
      /* undetermined atomic number */
      continue;
    }
    total_el += mol->atomic_number[j];
    float home_x = sp_mod(mol->pos[j*3]/rs_pixel_x,nx);
    float home_y = sp_mod(mol->pos[j*3+1]/rs_pixel_y,ny);
    float home_z = sp_mod(mol->pos[j*3+2]/rs_pixel_z,nz);
    for(int x = home_x-x_grid_radius;x<home_x+x_grid_radius;x++){
      for(int y = home_y-y_grid_radius;y<home_y+y_grid_radius;y++){
	for(int z = home_z-z_grid_radius;z<home_z+z_grid_radius;z++){
	  float ed = 0;
/*	  int z2 = 0;
	  int y2 = 0;
	  int x2 = 0;*/
	  for(int z2 = -1;z2<2;z2+=2){
	    for(int y2 = -1;y2<2;y2+=2){
	      for(int x2 = -1;x2<2;x2+=2){
		float dz = (home_z-(z+z2/2.0))*rs_pixel_z;
		float dy = (home_y-(y+y2/2.0))*rs_pixel_y;
		float dx = (home_x-(x+x2/2.0))*rs_pixel_x;
		float distance = sqrt(dz*dz+dy*dy+dx*dx);
		ed += electron_density(distance,mol->atomic_number[j])*rs_pixel_x*1e10*rs_pixel_y*1e10*rs_pixel_z*1e10;   
	      }
	    }
	  }
	  /* ed is the average of 27 points around the grid point */
	  ed *= rs_pixel_x*1e10*rs_pixel_y*1e10*rs_pixel_z*1e10/8;
	  /* Multiply by the voxel volume so that we have the number in electrons instead of electrons/a^3*/
/*	  if(!isnormal(ed)){
	    abort();
	  }*/
	  ed += sp_cabs(sp_image_get(rs,sp_mod(x,nx),sp_mod(y,ny),sp_mod(z,nz)));
	  sp_image_set(rs,sp_mod(x,nx),sp_mod(y,ny),sp_mod(z,nz),sp_cinit(ed,0));
	}
      }
    }
    fprintf(stderr,"%d done\n",j);
  }
  fprintf(stderr,"Total electrons - %d\n",total_el);
  sp_image_write(rs,"ed.vtk",0);  
  sp_image_write(rs,"ed.h5",sizeof(real));  
  Image * sf = sp_image_fft(rs);
  sp_image_free(rs);
  rs = sp_image_shift(sf);
  sp_image_free(sf);
  sf = rs;
  Diffraction_Pattern * res = malloc(sizeof(Diffraction_Pattern));
  res->HKL_list_size = nx*ny*nz;
  res->F = malloc(sizeof(Complex)*res->HKL_list_size);
  res->ints = malloc(sizeof(float)*res->HKL_list_size);
  res->HKL_list = malloc(sizeof(float)*3*res->HKL_list_size);
  int i = 0;
  float norm = 1.0;
  for(int x = 0;x<sp_image_x(sf);x++){
    for(int y = 0;y<sp_image_y(sf);y++){
      for(int z = 0;z<sp_image_z(sf);z++){
	res->F[i] = sp_cscale(sp_image_get(sf,x,y,z),norm);
	res->ints[i] = sp_cabs(res->F[i])*sp_cabs(res->F[i]);
	i++;
      }
    }
  }
  sp_image_free(sf);
  return res;
}

Diffraction_Pattern * compute_pattern_on_list(Molecule * mol, float * HKL_list, int HKL_list_size,float B,Experiment * exp,Options * opts){
  int timer = sp_timer_start();
  int i,j;
  float scattering_factor;
  float scattering_vector_length;
  float scattering_factor_cache[ELEMENTS];
  int is_element_in_molecule[ELEMENTS];
  Diffraction_Pattern * res = malloc(sizeof(Diffraction_Pattern));
  int HKL_list_start = 0;
  int HKL_list_end = 0;
  int points_per_percent;
  float * atom_ilumination = malloc(sizeof(float)*mol->natoms);
  get_my_loop_start_and_end(HKL_list_size,&HKL_list_start,&HKL_list_end);

  if(!atomsf_initialized){
    fill_ff_tables();
    atomsf_initialized = 1;
  }

  res->F = malloc(sizeof(Complex)*HKL_list_size);
  res->ints = malloc(sizeof(float)*HKL_list_size);
  res->HKL_list = malloc(sizeof(float)*3*HKL_list_size);
  memcpy(res->HKL_list,HKL_list,sizeof(float)*3*HKL_list_size);
  res->HKL_list_size = HKL_list_size;
  for(j = 0 ;j< ELEMENTS;j++){
    is_element_in_molecule[j] = 0;
  }
  for(j = 0 ;j< mol->natoms;j++){
    is_element_in_molecule[mol->atomic_number[j]] = 1;
    atom_ilumination[j] = ilumination_function(exp,&(mol->pos[j*3]));
  }

  points_per_percent = 1+(HKL_list_end-HKL_list_start)/100;
  for(i = HKL_list_start;i<HKL_list_end;i++){
#ifdef MPI    
    if(is_mpi_master()){
      if(i % points_per_percent == 0){
	fprintf(stderr,"%f percent done\n",(100.0*(i-HKL_list_start))/(HKL_list_end-HKL_list_start));
      }
    }
#else
      if(i % points_per_percent == 0){
	fprintf(stderr,"%f percent done\n",(100.0*(i-HKL_list_start))/(HKL_list_end-HKL_list_start));
      }

#endif

    sp_real(res->F[i]) = 0;
    sp_imag(res->F[i]) = 0;
    scattering_vector_length = sqrt(HKL_list[3*i]*HKL_list[3*i]+HKL_list[3*i+1]*HKL_list[3*i+1]+HKL_list[3*i+2]*HKL_list[3*i+2]);
    for(j = 0;j<ELEMENTS;j++){
      if(is_element_in_molecule[j]){
	scattering_factor_cache[j] = scatt_factor(scattering_vector_length,j,B);
      }
    }
    for(j = 0 ;j< mol->natoms;j++){
      if(!mol->atomic_number[j]){
	continue;
      }
      /* Multiply the scattering factor with the ilumination function (should it be the square root of it?)*/
      scattering_factor = scattering_factor_cache[mol->atomic_number[j]]*sqrt(atom_ilumination[j]);
/*      scattering_factor = 1;*/
      float tmp = 2*M_PI*(HKL_list[3*i]*-mol->pos[j*3]+HKL_list[3*i+1]*-mol->pos[j*3+1]+HKL_list[3*i+2]*-mol->pos[j*3+2]);
      if(!opts->delta_atoms){
	sp_real(res->F[i]) += scattering_factor*cos(tmp);
	sp_imag(res->F[i]) += scattering_factor*sin(tmp);
      }else{
	sp_real(res->F[i]) += cos(tmp);
	sp_imag(res->F[i]) += sin(tmp);
      }
    }
    res->ints[i] = sp_cabs(res->F[i])*sp_cabs(res->F[i]);
  }
  syncronize_patterns(res);
  printf("%g atoms.pixel/s\n",1.0e6*HKL_list_size*mol->natoms/sp_timer_stop(timer));
  return res;
}


Diffraction_Pattern * vector_compute_pattern_on_list(Molecule * mol, float * HKL_list, int HKL_list_size,float B,Experiment * exp,Options * opts){
  int timer = sp_timer_start();
  int i,j;
  float scattering_factor;
  float scattering_vector_length;
  float scattering_factor_cache[ELEMENTS];
  int is_element_in_molecule[ELEMENTS];
  Diffraction_Pattern * res = malloc(sizeof(Diffraction_Pattern));
  int HKL_list_start = 0;
  int HKL_list_end = 0;
  int points_per_percent;
  float * atom_ilumination = malloc(sizeof(float)*mol->natoms);
  get_my_loop_start_and_end(HKL_list_size,&HKL_list_start,&HKL_list_end);

  if(!atomsf_initialized){
    fill_ff_tables();
    atomsf_initialized = 1;
  }

  res->F = malloc(sizeof(Complex)*HKL_list_size);
  res->ints = malloc(sizeof(float)*HKL_list_size);
  res->HKL_list = malloc(sizeof(float)*3*HKL_list_size);
  memcpy(res->HKL_list,HKL_list,sizeof(float)*3*HKL_list_size);
  res->HKL_list_size = HKL_list_size;
  for(j = 0 ;j< ELEMENTS;j++){
    is_element_in_molecule[j] = 0;
  }
  for(j = 0 ;j< mol->natoms;j++){
    is_element_in_molecule[mol->atomic_number[j]] = 1;
    atom_ilumination[j] = sqrt(ilumination_function(exp,&(mol->pos[j*3])));
  }

  points_per_percent = 1+(HKL_list_end-HKL_list_start)/100;
  for(i = HKL_list_start;i<HKL_list_end;i++){
#ifdef MPI    
    if(is_mpi_master()){
      if(i % points_per_percent == 0){
	fprintf(stderr,"%f percent done\n",(100.0*(i-HKL_list_start))/(HKL_list_end-HKL_list_start));
      }
    }
#else
      if(i % points_per_percent == 0){
	fprintf(stderr,"%f percent done\n",(100.0*(i-HKL_list_start))/(HKL_list_end-HKL_list_start));
      }

#endif

    sp_real(res->F[i]) = 0;
    sp_imag(res->F[i]) = 0;
    scattering_vector_length = sqrt(HKL_list[3*i]*HKL_list[3*i]+HKL_list[3*i+1]*HKL_list[3*i+1]+HKL_list[3*i+2]*HKL_list[3*i+2]);
    for(j = 0;j<ELEMENTS;j++){
      if(is_element_in_molecule[j]){
	scattering_factor_cache[j] = scatt_factor(scattering_vector_length,j,B);
      }
    }
    for(j = 0 ;j< 4*(mol->natoms/4);j+=4){
      
      /* Multiply the scattering factor with the ilumination function (should it be the square root of it?)*/
      v4sf sf = {scattering_factor_cache[mol->atomic_number[j]]*atom_ilumination[j],
				    scattering_factor_cache[mol->atomic_number[j+1]]*atom_ilumination[j+1],
				    scattering_factor_cache[mol->atomic_number[j+2]]*atom_ilumination[j+2],
				    scattering_factor_cache[mol->atomic_number[j+3]]*atom_ilumination[j+3]};

      float tmp[4] = {2*M_PI*(HKL_list[3*i]*-mol->pos[j*3]+HKL_list[3*i+1]*-mol->pos[j*3+1]+HKL_list[3*i+2]*-mol->pos[j*3+2]),
		      2*M_PI*(HKL_list[3*i]*-mol->pos[(j+1)*3]+HKL_list[3*i+1]*-mol->pos[(j+1)*3+1]+HKL_list[3*i+2]*-mol->pos[(j+1)*3+2]),
		      2*M_PI*(HKL_list[3*i]*-mol->pos[(j+2)*3]+HKL_list[3*i+1]*-mol->pos[(j+2)*3+1]+HKL_list[3*i+2]*-mol->pos[(j+2)*3+2]),
		      2*M_PI*(HKL_list[3*i]*-mol->pos[(j+3)*3]+HKL_list[3*i+1]*-mol->pos[(j+3)*3+1]+HKL_list[3*i+2]*-mol->pos[(j+3)*3+2])};
      v4sf phase = __builtin_ia32_loadups(tmp);
      v4sf sin_phase;
      v4sf cos_phase;
      sincos_ps(phase,&sin_phase,&cos_phase);
      if(!opts->delta_atoms){
	sin_phase = __builtin_ia32_mulps(sin_phase,sf);
	cos_phase = __builtin_ia32_mulps(cos_phase,sf);
      }
      __builtin_ia32_storeups(tmp,cos_phase);
      float sum = 0;
      for(int ii = 0;ii<4;ii++){
	sum += tmp[ii];
      }
      sp_real(res->F[i]) += sum;
      __builtin_ia32_storeups(tmp,sin_phase);
      sum = 0;
      for(int ii = 0;ii<4;ii++){
	sum += tmp[ii];
      }
      sp_imag(res->F[i]) += sum;
    }
    for(;j< mol->natoms;j++){
      /* Multiply the scattering factor with the ilumination function (should it be the square root of it?)*/
      scattering_factor = scattering_factor_cache[mol->atomic_number[j]]*sqrt(atom_ilumination[j]);
/*      scattering_factor = 1;*/
      float tmp = 2*M_PI*(HKL_list[3*i]*-mol->pos[j*3]+HKL_list[3*i+1]*-mol->pos[j*3+1]+HKL_list[3*i+2]*-mol->pos[j*3+2]);
      sp_real(res->F[i]) += scattering_factor*cos(tmp);
      sp_imag(res->F[i]) += scattering_factor*sin(tmp);
    }
    res->ints[i] = sp_cabs(res->F[i])*sp_cabs(res->F[i]);
  }
  syncronize_patterns(res);
  printf("%g atoms.pixel/s\n",1.0e6*HKL_list_size*mol->natoms/sp_timer_stop(timer));
  return res;
}


Diffraction_Pattern * compute_fresnel_pattern_on_list(Molecule * mol, float * HKL_list, int HKL_list_size,float B,Experiment * exp){
  int i,j;
  double scattering_factor;
  double scattering_vector_length;
  double scattering_factor_cache[ELEMENTS];
  int is_element_in_molecule[ELEMENTS];
  Diffraction_Pattern * res = malloc(sizeof(Diffraction_Pattern));
  int HKL_list_start = 0;
  int HKL_list_end = 0;
  int points_per_percent;
  double * atom_ilumination = malloc(sizeof(double)*mol->natoms);
  double k;
  double distance = HKL_list[2]*exp->wavelength;
  get_my_loop_start_and_end(HKL_list_size,&HKL_list_start,&HKL_list_end);

  if(!atomsf_initialized){
    fill_ff_tables();
    atomsf_initialized = 1;
  }

  res->F = malloc(sizeof(Complex)*HKL_list_size);
  res->ints = malloc(sizeof(float)*HKL_list_size);
  res->HKL_list = malloc(sizeof(float)*3*HKL_list_size);
  memcpy(res->HKL_list,HKL_list,sizeof(float)*3*HKL_list_size);
  res->HKL_list_size = HKL_list_size;
  for(j = 0 ;j< ELEMENTS;j++){
    is_element_in_molecule[j] = 0;
  }
  for(j = 0 ;j< mol->natoms;j++){
    is_element_in_molecule[mol->atomic_number[j]] = 1;
    atom_ilumination[j] = ilumination_function(exp,&(mol->pos[j*3]));
  }

  points_per_percent = 1+(HKL_list_end-HKL_list_start)/100;
  for(i = HKL_list_start;i<HKL_list_end;i++){
#ifdef MPI    
    if(is_mpi_master()){
      if(i % points_per_percent == 0){
	fprintf(stderr,"%f percent done\n",(100.0*(i-HKL_list_start))/(HKL_list_end-HKL_list_start));
      }
    }
#else
      if(i % points_per_percent == 0){
	fprintf(stderr,"%f percent done\n",(100.0*(i-HKL_list_start))/(HKL_list_end-HKL_list_start));
      }

#endif

    sp_real(res->F[i]) = 0;
    sp_imag(res->F[i]) = 0;
    scattering_vector_length = sqrt(HKL_list[3*i]*HKL_list[3*i]+HKL_list[3*i+1]*HKL_list[3*i+1]+HKL_list[3*i+2]*HKL_list[3*i+2]);
    for(j = 0;j<ELEMENTS;j++){
      if(is_element_in_molecule[j]){
	scattering_factor_cache[j] = scatt_factor(scattering_vector_length,j,B);
      }
    }
    k = 2*M_PI/exp->wavelength;
    for(j = 0 ;j< mol->natoms;j++){
      if(!mol->atomic_number[j]){
	continue;
      }
      /* Multiply the scattering factor with the ilumination function (should it be the square root of it?)*/
      scattering_factor = scattering_factor_cache[mol->atomic_number[j]]*sqrt(atom_ilumination[j]);
      /*      scattering_factor = 1;*/
      

      sp_real(res->F[i]) += scattering_factor*cos(2*M_PI*(HKL_list[3*i]*mol->pos[j*3]+HKL_list[3*i+1]*mol->pos[j*3+1]+HKL_list[3*i+2]*mol->pos[j*3+2])+
	M_PI/(exp->wavelength*distance)*(mol->pos[j*3]*mol->pos[j*3]+mol->pos[j*3+1]*mol->pos[j*3+1]));
      sp_imag(res->F[i]) += scattering_factor*sin(2*M_PI*(HKL_list[3*i]*mol->pos[j*3]+HKL_list[3*i+1]*mol->pos[j*3+1]+HKL_list[3*i+2]*mol->pos[j*3+2])
						  -M_PI/(exp->wavelength*distance)*(mol->pos[j*3]*mol->pos[j*3]+mol->pos[j*3+1]*mol->pos[j*3+1]));

      /*
      sp_real(res->F[i]) += scattering_factor*cos(k/(2*distance)*((exp->wavelength*HKL_list[3*i]-mol->pos[j*3])*
												     (exp->wavelength*HKL_list[3*i]-mol->pos[j*3])+
												     (exp->wavelength*HKL_list[3*i+1]-mol->pos[j*3+1])*
												     (exp->wavelength*HKL_list[3*i+1]-mol->pos[j*3+1])));
      sp_real(res->F[i]) += scattering_factor*sin(k/(2*distance)*((exp->wavelength*HKL_list[3*i]-mol->pos[j*3])*
												     (exp->wavelength*HKL_list[3*i]-mol->pos[j*3])+
												     (exp->wavelength*HKL_list[3*i+1]-mol->pos[j*3+1])*
												     (exp->wavelength*HKL_list[3*i+1]-mol->pos[j*3+1])));
      */
    }
    res->ints[i] = sp_cabs(res->F[i])*sp_cabs(res->F[i]);
  }
  syncronize_patterns(res);
  return res;
}


double box_fourier_transform(Box box, float h, float k, float l){
  if(box.alpha != 90 ||
     box.beta != 90 ||
     box.gamma != 90){
    sp_error_fatal("Cannot handle non rectangular boxes\n");
  }
  return (sin(PI*box.a*h)/(PI*box.a*h))*(sin(PI*box.b*k)/(PI*box.b*k))*(sin(PI*box.c*l)/(PI*box.c*l));
}


Diffraction_Pattern * compute_box_on_list(Box box, float * HKL_list, int HKL_list_size){  
  int i;
  Diffraction_Pattern * res = malloc(sizeof(Diffraction_Pattern));
  res->F = malloc(sizeof(Complex)*HKL_list_size);
  res->ints = malloc(sizeof(float)*HKL_list_size);
  res->HKL_list = malloc(sizeof(float)*3*HKL_list_size);
  memcpy(res->HKL_list,HKL_list,sizeof(float)*3*HKL_list_size);
  res->HKL_list_size = HKL_list_size;

  for(i = 0;i<HKL_list_size;i++){
    sp_real(res->F[i]) = 0;
    sp_imag(res->F[i]) = 0;
    /*    scattering_vector_length = sqrt(HKL_list[3*i]*HKL_list[3*i]+HKL_list[3*i+1]*HKL_list[3*i+1]+HKL_list[3*i+2]*HKL_list[3*i+2]); */
    sp_real(res->F[i]) = box_fourier_transform(box,HKL_list[3*i],HKL_list[3*i+1],HKL_list[3*i+2]);
    res->ints[i] = sp_cabs(res->F[i])*sp_cabs(res->F[i]);
  }
  return res;
}



Diffraction_Pattern * load_pattern_from_file(CCD * det,char * filename, 
					     float * HKL_list, int HKL_list_size){
  Diffraction_Pattern * res = malloc(sizeof(Diffraction_Pattern));

  res->HKL_list = malloc(sizeof(float)*3*HKL_list_size);
  memcpy(res->HKL_list,HKL_list,sizeof(float)*3*HKL_list_size);
  res->HKL_list_size = HKL_list_size;

  res->ints = read_VTK_to_array(det->nx,det->ny,det->nz,filename);
  res->F = 0;
  return res;
}

void calculate_thomson_correction(CCD * det){
  int x,y,z;
  double px,py;
  int nx = det->nx;
  int ny = det->ny;
  int nz = det->nz;
  double r; /* distance to scatterer */
  int index;
  double r0 = 2.81794e-15; /* classical electron radius = e^2/(m*c^2)*/
  /* For the moment we're considering vertical polarization */
  double polarization_factor = 1;
  det->thomson_correction = malloc(sizeof(float)*nx*ny*nz);
  index = 0;
  for(x = 0;x<det->nx;x++){
    for(y = 0;y<det->ny;y++){
      for(z = 0;z<det->nz;z++){
	px = ((x-(nx-1.0)/2.0)/nx)*det->width/2;
	py = (((ny-1.0)/2.0-y)/ny)*det->height/2;
	r = sqrt(det->distance*det->distance+px*px+py*py);
	det->thomson_correction[index++] = (r0*r0)*polarization_factor;
      }
    }    
  }              
}

void calculate_pixel_solid_angle(CCD * det){
  /* 
     the effective pixel area is calculated by
     the projection of the pixel on a sphere centered
     on the scatterer and going through the center of the pixel
  */
  int i;
  int x,y,z;
  double px,py;
  double corners[4][2];
  double corner_distance[4];
  double projected_pixel_sides[4];
  int nx = det->nx;
  int ny = det->ny;
  int nz = det->nz;
  double r; /* distance to scatterer */
  int index;
  /* For the moment we're considering vertical polarization */
  det->solid_angle = malloc(sizeof(float)*nx*ny*nz);
  index = 0;
  for(x = 0;x<det->nx;x++){
    for(y = 0;y<det->ny;y++){
      for(z = 0;z<det->nz;z++){
	if(det->spherical){
	  r = det->distance;
	  det->solid_angle[index++] = det->pixel_width*det->pixel_height/(r*r);
	}else{
	  px = ((x-(nx-1.0)/2.0)/nx)*det->width/2;
	  py = (((ny-1.0)/2.0-y)/ny)*det->height/2;
	  r = sqrt(det->distance*det->distance+px*px+py*py);
	  /* top left */
	  corners[0][0] = px-det->pixel_width/2;
	  corners[0][1] = py+det->pixel_height/2;
	  
	  corner_distance[0] = sqrt(det->distance*det->distance+corners[0][0]*corners[0][0]+corners[0][1]*corners[0][1]);
	  /* top right */
	  corners[1][0] = px+det->pixel_width/2;
	  corners[1][1] = py+det->pixel_height/2;
	  corner_distance[1] = sqrt(det->distance*det->distance+corners[1][0]*corners[1][0]+corners[1][1]*corners[1][1]);
	  /* bottom right */
	  corners[2][0] = px+det->pixel_width/2;
	  corners[2][1] = py-det->pixel_height/2;
	  corner_distance[2] = sqrt(det->distance*det->distance+corners[2][0]*corners[2][0]+corners[2][1]*corners[2][1]);
	  /* bottom left */
	  corners[3][0] = px-det->pixel_width/2;
	  corners[3][1] = py-det->pixel_height/2;
	  corner_distance[3] = sqrt(det->distance*det->distance+corners[3][0]*corners[3][0]+corners[3][1]*corners[3][1]);
	  /* project on plane*/
	  if(!det->spherical){
	    for(i = 0;i<4;i++){
	      corners[i][0] *= r/corner_distance[i];
	      corners[i][1] *= r/corner_distance[i];	
	    }
	  }
	  /* top */
	  projected_pixel_sides[0] = sqrt((corners[0][0]-corners[1][0])*(corners[0][0]-corners[1][0])+(corners[0][1]-corners[1][1])*(corners[0][1]-corners[1][1]));
	  /* left */
	  projected_pixel_sides[1] = sqrt((corners[0][0]-corners[3][0])*(corners[0][0]-corners[3][0])+(corners[0][1]-corners[3][1])*(corners[0][1]-corners[3][1]));
	  det->solid_angle[index++] = projected_pixel_sides[0]*projected_pixel_sides[1]/(r*r);
	}
      }
    }
  }
}                    


void calculate_photons_per_pixel(Diffraction_Pattern * pattern, Options * opts){
  CCD * det = opts->detector;
  Experiment * experiment = opts->experiment;
  int i;
  if(!det->thomson_correction){
    calculate_thomson_correction(det);
  }
  if(!det->solid_angle){
    calculate_pixel_solid_angle(det);
  }
  det->photons_per_pixel = malloc(sizeof(float)*det->nx*det->ny*det->nz);
  for(i = 0;i<det->nx*det->ny*det->nz;i++){
    det->photons_per_pixel[i] = det->thomson_correction[i] * det->solid_angle[i] * pattern->ints[i] * experiment->beam_intensity;
  }
}


void write_hkl_grid(float * list, Molecule * mol,CCD * det){
  Image * hkl_grid;
  int x,y,z;
  int i;
  int index;
  double min_x = FLT_MAX;
  double min_y = FLT_MAX; 
  double min_z = FLT_MAX; 
  double max_x = -FLT_MAX;
  double max_y = -FLT_MAX;
  double max_z = -FLT_MAX;
  double max_dim;
#ifdef MPI
  if(!is_mpi_master()){
    return;
  }
#endif
  for(i =0 ;i<mol->natoms;i++){
    if(mol->pos[i*3] < min_x){
      min_x = mol->pos[i*3];
    }
    if(mol->pos[i*3] > max_x){
      max_x = mol->pos[i*3];
    }
    if(mol->pos[i*3+1] < min_y){
      min_y = mol->pos[i*3+1];
    }
    if(mol->pos[i*3+1] > max_y){
      max_y = mol->pos[i*3+1];
    }
    if(mol->pos[i*3+2] < min_z){
      min_z = mol->pos[i*3+2];
    }
    if(mol->pos[i*3+2] > max_z){
      max_z = mol->pos[i*3+2];
    }
  }
  
  max_dim = sp_max(max_x-min_x,max_y-min_y);
  max_dim = sp_max(max_z-min_z,max_dim);
  hkl_grid = sp_image_alloc(det->nx/det->binning_x,det->ny/det->binning_y,det->nz/det->binning_z);
  hkl_grid->detector->image_center[0] = sp_image_x(hkl_grid)/2;
  hkl_grid->detector->image_center[1] = sp_image_y(hkl_grid)/2;
  hkl_grid->detector->image_center[2] = sp_image_z(hkl_grid)/2;
  i = 0;
  index = 0;
  for(x = 0;x<sp_image_x(hkl_grid);x++){
    for(y = 0;y<sp_image_y(hkl_grid);y++){
      for(z = 0;z<sp_image_z(hkl_grid);z++){
	if(fabs(x-hkl_grid->detector->image_center[0]) > fabs(y-hkl_grid->detector->image_center[1])){
	  sp_real(hkl_grid->image->data[i]) = list[index]*max_dim;
	}else{
	  sp_real(hkl_grid->image->data[i]) = list[index+1]*max_dim;
	}
	hkl_grid->mask->data[i] = 1;
	index += 3*det->binning_x;
	i++;
      }
    }
    index += 3*(det->binning_x-1)*(det->ny);
  }
  sp_image_write(hkl_grid,"hkl_grid.h5",sizeof(real));
  sp_image_free(hkl_grid);
}
