#include "crystal.h"
#include "cufft.h"


__global__ void CUDA_add_shifted_cell(cufftComplex * cF, float * I,cufftComplex * F, const float * HKL_list,const int hkl_size,const float dx,const float dy, const float dz);


void calculate_pattern_from_crystal_cuda(float * d_I, cufftComplex * d_F,float * d_HKL_list, int HKL_list_size,Options * opts){
  dim3 threads_per_block(8,8);
  const int HKL_side = (sqrt(HKL_list_size)+1);
  dim3 number_of_blocks( (HKL_side+threads_per_block.x-1)/threads_per_block.x,
			 (HKL_side+threads_per_block.y-1)/threads_per_block.y );

  //  int threads_per_block = 64;
  //  int number_of_blocks = (HKL_list_size+threads_per_block-1)/threads_per_block;
  float cell[9];
  cufftComplex * d_cF;
  cutilSafeCall(cudaMalloc((void **)&d_cF,sizeof(cufftComplex)*HKL_list_size));
  cutilSafeCall(cudaMemset(d_cF,0,sizeof(cufftComplex)*HKL_list_size));
  crystal_cell_matrix(opts,cell);
  float dr[3];
  int k =0;
  for(int a = 0;a<opts->crystal_size[0];a++){
    for(int b = 0;b<opts->crystal_size[1];b++){
      for(int c = 0;c<opts->crystal_size[2];c++){
	dr[0] = cell[0]*a+cell[3]*b+cell[6]*c;
	dr[1] = cell[1]*a+cell[4]*b+cell[7]*c;
	dr[2] = cell[2]*a+cell[5]*b+cell[8]*c;
	printf("%f%% done\n",
	       100.0*k/(opts->crystal_size[0]*opts->crystal_size[1]*opts->crystal_size[2]));
	k++;
	CUDA_add_shifted_cell<<<number_of_blocks, threads_per_block>>>(d_cF,d_I,d_F,d_HKL_list,HKL_list_size,dr[0],dr[1],dr[2]);
	cudaThreadSynchronize();
      }
    }
  }
  cudaFree(d_cF);
}

__global__ void CUDA_add_shifted_cell(cufftComplex * cF, float * I,cufftComplex * F, const float * HKL_list,const int hkl_size,const float dx,const float dy, const float dz){
  const int i = ((blockIdx.y*blockDim.y + threadIdx.y)*gridDim.x + blockIdx.x )*blockDim.x + threadIdx.x;
  //  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i<hkl_size){
    /* use atan2 instead of carg */
    float phi = 2*3.14159265F*(HKL_list[3*i]*-dx+HKL_list[3*i+1]*-dy+HKL_list[3*i+2]*-dz) + atan2(F[i].y,F[i].x);
    float amp = sqrt(F[i].x*F[i].x+F[i].y*F[i].y);
    cF[i].x += cos(phi)*amp;
    cF[i].y += sin(phi)*amp;
    I[i] = cF[i].x*cF[i].x+cF[i].y*cF[i].y;
  }
}
