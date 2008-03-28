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



#include "diffraction.h"

#ifdef MPI

#include "mpi.h"

#endif

void get_my_loop_start_and_end(int size, int * start, int * end){
#ifndef MPI
  *start = 0;
  *end = size;
#else
  int np;
  int myid;
  MPI_Comm_size(MPI_COMM_WORLD,&np);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  if(myid != np-1){
    *start = myid*size/np;
    *end = (myid+1)*size/np;
  }else{
    /* last node takes a bit of extra work due to rounding */
    *start = myid*size/np;
    *end = size;
  } 
#endif
}

void syncronize_patterns(Diffraction_Pattern * pat){
#ifndef MPI
  return;
#else
  if(is_mpi_master()){
    mpi_receive_pattern(pat);
  }else{
    mpi_send_pattern(pat);
  }
#endif
}

void sum_patterns(Diffraction_Pattern * pat){
#ifndef MPI
  return;
#else
  if(is_mpi_master()){
    mpi_sum_receive_pattern(pat);
  }else{
    mpi_sum_send_pattern(pat);
  }
#endif
}


#ifdef MPI

int is_mpi_master(){
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  return (myid == 0);
}

void get_id_loop_start_and_end(int id, int size, int * start, int * end){
  int np;
  MPI_Comm_size(MPI_COMM_WORLD,&np);
  if(id != np-1){
    *start = id*size/np;
    *end = (id+1)*size/np;
  }else{
    /* last node takes a bit of extra work due to rounding */
    *start = id*size/np;
    *end = size;
  } 
}

int mpi_receive_pattern(Diffraction_Pattern * pat){
  int np;
  int i,j;
  int start,end;
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD,&np);
  MPI_Datatype mpi_datatype;
  if(sizeof(Complex) == 8){
    mpi_datatype = MPI_FLOAT;
  }else if(sizeof(Complex) == 16){
    mpi_datatype = MPI_DOUBLE;
  }
  for(i = 1;i<np;i++){
    get_id_loop_start_and_end(i,pat->HKL_list_size,&start,&end);
    /*count is (end-start)*2 because we're sending fftw_complex not doubles */
    MPI_Recv(&(pat->F[start]),(end-start)*2,mpi_datatype,i,0,MPI_COMM_WORLD,&status);
    for(j = start;j<end;j++){       
      pat->ints[j] = sp_cabs(pat->F[j])*sp_cabs(pat->F[j]);
    }
  }  
  return 0;
}

int mpi_sum_receive_pattern(Diffraction_Pattern * pat){
  int np;
  int i,j;
  int start,end;
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD,&np);
  Complex * buffer = malloc(sizeof(Complex)*pat->HKL_list_size);
  MPI_Datatype mpi_datatype;
  if(sizeof(Complex) == 8){
    mpi_datatype = MPI_FLOAT;
  }else if(sizeof(Complex) == 16){
    mpi_datatype = MPI_DOUBLE;
  }
  for(i = 1;i<np;i++){
    MPI_Recv(buffer,pat->HKL_list_size*2,mpi_datatype,i,0,MPI_COMM_WORLD,&status);
    for(j = 0;j<pat->HKL_list_size;j++){       
      sp_real(pat->F[j]) += sp_real(buffer[j]);
      sp_imag(pat->F[j]) += sp_imag(buffer[j]);
    }
  }  
  free(buffer);
  for(j = 0;j<pat->HKL_list_size;j++){       
    pat->ints[j] = sp_cabs(pat->F[j])*sp_cabs(pat->F[j]);
  }
  return 0;
}

int mpi_send_pattern(Diffraction_Pattern * pat){
  int start,end;
  get_my_loop_start_and_end(pat->HKL_list_size,&start,&end);
  /*count is (end-start)*2 because we're sending fftw_complex not doubles */
  if(sizeof(Complex) == 8){
    MPI_Send (&(pat->F[start]),(end-start)*2,MPI_FLOAT,0,0,MPI_COMM_WORLD);
  }else if(sizeof(Complex) == 16){
    MPI_Send (&(pat->F[start]),(end-start)*2,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
  }
  return 0;    
}

int mpi_sum_send_pattern(Diffraction_Pattern * pat){
  int start,end;
  get_my_loop_start_and_end(pat->HKL_list_size,&start,&end);
  /*count is (end-start)*2 because we're sending fftw_complex not doubles */
  if(sizeof(Complex) == 8){
    MPI_Send (pat->F,pat->HKL_list_size*2,MPI_FLOAT,0,0,MPI_COMM_WORLD);
  }else if(sizeof(Complex) == 16){
    MPI_Send (pat->F,pat->HKL_list_size*2,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
  }
  return 0;    
}
#endif
