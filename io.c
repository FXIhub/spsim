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
#include "config.h"
#include "diffraction.h"

#include "io.h"
#include "mpi.h"

void write_2D_array_to_vtk(float * F,int nx, int ny, char * filename){
  FILE * f;
  int i;
#ifdef MPI  
  if(!is_mpi_master()){
    return;
  }
#endif
  f = fopen(filename,"w");
  if(!f){
    perror("Bad file in write_vtk!");
    abort();
  }
  fprintf(f,"# vtk DataFile Version 2.0\n");
  fprintf(f,"Generated by image_util write_vtk()\n");
  fprintf(f,"ASCII\n");
  fprintf(f,"DATASET STRUCTURED_POINTS\n");
  fprintf(f,"DIMENSIONS %d %d 1\n",ny,nx);
  fprintf(f,"ORIGIN 0 0 0\n");
  fprintf(f,"SPACING 1 1 1\n");
  fprintf(f,"POINT_DATA %d\n",nx*ny);
  fprintf(f,"SCALARS intensity float 1\n");
  fprintf(f,"LOOKUP_TABLE default\n");
  fprintf(f,"%6e",F[0]);
  for(i = 1;i<nx*ny;i++){
    fprintf(f," %e",F[i]);
  }
  fprintf(f,"\n");
  fflush(f);
  fclose(f);
}


float * read_VTK_to_array(int nx, int ny, char * filename){
  FILE * fp;
  int i;
  int buffer_size = 1024;
  char * buffer = malloc(sizeof(char)*buffer_size);
  char * start = buffer;
  char * end;
  float * res = malloc(sizeof(float)*nx*ny);
  fp = fopen(filename,"r");
  if(!fp){
    perror("Bad file in read_vtk!");
    abort();
  }
  for(i = 0;i<10;i++){
    fgets(buffer,1024,fp);    
  }
  buffer[buffer_size-2] = 1;
  while(buffer[buffer_size-2]){
    buffer[buffer_size-2] = 0;
    fgets(start,buffer_size-(start-buffer),fp);
    if(buffer[buffer_size-2]){
      buffer_size *= 2;
      buffer = realloc(buffer,sizeof(char)*buffer_size);
      start = &buffer[buffer_size/2-1];
      buffer[buffer_size-2] = 1;
    }
  }
  start = NULL;
  end = buffer;
  i = 0;
  while(start != end){
    start = end;
    res[i++] = strtod(start,&end);
  }
  fclose(fp);
  return res;
}
