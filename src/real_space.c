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
#include <spimage.h>
#include <time.h>
#include "config.h"
#include "diffraction.h"
#include "mpi_comm.h"
#include "io.h"

int descend_complex_compare(const void * pa,const void * pb);

Image * calculate_noiseless_real_space(Options * opts, Diffraction_Pattern * pattern){
  Image * real_space = sp_image_alloc(opts->detector->nx/opts->detector->binning_x,
				      opts->detector->ny/opts->detector->binning_y,
				      opts->detector->nz/opts->detector->binning_z);
  real_space->phased = 1;
  real_space->shifted = 0;
  int i = 0;
  for(int z = 0;z<sp_image_z(real_space);z++){
    for(int y = 0;y<sp_image_y(real_space);y++){
      for(int x = 0;x<sp_image_x(real_space);x++){
	Complex phased = pattern->F[(x*opts->detector->binning_x) + 
				    (y*opts->detector->binning_y) * opts->detector->nx +
				    (z*opts->detector->binning_z) * opts->detector->nx * opts->detector->ny];
/*	sp_cscale(phased,1.0/sp_cabs(phased));
	sp_cscale(phased,opts->detector->noiseless_output[i++]);*/
	if(sp_cabs(phased)){
	  sp_image_set(real_space,x,y,z,sp_cscale(phased,sqrt(opts->detector->noiseless_output[i])/sp_cabs(phased)));
	}else{
	  sp_image_set(real_space,x,y,z,phased);
	}
	sp_i3matrix_set(real_space->mask,x,y,z,1);
	i++;
      }
    }
  }
  sp_image_write(real_space,"unshifted_pattern.vtk",0);
  real_space->detector->image_center[0] = (sp_image_x(real_space))/2;
  real_space->detector->image_center[1] = (sp_image_y(real_space))/2;
  real_space->detector->image_center[2] = (sp_image_z(real_space))/2;
  Image * rs = sp_image_shift(real_space);
  sp_image_free(real_space);
  sp_image_write(rs,"shifted_pattern.vtk",0);
  real_space = sp_image_ifft(rs);
  sp_image_free(rs);
  /* Move mass to the middle of the image */
  sp_image_write(real_space,"real_space_uncentered.vtk",0);
  sp_image_write(real_space,"real_space_uncentered.h5",0);
  sp_vector * center = sp_image_center_of_mass(real_space);
  /*FM: This is buggy the center of mass should take wrap around into account*/
  sp_image_translate(real_space,
		     sp_image_x(real_space)/2-center->data[0],
		     sp_image_y(real_space)/2-center->data[1],
		     sp_image_z(real_space)/2-center->data[2],
		     SP_TRANSLATE_WRAP_AROUND);
  center = sp_image_center_of_mass(real_space);
  sp_image_translate(real_space,
		     sp_image_x(real_space)/2-center->data[0],
		     sp_image_y(real_space)/2-center->data[1],
		     sp_image_z(real_space)/2-center->data[2],
		     SP_TRANSLATE_WRAP_AROUND);
  sp_image_write(real_space,"real_space_centered.vtk",0);
  rs = real_space;
/*  real_space->shifted = 1;
  real_space->detector->image_center[0] = 0;
  real_space->detector->image_center[1] = 0;
  real_space->detector->image_center[2] = 0;
  rs = sp_image_shift(real_space);*/
  for(int x = 0;x<sp_image_x(real_space);x++){
    for(int y = 0;y<sp_image_y(real_space);y++){
      for(int z = 0;z<sp_image_z(real_space);z++){
	real size = (sp_image_size(rs));
	sp_image_set(rs,x,y,z,sp_cinit(sp_cabs(sp_image_get(rs,x,y,z))/size,0));
	sp_i3matrix_set(rs->mask,x,y,z,1);
      }
    }
  }
  return rs;
}

void write_density_histogram(Image * a){
  qsort(a->image->data,sp_image_size(a),sizeof(Complex),descend_complex_compare);  
  FILE * fp = fopen("density_histogram.log","w");
  fprintf(fp,"# Spsim version %s compiled on %s at %s\n",VERSION,__DATE__,__TIME__);
  time_t t ;
  fprintf(fp,"# Date: %s\n",ctime(&t));
  fprintf(fp,"@TYPE xy\n");
  fprintf(fp,"@ view 0.15, 0.15, 0.75, 0.8\n");
  fprintf(fp,"@ legend on\n");
  fprintf(fp,"@ legend box on\n");
  fprintf(fp,"@ legend loctype view\n");
  fprintf(fp,"@ legend 0.78, 0.8\n");
  fprintf(fp,"@ legend length 2\n");
  fprintf(fp,"@ s1 legend \"Density \\N\"\n");
  for(int i = 0;i<sp_image_size(a);i++){  
    fprintf(fp,"%f\t%f\n",(real)i/sp_image_size(a),sp_cabs(a->image->data[i]));
  }
  fclose(fp);
}


int descend_complex_compare(const void * pa,const void * pb){
  Complex a,b;
  a = *((Complex *)pa);
  b = *((Complex *)pb);
  if(sp_cabs(a) < sp_cabs(b)){
    return 1;
  }else if(sp_cabs(a) == sp_cabs(b)){
    return 0;
  }else{
    return -1;
  }
}

