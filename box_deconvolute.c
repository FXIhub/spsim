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
#include <spimage.h>
#include "config.h"
#include "diffraction.h"
#include "molecule.h"
#include "io.h"
#include "mpi.h"
#include "noise.h"
#include "amplification.h"
#include "real_space.h"


typedef struct{
  double a;
  double b;
  double c;
  double alpha;
  double beta;
  double gamma;   
}Box;

 pattern = compute_box_on_list(box,HKL_list,HKL_list_size);

int main(int argc, char ** argv){
  Options * opts = set_defaults();
  Diffraction_Pattern * pattern;
  float * HKL_list;
  int HKL_list_size = 0;
  int i;
  int is3D = 0;
  Image * input;
  Image * output;
  Box box;
  read_options_file("spsim.conf",opts);
  if(opts->detector->nz > 1){
    is3D = 1;
  }else{
    is3D = 0;
  }
  if(argc < 4){
    sp_error_fatal("Usage: box_deconvolute <image> <a b c> [alpha beta gamma]\n");
  }
  if(argc > 4){
    sp_error_fatal("Sorry non rectangular boxes are not supported yet\n");
  }

  input = sp_image_read(argv[1]);

  if(is3D){
    HKL_list = get_HKL_list_for_3d_detector(opts->detector,opts->experiment, &HKL_list_size);
    apply_orientation_to_HKL_list(&HKL_list,&HKL_list_size,opts);
  }else{
    HKL_list = get_HKL_list_for_detector(opts->detector,opts->experiment, &HKL_list_size);
    apply_orientation_to_HKL_list(&HKL_list,&HKL_list_size,opts);
  }
  pattern = compute_box_on_list(box,HKL_list,HKL_list_size);
  
  gaussian_blur_pattern(opts,pattern);

  write_3D_array_to_vtk(pattern->ints,opts->detector->nx,opts->detector->ny,opts->detector->nz,
			"box_scattering_factor.vtk");  
  write_3D_array_as_structured_grid_to_vtk(pattern->ints,opts->detector,HKL_list,HKL_list_size,opts->n_patterns,"ewald.vtk");  
  calculate_thomson_correction(opts->detector);
  calculate_pixel_solid_angle(opts->detector);
  write_3D_array_to_vtk(opts->detector->thomson_correction,opts->detector->nx,opts->detector->ny,
			opts->detector->nz,"box_thomson_correction.vtk");
  write_3D_array_to_vtk(opts->detector->solid_angle,opts->detector->nx,opts->detector->ny,
			opts->detector->nz,"box_solid_angle.vtk");
  calculate_photons_per_pixel(pattern,opts->detector,opts->experiment);
  write_3D_array_to_vtk(opts->detector->photons_per_pixel,opts->detector->nx,opts->detector->ny,
			opts->detector->nz,"box_pattern.vtk");
  generate_poisson_noise(opts->detector);
  write_3D_array_to_vtk(opts->detector->photon_count,opts->detector->nx,opts->detector->ny,
			opts->detector->nz,"box_photon_count.vtk");
  calculate_electrons_per_pixel(opts->detector,opts->experiment);
  write_3D_array_to_vtk(opts->detector->electrons_per_pixel,opts->detector->nx,opts->detector->ny,opts->detector->nz,"box_electrons_per_pixel.vtk");
  calculate_real_detector_output(opts->detector,opts->experiment);
  output = sp_image_alloc(opts->detector->nx/opts->detector->binning_x,opts->detector->ny/opts->detector->binning_y,opts->detector->nz/opts->detector->binning_z);
  i = 0;
  for(int x = 0;x<sp_image_x(output);x++){
    for(int y = 0;y<sp_image_y(output);y++){
      for(int z = 0;z<sp_image_z(output);z++){
	sp_image_set(output,x,y,z,sp_cinit(opts->detector->real_output[i++],0));
	sp_i3matrix_set(output->mask,x,y,z,1);
      }
    }
  }
  sp_image_write(output,"box_real_output.h5",sizeof(real));
  sp_image_write(output,"box_real_output.vtk",0);

  calculate_noiseless_detector_output(opts->detector,opts->experiment);
  noiseless = sp_image_alloc(opts->detector->nx/opts->detector->binning_x,opts->detector->ny/opts->detector->binning_y,
			     opts->detector->nz/opts->detector->binning_z);
  i = 0;
  for(int x = 0;x<sp_image_x(noiseless);x++){
    for(int y = 0;y<sp_image_y(noiseless);y++){
      for(int z = 0;z<sp_image_z(noiseless);z++){
	sp_image_set(noiseless,x,y,z,sp_cinit(opts->detector->noiseless_output[i++],0));
	sp_i3matrix_set(noiseless->mask,x,y,z,1);
      }
    }
  }
  sp_image_write(noiseless,"box_noiseless_output.h5",sizeof(real));
  sp_image_write(noiseless,"box_noiseless_output.vtk",0);
  


  if(pattern->F){
    Image * rs = calculate_noiseless_real_space(opts,pattern);
    
    sp_image_write(rs,"box_real_space.h5",sizeof(real));
    sp_image_write(rs,"box_real_space.vtk",0);

    write_density_histogram(rs);
  }
  return 0;
}
