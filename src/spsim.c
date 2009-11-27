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
#include "mpi_comm.h"
#include "noise.h"
#include "amplification.h"
#include "real_space.h"

void gaussian_blur_real_space(Options * opts,Diffraction_Pattern * pattern){
  Image * real_space = calculate_noiseless_real_space(opts,pattern);
  float radius = opts->detector->real_space_blurring;
  //  float d_c = sqrt(opts->detector->nx*opts->detector->nx/4+opts->detector->ny*opts->detector->ny/4+opts->detector->nz*opts->detector->nz/4);
  if(!radius){
    return;
  }
  /*  radius *= d_c; */
  /*  f(x,y) = 1/sqrt(2*M_PI*radius) * exp(-(x^2+y^2+dz^2)/(2*radius^2)) */
  
  int i = 0;
  for(int x = 0;x<opts->detector->nx;x++){
    float dx = fabs(x-((opts->detector->nx-1)/2));
    for(int y = 0;y<opts->detector->ny;y++){
    float dy = fabs(y-((opts->detector->ny-1)/2));
      for(int z = 0;z<opts->detector->nz;z++){
	float dz = fabs(z-((opts->detector->nz-1)/2));
	float factor = 1/sqrt(2*M_PI*radius) * exp(-(dx*dx+dy*dy+dz*dz)/(2*radius*radius));
	sp_real(real_space->image->data[i]) *= factor;
	sp_imag(real_space->image->data[i]) *= factor;
	i++;
      }
    }
  }
  sp_image_write(real_space,"blurred_real_space.vtk",0);
  Image * tmp = sp_image_shift(sp_image_fft(real_space));
  sp_image_write(tmp,"blurred_real_space_pattern.vtk",0);
  for(i = 0;i<sp_image_size(real_space);i++){
    pattern->F[i] = tmp->image->data[i];
    pattern->ints[i] = sp_cabs(pattern->F[i])*sp_cabs(pattern->F[i]);
  }
}

void gaussian_blur_pattern(Options * opts,Diffraction_Pattern * pattern){
  if(!opts->detector->gaussian_blurring){
    return;
  }
  float radius = opts->detector->gaussian_blurring;
  float d_c = sqrt(opts->detector->nx*opts->detector->nx/4+opts->detector->ny*opts->detector->ny/4+opts->detector->nz*opts->detector->nz/4);
  radius *= d_c;
  /*  f(x,y) = 1/sqrt(2*M_PI*radius) * exp(-(x^2+y^2+dz^2)/(2*radius^2)) */
  
  int i = 0;
  for(int x = 0;x<opts->detector->nx;x++){
    float dx = fabs(x-((opts->detector->nx-1)/2));
    for(int y = 0;y<opts->detector->ny;y++){
    float dy = fabs(y-((opts->detector->ny-1)/2));
      for(int z = 0;z<opts->detector->nz;z++){
	float dz = fabs(z-((opts->detector->nz-1)/2));
	float factor = 1/sqrt(2*M_PI*radius) * exp(-(dx*dx+dy*dy+dz*dz)/(2*radius*radius));
	sp_real(pattern->F[i]) *= factor;
	sp_imag(pattern->F[i]) *= factor;
	pattern->ints[i] = sp_cabs(pattern->F[i])*sp_cabs(pattern->F[i]);
	i++;
      }
    }
  }
}

int main(int argc, char ** argv){
  Options * opts = set_defaults();
  Diffraction_Pattern * pattern;
  float * HKL_list;
  int HKL_list_size = 0;
  int i;
  int is3D = 0;
  Image * noiseless;
  Image * output;
  Molecule * mol = NULL;
#ifdef MPI
  MPI_Init(&argc, &argv);
#endif
  read_options_file("spsim.conf",opts);
  write_options_file("spsim.confout", opts);
  if(opts->detector->nz > 1){
    is3D = 1;
  }else{
    is3D = 0;
  }
  if(opts->input_type == CHEM_FORMULA){
    mol = get_Molecule_from_formula(opts->chem_formula,opts);
    write_pdb_from_mol("molout.pdb",mol);
  }else if(opts->input_type == PDB){
    mol = get_Molecule_from_pdb(opts->pdb_filename);
  }else{
    fprintf(stderr,"input_type not specified!\n");
    exit(0);
  }

  if(opts->sf_filename[0]){
    pattern = load_pattern_from_file(opts->detector,opts->sf_filename,HKL_list,HKL_list_size);
  }else{
    if(is3D){
      HKL_list = get_HKL_list_for_3d_detector(opts->detector,opts->experiment, &HKL_list_size);
      apply_orientation_to_HKL_list(&HKL_list,&HKL_list_size,opts);
      if(opts->use_fft_for_sf){
	pattern = compute_pattern_by_fft(mol,opts->detector,opts->experiment,opts->b_factor);
      }else if(opts->use_nfft_for_sf){
	pattern = compute_pattern_by_nfft(mol,opts->detector,opts->experiment,opts->b_factor,HKL_list,opts);
      }else{
#ifdef _USE_CUDA
	if(opts->use_cuda){
	  pattern = cuda_compute_pattern_on_list(mol,HKL_list,HKL_list_size,opts->b_factor,opts->experiment,opts);	
	}else
#endif
	if(opts->vectorize){
	  pattern = vector_compute_pattern_on_list(mol,HKL_list,HKL_list_size,opts->b_factor,opts->experiment,opts);	
	}else{
	  pattern = compute_pattern_on_list(mol,HKL_list,HKL_list_size,opts->b_factor,opts->experiment,opts);
	}
      }
    }else{
      HKL_list = get_HKL_list_for_detector(opts->detector,opts->experiment, &HKL_list_size);
      apply_orientation_to_HKL_list(&HKL_list,&HKL_list_size,opts);
      if(opts->use_fft_for_sf){
	fprintf(stderr,"Cannot use fft for 2D pattern calculation!\n");
	abort();
      }else if(opts->use_nfft_for_sf){
	pattern = compute_pattern_on_list_by_nfft(mol,HKL_list,HKL_list_size,opts->detector,opts->b_factor,opts);
      }else{
#ifdef _USE_CUDA
	if(opts->use_cuda){
	  pattern = cuda_compute_pattern_on_list(mol,HKL_list,HKL_list_size,opts->b_factor,opts->experiment,opts);	
	}else
#endif
	if(opts->vectorize){
	  pattern = vector_compute_pattern_on_list(mol,HKL_list,HKL_list_size,opts->b_factor,opts->experiment,opts);	
	}else{
	  pattern = compute_pattern_on_list(mol,HKL_list,HKL_list_size,opts->b_factor,opts->experiment,opts);
	}
      }
    }
  }
#ifdef MPI
  if(!is_mpi_master()){
    MPI_Finalize();
    return 0;
  }
#endif

  gaussian_blur_pattern(opts,pattern);
  //  gaussian_blur_real_space(opts,pattern);

  Image * scattering_amp = sp_image_alloc(opts->detector->nx,opts->detector->ny,opts->detector->nz);

  i = 0;
  for(int x = 0;x<sp_image_x(scattering_amp);x++){
    for(int y = 0;y<sp_image_y(scattering_amp);y++){
      for(int z = 0;z<sp_image_z(scattering_amp);z++){
	sp_image_set(scattering_amp,x,y,z,pattern->F[i++]);
	sp_i3matrix_set(scattering_amp->mask,x,y,z,1);
      }
    }
  }

  Image * scattering_int = sp_image_alloc(opts->detector->nx,opts->detector->ny,opts->detector->nz);

  i = 0;
  for(int x = 0;x<sp_image_x(scattering_int);x++){
    for(int y = 0;y<sp_image_y(scattering_int);y++){
      for(int z = 0;z<sp_image_z(scattering_int);z++){
	sp_image_set(scattering_int,x,y,z,sp_cinit(pattern->ints[i++],0));
	sp_i3matrix_set(scattering_int->mask,x,y,z,1);
      }
    }
  }

  scattering_amp->phased = 1;
  /*  sp_image_write(scattering_amp,"scattering_amp.h5",0);*/
  sp_image_write(scattering_int,"scattering_int.h5",0);
  if(opts->fast_exit){
#ifdef MPI
    MPI_Finalize();
#endif
    return 0;
  }
  write_3D_array_to_vtk(pattern->ints,opts->detector->nx,opts->detector->ny,opts->detector->nz,"scattering_factor.vtk");  

  write_3D_array_as_structured_grid_to_vtk(pattern->ints,opts->detector,HKL_list,HKL_list_size,opts->n_patterns,"ewald.vtk");  
  calculate_thomson_correction(opts->detector);
  calculate_pixel_solid_angle(opts->detector);
  /*  write_3D_array_to_vtk(opts->detector->thomson_correction,opts->detector->nx,opts->detector->ny,
			opts->detector->nz,"thomson_correction.vtk");*/
  /*  write_3D_array_to_vtk(opts->detector->solid_angle,opts->detector->nx,opts->detector->ny,
      opts->detector->nz,"solid_angle.vtk");*/
  calculate_photons_per_pixel(pattern,opts->detector,opts->experiment);
  /*  write_3D_array_to_vtk(opts->detector->photons_per_pixel,opts->detector->nx,opts->detector->ny,
      opts->detector->nz,"pattern.vtk");*/
  //generate_gaussian_noise(opts->detector);
  generate_poisson_noise(opts->detector);
  /*  write_3D_array_to_vtk(opts->detector->photon_count,opts->detector->nx,opts->detector->ny,
      opts->detector->nz,"photon_count.vtk");*/

  output = sp_image_alloc(opts->detector->nx/opts->detector->binning_x,opts->detector->ny/opts->detector->binning_y,opts->detector->nz/opts->detector->binning_z);

  i = 0;
  for(int x = 0;x<sp_image_x(output);x++){
    for(int y = 0;y<sp_image_y(output);y++){
      for(int z = 0;z<sp_image_z(output);z++){
	sp_image_set(output,x,y,z,sp_cinit(opts->detector->photon_count[i++],0));
	sp_i3matrix_set(output->mask,x,y,z,1);
      }
    }
  }
  output->detector->wavelength = opts->experiment->wavelength;
  output->detector->pixel_size[0] = opts->detector->pixel_width*opts->detector->binning_x;
  output->detector->pixel_size[1] = opts->detector->pixel_height*opts->detector->binning_y;
  output->detector->detector_distance = opts->detector->distance;
  sp_image_write(output,"photon_out.h5",sizeof(real));
  /*      sp_image_write(output,"photon_out.vtk",0);*/

  calculate_electrons_per_pixel(opts->detector,opts->experiment);
  write_3D_array_to_vtk(opts->detector->electrons_per_pixel,opts->detector->nx,opts->detector->ny,opts->detector->nz,"electrons_per_pixel.vtk");
  calculate_real_detector_output(opts->detector,opts->experiment);

  i = 0;
  for(int x = 0;x<sp_image_x(output);x++){
    for(int y = 0;y<sp_image_y(output);y++){
      for(int z = 0;z<sp_image_z(output);z++){
	sp_image_set(output,x,y,z,sp_cinit(opts->detector->real_output[i++],0));
	sp_i3matrix_set(output->mask,x,y,z,1);
      }
    }
  }
  output->detector->wavelength = opts->experiment->wavelength;
  output->detector->pixel_size[0] = opts->detector->pixel_width*opts->detector->binning_x;
  output->detector->pixel_size[1] = opts->detector->pixel_height*opts->detector->binning_y;
  output->detector->detector_distance = opts->detector->distance;

  /*  sp_image_write(output,"real_output.h5",sizeof(real));
      sp_image_write(output,"real_output.vtk",0);*/


  calculate_noiseless_detector_output(opts->detector,opts->experiment);
  noiseless = sp_image_alloc(opts->detector->nx/opts->detector->binning_x,opts->detector->ny/opts->detector->binning_y,
			     opts->detector->nz/opts->detector->binning_z);
  noiseless->detector->wavelength = opts->experiment->wavelength;
  noiseless->detector->pixel_size[0] = opts->detector->pixel_width*opts->detector->binning_x;
  noiseless->detector->pixel_size[1] = opts->detector->pixel_height*opts->detector->binning_y;
  noiseless->detector->detector_distance = opts->detector->distance;
  i = 0;
  //  for(int u = 0; u < opts->n_patterns; u++) {
  //printf("pattern %i\n",u);
  for(int x = 0;x<sp_image_x(noiseless);x++){
    for(int y = 0;y<sp_image_y(noiseless);y++){
      for(int z = 0;z<sp_image_z(noiseless);z++){
	//printf("%i %i %i\n",x,y,z);
	sp_image_set(noiseless,x,y,z,sp_cinit(opts->detector->noiseless_output[i++],0));
	//printf("set image to %g\n",opts->detector->noiseless_output[i-1]);
	sp_i3matrix_set(noiseless->mask,x,y,z,1);
	//printf("set mask\n");
      }
    }
  }
  //char buffer_u[1024];
  //printf("init buffer\n");
  //sprintf(buffer_u,"noiseless%i.h5",u);
  //sp_image_write(noiseless,buffer_u,sizeof(real));
  //printf("write h5\n");
  //sprintf(buffer_u,"noiseless%i.vtk",u);
  //sp_image_write(noiseless,buffer_u,0);
  //printf("write vtk\n");
  sp_image_write(noiseless,"noiseless_output.h5",sizeof(real));
  //  sp_image_write(noiseless,"noiseless_output.vtk",0);
  //}
  

  if(pattern->F){
    Image * rs = calculate_noiseless_real_space(opts,pattern);
    
    /*    sp_image_write(rs,"real_space.h5",sizeof(real));
	  sp_image_write(rs,"real_space.vtk",0);*/

    write_density_histogram(rs);
  }
#ifdef MPI
  MPI_Finalize();
#endif
  return 0;
}
