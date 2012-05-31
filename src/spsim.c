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
#include <sys/types.h>
#include <unistd.h>

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
  Molecule * mol = NULL;
  SpRotation ** rot = NULL;
  char buffer[1024];

#ifdef MPI
  MPI_Init(&argc, &argv);
#endif
  read_options_file("spsim.conf",opts);
  if(opts->random_seed < 0){
    opts->random_seed = getpid();
  }
  sp_srand(opts->random_seed);
  write_options_file("spsim.confout", opts);
  if(opts->detector->nz > 1){
    is3D = 1;
  }else{
    is3D = 0;
  }
  opts->detector->n_patterns = opts->n_patterns;

  if(opts->input_type == CHEM_FORMULA){
    mol = get_Molecule_from_formula(opts->chem_formula,opts);
    write_pdb_from_mol("molout.pdb",mol);
  }else if(opts->input_type == PDB){
    mol = get_Molecule_from_pdb(opts->pdb_filename);
  }else{
    fprintf(stderr,"input_type not specified!\n");
    exit(0);
  }

  rot = malloc(sizeof(SpRotation *)*opts->n_patterns);
  for(int n=0;n<opts->n_patterns;n++){
    if(opts->sf_filename[0]){
      pattern = load_pattern_from_file(opts->detector,opts->sf_filename,HKL_list,HKL_list_size);
    }else{
      if(is3D){
	HKL_list = get_HKL_list_for_3d_detector(opts->detector,opts->experiment, &HKL_list_size);
	rot[n] = apply_orientation_to_HKL_list(&HKL_list,&HKL_list_size,opts);
	if(opts->use_fft_for_sf){
	  pattern = compute_pattern_by_fft(mol,opts->detector,opts->experiment,opts->b_factor);
	}else if(opts->use_nfft_for_sf){
	  pattern = compute_pattern_by_nfft(mol,opts->detector,opts->experiment,opts->b_factor,HKL_list,opts);
	}else{
#ifdef _USE_CUDA
	  if(opts->use_cuda){
	    if(opts->use_cuda == 2){
	      pattern = cuda_compute_pattern_on_list(mol,HKL_list,HKL_list_size,opts->b_factor,opts->experiment,opts);	
	    }else{
	      pattern = cuda_compute_pattern_on_list2(mol,HKL_list,HKL_list_size,opts->b_factor,opts->experiment,opts);	
	    }
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
	rot[n] = apply_orientation_to_HKL_list(&HKL_list,&HKL_list_size,opts);
	if(opts->use_fft_for_sf){
	  fprintf(stderr,"Cannot use fft for 2D pattern calculation!\n");
	  abort();
	}else if(opts->use_nfft_for_sf){
	  pattern = compute_pattern_on_list_by_nfft(mol,HKL_list,HKL_list_size,opts->detector,opts->b_factor,opts);
	}else{
#ifdef _USE_CUDA
	  if(opts->use_cuda){
	    if(opts->use_cuda == 2){
	      pattern = cuda_compute_pattern_on_list(mol,HKL_list,HKL_list_size,opts->b_factor,opts->experiment,opts);	
	    }else{
	      pattern = cuda_compute_pattern_on_list2(mol,HKL_list,HKL_list_size,opts->b_factor,opts->experiment,opts);	
	    }
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
    
    if(opts->output_scatt_int){
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
      
      sprintf(buffer,"scattering_int-%05d.h5",n);
      sp_image_write(scattering_int,buffer,0);
      sp_image_free(scattering_int);
    }
    if(opts->fast_exit){
#ifdef MPI
      MPI_Finalize();
#endif
      return 0;
    }
    if(opts->output_sf_vtk){
      sprintf(buffer,"scattering_factor-%05d.vtk",n);
      write_3D_array_to_vtk(pattern->ints,opts->detector->nx,opts->detector->ny,opts->detector->nz,buffer);
    }
    if(opts->output_ewald_vtk){
      sprintf(buffer,"ewald-%05d.vtk",n);
      write_3D_array_as_structured_grid_to_vtk(pattern->ints,opts->detector,HKL_list,HKL_list_size,opts->n_patterns,buffer);  
    }
    if(opts->output_intensities){
      Image * ints = sp_image_alloc(opts->detector->nx,opts->detector->ny,opts->detector->nz);    
      i = 0;
      for(int x = 0;x<sp_image_x(ints);x++){
	for(int y = 0;y<sp_image_y(ints);y++){
	  for(int z = 0;z<sp_image_z(ints);z++){
	    sp_image_set(ints,x,y,z,sp_cinit(pattern->ints[sp_image_size(ints)*n+i++],0));
	    sp_i3matrix_set(ints->mask,x,y,z,1);
	  }
	}
      }
      ints->detector->wavelength = opts->experiment->wavelength;
      ints->detector->pixel_size[0] = opts->detector->pixel_width*opts->detector->binning_x;
      ints->detector->pixel_size[1] = opts->detector->pixel_height*opts->detector->binning_y;
      ints->detector->detector_distance = opts->detector->distance;
      if(rot){
	ints->detector->orientation = rot[n];
      }
      char buffer[1024];
      sprintf(buffer,"intensities-%05d.h5",n);
      sp_image_write(ints,buffer,sizeof(real));    
      sp_image_free(ints);
    }
    
    calculate_thomson_correction(opts->detector);
    calculate_pixel_solid_angle(opts->detector);
    calculate_photons_per_pixel(pattern,opts);
    
    
    if(opts->output_noiseless_photons){
      i = 0;
      Image * output = sp_image_alloc(opts->detector->nx,opts->detector->ny,opts->detector->nz);
      
      for(int x = 0;x<sp_image_x(output);x++){
	for(int y = 0;y<sp_image_y(output);y++){
	  for(int z = 0;z<sp_image_z(output);z++){
	    sp_image_set(output,x,y,z,sp_cinit(opts->detector->photons_per_pixel[i++],0));
	    sp_i3matrix_set(output->mask,x,y,z,1);
	  }
	}
      }
      output->detector->wavelength = opts->experiment->wavelength;
      output->detector->pixel_size[0] = opts->detector->pixel_width*opts->detector->binning_x;
      output->detector->pixel_size[1] = opts->detector->pixel_height*opts->detector->binning_y;
      output->detector->detector_distance = opts->detector->distance;
      if(rot){
	output->detector->orientation = rot[n];
      }
      sprintf(buffer,"noiseless_photon_out-%05d.cxi",n);
      sp_image_write(output,buffer,sizeof(real));
      sp_image_free(output);
    }
    
    generate_poisson_noise(opts->detector);
    
    if(opts->output_photons){
      i = 0;
      Image * output = sp_image_alloc(opts->detector->nx,opts->detector->ny,opts->detector->nz);
      
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
      if(rot){
	output->detector->orientation = rot[n];
      }
      sprintf(buffer,"photon_out-%05d.cxi",n);
      sp_image_write(output,buffer,sizeof(real));
      sp_image_free(output);
    }

    calculate_electrons_per_pixel(opts->detector,opts->experiment);
    if(opts->output_electron_pixel_vtk){
      sprintf(buffer,"electrons_per_pixel-%05d.vtk",n);
      write_3D_array_to_vtk(opts->detector->electrons_per_pixel,opts->detector->nx,opts->detector->ny,opts->detector->nz,buffer);
    }

    calculate_real_detector_output(opts->detector,opts->experiment);
    if(opts->output_count){      
      Image * output = sp_image_alloc(opts->detector->nx/opts->detector->binning_x,opts->detector->ny/opts->detector->binning_y,opts->detector->nz/opts->detector->binning_z);
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
      sprintf(buffer,"output-%05d.h5",n);
      sp_image_write(output,buffer,sizeof(real));
      sp_image_free(output);
    }
    
    
    calculate_noiseless_detector_output(opts->detector,opts->experiment);

    if(opts->output_noiseless_count){
      Image * noiseless = sp_image_alloc(opts->detector->nx/opts->detector->binning_x,opts->detector->ny/opts->detector->binning_y,opts->detector->nz/opts->detector->binning_z);
      noiseless->detector->wavelength = opts->experiment->wavelength;
      noiseless->detector->pixel_size[0] = opts->detector->pixel_width*opts->detector->binning_x;
      noiseless->detector->pixel_size[1] = opts->detector->pixel_height*opts->detector->binning_y;
      noiseless->detector->detector_distance = opts->detector->distance;
      if(rot){
	noiseless->detector->orientation = rot[n];
      }
      i = 0;
      for(int x = 0;x<sp_image_x(noiseless);x++){
	for(int y = 0;y<sp_image_y(noiseless);y++){
	  for(int z = 0;z<sp_image_z(noiseless);z++){
	    sp_image_set(noiseless,x,y,z,sp_cinit(opts->detector->noiseless_output[i++],0));
	    sp_i3matrix_set(noiseless->mask,x,y,z,1);
	  }
	}
      }
      sprintf(buffer,"noiseless_output-%05d.h5",n);
      sp_image_write(noiseless,buffer,sizeof(real));
      sp_image_free(noiseless);
    }
    if(opts->output_realspace_histogram && pattern->F){
      Image * rs = calculate_noiseless_real_space(opts,pattern);      
      write_density_histogram(rs);
      sp_image_free(rs);
    }
#ifdef MPI
    MPI_Finalize();
#endif
    free(HKL_list);
    free(pattern->F);
    free(pattern->ints);
    free(pattern->HKL_list);
    free(pattern);
    free(opts->detector->photons_per_pixel);
    free(opts->detector->thomson_correction);
    free(opts->detector->solid_angle);
    free(opts->detector->photon_count);
    free(opts->detector->electrons_per_pixel);
    free(opts->detector->real_output);
    free(opts->detector->noiseless_output);
  }
  free(mol->atomic_number);
  free(mol->pos);  
  free(opts->chem_formula);
  free(opts->experiment);
  free(opts->detector);
  free(opts);
  free(rot);

  return 0;
}
  
