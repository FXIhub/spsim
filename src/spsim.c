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
#include <sys/types.h>
#include <unistd.h>

#include "spsim.h"
#include "config.h"

void gaussian_blur_real_space(Options * opts,Diffraction_Pattern * pattern){
  float radius = opts->detector->real_space_blurring;
  if(!radius){
    return;
  }
  Image * real_space = calculate_noiseless_real_space(opts,pattern);
  
  int i = 0;
  for(int z = 0;z<opts->detector->nz;z++){
    float dz = abs(z-((opts->detector->nz-1)/2));
    for(int y = 0;y<opts->detector->ny;y++){
      float dy = abs(y-((opts->detector->ny-1)/2));
      for(int x = 0;x<opts->detector->nx;x++){
	float dx = abs(x-((opts->detector->nx-1)/2));
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
  float d_c = sqrt(opts->detector->nx*opts->detector->nx/4.+opts->detector->ny*opts->detector->ny/4.+opts->detector->nz*opts->detector->nz/4.);
  radius *= d_c;
  /*  f(x,y) = 1/sqrt(2*M_PI*radius) * exp(-(x^2+y^2+dz^2)/(2*radius^2)) */
  
  int i = 0;
  for(int z = 0;z<opts->detector->nz;z++){
    float dz = abs(z-((opts->detector->nz-1)/2));
    for(int y = 0;y<opts->detector->ny;y++){
      float dy = abs(y-((opts->detector->ny-1)/2));
      for(int x = 0;x<opts->detector->nx;x++){
	float dx = abs(x-((opts->detector->nx-1)/2));
        float factor = 1/sqrt(2*M_PI*radius) * exp(-(dx*dx+dy*dy+dz*dz)/(2*radius*radius));
        sp_real(pattern->F[i]) *= factor;
        sp_imag(pattern->F[i]) *= factor;
        pattern->ints[i] = sp_cabs(pattern->F[i])*sp_cabs(pattern->F[i]);
        i++;
      }
    }
  }
}


Molecule * get_molecule(Options * opts){
  Molecule * mol = NULL;
  if(opts->input_type == CHEM_FORMULA){
    mol = get_Molecule_from_formula(opts->chem_formula,opts);
    write_pdb_from_mol("molout.pdb",mol);    
  }else if(opts->input_type == PDB){
    mol = get_Molecule_from_pdb(opts->pdb_filename);
  }else{
    fprintf(stderr,"input_type not specified!\n");
    exit(0);
  }
  return mol;
}

float * get_HKL_list(Options * opts, int * HKL_list_size){
  float * HKL_list = NULL;
  int is3D = (opts->detector->nz > 1);
  if(is3D){
    HKL_list = get_HKL_list_for_3d_detector(opts->detector,opts->experiment, HKL_list_size);
  }else{
    HKL_list = get_HKL_list_for_detector(opts->detector,opts->experiment, HKL_list_size);
  }

  return HKL_list;
}

Diffraction_Pattern * compute_sf(Molecule * mol, float * HKL_list, int HKL_list_size, Options * opts){
  Diffraction_Pattern * pattern = NULL;
  int is3D = (opts->detector->nz > 1);
  if(opts->use_fft_for_sf){
    if(is3D == 0){
      fprintf(stderr,"Cannot use fft for 2D pattern calculation!\n");
      abort();
    }
    pattern = compute_pattern_by_fft(mol,opts->detector,opts->experiment,opts->b_factor);
  }else if(opts->use_nfft_for_sf){
#ifdef NFFT_SUPPORT    
    pattern = compute_pattern_by_nfft(mol,opts->detector,opts->experiment,opts->b_factor,HKL_list,opts);
#else
    fprintf(stderr,"spsim built without NFFT support!\n");
#endif
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
  pattern->rot = NULL;
  return pattern;
}

void output_array(char * basename, int index, float * array, SpRotation * rot, Options * opts){
  Image * img = make_image(array, rot, opts);
  char buffer[1024];
  sprintf(buffer,"%s-%05d.cxi", basename, index);
  sp_image_write(img,buffer,0);    
  sp_image_free(img);
}

Image * make_image(float * array, SpRotation * rot, Options * opts){
  Image * img = sp_image_alloc(opts->detector->nx,opts->detector->ny,opts->detector->nz);    
  int i = 0;
  for(int z = 0;z<sp_image_z(img);z++){
    for(int y = 0;y<sp_image_y(img);y++){
      for(int x = 0;x<sp_image_x(img);x++){
	sp_image_set(img,x,y,z,sp_cinit(array[i++],0));
	sp_i3matrix_set(img->mask,x,y,z,1);
      }
    }
  }
  img->detector->wavelength = opts->experiment->wavelength;
  img->detector->pixel_size[0] = opts->detector->pixel_width*opts->detector->binning_x;
  img->detector->pixel_size[1] = opts->detector->pixel_height*opts->detector->binning_y;
  img->detector->pixel_size[2] = opts->detector->pixel_depth*opts->detector->binning_z;
  img->detector->detector_distance = opts->detector->distance;
  if(rot){
    img->detector->orientation = rot;
  }
  return img;
}

void output_carray(char * basename, int index, Complex * array, SpRotation * rot, Options * opts){
  Image * img = make_cimage(array, rot, opts);
  char buffer[1024];
  sprintf(buffer,"%s-%05d.cxi", basename, index);
  sp_image_write(img,buffer,0);    
  sp_image_free(img);
}

void array_to_image(float * array, Image * img) {
  int i = 0;
  for(int z = 0;z<sp_image_z(img);z++){
    for(int y = 0;y<sp_image_y(img);y++){
      for(int x = 0;x<sp_image_x(img);x++){
	sp_image_set(img,x,y,z,sp_cinit(array[i++],0));
	sp_i3matrix_set(img->mask,x,y,z,1);
      }
    }
  }
}

void iarray_to_image(int * array, Image * img) {
  int i = 0;
  for(int z = 0;z<sp_image_z(img);z++){
    for(int y = 0;y<sp_image_y(img);y++){
      for(int x = 0;x<sp_image_x(img);x++){
	sp_image_set(img,x,y,z,sp_cinit((float) array[i++],0));
	sp_i3matrix_set(img->mask,x,y,z,1);
      }
    }
  }
}

Image * make_cimage(Complex * array, SpRotation * rot, Options * opts){
  Image * img = sp_image_alloc(opts->detector->nx,opts->detector->ny,opts->detector->nz);    
  int i = 0;
  for(int z = 0;z<sp_image_z(img);z++){
    for(int y = 0;y<sp_image_y(img);y++){
      for(int x = 0;x<sp_image_x(img);x++){
	sp_image_set(img,x,y,z,array[i++]);
	sp_i3matrix_set(img->mask,x,y,z,1);
      }
    }
  }
  img->detector->wavelength = opts->experiment->wavelength;
  img->detector->pixel_size[0] = opts->detector->pixel_width*opts->detector->binning_x;
  img->detector->pixel_size[1] = opts->detector->pixel_height*opts->detector->binning_y;
  img->detector->pixel_size[2] = opts->detector->pixel_depth*opts->detector->binning_z;
  img->detector->detector_distance = opts->detector->distance;
  if(rot){
    img->detector->orientation = rot;
  }
  return img;
}

void output_files( Diffraction_Pattern * pattern, Options * opts, int index){
  if(opts->output_scatt_int || opts->output_intensities){
    output_array("scattering_int", index, pattern->ints, pattern->rot, opts);
  }

  if(opts->output_sf_vtk){
    char buffer[1024];
    sprintf(buffer,"scattering_factor-%05d.vtk",index);
    write_3D_array_to_vtk(pattern->ints,opts->detector->nx,opts->detector->ny,opts->detector->nz,buffer);
  }
    
  if(opts->output_noiseless_photons){
    output_array("noiseless_photon_out", index, opts->detector->photons_per_pixel, pattern->rot, opts);
  }

  if(opts->output_scattering_factors){
    output_carray("scattering_factors", index, pattern->F, pattern->rot, opts);
  }

  if(opts->output_real_space){ 
    Image * real_space = calculate_noiseless_real_space(opts,pattern);
    char buffer[1024];     
    sprintf(buffer,"real_space-%05d.cxi",index);
    sp_image_write(real_space,buffer,0);
    sp_image_free(real_space);
  }
    
  if(opts->output_photons){
    output_array("photon_out", index, opts->detector->photon_count, pattern->rot, opts);
  }

  if(opts->output_count){      
    output_array("counts_out", index, opts->detector->real_output, pattern->rot, opts);    
  }


  if(opts->output_noiseless_count){
    output_array("noiseless_counts_out", index, opts->detector->noiseless_output, pattern->rot, opts);    
  }

  if(opts->output_realspace_histogram && pattern->F){
    Image * rs = calculate_noiseless_real_space(opts,pattern);      
    write_density_histogram(rs);
    sp_image_free(rs);
  }
}

Diffraction_Pattern * simulate_shot(Molecule * mol, Options * opts){  
  SpRotation * rot = NULL;
  Diffraction_Pattern * pattern = NULL;
  float * HKL_list;

  //printf("mol->pos[13700]=%e",mol->pos[13700]);
  //float foo_x0,foo_x1;
  //float foo_y0,foo_y1;
  //float foo_z0,foo_z1;
  //long i = 0;

  /*
  foo_x0 = 10000000000.;
  foo_x1 = -10000000000.;
  foo_y0 = 10000000000.;
  foo_y1 = -10000000000.;
  foo_z0 = 10000000000.;
  foo_z1 = -10000000000.;
  for (i=0;i<mol->natoms*3;i+=3){
    if (mol->pos[i+2] < foo_x0) {
      foo_x0 = mol->pos[i+2];
    }
    if (mol->pos[i+1] < foo_y0) {
      foo_y0 = mol->pos[i+1];
    }
    if (mol->pos[i+0] < foo_z0) {
      foo_z0 = mol->pos[i+0];
    }
    if (mol->pos[i+2] > foo_x1) {
      foo_x1 = mol->pos[i+2];
    }
    if (mol->pos[i+1] > foo_y1) {
      foo_y1 = mol->pos[i+1];
    }
    if (mol->pos[i+0] > foo_z1) {
      foo_z1 = mol->pos[i+0];
    }
  }
  */
  //printf("Lx=%e, Ly=%e, Lz=%e\n",foo_x1-foo_x0,foo_y1-foo_y0,foo_z1-foo_z0);

  int HKL_list_size = 0;  
  if(opts->sf_filename[0]){
    pattern = load_pattern_from_file(opts->detector,opts->sf_filename,HKL_list,HKL_list_size);
  }else{
    HKL_list = get_HKL_list(opts, &HKL_list_size);
    rot = apply_orientation_to_HKL_list(&HKL_list,&HKL_list_size,opts);
    pattern = compute_sf(mol, HKL_list, HKL_list_size, opts);
    pattern->rot = rot;
  }
  
  // Blur the pattern with a gaussian, if you feel like it
  gaussian_blur_pattern(opts,pattern);    
  calculate_thomson_correction(opts->detector);
  calculate_pixel_solid_angle(opts->detector);
  calculate_photons_per_pixel(pattern,opts);
    
  generate_poisson_noise(opts->detector);
    
  calculate_electrons_per_pixel(opts->detector,opts->experiment);
  calculate_real_detector_output(opts->detector,opts->experiment);
  calculate_noiseless_detector_output(opts->detector,opts->experiment);
  
  free(HKL_list);
  return pattern;
}

void free_stuff(Diffraction_Pattern * pattern, Options * opts){
  free_diffraction_pattern(pattern);
  free_output_in_options(opts);
}

void free_diffraction_pattern(Diffraction_Pattern * pattern) {
  free(pattern->F);
  free(pattern->ints);
  free(pattern->HKL_list);
  if(pattern->rot){
    free(pattern->rot);
  }
  free(pattern);
}

void free_output_in_options(Options * opts){
  free(opts->detector->photons_per_pixel);
  free(opts->detector->thomson_correction);
  free(opts->detector->solid_angle);
  free(opts->detector->photon_count);
  free(opts->detector->electrons_per_pixel);
  free(opts->detector->real_output);
  free(opts->detector->noiseless_output);  
}

#ifndef _LIBRARY
int main(int argc, char ** argv){
  Options * opts = set_defaults();
  Molecule * mol = NULL;
  
#ifdef MPI
  MPI_Init(&argc, &argv);
#endif
  read_options_file("spsim.conf",opts);
  if(opts->random_seed < 0){
    opts->random_seed = getpid();
  }
  sp_srand(opts->random_seed);
  write_options_file("spsim.confout", opts);

  opts->detector->n_patterns = opts->n_patterns;
  
  mol = get_molecule(opts);

  for(int n=0;n<opts->n_patterns;n++){
    Diffraction_Pattern * pattern = simulate_shot(mol, opts);
    output_files(pattern, opts, n);
    free_stuff(pattern, opts);
  }
#ifdef MPI
  if(!is_mpi_master()){
    MPI_Finalize();
    return 0;
  }
#endif
  free(mol->atomic_number);
  free(mol->pos);  
  free(opts->chem_formula);
  free(opts->experiment);
  free(opts->detector);
  free(opts);

  return 0;
}
#endif /* _LIBRARY */
  
