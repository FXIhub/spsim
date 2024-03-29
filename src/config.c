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


#include <string.h>
#include <libconfig.h>
#include <stdlib.h>
#include <math.h>

#include "config.h"
#include "mpi_comm.h"


Options * set_defaults(){
  Options * opt = calloc(1,sizeof(Options));
  opt->chem_formula = calloc(1,sizeof(Chem_Formula));
  opt->experiment = calloc(1,sizeof(Experiment));
  opt->experiment->polarization = POLARIZATION_IGNORE;
  opt->detector = calloc(1,sizeof(CCD));
  opt->detector->binning_x = 1;
  opt->detector->binning_y = 1;
  opt->detector->binning_z = 1;
  opt->detector->center_x = 0;
  opt->detector->center_y = 0;
  opt->detector->center_z = 0;
  opt->box_type = BOX_SPHERICAL;
  opt->box_dimension = 1e-9;
  opt->use_fft_for_sf = 0;
  opt->use_nfft_for_sf = 0;
  opt->b_factor = 0;
  opt->n_patterns = 1;
  opt->origin_to_com = 0;
  opt->euler_orientation[0]= 0;
  opt->euler_orientation[1]= 0;
  opt->euler_orientation[2]= 0;
  opt->random_orientation = 0;
  opt->vectorize = 1;
#ifdef _USE_CUDA
  opt->use_cuda = 1;
#else
  opt->use_cuda = 0;
#endif
  opt->delta_atoms = 0;
  opt->fast_exit = 0;
  opt->crystal_size[0] = 1;
  opt->crystal_size[1] = 1;
  opt->crystal_size[2] = 1;
  opt->crystal_cell[0] = 0;
  opt->crystal_cell[1] = 0;
  opt->crystal_cell[2] = 0;
  opt->crystal_cell[3] = 90;
  opt->crystal_cell[4] = 90;
  opt->crystal_cell[5] = 90;
  opt->wavelength_samples = 5;
  opt->random_seed = -1;
  opt->output_sf_vtk = 0;
  opt->output_scatt_int = 0;
  opt->output_ewald_vtk = 0;
  opt->output_intensities = 0;
  opt->output_noiseless_photons = 1;
  opt->output_photons = 1;
  opt->output_electron_pixel_vtk = 0;
  opt->output_noiseless_count = 0;
  opt->output_count = 0;
  opt->output_realspace_histogram = 0;
  opt->output_scattering_factors = 0;
  opt->output_solid_angles = 0;
  opt->output_real_space = 0;
  opt->verbosity_level = 1;
  return opt;
}

void read_options_file(char * filename, Options * res){
  config_t config;
  const char * tmp;
  config_setting_t * root;
  config_setting_t * chem_formula;
  int i;
  char buffer[1024];

  config_init(&config);
  root = config_root_setting(&config);
  if(!config_read_file(&config,filename)){
    fprintf(stderr,"Error parsing %s on line %d:\n,%s\n",
	    filename,config_error_line(&config),config_error_text(&config));
    exit(1);
  }

  if((tmp = config_lookup_string(&config,"box_type"))){
    if(strcmp(tmp,"spherical") == 0){
      res->box_type = BOX_SPHERICAL;
    }else if(strcmp(tmp,"parallelepipedic") == 0){
      res->box_type = BOX_PARALLEL;
    }
  }

  if(config_lookup(&config,"box_dimension")){
    res->box_dimension = config_lookup_float(&config,"box_dimension");
  }

  if((tmp = config_lookup_string(&config,"input_type"))){
    if(strcmp(tmp,"chemical_formula") == 0){
      res->input_type = CHEM_FORMULA;
      res->chem_formula->n_elements = 0; 
      chem_formula = config_setting_get_member(root, "chemical_formula");
      for(i = 0;i<100;i++){
	sprintf(buffer,"chemical_formula/atom%d",i);
	if(config_lookup_int(&config,buffer)){
	  res->chem_formula->n_elements++;
	  res->chem_formula->atomic_number = realloc(res->chem_formula->atomic_number,sizeof(int)*res->chem_formula->n_elements);
	  res->chem_formula->atomic_number[res->chem_formula->n_elements-1] = i;
	  res->chem_formula->quantity = realloc(res->chem_formula->quantity,sizeof(int)*res->chem_formula->n_elements);
	  res->chem_formula->quantity[res->chem_formula->n_elements-1] = config_lookup_int(&config,buffer);
	}
      }      
    }else if(strcmp(tmp,"pdb") == 0){
      res->input_type = PDB;     
      if(config_lookup(&config,"pdb_filename")){
	strcpy(res->pdb_filename,config_lookup_string(&config,"pdb_filename"));
      }else{ 
	fprintf(stderr,"Error: You need to specify \"pdb_filename\" when you use input_type = \"pdb\"\n");
	exit(0);
      }
    }
  }
  if(config_lookup(&config,"detector_distance")){
    res->detector->distance = config_lookup_float(&config,"detector_distance");
  }
  if(config_lookup(&config,"detector_center_x")){
    res->detector->center_x = config_lookup_float(&config,"detector_center_x");
  }
  if(config_lookup(&config,"detector_center_y")){
    res->detector->center_y = config_lookup_float(&config,"detector_center_y");
  }
  if(config_lookup(&config,"detector_center_z")){
    res->detector->center_z = config_lookup_float(&config,"detector_center_z");
  }
  if(config_lookup(&config,"detector_width")){
    res->detector->width = config_lookup_float(&config,"detector_width");
  }
  if(config_lookup(&config,"detector_height")){
    res->detector->height = config_lookup_float(&config,"detector_height");
  }
  if(config_lookup(&config,"detector_depth")){
    res->detector->depth = config_lookup_float(&config,"detector_depth");
  }
  if(config_lookup(&config,"detector_pixel_width")){
    res->detector->pixel_width = config_lookup_float(&config,"detector_pixel_width");
  }
  if(config_lookup(&config,"detector_pixel_height")){
    res->detector->pixel_height = config_lookup_float(&config,"detector_pixel_height");
  }
  if(config_lookup(&config,"detector_pixel_depth")){
    res->detector->pixel_depth = config_lookup_float(&config,"detector_pixel_depth");
  }
  if(config_lookup(&config,"detector_quantum_efficiency")){
    res->detector->quantum_efficiency = config_lookup_float(&config,"detector_quantum_efficiency");
  }
  if(config_lookup(&config,"detector_electron_hole_production_energy")){
    res->detector->electron_hole_production_energy = config_lookup_float(&config,"detector_electron_hole_production_energy");
  }
  if(config_lookup(&config,"detector_readout_noise")){
    res->detector->readout_noise = config_lookup_float(&config,"detector_readout_noise");
  }
  if(config_lookup(&config,"detector_dark_current")){
    res->detector->dark_current = config_lookup_float(&config,"detector_dark_current");
  }
  if(config_lookup(&config,"detector_linear_full_well")){
    res->detector->linear_full_well = config_lookup_float(&config,"detector_linear_full_well");
  }
  if(config_lookup(&config,"detector_binning")){
    res->detector->binning_x = config_lookup_int(&config,"detector_binning");
    res->detector->binning_y = config_lookup_int(&config,"detector_binning");
    res->detector->binning_z = 1;
  }
  if(config_lookup(&config,"detector_binning_x")){
    res->detector->binning_x = config_lookup_int(&config,"detector_binning_x");
  }
  if(config_lookup(&config,"detector_binning_y")){
    res->detector->binning_y = config_lookup_int(&config,"detector_binning_y");
  }
  if(config_lookup(&config,"detector_binning_z")){
    res->detector->binning_z = config_lookup_int(&config,"detector_binning_z");
  }
  if(config_lookup(&config,"detector_maximum_value")){
    res->detector->maximum_value = config_lookup_float(&config,"detector_maximum_value");
  }

  if(config_lookup(&config,"experiment_wavelength")){
    res->experiment->wavelength = config_lookup_float(&config,"experiment_wavelength");
  }
  if(config_lookup(&config,"experiment_photon_energy")){
    res->experiment->photon_energy = config_lookup_float(&config,"experiment_photon_energy");
  }
  if(config_lookup(&config,"experiment_bandwidth")){
    res->experiment->bandwidth = config_lookup_float(&config,"experiment_bandwidth");
  }
  if(config_lookup(&config,"experiment_exposure_time")){
    res->experiment->exposure_time = config_lookup_float(&config,"experiment_exposure_time");
  }
  if(config_lookup(&config,"experiment_beam_intensity")){
    res->experiment->beam_intensity = config_lookup_float(&config,"experiment_beam_intensity");
  }
  if(config_lookup(&config,"experiment_beam_center_x")){
    res->experiment->beam_center_x = config_lookup_float(&config,"experiment_beam_center_x");
  }
  if(config_lookup(&config,"experiment_beam_center_y")){
    res->experiment->beam_center_y = config_lookup_float(&config,"experiment_beam_center_y");
  }
  if(config_lookup(&config,"experiment_beam_fwhm")){
    res->experiment->beam_fwhm = config_lookup_float(&config,"experiment_beam_fwhm");
  }
  if(config_lookup(&config,"experiment_beam_energy")){
    res->experiment->beam_energy = config_lookup_float(&config,"experiment_beam_energy");
  }
  if(config_lookup(&config,"experiment_focal_diameter")){
    res->experiment->focal_diameter = config_lookup_float(&config,"experiment_focal_diameter");
  }

  if((tmp = config_lookup_string(&config,"experiment_polarization"))){
    if(strcmp(tmp,"ignore") == 0){
      res->experiment->polarization = POLARIZATION_IGNORE;
    }else if(strcmp(tmp,"horizontal") == 0){
      res->experiment->polarization = POLARIZATION_HORIZONTAL;
    }else if(strcmp(tmp,"vertical") == 0){
      res->experiment->polarization = POLARIZATION_VERTICAL;
    }else if(strcmp(tmp,"unpolarized") == 0){
      res->experiment->polarization = POLARIZATION_UNPOLARIZED;
    }else{
      fprintf(stderr,"Warning: polarization is set to an undefined option!\n");
    }
  }
  
  if(config_lookup_string(&config,"precalculated_sf")){
    strcpy(res->sf_filename,config_lookup_string(&config,"precalculated_sf"));
  }

  if(config_lookup_string(&config,"hkl_grid")){
    strcpy(res->hkl_grid_filename,config_lookup_string(&config,"hkl_grid"));
  }

  if(config_lookup(&config,"detector_spherical")){
    res->detector->spherical = config_lookup_int(&config,"detector_spherical");
  }

  if(config_lookup(&config,"detector_gaussian_blurring")){
    res->detector->gaussian_blurring = config_lookup_float(&config,"detector_gaussian_blurring");
  }

  if(config_lookup(&config,"real_space_blurring")){
    res->detector->real_space_blurring = config_lookup_float(&config,"real_space_blurring");
  }

  if(config_lookup(&config,"use_fft_for_sf")){
    res->use_fft_for_sf = config_lookup_int(&config,"use_fft_for_sf");
  }
  if(config_lookup(&config,"use_nfft_for_sf")){
    res->use_nfft_for_sf = config_lookup_int(&config,"use_nfft_for_sf");
  }
  if(config_lookup(&config,"b_factor")){
    res->b_factor = config_lookup_float(&config,"b_factor");
  }
  if(config_lookup(&config,"origin_to_com")){
    res->origin_to_com = config_lookup_int(&config,"origin_to_com");
  }
  if(config_lookup(&config,"phi")){
    res->euler_orientation[0] = config_lookup_float(&config,"phi");
  }
  if(config_lookup(&config,"theta")){
    res->euler_orientation[1] = config_lookup_float(&config,"theta");
  }
  if(config_lookup(&config,"psi")){
    res->euler_orientation[2] = config_lookup_float(&config,"psi");
  }
  if(config_lookup(&config,"random_orientation")){
    res->random_orientation = config_lookup_int(&config,"random_orientation");
  }
  if(config_lookup(&config,"number_of_patterns")){
    res->n_patterns = config_lookup_int(&config,"number_of_patterns");
  }
  if(config_lookup(&config,"vectorize")){
    res->vectorize = config_lookup_int(&config,"vectorize");
  }
  if(config_lookup(&config,"delta_atoms")){
    res->delta_atoms = config_lookup_int(&config,"delta_atoms");
  }
  if(config_lookup(&config,"fast_exit")){
    res->fast_exit = config_lookup_int(&config,"fast_exit");
  }
  if(config_lookup(&config,"crystal_size_a")){
    res->crystal_size[0] = config_lookup_int(&config,"crystal_size_a");
  }
  if(config_lookup(&config,"crystal_size_b")){
    res->crystal_size[1] = config_lookup_int(&config,"crystal_size_b");
  }
  if(config_lookup(&config,"crystal_size_c")){
    res->crystal_size[2] = config_lookup_int(&config,"crystal_size_c");
  }
  if(config_lookup(&config,"crystal_cell_a")){
    res->crystal_cell[0] = config_lookup_float(&config,"crystal_cell_a");
  }
  if(config_lookup(&config,"crystal_cell_b")){
    res->crystal_cell[1] = config_lookup_float(&config,"crystal_cell_b");
  }
  if(config_lookup(&config,"crystal_cell_c")){
    res->crystal_cell[2] = config_lookup_float(&config,"crystal_cell_c");
  }
  if(config_lookup(&config,"crystal_cell_alpha")){
    res->crystal_cell[3] = config_lookup_float(&config,"crystal_cell_alpha");
  }
  if(config_lookup(&config,"crystal_cell_beta")){
    res->crystal_cell[4] = config_lookup_float(&config,"crystal_cell_beta");
  }
  if(config_lookup(&config,"crystal_cell_gamma")){
    res->crystal_cell[5] = config_lookup_float(&config,"crystal_cell_gamma");
  }

  if(config_lookup(&config,"use_cuda")){
    res->use_cuda = config_lookup_int(&config,"use_cuda");
  }
#ifndef _USE_CUDA
  if(res->use_cuda){
    res->use_cuda = 0;
    fprintf(stderr,"Warning: use_cuda set to true but spsim was not compiled with CUDA support!\n");
    fprintf(stderr,"Warning: Setting use_cuda to false.\n");
  }
#else
  if(res->use_cuda){
    int deviceCount;
    cudaError_t cudaResultCode = cudaGetDeviceCount(&deviceCount);
    if(cudaResultCode != cudaSuccess){
      res->use_cuda = 0;
      fprintf(stderr,"Warning: use_cuda set to true but got \"%s\" while initializing CUDA!\n",cudaGetErrorString(cudaResultCode));
      fprintf(stderr,"Warning: Setting use_cuda to false.\n");
    }else if(deviceCount == 0){
      res->use_cuda = 0;
      fprintf(stderr,"Warning: use_cuda set to true but did not find any CUDA devices!\n");
      fprintf(stderr,"Warning: Setting use_cuda to false.\n");
    }
  }
#endif
  if(config_lookup(&config,"wavelength_samples")){
    res->wavelength_samples = config_lookup_int(&config,"wavelength_samples");
  }
  if(config_lookup(&config,"random_seed")){
    res->random_seed = config_lookup_int(&config,"random_seed");
  }

  if(config_lookup(&config,"output_noiseless_photons")){
    res->output_noiseless_photons = config_lookup_int(&config,"output_noiseless_photons");
  }
  if(config_lookup(&config,"output_photons")){
    res->output_photons = config_lookup_int(&config,"output_photons");
  }
  if(config_lookup(&config,"output_noiseless_count")){
    res->output_noiseless_count = config_lookup_int(&config,"output_noiseless_count");
  }
  if(config_lookup(&config,"output_count")){
    res->output_count = config_lookup_int(&config,"output_count");
  }
  if(config_lookup(&config,"output_scattering_factors")){
    res->output_scattering_factors = config_lookup_int(&config,"output_scattering_factors");
  }
  if(config_lookup(&config,"output_solid_angles")){
    res->output_solid_angles = config_lookup_int(&config,"output_solid_angles");
  }
  if(config_lookup(&config,"output_real_space")){
    res->output_real_space = config_lookup_int(&config,"output_real_space");
  }
  res->detector->nx = rint(res->detector->width/res->detector->pixel_width);
  res->detector->ny = rint(res->detector->height/res->detector->pixel_height);
  if(res->detector->pixel_depth){
    res->detector->nz = rint(res->detector->depth/res->detector->pixel_depth);
  }else{
    res->detector->nz = 1;
    res->detector->binning_z = 1;
  }
  if((res->detector->nz > 1) && (res->experiment->polarization != POLARIZATION_IGNORE)){
    fprintf(stderr,"Warning: For a 3D detector polarization effects (Thomson correction) have to be neglected (polarization = \"ignore\")!\n");
    res->experiment->polarization = POLARIZATION_IGNORE;
  }

  if(res->experiment->photon_energy){
    float lambda = 1.240e-6/res->experiment->photon_energy;
    if(res->experiment->wavelength){
      if((lambda-res->experiment->wavelength)/(lambda+res->experiment->wavelength) > 0.01){
	fprintf(stderr,"Warning: experiment_wavelength does not agree with experiment_photon_energy!\n");
      }
    }else{
      res->experiment->wavelength = lambda;
    }
  }else{
    if(res->experiment->wavelength){
      float eV = 1.240e-6/res->experiment->wavelength;
      res->experiment->photon_energy = eV;
    }else{
      fprintf(stderr,"Warning: you need to specify the experiment_wavelength!\n");
    }
  }
  if(!res->experiment->beam_intensity){
    if(res->experiment->focal_diameter && res->experiment->beam_energy){
      float eV = 1.240e-6/res->experiment->wavelength;
      float nphotons = res->experiment->beam_energy*6.24150974e18/eV;
      float area = res->experiment->focal_diameter*res->experiment->focal_diameter;
      res->experiment->beam_intensity = nphotons/area;
    }
  }
  if(config_lookup(&config,"verbosity_level")){
    res->verbosity_level = config_lookup_int(&config,"verbosity_level");
  }
}





void write_options_file(char * filename, Options * res){
  config_t config;
  config_setting_t * root;
  config_setting_t * s;
  config_setting_t * chem_formula;
  int i;
  char buffer[1024];
#ifdef MPI
  if(!is_mpi_master()){
    return;
  }
#endif  
  config_init(&config);
  root = config_root_setting(&config);
  if(res->input_type == CHEM_FORMULA){
    s = config_setting_add(root,"input_type",CONFIG_TYPE_STRING);
    config_setting_set_string(s,"chemical_formula");
    chem_formula = config_setting_add(root,"chemical_formula",CONFIG_TYPE_GROUP);
    for(i = 0;i<res->chem_formula->n_elements;i++){
      sprintf(buffer,"atom%d",res->chem_formula->atomic_number[i]);
      s = config_setting_add(chem_formula,buffer,CONFIG_TYPE_INT);
      config_setting_set_int(s,res->chem_formula->quantity[i]);      
    }    
  }else if(res->input_type == PDB){
    s = config_setting_add(root,"input_type",CONFIG_TYPE_STRING);
    config_setting_set_string(s,"pdb");
    s = config_setting_add(root,"pdb_filename",CONFIG_TYPE_STRING);
    config_setting_set_string(s,res->pdb_filename);
  }else if(!res->input_type){
    s = config_setting_add(root,"input_type",CONFIG_TYPE_STRING);
    config_setting_set_string(s,"none");
  }

  if(res->box_type == BOX_SPHERICAL){
    s = config_setting_add(root,"box_type",CONFIG_TYPE_STRING);
    config_setting_set_string(s,"spherical");
  }else if(res->box_type == BOX_PARALLEL){
    s = config_setting_add(root,"box_type",CONFIG_TYPE_STRING);
    config_setting_set_string(s,"parallelepipedic");
  }
  
  s = config_setting_add(root,"detector_distance",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->distance);
  s = config_setting_add(root,"detector_width",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->width);
  s = config_setting_add(root,"detector_height",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->height);
  s = config_setting_add(root,"detector_depth",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->depth);
  s = config_setting_add(root,"detector_center_x",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->center_x);
  s = config_setting_add(root,"detector_center_y",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->center_y);
  s = config_setting_add(root,"detector_center_z",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->center_z);
  s = config_setting_add(root,"detector_pixel_width",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->pixel_width);
  s = config_setting_add(root,"detector_pixel_height",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->pixel_height);
  s = config_setting_add(root,"detector_pixel_depth",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->pixel_depth);
  s = config_setting_add(root,"detector_quantum_efficiency",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->quantum_efficiency);
  s = config_setting_add(root,"detector_electron_hole_production_energy",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->electron_hole_production_energy);
  s = config_setting_add(root,"detector_readout_noise",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->readout_noise);
  s = config_setting_add(root,"detector_dark_current",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->dark_current);
  s = config_setting_add(root,"detector_linear_full_well",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->linear_full_well);
  s = config_setting_add(root,"detector_binning_x",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->detector->binning_x);
  s = config_setting_add(root,"detector_binning_y",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->detector->binning_y);
  s = config_setting_add(root,"detector_binning_z",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->detector->binning_z);
  s = config_setting_add(root,"detector_maximum_value",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->maximum_value);



  s = config_setting_add(root,"experiment_wavelength",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->experiment->wavelength);
  s = config_setting_add(root,"experiment_bandwidth",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->experiment->bandwidth);
  s = config_setting_add(root,"experiment_exposure_time",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->experiment->exposure_time);
  s = config_setting_add(root,"experiment_beam_intensity",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->experiment->beam_intensity);

  s = config_setting_add(root,"experiment_beam_center_x",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->experiment->beam_center_x);
  s = config_setting_add(root,"experiment_beam_center_y",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->experiment->beam_center_y);
  s = config_setting_add(root,"experiment_beam_fwhm",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->experiment->beam_fwhm);
  s = config_setting_add(root,"experiment_beam_energy",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->experiment->beam_energy);
  s = config_setting_add(root,"experiment_photon_energy",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->experiment->photon_energy);
  s = config_setting_add(root,"experiment_focal_diameter",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->experiment->focal_diameter);


  s = config_setting_add(root,"box_dimension",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->box_dimension);
  
  s = config_setting_add(root,"detector_spherical",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->detector->spherical);

  s = config_setting_add(root,"detector_gaussian_blurring",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->gaussian_blurring);

  s = config_setting_add(root,"real_sapce_blurring",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->real_space_blurring);

  s = config_setting_add(root,"use_fft_for_sf",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->use_fft_for_sf);

  s = config_setting_add(root,"use_nfft_for_sf",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->use_nfft_for_sf);

  s = config_setting_add(root,"b_factor",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->b_factor);

  s = config_setting_add(root,"origin_to_com",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->origin_to_com);
  
  s = config_setting_add(root,"phi",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->euler_orientation[0]);
  s = config_setting_add(root,"theta",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->euler_orientation[1]);
  s = config_setting_add(root,"psi",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->euler_orientation[2]);

  s = config_setting_add(root,"random_orientation",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->random_orientation);
  s = config_setting_add(root,"number_of_patterns",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->n_patterns);
  s = config_setting_add(root,"vectorize",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->vectorize);
  s = config_setting_add(root,"delta_atoms",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->delta_atoms);
  s = config_setting_add(root,"fast_exit",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->fast_exit);

  s = config_setting_add(root,"crystal_size_a",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->crystal_size[0]);
  s = config_setting_add(root,"crystal_size_b",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->crystal_size[1]);
  s = config_setting_add(root,"crystal_size_c",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->crystal_size[2]);

  s = config_setting_add(root,"crystal_cell_a",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->crystal_cell[0]);
  s = config_setting_add(root,"crystal_cell_b",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->crystal_cell[1]);
  s = config_setting_add(root,"crystal_cell_c",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->crystal_cell[2]);
  s = config_setting_add(root,"crystal_cell_alpha",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->crystal_cell[3]);
  s = config_setting_add(root,"crystal_cell_beta",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->crystal_cell[4]);
  s = config_setting_add(root,"crystal_cell_gamma",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->crystal_cell[5]);

  if(res->experiment->polarization == POLARIZATION_IGNORE){
    s = config_setting_add(root,"experiment_polarization",CONFIG_TYPE_STRING);
    config_setting_set_string(s,"ignore");
  }else if(res->experiment->polarization == POLARIZATION_VERTICAL){
    s = config_setting_add(root,"experiment_polarization",CONFIG_TYPE_STRING);
    config_setting_set_string(s,"vertical");
  }else if(res->experiment->polarization == POLARIZATION_HORIZONTAL){
    s = config_setting_add(root,"experiment_polarization",CONFIG_TYPE_STRING);
    config_setting_set_string(s,"horizontal");
  }else if(res->experiment->polarization == POLARIZATION_UNPOLARIZED){
    s = config_setting_add(root,"experiment_polarization",CONFIG_TYPE_STRING);
    config_setting_set_string(s,"unpolarized");
  }
  
  s = config_setting_add(root,"use_cuda",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->use_cuda);
  s = config_setting_add(root,"wavelength_samples",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->wavelength_samples);
  s = config_setting_add(root,"random_seed",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->random_seed);
  s = config_setting_add(root,"output_photons",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->output_photons);
  s = config_setting_add(root,"output_noiseless_photons",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->output_noiseless_photons);
  s = config_setting_add(root,"output_count",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->output_count);
  s = config_setting_add(root,"output_noiseless_count",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->output_noiseless_count);
  s = config_setting_add(root,"output_solid_angles",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->output_solid_angles);
  s = config_setting_add(root,"output_scattering_factors",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->output_scattering_factors);
  s = config_setting_add(root,"output_real_space",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->output_real_space);
  s = config_setting_add(root,"verbosity_level",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->verbosity_level);

  if(res->sf_filename[0]){
    s = config_setting_add(root,"precalculated_sf",CONFIG_TYPE_STRING);
    config_setting_set_string(s,res->sf_filename);
  }
  config_write_file(&config,filename);
}
