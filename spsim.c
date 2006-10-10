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

int main(int argc, char ** argv){
  Options * opts = set_defaults();
  Diffraction_Pattern * pattern;
  /* 20x100x100 nm box */
  float box[3] = {20e-9,100e-9,100e-9};
  float * HKL_list;
  int HKL_list_size = 0;
  int i;
  Image * noiseless;
  Molecule * mol = NULL;
#ifdef MPI
  MPI_Init(&argc, &argv);
#endif
  read_options_file("spsim.conf",opts);
  write_options_file("spsim.confout", opts);
  if(opts->input_type == CHEM_FORMULA){
    mol = get_Molecule_from_formula(opts->chem_formula,box);
  }else if(opts->input_type == PDB){
    mol = get_Molecule_from_pdb(opts->pdb_filename);
  }else{
    fprintf(stderr,"input_type not specified!\n");
    exit(0);
  }
  HKL_list = get_HKL_list_for_detector(opts->detector,opts->experiment, &HKL_list_size);

  write_hkl_grid(HKL_list,mol,opts->detector);

  if(opts->sf_filename[0]){
    pattern = load_pattern_from_file(opts->detector,opts->sf_filename,HKL_list,HKL_list_size);
  }else{
    pattern = compute_pattern_on_list(mol,HKL_list,HKL_list_size);
  }
#ifdef MPI
  if(!is_mpi_master()){
    return 0;
  }
#endif
  write_2D_array_to_vtk(pattern->I,opts->detector->nx,opts->detector->ny,"scattering_factor.vtk");
  calculate_thomson_correction(opts->detector);
  calculate_pixel_solid_angle(opts->detector);
  write_2D_array_to_vtk(opts->detector->thomson_correction,opts->detector->nx,opts->detector->ny,"thomson_correction.vtk");
  write_2D_array_to_vtk(opts->detector->solid_angle,opts->detector->nx,opts->detector->ny,"solid_angle.vtk");
  calculate_photons_per_pixel(pattern,opts->detector,opts->experiment);
  write_2D_array_to_vtk(opts->detector->photons_per_pixel,opts->detector->nx,opts->detector->ny,"pattern.vtk");
  generate_poisson_noise(opts->detector);
  write_2D_array_to_vtk(opts->detector->photon_count,opts->detector->nx,opts->detector->ny,"photon_count.vtk");
  calculate_electrons_per_pixel(opts->detector,opts->experiment);
  write_2D_array_to_vtk(opts->detector->electrons_per_pixel,opts->detector->nx,opts->detector->ny,"electrons_per_pixel.vtk");
  calculate_real_detector_output(opts->detector,opts->experiment);
  write_2D_array_to_vtk(opts->detector->real_output,opts->detector->nx/opts->detector->binning,
			opts->detector->ny/opts->detector->binning,"real_output.vtk");
  calculate_noiseless_detector_output(opts->detector,opts->experiment);
  write_2D_array_to_vtk(opts->detector->noiseless_output,opts->detector->nx/opts->detector->binning,
			opts->detector->ny/opts->detector->binning,"noiseless_output.vtk");
  noiseless = create_new_img(opts->detector->nx/opts->detector->binning,opts->detector->ny/opts->detector->binning);
  for(i = 0;i<TSIZE(noiseless);i++){
    noiseless->image[i] = opts->detector->noiseless_output[i];
    noiseless->mask[i] = 1;
  }
  write_img(noiseless,"noiseless_output.h5",sizeof(real));
#ifdef MPI
  MPI_Finalize();
#endif
  return 0;
}
