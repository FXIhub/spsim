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


#ifndef _SPSIM_H_
#define _SPSIM_H_ 1

#include "config.h"
#include "diffraction.h"
#include "molecule.h"
#include "io.h"
#include "mpi_comm.h"
#include "noise.h"
#include "amplification.h"
#include "real_space.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

void gaussian_blur_real_space(Options * opts,Diffraction_Pattern * pattern);
void gaussian_blur_pattern(Options * opts,Diffraction_Pattern * pattern);
Molecule * get_molecule(Options * opts);
float * get_HKL_list(Options * opts, int * HKL_list_size);
Diffraction_Pattern * compute_sf(Molecule * mol, float * HKL_list, int HKL_list_size, Options * opts);
void output_files( Diffraction_Pattern * pattern, Options * opts, int index);
Diffraction_Pattern * simulate_shot(Molecule * mol, Options * opts);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */


#endif 
