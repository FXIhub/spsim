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
#ifndef _CRYSTAL_H_
#define _CRYSTAL_H_ 1

#include "config.h"
#include <spimage.h>

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

typedef enum{CrystalA,CrystalB,CrystalC}CrystalAxis;
void crystal_axis_to_cartesian(Options * opts,CrystalAxis axis,float * x,float * y,float *z);
void crystal_cell_matrix(Options * opts,float * matrix);
  void calculate_pattern_from_crystal(float * I,Complex * F,float * HKL_list, int HKL_list_size,Options * opts);
  void calculate_pattern_from_crystal_cuda(float * d_I, cufftComplex * d_F,float * d_HKL_list, int HKL_list_size,Options * opts);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
