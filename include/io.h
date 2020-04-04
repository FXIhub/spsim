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

#include "config.h"

void write_2D_array_to_vtk(float * f,int nx, int ny, char * filename);
void write_3D_array_to_vtk(float * f,int nx, int ny,int nz, char * filename);
float * read_VTK_to_array(int nx, int ny, int nz,char * filename);
void write_3D_array_as_structured_grid_to_vtk(float * F,CCD * det, float * HKL_list, int HKL_list_size, int n_patterns,char * filename);
