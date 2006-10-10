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



#include "diffraction.h"

#ifdef MPI
#include <mpi.h>
#endif

void get_my_loop_start_and_end(int size, int * start, int * end);
void syncronize_patterns(Diffraction_Pattern * pat);

#ifdef MPI

int is_mpi_master();
void get_id_loop_start_and_end(int id, int size, int * start, int * end);
int mpi_receive_pattern(Diffraction_Pattern * pat);
int mpi_send_pattern(Diffraction_Pattern * pat);

#endif
