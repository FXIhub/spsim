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

#include "spimage.h"

Molecule * get_Molecule_from_formula(Chem_Formula * form, Options * opts);
Molecule * get_Molecule_from_pdb(char * filename);
void write_pdb_from_mol(char *filename,Molecule * mol);

Molecule * alloc_mol();
Molecule * make_mol(Image * atomic_number, Image * pos);
void add_atom_to_mol(Molecule * mol, int atomic_number, float x, float y, float z);
void free_mol(Molecule * mol);

void origin_to_center_of_mass(Molecule * mol);
