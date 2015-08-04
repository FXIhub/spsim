# SPSIM PYTHON MODULE

# Wrapped C funcitons
from spsim_pybackend import *

import numpy

# Convenience functions
def get_molecule_from_opts(opts):
    mol = get_molecule(opts)
    return mol
    
def get_molecule_from_atoms(atomic_numbers = None, atomic_positions = None):
    mol = alloc_mol()
    for j,(pos0,pos1,pos2) in zip(atomic_numbers, atomic_positions):
        add_atom_to_mol(mol, int(j), float(pos0), float(pos1), float(pos2))
    return mol

def get_atoms_from_molecule(mol):
    pos_img = sp_image_alloc(mol.natoms, 3, 1)
    array_to_image(mol.pos, pos_img)
    temp = numpy.array(pos_img.image.real, dtype=numpy.float64)
    temp = temp.reshape((3,mol.natoms))
    pos = numpy.zeros(shape=(mol.natoms,3), dtype=numpy.float64)
    pos[:,0] = temp[0,:]
    pos[:,1] = temp[1,:]
    pos[:,2] = temp[2,:]
    anum_img = sp_image_alloc(mol.natoms, 1, 1)
    iarray_to_image(mol.atomic_number, anum_img)
    anum = numpy.array(anum_img.image.real, dtype=numpy.int32)
    anum = anum.reshape(anum.size)
    sp_image_free(pos_img)
    sp_image_free(anum_img)
    return anum, pos
