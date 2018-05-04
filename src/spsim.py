# SPSIM PYTHON MODULE

__version__ = "0.1.0"

# Wrapped C funcitons
from spsim_pybackend import *

import numpy

try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen
from io import StringIO
from io import BytesIO
import gzip

# Convenience functions
def get_molecule_from_opts(opts):
    mol = get_molecule(opts)
    return mol
    
def get_molecule_from_atoms(atomic_numbers = None, atomic_positions = None):
    mol = alloc_mol()
    for j,(x,y,z) in zip(atomic_numbers, atomic_positions):
        add_atom_to_mol(mol, int(j), float(x), float(y), float(z))
    return mol

def get_atoms_from_molecule(mol):
    pos_img = sp_image_alloc(3, mol.natoms, 1)
    array_to_image(mol.pos, pos_img)
    pos = numpy.array(pos_img.image.real, dtype=numpy.float64)
    pos = pos.reshape((mol.natoms,3))
    #print "from molecule: Lz=%e Ly=%e Lx=%e" % (pos[:,1].max()-pos[:,1].min(),pos[:,2].max()-pos[:,2].min())
    anum_img = sp_image_alloc(mol.natoms, 1, 1)
    iarray_to_image(mol.atomic_number, anum_img)
    anum = numpy.array(anum_img.image.real, dtype=numpy.int32)
    anum = anum.reshape(anum.size)
    sp_image_free(pos_img)
    sp_image_free(anum_img)
    return anum, pos

def fetch_pdb(pdb_id):
    url = "http://www.rcsb.org/pdb/files/%s.pdb.gz" % str(pdb_id)
    filename = "./%s.pdb" % str(pdb_id)
    response = urlopen(url)
    compressedFile = BytesIO()
    compressedFile.write(response.read())
    compressedFile.seek(0)
    decompressedFile = gzip.GzipFile(fileobj=compressedFile, mode='rb')
    with open(filename, 'wb') as outfile:
        outfile.write(decompressedFile.read())
    return filename
