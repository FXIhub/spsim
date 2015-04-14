%module spsim_pybackend
%{
  /* Includes the header in the wrapper code */
#include "../include/spsim.h"
  /*
#include "../include/box.h"
#include "../include/config.h"
#include "../include/crystal.h"
#include "../include/diffraction.h"
#include "../include/fftw3.h"
#include "../include/io.h"
#include "../include/molecule.h"
#include "../include/mpi_comm.h"
#include "../include/noise.h"
#include "../include/real_space.h"
  */
  %}

%include "../include/spsim.h"
