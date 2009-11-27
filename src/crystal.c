#include <math.h>
#include "crystal.h"

void crystal_axis_to_cartesian(Options * opts,CrystalAxis axis,int * x,int * y,int *z){
  const float pi = 3.14159265;
  const float grad_to_rad = pi/180.0;
  if(axis == CrystalA){
    *x = opts->crystal_cell[0];
    *y = 0;
    *z = 0;
  }else if(axis == CrystalB){
    *x = cos(opts->crystal_cell[5]*grad_to_rad)*opts->crystal_cell[1];
    *y = sin(opts->crystal_cell[5]*grad_to_rad)*opts->crystal_cell[1];
    *z = 0;    
  }else if(axis == CrystalC){
    *x  = opts->crystal_cell[2]*cos(opts->crystal_cell[4]*grad_to_rad);
    *y = opts->crystal_cell[2]*(cos(opts->crystal_cell[3]*grad_to_rad)-cos(opts->crystal_cell[4]*grad_to_rad)*cos(opts->crystal_cell[5]*grad_to_rad))
      /sin(opts->crystal_cell[5]*grad_to_rad);
    *z = sqrt(opts->crystal_cell[2]*opts->crystal_cell[2]-(*x)*(*x)-(*y)*(*y));
  }
}
