#include <math.h>
#include "crystal.h"

void crystal_axis_to_cartesian(Options * opts,CrystalAxis axis,float * x,float * y,float *z){
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

void crystal_cell_matrix(Options * opts,float * matrix){
  const float pi = 3.14159265;
  const float grad_to_rad = pi/180.0;

  matrix[0] = opts->crystal_cell[0];
  matrix[1] = 0;
  matrix[2] = 0;
  matrix[3] = cos(opts->crystal_cell[5]*grad_to_rad)*opts->crystal_cell[1];
  matrix[4] = sin(opts->crystal_cell[5]*grad_to_rad)*opts->crystal_cell[1];
  matrix[5] = 0;    
  matrix[6]  = opts->crystal_cell[2]*cos(opts->crystal_cell[4]*grad_to_rad);
  matrix[7] = opts->crystal_cell[2]*(cos(opts->crystal_cell[3]*grad_to_rad)-cos(opts->crystal_cell[4]*grad_to_rad)*cos(opts->crystal_cell[5]*grad_to_rad))
    /sin(opts->crystal_cell[5]*grad_to_rad);
  matrix[8] = sqrt(opts->crystal_cell[2]*opts->crystal_cell[2]-(matrix[6])*(matrix[6])-(matrix[7])*(matrix[7]));
}

Complex * calculate_pattern_from_crystal(Complex * F,float * HKL_list, int HKL_list_size,Options * opts){
  float cell[9];
  Complex * cF = (Complex *)malloc(sizeof(Complex)*HKL_list_size);
  crystal_cell_matrix(opts,cell);
  for(int i =0 ;i<HKL_list_size;i++){
    cF[i] = sp_cinit(0,0);
  }
  float dr[3];
  for(int a = 0;a<opts->crystal_size[0];a++){
    for(int b = 0;b<opts->crystal_size[1];b++){
      for(int c = 0;c<opts->crystal_size[2];c++){
	dr[0] = cell[0]*a+cell[3]*b+cell[6]*c;
	dr[1] = cell[1]*a+cell[4]*b+cell[7]*c;
	dr[2] = cell[2]*a+cell[5]*b+cell[8]*c;
	printf(".\n");
	for(int i =0 ;i<HKL_list_size;i++){	  
	  float phi = 2*M_PI*(HKL_list[3*i]*-dr[0]+HKL_list[3*i+1]*-dr[1]+HKL_list[3*i+2]*-dr[2]) + sp_carg(F[i]);
	  float amp = sp_cabs(F[i]);
	  sp_real(cF[i]) += cos(phi)*amp;
	  sp_imag(cF[i]) += sin(phi)*amp;
	}
      }
    }
  }
  return cF;
}
