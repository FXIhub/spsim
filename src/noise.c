#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <spimage.h>
#include "config.h"

static gsl_rng * r = NULL;

void init_random_generator(){
  const gsl_rng_type * T;
  T = gsl_rng_default;
  if(!r){
    r = gsl_rng_alloc (T);    
  }
}


int get_poisson_random_number(double L){
  return gsl_ran_poisson(r,L);
}


void generate_poisson_noise(CCD * det){
  int i;
  init_random_generator();
  if(!det->photons_per_pixel){
    fprintf(stderr,"Please calculate photons per pixel first\n");
    return;
  }
  det->photon_count = (float *)malloc(sizeof(float)*det->nx*det->ny*det->nz);
  for(i = 0;i<det->nx*det->ny*det->nz;i++){
    det->photon_count[i] = get_poisson_random_number(det->photons_per_pixel[i]*det->quantum_efficiency);
  }  
}

void generate_gaussian_noise(CCD * det){
  int i;
  if(!det->photons_per_pixel){
    fprintf(stderr,"Please calculate photons per pixel first\n");
    return;
  }
  det->photon_count = (float *)malloc(sizeof(float)*det->nx*det->ny*det->nz);
  for(i = 0;i<det->nx*det->ny*det->nz;i++){
    det->photon_count[i] = ((p_drand48()-0.5))*sqrt(det->photons_per_pixel[i]*det->quantum_efficiency)+det->photons_per_pixel[i]*det->quantum_efficiency;
  }  
}
