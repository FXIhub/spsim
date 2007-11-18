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


#include <time.h>                      // define time()
#include "randomc.h"                   // define classes for random number generators
#include "mersenne.cpp"                // code for random number generator
#define RANDOM_GENERATOR TRandomMersenne // define which random number generator to use
#include "stocc.h"                     // define random library classes
#include "stoc1.cpp"                   // random library source code
#include "stoc2.cpp"                   // random library source code
#include "userintf.cpp"                // define system specific user interface

#include "config.h"
#include "noise.h"

static StochasticLib2 * sto = NULL;

void init_random_generator(){
  if(!sto){
    int seed = time(0);                // random seed
    sto = new StochasticLib2(seed);   // make instance of random library
  }
}

int get_poisson_random_number(double L){
  return sto->Poisson(L);  
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
