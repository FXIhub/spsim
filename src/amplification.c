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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <spimage.h>
#include "config.h"
#include "amplification.h"



/* mean m, standard deviation s */
/* gaussian distribution generator */
float box_muller(float m, float s){        
  float x1, x2, w, y1;
  static float y2;
  static int use_last = 0;
  
  if (use_last){        /* use value from previous call */
    y1 = y2;
    use_last = 0;
  }  else {
    do {
      x1 = 2.0 * p_drand48() - 1.0;
      x2 = 2.0 * p_drand48() - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;
    use_last = 1;
  }  
  return( m + y1 * s );
}


void calculate_electrons_per_pixel(CCD * det, Experiment * experiment){
  int i;
  double h_c = 1.98644521e-25; /* Planck's constant * the speed of light */
  double photon_energy = h_c/experiment->wavelength; /* photon energy in joules */
  float electrons_per_photon = 1;
  if(det->electron_hole_production_energy){
    electrons_per_photon = photon_energy/det->electron_hole_production_energy;
  }
  if(!det->photon_count){
    fprintf(stderr,"Calculate photon_count first\n");
    return;
  }
  det->electrons_per_pixel = malloc(sizeof(float)*det->nx*det->ny*det->nz);
  for(i = 0;i<det->nx*det->ny*det->nz;i++){
    det->electrons_per_pixel[i] =  det->photon_count[i]*electrons_per_photon;
  }
}


/* This functions calculates the real detector output taking into account
 readout noise, thermal noise, linear full well and maximum detector value */

/* It doesn't take into account electron spill */
void calculate_real_detector_output(CCD * det, Experiment * experiment){
  int i;
  float ADC_constant = 1;
  if(det->linear_full_well && det->maximum_value){
    ADC_constant = det->maximum_value/det->linear_full_well;
  }
  int x,y,z,xi,yi,zi;
  if(!det->electrons_per_pixel){
    fprintf(stderr,"Calculate electrons_per_pixel first\n");
    return;
  }
  det->real_output = malloc(sizeof(float)*(det->nx/det->binning_x)*(det->ny/det->binning_y)*(det->nz/det->binning_z));
  i = 0;
  for(z = 0;z<det->nz/det->binning_z;z++){
    for(y = 0;y<det->ny/det->binning_y;y++){
      for(x = 0;x<det->nx/det->binning_x;x++){
	det->real_output[i] = 0;
	for(zi = 0;zi<det->binning_z;zi++){
	  for(yi = 0;yi<det->binning_y;yi++){
	    for(xi = 0;xi<det->binning_x;xi++){
	      det->real_output[i] +=  (det->electrons_per_pixel[(z*det->binning_z+zi)*(det->ny)*(det->nx)+
								(y*det->binning_y+yi)*(det->nx)+
				                                (x*det->binning_x+xi)]+
				       det->dark_current*experiment->exposure_time);
	    }
	  }
	}
	det->real_output[i] +=+box_muller(0,det->readout_noise);
	if(det->real_output[i] < 0){
	  det->real_output[i] = 0;
	}
	det->real_output[i] *= ADC_constant;
	i++;
      }
    }
  }

}


/* This function calculates the noiseless detector output directly from
 the photons_per_pixel value */
void calculate_noiseless_detector_output(CCD * det, Experiment * experiment){
  int i;
  double h_c = 1.98644521e-25; /* Planck's constant * the speed of light */
  double photon_energy = h_c/experiment->wavelength; /* photon energy in joules */
  float electrons_per_photon = photon_energy/det->electron_hole_production_energy;
  float ADC_constant = 1;
  if(det->linear_full_well && det->maximum_value){
    ADC_constant = det->maximum_value/det->linear_full_well;
  }
  int x,y,z,xi,yi,zi;
  if(!det->photons_per_pixel){
    fprintf(stderr,"Calculate photons_per_pixel first\n");
    return;
  }
  det->noiseless_output = malloc(sizeof(float)*(det->nx/det->binning_x)*(det->ny/det->binning_y)*(det->nz/det->binning_z));
  i = 0;
  for(z = 0;z<det->nz/det->binning_z;z++){
    for(y = 0;y<det->ny/det->binning_y;y++){
      for(x = 0;x<det->nx/det->binning_x;x++){    
	det->noiseless_output[i] = 0;
	for(zi = 0;zi<det->binning_z;zi++){
	  for(yi = 0;yi<det->binning_y;yi++){
	    for(xi = 0;xi<det->binning_x;xi++){
	      det->noiseless_output[i] +=  (det->photons_per_pixel[(z*det->binning_z+zi)*det->ny*det->nx+
								   (y*det->binning_y+yi)*det->nx+
								   (x*det->binning_x+xi)]
					    *det->quantum_efficiency*electrons_per_photon*ADC_constant);
	    }
	  }
	}
	i++;
      }
    }
  }
}
