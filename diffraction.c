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

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <spimage.h>
#include "config.h"
#include "diffraction.h"
#include "mpi.h"
#include "io.h"


#ifndef M_PI
#define M_PI     3.1415926535897932384626433832795029L
#endif

#define ELEMENTS 100
static float atomsf[ELEMENTS][9] = 
#include "atomsf.cdata"

static int atomsf_initialized = 0;


/*
  Use atomsf.lib to obtain details for all atoms in a problem.
*/

/*
static void write_ff_tables(){
  int i;
  FILE * fp = fopen("atomsf.data","w");
  fprintf(fp," {\n");
  i = 0;
  fprintf(fp,"{%e,%e,%e,%e,%e,%e,%e,%e,%e}\n",atomsf[i][0],
	  atomsf[i][1],atomsf[i][2],atomsf[i][3],
	  atomsf[i][4],atomsf[i][5],atomsf[i][6],
	  atomsf[i][7],atomsf[i][8]);
  for(i =1;i<ELEMENTS;i++){
    fprintf(fp,",{%e,%e,%e,%e,%e,%e,%e,%e,%e}\n",atomsf[i][0],
	    atomsf[i][1],atomsf[i][2],atomsf[i][3],
	    atomsf[i][4],atomsf[i][5],atomsf[i][6],
	    atomsf[i][7],atomsf[i][8]);
  };
  fprintf(fp,"};\n");   
  fclose(fp);
}
*/

static void fill_ff_tables(){
  FILE * ff;
  char line[1024];
  char * p;
  int z;
  int i;
  ff = fopen("atomsf.lib","r");
  if(!ff){
/*    fprintf(stderr,"Info: atomsf.lib not found in current directory, using internal atomsf.lib\n");*/
    return;
  }
  
  /* Init atomsf array */
  for(i = 0;i<ELEMENTS;i++){
    atomsf[i][0] = -1;
  }
  /* skip the comments */
  for(fgets(line,1024,ff);strstr(line,"AD") == line;fgets(line,1024,ff));

  /* we're at the beginning of an atom record */
  while(line[0]){
    fgets(line,1024,ff);
    /* get Z */
    sscanf(line,"%d",&z);
    if(atomsf[z][0] == -1){
      p = line+23;
      sscanf(p,"%f",&atomsf[z][8]);
      fgets(line,1024,ff);
      p = line;
      sscanf(p,"%f",&atomsf[z][0]);
      p+= 16; 
      sscanf(p,"%f",&atomsf[z][1]);
      p+= 16; 
      sscanf(p,"%f",&atomsf[z][2]);
      p+= 16; 
      sscanf(p,"%f",&atomsf[z][3]);
      fgets(line,1024,ff);
      p = line;
      sscanf(p,"%f",&atomsf[z][4]);
      p+= 16; 
      sscanf(p,"%f",&atomsf[z][5]);
      p+= 16; 
      sscanf(p,"%f",&atomsf[z][6]);
      p+= 16; 
      sscanf(p,"%f",&atomsf[z][7]);
      /* get the last line of the atom record */
      fgets(line,1024,ff);
      /* get the first line of the next atom record */
      line[0] = 0;
      fgets(line,1024,ff);      
    }else{
      /* skip record */
      fgets(line,1024,ff);      
      fgets(line,1024,ff);      
      fgets(line,1024,ff);      
      line[0] = 0;
      fgets(line,1024,ff);            
    }
  }
}



/* d should be the size of the scattering vector |H| */
/* d is in m^-1 but we need it in A^-1 so divide  by 1^10*/
static float  scatt_factor(float d,int Z){
  float res = 0;
  int i;  
  d *= 1e-10;
  /* the 0.25 is there because the formula is (|H|/2)^2 */
  for(i = 0;i<4;i++){
    res+= atomsf[Z][i]*exp((-atomsf[Z][i+4])*d*d*0.25);
  }                
  res += atomsf[Z][8];
  return res;    
}


float * get_HKL_list_for_detector(CCD * det, Experiment * exp,int * HKL_list_size){
  /* number of pixels */
  int nx, ny;
  /* pixel index */
  int x,y;
  /* physical location of pixel*/
  double px,py;
  /* reciprocal coordinates */
  double rx,ry;
  double real_to_reciprocal = 1.0/(det->distance*exp->wavelength);
  double ewald_radius = 1.0/exp->wavelength;
  double distance_to_ewald_sphere_center;

  float * HKL_list;
  int index = 0;
  
  nx = det->nx;
  ny = det->ny;

  HKL_list = malloc(sizeof(float)*nx*ny*3);
  for(x = 0;x<nx;x++){
    for(y = 0;y<ny;y++){
      /* 
	 Calculate the pixel coordinates in reciprocal space 	 
	 by dividing the physical position by detector_distance*wavelength.
	 
	 CCD center at (nx-1)/2,(ny-1)/2

	 Upper left corner of the detector with negative x and positive y
      */
      px = ((x-(nx-1.0)/2.0)/nx)*det->width;
      py = (((ny-1.0)/2.0-y)/ny)*det->height;

      rx = px*real_to_reciprocal;
      ry = py*real_to_reciprocal;
      
      /* Project pixel into Ewald sphere. */
      if(!det->spherical){
	distance_to_ewald_sphere_center = sqrt(rx*rx+ry*ry+ewald_radius*ewald_radius);
	HKL_list[index++] = rx * ewald_radius/distance_to_ewald_sphere_center;
	HKL_list[index++] = ry * ewald_radius/distance_to_ewald_sphere_center;
	HKL_list[index++] = ewald_radius-(ewald_radius * ewald_radius/distance_to_ewald_sphere_center);
      }else{
	HKL_list[index++] = rx;
	HKL_list[index++] = ry;
	HKL_list[index++] = 0;
      }
	 
    }
  }
  *HKL_list_size = nx*ny;
  return HKL_list;
}

Diffraction_Pattern * compute_pattern_on_list(Molecule * mol, float * HKL_list, int HKL_list_size){
  int i,j;
  double scattering_factor;
  double scattering_vector_length;
  double scattering_factor_cache[ELEMENTS];
  int is_element_in_molecule[ELEMENTS];
  Diffraction_Pattern * res = malloc(sizeof(Diffraction_Pattern));
  int HKL_list_start = 0;
  int HKL_list_end = 0;
  int points_per_percent;
  get_my_loop_start_and_end(HKL_list_size,&HKL_list_start,&HKL_list_end);

  if(!atomsf_initialized){
    fill_ff_tables();
    atomsf_initialized = 1;
  }

  res->F = malloc(sizeof(fftw_complex)*HKL_list_size);
  res->ints = malloc(sizeof(float)*HKL_list_size);
  res->HKL_list = malloc(sizeof(float)*3*HKL_list_size);
  memcpy(res->HKL_list,HKL_list,sizeof(float)*3*HKL_list_size);
  res->HKL_list_size = HKL_list_size;
  for(j = 0 ;j< ELEMENTS;j++){
    is_element_in_molecule[j] = 0;
  }
  for(j = 0 ;j< mol->natoms;j++){
    is_element_in_molecule[mol->atomic_number[j]] = 1;
  }

  points_per_percent = (HKL_list_end-HKL_list_start)/100;
  for(i = HKL_list_start;i<HKL_list_end;i++){
#ifdef MPI    
    if(is_mpi_master()){
      if(i % points_per_percent == 0){
	fprintf(stderr,"%f percent done\n",(100.0*(i-HKL_list_start))/(HKL_list_end-HKL_list_start));
      }
    }
#else
      if(i % points_per_percent == 0){
	fprintf(stderr,"%f percent done\n",(100.0*(i-HKL_list_start))/(HKL_list_end-HKL_list_start));
      }

#endif

    res->F[i] = 0;
    res->F[i] = 0;
    scattering_vector_length = sqrt(HKL_list[3*i]*HKL_list[3*i]+HKL_list[3*i+1]*HKL_list[3*i+1]+HKL_list[3*i+2]*HKL_list[3*i+2]);
    for(j = 0;j<ELEMENTS;j++){
      if(is_element_in_molecule[j]){
	scattering_factor_cache[j] = scatt_factor(scattering_vector_length,j);
      }
    }
    for(j = 0 ;j< mol->natoms;j++){
      scattering_factor = scattering_factor_cache[mol->atomic_number[j]];
      res->F[i] += scattering_factor*cos(2*M_PI*(HKL_list[3*i]*mol->pos[j*3]+HKL_list[3*i+1]*mol->pos[j*3+1]+HKL_list[3*i+2]*mol->pos[j*3+2]));
      res->F[i] += I*scattering_factor*sin(2*M_PI*(HKL_list[3*i]*mol->pos[j*3]+HKL_list[3*i+1]*mol->pos[j*3+1]+HKL_list[3*i+2]*mol->pos[j*3+2]));
    }
    res->ints[i] = cabsr(res->F[i])*cabsr(res->F[i]);
  }
  syncronize_patterns(res);
  return res;
}



Diffraction_Pattern * load_pattern_from_file(CCD * det,char * filename, 
					     float * HKL_list, int HKL_list_size){
  Diffraction_Pattern * res = malloc(sizeof(Diffraction_Pattern));

  res->HKL_list = malloc(sizeof(float)*3*HKL_list_size);
  memcpy(res->HKL_list,HKL_list,sizeof(float)*3*HKL_list_size);
  res->HKL_list_size = HKL_list_size;

  res->ints = read_VTK_to_array(det->nx,det->ny,filename);
  res->F = 0;
  return res;
}

void calculate_thomson_correction(CCD * det){
  int x,y;
  double px,py;
  int nx = det->nx;
  int ny = det->ny;
  double r; /* distance to scatterer */
  int index;
  double r0 = 2.81794e-15; /* classical electron radius = e^2/(m*c^2)*/
  /* For the moment we're considering vertical polarization */
  double polarization_factor = 1;
  det->thomson_correction = malloc(sizeof(float)*nx*ny);
  index = 0;
  for(x = 0;x<det->nx;x++){
    for(y = 0;y<det->ny;y++){
      px = ((x-(nx-1.0)/2.0)/nx)*det->width/2;
      py = (((ny-1.0)/2.0-y)/ny)*det->height/2;
      r = sqrt(det->distance*det->distance+px*px+py*py);
      det->thomson_correction[index++] = (r0*r0)*polarization_factor;
    }
  }                  
}

void calculate_pixel_solid_angle(CCD * det){
  /* 
     the effective pixel area is calculated by
     the projection of the pixel on a sphere centered
     on the scatterer and going through the center of the pixel
  */
  int i;
  int x,y;
  double px,py;
  double corners[4][2];
  double corner_distance[4];
  double projected_pixel_sides[4];
  int nx = det->nx;
  int ny = det->ny;
  double r; /* distance to scatterer */
  int index;
  /* For the moment we're considering vertical polarization */
  det->solid_angle = malloc(sizeof(float)*nx*ny);
  index = 0;
  for(x = 0;x<det->nx;x++){
    for(y = 0;y<det->ny;y++){
      if(det->spherical){
	r = det->distance;
	det->solid_angle[index++] = det->pixel_width*det->pixel_height/(r*r);
      }else{
	px = ((x-(nx-1.0)/2.0)/nx)*det->width/2;
	py = (((ny-1.0)/2.0-y)/ny)*det->height/2;
	r = sqrt(det->distance*det->distance+px*px+py*py);
	/* top left */
	corners[0][0] = px-det->pixel_width/2;
	corners[0][1] = py+det->pixel_height/2;
	
	corner_distance[0] = sqrt(det->distance*det->distance+corners[0][0]*corners[0][0]+corners[0][1]*corners[0][1]);
	/* top right */
	corners[1][0] = px+det->pixel_width/2;
	corners[1][1] = py+det->pixel_height/2;
	corner_distance[1] = sqrt(det->distance*det->distance+corners[1][0]*corners[1][0]+corners[1][1]*corners[1][1]);
	/* bottom right */
	corners[2][0] = px+det->pixel_width/2;
	corners[2][1] = py-det->pixel_height/2;
	corner_distance[2] = sqrt(det->distance*det->distance+corners[2][0]*corners[2][0]+corners[2][1]*corners[2][1]);
	/* bottom left */
	corners[3][0] = px-det->pixel_width/2;
	corners[3][1] = py-det->pixel_height/2;
	corner_distance[3] = sqrt(det->distance*det->distance+corners[3][0]*corners[3][0]+corners[3][1]*corners[3][1]);
	/* project on plane*/
	if(!det->spherical){
	  for(i = 0;i<4;i++){
	    corners[i][0] *= r/corner_distance[i];
	    corners[i][1] *= r/corner_distance[i];	
	  }
	}
	/* top */
	projected_pixel_sides[0] = sqrt((corners[0][0]-corners[1][0])*(corners[0][0]-corners[1][0])+(corners[0][1]-corners[1][1])*(corners[0][1]-corners[1][1]));
	/* left */
	projected_pixel_sides[1] = sqrt((corners[0][0]-corners[3][0])*(corners[0][0]-corners[3][0])+(corners[0][1]-corners[3][1])*(corners[0][1]-corners[3][1]));
	det->solid_angle[index++] = projected_pixel_sides[0]*projected_pixel_sides[1]/(r*r);
      }
    }
  }
}                    


void calculate_photons_per_pixel(Diffraction_Pattern * pattern, CCD * det, Experiment * experiment){
  int i;
  if(!det->thomson_correction){
    calculate_thomson_correction(det);
  }
  if(!det->solid_angle){
    calculate_pixel_solid_angle(det);
  }
  det->photons_per_pixel = malloc(sizeof(float)*det->nx*det->ny);
  for(i = 0;i<det->nx*det->ny;i++){
    det->photons_per_pixel[i] = det->thomson_correction[i] * det->solid_angle[i] * pattern->ints[i] * experiment->beam_intensity;
  }
}


void write_hkl_grid(float * list, Molecule * mol,CCD * det){
  Image * hkl_grid;
  int x,y;
  int i;
  int index;
  double min_x = FLT_MAX;
  double min_y = FLT_MAX; 
  double max_x = -FLT_MAX;
  double max_y = -FLT_MAX;
  double max_dim;
#ifdef MPI
  if(!is_mpi_master()){
    return;
  }
#endif
  for(i =0 ;i<mol->natoms;i++){
    if(mol->pos[i*3] < min_x){
      min_x = mol->pos[i*3];
    }
    if(mol->pos[i*3] > max_x){
      max_x = mol->pos[i*3];
    }
    if(mol->pos[i*3+1] < min_y){
      min_y = mol->pos[i*3+1];
    }
    if(mol->pos[i*3+1] > max_y){
      max_y = mol->pos[i*3+1];
    }
  }
  max_dim = sp_max(max_x-min_x,max_y-min_y);
  hkl_grid = sp_image_alloc(det->nx/det->binning,det->ny/det->binning);
  hkl_grid->detector->image_center[0] = sp_image_width(hkl_grid)/2;
  hkl_grid->detector->image_center[1] = sp_image_height(hkl_grid)/2;
  i = 0;
  index = 0;
  for(x = 0;x<sp_image_width(hkl_grid);x++){
    for(y = 0;y<sp_image_height(hkl_grid);y++){
      if(fabs(x-hkl_grid->detector->image_center[0]) > fabs(y-hkl_grid->detector->image_center[1])){
	hkl_grid->image->data[i] = list[index]*max_dim;
      }else{
	hkl_grid->image->data[i] = list[index+1]*max_dim;
      }
      hkl_grid->mask->data[i] = 1;
      index += 3*det->binning;
      i++;
    }
    index += 3*(det->binning-1)*(det->ny);
  }
  sp_image_write(hkl_grid,"hkl_grid.h5",sizeof(real));
}
