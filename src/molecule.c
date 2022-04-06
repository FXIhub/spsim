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
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "config.h"
#include "molecule.h"
#include "mpi_comm.h"

static char legal_atom_names[212];
static void get_legal_atom_names();
static int getZfromSymbol(char * symbol);

Molecule * get_Molecule_from_formula(Chem_Formula * form, Options * opts){
  int i,j,d;
  Molecule * res = malloc(sizeof(Molecule));
  float distance;
  res->atomic_number = NULL;
  res->pos = NULL;
  res->natoms = 0;
  for(i = 0;i<form->n_elements;i++){
    res->natoms += form->quantity[i];
    res->atomic_number = realloc(res->atomic_number,sizeof(int)*res->natoms);
    res->pos = realloc(res->pos,sizeof(float)*res->natoms*3);
    for(j = res->natoms-form->quantity[i];j<res->natoms;j++){
      res->atomic_number[j] = form->atomic_number[i];
      if(opts->box_type == BOX_SPHERICAL){
      do{
        distance = 0;
        for(d = 0;d<3;d++){
          res->pos[3*j+d] = (p_drand48()-0.5)*opts->box_dimension;
          distance += res->pos[3*j+d]*res->pos[3*j+d];
        }
        distance = sqrt(distance);
      }while(distance > opts->box_dimension/2.0);
      }else if(opts->box_type == BOX_PARALLEL){
      for(d = 0;d<3;d++){
        res->pos[3*j+d] = p_drand48()*opts->box_dimension;
      }      
      }
    }
  }
  return res;  
}

Molecule * get_Molecule_from_pdb(char * filename){
  FILE * fp = fopen(filename,"r");
  char buffer[1024];
  int maxnatoms = 1024;
  char      gid[10], aid[10] ,atid[3];
  char  element_buffer[2];
  char       tmp_char;
  char       *next_field_start;
  float t4,t5;
  float x,y,z;
  int Z;
  int total_atomic_number = 0;
  float *rot_mats, *trans_vecs;
  int elem_num, dim_num;
  int num_symm = 0, max_num_symm = 1;
  Molecule * res = malloc(sizeof(Molecule));
  if(!fp){
    perror("Error reading PDB\n");
    abort();
  }
  res->atomic_number = malloc(sizeof(int)*maxnatoms);
  res->pos = malloc(3*sizeof(float)*maxnatoms);
  res->natoms = 0;
  rot_mats = malloc(9*max_num_symm*sizeof(float));
  trans_vecs = malloc(3*max_num_symm*sizeof(float));

  get_legal_atom_names();

  // Loop over lines in file
  while (fgets(buffer, 1024, fp) != NULL)  {
    // Only process atom, hetatm and biomt lines
    if ((int)strlen(buffer) <= 1) {
      continue ;
    }
    else if ((strstr(buffer, "ATOM") == buffer) || (strstr(buffer, "HETATM") == buffer)) {
      // Increase arrays if needed
      if (res->natoms >= maxnatoms) {
      /* increase array sizes */
      maxnatoms *= 2;
      //printf("%i\n", maxnatoms);
      res->atomic_number = realloc(res->atomic_number,sizeof(int)*maxnatoms);
      res->pos = realloc(res->pos,3*sizeof(float)*maxnatoms);
      }
      
      /**************************************************
                 Original code had:
                 sscanf(buffer, "%*s %d %s %s %d %g %g %g %g %g", 
                        &atno[Nin], aid, gid, 
                        &groupno[Nin], &ud->t1, &ud->t2, &ud->t3, &ud->t4, &ud->t5) ;

                 Filipe Maia: Extract the fixed width fields from 
                 the pdb, according to 
                 http://www.ccp4.ac.uk/dist/html/pdbformat.html
                 (which he believes is exactly the same as the formal 
                 format description).  
                 The atof/atoi should be changed to strtof/strtol to 
                 detect input errors.
      **************************************************/

      tmp_char = buffer[11];            /* limit the field*/
      buffer[11] = 0;
      next_field_start = &(buffer[6]);
      /*             res->pdb_atom_number[res->natoms] = atoi(next_field_start);*/
      buffer[11] = tmp_char;            /* advance field pointer */
      next_field_start+= 6;
      tmp_char = buffer[16];
      buffer[16] = 0;

      /* Next: skip over alternate location indicator */
      next_field_start+= 5;
      tmp_char = buffer[20];
      buffer[20] = 0;
      strncpy(gid,next_field_start,4);
      gid[4] = 0;
      
      /* Next: skip over chain identifier */
      next_field_start+= 5;
      tmp_char = buffer[26];
      buffer[26] = 0;
      /*             res->pdb_groupno[res->natoms] = atoi(next_field_start);*/
      buffer[26] = tmp_char;
      
      /* Start retrieving coordinates*/
      next_field_start += 8;
      tmp_char = buffer[38];
      buffer[38] = 0;
      x = (float)atof(next_field_start);
      buffer[38] = tmp_char;
      next_field_start += 8;
      tmp_char = buffer[46];
      buffer[46] = 0;
      y = (float)atof(next_field_start);
      buffer[46] = tmp_char;
      next_field_start += 8;
      tmp_char = buffer[54];
      buffer[54] = 0;
      z = (float)atof(next_field_start);
      
      /* Now retrieve Pdb->Occupancy and B */
      buffer[54] = tmp_char;
      next_field_start += 8;
      tmp_char = buffer[60];
      buffer[60] = 0;
      t4 = (float)atof(next_field_start);
      buffer[60] = tmp_char;
      next_field_start += 6;
      tmp_char = buffer[66];
      buffer[66] = 0;
      t5 = (float)atof(next_field_start);
      buffer[66] = tmp_char;
      
      /*               strncpy(pdb->atid[pdb->Nin], aid, 2) ; 
                   strncpy(pdb->fullatid[pdb->Nin], aid, 5) ; 
                   pdb->fullatid[pdb->Nin][4] = 0;
                   pdb->atid[pdb->Nin][2] = 0;*/
      /*             if(isspace(pdb->atid[pdb->Nin][0]) ||
                   isdigit(pdb->atid[pdb->Nin][0])){*/
      /* left justify */
      /*             pdb->atid[pdb->Nin][0] = pdb->atid[pdb->Nin][1];      
                   pdb->atid[pdb->Nin][1] = ' ';
                   }
                   strncpy(pdb->groupid[pdb->Nin], gid, 3) ; 
                   pdb->groupid[pdb->Nin][3] = 0;*/

      /* convert to meters */
      // right handed coordinate system
      res->pos[res->natoms*3+0] = x*1e-10 ;
      res->pos[res->natoms*3+1] = y*1e-10 ;
      res->pos[res->natoms*3+2] = z*1e-10 ;

      /* Retrieve element entry */
      next_field_start = &(buffer[77-1]);
      strncpy(element_buffer, next_field_start, sizeof(element_buffer));
      if(isspace(element_buffer[0])){
      element_buffer[0] = element_buffer[1];            /* left justify */
      element_buffer[1] = ' ';
      }
      Z = getZfromSymbol(element_buffer);
      if(!Z){
      fprintf(stderr,"WARNING: Null atom at line '%s'. Skipping\n",buffer);
      res->natoms--;
      }
      total_atomic_number += Z;
      res->atomic_number[res->natoms] = Z;
      res->natoms++;
    }
    else if (strstr(buffer, "REMARK 350   BIOMT") == buffer) {
      tmp_char = buffer[19];
      buffer[19] = 0;
      dim_num = atoi(&buffer[18]) - 1;
      buffer[19] = tmp_char;

      tmp_char = buffer[24];
      buffer[24] = 0;
      elem_num = atoi(&buffer[20]) - 2;
      buffer[24] = tmp_char;
      if (elem_num == -1) // Identity operation
        continue ;

      if (elem_num >= max_num_symm) { // Expand array if needed
        /* increase array sizes */
        max_num_symm *= 2;
        //printf("%i\n", maxnatoms);
        rot_mats = realloc(rot_mats, 9*max_num_symm*sizeof(float));
        trans_vecs = realloc(trans_vecs ,3*max_num_symm*sizeof(float));
      }
      if (elem_num >= num_symm)
        num_symm = elem_num + 1;

      /* Read transformation elements */
      tmp_char = buffer[33];
      buffer[33] = 0;
      rot_mats[elem_num*9 + dim_num*3 + 0] = atof(&buffer[24]);
      buffer[33] = tmp_char;
      
      tmp_char = buffer[43];
      buffer[43] = 0;
      rot_mats[elem_num*9 + dim_num*3 + 1] = atof(&buffer[34]);
      buffer[43] = tmp_char;
      
      tmp_char = buffer[53];
      buffer[53] = 0;
      rot_mats[elem_num*9 + dim_num*3 + 2] = atof(&buffer[44]);
      buffer[53] = tmp_char;
      
      tmp_char = buffer[68];
      buffer[68] = 0;
      trans_vecs[elem_num*3 + dim_num] = atof(&buffer[59]);
      buffer[68] = tmp_char;
    }                  /* skip anything else */
  }
  fprintf(stderr, "Read %d atoms with %d electrons\n",res->natoms,total_atomic_number);
  fprintf(stderr, "Read %d symmetry transformations\n", num_symm);
  //printf("Read %d atoms with %d electrons\n",res->natoms,total_atomic_number);
  fclose(fp);
  res->atomic_number = realloc(res->atomic_number,sizeof(int)*res->natoms);
  res->pos = realloc(res->pos,3*sizeof(float)*res->natoms);  
  //printf("%d atoms\n", res->natoms);
  free(rot_mats);
  free(trans_vecs);
  return res;  
}

static void add_symmetry_equivalents(Molecule *mol, float *rot_mats, float *trans_vecs, int num_symm) {
}

/* Return the atomic number of a 2 char symbol or 0 on error */

static int getZfromSymbol(char * symbol)
{
  int z;

  get_legal_atom_names() ;

  if(isalpha(symbol[0])){
    symbol[0] = toupper(symbol[0]);
  }else if(isalpha(symbol[1])){
    symbol[0] = toupper(symbol[1]);
    symbol[1] = ' ';
  }else{
    return 0;
  }
  if(isalpha(symbol[1])){
    symbol[1] = toupper(symbol[1]);
  }else{
    symbol[1] = ' ';
  }
  symbol[2] = 0;
  for (z = 0; z < 106; z++){
    if ((symbol[0] == legal_atom_names[z*2+0]) &&
      (symbol[1] == legal_atom_names[z*2+1])) {
      return z+1;     
    }
  }
  /* If we didn't match any try to match only the first letter */
  symbol[1] = ' ';
  for (z = 0; z < 106; z++){
    if ((symbol[0] == legal_atom_names[z*2+0]) &&
      (symbol[1] == legal_atom_names[z*2+1])) {
      return z+1;     
    }
  }
  /* if we're here print the error */
  fprintf(stderr, "Could not find element for atom named '%c%c'\n", symbol[0], symbol[1]);
  return 0;
}



static void get_legal_atom_names()
{
     /* The full periodic table as 2 character identifiers, left-justified */
     
     strcpy(&legal_atom_names[0], "H HELIBEB C N O F NE");
     strcpy(&(legal_atom_names[20]), "NAMGALSIP S CLARK CA"); 

     strcpy(&legal_atom_names[40], "SCTIV CRMNFECONICUZN") ; 
     strcpy(&legal_atom_names[60], "GAGEASSEBRKRRBSRY ZR") ; 
     strcpy(&legal_atom_names[80], "NBMOTCRURHPDAGCDINSN") ; 
     strcpy(&legal_atom_names[100], "SBTEI XECSBALACEPRND") ; 
     strcpy(&legal_atom_names[120], "PMSMEUGDTBDYHOERTMYB") ; 
     strcpy(&legal_atom_names[140], "LUHFTAW REOSIRPTAUHG") ;
     strcpy(&legal_atom_names[160], "TLPBBIPOATRNFRRAACTH") ; 
     strcpy(&legal_atom_names[180], "PAU NPPUAMCMBKCFESFM") ; 
     strcpy(&legal_atom_names[200], "MDNOLRRFHA") ;

     return;
}


void    write_pdb_from_mol(char *filename,Molecule * mol){
  FILE     *fpout;
  int      i ;
#ifdef MPI
  if(!is_mpi_master()){
    return ;
  }
#endif  
  get_legal_atom_names();
  
  if ((fpout = fopen(filename, "w")) == NULL){
    perror("Cannot open file\n");
    abort();
  }

  for (i = 0; i <  mol->natoms; i++) {
    fprintf(fpout,"ATOM  %5d  %.2s      A   1    %8.3f%8.3f%8.3f\n",
          i%99999,
          &legal_atom_names[2*(mol->atomic_number[i]-1)],
          mol->pos[i*3+0]*1e10,
          mol->pos[i*3+1]*1e10,
          mol->pos[i*3+2]*1e10);
  }
  fprintf(fpout, "END\n") ;
  fclose(fpout);
}

Molecule * alloc_mol() {
  Molecule * mol = malloc(sizeof(Molecule));
  mol->natoms = 0;
  mol->atomic_number = NULL;
  mol->pos = NULL;
  return mol;
}

void add_atom_to_mol(Molecule * mol, int atomic_number, float x, float y, float z) {
  int i = mol->natoms;
  mol->natoms++;
  mol->atomic_number = realloc(mol->atomic_number,sizeof(int)*mol->natoms);
  mol->pos = realloc(mol->pos,sizeof(float)*mol->natoms*3);
  mol->atomic_number[i] = atomic_number;
  mol->pos[i*3+0] = x;
  mol->pos[i*3+1] = y;
  mol->pos[i*3+2] = z;
}

void free_mol(Molecule * mol) {
  free(mol->atomic_number);
  free(mol->pos);
  free(mol);
}

Molecule * make_mol(Image * atomic_number, Image * pos) {
  int i, N;
  Molecule * mol;
  if (pos->image->x != 3) {
    fprintf(stderr,"Array of positions must have 3 times N dimensions, where N denotes the number of atoms.\n");
    return NULL;
  } 
  if (atomic_number->image->x != pos->image->y) {
    fprintf(stderr,"Array of atomic numbers and atomic positions do not have compatible shapes.\n");
    return NULL;
  }
  N = atomic_number->image->x;
  mol = alloc_mol();
  mol->natoms = N;
  mol->atomic_number = realloc(mol->atomic_number,sizeof(int)*N);
  mol->pos = realloc(mol->pos,sizeof(float)*N*3);
  for (i = 0; i < N; i++) {
    mol->atomic_number[i] = (int) atomic_number->image->data[i].re;
    mol->pos[i*3+0] = pos->image->data[i*3+0].re;
    mol->pos[i*3+1] = pos->image->data[i*3+1].re;
    mol->pos[i*3+2] = pos->image->data[i*3+2].re;
  }
  return mol;
}

void origin_to_center_of_mass(Molecule * mol) {
  int i;
  float com[3] = {0., 0., 0.};
  int N = 0;
  // Calculate center of mass (of electron density, not real mass)
  for (i=0; i<mol->natoms; i++) {
    N += mol->atomic_number[i];
    com[0] += mol->pos[i*3+0] * mol->atomic_number[i];
    com[1] += mol->pos[i*3+1] * mol->atomic_number[i];
    com[2] += mol->pos[i*3+2] * mol->atomic_number[i];
  }
  com[0] /= N;
  com[1] /= N;
  com[2] /= N;
  // Now move origin to center of mass
  for (i=0; i<mol->natoms; i++){ 
    mol->pos[i*3+0] -= com[0];
    mol->pos[i*3+1] -= com[1];
    mol->pos[i*3+2] -= com[2];
  }
}
