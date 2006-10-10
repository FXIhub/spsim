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
#include "config.h"
#include "molecule.h"


static char legal_atom_names[106][2] ;
static void get_legal_atom_names();
static int getZfromSymbol(char * symbol);

Molecule * get_Molecule_from_formula(Chem_Formula * form, float * box){
  int i,j,d;
  Molecule * res = malloc(sizeof(Molecule));
  res->atomic_number = NULL;
  res->pos = NULL;
  res->natoms = 0;
  for(d = 0;d<3;d++){
    res->box[d] = box[d];
  }
  for(i = 0;i<form->n_elements;i++){
    res->natoms += form->quantity[i];
    res->atomic_number = realloc(res->atomic_number,sizeof(int)*res->natoms);
    res->pos = realloc(res->pos,sizeof(float)*res->natoms*3);
    for(j = res->natoms-form->quantity[i];j<res->natoms;j++){
      res->atomic_number[j] = form->atomic_number[i];
      for(d = 0;d<3;d++){
	res->pos[3*j+d] = drand48();
      }
    }
  }
  return res;  
}

Molecule * get_Molecule_from_pdb(char * filename){
  FILE * fp = fopen(filename,"r");
  char buffer[1024];
  int maxnatoms = 1024;
  char	gid[10], aid[10] ,atid[3];
  char       tmp_char;
  char       *next_field_start;
  float t1,t2,t3,t4,t5;
  Molecule * res = malloc(sizeof(Molecule));
  if(!fp){
    perror("Error reading PDB\n");
    abort();
  }
  res->atomic_number = malloc(sizeof(int)*maxnatoms);
  res->pos = malloc(3*sizeof(float)*maxnatoms);
  res->natoms = 0;

  get_legal_atom_names();

  while (fgets(buffer, 1024, fp) != NULL)  {
    if (((int)strlen(buffer) > 1) &&
	((strstr(buffer, "ATOM") == buffer) ||
	 (strstr(buffer, "HETATM") == buffer))) {
      
      if (res->natoms >= maxnatoms) {
	/* increase array sizes */
	maxnatoms *= 2;
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

      tmp_char = buffer[11];		/* limit the field*/
      buffer[11] = 0;
      next_field_start = &(buffer[6]);
      /*	       res->pdb_atom_number[res->natoms] = atoi(next_field_start);*/
      buffer[11] = tmp_char;		/* advance field pointer */
      next_field_start+= 6;
      tmp_char = buffer[16];
      buffer[16] = 0;
      
      /* Next copy Remoteness indicator and Branch 
	 designator besides the Chemical symbol */
      strncpy(aid,next_field_start,5);	       
      aid[4] = 0;
      buffer[16] = tmp_char;
      
      strncpy(atid, aid, 2) ; 
      if(isspace(atid[0]) || isdigit(atid[0])){
	atid[0] = atid[1];		/* left justify */
	atid[1] = ' ';
      }
      res->atomic_number[res->natoms] = getZfromSymbol(atid);
      
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
      /*	       res->pdb_groupno[res->natoms] = atoi(next_field_start);*/
      buffer[26] = tmp_char;
      
      /* Start retrieving coordinates*/
      next_field_start += 8;
      tmp_char = buffer[38];
      buffer[38] = 0;
      t1 = (float)atof(next_field_start);
      buffer[38] = tmp_char;
      next_field_start += 8;
      tmp_char = buffer[38];
      buffer[46] = 0;
      t2 = (float)atof(next_field_start);
      buffer[46] = tmp_char;
      next_field_start += 8;
      tmp_char = buffer[54];
      buffer[54] = 0;
      t3 = (float)atof(next_field_start);
      
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
      /*	       if(isspace(pdb->atid[pdb->Nin][0]) ||
		       isdigit(pdb->atid[pdb->Nin][0])){*/
      /* left justify */
      /*		 pdb->atid[pdb->Nin][0] = pdb->atid[pdb->Nin][1];	
			 pdb->atid[pdb->Nin][1] = ' ';
			 }
			 strncpy(pdb->groupid[pdb->Nin], gid, 3) ; 
			 pdb->groupid[pdb->Nin][3] = 0;*/

      /* convert to meters */
      res->pos[res->natoms*3] = t1*1e-10 ;
      res->pos[res->natoms*3+1] = t2*1e-10 ;
      res->pos[res->natoms*3+2] = t3*1e-10 ;
      res->natoms++;
    }			/* skip anything else */
  }
  fclose(fp);
  res->atomic_number = realloc(res->atomic_number,sizeof(int)*res->natoms);
  res->pos = realloc(res->pos,3*sizeof(float)*res->natoms);  
  return res;  
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
  if(isalpha(symbol[1]))
    symbol[1] = toupper(symbol[1]);
  else
    symbol[1] = ' ';
  symbol[2] = 0;
  for (z = 0; z < 106; z++) 
    if ((symbol[0] == legal_atom_names[z][0]) &&
	(symbol[1] == legal_atom_names[z][1])) 
      return z+1;     
  return 0;
}



static void get_legal_atom_names()
{
     /* The full periodic table as 2 character identifiers, left-justified */
     
     strcpy(legal_atom_names [0], "H HELIBEB C N O F NE") ;
     strcpy(legal_atom_names[10], "NAMGALSIP S CLARK CA") ; 
     strcpy(legal_atom_names[20], "SCTIV CRMNFECONICUZN") ; 
     strcpy(legal_atom_names[30], "GAGEASSEBRKRRBSRY ZR") ; 
     strcpy(legal_atom_names[40], "NBMOTCRURHPDAGCDINSN") ; 
     strcpy(legal_atom_names[50], "SBTEI XECSBALACEPRND") ; 
     strcpy(legal_atom_names[60], "PMSMEUGDTBDYHOERTMYB") ; 
     strcpy(legal_atom_names[70], "LUHFTAW REOSIRPTAUHG") ; 
     strcpy(legal_atom_names[80], "TLPBBIPOATRNFRRAACTH") ; 
     strcpy(legal_atom_names[90], "PAU NPPUAMCMBKCFESFM") ; 
     strcpy(legal_atom_names[100], "MDNOLRRFHA") ;

     return;
}
