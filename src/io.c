#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>

#include "prototypes.h"
#include "morton.h"
#include "constant.h"

float assign_grid_field(int field,struct CELL *cell){

/**
  * This function return the appropriate field, depending of a given ID.
  * see also parameters.c for the link between the ID and the input given in param.output.
  */

  float res;


  switch(field){
  case 4:
    res=cell->rank;
    break;
  case 5:
    res=cell->gdata.d;
    break;
  case 6:
    res=cell->cicdens;
    break;
  case 7:
    res=cell->gdata.p;
    break;
  case 8:
    res=cell->gdata.res;
    break;
  case 9:
    res=cell->gdata.f[0];
    break;
  case 10:
    res=cell->gdata.f[1];
    break;
  case 11:
    res=cell->gdata.f[2];
    break;
  case 12:
    res=cell->marked;
    break;

#ifdef WHYDRO2
  case 13:
    res=cell->field.d;
    break;
  case 14:
    res=cell->field.u;
    break;
  case 15:
    res=cell->field.v;
    break;
  case 16:
    res=cell->field.w;
    break;
  case 17:
    res=cell->field.p;
    break;
#endif // WHYDRO2

#ifdef WRAD
  case 18:
    res=cell->rfield.e[0];
    break;
  case 19:
    res=cell->rfield.fx[0];
    break;
  case 20:
    res=cell->rfield.fy[0];
    break;
  case 21:
    res=cell->rfield.fz[0];
    break;
  case 22:
    res=cell->rfield.e[1];
    break;
  case 23:
    res=cell->rfield.fx[1];
    break;
  case 24:
    res=cell->rfield.fy[1];
    break;
  case 25:
    res=cell->rfield.fz[1];
    break;
  case 26:
#ifdef SUPERNOVAE
//    res=cell->rfield.snfb;
#endif // SUPERNOVAE
    break;


  case 27:
    res=cell->rfield.e[2];
    break;
  case 28:
    res=cell->rfield.fx[2];
    break;
  case 29:
    res=cell->rfield.fy[2];
    break;
  case 30:
    res=cell->rfield.fz[2];
    break;

  case 31:
    res=0;
    int igrp;
    for(igrp=0;igrp<NGRP;igrp++) res+=cell->rfield.src[igrp];
    break;
#ifdef WCHEM
  case 32:
    res=cell->rfield.nhplus/cell->rfield.nh;
    break;
#ifdef WRADHYD
  case 33:
    res=cell->field.dX;
    break;
#ifdef HELIUM
  case 34:
    res=cell->field.dXHE;
    break;
  case 35:
    res=cell->field.dXXHE;
    break;
  case 36:
    res=cell->field.xHE;
    break;
  case 37:
    res=cell->field.xxHE;
    break;
#endif // HELIUM
#endif // WRADHYD
  case 38:
    res=cell->rfield.temp;
    break;
#endif // WCHEM
  case 39:
    res=cell->z_first_xion;
    break;
  case 40:
    res=cell->z_last_xion;
    break;
#endif // WRAD
  }
  return res;
}


// =====================================================================

void dumpalloct_serial(char folder[],REAL tsim, struct PARAM *param, struct CPU *cpu, int ldump){


  /**
   * This function dump the output data with STDIO
   * only the most reffined cell are dumped
   *
   * format:
   *   - N_fields time N_procs files corresponding to the defined field in param.output
   *     int n : the number of dumped cells
   *     followed by n times:
   *     float field : value of the corresponding field
   *
   *     TODO: for the flux and all the fields depending of NGRP, dump in 1 field instead of NGRP field
   *
   */

  const int debug=0;

  int i;
  int n_field=0;
  int n_cell=0;
  float xmin,xmax,ymin,ymax,zmin,zmax;

  xmin=2;
  xmax=-1;
  ymin=2;
  ymax=-1;
  zmin=2;
  zmax=-1;

  // Opening all the fields files
  if(debug) printf("Allocating %d file pointers \n", param->out_grid->n_field);

  FILE **f_dat;
  f_dat=(FILE **)malloc(param->out_grid->n_field*sizeof(FILE *));


  for(i=0;i<param->out_grid->n_field_tot;i++){
    if(param->out_grid->field_id[i]){

      char folder_field[128];
      sprintf(folder_field,"%sgrid_%s/",folder, param->out_grid->field_name[i]);
      mkdir(folder_field, 0755);
      char dat_name[256];
      sprintf(dat_name,"%s%s.%05d.p%05d",folder_field,param->out_grid->field_name[i],cpu->ndumps,cpu->rank);

      if(debug) printf("Openning : %s",dat_name);

      f_dat[n_field]=fopen(dat_name,"wb");
      if(f_dat[n_field] == NULL){
	printf("Cannot open %s\n", dat_name);
	abort();
      }

      fwrite(&n_cell,sizeof(int),1,f_dat[n_field]);
      fwrite(&tsim,sizeof(float),1,f_dat[n_field]);
      fwrite(&xmin,sizeof(float),1,f_dat[n_field]);
      fwrite(&xmax,sizeof(float),1,f_dat[n_field]);
      fwrite(&ymin,sizeof(float),1,f_dat[n_field]);
      fwrite(&ymax,sizeof(float),1,f_dat[n_field]);
      fwrite(&zmin,sizeof(float),1,f_dat[n_field]);
      fwrite(&zmax,sizeof(float),1,f_dat[n_field]);

      n_field++;
    }
  }

  if(debug) printf("Files open, let's write");

  // writing the data
  unsigned int level;
  for(level=param->lcoarse;level<=ldump;level++){
    double dx=POW(0.5,level); // oct size

    struct CELL *cell=&(cpu->grid[cpu->firstcell[level]]);
    unsigned long nidx=cpu->ncell[level];
    unsigned long idx;

    for(idx=0;idx<nidx;idx++){

      if((cell->child==0)||(level==ldump)){
	n_cell++;
	
	REAL x[3];
	key2pos(cell->key,x,&level);
	int ii=0;
	for (i =0;i<param->out_grid->n_field_tot; i++){
	  
	  if(param->out_grid->field_id[i]){
	    float dat;
	    if(i==0){
	      dat=x[0];
	    }
	    else if(i==1){
	      dat=x[1];
	    }
	    else if(i==2){
	      dat=x[2];
	    }
	    else if(i==3){
	      dat=level;
	    }
	    else{
	      dat = (float)assign_grid_field(i,cell);
	    }
	    fwrite(&dat,sizeof(float),1,f_dat[ii]);
	    ii++;

	  }
	}
	


	
	// update file boundaries
	if(x[0]<xmin) xmin=x[0];
	if(x[1]<ymin) ymin=x[1];
	if(x[2]<zmin) zmin=x[2];
	
	if(x[0]+dx>xmax) xmax=x[0]+dx;
	if(x[1]+dx>ymax) ymax=x[1]+dx;
	if(x[2]+dx>zmax) zmax=x[2]+dx;
      }

      cell++;
    }
  }


  if(debug) printf("Write OK, header update and close");
  
  // write n_cells and close the fields files
  n_field=0;
  for(i=0;i<param->out_grid->n_field_tot;i++){
    if(param->out_grid->field_id[i]){
      rewind(f_dat[n_field]);
      fwrite(&n_cell,sizeof(int),1,f_dat[n_field]);
      fwrite(&tsim,sizeof(float),1,f_dat[n_field]);
      fwrite(&xmin,sizeof(float),1,f_dat[n_field]);
      fwrite(&xmax,sizeof(float),1,f_dat[n_field]);
      fwrite(&ymin,sizeof(float),1,f_dat[n_field]);
      fwrite(&ymax,sizeof(float),1,f_dat[n_field]);
      fwrite(&zmin,sizeof(float),1,f_dat[n_field]);
      fwrite(&zmax,sizeof(float),1,f_dat[n_field]);
      fclose(f_dat[n_field]);
      n_field++;
    }
  }

  free(f_dat);
}
