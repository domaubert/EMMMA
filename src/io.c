#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>

#include "prototypes.h"
#include "morton.h"
#include "constant.h"

// =========================================================
float assign_part_field(int field,struct PART *curp){

/**
  * return the appropriate particle field, depending of a given ID.
  */

  float res=0;
  switch(field){
  case 0:
    res=(float)curp->x[0];
    break;
  case 1:
    res=(float)curp->x[1];
    break;
  case 2:
    res=(float)curp->x[2];
    break;

  case 3:
    res=(float)curp->v[0];
    break;
  case 4:
    res=(float)curp->v[1];
    break;
  case 5:
    res=(float)curp->v[2];
    break;

  /* case 6: */
  /*   res=(float)curp->fx; */
  /*   break; */
  /* case 7: */
  /*   res=(float)curp->fy; */
  /*   break; */
  /* case 8: */
  /*   res=(float)curp->fz; */
  /*   break; */
    
  case 9:
    res=(float)(curp->idx);
    break;
    
#ifdef STARS
  case 10:
    res=(float)(curp->isStar);
    break;
#endif // STARS
    
  /* case 11: */
  /*   res=(float)(curp->epot); */
  /*   break; */
  /* case 12: */
  /*   res=(float)(curp->ekin); */
  /*   break; */
    
  case 13:
    res=(float)(curp->mass);
    break;
    
#ifdef STARS
    
  case 14:
    if(curp->isStar) {
      res=(float)(curp->age);
    }
    break;
#endif // STARS
  }
  return res;
}

// =========================================================
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
      sprintf(folder_field,"%s/%05d",folder, cpu->ndumps);
      mkdir(folder_field, 0755);

      sprintf(folder_field,"%s/%05d/grid_%s/",folder, cpu->ndumps,param->out_grid->field_name[i]);
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
	  
	  if(param->out_part->field_id[i]){
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

//void dumpalloct_serial(char folder[],REAL tsim, struct PARAM *param, struct CPU *cpu, int ldump){

void dumppart_serial(char folder[],REAL tsim, struct PARAM *param, struct CPU *cpu, int ldump){
  //(struct RUNPARAMS *param, struct OCT **firstoct,char filename[], int levelcoarse, int levelmax, REAL tsim, struct CPUINFO *cpu){

  const int debug =0;

  FILE *fp = NULL;
  int ipart=0;
  unsigned int level;
  float tsimf=tsim;

  int npart=0;
  int first=-1;
  int n_cell=0;

  float part_xmin=2;
  float part_xmax=-1;
  float part_ymin=2;
  float part_ymax=-1;
  float part_zmin=2;
  float part_zmax=-1;

#ifdef STARS
  int nstar=0;
  float star_xmin=2;
  float star_xmax=-1;
  float star_ymin=2;
  float star_ymax=-1;
  float star_zmin=2;
  float star_zmax=-1;
#endif

  FILE **f_part=(FILE **)malloc(param->out_part->n_field*sizeof(FILE *));

#ifdef STARS
  FILE **f_star=(FILE **)malloc(param->out_part->n_field*sizeof(FILE *));
#endif // STARS


  int n_field=0;
  int i;
  for(i=0;i<param->out_part->n_field_tot;i++){
    if(param->out_part->field_id[i]){

      if(first==-1) first=i; // looking for the first non nil field

#ifdef STARS
      char filenamestar[128];
      char folder_field_star[128];
      sprintf(folder_field_star,"%s/%05d/star_%s/",folder,cpu->ndumps,param->out_part->field_name[i]);
      mkdir(folder_field_star, 0755);
      sprintf(filenamestar,"%s%s.%05d.p%05d",folder_field_star,param->out_part->field_name[i],cpu->ndumps,cpu->rank);


      if(debug) printf("openning %s at %p\n",filenamestar, f_star[n_field]);
      f_star[n_field]=fopen(filenamestar,"wb");
      if(f_star[n_field] == NULL) {
        printf("Cannot open %s\n", filenamestar);
        abort();
      }
      fwrite(&nstar,1,sizeof(int)  ,f_star[n_field]);
      fwrite(&tsimf,1,sizeof(float),f_star[n_field]);
      fwrite(&star_xmin,sizeof(float),1,f_star[n_field]);
      fwrite(&star_xmax,sizeof(float),1,f_star[n_field]);
      fwrite(&star_ymin,sizeof(float),1,f_star[n_field]);
      fwrite(&star_ymax,sizeof(float),1,f_star[n_field]);
      fwrite(&star_zmin,sizeof(float),1,f_star[n_field]);
      fwrite(&star_zmax,sizeof(float),1,f_star[n_field]);
#endif

      char filenamepart[128];
      char folder_field_part[128];
      sprintf(folder_field_part,"%s/%05d/part_%s/",folder,cpu->ndumps,param->out_part->field_name[i]);
      mkdir(folder_field_part, 0755);
      sprintf(filenamepart,"%s%s.%05d.p%05d",folder_field_part,param->out_part->field_name[i],cpu->ndumps,cpu->rank);

      if(debug) printf("openning %s at %p\n",filenamepart,f_part[n_field]);
      f_part[n_field]=fopen(filenamepart,"wb");
      if(f_part[n_field] == NULL){
        printf("Cannot open %s\n", filenamepart);
        abort();
      }

      fwrite(&npart,1,sizeof(int)  ,f_part[n_field]);
      fwrite(&tsimf,1,sizeof(float),f_part[n_field]);
      fwrite(&part_xmin,sizeof(float),1,f_part[n_field]);
      fwrite(&part_xmax,sizeof(float),1,f_part[n_field]);
      fwrite(&part_ymin,sizeof(float),1,f_part[n_field]);
      fwrite(&part_ymax,sizeof(float),1,f_part[n_field]);
      fwrite(&part_zmin,sizeof(float),1,f_part[n_field]);
      fwrite(&part_zmax,sizeof(float),1,f_part[n_field]);

      n_field++;
    }
  }

  if(debug) printf("opening file OK\n");


  // writing the data
  for(level=param->lcoarse;level<=ldump;level++){
    double dx=POW(0.5,level); // oct size

    struct CELL *cell=&(cpu->grid[cpu->firstcell[level]]);
    unsigned long nidx=cpu->ncell[level];
    unsigned long idx;

    for(idx=0;idx<nidx;idx++){

      if((cell->child==0)||(level==ldump)){
	n_cell++;

	struct PART *p=NULL;
	p=getpart(&(cell->key),level,cpu); // returns the first part of the cell
	if(p!=NULL){
	  ipart++;
	  // ========================================================
	  int ii=0;
	  for (i=0;i<param->out_part->n_field_tot; i++){
	    if(param->out_part->field_id[i]){

	      //    if(debug) printf("field_id=%d\n",param->out_part->field_id[i]);

#ifdef STARS
	      if(p->isStar) 	{	
		fp=f_star[ii];	if(i==first) nstar++;	
	      }
	      else{	
		fp=f_part[ii];	if(i==first) npart++;	
	      }
#else
	      fp=f_part[ii];	if(i==first) npart++;	
#endif

	      float dat = (float)assign_part_field(i,p);
	      fwrite(&dat,sizeof(float),1,fp);
	      ii++;
	    }
	  }
	  // ========================================================

	  
	  while(p<(cpu->part+cpu->nparttotal-1)){
	    p++;
	    if(p->key!=cell->key){
	      break; // end of current cell particle stream
	    }
	    else{
	      ipart++;
	      // ========================================================
	      int ii=0;
	      for (i=0;i<param->out_part->n_field_tot; i++){
		if(param->out_part->field_id[i]){
		  
		  //    if(debug) printf("field_id=%d\n",param->out_part->field_id[i]);
		  
#ifdef STARS
		  if(p->isStar) 	{	
		    fp=f_star[ii];	if(i==first) nstar++;	
		  }
		  else{	
		    fp=f_part[ii];	if(i==first) npart++;	
		  }
#else
		    fp=f_part[ii];	if(i==first) npart++;	
#endif		  
		  float dat = (float)assign_part_field(i,p);
		  fwrite(&dat,sizeof(float),1,fp);
		  ii++;
		}
	      }
	      // ========================================================
	    }
	  }
	}
	
	REAL x[3];
	key2pos(cell->key,x,&level);

	// update file boundaries
#ifdef STARS
	if(x[0]<star_xmin) star_xmin=x[0];
	if(x[1]<star_ymin) star_ymin=x[1];
	if(x[2]<star_zmin) star_zmin=x[2];
	
	if(x[0]+dx>star_xmax) star_xmax=x[0]+dx;
	if(x[1]+dx>star_ymax) star_ymax=x[1]+dx;
	if(x[2]+dx>star_zmax) star_zmax=x[2]+dx;
#endif

	if(x[0]<part_xmin) part_xmin=x[0];
	if(x[1]<part_ymin) part_ymin=x[1];
	if(x[2]<part_zmin) part_zmin=x[2];
	
	if(x[0]+dx>part_xmax) part_xmax=x[0]+dx;
	if(x[1]+dx>part_ymax) part_ymax=x[1]+dx;
	if(x[2]+dx>part_zmax) part_zmax=x[2]+dx;

      }

      cell++;
    }
    printf("ipart=%d\n",ipart);
  }

  if (debug) printf("writing OK\n");

  n_field=0;
  for(i=0;i<param->out_part->n_field_tot;i++){
    if(param->out_part->field_id[i]){

      rewind(f_part[n_field]);
      fwrite(&npart,1,sizeof(int)  ,f_part[n_field]);
      fwrite(&tsimf,1,sizeof(float),f_part[n_field]);
      fwrite(&part_xmin,sizeof(float),1,f_part[n_field]);
      fwrite(&part_xmax,sizeof(float),1,f_part[n_field]);
      fwrite(&part_ymin,sizeof(float),1,f_part[n_field]);
      fwrite(&part_ymax,sizeof(float),1,f_part[n_field]);
      fwrite(&part_zmin,sizeof(float),1,f_part[n_field]);
      fwrite(&part_zmax,sizeof(float),1,f_part[n_field]);
      fclose(f_part[n_field]);

#ifdef STARS
      rewind(f_star[n_field]);
      fwrite(&nstar,1,sizeof(int)  ,f_star[n_field]);
      fwrite(&tsimf,1,sizeof(float),f_star[n_field]);
      fwrite(&star_xmin,sizeof(float),1,f_star[n_field]);
      fwrite(&star_xmax,sizeof(float),1,f_star[n_field]);
      fwrite(&star_ymin,sizeof(float),1,f_star[n_field]);
      fwrite(&star_ymax,sizeof(float),1,f_star[n_field]);
      fwrite(&star_zmin,sizeof(float),1,f_star[n_field]);
      fwrite(&star_zmax,sizeof(float),1,f_star[n_field]);
      fclose(f_star[n_field]);
#endif
      n_field++;
    }
  }

  if (debug) printf("closing  OK\n");
  printf("wrote %d particles (%d expected)\n",ipart,npart);

  free(f_part);
#ifdef STARS
  free(f_star);
#endif
}
