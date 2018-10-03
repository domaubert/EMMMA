
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "constant.h"
#include "prototypes.h"
#include "morton.h"
#include "cic.h"


int countpart(struct CELL *cell, struct CPU *cpu, struct PARAM *param, unsigned int level){
  unsigned int npart=0;

  if(!cell->child){
    struct PART *p=NULL;
    p=getpart(&(cell->key),level,cpu); // returns the first part of the cell
    
    if(p!=NULL){
      npart++;
      while(p<(cpu->part+cpu->nparttotal-1)){
	p++;
	if(p->key!=cell->key){
	  break; // end of current cell particle stream
	}
	else{
	  npart++;
	}
      }
    }
  }
  else{
    unsigned long icell;
    unsigned long root=cell->key<<3; // root key for all the childs
    for(icell=0;icell<8;icell++){
      unsigned long key=(root|icell); //key for a child
      struct CELL *lcell;
      lcell=getcell(&key,level+1,cpu);
      npart+=countpart(lcell,cpu,param,level+1);
    }
  }
  return npart;
}

//=============================================================
void L_accelpart(unsigned int level, struct CPU *cpu, struct PARAM *param,int is){
  
  struct CELL *cell=&(cpu->grid[cpu->firstcell[level]]);
  struct CELL *mumcell;
  unsigned long nidx=cpu->ncell[level];
  unsigned long idx;
  unsigned long neikey[27];
  unsigned long neikeymum[27];
  unsigned long keymum;

  for(idx=0;idx<nidx;idx++){
    nei27(cell,neikey);

    // mum treatment in case of particles at boundary
    keymum=(cell->key>>3);
    mumcell=getcell(&keymum,level-1,cpu);
    nei27(mumcell,neikeymum);

    // FIRST THE EASY CASE : contribution to level L cells
    if(!cell->child){
      struct PART *p=NULL;
      p=getpart(&(cell->key),level,cpu); // returns the first part of the cell
      if(p!=NULL){
	
	if((p->level==level)||(p->is==0)||(is==-1)){
	  cell2part(p,cell,cpu,neikey,level,mumcell,neikeymum,cpu->adt[level]); 
	}
	while(p<(cpu->part+cpu->nparttotal-1)){
	  p++;
	  if(p->key!=cell->key){
	    break; // end of current cell particle stream
	  }
	  else{
	    if((p->level==level)||(p->is==0)||(is==-1)){
	      cell2part(p,cell,cpu,neikey,level,mumcell,neikeymum, cpu->adt[level]); //here 0.5 because of half timestep for midpoint rule
	    }
	  }
	}
      }
    }
    cell++;
  }
}

//==========================================================================
//==========================================================================
//==========================================================================

REAL L_comptstep(unsigned int level, struct CPU *cpu, struct PARAM *param){
  
  struct CELL *cell=&(cpu->grid[cpu->firstcell[level]]);
  unsigned long nidx=cpu->ncell[level];
  unsigned long idx;
  int ic;
  REAL va;
  REAL dtnew=1e9;
  REAL dxcur=1./(1<<level);
  REAL dtlev;

  for(idx=0;idx<nidx;idx++){
    // FIRST THE EASY CASE : contribution to level L cells
    if(!cell->child){
      struct PART *p=NULL;
      p=getpart(&(cell->key),level,cpu); // returns the first part of the cell
      if(p!=NULL){
	
	// === do stuff
	va=0.;
	for(ic=0;ic<3;ic++){
	  va+=p->v[ic]*p->v[ic];
	}
	
	if(va>0){
	  dtlev=(FRACDX*dxcur)/SQRT(va);
	}
	else{
	  dtlev=1e9;
	}

	if(dtnew>dtlev) dtnew=dtlev;
	// === END do stuff


	while(p<(cpu->part+cpu->nparttotal-1)){
	  p++;
	  if(p->key!=cell->key){
	    break; // end of current cell particle stream
	  }
	  else{
	    // === do stuff
	    va=0.;
	    for(ic=0;ic<3;ic++){
	      va+=p->v[ic]*p->v[ic];
	    }
	    
	    if(va>0){
	      dtlev=(FRACDX*dxcur)/SQRT(va);
	    }
	    else{
	      dtlev=1e9;
	    }
	    
	    if(dtnew>dtlev) dtnew=dtlev;
	    // === END do stuff
	  }
	}
      }
    }
    cell++;
  }
  return dtnew;
}

//=============================================================================================
//=============================================================================================

void L_movepart(unsigned int level, struct CPU *cpu, struct PARAM *param, int is){
  
  struct CELL *cell=&(cpu->grid[cpu->firstcell[level]]);
  unsigned long nidx=cpu->ncell[level];
  unsigned long idx;
  REAL dt=cpu->adt[level];
  int ic;
  REAL mdisp=0.,disp;
  REAL dxcur=1./(1<<level);

  for(idx=0;idx<nidx;idx++){
    // FIRST THE EASY CASE : contribution to level L cells
    if(!cell->child){
      struct PART *p=NULL;
      p=getpart(&(cell->key),level,cpu); // returns the first part of the cell
      if(p!=NULL){
	if((p->level==level)||(p->is==0)){
	  disp=0.;
	  for(ic=0;ic<3;ic++){
	    p->x[ic]+=p->v[ic]*dt;
	    disp+=p->v[ic]*dt*p->v[ic]*dt;

	    // PERIODIC BOUNDARY CONDITIONS
	    p->x[ic]+=(p->x[ic]<0?1.0:(p->x[ic]>=1.?-1.:0));
	    if(p->x[ic]==1.0) p->x[ic]=0.;
	  }

	  // New Key
	  p->newkey=pos2key(p->x,level);
	  

	  if(disp>mdisp) mdisp=disp;
	}
	
	while(p<(cpu->part+cpu->nparttotal-1)){
	  p++;
	  if(p->key!=cell->key){
	    break; // end of current cell particle stream
	  }
	  else{
	    if((p->level==level)||(p->is==0)){
	      disp=0.;
	      for(ic=0;ic<3;ic++){
		p->x[ic]+=p->v[ic]*dt;
		disp+=p->v[ic]*dt*p->v[ic]*dt;

		// PERIODIC BOUNDARY CONDITIONS
		p->x[ic]+=(p->x[ic]<0?1.0:(p->x[ic]>=1.?-1.0:0));
		if(p->x[ic]==1.0) p->x[ic]=0.;
	      }

	      // New Key
	      p->newkey=pos2key(p->x,level);
	      
	      if(disp>mdisp) mdisp=disp;
	    }
	  }
	}
      }
    }
    cell++;
  }

  REAL mmdisp;
  mmdisp=mdisp;
/* #ifdef WMPI */
/*   MPI_Allreduce(&mdisp,&mmdisp,1,MPI_REEL,MPI_MAX,cpu->comm); */
/*   //  mdisp=mmdisp; */
/* #else */
  /* mmdisp=mdisp; */
/* #endif */

  if(cpu->rank==RANK_DISP) printf("level=%d maxdisp=%e or %e dx\n",level,SQRT(mmdisp),SQRT(mmdisp)/dxcur);

}

// ===============================================================================================================
// ===============================================================================================================

void L_levpart(unsigned int level, struct CPU *cpu, struct PARAM *param,int is){
  
  struct CELL *cell=&(cpu->grid[cpu->firstcell[level]]);
  unsigned long nidx=cpu->ncell[level];
  unsigned long idx;

  for(idx=0;idx<nidx;idx++){
    if(!cell->child){
      struct PART *p=NULL;
      p=getpart(&(cell->key),level,cpu); // returns the first part of the cell
      if(p!=NULL){

	p->level=level;
	p->is=is;

	while(p<(cpu->part+cpu->nparttotal-1)){
	  p++;
	  if(p->key!=cell->key){
	    break; // end of current cell particle stream
	  }
	  else{
	    p->level=level;
	    p->is=is;
	  }
	}
      }
    }
    cell++;
  }
}

// ================================================================================
// ================================================================================
void L_resetispart(unsigned int level, struct CPU *cpu, struct PARAM *param){
  
  struct CELL *cell=&(cpu->grid[cpu->firstcell[level]]);
  unsigned long nidx=cpu->ncell[level];
  unsigned long idx;

  for(idx=0;idx<nidx;idx++){
    if(!cell->child){
      struct PART *p=NULL;
      p=getpart(&(cell->key),level,cpu); // returns the first part of the cell
      if(p!=NULL){

	p->is=0;

	while(p<(cpu->part+cpu->nparttotal-1)){
	  p++;
	  if(p->key!=cell->key){
	    break; // end of current cell particle stream
	  }
	  else{
	    p->is=0;
	  }
	}
      }
    }
    cell++;
  }
}
