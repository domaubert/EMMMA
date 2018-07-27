#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "constant.h"
#include "prototypes.h"
#include "morton.h"

// CIC Assign ===================================
void part2cell(struct PART *part, struct CELL *cell, struct CPU *cpu, unsigned long *neikey, int level){
  REAL dx=1./(1<<level);
  REAL dpart=part->mass/(dx*dx*dx);
  int inei;
  for(inei=0;inei<27;inei++){ // Brute force, we look for the 27 neighbors
    struct CELL *lcell;
    lcell=getcell(neikey+inei,level,cpu);
    if(lcell){
      unsigned int lloc;
      REAL xc[3];
      key2cen(neikey[inei],xc,&lloc);
      
      REAL w=1.;
      int i;
      for(i=0;i<3;i++){
	REAL dpc=part->x[i]-xc[i];
	dpc=(dpc>0.5?1.-dpc:(dpc<-0.5?1.+dpc:dpc)); // dealing with boundary conditions

	REAL d=dpc/dx;
	w*=(d>=1.?0.:(d<=-1.?0.:(d<0.?1.+d:1.-d)));
      }
      lcell->cicdens+=w*dpart;
    }
  }

}



// CIC Assign ===================================
void cell2part(struct PART *part, struct CELL *cell, struct CPU *cpu, unsigned long *neikey, int level, struct CELL *mumcell, unsigned long *neikeymum, REAL dt){

  REAL dx=1./(1<<level);
  REAL accel[3]={0.,0.,0.};
  int ic; // for directions
  int lowres=0;
  int inei;
  for(inei=0;inei<27;inei++){ // Brute force, we look for the 27 neighbors
    struct CELL *lcell;
    lcell=getcell(neikey+inei,level,cpu);
    if(lcell){
      unsigned int lloc;
      REAL xc[3];
      key2cen(neikey[inei],xc,&lloc);
      
      REAL w=1.;
      int i;
      for(i=0;i<3;i++){
	REAL d=(part->x[i]-xc[i])/dx;
	w*=(d>=1.?0.:(d<=-1.?0.:(d<0.?1.+d:1.-d)));
      }
      for(ic=0;ic<3;ic++) accel[ic]+=lcell->gdata.f[ic]*w;
    }
    else{
      lowres=1;
      break;
    }
  }

  // boundary particle, let's start all over at lower res
  if(lowres){
    for(ic=0;ic<3;ic++) accel[ic]=0.;
    dx=1./(1<<(level-1));
    int inei;
    for(inei=0;inei<27;inei++){ // Brute force, we look for the 27 neighbors
      struct CELL *lcell;
      lcell=getcell(neikeymum+inei,level-1,cpu);
      if(lcell){
	unsigned int lloc;
	REAL xc[3];
	key2cen(neikeymum[inei],xc,&lloc);
	
	REAL w=1.;
	int i;
	for(i=0;i<3;i++){
	  REAL d=(part->x[i]-xc[i])/dx;
	  w*=(d>=1.?0.:(d<=-1.?0.:(d<0.?1.+d:1.-d)));
	}
	for(ic=0;ic<3;ic++) accel[ic]+=lcell->gdata.f[ic]*w;
      }
      else{
	printf("ISSUE with force boundary... abort\n");
	abort();
      }
    }
  }


  // at this stage we have the boundaries, let's accelerate the particle

  for(ic=0;ic<3;ic++) part->v[ic]+=-accel[ic]*dt;

}


//==================================================================

void part2cell_13(struct PART *part, struct CELL *cell, struct CPU *cpu,  unsigned long *neikey, int level){
  REAL dx=1./(1<<level);
  REAL dpart=part->mass/(dx*dx*dx);
  int inei;
  for(inei=13;inei<14;inei++){ // Brute force, we look for the 27 neighbors
    struct CELL *lcell;
    lcell=getcell(neikey+inei,level,cpu);
    if(lcell){
      unsigned int lloc;
      REAL xc[3];
      key2cen(neikey[inei],xc,&lloc);
      
      REAL w=1.;
      int i;
      for(i=0;i<3;i++){
	REAL d=(part->x[i]-xc[i])/dx;
	w*=(d>=1.?0.:(d<=-1.?0.:(d<0.?1.+d:1.-d)));
      }
      lcell->cicdens+=w*dpart;
    }
  }
}

//==========================================
//==========================================

void cic_child(struct CELL *cellorg, unsigned int levelorg,  unsigned long *neikeyorg, struct CELL *cellpart, unsigned int level, struct CPU *cpu){

  // here cellorg is the cell that should be CIC assigned
  // here cellpart is the mother cell that is being explored for particles in her leaves
  // cellpart level is given by level

  unsigned long icell;

  for(icell=0;icell<8;icell++){
    unsigned long newkey;
    newkey=(cellpart->key<<3)|icell; // getting the daughter's key
    struct CELL *newcell;
    newcell=getcell(&newkey,level+1,cpu);

    if(!newcell->child){
      // it's a leaf lets look for particles

      // === PARTICLE STREAM
      struct PART *p=NULL;
      p=getpart(&(newcell->key),level+1,cpu); // returns the first part of the cell
      if(p!=NULL){
	part2cell(p,cellorg,cpu,neikeyorg,levelorg);
	while(p<(cpu->part+cpu->nparttotal-1)){
	  p++;
	  if(p->key!=newcell->key){
	    break; // end of current cell particle stream
	  }
	  else{
	    part2cell(p,cellorg,cpu,neikeyorg,levelorg);
	  }
	}
      }
      // === END PARTICLE STREAM
    }
    else{
      cic_child(cellorg,levelorg,neikeyorg,newcell,level+1,cpu);
    }
  }
}

//=============================================================
void cic(unsigned int level, struct CPU *cpu, struct PARAM *param){
  
  struct CELL *cell=&(cpu->grid[cpu->firstcell[level]]);
  unsigned long nidx=cpu->ncell[level];
  unsigned long idx;
  unsigned long neikey[27];

  // first cleaning the density et set them to zero
  for(idx=0;idx<nidx;idx++){
    cell->cicdens=0.;
    cell++;
  }


  // let's go
  cell=&(cpu->grid[cpu->firstcell[level]]); //reset to first cell

  for(idx=0;idx<nidx;idx++){
    nei27(cell,neikey);
    // FIRST THE EASY CASE : contribution to level L cells
    if(cell->child){
      // the cell is split let's look down
      cic_child(cell,level,neikey,cell,level,cpu);
    }
    else{
      struct PART *p=NULL;
      p=getpart(&(cell->key),level,cpu); // returns the first part of the cell
      if(p!=NULL){
	part2cell(p,cell,cpu,neikey,level);
	while(p<(cpu->part+cpu->nparttotal-1)){
	  p++;
	  if(p->key!=cell->key){
	    break; // end of current cell particle stream
	  }
	  else{
	    part2cell(p,cell,cpu,neikey,level);
	  }
	}
      }
    }

    if(level>param->lcoarse){
      // SECOND let's exlore L-1 neighbors
      unsigned long mumkey;
      mumkey=(cell->key>>3);
      struct CELL *mumcell;
      mumcell=getcell(&mumkey,level-1,cpu);
      unsigned long neikeymum[27];
      nei27(mumcell,neikeymum);
      int inei;
      for(inei=0;inei<27;inei++){
	struct CELL *lcell;
	lcell=getcell(neikeymum+inei,level-1,cpu);
	if(lcell){
	  if(!lcell->child){
	    struct PART *p=NULL;
	    p=getpart(&(lcell->key),level-1,cpu); // returns the first part of the cell
	    if(p!=NULL){
	      part2cell_13(p,cell,cpu,neikey,level);
	      while(p<(cpu->part+cpu->nparttotal-1)){
		p++;
		if(p->key!=lcell->key){
		  break; // end of current cell particle stream
		}
		else{
		  part2cell_13(p,cell,cpu,neikey,level);
		}
	      }
	    }
	  }
	}
      }
    }


    cell++;
  }
}


