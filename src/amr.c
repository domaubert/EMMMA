#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <string.h>

#include "prototypes.h"
#include "constant.h"
#include "morton.h"
#include "particle.h"
//=====================================================

// SPHERICAL INNER REGION

/* int crit_phy(struct CELL *cell, struct CPU *cpu, struct PARAM *param, unsigned int level){ */
/*   unsigned int levelloc; */
/*   REAL x[3]; */

/*   key2cen(cell->key,x,&levelloc); */
/*   if(POW(x[0]-0.5,2)+POW(x[1]-0.5,2)+POW(x[2]-0.5,2)<POW(0.15,2)) return 1; */
/*   return 0; */
/* } */


// LAGRANGIAN CRITERION =================================

int crit_phy(struct CELL *cell, struct CPU *cpu, struct PARAM *param, unsigned int level){
  int npart=countpart(cell,cpu,param,level);
  if(npart>param->amrthresh0) return 1;
  return 0;
}


//=====================================================
void mark_phy(unsigned int level, struct CPU *cpu, struct PARAM *param, int ismooth){

  struct CELL *cell=&(cpu->grid[cpu->firstcell[level]]);
  unsigned long nidx=cpu->ncell[level];
  unsigned long idx;
  int nm=0;
  int marker=3+3*ismooth;

  for(idx=0;idx<nidx;idx++){
    if((cell->marked==0)&&(crit_phy(cell,cpu,param,level))){
      cell->marked=marker;
      nm++;
    }
    cell++;
  }

  printf("%d/%ld PHY marked  on level %d\n",nm,nidx,level);
}


//=====================================================
void mark_child(unsigned int level, struct CPU *cpu, struct PARAM *param, int ismooth){

  struct CELL *cell=&(cpu->grid[cpu->firstcell[level]]);
  unsigned long nidx=cpu->ncell[level];
  unsigned long idx;
  int nm=0;
  int marker=1+3*ismooth;

  for(idx=0;idx<nidx;idx++){
    if((cell->marked==0)&&(cell->child)){
      unsigned long icell;
      unsigned long root=cell->key<<3; // root key for all the childs
      int smark=0;
      for(icell=0;icell<8;icell++){
	unsigned long key=(root|icell); //key for a child
	struct CELL *lcell;
	lcell=getcell(&key,level+1,cpu);
	smark+=((lcell->child)||(lcell->marked));
      }
      if(smark){
	cell->marked=marker;
	nm++;
      }
    }
    cell++;
  }

  printf("%d/%ld CHI marked  on level %d\n",nm,nidx,level);
}


//=====================================================
void mark_nei(unsigned int level, struct CPU *cpu, struct PARAM *param, int ismooth){

  struct CELL *cell=&(cpu->grid[cpu->firstcell[level]]);
  unsigned long nidx=cpu->ncell[level];
  unsigned long idx;
  int nm=0;
  int marker=2+3*ismooth;

  for(idx=0;idx<nidx;idx++){
    if((cell->marked>0)&&(cell->marked<marker)){
      unsigned long neikey[27];
      nei27(cell,neikey);
      int inei;
      for(inei=0;inei<27;inei++){
	struct CELL *lcell;
	if(inei==13) continue; //skip central cell
	lcell=getcell(neikey+inei,level,cpu);
	if((lcell!=NULL)&&(lcell->marked==0)){
	  lcell->marked=marker;
	  nm++;
      	}
      }
    }
    cell++;
  }

  printf("%d/%ld NEI marked  on level %d\n",nm,nidx,level);
}


//=============================================================
void clean_mark(unsigned int level, struct CPU *cpu, struct PARAM *param){
  
  struct CELL *cell=&(cpu->grid[cpu->firstcell[level]]);
  unsigned long nidx=cpu->ncell[level];
  unsigned long idx;

  // first cleaning the density et set them to zero
  for(idx=0;idx<nidx;idx++){
    cell->marked=0.;
    cell++;
  }
}


//=====================================================
void create_cell(unsigned int level, struct CPU *cpu, struct PARAM *param){

  struct CELL *cell=&(cpu->grid[cpu->firstcell[level]]);
  unsigned long nidx=cpu->ncell[level];
  unsigned long idx;
  int np=0;

  for(idx=0;idx<nidx;idx++){

    // let's create new cells
    if((cell->marked)&&(!cell->child)){
      cell->child=1; // the current cell has a child
      unsigned long icell;
      for(icell=0;icell<8;icell++){
	struct CELL *newcell=&cpu->grid[cpu->ncelltotal];
	newcell->key=(cell->key<<3)|icell;
	newcell->child=0;
	newcell->marked=0;

	// Straight injection
	newcell->gdata.d=cell->gdata.d;
	newcell->gdata.p=cell->gdata.p; // for testing purpose
	
	// global update
	cpu->ncell[level+1]++;
	cpu->ncelltotal++;

      }

#ifdef PIC
      // === PARTICLE STREAM
      struct PART *p=NULL;
      p=getpart(&(cell->key),level,cpu); // returns the first part of the cell

      if(p!=NULL){
	REAL dx=1./(1<<(level+1));

	// we go brute force to assign a new key to the particles
	unsigned int x=(unsigned int)(p->x[0]/dx);
	unsigned int y=(unsigned int)(p->x[1]/dx);
	unsigned int z=(unsigned int)(p->x[2]/dx);
	unsigned long key;
	LC2M(&key,x,y,z,level+1);
	p->newkey=key;
	p->level+=1;
	np++;

	while(p<(cpu->part+cpu->nparttotal-1)){
	  p++;
	  if(p->key!=cell->key){
	    break; // end of current cell particle stream
	  }
	  else{
	    unsigned int x=(unsigned int)(p->x[0]/dx);
	    unsigned int y=(unsigned int)(p->x[1]/dx);
	    unsigned int z=(unsigned int)(p->x[2]/dx);
	    unsigned long key;
	    LC2M(&key,x,y,z,level+1);
	    p->newkey=key;
	    p->level+=1;
	    np++;
	  }
	}
      }
      // === END PARTICLE STREAM
#endif // PIC

    }
    cell++;
  }
#ifdef PIC
  printf("%d particles to be refined\n",np);
#endif
}

//===================================================================================

void destroy_cell(unsigned int level, struct CPU *cpu, struct PARAM *param){

  struct CELL *cell=&(cpu->grid[cpu->firstcell[level]]);
  unsigned long nidx=cpu->ncell[level];
  unsigned long idx;
  int np=0;
  unsigned long KEYMAX=0;
  KEYMAX=~(KEYMAX&0);

  for(idx=0;idx<nidx;idx++){
    // let's destroy cells 
    if((cell->marked==0)&&(cell->child)){
      
      //flag destroying child cells
      unsigned long icell;
      for(icell=0;icell<8;icell++){
	unsigned long key;
	struct CELL *newcell=&cpu->grid[cpu->ncelltotal];
	key=(cell->key<<3)|icell;
	newcell=getcell(&key,level+1,cpu);
	newcell->desflag=1;

#ifdef PIC
	// === PARTICLE STREAM
	struct PART *p=NULL;
	p=getpart(&(newcell->key),level+1,cpu); // returns the first part of the cell
	if(p!=NULL){
	  p->newkey=cell->key;
	  p->level-=1;
	  np++;

	  while(p<(cpu->part+cpu->nparttotal-1)){
	    p++;
	    if(p->key!=newcell->key){
	      break; // end of current cell particle stream
	    }
	    else{
	      p->newkey=cell->key;
	      p->level-=1;
	      np++;
	    }
	  }
	}
      // ==== END PARTICLE STREAM
#endif // PIC

      }

      cell->child=0; // the current cell has no child
    }
    cell++;
  }

#ifdef PIC
  printf("%d particles to be derefined\n",np);
#endif

  // once done, we rescan all the level+1 cell to nullify them and change their keys
  cell=&(cpu->grid[cpu->firstcell[level+1]]);
  nidx=cpu->ncell[level+1];
  for(idx=0;idx<nidx;idx++){
    if(cell->desflag){
	memset(cell,0,sizeof(struct CELL));
	cell->key=KEYMAX;

	// global update
	cpu->ncell[level+1]--;
	cpu->ncelltotal--;
    }
    cell++;
  }

}
