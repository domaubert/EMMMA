#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>

#include "prototypes.h"
#include "morton.h"
#include "constant.h"

// =========================================================
void allocation(struct CPU *cpu,struct PARAM *param){
  unsigned long KEYMAX=0;
  unsigned long footprint;
  KEYMAX=~(KEYMAX&0);

  cpu->grid = (struct CELL*)calloc(param->ngridmax,sizeof(struct CELL));
  footprint = param->ngridmax*sizeof(struct CELL)/(1024*1024);
  int i;
  for(i=0;i<param->ngridmax;i++) cpu->grid[i].key=KEYMAX; // INIT TO HUGE KEY


  cpu->firstcell=(unsigned long *)calloc(param->lmax+1,sizeof(unsigned long));
  cpu->ncell=(unsigned long *)calloc(param->lmax+1,sizeof(unsigned long));
  cpu->ncelltotal=0;

#ifdef PIC
  cpu->part=(struct PART*)calloc(param->npartmax,sizeof(struct PART));
  footprint += param->npartmax*sizeof(struct PART)/(1024*1024);
  for(i=0;i<param->npartmax;i++){ 
    cpu->part[i].key=KEYMAX; // INIT TO HUGE KEY
    cpu->part[i].newkey=KEYMAX; // INIT TO HUGE KEY
  }
  cpu->firstpart=(unsigned long *)calloc(param->lmax+1,sizeof(unsigned long));
  cpu->npart=(unsigned long *)calloc(param->lmax+1,sizeof(unsigned long));
  cpu->nparttotal=0;
#endif


  if(cpu->rank==RANK_DISP) printf("allocated %ld MB/proc (total = %ld GB)\n",footprint,footprint*cpu->nproc/1024);

 /* { */
 /*      unsigned long ktest; */
 /*      ktest=2788426; */
 /*      struct CELL *newcell=getcell(&ktest,7,cpu); */
 /*      printf("allocktest= %p\n",newcell); */
 /*    } */

}

#ifdef PIC
// =========================================================
void build_init_part(struct CPU* cpu, struct PARAM *param){
  int i;
  for(i=0;i<param->npartmax;i++){
    unsigned int x,y,z;
    unsigned long key;
    int level=param->lcoarse;
    REAL dx=1./(1<<level);
    if(cpu->part[i].mass){
      x=(int)(cpu->part[i].x[0]/dx);
      y=(int)(cpu->part[i].x[1]/dx);
      z=(int)(cpu->part[i].x[2]/dx);
      LC2M(&key,x,y,z,level);
    
      cpu->part[i].key=key; 
      cpu->part[i].newkey=key; 
      //cpu->part[i].mass=1./param->npartmax;
      cpu->npart[level]++; 
      cpu->nparttotal++; 
    }
  }
}
#endif

// =========================================================
unsigned long code2rank(unsigned long code, int level, struct CPU *cpu){
  unsigned long ntotlevel;
  unsigned long rank;
  ntotlevel=(1<<(3*level));
  
  if(ntotlevel<(8*cpu->nproc)){
    // for very coarse levels we force the attachment to this proc 
    rank=cpu->rank;
  }
  else{
    unsigned long nperrank;
    nperrank=ntotlevel/cpu->nproc;
    rank=code/nperrank;
    //    printf("%ld %ld\n",rank,nperrank);
  }
  return rank;
}


// =========================================================
void build_init_grid(struct CPU *cpu,struct PARAM *param){
  int level;
  int idx=0;
  unsigned long code,code0;
  unsigned long crank=cpu->rank;

  for(level=0;level<=param->lcoarse;level++){
    unsigned int n1=1<<level;
    unsigned int i,j,k;
    unsigned long nlev=0;
    REAL dx=1./(1<<level);
#ifdef WMPI
    LC2M(&code0,0,0,0,level);
#endif
    // cartesian sweep
    for(k=0;k<n1;k++){
      for(j=0;j<n1;j++){
	for(i=0;i<n1;i++){
	  LC2M(&code,i,j,k,level);
#ifdef WMPI
	  // computing local rank
	  crank=code2rank(code-code0,level,cpu);
#endif
	  if(crank==cpu->rank){
	    cpu->grid[idx].key=code;
	    cpu->grid[idx].rank=cpu->rank;
	    if(level<param->lcoarse) cpu->grid[idx].child=1;
	    cpu->grid[idx].gdata.d=0.;
	    cpu->grid[idx].gdata.p=0.;
	  
	    nlev++;
	    idx++;
	  }
	  else{
	    //we scan the 27 neighbors of this cell
	    // if one of them contains the current rank as a direct neighbor, it is being added
	    unsigned long neikey[27]; 
	    unsigned long cranknei; 
	    struct CELL dummy;
	    dummy.key=code;
	    nei27(&dummy,neikey);
	    int inei;
	    for(inei=0;inei<27;inei++){	  
	 
	  cranknei=code2rank(neikey[inei]-code0,level,cpu);
	  if(cranknei==cpu->rank){
	  cpu->grid[idx].key=code;
	  cpu->grid[idx].rank=crank;
	  if(level<param->lcoarse) cpu->grid[idx].child=1; // QUESTION DO WE NEED THIS?
	  cpu->grid[idx].gdata.d=0.;
	  cpu->grid[idx].gdata.p=0.;
	  
	  nlev++;
	  idx++;
	  break;
	}
	}
	}

	}
      }
    }

    cpu->ncell[level]=nlev;
    cpu->ncelltotal+=nlev;
    if(cpu->rank==RANK_DISP){
      printf("%ld on rank %ld for level=%d\n",cpu->ncell[level],cpu->rank,level);
    }
  }
}
