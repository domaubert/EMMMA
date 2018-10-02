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
  KEYMAX=~(KEYMAX&0);

  cpu->grid=(struct CELL*)calloc(param->ngridmax,sizeof(struct CELL));
  printf("allocated %ld MB\n",param->ngridmax*sizeof(struct CELL)/(1024*1024));
  int i;
  for(i=0;i<param->ngridmax;i++) cpu->grid[i].key=KEYMAX; // INIT TO HUGE KEY


  cpu->firstcell=(unsigned long *)calloc(param->lmax+1,sizeof(unsigned long));
  cpu->ncell=(unsigned long *)calloc(param->lmax+1,sizeof(unsigned long));
  cpu->ncelltotal=0;

#ifdef PIC
  cpu->part=(struct PART*)calloc(param->npartmax,sizeof(struct PART));
  for(i=0;i<param->npartmax;i++){ 
    cpu->part[i].key=KEYMAX; // INIT TO HUGE KEY
    cpu->part[i].newkey=KEYMAX; // INIT TO HUGE KEY
  }
  cpu->firstpart=(unsigned long *)calloc(param->lmax+1,sizeof(unsigned long));
  cpu->npart=(unsigned long *)calloc(param->lmax+1,sizeof(unsigned long));
  cpu->nparttotal=0;
#endif
 {
      unsigned long ktest;
      ktest=2788426;
      struct CELL *newcell=getcell(&ktest,7,cpu);
      printf("allocktest= %p\n",newcell);
    }

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
void build_init_grid(struct CPU *cpu,struct PARAM *param){
  int level;
  int idx=0;
  unsigned long code;

  for(level=0;level<=param->lcoarse;level++){
    unsigned int n1=1<<level;
    unsigned int i,j,k;
    unsigned long nlev=0;
    REAL dx=1./(1<<level);
    // cartesian sweep
    for(k=0;k<n1;k++){
      for(j=0;j<n1;j++){
	for(i=0;i<n1;i++){
	  LC2M(&code,i,j,k,level);
	  cpu->grid[idx].key=code;
	  if(level<param->lcoarse) cpu->grid[idx].child=1;

	  cpu->grid[idx].gdata.d=0.;
	  cpu->grid[idx].gdata.p=0.;

	  nlev++;
	  idx++;
	}
      }
    }

    cpu->ncell[level]=nlev;
    cpu->ncelltotal+=nlev;
    printf("%ld\n",cpu->ncell[level]);
  }
}
