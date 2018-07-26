#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>

#include "prototypes.h"
#include "morton.h"
#include "constant.h"

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

}

#ifdef PIC
// =========================================================
void build_init_part(struct CPU* cpu, struct PARAM *param){
  int i;
  for(i=0;i<param->npartmax;i++){

/*     for(int p=0;p<3;p++){ */
/*       cpu->part[i].x[p]=(rand()*1.0)/RAND_MAX*0.6+0.2; */
/*     } */
    
    unsigned int x,y,z;
    unsigned long key;
    int level=param->lcoarse;
    REAL dx=1./(1<<level);
    if(cpu->part[i].mass){
      x=(int)(cpu->part[i].x[0]/dx);
      y=(int)(cpu->part[i].x[1]/dx);
      z=(int)(cpu->part[i].x[2]/dx);
      LC2M(&key,x,y,z,level);

    // depth traversal
/*     struct CELL *lcell; */
/*     lcell=getcell(&key,level,cpu); */
/*     level++; */
/*     while(lcell->child==1){ */
/*       REAL dx=1./(1<<level); */
/*       x=(int)(cpu->part[i].x[0]/dx); */
/*       y=(int)(cpu->part[i].x[1]/dx); */
/*       z=(int)(cpu->part[i].x[2]/dx); */
/*       LC2M(&key,x,y,z,level); */
      
/*       struct CELL *lcell; */
/*       lcell=getcell(&key,level,cpu); */
/*       level++; */
/*     } */

    
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

	  cpu->grid[idx].x=i*dx;
	  cpu->grid[idx].y=j*dx;
	  cpu->grid[idx].z=k*dx;

	  cpu->grid[idx].gdata.d=0.;
	  cpu->grid[idx].gdata.p=0.;

#if 0 // TO GENERATE CENTRAL BLOP
	  if((level==param->lcoarse)&&(POW(cpu->grid[idx].x+0.5*dx-0.5,2)+POW(cpu->grid[idx].y+0.5*dx-0.5,2)+POW(cpu->grid[idx].z+0.5*dx-0.5,2)<POW(0.2,2))){
	    cpu->grid[idx].gdata.d=1.;
	    cpu->grid[idx].gdata.p=1.;
	  }
	  else{
	  }
#endif
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
