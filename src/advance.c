#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prototypes.h"
#include "constant.h"
#include "morton.h"
#include "io.h"
#include "amr.h"
#include "poisson.h"
#include "cic.h"
#include "parameters.h"
#include "ic.h"
#include "friedmann.h"
#include "particle.h"

//void Advance_level(unsigned int level, REAL *adt, REAL *aaexp, REAL *atime, int *ndt, struct CPU *cpu, struct PARAM *param){
REAL Advance_level(unsigned int level, struct CPU *cpu, struct PARAM *param){
  
  
  int nsub=param->nsubcycles;
  REAL dt=0.; // local time (origin first call);
  int is=0;
  int ll;
  REAL dtfine;
  long int nparttotal;
  long int ncelltotal;

  printf("Entering Level %d\n",level);

  //--------------------- inner loop start
  do{
    
    if(cpu->rank==RANK_DISP){
      printf("----\n");
      printf("subscyle #%d subt=%e(->%e) nsub=%d aexp=%e time=%e\n",is,dt,cpu->adt[level-1],nsub,cpu->aaexp[level],cpu->atime[level]);
    }
    
    // compute local expansion factor
    
    
#if 1
    {
      unsigned long ktest;
      ktest=3932848;
      struct CELL *newcell=getcell(&ktest,7,cpu);
      printf("ktest= %p\n",newcell);
    }

    // REFINE/DESTROY ===========================================
    if((level<param->lmax)&&(param->lmax!=param->lcoarse)) {
      //if(is==0){ // we only refine when entering the level
      	destroy_cell(level,cpu,param);
	reorg(cpu,param); //reorganising the grid
	
	create_cell(level,cpu,param);
	reorg(cpu,param); //reorganising the grid

	amr_update_key_part(level,cpu,param);  //update particles key
	reorgpart(cpu,param); // reorganising the grid
	//}
    }

    
    //===================
    {
      unsigned long ktest;
      ktest=3932848;
      struct CELL *newcell=getcell(&ktest,7,cpu);
      printf("ktest= %p\n",newcell);
    }
    
    printf("CHECK POST REFINE\n");
    nparttotal=0.;
    ncelltotal=0.;
    for(ll=param->lcoarse;ll<=param->lmax;ll++){
      if(cpu->npart[ll]==0) continue;
      printf("ll=%d FP=%lu key=%lu npart=%lu ncell=%lu\n",ll,cpu->firstpart[ll],cpu->part[cpu->firstpart[ll]].key,cpu->npart[ll],cpu->ncell[ll]); 
      nparttotal+=cpu->npart[ll];
      ncelltotal+=cpu->ncell[ll];
    }
    printf("total number of particles across levels =%ld // of cells =%lu\n",nparttotal,ncelltotal);


    // END REFINE/DESTROY========================================
#endif


    // COMPUTING TIMESTEPS ======================================
    REAL dtnew;
    dtnew=param->dt;
    
    REAL dtpic;
    dtpic=L_comptstep(level,cpu,param);
    dtnew=(dtpic<dtnew?dtpic:dtnew);
    cpu->adt[level]=dtnew;

   

    



    // CIC ======================================================
    printf("CIC on level %d\n",level);
    cic(level,cpu,param);

    // POISSON SOLVER ==========================================
    FillDens(level,cpu,param);
    PoissonSolver(level,param,cpu);
    PoissonForce(level,cpu,param);


    // RECURSIVE CALL =========================================
    if(level<param->lmax){
      if(cpu->ncell[level+1]){
	dtfine=Advance_level(level+1,cpu,param);    

	// synchronization of levels
	cpu->adt[level]=dtfine;

      }

    }
    
    if(level==param->lcoarse) cpu->adt[level-1]=cpu->adt[level];    // to prevent lcoarse subcycle
    // MOVING PARTICLES ========================================
    L_accelpart(level,cpu,param,is);
    L_movepart(level,cpu,param,is);
    update_key_part(level,cpu,param);  //update particles key
    reorgpart(cpu,param); // we sort them again

    // MARKING CURRENT LEVEL ===================================
    
    clean_mark(level,cpu,param);
    
    if(level<param->lmax){
      int ns;
      for(ns=0;ns<param->nsmooth;ns++){
	mark_child(level,cpu,param,ns);
	mark_nei(level,cpu,param,ns);
	mark_phy(level,cpu,param,ns);
      }
    }
      
    
    // BOOKKEEPING ==============================================
    dt+=cpu->adt[level]; // advance local time
    cpu->atime[level]+=cpu->adt[level]; // advance local time
    cpu->aaexp[level]=interp_aexp(cpu->atime[level],(double *)param->cosmo->tab_aexp,(double *)param->cosmo->tab_ttilde); // update local expansion factor
    is++;


    printf("CHECK POST END\n");
    nparttotal=0.;
    ncelltotal=0.;
    for(ll=param->lcoarse;ll<=param->lmax;ll++){
      if(cpu->npart[ll]==0) continue;
      printf("ll=%d FP=%lu key=%lu npart=%lu ncell=%lu\n",ll,cpu->firstpart[ll],cpu->part[cpu->firstpart[ll]].key,cpu->npart[ll],cpu->ncell[ll]); 
      nparttotal+=cpu->npart[ll];
      ncelltotal+=cpu->ncell[ll];
    }
    printf("total number of particles across levels =%ld // of cells =%lu\n",nparttotal,ncelltotal);


  }while((dt<cpu->adt[level-1])&&(is<nsub));


  printf("Exiting Level %d\n",level);

  return dt;


}
