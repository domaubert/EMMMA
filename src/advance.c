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

void Advance_level(unsigned int level, REAL *adt, REAL *aaexp, REAL *atime, int *ndt, struct CPU *cpu, struct PARAM *param){
  
  
  int nsub=param->nsubcycles;
  REAL dt=0.; // local time (origin first call);
  int is=0;
  
  //--------------------- inner loop start
  do{
    
    if(cpu->rank==RANK_DISP){
      printf("----\n");
      printf("subscyle #%d subt=%e nsub=%d ndt=%d\n",is,dt,nsub,ndt[level-1]);
    }
    
    // compute local expansion factor
    
    
#if 0
    // REFINE/DESTROY ===========================================
    if(ndt[level]==0){ // we only refine when entering the level
      if(level>param->lcoarse) {
	destroy_cell(level,cpu,param);
	reorg(cpu,param); //reorganising the grid
      }
      
      create_cell(level,cpu,param);
      reorg(cpu,param); //reorganising the grid
      
#ifdef PIC
      amr_update_key_part(level,cpu,param);  //update particles key
      reorgpart(cpu,param); // reorganising the grid
#endif
    }
    for(int ll=param->lcoarse;ll<=param->lmax;ll++){
      printf("ll=%d FP=%lu key=%lu npart=%lu\n",ll,cpu.firstpart[ll],cpu.part[cpu.firstpart[ll]].key,cpu.npart[ll]); 
    } 
#endif


    // COMPUTING TIMESTEPS ======================================
    REAL dtnew;
    dtnew=param->dt;
    
    REAL dtpic;
    dtpic=L_comptstep(level,cpu,param);
    printf("dtpic=%e dtnew=%e\n",dtpic,dtnew);
    dtnew=(dtpic<dtnew?dtpic:dtnew);
    adt[level]=dtnew;
    
    // CIC ======================================================
    cic(level,cpu,param);
    
    // POISSON SOLVER ==========================================
    FillDens(level,cpu,param);
    PoissonSolver(level,param,cpu,aaexp[level]);
    PoissonForce(level,cpu,param,aaexp[level]);


    // RECURSIVE CALL =========================================
    if(level<param->lmax){
      if(cpu->ncell[level+1]){
	Advance_level(level+1,adt,aaexp,atime,ndt,cpu,param);
      }

    }
    

    // MOVING PARTICLES ========================================
    L_accelpart(level,cpu,param,adt,is);
    L_movepart(level,cpu,param,adt,is);
    update_key_part(level,cpu,param);  //update particles key
    reorgpart(cpu,param); // we sort them again

    // MARKING CURRENT LEVEL ===================================
    int ns;
    for(ns=0;ns<param->nsmooth;ns++){
      mark_child(level,cpu,param,ns);
      mark_nei(level,cpu,param,ns);
      mark_phy(level,cpu,param,ns);
    }
    
    // BOOKKEEPING ==============================================
    dt+=adt[level]; // advance local time
    atime[level]+=adt[level]; // advance local time
    aaexp[level]=interp_aexp(atime[level],(double *)param->cosmo->tab_aexp,(double *)param->cosmo->tab_ttilde); // update local expansion factor
    ndt[level]++;
    is++;


  }while((dt<adt[level-1])&&(is<nsub));


}
