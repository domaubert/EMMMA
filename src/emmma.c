
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#include "constant.h"
#include "prototypes.h"
#include "allocation.h"
#include "morton.h"
#include "io.h"
#include "amr.h"
#include "poisson.h"
#include "cic.h"
#include "parameters.h"
#include "ic.h"
#include "friedmann.h"
#include "advance.h"

#define LEVCOARSE 5
#define LEVMAX 6
#define NGRIDMAX 30000000
#define NPARTMAX 100000


int main(int argc, char *argv[]){
  // ================= ZE DATA
  struct CPU cpu;
  struct PARAM param;
  struct RUN run;
  param.run=&run;

  // ================= COSMO TABS
  double tab_aexp[NCOSMOTAB];
  double tab_ttilde[NCOSMOTAB];
  double tab_t[NCOSMOTAB];
  REAL amax;
  struct COSMOPARAM cosmo;
  param.cosmo=&cosmo;
  // ================= OUTPUTS MANAGEMENT
  struct OUTPUTPARAM out_grid;
  param.out_grid=&out_grid;

  struct OUTPUTPARAM out_part;
  param.out_part=&out_part;


  struct PHYSICAL_STATE physical_state; 
  param.physical_state = &physical_state; 


  //=========== some initial calls =============
  GetParameters(argv[1],&param); // reading the parameters file
  strcpy(param.paramrunfile,argv[1]);
  //in  cosmo case tmax is understood as a maximal expansion factor
  amax=param.tmax;


  char gridoutput[512];
  strcpy(gridoutput,param.paramrunfile);
  strcat(gridoutput,".grid_output");
  readOutputParam_grid(gridoutput, &param);
  char partoutput[512];
  strcpy(partoutput,param.paramrunfile);
  strcat(partoutput,".part_output");
  readOutputParam_part(partoutput, &param);


  cpu.rank=0;
  cpu.nproc=1;

  
  // =================== Starting Banner

  if(cpu.rank==RANK_DISP){
    printf("================================\n");
    printf("            EMMA V1.2           \n");
    printf("      Engines Are Running on    \n");
    printf("             %d process         \n",cpu.nproc);
    printf("================================\n");

    copy_file(param.paramrunfile, "data/param.run");
  }



  //==================================   reading outputlist
  char outputlist[512];
  strcpy(outputlist,param.paramrunfile);
  strcat(outputlist,".list_aexp");
  FILE *foutputs;

  run.aexpdump=0;
  float tempa;
  if((foutputs=fopen(outputlist,"r"))!=NULL){
    int dump = fscanf(foutputs,"%e",&tempa);
    run.aexpdump=tempa;
    if(cpu.rank==RANK_DISP){
      printf("Reading outputs from %s : first dump at aexp=%e\n",outputlist,run.aexpdump);
    }
  }
  else{
    if(cpu.rank==RANK_DISP)
      printf("WARNING NOT OUTPUT LIST FOUND !! \n");
  }

  // allocations ===================
  allocation(&cpu,&param);

  int ncellscoarse = POW(2,3*param.lcoarse)/8; // number of cells before refinement
  int ncellsmax    = (param.lmax>param.lcoarse?3:1) * ncellscoarse; 		 // max number of cells after refinement
  int lbfg = 2; 				 // load balancing factor for the grid
  int noon = (ncellsmax * lbfg) /cpu.nproc;	 // number of octs needed
  if (param.ngridmax < noon && cpu.rank==RANK_DISP ) {
	printf("\n");
	printf("YOU MAY NEED MORE MEMORY SPACE TO COMPUTE THE GRID\n");
	printf("%d oct allocated per processor \n",param.ngridmax);
        printf("%d oct approximately needed\n", noon);
	printf("\n");
  }


  // ================================  building the initial coarse grid
  build_init_grid(&cpu,&param);
  reorg(&cpu,&param); //reorganising the grid
  
  if(cpu.rank==RANK_DISP) printf("init grid done\n");


  //=================================  building the array of timesteps
  int level;
  REAL *adt;
  adt=(REAL *)malloc(sizeof(REAL)*(param.lmax+1));
  for(level=0;level<=param.lmax;level++) adt[level]=param.dt;

  REAL *aaexp;
  aaexp=(REAL *)malloc(sizeof(REAL)*(param.lmax+1));
  REAL *atime;
  atime=(REAL *)malloc(sizeof(REAL)*(param.lmax+1));

#ifdef COARSERAD
  REAL *adt_rad;
  adt_rad=(REAL *)malloc(sizeof(REAL)*(param.lmax+1));
  for(level=0;level<=param.lmax;level++) adt_rad[level]=param.dt;
#endif

  int *ndt;
  ndt=(int *)malloc(sizeof(int)*(param.lmax+1));

  // INITIALISATION FROM INITIAL CONDITIONS =========================
  REAL munit;
  REAL ainit,tinit,tsim;
  unsigned long npart;
  int nstepstart;

  if(param.nrestart==0){
    if(cpu.rank==RANK_DISP) printf("==> starting part\n");

    nstepstart=0;

#ifdef SPLIT
    read_split_grafic_part(&cpu, &munit, &ainit, &npart, &param, param.lcoarse);
    
#ifdef WMPI
    long ntotsplit=npart;
    MPI_Allreduce(MPI_IN_PLACE,&ntotsplit,1,MPI_LONG,MPI_SUM,cpu.comm);
    
    if(cpu.rank==RANK_DISP) printf("found %ld particles among the split ICs files\n",ntotsplit);
      
#endif
#else
      read_grafic_part(&cpu, &munit, &ainit, &npart, &param, param.lcoarse);
#endif
      
      // once particles have been read from grafic, we sort and reoganize part data
      build_init_part(&cpu,&param);
      reorgpart(&cpu,&param);
      tinit=ainit;
      tsim=tinit;

      printf("%d %d\n",param.lcoarse,param.lmax);

      for(level=param.lcoarse;level<=param.lmax;level++){
	printf("level=%d FP=%lu key=%lu npart=%lu\n",level,cpu.firstpart[level],cpu.part[cpu.firstpart[level]].key,cpu.npart[level]); 
      } 


  }
  else{

#if 0 // DISABLE RESTART COMPILATION
    
// dealing with restart
    char filename[256];
    if(cpu.rank==RANK_DISP) printf("Restarting from snap #%d\n", param.nrestart);
#ifdef PIC
    sprintf(filename,"data/bkp/part.%05d.p%05d",param.nrestart,cpu.rank);
    restore_part(filename,&tsim,&param,&cpu);
#endif

    sprintf(filename,"data/bkp/grid.%05d.p%05d",param.nrestart,cpu.rank);
    freeoct=restore_amr(filename,&tsim,&tinit,&nstepstart,&ndumps,&param,&cpu,part,adt);


    nstepstart+=1.; // next timestep is n+1
    ndumps+=1.;    // next timestep is n+1


    if(cpu.rank==RANK_DISP){
      printf(" ... Restarting from file #%d with nstep=%d tsim=%e ndumps=%d\n",param.nrestart,nstepstart,tsim,ndumps);
    }


    // prepare the next in aexplist
    if(run.aexpdump){
      while(run.aexpdump<=tsim){
	  if(fscanf(foutputs,"%e",&tempa)==EOF){
	    run.aexpdump=0;
	    break;
	  }
	  else{
	    run.aexpdump=tempa;
	  }
      }
    }
    if(cpu.rank==RANK_DISP){
      printf("Next dump in the list at aexp=%e\n",run.aexpdump);
    }

    // temporal boundaries of the full run
    ainit=tinit;

    //==================================== END Restart =================================================
#endif // ENABLE RESTART COMPILATION
  }


  // ================== COMPUTATION OF FRIEDMANN TABLES
  REAL treal,treal0,trealBB;
  REAL aexp;
  REAL tmax;
  // we compute the friedmann tables
  aexp=tsim;
  
  // at this stage we have to compute the conformal time
  tsim=-0.5*SQRT(cosmo.om)*integ_da_dt_tilde(aexp,1.0,cosmo.om,cosmo.ov,1e-8);
  
  // real times in units of 1./H0
  treal=-integ_da_dt(aexp,1.0,cosmo.om,cosmo.ov,1e-8);
  trealBB=-integ_da_dt(1e-5,1.0,cosmo.om,cosmo.ov,1e-8);
  treal0=treal;
  
  
  // interpolation table
  if(cpu.rank==RANK_DISP) printf("computing friedmann tables with ainit=%e amax=%e\n",ainit,amax);
  compute_friedmann(ainit*0.95,amax,NCOSMOTAB,cosmo.om,cosmo.ov,tab_aexp,tab_ttilde,tab_t);
  
  tmax=-0.5*SQRT(cosmo.om)*integ_da_dt_tilde(amax,1.0+1e-6,cosmo.om,cosmo.ov,1e-8);
  if(cpu.rank==RANK_DISP) printf("tmax=%e treal=%e\n",tmax,treal);
  cosmo.tab_aexp=(REAL *)tab_aexp;
  cosmo.tab_ttilde=(REAL *)tab_ttilde;
  
  param.time_max=tmax;
  printf("tmax=%e\n",tmax);
  mkdir("data/", 0755);
  //if(cpu.rank==RANK_DISP) dumpHeader(&param,&cpu,argv[1]);

  // test if each cpu will have at least one oct in the minimum level of multigrid
  int Lmin = 1+(int)(log(cpu.nproc)/(3*log(2.)));


  if( param.mgridlmin>0 && param.mgridlmin < Lmin ){
    param.mgridlmin = Lmin;
    if(cpu.rank==RANK_DISP){
      printf("Conflict between mgridlmin and ncpu : mgridlmin set to %d\n",param.mgridlmin);
    }
  }

  //================================================================================
  //
  //          AT THIS STAGE THE INITIAL SETUP HAS BEEN COMPLETED
  //
  //================================================================================

  /*  // CIC */
  /* cic(param.lcoarse,&cpu,&param); */
  


  /* // POISSON SOLVER */
  /* FillDens(param.lcoarse,&cpu,&param); */
  /* PoissonSolver(param.lcoarse,&param,&cpu,1.0); */
  
  
  /* //DUMP */
  /* dumpalloct_serial("./data/",1.0,&param, &cpu,param.lmax); */
  


#if 1 // Balise


    //==================================== MAIN LOOP ================================================
    //===============================================================================================
    // Loop over time

  int nsteps;
  for(nsteps=nstepstart;(nsteps<=param.nsteps)*(tsim<tmax);nsteps++){

    //cpu.nsteps=nsteps;
    aexp=interp_aexp(tsim,(double *)cosmo.tab_aexp,(double *)cosmo.tab_ttilde);
    cosmo.aexp=aexp;
    cosmo.tsim=tsim;
    if(cpu.rank==RANK_DISP) printf("\n============== STEP %d aexp=%e z=%lf tconf=%e tmax=%e================\n",nsteps,cosmo.aexp,1./cosmo.aexp-1.,tsim,tmax);


      // Resetting the timesteps
    
      for(level=0;level<=param.lmax;level++){
	ndt[level]=0;
	aaexp[level]=aexp;
	atime[level]=tsim;
      }


      // ---------------  Recursive Call
      level=param.lcoarse;
      //Advance_level(param.lcoarse,adt,&cpu,&param,ndt,nsteps,tsim);


      // COMPUTING TIMESTEPS ======================================
      REAL dtnew;
      dtnew=param.dt;

      REAL dtpic;
      dtpic=L_comptstep(level,&cpu,&param);
      printf("dtpic=%e dtnew=%e\n",dtpic,dtnew);
      dtnew=(dtpic<dtnew?dtpic:dtnew);
      adt[level]=dtnew;
      
      // CIC ======================================================
      cic(level,&cpu,&param);

      // POISSON SOLVER ==========================================
      FillDens(level,&cpu,&param);
      PoissonSolver(level,&param,&cpu,aaexp[level]);
      PoissonForce(level,&cpu,&param,aaexp[level]);


      
      // MOVING PARTICLES ========================================
      int is=0;
      L_accelpart(level,&cpu,&param,adt,is);
      L_movepart(level,&cpu,&param,adt,is);
      update_key_part(level,&cpu,&param);  //update particles key
      reorgpart(&cpu,&param); // we sort them again


      //==================================== timestep completed
      REAL dt=adt[param.lcoarse];
      tsim+=dt;
      

      //=================     
      REAL aexpdump=interp_aexp(tsim,(double *)cosmo.tab_aexp,(double *)cosmo.tab_ttilde);
      if(nsteps%10==0){
	dumpalloct_serial("./data/",aexpdump,&param, &cpu,param.lmax);
	dumppart_serial("./data/",aexpdump,&param, &cpu,param.lmax);
      }


#if 0
      // ==================================== dump
      cond1 = nsteps%param.ndumps==0;
      cond2 = 0;
      cond3 = tsim+adt[levelcoarse-1]>=tmax;
      cond4 = 0;

      if(run.aexpdump){
	// dumpfile at specific outputs
	cond4=cosmo.aexp>run.aexpdump;
	if(cond4){
	  if(fscanf(foutputs,"%e",&tempa)==EOF){
	    run.aexpdump=0;
	  }
	  else{
	    run.aexpdump=tempa;
	    if(cpu.rank==RANK_DISP){
	      printf("next output aexp=%e\n",run.aexpdump);
	    }
	  }
	}
      }

      if (param.dt_dump){
	cond1=0;
	int offset=0;
	
	if (nsteps==0) offset = (int)(param.cosmo->tphy/param.dt_dump);
	REAL a=param.cosmo->tphy;
	REAL b=(int)(ndumps+offset)*param.dt_dump;
	cond2=a>b;
	if(cpu.rank==RANK_DISP)printf("t=%.2e yrs next dump at %.2e yrs\n",a,b+(a>b)*param.dt_dump);
      }

      if(cond1||cond2||cond3||cond4){

	int fdump=FDUMP;
	if(cpu.nproc>fdump){
	  // dumping fields only
	  int idump;
	  for(idump=0;idump<fdump;idump++){
	    if(cpu.rank==RANK_DISP) printf("Dump batch # %d/%d\n",idump,fdump-1);
	    if(cpu.rank%fdump==idump) dumpIO(tsim+adt[param.lcoarse],&param,&cpu,adt,0);
	    sleep(1);
	  }
	}
	else{
	  dumpIO(tsim+adt[param.lcoarse],&param,&cpu,adt,0);
	}
	ndumps++;
      }

      dumpStepInfo(&param, &cpu,nsteps,adt[param.lcoarse],(float)tsim);
#endif



    }// END main loop

    
#if 0
    // ========================================= DUMP LAST SNAP
    ndumps-=1;
    int fdump=FDUMP;
    if(cpu.nproc>fdump){
      // dumping fields only
      int idump;
      for(idump=0;idump<fdump;idump++){
	if(cpu.rank==RANK_DISP && FDUMP>1) printf("Dump batch # %d/%d\n",idump,fdump-1);
	if(cpu.rank%fdump==idump) dumpIO(tsim,&param,&cpu,adt,1);
	sleep(1);
#ifdef WMPI
	MPI_Barrier(cpu.comm);
#endif
      }
    }
    else{
      dumpIO(tsim,&param,&cpu,adt,1);
    }
    // ===========================================END  DUMP LAST SNAP
#endif

#endif // BALISE COMPIL


    printf("Done\n");
    
    return 0;
}



#if 0
// ======================================================

int main(){


  struct CPU cpu;
  struct PARAM param;
  struct OUTPUTPARAM out_grid;
  int level;

  // parameters ===================
  param.ngridmax=NGRIDMAX;
  param.npartmax=NPARTMAX;
  param.lmax=LEVMAX;
  param.lcoarse=LEVCOARSE;
  param.nsmooth=2;
  param.out_grid=&out_grid;

  param.out_grid->n_field=8; // x,y,z;
  param.out_grid->n_field_tot=  param.out_grid->n_field; 
  param.out_grid->field_name[0]="x";
  param.out_grid->field_name[1]="y";
  param.out_grid->field_name[2]="z";
  param.out_grid->field_name[3]="l";
  param.out_grid->field_name[4]="cpu";
  param.out_grid->field_name[5]="gdata_d";
  param.out_grid->field_name[6]="cicdens";
  param.out_grid->field_name[7]="gdata_p";

  for(int ifield=0;ifield<param.out_grid->n_field;ifield++) param.out_grid->field_id[ifield]=1;

  param.mgridlmin=2;
  param.niter=100;
  param.nrelax=5;
  param.poissonacc=1e-5;
  param.nvcycles=3;
  // cpu parameters ================
  cpu.ndumps=0;
  cpu.rank=0;


  // allocations ===================
  allocation(&cpu,&param);

  // ====== building the initial coarse grid
  build_init_grid(&cpu,&param);
  reorg(&cpu,&param); //reorganising the grid

  // =============

#ifdef PIC
  // ====== building the initial part distrib
  build_init_part(&cpu,&param);
  reorgpart(&cpu,&param);
  
/*   int i; */
/*   for(i=0;i<param.npartmax;i++){ */
/*     printf("%f %lu // ",cpu.part[i].x[0],cpu.part[i].key); */
/*   } */

/*   for(int level=0;level<=param.lcoarse;level++){ */
/*     printf("level=%d FP=%lu key=%lu x=%f\n",level,cpu.firstpart[level],cpu.part[cpu.firstpart[level]].key,cpu.part[cpu.firstpart[level]].x[0]); */
/*   } */
#endif


#if 0
  // CIC
  cic(param.lcoarse,&cpu,&param);
#endif


#if 0
/*   struct CELL *lcell; */
/*   unsigned long kk=103; */

  for(level=0;level<=param.lcoarse;level++){
    printf("level=%d FO=%lu key=%lu ncell=%lu\n",level,cpu.firstcell[level],cpu.grid[cpu.firstcell[level]].key,cpu.ncell[level]);
  }
  printf("ncellotal=%ld\n",cpu.ncelltotal);

  PoissonSolver(param.lcoarse,&param,&cpu,1.0);
#endif

# if 1
  // ASSUMING level;
  level=param.lcoarse;

  // MARKING CURRENT LEVEL ===================================
  int ns;
  for(ns=0;ns<param.nsmooth;ns++){
    mark_child(level,&cpu,&param,ns);
    mark_nei(level,&cpu,&param,ns);
    mark_phy(level,&cpu,&param,ns);
  }
  // ==========================================================


  // REFINE/DESTROY ===========================================
  if(level>param.lcoarse) {
    destroy_cell(level,&cpu,&param);
    reorg(&cpu,&param); //reorganising the grid
  }

  create_cell(level,&cpu,&param);
  reorg(&cpu,&param); //reorganising the grid

#ifdef PIC
  amr_update_key_part(level,&cpu,&param);  //update particles key
  reorgpart(&cpu,&param);
#endif

  for(int level=param.lcoarse;level<=param.lmax;level++){
    printf("level=%d FP=%lu key=%lu npart=%lu\n",level,cpu.firstpart[level],cpu.part[cpu.firstpart[level]].key,cpu.npart[level]); 
  } 

#if 1
  // CIC =================================================
  cic(param.lcoarse,&cpu,&param);
  cic(param.lcoarse+1,&cpu,&param);
#endif

#endif
  // ==========================================================

#if 0
  for(level=0;level<=param.lmax;level++){
    printf("level=%d FO=%lu key=%lu ncell=%lu\n",level,cpu.firstcell[level],cpu.grid[cpu.firstcell[level]].key,cpu.ncell[level]);
  }

  PoissonSolver(param.lcoarse+1,&param,&cpu,1.0);

#endif

  //
  
  dumpalloct_serial("./data/",1.0,&param, &cpu,param.lmax);

  printf("Done\n");

  return 0;

}


#endif
