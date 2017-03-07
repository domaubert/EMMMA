#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <sys/stat.h>
#include <math.h>

#include "prototypes.h"
#include "io.h" // assign_grid_field
#include "constant.h"

#if defined(UVBKG) || defined(STARS_TO_UVBKG)
#include "src_utils.h" //setUVBKG
#endif

int copy_file(char const * const source, char const * const destination) {
/**
  * copy a file name 'source' into 'destination'
  */
  FILE* fSrc;
  FILE* fDest;
  char buffer[512];
  int NbLus;

  if ((fSrc=fopen(source,"rb"))==NULL) {return 1;}

  if ((fDest=fopen(destination,"wb"))==NULL){
      fclose(fSrc);
      return 2;
  }

  while ((NbLus = fread(buffer, 1, 512, fSrc)) != 0)
      fwrite(buffer, 1, NbLus, fDest);

  fclose(fDest);
  fclose(fSrc);

  return 0;
}

void copy_param(const char *folder){
/**
  * copy all parameters files into each sudfolders
  */

  char param[128];
  char param_src[128];
  char param_dest[128];

  sprintf(param,"param.run");
  sprintf(param_src,"%s%s","data/",param);
  sprintf(param_dest,"%s%s",folder,param);
  copy_file(param_src, param_dest);

  sprintf(param,"src/param.mk");
  sprintf(param_src,"%s%s","data/",param);
  sprintf(param_dest,"%s%s",folder,param);
  copy_file(param_src, param_dest);

  sprintf(param,"param.info");
  sprintf(param_src,"%s%s","data/",param);
  sprintf(param_dest,"%s%s",folder,param);
  copy_file(param_src, param_dest);

#ifdef GPUAXL
  sprintf(param,"param.avg.gpu");
#else
  sprintf(param,"param.avg.cpu");
#endif

  sprintf(param_src,"%s%s","data/",param);
  sprintf(param_dest,"%s%s",folder,param);
  copy_file(param_src, param_dest);

  sprintf(param,"param.h");
  sprintf(param_src,"%s%s","data/",param);
  sprintf(param_dest,"%s%s",folder,param);
  copy_file(param_src, param_dest);
#ifdef ALLOCT
  sprintf(param,"param.output");
  sprintf(param_src,"%s%s","data/",param);
  sprintf(param_dest,"%s%s",folder,param);
  copy_file(param_src, param_dest);
  #endif // ALLOCT
}

void printFileOnScreen(char *filename){
/**
  * print file on screen
  */

  FILE* fp=fopen(filename,"r");
  if(fp == NULL){
    printf("Cannot open %s : you may want to link it in the current directory\n", filename);
  }
  char ch;
  while((ch=fgetc(fp))!=EOF) printf("%c",ch);
}

void dumpFile(char *filename_in, char *filename_out){
/**
  * copy filename_in into filename_out and print it on screen
  */
  int fileok=1;
  FILE *fps[2] = {stdout, NULL};
  fps[1]=fopen(filename_out,"w");
  if(fps[1] == NULL) {
    printf("Cannot open %s\n", filename_out);
    fileok=0;
  }

  FILE* buf=NULL;
  buf=fopen(filename_in,"r");
  if(buf == NULL){
    printf("Cannot open %s : you may want to link it in the current directory\n", filename_in);
    fileok=0;
  }

  int i;
  if(fileok){
    for(i=0;i<2;i++){
      FILE *fp = fps[i];
      char ch;
      fseek(buf,0,SEEK_SET);
      while((ch=fgetc(buf))!=EOF) fprintf(fp,"%c",ch);
    }

    fclose(fps[1]);
    fclose(buf);
  }
}

void readOutputParam_grid(char *fparam, struct PARAM *param){
  int debug=0;

  char *field_name [] ={
    // The field order has to be the same as in param.output for consistency
    "x",
    "y",
    "z",
    "l",
    "cpu",
    "gdata_d",
    "density",
    "gdata.p",
    "res",
    "f0",
    "f1",
    "f2",
    "marked",
    "field_d",
    "field_u",
    "field_v",
    "field_w",
    "field_p",
    "rfield_e0",
    "rfield_fx0",
    "rfield_fy0",
    "rfield_fz0",
    "rfield_e1",
    "rfield_fx1",
    "rfield_fy1",
    "rfield_fz1",
    "rfield_snfb",
    "rfield_e2",
    "rfield_fx2",
    "rfield_fy2",
    "rfield_fz2",
    "rfield_src",
    "xion",
    "field_dX",
    "field_dXHE",
    "field_dXXHE",
    "field_xHE",
    "field_xxHE",
    "rfield_temp",
    "z_first_xion",
    "z_last_xion"
  };

  int n_field=0;
  int n_field_tot=0;
  param->out_grid->n_field_movie=0;

  FILE *f=fopen(fparam,"r");
  if(f==NULL){
      printf("ERROR : cannot open the parameter file (%s given), please check\n",fparam);
      abort();
  }

  if (debug) printf("opening OK\n");

  char stream[256];

  while (fscanf(f,"%s",stream)!= EOF) {
    param->out_grid->field_name[n_field_tot]= field_name[n_field_tot];
    if( !strncmp(stream, "#" ,1)){
      // skip line if first char == #
      char* line=NULL;
      size_t len=0;
      size_t status= getline(&line,&len,f);
      if (debug) printf("%s\t",line);
      free(line);
      continue;
    }

#ifdef SINGLEPRECISION
    //char* type= "%d %d %d %e %e\n";
char* type= "%d %d %d\n";
#else
//char* type= "%d %d %d %le %le\n";
char* type= "%d %d %d\n";
#endif // SINGLEPRECISION

    size_t status=fscanf(f, type,
          &(param->out_grid->field_state_grid[n_field_tot]),
          &(param->out_grid->field_state_movie[n_field_tot]),
          &(param->out_grid->field_state_stat[n_field_tot])
			 //          &(param->physical_state->field[n_field_tot].bin_min),
			 //&(param->physical_state->field[n_field_tot].bin_max)
          );

    if (debug) printf("%s\t",stream);
    //if (debug) printf("grid=%d stat=%d\t",param->out_grid->field_state_grid[n_field_tot],  param->out_grid->field_state_stat[n_field_tot]);
    //if (debug) printf("bin_min=%e bin_max=%e\n", param->physical_state->field[n_field_tot].bin_min,  param->physical_state->field[n_field_tot].bin_max);


    if (param->out_grid->field_state_movie[n_field_tot]) param->out_grid->n_field_movie++;

    if (param->out_grid->field_state_grid[n_field_tot]){
      param->out_grid->field_id[n_field_tot]=1;

      n_field++;
    }else{
      param->out_grid->field_id[n_field_tot]=0;
    }
    n_field_tot++;
  }

  if (debug) printf("n_field_tot=%d\n",n_field_tot);

  fclose(f);

  if (debug) printf("read OK\n");

  param->out_grid->n_field=n_field;
  param->out_grid->n_field_tot=n_field_tot;

  if(debug){
    int i;
    for (i=0;i<n_field_tot; i++){
      if (param->out_grid->field_state_movie[i])
      printf("%d\t%s\n",param->out_grid->field_id[i], param->out_grid->field_name[i]);
    }
  abort();
  }

}

void readOutputParam_part(char *fparam, struct PARAM *param){

  int debug=0;

  char *field_name [] ={
    // The field order has to be the same as in param.output for consistency
    "x",
    "y",
    "z",
    "vx",
    "vy",
    "vz",
    "fx",
    "fy",
    "fz",
    "idx",
    "isStar",
    "epot",
    "ekin",
    "mass",
    "age"
  };

  int n_field=0;
  int n_field_tot=0;

  FILE *f=NULL;
  f=fopen(fparam,"r");
  if(f==NULL){
      printf("ERROR : cannot open the parameter file (%s given), please check\n",fparam);
      abort();
  }
  char stream[256];
  int state;

  while (fscanf(f,"%s",stream)!= EOF) {
  param->out_part->field_name[n_field_tot]= field_name[n_field_tot];
  if( !strncmp(stream, "#" ,1)){
    // skip line if first char == #
    char* line=NULL;
    size_t len=0;
    size_t status= getline(&line,&len,f);
    if (debug) printf("%s\t",line);
    free(line);
    continue;
  }

  size_t status=fscanf(f, "%d\n",&state);
    if (state){
      param->out_part->field_id[n_field_tot] = 1;
      param->out_part->field_name[n_field_tot] = field_name[n_field_tot];
      n_field++;
    }else{
      param->out_part->field_id[n_field_tot] = 0;
    }
    n_field_tot++;
  }
  fclose(f);


  param->out_part->n_field=n_field;
  param->out_part->n_field_tot=n_field_tot;

  if(debug){
    int i;
    for (i=0;i<n_field_tot; i++){
      if (param->out_part->field_id[i])
      printf("%d\t%s\n",param->out_part->field_id[i], param->out_part->field_name[i]);
    }
    abort();
  }

}


void ReadParameters(char *fparam, struct PARAM *param){
  int debug=0;
  FILE *buf=NULL;
  char stream[256];
  size_t rstat;
  double dummyf;
  char RF[]="%s %lf";

  buf=fopen(fparam,"r");
  if(buf==NULL)
    {
      printf("ERROR : cannot open the parameter file (%s given), please check\n",fparam);
      abort();
    }
  else
    {

      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->ngridmax);  if (debug) printf("param->ngridmax=%d\n", param->ngridmax);
      rstat=fscanf(buf,"%s %d",stream,&param->npartmax);  if (debug) printf("param->npartmax=%d\n", param->npartmax);
      rstat=fscanf(buf,"%s %d",stream,&param->nbuff);     if (debug) printf("param->nbuff=%d\n", param->nbuff);

      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->ndumps);    if (debug) printf("param->ndumps=%d\n", param->ndumps);
      rstat=fscanf(buf,RF,stream,&param->dt_dump);        if (debug) printf("param->dt_dump=%e\n", param->dt_dump);
      rstat=fscanf(buf,"%s %d",stream,&param->nsteps);    if (debug) printf("param->nsteps=%d\n", param->nsteps);
      rstat=fscanf(buf,RF,stream,&dummyf);param->dt=(REAL)dummyf;   if (debug) printf("param->dt=%e\n", param->dt);
      rstat=fscanf(buf,RF,stream,&dummyf);param->tmax=(REAL)dummyf; if (debug) printf("param->tmax=%e\n", param->tmax);

      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->lcoarse); if (debug) printf("param->lcoarse=%d\n", param->lcoarse);
      rstat=fscanf(buf,"%s %d",stream,&param->lmax);    if (debug) printf("param->lmax=%d\n", param->lmax);
      rstat=fscanf(buf,RF,stream,&dummyf);param->amrthresh0=(REAL)dummyf; if (debug) printf("param->amrthresh=%e\n", param->amrthresh0);
      rstat=fscanf(buf,"%s %d",stream,&param->nsmooth); if (debug) printf("param->nsmooth=%d\n", param->nsmooth);


      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->DM_res);  if (debug) printf("param->DM_res=%d\n", param->DM_res);
      rstat=fscanf(buf,RF,stream,&dummyf);param->dx_res=(REAL)dummyf; if (debug) printf("param->dx_res=%e\n", param->dx_res);

      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->niter);   if (debug) printf("param->niter=%d\n", param->niter);
      rstat=fscanf(buf,RF,stream,&dummyf);param->poissonacc=(REAL)dummyf; if (debug) printf("param->poissonacc=%e\n", param->poissonacc);
      rstat=fscanf(buf,"%s %d",stream,&param->mgridlmin); if (debug) printf("param->mgridlmin=%d\n", param->mgridlmin);
      if(param->mgridlmin<0){
        param->mgridlmin=param->lcoarse-param->lcoarse;
      }

      rstat=fscanf(buf,"%s %d",stream,&param->nvcycles);  if (debug) printf("param->nvcycles=%d\n", param->nvcycles);
      rstat=fscanf(buf,"%s %d",stream,&param->nrelax);    if (debug) printf("param->nrelax=%d\n", param->nrelax);

      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->nrestart);  if (debug) printf("param->nrestart=%d\n", param->nrestart);

      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->gstride);    if (debug) printf("param->gstride=%d\n", param->gstride);
      rstat=fscanf(buf,"%s %d",stream,&param->hstride);    if (debug) printf("param->hstride=%d\n", param->hstride);
      rstat=fscanf(buf,"%s %d",stream,&param->nsubcycles); if (debug) printf("param->nsubcycles=%d\n", param->nsubcycles);
      rstat=fscanf(buf,"%s %d",stream,&param->nthread);    if (debug) printf("param->nthread=%d\n", param->nthread);
      rstat=fscanf(buf,"%s %d",stream,&param->nstream);    if (debug) printf("param->nstream=%d\n", param->nstream);
      rstat=fscanf(buf,"%s %d",stream,&param->ompthread);  if (debug) printf("param->ompthread=%d\n", param->ompthread);

      rstat=fscanf(buf,"%s",stream);
#ifdef WRAD
      rstat=fscanf(buf,RF,stream,&dummyf);param->clight=(REAL)dummyf;param->clightorg=(REAL)dummyf; if (debug) printf("param->clight=%e\n", param->clight);
      rstat=fscanf(buf,RF,stream,&dummyf);param->denthresh=(REAL)dummyf;  if (debug) printf("param->denthresh=%e\n", param->denthresh);
      rstat=fscanf(buf,RF,stream,&dummyf);param->tmpthresh=(REAL)dummyf;  if (debug) printf("param->tmpthresh=%e\n", param->tmpthresh);
      rstat=fscanf(buf,RF,stream,&dummyf);param->srcint=(REAL)dummyf;     if (debug) printf("param->srcint=%e\n", param->srcint);
      rstat=fscanf(buf,RF,stream,&dummyf);param->fesc=(REAL)dummyf;       if (debug) printf("param->fesc=%e\n", param->fesc);

      param->srcint*=param->fesc;

      char filename[256];
      rstat=fscanf(buf,"%s %s",stream, filename);
      sprintf(param->atomic.path,"./SRC/src/atomic_data/%s",filename);   if (debug) printf("param->atomic.path=%s\n", param->atomic.path);
      param->fudgecool=1.0;
      param->ncvgcool=0;
#else
      int i;
	for (i=0; i<5; i++)	rstat=fscanf(buf,RF,stream,&dummyf);
                      rstat=fscanf(buf,"%s %s",stream, stream);
#endif

      rstat=fscanf(buf,"%s",stream);
#ifdef STARS
      rstat=fscanf(buf,RF,stream,&dummyf);param->stars->overdensity_cond=(REAL)dummyf;   if (debug) printf("param->stars->overdensity_cond=%e\n", param->stars->overdensity_cond);
      rstat=fscanf(buf,RF,stream,&dummyf);param->stars->density_cond=(REAL)dummyf;       if (debug) printf("param->stars->density_cond=%e\n", param->stars->density_cond);
      rstat=fscanf(buf,RF,stream,&dummyf);param->stars->efficiency=(REAL)dummyf;         if (debug) printf("param->stars->efficiency=%e\n", param->stars->efficiency);
      rstat=fscanf(buf,RF,stream,&dummyf);param->stars->tlife=(REAL)dummyf;              if (debug) printf("param->stars->tlife=%e\n", param->stars->tlife);
      rstat=fscanf(buf,RF,stream,&dummyf);param->stars->mass_res=(REAL)dummyf;           if (debug) printf("param->stars->mass_res=%e\n", param->stars->mass_res);

#else
	for (i=0; i<5; i++)	rstat=fscanf(buf,RF,stream,&dummyf);
#endif

      rstat=fscanf(buf,"%s",stream);
#ifdef SUPERNOVAE
      rstat=fscanf(buf,RF,stream,&dummyf);param->sn->feedback_eff	=(REAL)dummyf;        if (debug) printf("param->sn->feedback_eff=%e\n", param->sn->feedback_eff);
      rstat=fscanf(buf,RF,stream,&dummyf);param->sn->feedback_frac	=(REAL)dummyf;      if (debug) printf("param->sn->feedback_frac=%e\n", param->sn->feedback_frac);
      rstat=fscanf(buf,RF,stream,&dummyf);param->sn->ejecta_proportion	=(REAL)dummyf;  if (debug) printf("param->sn->ejecta_proportion=%e\n",param->sn->ejecta_proportion);
      rstat=fscanf(buf,RF,stream,&dummyf);param->sn->sn_egy	=(REAL)dummyf;              if (debug) printf("param->sn->sn_egy=%e\n", param->sn->sn_egy);
      rstat=fscanf(buf,RF,stream,&dummyf);param->sn->tlife	=(REAL)dummyf;              if (debug) printf("param->sn->tlife=%e\n", param->sn->tlife);
      rstat=fscanf(buf,RF,stream,&dummyf);param->sn->load_factor =(REAL)dummyf;         if (debug) printf("param->sn->load_factor=%e\n", param->sn->load_factor);

#else
	for (i=0; i<5; i++)	rstat=fscanf(buf,RF,stream,&dummyf);
#endif

#ifdef MOVIE
      rstat=fscanf(buf,"%s",stream);
      rstat=fscanf(buf,"%s %d",stream,&param->movie->lmap);                    if (debug) printf("param->movie->lmap=%d\n", param->movie->lmap);
      if (param->movie->lmap>param->lmax) param->movie->lmap=param->lmax;
      rstat=fscanf(buf,"%s %s",stream,param->movie->mode_str);                if (debug) printf("param->movie->mode_str=%s\n", param->movie->mode_str);

      rstat=fscanf(buf,RF,stream,&param->movie->xmin);    if (debug) printf("param->movie->xmin=%e\n", param->movie->xmin);
      rstat=fscanf(buf,RF,stream,&param->movie->xmax);    if (debug) printf("param->movie->xmax=%e\n", param->movie->xmax);
      rstat=fscanf(buf,RF,stream,&param->movie->ymin);    if (debug) printf("param->movie->ymin=%e\n", param->movie->ymin);
      rstat=fscanf(buf,RF,stream,&param->movie->ymax);    if (debug) printf("param->movie->ymax=%e\n", param->movie->ymax);
      rstat=fscanf(buf,RF,stream,&param->movie->zmin);    if (debug) printf("param->movie->zmin=%e\n", param->movie->zmin);
      rstat=fscanf(buf,RF,stream,&param->movie->zmax);    if (debug) printf("param->movie->zmax=%e\n", param->movie->zmax);
#endif
      fclose(buf);
    }

 if (debug) abort();
}

// =============================== =============================== ===============================
// =============================== =============================== ===============================

void GetParameters(char *fparam, struct PARAM *param){

  ReadParameters(fparam, param);


  // ====================== some checks

  // stencil/streams conformity
#ifdef GPUAXL
  if(param->hstride<(param->nthread*param->nstream)){
    printf(" Stream Thread granulosity too high : nt=%d ns=%d stencil=%d\n",param->hstride,param->nthread,param->nstream);
    abort();
  }
#endif

#ifdef STARS
    param->stars->n		= 0;
#ifdef AGN
    param->stars->nagn		= 0;
#endif
#endif

#ifdef SRCINT
  param->srcint*=SRCINT;
#endif

#if defined(UVBKG) || defined(STARS_TO_UVBKG)
  setUVBKG(param, "SRC/src/phys_data/uvbkg.dat");
#endif // UVBKG


#ifdef WRAD
  readAtomic(param);
#endif // WRAD

#ifdef SUPERNOVAE
  //read_egy_loss(param);
  //read_mass_loss(param);
#endif // SUPERNOVAE


}

// ========================================================================================================

void dumpInfo(char *filename_info, struct PARAM *param, struct CPU *cpu){
/**
  * Write a file containing information about the simulation like the cosmology or the resolution.
  */
  FILE *fps[2] = {stdout, NULL};

  fps[1]=fopen(filename_info,"w");
  if(fps[1] == NULL) printf("Cannot open %s\n", filename_info);

  char* int_format = "%-24s%d\n";
  char* float_format = "%-24s%.2f\n";
  char* real_format = "%-24s%e\n";

  int i;
  for(i=0;i<2;i++){
    FILE *fp = fps[i];

    fprintf(fp, int_format,"nproc",(cpu->nproc)); 		// number of processor
    fprintf(fp, float_format,"box_size_hm1_Mpc",(param->unit.unit_l/PARSEC/1e6*param->cosmo->H0/100));
    fprintf(fp, int_format,"level_min",(param->lcoarse) );
    fprintf(fp, int_format,"level_max",(param->lmax) );

  fprintf(fp,"##=Unit_code->SI====================\n" );

  fprintf(fp, real_format,"unit_l",(param->unit.unit_l) );		// comoving length size of the box [meters]
  fprintf(fp, real_format,"unit_v",(param->unit.unit_v) );		// unit velocity
  fprintf(fp, real_format,"unit_t",(param->unit.unit_t) );		// unit time [seconds]
  fprintf(fp, real_format,"unit_N",(param->unit.unit_N) );		// unit number [moles typically]
  fprintf(fp, real_format,"unit_mass",(param->unit.unit_mass) );	// unit mass [in kg, total mass is equal to one in unit codes]
  //  fprintf(fp,"\n");

  fprintf(fp,"##=Cosmology========================\n" );
    fprintf(fp, real_format,"om",(param->cosmo->om) );			// Omega matter
    fprintf(fp, real_format,"ov",(param->cosmo->ov) );			// Omega vacuum
    fprintf(fp, real_format,"ob",(param->cosmo->ob) );			// Omega baryon
    fprintf(fp, real_format,"H0",(param->cosmo->H0) );			// Hubble constant
  //  fprintf(fp,"\n");

#ifdef PIC
    fprintf(fp,"##=Mass_resolution_(Mo)=============\n" );
    //REAL mass_res_DM =  (1.- param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*(param->lcoarse+param->DM_res))*param->unit.unit_mass/SOLAR_MASS;

    // below hack to get the correct result even in the case of SP calculation
    double munpercell=param->cosmo->om*(3.*POW(param->cosmo->H0*1e3/PARSEC/1e6,2)/(8.*M_PI*NEWTON_G))*POW(param->unit.unit_l/POW(2.0,param->lcoarse),3);
    REAL mass_res_DM=(REAL)((1.- param->cosmo->ob/param->cosmo->om)*munpercell/POW(2.0,3.*param->DM_res)/SOLAR_MASS);

    fprintf(fp, real_format,"mass_res_DM",mass_res_DM );
#ifdef STARS
    REAL res = param->stars->mass_res;
    if(res>100){
        fprintf(fp, real_format,"mass_res_star",param->stars->mass_res);
    }else{

      if(res>=0){
        /* REAL mstars_level=(param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*(param->lcoarse+param->stars->mass_res)); */
        /* REAL mass_res_star = mstars_level * param->unit.unit_mass /SOLAR_MASS; */

	double mstars_level=(param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*(param->stars->mass_res));
	REAL mass_res_star =(REAL)( mstars_level * munpercell/SOLAR_MASS);


        fprintf(fp, real_format,"mass_res_star",mass_res_star);
      }else{
        int level;
        for(level=param->lcoarse;level<=param->lmax;level++){
          REAL mlevel=level-1;
          REAL restmp=-param->stars->mass_res;
          REAL mstars_level=(param->cosmo->ob/param->cosmo->om) * POW(2.0,-3.0*(mlevel+restmp));
          REAL mass_res_star = mstars_level * param->unit.unit_mass /SOLAR_MASS;
          char mlev[128];
          sprintf(mlev,"mass_star_L%d",level);
          fprintf(fp, real_format,mlev,mass_res_star );
        }
      }
    }
#endif // STARS
#endif // PIC

    fprintf(fp,"##=Spatial_resolution_(Kpc)=========\n" );
    int level;
    for(level=param->lcoarse;level<=param->lmax;level++){
      REAL dx = param->unit.unit_l * POW(2,-level) /1e3/PARSEC;
      char dxlev[128];
      sprintf(dxlev,"dx_L%d",level);
      fprintf(fp, real_format,dxlev,dx);
    }

  }
  fclose(fps[1]);
}

void dumpHeader(struct PARAM *param, struct CPU *cpu,char *fparam){
/**
  * Dump on screen and into files, a set of parameters
  */

  printf("\n");
  printf("--------------------------------------------------------------\n");
  dumpInfo("data/param.info", param, cpu);
  printf("\n");
  printf("--------------------------------------------------------------\n");
  printFileOnScreen("SRC/param.mk");
  printf("\n");
  printf("--------------------------------------------------------------\n");
  printFileOnScreen("SRC/src/param.h");
  printf("\n");
  printf("--------------------------------------------------------------\n");
  printFileOnScreen(param->paramrunfile);
  printf("\n");
  printf("--------------------------------------------------------------\n");
#ifdef ALLOCT
  char partoutput[512];
  strcpy(partoutput,param->paramrunfile);
  strcat(partoutput,".part_output");
  printFileOnScreen(partoutput);
  printf("\n");
  printf("--------------------------------------------------------------\n");
  char gridoutput[512];
  strcpy(gridoutput,param->paramrunfile);
  strcat(gridoutput,".grid_output");
  printf("\n");
  printf("--------------------------------------------------------------\n");
#endif // ALLOCT


  REAL threshold=param->amrthresh0;
#ifndef ZOOM
  threshold*=POW(2.0,-3.0*param->lcoarse);
#else
  threshold*=POW(2.0,-3.0*param->lmaxzoom);
#endif
  if(cpu->rank==RANK_DISP)
    printf("amrthresh : maximum number of part in a cell before refinement : %d -> compute density thresold of %e \n ", (int)param->amrthresh0, threshold);


  printf("\n");
}

