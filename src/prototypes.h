#include "constant.h"


#ifdef SINGLEPRECISION
// SINGLE PRECISION CASE

typedef float REAL;
#ifdef WMPI
#define MPI_REEL MPI_FLOAT
#endif

#define POW(A,B) powf(A,B)
#define SQRT(A) sqrtf(A)
#define EXP(A) expf(A)
#define FMIN(A,B) fminf(A,B)
#define FMAX(A,B) fmaxf(A,B)
#define FABS(A) fabsf(A)
#define CUDPP_REAL CUDPP_FLOAT

#else
// DOUBLE PRECISION CASE (BY DEFAULT)

typedef double REAL;
#ifdef WMPI
#define MPI_REEL MPI_DOUBLE
#endif

#define POW(A,B) pow(A,B)
#define SQRT(A) sqrt(A)
#define EXP(A) exp(A)
#define FMIN(A,B) fmin(A,B)
#define FMAX(A,B) fmax(A,B)
#define FABS(A) fabs(A)
#define CUDPP_REAL CUDPP_DOUBLE

#endif // SINGLEPRECISION

// ===================================


struct UNITS{
  double unit_l;///< comoving length size of the box [meters]
  double unit_v;///< unit velocity
  double unit_t;///< unit time [seconds]
  double unit_n;///< unit number [moles typically]
  double unit_mass;///< unit mass [in kg, total mass is equal to one in unit codes]
  double unit_d;///< density unit [typically Omegam*rhoc in kg/m3]
  double unit_N;///< number density unit [typically Omegam*rhoc/mp in 1./m3]
};

struct SCALE{
  REAL l;
  REAL v;
  REAL t;
  REAL d;
  REAL p;
  REAL E;
  REAL n;
  REAL mass;
  REAL N;
};


struct COSMOPARAM{
  REAL aexp;
  REAL om;
  REAL ov;
  REAL ob;
  REAL H0;
  REAL *tab_aexp;
  REAL *tab_ttilde;
  REAL tsim;
  REAL unit_l;
  REAL tphy;
};


#define OUTPUTPARAM_n_field_max (50)
struct OUTPUTPARAM{
  int n_field;
  int n_field_tot;

  char *field_name[OUTPUTPARAM_n_field_max];
  int field_id[OUTPUTPARAM_n_field_max];

  int field_state_grid[OUTPUTPARAM_n_field_max];
  int field_state_stat[OUTPUTPARAM_n_field_max];
  int field_state_movie[OUTPUTPARAM_n_field_max];

  int n_field_movie;

};

struct GDATA{
  REAL d;
  REAL p;
  REAL pnew;
  REAL res;
  REAL f[3];
};


struct PART{
  unsigned long key;
  unsigned long newkey;
  unsigned long idx;
  unsigned int level;
  unsigned int is;
  REAL x[3];
  REAL v[3];
  REAL mass;
  
};

struct CELL{
  unsigned long key;
  unsigned int rank;
  struct GDATA gdata;
#ifdef PIC
  REAL cicdens;
#endif

  char marked;
  char child;
  char desflag; // for destruction
};

struct CPU{
  struct CELL *grid;
  unsigned long *firstcell;
  unsigned long *ncell;
  unsigned long ncelltotal;


#ifdef PIC
  struct PART *part;
  unsigned long *firstpart;
  unsigned long *npart;
  unsigned long nparttotal;
#endif
  
  unsigned int ndumps;
  unsigned int rank;
  unsigned int nproc;
  unsigned long key_coarse_min;
  unsigned long key_coarse_max;

  REAL *adt; // array of dts (1 per level)
  REAL *aaexp; // array of expansion factors (1 per level)
  REAL *atime; // array of times (1 per level)
};

struct FIELD_INFO{
  double min;
  double max;
  double mean;
  double sigma;

  double bin_min;
  double bin_max;
  double bins_edges[N_BIN_PDF+1];
  double pdf[N_BIN_PDF];
};


struct PHYSICAL_STATE{

  int n_field;
  struct FIELD_INFO field[OUTPUTPARAM_n_field_max];

  double sfr;
  double v;

  double src;
  double mstar;
  double mstar_sfr;
  double t;

  int max_level;
  int Nsn;

  REAL dt_ff;
  REAL dt_hydro;
  REAL dt_cosmo;
  REAL dt_pic;
  REAL dt_rad;
};


struct RUN{
  int nstepstart;
  REAL aexpdump;
};



struct PARAM{
  unsigned int lmax;
  unsigned int lcoarse;
  unsigned int ngridmax;
  unsigned int npartmax;
  unsigned int nsmooth;
  unsigned int nbuff; ///< the mpi buffer size

  REAL dt; ///< the timsestep
  REAL tmax; ///< the simulation stops at tmax : corresponds to amax in cosmo
  REAL time_max; ///< for cosmo only : contains the time equivalent to amax (contained in tmax, yeah its obfuscated)


  int ndumps; ///< the frequency of outputs
  REAL dt_dump; ///< the physical time between 2 dumps in years
  int nsteps; ///< the maximal number of timesteps

  int nrestart; ///< the restart snapshot
  int nsubcycles; ///< number of subcyles in AMR advance procedure

  int mgridlmin; 
  unsigned int niter;
  unsigned int nrelax;
  unsigned int nvcycles;
  REAL poissonacc;

  struct OUTPUTPARAM *out_grid;
  struct OUTPUTPARAM *out_part;

  REAL amrthresh0;
  REAL amrthresh; ///< the refinement criterion (refine if mcell>amrthresh)

  int DM_res; ///< resolution of dark matter particle (equivalent level of lcoarse + DM_res)
  REAL dx_res; ///< maximum spatial resolution before blocking AMR in Parsec

  int gstride; ///< the size of the stencil for vector based computations (gravity) // OBSOLETE
  int hstride; ///< the size of the stencil for vector based computations (hydro) // OBSOLETE
  int nthread; ///< number of GPU threads
  int nstream; ///< number of GPU streams
  int ompthread; ///< numberf of OMP threads

  char paramrunfile[256]; // contains the parameter file name used on command line

  struct COSMOPARAM *cosmo; ///< the cosmological parameters
  struct RUN *run;// current run parameters

  struct PHYSICAL_STATE *physical_state;

  struct UNITS unit; ///< contains the units
  struct SCALE scale; ///< contains the scaling factor for units convertion (function of aexp)

};
