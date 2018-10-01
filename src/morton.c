#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prototypes.h"
#include "constant.h"

// ===================================

unsigned long bitscat(unsigned int a)
{
  unsigned long x;
  x=a;
  x = x& 0x1fffff;
  x = (x | x << 32) & 0x1f00000000ffff;
  x = (x | x << 16) & 0x1f0000ff0000ff;
  x = (x | x << 8) & 0x100f00f00f00f00f;
  x = (x | x << 4) & 0x10c30c30c30c30c3;
  x = (x | x << 2) & 0x1249249249249249;


  return x;
}


unsigned int bitgat(unsigned long x)
{

  x = x& 0x1249249249249249;;
  x = (x ^ (x >> 2)) & 0x10c30c30c30c30c3;
  x = (x ^ (x >> 4)) & 0x100f00f00f00f00f;
  x = (x ^ (x >> 8)) & 0x1f0000ff0000ff;
  x = (x ^ (x >> 16)) & 0x1f00000000ffff;
  x = (x ^ (x >> 32)) & 0x1fffff;

  return (unsigned int) x;
}



void C2M(unsigned long *c, unsigned int x,unsigned int y, unsigned int z){
  
  *c=bitscat(x)| (bitscat(y)<<1)|bitscat(z)<<2;

}


void LC2M(unsigned long *c, unsigned int x,unsigned int y, unsigned int z, unsigned int level){
  
  *c=bitscat(x)| (bitscat(y)<<1)|bitscat(z)<<2;
  *c|=(1<<(3*level));
}

int get_level(unsigned long c){
  int level;
  for (level=0; c!=1; c>>=3, level++);
  return level;
}

void M2C(unsigned long c, unsigned int *x,unsigned int *y, unsigned int *z, unsigned int *level){
  *level=get_level(c);
  c &= ~(1 << (3*(*level)));
  *x= bitgat(c);
  *y=bitgat(c>>1);
  *z=bitgat(c>>2);

}

void assign_level(int level, unsigned long *c){
  *c=*c|(1<<(3*level));
}


int comp(const void * a, const void * b){

  struct CELL *oa;
  struct CELL *ob;

  oa=(struct CELL *) a;
  ob=(struct CELL *) b;

  int res=0;
  if(oa->key>ob->key){
    res=1;
  }
  else if(oa->key<ob->key){
    res=-1;
  }
  return res;
}

#ifdef PIC
// ===========================================
int comppart(const void * a, const void * b){

  struct PART *oa;
  struct PART *ob;

  oa=(struct PART *) a;
  ob=(struct PART *) b;

  int res=0;
  if(oa->key>ob->key){
    res=1;
  }
  else if(oa->key<ob->key){
    res=-1;
  }
  return res;
}
#endif

unsigned long floorkey(struct CELL *grid,unsigned long target, unsigned long imin, unsigned long imax){
  while(imin!=imax){
    unsigned long mid =(imin+imax)/2;
    if(grid[mid].key<target){
      imin=mid+1;
    }
    else{
      imax=mid;
    }
  }

  return imax;
}

#ifdef PIC
unsigned long floorkeypart(struct PART *part,unsigned long target, unsigned long imin, unsigned long imax){
  while(imin!=imax){
    unsigned long mid =(imin+imax)/2;
    if(part[mid].key<target){
      imin=mid+1;
    }
    else{
      imax=mid;
    }
  }

  return imax;
}
#endif

//=========================
//=========================

void key2pos(unsigned long key, REAL *x, unsigned int *level){
  unsigned int xi;
  unsigned int yi;
  unsigned int zi;

  M2C(key,&xi,&yi,&zi,level);
  REAL dx=1./(1<<(*level));

  x[0]=(xi)*dx;
  x[1]=(yi)*dx;
  x[2]=(zi)*dx;
}

//============================================================
void key2cen(unsigned long key, REAL *x, unsigned int *level){
  unsigned int xi;
  unsigned int yi;
  unsigned int zi;

  M2C(key,&xi,&yi,&zi,level);
  //printf("%ld %d %d %d\n",key,xi,yi,zi);
  REAL dx=1./(1<<(*level));

  x[0]=(xi+0.5)*dx;
  x[1]=(yi+0.5)*dx;
  x[2]=(zi+0.5)*dx;
}
//============================================================
unsigned long pos2key(REAL *pos, unsigned int level){
  unsigned long key;
  REAL dxcur=1./(1<<(level));

  unsigned int x=(unsigned int)(pos[0]/dxcur);
  unsigned int y=(unsigned int)(pos[1]/dxcur);
  unsigned int z=(unsigned int)(pos[2]/dxcur);
  LC2M(&key,x,y,z,level);

  return key;
}




// ================================================================================
void reorg(struct CPU *cpu, struct PARAM *param){

  int level;
  qsort(cpu->grid,param->ngridmax,sizeof(struct CELL),comp); // sorting the keys

  // lookup table for the first cell/levell
  for(level=0;level<=param->lmax;level++){
    unsigned long imin=0;
    unsigned long imax=cpu->ncelltotal;
    unsigned long target;
    LC2M(&target,0,0,0,level);
    cpu->firstcell[level]=floorkey(cpu->grid,target,imin,imax);
  }
}


#ifdef PIC
// ================================================================================
void reorgpart(struct CPU *cpu, struct PARAM *param){

  int level;
  qsort(cpu->part,param->npartmax,sizeof(struct PART),comppart); // sorting the keys

  // lookup table for the first part/level
#if 1
  for(level=param->lcoarse;level<=param->lmax;level++){
    unsigned long imin=0;
    unsigned long imax=cpu->nparttotal;
    unsigned long target;
    LC2M(&target,0,0,0,level);
    unsigned long fp=floorkeypart(cpu->part,target,imin,imax);

    cpu->firstpart[level]=fp;
    
    /* while((cpu->part[fp].key>=target)&&(fp>0)){  */
    /*   fp--;  */
    /* }  */

    /* if(fp==0){ */
    /*   cpu->firstpart[level]=0; */
    /* } */
    /* else{ */
    /*   cpu->firstpart[level]=fp+1; */
    /* } */
      
  }
#endif
}
#endif

// ================================================================================
struct CELL * getcell(unsigned long *key,int level,struct CPU *cpu){
  struct CELL *cell;
  cell=(struct CELL *)bsearch(key,cpu->grid+cpu->firstcell[level],cpu->ncell[level],sizeof(struct CELL),comp);
  return cell;
}


#ifdef PIC
// ================================================================================
struct PART * getpart(unsigned long *key,int level,struct CPU *cpu){
  struct PART *part;
  part=(struct PART *)bsearch(key,cpu->part+cpu->firstpart[level],cpu->npart[level],sizeof(struct PART),comppart);

  // rewind in case of multiple occurences
  if((part)&&(part>cpu->part)){
    while((part-1)->key==*key){
      part--;
      if(part==cpu->part) break;
    }
  }
  return part;
}
#endif


#ifdef PIC
//====================================================================================
//====================================================================================

void amr_update_key_part(unsigned int level, struct CPU *cpu, struct PARAM *param)
{
  // post refinement/derefinement key regularisation

  struct PART *part=&(cpu->part[cpu->firstpart[level]]);
  unsigned long nidx=cpu->npart[level];
  int nplp1=0;
  int nplm1=0;
  int idx;
  for(idx=0;idx<nidx;idx++){
    if(part->key<part->newkey){
      part->key=part->newkey;
      nplp1++;
    }
    else if(part->key>part->newkey){
      part->key=part->newkey;
      nplm1++;
    }
    part++;
  }
  cpu->npart[level+1]+=nplp1;
  cpu->npart[level-1]+=nplm1;
  cpu->npart[level]=cpu->npart[level]-(nplm1+nplp1);

  printf("AMR SHUFFLE adding %d part to level %d and removing %d part from level %d\n",nplp1,level+1,nplm1,level-1);
}


//====================================================================================
void update_key_part(unsigned int level, struct CPU *cpu, struct PARAM *param)
{
  // post refinement/derefinement key regularisation

  struct PART *part=&(cpu->part[cpu->firstpart[level]]);
  unsigned long nidx=cpu->npart[level];
  int nplp1=0;
  int nplm1=0;
  int idx;
  for(idx=0;idx<nidx;idx++){

    if(part->key!=part->newkey){
      // we should check if the cell with newkey exists !
      struct CELL *lcell;
      lcell=getcell(&part->newkey,level,cpu);
      if(lcell==NULL){
	// we should go LOWRES because it does not exist
	part->newkey=pos2key(part->x,level-1);
	nplm1++;
      }
      else{
	if(lcell->child){
	  // it exists it's refined, let's compute the key again
	  part->newkey=pos2key(part->x,level+1);
	  nplp1++;
	}
      }

      part->key=part->newkey;
      

    }

    part++;
  }

  cpu->npart[level+1]+=nplp1;
  cpu->npart[level-1]+=nplm1;
  cpu->npart[level]=cpu->npart[level]-nplm1-nplp1;

  printf("MOTION SHUFFLE adding %d part to level %d and adding %d part to level %d\n",nplp1,level+1,nplm1,level-1);


}

#endif


// ==============================================================================

void nei6(struct CELL *cell, unsigned long neikey[]){
  unsigned int i;
  unsigned int j;
  unsigned int k;
  unsigned int level;

  M2C(cell->key,&i,&j,&k,&level); // get the position

  unsigned int dmax=(1<<level)-1; //2**level-1
  unsigned int d;

  d=((int)i-1<0?dmax:i-1);
  LC2M(neikey,d,j,k,level);

  d=(i+1>dmax?0:i+1);
  LC2M(neikey+1,d,j,k,level);

  d=((int)j-1<0?dmax:j-1);
  LC2M(neikey+2,i,d,k,level);

  d=(j+1>dmax?0:j+1);
  LC2M(neikey+3,i,d,k,level);

  d=((int)k-1<0?dmax:k-1);
  LC2M(neikey+4,i,j,d,level);

  d=(k+1>dmax?0:k+1);
  LC2M(neikey+5,i,j,d,level);
}

// =============================================

void nei27(struct CELL *cell, unsigned long neikey[]){
  unsigned int i;
  unsigned int j;
  unsigned int k;
  unsigned int level;

  M2C(cell->key,&i,&j,&k,&level); // get the position

  unsigned int dmax=(1<<level)-1; //2**level-1
  int di;
  int dj;
  int dk;

  int count=0;
  for(dk=-1;dk<=1;dk++)
    {
      for(dj=-1;dj<=1;dj++)
	{
	  for(di=-1;di<=1;di++)
	    {
	      unsigned int ii;
	      unsigned int jj;
	      unsigned int kk;

	      ii=((int)i+di<0?dmax:((int)i+di>dmax?0:i+di));
	      jj=((int)j+dj<0?dmax:((int)j+dj>dmax?0:j+dj));
	      kk=((int)k+dk<0?dmax:((int)k+dk>dmax?0:k+dk));
	      LC2M(neikey+count,ii,jj,kk,level);
	      count++;
	    }
	}
    }
}


