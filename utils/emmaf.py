import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
import subprocess
import os
from scipy.integrate import simps, quad

#==============
class Raw:
#==============
    def load(self,fname):
        f=open(fname,"rb")
        self.ncell=np.fromfile(f,dtype=np.int32,count=1)[0]
        self.time=np.fromfile(f,dtype=np.float32,count=1)[0]
        self.limits=np.fromfile(f,dtype=np.float32,count=6)
        self.data=np.fromfile(f,dtype=np.float32,count=self.ncell)
        f.close
          
#==============
    def __init__(self):
        self.ncell=[]
        self.time=[]
        self.limits=[]
        self.data=[]


#==============
class Grid:
    def load(self,fname):
        f=open(fname,"rb")
        self.ncell=np.fromfile(f,dtype=np.int32,count=3)
        ncelltot=self.ncell[0]*self.ncell[1]*self.ncell[2]
        self.time=np.fromfile(f,dtype=np.float32,count=1)[0]
        self.limits=np.fromfile(f,dtype=np.float32,count=6)
        self.data=np.fromfile(f,dtype=np.float32,count=ncelltot)
        f.close
        self.data=np.reshape(self.data,(self.ncell[0],self.ncell[1],self.ncell[2]),order='F')
    
    def __init__(self):
        self.ncell=[]
        self.time=[]
        self.limits=[]
        self.data=[]
        
        

#==============
class Rawpart:
    def load(self,fname):
        f=open(fname,"rb")
        self.npart=np.fromfile(f,dtype=np.int32,count=1)[0]
        self.time=np.fromfile(f,dtype=np.float32,count=1)[0]
        plims=np.fromfile(f,dtype=np.float32,count=6)[0]
        self.data=np.fromfile(f,dtype=np.float32,count=self.npart)
        f.close
                              
    def __init__(self):
        self.time=[]
        self.npart=[]
        self.data=[]

#=========
class Part:
    def load(self,directory,isnap,star=False):
        repsnap=directory+'{:05d}/'.format(isnap)
        if(star):
            repx=repsnap+'star_x'
            repy=repsnap+'star_y'
            repz=repsnap+'star_z'
            repm=repsnap+'star_mass'
            repa=repsnap+'star_age'
            
        else:
            repx=repsnap+'part_x'
            repy=repsnap+'part_y'
            repz=repsnap+'part_z'
            repm=repsnap+'part_mass'
            
        ncpu=len(os.listdir(repx))
        
        raw=Rawpart()
        npart=0
        for icpu in range(ncpu):
           
            fname=repx+'/x.{:05d}.p{:05d}'.format(isnap,icpu)
            f=open(fname,"rb")
            npartloc=np.fromfile(f,dtype=np.int32,count=1)[0]
            npart+=npartloc
        
        print('Found {:d} files in {!s} and {:d} particles'.format(ncpu,repx,npart))
        
      
        self.x=np.zeros(npart,dtype=np.float32)
        self.y=np.zeros(npart,dtype=np.float32)
        self.z=np.zeros(npart,dtype=np.float32)
        self.mass=np.zeros(npart,dtype=np.float32)
        self.age=np.zeros(npart,dtype=np.float32)
        self.npart=npart
        offset=0
        raw=Rawpart()
        for icpu in range(ncpu):
            fname=repx+'/x.{:05d}.p{:05d}'.format(isnap,icpu)
            raw.load(fname)
            self.x[offset:offset+raw.npart]=raw.data
            
            fname=repy+'/y.{:05d}.p{:05d}'.format(isnap,icpu)
            raw.load(fname)
            self.y[offset:offset+raw.npart]=raw.data
            
            fname=repz+'/z.{:05d}.p{:05d}'.format(isnap,icpu)
            raw.load(fname)
            self.z[offset:offset+raw.npart]=raw.data
            
            fname=repm+'/mass.{:05d}.p{:05d}'.format(isnap,icpu)
            raw.load(fname)
            self.mass[offset:offset+raw.npart]=raw.data
            
            if(star):
                fname=repa+'/age.{:05d}.p{:05d}'.format(isnap,icpu)
                raw.load(fname)
                self.age[offset:offset+raw.npart]=raw.data
            
            offset=offset+raw.npart
        
        self.time=raw.time
        
    def loadtag(self,fname):
        print('TAG LOAD reading '+fname)
        f=open(fname,"rb")
        dummy=np.fromfile(f,dtype=np.int32,count=5)
        self.tag=np.fromfile(f,dtype=np.int32,count=self.npart)
        f.close
        
    def loadden(self,fname):
        print('TAG LOAD reading '+fname)
        f=open(fname,"rb")
        dummy=np.fromfile(f,dtype=np.int32,count=1)
        self.den=np.fromfile(f,dtype=np.int32,count=self.npart)
        f.close

    
    def __init__(self):
        self.time=[]
        self.x=[]
        self.y=[]
        self.z=[]
        self.tag=[]
        self.age=[]
        self.mass=[]
        self.npart=[]
        self.den=[]

        
#================
class Cell:
    def load(self,directory,fieldname,isnap):
        repsnap=directory+'{:05d}/'.format(isnap)
        
        repx=repsnap+'grid_x'
        repy=repsnap+'grid_y'
        repz=repsnap+'grid_z'
        repl=repsnap+'grid_l'
        repf=repsnap+'grid_'+fieldname
            
        ncpu=len(os.listdir(repx))
        
        ncell=0
        for icpu in range(ncpu):
            fname=repx+'/x.{:05d}.p{:05d}'.format(isnap,icpu)
            f=open(fname,"rb")
            ncellloc=np.fromfile(f,dtype=np.int32,count=1)[0]
            ncell+=ncellloc
        
        print('Found {:d} files in {!s} and {:d} cells'.format(ncpu,repx,ncell))
        
      
        self.x=np.zeros(ncell,dtype=np.float32)
        self.y=np.zeros(ncell,dtype=np.float32)
        self.z=np.zeros(ncell,dtype=np.float32)
        self.level=np.zeros(ncell,dtype=np.float32)
        self.field=np.zeros(ncell,dtype=np.float32)
        self.ncell=ncell
        offset=0
        raw=Raw()
        for icpu in range(ncpu):
            fname=repx+'/x.{:05d}.p{:05d}'.format(isnap,icpu)
            raw.load(fname)
            self.x[offset:offset+raw.ncell]=raw.data
            
            fname=repy+'/y.{:05d}.p{:05d}'.format(isnap,icpu)
            raw.load(fname)
            self.y[offset:offset+raw.ncell]=raw.data
            
            fname=repz+'/z.{:05d}.p{:05d}'.format(isnap,icpu)
            raw.load(fname)
            self.z[offset:offset+raw.ncell]=raw.data
            
            fname=repl+'/l.{:05d}.p{:05d}'.format(isnap,icpu)
            raw.load(fname)
            self.level[offset:offset+raw.ncell]=raw.data
            
            fname=repf+'/'+fieldname+'.{:05d}.p{:05d}'.format(isnap,icpu)
            raw.load(fname)
            self.field[offset:offset+raw.ncell]=raw.data
            
            offset=offset+raw.ncell
        
        self.time=raw.time
     
    
    def __init__(self):
        self.time=[]
        self.x=[]
        self.y=[]
        self.z=[]
        self.level=[]
        self.data=[]
     


#==============

def getR200(rc,m,delta=200):
    s=np.argsort(rc)
    rc=rc[s]
    m=m[s]
    i=np.arange(rc.size)*m
    dens=i/(4./3.*np.pi*pow(rc+1e-6,3))
    ww=np.where(dens<=200)
    if ww[0].size>1:
        res=rc[ww[0][1]]
    else:
        res=np.zeros(1)

    return res


def getcdm(p,w):
    xc=np.average(p.x[w])
    yc=np.average(p.y[w])
    zc=np.average(p.z[w])
    return xc,yc,zc

def getdensest(p,w):
    iloc=p.den[w].argmax()
    xc=p.x[w][iloc]
    yc=p.y[w][iloc]
    zc=p.z[w][iloc]
    return xc,yc,zc

def getdist2cdm(p,w):
    xc,yc,zc=getcdm(p,w)
    rc=np.sqrt(pow(p.x[w]-xc,2)+pow(p.y[w]-yc,2)+pow(p.z[w]-zc,2))
    return rc

def radselect(p,x,y,z,r):
    x2=np.where(p.x-x>0.5,p.x-x-1.,px-x)
    y2=np.where(p.y-y>0.5,p.y-y-1.,py-y)
    z2=np.where(p.z-z>0.5,p.z-z-1.,pz-z)
    r2=np.sqrt(x2**2+y2**2+z2**2)
    return np.where(r2<r)


def field2grid(dirname,fieldname,cen,dcen,level,execname='./field2grid'):
    
    commandeformat=execname+' {!s} {!s} {:f} {:f} {:f} {:f} {:f} {:f} {:d}'
    commande=commandeformat.format(dirname,fieldname,cen[0],cen[1],cen[2],dcen[0],dcen[1],dcen[2],level)
    subprocess.check_call(commande,shell=True)
    field=Grid()
    field.load('dump.dat')
    return field
    
# ========================================================
def dtda(a,H0=70,Om=0.3,Ov=0.7):
    H0=H0*1e3/3.08567758e22 #s-1
    res=(1./H0)/np.sqrt(Om/a+Ov*a*a)
    return res

def a2t(a,H0=70,Om=0.3,Ov=0.7):
    res=quad(lambda x: dtda(x,H0=H0,Om=Om,Ov=Ov), 0, a)
    return res[0]   
    

    