#===============================================================================
# Architecture parameter
#===============================================================================
# ARCH = GPU

#===============================================================================
# Code parameters
#===============================================================================
DEFINES =
include param.mk

ifeq ($(ARCH),GPU)
	DEFINES += -DGPUAXL
endif

#===============================================================================
# Compilers
#===============================================================================

# main compiler
CC= mpicc
C_FLAGS=  -g -O2 -Wall -ftree-vectorize -ffast-math -fno-cx-limited-range -Wimplicit
C_LIBS= -lm #-fopenmp -lstdc++

#GPU compiler
NVCC= nvcc
NVCC_FLAGS= #-lstdc++ --ptxas-options=-v -arch=sm_35 #-g -G #--device-emulation
NVCC_LIBS=

#===============================================================================
# External libs
#===============================================================================

# HDF5
#I_DIR += -I/usr/include/hdf5/mpich
#LD_DIR += -L/usr/lib/x86_64-linux-gnu/hdf5/mpich
#LD_FLAGS += -lhdf5

# CUDA
# I_DIR +=  -I/workdir/observatoire/aubert/cudpp_src_2.0/include -I/usr/local/cuda-5.0/include
# LD_DIR += -L/workdir/observatoire/aubert/cudpp_src_2.0/lib -L/usr/local/cuda-5.0/lib64 -L/usr/lib/openmpi/lib/
# LD_FLAGS += -lcudart -lmpi -lopen-rte -lopen-pal -ldl -lnsl -lutil -ldl -lcudpp #-lcuda

I_DIR += -I/usr/include/
# LD_DIR += -L/
#LD_FLAGS += -lcuda

#MPI
I_DIR += -I/usr/include/mpi
LD_DIR += -L/usr/lib
LD_FLAGS += -lmpi

#===============================================================================
# Obj
#===============================================================================

C_OBJS= emmma.o allocation.o morton.o io.o amr.o poisson.o cic.o parameters.o ic.o friedmann.o advance.o particle.o communication.o

CUDA_OBJS= \
	interface.o \
	poisson_utils_gpu.o \
	hydro_utils_gpu.o \
	rad_utils_gpu.o \
	chem_utils_gpu.o \
	# cic_gpu.o

#===============================================================================
# Compile
#===============================================================================

OBJDIR = obj
SRCDIR = src

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(C_FLAGS) $(C_LIBS) $(DEFINES) $(I_DIR) -c $< -o $@ 

ifeq ($(ARCH),GPU)
$(OBJDIR)/%.o: $(SRCDIR)/%.cu
	$(NVCC) $(NVCC_FLAGS) $(DEFINES) $(I_DIR) -c $< -o $@
endif

#===============================================================================
# Link
#===============================================================================

ifeq ($(ARCH),GPU)
EXECUTABLE = emmmagpu
else
EXECUTABLE = emmmacpu
endif

OBJ=$(patsubst %,$(OBJDIR)/%,$(C_OBJS))
ifeq ($(ARCH),GPU)
OBJ += $(patsubst %,$(OBJDIR)/%,$(CUDA_OBJS))
endif

ifeq ($(ARCH),GPU)
LCC=$(NVCC)
else
LCC=$(CC)
endif

all:$(OBJ)
	$(LCC) $(OBJ) $(LD_DIR) $(C_FLAGS) $(C_LIBS) -o $(EXECUTABLE) $(LD_FLAGS) 

# oct2grid:
# 	$(CC) $(DEFINESGLOB) $(C_LIBS) $(C_FLAGS) utils/oct2grid.c -o utils/oct2grid -lm
# field2grid:
# 	$(CC) $(DEFINESGLOB) $(C_LIBS) $(C_FLAGS) utils/field2grid.c -o utils/field2grid -lm
# alloct:
# 	$(CC) $(DEFINESGLOB) $(C_LIBS) $(C_FLAGS) -o utils/alloct utils/alloct.c -lm
# part2cic:
# 	$(CC) $(DEFINESGLOB) $(C_LIBS) $(C_FLAGS) utils/part2cic.c -o utils/part2cic -lm
# cube2silo:
# 	$(CC) $(DEFINESGLOB) $(C_LIBS) $(C_FLAGS) -o utils/cube2silo utils/cube2silo.c utils/libsilo.a -lm
# part2silo:
# 	$(CC) $(DEFINESGLOB) $(C_LIBS) $(C_FLAGS) -o utils/part2silo utils/part2silo.c utils/libsilo.a -lm

clean:
	rm -f $(OBJDIR)/*.o *.cudafe1.* *.cudafe2.* *.hash *.ptx *fatbin.c *.cubin *.cpp* $(EXECUTABLE) *~
