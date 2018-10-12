#=======================================================================
# Makefile for Incompact3D modified by Ricardo
#=======================================================================

# Choose pre-processing options
#   -DDOUBLE_PREC - use double-precision
#   -DSAVE_SINGLE - Save 3D data in single-precision
#   -DDEBG        - debuggin incompact3d.f90
#   -DIBM         - enable IBM calls
#   -DPOST        - enable statistics processing
#   -DVISU        - enable visu.f90
#   -DVISUEXTRA   - enable extra options visu.f90
#   -DFORCES      - enable lift and drag computing over solid body
#   -DELES        - enable explicit LES modelling
#   -DSTRETCHING  - enable mesh stretching in y direction
# generate a Git version string
GIT_VERSION := $(shell git describe --tag --long --always)

#FLOW_TYPE = #Lock-exchange# TGV# Channel-flow# Cylinder

DEFS = -DVISU -DVISUEXTRA -DDOUBLE_PREC -DVERSION=\"$(GIT_VERSION)\"

LCL = local# local,lad,sdu,archer
IVER = 17# 15,16,17,18
CMP = gcc# intel,gcc
FFT = generic# mkl,generic,fftw3

#######Minimum defs###########
ifeq ($(FLOW_TYPE),Channel-flow)
DEFS2 = -DSTRETCHING -DPOST
else ifeq ($(FLOW_TYPE),Cylinder)
DEFS2 = -DIBM -DFORCES
else ifeq ($(FLOW_TYPE),Lock-exchange)
DEFS2 = -DPOST
else ifeq ($(FLOW_TYPE),Periodic-hill)
DEFS2 = -DIBM -DSTRETCHING -DPOST
else ifeq ($(FLOW_TYPE),TGV)
DEFS2 = -DPOST
endif

#######CMP settings###########
ifeq ($(CMP),intel)
FC = mpiifort
FFLAGS = -fpp -O3 -xHost -heap-arrays -shared-intel -mcmodel=large -safe-cray-ptr -g -traceback
##debuggin test: -check all -check bounds -chintel eck uninit -gen-interfaces -warn interfaces
else ifeq ($(CMP),gcc)
FC = mpif90
#FFLAGS = -O3 -funroll-loops -floop-optimize -g -Warray-bounds -fcray-pointer -x f95-cpp-input
FFLAGS = -cpp -O3 -funroll-loops -floop-optimize -g -Warray-bounds -fcray-pointer -fbacktrace -march=native -ffree-line-length-none
endif

### List of files for the main code
SRC = decomp_2d.f90 glassman.f90 fft_$(FFT).f90 module_param.f90 io.f90 variables.f90 poisson.f90 schemes.f90 BC-$(FLOW_TYPE).f90 implicit.f90 convdiff.f90 navier.f90 derive.f90 parameters.f90 tools.f90 visu.f90 paraview.f90 genepsi3d.f90 filter.f90 les_models.f90 incompact3d.f90

### List of files for the post-processing code
PSRC = decomp_2d.f90 module_param.f90 io.f90 variables.f90 schemes.f90 derive.f90 BC-$(FLOW_TYPE).f90 parameters.f90 tools.f90 visu.f90 paraview.f90 post.f90

######MKL INSTALL PATH######
ifeq ($(LCL),local)
  MKLROOT=/opt/intel/mkl

else ifeq ($(LCL),lad)
ifeq ($(IVER),17)
  MKLROOT=/usr/local/Intel_Cluster_Studio_XE_2017/parallel_studio_xe_2017/mkl
else ifeq ($(IVER),16)
  MKLROOT=/usr/local/Intel_Cluster_Studio_XE_2016/parallel_studio_xe_2016_update3/mkl
else ifeq ($(IVER),13)
  MKLROOT=/usr/local/Intel_Cluster_Studio_XE_2013/parallel_studio_xe_2013/l_ics_2013.0.028/mkl
endif
else ifeq ($(LCL),sdu)
ifeq ($(IVER),17)
  MKLROOT=/opt/intel/parallel_studio_xe_2017/mkl
else ifeq ($(IVER),16)
  MKLROOT=/opt/intel/parallel_studio_xe_2016/mkl
endif
endif

#######FFT settings##########
ifeq ($(FFT),mkl)
  SRC := mkl_dfti.f90 $(SRC)
  LIBFFT=-Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread
  INC=-I$(MKLROOT)/include
  MKL_MOD=mkl_mod
  MKL_DFTI=mkl_dfti
else ifeq ($(FFT),fftw3)
  #FFTW3_PATH=/usr 
  #FFTW3_PATH=/usr/lib64
  FFTW3_PATH=/usr/local/Cellar/fftw/3.3.7_1
  INC=-I$(FFTW3_PATH)/include
  LIBFFT=-L$(FFTW3_PATH) -lfftw3 -lfftw3f
else ifeq ($(FFT),fftw3_f03)
  FFTW3_PATH=/usr                                #ubuntu # apt install libfftw3-dev
  #FFTW3_PATH=/usr/lib64                         #fedora # dnf install fftw fftw-devel
  #FFTW3_PATH=/usr/local/Cellar/fftw/3.3.7_1     #macOS  # brew install fftw
  INC=-I$(FFTW3_PATH)/include
  LIBFFT=-L$(FFTW3_PATH)/lib -lfftw3 -lfftw3f
else ifeq ($(FFT),generic)
  INC=
  LIBFFT=
endif

#######OPTIONS settings###########
#ifneq (,$(findstring DSHM,$(OPTIONS)))
#SRC := FreeIPC.f90 $(SRC)
#OBJ = $(SRC:.f90=.o) alloc_shm.o FreeIPC_c.o
#else
OBJ = $(SRC:.f90=.o)
#endif

#-----------------------------------------------------------------------
# Normally no need to change anything below

all: incompact3d

#alloc_shm.o: alloc_shm.c
#	$(CC) $(CFLAGS) -c $<

#FreeIPC_c.o: FreeIPC_c.c
#	$(CC) $(CFLAGS) -c $<

incompact3d : $(OBJ)
	$(FC) -O3 -o $@ $(OBJ) $(LIBFFT)

%.o : %.f90
	$(FC) $(FFLAGS) $(DEFS) $(DEFS2) $(INC) -c $<

.PHONY: post
post:
	$(FC) $(FFLAGS) $(DEFS) $(DEFS2) post.f90 -c
	$(FC) $(FFLAGS) -o $@ $(PSRC:.f90=.o)

.PHONY: clean
clean:
	rm -f *.o *.mod incompact3d post
.PHONY: cleanall
cleanall: clean
	rm -f *~ \#*\# out/* data/* stats/* planes/* *.xdmf *.log *.out nodefile core sauve*
