#=======================================================================
# Makefile for Incompact3D modified by Ricardo
#=======================================================================
# Choose pre-processing options
#   -DDOUBLE_PREC - use double-precision
#   -DSAVE_SINGLE - Save 3D data in single-precision
#   -DDEBG        - debuggin incompact3d.f90
#   -DVISU        - enable visu.f90
#   -DVISUEXTRA   - enable extra options visu.f90
#   -DFORCES      - enable lift and drag computing over solid body
# generate a Git version string
GIT_VERSION := $(shell git describe --tag --long --always)

#######Select Flow Type#######
# FLOW_TYPE = Lock-exchange
# FLOW_TYPE = TGV
  FLOW_TYPE = Channel-flow
# FLOW_TYPE = Periodic-hill
# FLOW_TYPE = Cylinder
# FLOW_TYPE = dbg-schemes

DEFS = -DDOUBLE_PREC -DVERSION=\"$(GIT_VERSION)\"

LCL = local# local,lad,sdu,archer
IVER = 17# 15,16,17,18
CMP = gcc# intel,gcc
FFT = generic# mkl,generic,fftw3

#######Minimum defs###########
TWOD = 0
ifeq ($(FLOW_TYPE),Channel-flow)
DEFS2 = -DSTRETCHING -DPOST 
else ifeq ($(FLOW_TYPE),Cylinder)
DEFS2 = -DFORCES
else ifeq ($(FLOW_TYPE),Lock-exchange)
DEFS2 = -DPOST
else ifeq ($(FLOW_TYPE),Periodic-hill)
DEFS2 = -DSTRETCHING -DPOST
else ifeq ($(FLOW_TYPE),TGV)
DEFS2 = -DPOST
endif
ifeq ($(TWOD), 1)
DEFS2 += -DTWOD
endif
DEFS2 += -DVISU

#######CMP settings###########
ifeq ($(CMP),intel)
FC = mpiifort
#FFLAGS = -fpp -O3 -xHost -heap-arrays -shared-intel -mcmodel=large -safe-cray-ptr -g -traceback
FFLAGS = -fpp -O3 -xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2 -ipo -fp-model fast=2 -mcmodel=large -safe-cray-ptr
##debuggin test: -check all -check bounds -chintel eck uninit -gen-interfaces -warn interfaces
else ifeq ($(CMP),gcc)
FC = mpif90
#FFLAGS = -O3 -funroll-loops -floop-optimize -g -Warray-bounds -fcray-pointer -x f95-cpp-input
FFLAGS = -cpp  -funroll-loops -floop-optimize -g -Warray-bounds -fcray-pointer -fbacktrace -ffree-line-length-none
endif


MODDIR = ./mod  
DECOMPDIR = ./decomp2d
SRCDIR = ./src

### List of files for the main code
SRCDECOMP = $(DECOMPDIR)/decomp_2d.f90 $(DECOMPDIR)/glassman.f90 $(DECOMPDIR)/fft_$(FFT).f90 $(DECOMPDIR)/module_param.f90 $(DECOMPDIR)/io.f90 
OBJDECOMP = $(SRCDECOMP:%.f90=%.o)
SRC = $(SRCDIR)/variables.f90 $(SRCDIR)/poisson.f90 $(SRCDIR)/schemes.f90 $(SRCDIR)/derive.f90 $(SRCDIR)/parameters.f90 $(SRCDIR)/*.f90
OBJ = $(SRC:%.f90=%.o)
SRC = $(SRCDIR)/variables.f90 $(SRCDIR)/poisson.f90 $(SRCDIR)/schemes.f90 $(SRCDIR)/BC-TGV.f90 $(SRCDIR)/BC-Channel-flow.f90 $(SRCDIR)/BC-Cylinder.f90 $(SRCDIR)/BC-dbg-schemes.f90 $(SRCDIR)/case.f90 $(SRCDIR)/transeq.f90 $(SRCDIR)/forces.f90 $(SRCDIR)/navier.f90 $(SRCDIR)/derive.f90 $(SRCDIR)/filters.f90 $(SRCDIR)/parameters.f90 $(SRCDIR)/tools.f90 $(SRCDIR)/visu.f90 $(SRCDIR)/paraview.f90 $(SRCDIR)/genepsi3d.f90 $(SRCDIR)/les_models.f90 $(SRCDIR)/incompact3d.f90


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
OPT = -I$(SRCDIR) -I$(DECOMPDIR) $(FFLAGS)
LINKOPT = $(FFLAGS)
#-----------------------------------------------------------------------
# Normally no need to change anything below

all: incompact3d

incompact3d : $(OBJDECOMP) $(OBJ)
	$(FC) -o $@ $(LINKOPT) $(OBJDECOMP) $(OBJ) $(LIBFFT)

$(OBJDECOMP):$(DECOMPDIR)%.o : $(DECOMPDIR)%.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(DEFS2) $(INC) -c $<
	mv $(@F) ${DECOMPDIR}
	#mv *.mod ${DECOMPDIR}


$(OBJ):$(SRCDIR)%.o : $(SRCDIR)%.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(DEFS2) $(INC) -c $<
	mv $(@F) ${SRCDIR}
	#mv *.mod ${SRCDIR}

## This %.o : %.f90 doesn't appear to be called...
%.o : %.f90
	$(FC) $(FFLAGS) $(DEFS) $(DEFS2) $(INC) -c $<

.PHONY: post
post:
	$(FC) $(FFLAGS) $(DEFS) $(DEFS2) post.f90 -c
	$(FC) $(FFLAGS) -o $@ $(PSRC:.f90=.o)

.PHONY: clean

clean:
	rm -f $(DECOMPDIR)/*.o $(DECOMPDIR)/*.mod 
	rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod 
	rm -f *.o *.mod incompact3d post

.PHONY: cleanall
cleanall: clean
	rm -f *~ \#*\# out/* data/* stats/* planes/* *.xdmf *.log *.out nodefile core sauve*
