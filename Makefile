#=======================================================================
# Makefile for Xcompact3D
#=======================================================================
# Choose pre-processing options
#   -DDOUBLE_PREC - use double-precision
#   -DSAVE_SINGLE - Save 3D data in single-precision
#   -DDEBG        - debuggin xcompact3d.f90
# generate a Git version string
GIT_VERSION := $(shell git describe --tag --long --always)

DEFS = -DOVERWRITE -DSAVE_SINGLE -DDOUBLE_PREC -DVERSION=\"$(GIT_VERSION)\"

LCL = local# local,lad,sdu,archer
IVER = 17# 15,16,17,18
CMP = gcc# intel,gcc
FFT = generic# generic,fftw3,mkl

#######CMP settings###########
ifeq ($(CMP),intel)
CC = mpiicc
FC = mpiifort
#FFLAGS = -fpp -O3 -xHost -heap-arrays -shared-intel -mcmodel=large -safe-cray-ptr -g -traceback
FFLAGS = -fpp -O3 -xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2 -ipo -fp-model fast=2 -mcmodel=large -safe-cray-ptr -I$(MPI_ROOT)/lib
##debuggin test: -check all -check bounds -chintel eck uninit -gen-interfaces -warn interfaces
else ifeq ($(CMP),gcc)
CC = mpicc
FC = mpif90
#FFLAGS = -O3 -funroll-loops -floop-optimize -g -Warray-bounds -fcray-pointer -x f95-cpp-input
FFLAGS = -cpp -O0 -funroll-loops -floop-optimize -g -fno-inline -finstrument-functions -Warray-bounds -fcray-pointer -fbacktrace -ffree-line-length-none -fallow-argument-mismatch -msse4.1 -msse4.2 -mavx2
#-ffpe-trap=invalid,zero
else ifeq ($(CMP),nagfor)
FC = mpinagfor
FFLAGS = -fpp
else ifeq ($(CMP),cray)
FC = ftn
FFLAGS = -eF -g -O3 -N 1023
endif


MODDIR = ./mod
DECOMPDIR = ./decomp2d
SRCDIR = ./src
TURBDIR = ./src

### List of files for the main code
SRCSIG = $(SRCDIR)/signal.c
OBJSIG = $(SRCSIG:%.c=%.o)
SRCDECOMP = $(DECOMPDIR)/decomp_2d.f90 $(DECOMPDIR)/glassman.f90 $(DECOMPDIR)/fft_$(FFT).f90 $(DECOMPDIR)/io.f90
OBJDECOMP = $(SRCDECOMP:%.f90=%.o)
SRC = $(SRCDIR)/module_param.f90 $(SRCDIR)/variables.f90 $(SRCDIR)/poisson.f90 $(SRCDIR)/derive.f90 $(SRCDIR)/implicit.f90 $(SRCDIR)/schemes.f90 $(SRCDIR)/parameters.f90 $(SRCDIR)/*.f90
OBJ = $(SRC:%.f90=%.o)
SRC = $(SRCDIR)/module_param.f90 $(SRCDIR)/variables.f90 $(SRCDIR)/poisson.f90 $(SRCDIR)/ibm.f90 $(SRCDIR)/derive.f90 $(SRCDIR)/implicit.f90 $(SRCDIR)/schemes.f90 $(SRCDIR)/forces.f90 $(SRCDIR)/probes.f90 $(SRCDIR)/navier.f90 $(SRCDIR)/tools.f90 $(SRCDIR)/visu.f90 $(SRCDIR)/BC-TBL.f90 $(SRCDIR)/BC-ABL.f90 $(SRCDIR)/les_models.f90 $(SRCDIR)/BC-Lock-exchange.f90 $(SRCDIR)/time_integrators.f90 $(SRCDIR)/filters.f90 $(SRCDIR)/parameters.f90 $(SRCDIR)/BC-User.f90 $(SRCDIR)/BC-TGV.f90 $(SRCDIR)/BC-Channel-flow.f90 $(SRCDIR)/BC-Cavity.f90 $(SRCDIR)/BC-Periodic-hill.f90 $(SRCDIR)/BC-Cylinder.f90 $(SRCDIR)/BC-Mixing-layer.f90 $(SRCDIR)/BC-Jet.f90 $(SRCDIR)/BC-dbg-schemes.f90 $(SRCDIR)/BC-Uniform.f90 $(TURBDIR)/constants.f90 $(TURBDIR)/acl_utils.f90 $(TURBDIR)/airfoils.f90 $(TURBDIR)/dynstall.f90 $(TURBDIR)/dynstall_legacy.f90 $(TURBDIR)/acl_elem.f90 $(TURBDIR)/acl_controller.f90 $(TURBDIR)/acl_turb.f90 $(TURBDIR)/acl_out.f90 $(TURBDIR)/acl_farm_controller.f90 $(TURBDIR)/acl_model.f90 $(TURBDIR)/acl_source.f90 $(TURBDIR)/adm.f90 $(TURBDIR)/turbine.f90 $(SRCDIR)/statistics.f90 $(SRCDIR)/case.f90 $(SRCDIR)/transeq.f90 $(SRCDIR)/genepsi3d.f90 $(SRCDIR)/xcompact3d.f90


#######FFT settings##########
ifeq ($(FFT),fftw3)
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
else ifeq ($(FFT),mkl)
  SRCDECOMP := $(DECOMPDIR)/mkl_dfti.f90 $(SRCDECOMP)
  LIBFFT=-Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread
	INC=-I$(MKLROOT)/include
endif

#######OPTIONS settings###########
OPT = -I$(SRCDIR) -I$(DECOMPDIR) $(FFLAGS)
LINKOPT = $(FFLAGS)
#-----------------------------------------------------------------------
# Normally no need to change anything below

all: xcompact3d

xcompact3d : $(OBJDECOMP) $(OBJ) $(OBJSIG)
	$(FC) -o $@ $(LINKOPT) $(OBJDECOMP) $(OBJ) $(OBJSIG) $(LIBFFT)

$(OBJDECOMP):$(DECOMPDIR)%.o : $(DECOMPDIR)%.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(DEFS2) $(INC) -c $<
	mv $(@F) ${DECOMPDIR}

$(OBJ):$(SRCDIR)%.o : $(SRCDIR)%.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(DEFS2) $(INC) -c $<
	mv $(@F) ${SRCDIR}

$(OBJSIG):$(SRCDIR)%.o : $(SRCDIR)%.c
	$(CC) $(DEFS) $(DEFS2) -c $<
	mv $(@F) ${SRCDIR}

#
# Explicit dependencies allow parallel make
#
# Decomp 2D
#
$(DECOMPDIR)/glassman.o: $(DECOMPDIR)/decomp_2d.o
$(DECOMPDIR)/io.o: $(DECOMPDIR)/decomp_2d.o
$(DECOMPDIR)/fft_$(FFT).o: $(DECOMPDIR)/glassman.o
#
# X3D
#
$(SRCDIR)/acl_controller.o: $(DECOMPDIR)/decomp_2d.o
$(SRCDIR)/acl_elem.o: $(SRCDIR)/dynstall_legacy.o
$(SRCDIR)/acl_farm_controller.o: $(SRCDIR)/acl_turb.o
$(SRCDIR)/acl_model.o: $(SRCDIR)/acl_elem.o $(SRCDIR)/acl_out.o
$(SRCDIR)/acl_out.o: $(SRCDIR)/acl_turb.o
$(SRCDIR)/acl_source.o: $(SRCDIR)/acl_model.o $(SRCDIR)/acl_utils.o
$(SRCDIR)/acl_turb.o: $(SRCDIR)/module_param.o $(SRCDIR)/acl_elem.o
$(SRCDIR)/acl_utils.o: $(DECOMPDIR)/decomp_2d.o $(SRCDIR)/constants.o
$(SRCDIR)/adm.o: $(SRCDIR)/airfoils.o $(SRCDIR)/variables.o
$(SRCDIR)/airfoils.o: $(SRCDIR)/acl_utils.o
$(SRCDIR)/BC-ABL.o: $(SRCDIR)/visu.o
$(SRCDIR)/BC-Channel-flow.o: $(SRCDIR)/visu.o
$(SRCDIR)/BC-Cylinder.o: $(SRCDIR)/visu.o
$(SRCDIR)/BC-dbg-schemes.o: $(SRCDIR)/visu.o
$(SRCDIR)/BC-Jet.o: $(SRCDIR)/visu.o
$(SRCDIR)/BC-Lock-exchange.o: $(SRCDIR)/visu.o
$(SRCDIR)/BC-Mixing-layer.o: $(SRCDIR)/visu.o
$(SRCDIR)/BC-Periodic-hill.o: $(SRCDIR)/visu.o
$(SRCDIR)/BC-TBL.o: $(SRCDIR)/visu.o
$(SRCDIR)/BC-TGV.o: $(SRCDIR)/visu.o
$(SRCDIR)/BC-Uniform.o: $(SRCDIR)/visu.o
$(SRCDIR)/BC-User.o: $(SRCDIR)/visu.o
$(SRCDIR)/BC-Cavity.o: $(SRCDIR)/visu.o
$(SRCDIR)/case.o: $(SRCDIR)/BC-Lock-exchange.o $(SRCDIR)/BC-Channel-flow.o $(SRCDIR)/acl_turb.o $(SRCDIR)/visu.o $(SRCDIR)/BC-Cavity.o
$(SRCDIR)/constants.o: $(DECOMPDIR)/decomp_2d.o
$(SRCDIR)/derive.o: $(SRCDIR)/ibm.o
$(SRCDIR)/dynstall.o: $(SRCDIR)/airfoils.o
$(SRCDIR)/dynstall_legacy.o: $(SRCDIR)/airfoils.o
$(SRCDIR)/filters.o: $(SRCDIR)/ibm.o
$(SRCDIR)/forces.o: $(SRCDIR)/module_param.o $(DECOMPDIR)/io.o
$(SRCDIR)/genepsi3d.o: $(SRCDIR)/BC-Cylinder.o $(SRCDIR)/BC-Periodic-hill.o
$(SRCDIR)/ibm.o: $(SRCDIR)/module_param.o
$(SRCDIR)/implicit.o: $(SRCDIR)/variables.o
$(SRCDIR)/les_models.o: $(SRCDIR)/BC-ABL.o $(SRCDIR)/tools.o
$(SRCDIR)/module_param.o: $(DECOMPDIR)/decomp_2d.o
$(SRCDIR)/navier.o: $(SRCDIR)/forces.o $(SRCDIR)/poisson.o
$(SRCDIR)/parameters.o: $(SRCDIR)/BC-Lock-exchange.o $(SRCDIR)/visu.o
$(SRCDIR)/poisson.o: $(SRCDIR)/variables.o
$(SRCDIR)/probes.o: $(SRCDIR)/variables.o
$(SRCDIR)/schemes.o: $(SRCDIR)/implicit.o
$(SRCDIR)/statistics.o: $(SRCDIR)/variables.o $(SRCDIR)/tools.o
$(SRCDIR)/time_integrators.o: $(SRCDIR)/implicit.o $(SRCDIR)/navier.o
$(SRCDIR)/tools.o: $(SRCDIR)/navier.o
$(SRCDIR)/transeq.o: $(SRCDIR)/case.o $(SRCDIR)/les_models.o $(SRCDIR)/probes.o
$(SRCDIR)/turbine.o: $(SRCDIR)/acl_model.o $(SRCDIR)/acl_source.o $(SRCDIR)/adm.o
$(SRCDIR)/variables.o: $(SRCDIR)/module_param.o
$(SRCDIR)/visu.o: $(SRCDIR)/tools.o
$(SRCDIR)/xcompact3d.o: $(SRCDIR)/transeq.o

.PHONY: post
post:
	$(FC) $(FFLAGS) $(DEFS) $(DEFS2) post.f90 -c
	$(FC) $(FFLAGS) -o $@ $(PSRC:.f90=.o)

.PHONY: clean


clean:
	rm -f $(DECOMPDIR)/*.o $(DECOMPDIR)/*.mod
	rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod
	rm -f *.o *.mod xcompact3d post

.PHONY: cleanall
cleanall: clean
	rm -f *~ \#*\# out/* data/* stats/* planes/* *.xdmf *.log *.out nodefile core sauve*
