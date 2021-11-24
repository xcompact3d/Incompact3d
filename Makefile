#=======================================================================
# Makefile for Xcompact3D
#=======================================================================
# Choose pre-processing options
#   -DDOUBLE_PREC - use double-precision
#   -DSAVE_SINGLE - Save 3D data in single-precision
#   -DDEBG        - debuggin xcompact3d.f90
# generate a Git version string
GIT_VERSION := $(shell git describe --tag --long --always)

DEFS = -DDOUBLE_PREC -DVERSION=\"$(GIT_VERSION)\"

LCL = local# local,lad,sdu,archer
IVER = 17# 15,16,17,18
CMP = gcc# intel,gcc
FFT = generic# generic,fftw3,mkl

#######CMP settings###########
ifeq ($(CMP),intel)
FC = mpiifort
#FFLAGS = -fpp -O3 -xHost -heap-arrays -shared-intel -mcmodel=large -safe-cray-ptr -g -traceback
FFLAGS = -fpp -O3 -xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2 -ipo -fp-model fast=2 -mcmodel=large -safe-cray-ptr -I$(MPI_ROOT)/lib
##debuggin test: -check all -check bounds -chintel eck uninit -gen-interfaces -warn interfaces
else ifeq ($(CMP),gcc)
FC = mpif90
#FFLAGS = -O3 -funroll-loops -floop-optimize -g -Warray-bounds -fcray-pointer -x f95-cpp-input
FFLAGS = -cpp -O3 -funroll-loops -floop-optimize -g -Warray-bounds -fcray-pointer -fbacktrace -ffree-line-length-none
#-ffpe-trap=invalid,zero
ifeq "$(shell expr `gfortran -dumpversion | cut -f1 -d.` \>= 10)" "1"
FFLAGS += -fallow-argument-mismatch
endif
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
SRCDECOMP = $(DECOMPDIR)/decomp_2d.f90 $(DECOMPDIR)/glassman.f90 $(DECOMPDIR)/fft_$(FFT).f90 $(DECOMPDIR)/io.f90
OBJDECOMP = $(SRCDECOMP:%.f90=%.o)
SRC = $(SRCDIR)/module_param.f90 $(SRCDIR)/variables.f90 $(SRCDIR)/poisson.f90 $(SRCDIR)/derive.f90 $(SRCDIR)/implicit.f90 $(SRCDIR)/schemes.f90 $(SRCDIR)/parameters.f90 $(SRCDIR)/*.f90
OBJ = $(SRC:%.f90=%.o)
SRC = $(SRCDIR)/module_param.f90 $(SRCDIR)/variables.f90 $(SRCDIR)/BC-dbg-schemes.f90 $(SRCDIR)/poisson.f90 $(SRCDIR)/ibm.f90 $(SRCDIR)/derive.f90 $(SRCDIR)/implicit.f90 $(SRCDIR)/schemes.f90 $(SRCDIR)/forces.f90 $(SRCDIR)/probes.f90 $(SRCDIR)/navier.f90 $(SRCDIR)/tools.f90 $(SRCDIR)/visu.f90 $(SRCDIR)/BC-TBL.f90 $(SRCDIR)/BC-ABL.f90 $(SRCDIR)/les_models.f90 $(SRCDIR)/BC-Lock-exchange.f90 $(SRCDIR)/time_integrators.f90 $(SRCDIR)/filters.f90 $(SRCDIR)/parameters.f90 $(SRCDIR)/BC-User.f90 $(SRCDIR)/BC-TGV.f90 $(SRCDIR)/BC-Channel-flow.f90 $(SRCDIR)/BC-Periodic-hill.f90 $(SRCDIR)/BC-Cylinder.f90 $(SRCDIR)/BC-Mixing-layer.f90 $(SRCDIR)/BC-Jet.f90 $(SRCDIR)/BC-Sandbox.f90 $(SRCDIR)/BC-Uniform.f90 $(TURBDIR)/constants.f90 $(TURBDIR)/acl_utils.f90 $(TURBDIR)/airfoils.f90 $(TURBDIR)/dynstall.f90 $(TURBDIR)/dynstall_legacy.f90 $(TURBDIR)/acl_elem.f90 $(TURBDIR)/acl_controller.f90 $(TURBDIR)/acl_turb.f90 $(TURBDIR)/acl_out.f90 $(TURBDIR)/acl_farm_controller.f90 $(TURBDIR)/acl_model.f90 $(TURBDIR)/acl_source.f90 $(TURBDIR)/adm.f90 $(TURBDIR)/turbine.f90 $(SRCDIR)/statistics.f90 $(SRCDIR)/case.f90 $(SRCDIR)/transeq.f90 $(SRCDIR)/genepsi3d.f90 $(SRCDIR)/xcompact3d.f90


#######FFT settings##########
ifeq ($(FFT),fftw3)
  #FFTW3_PATH=/usr
  #FFTW3_PATH=/usr/lib64
  FFTW3_PATH=/usr/local/Cellar/fftw/3.3.7_1
  INC:=-I$(FFTW3_PATH)/include
  LIBFFT=-L$(FFTW3_PATH) -lfftw3 -lfftw3f
else ifeq ($(FFT),fftw3_f03)
  FFTW3_PATH=/usr                                #ubuntu # apt install libfftw3-dev
  #FFTW3_PATH=/usr/lib64                         #fedora # dnf install fftw fftw-devel
  #FFTW3_PATH=/usr/local/Cellar/fftw/3.3.7_1     #macOS  # brew install fftw
  INC:=-I$(FFTW3_PATH)/include
  LIBFFT=-L$(FFTW3_PATH)/lib -lfftw3 -lfftw3f
else ifeq ($(FFT),generic)
  INC:=
  LIBFFT=
else ifeq ($(FFT),mkl)
  SRCDECOMP := $(DECOMPDIR)/mkl_dfti.f90 $(SRCDECOMP)
  LIBFFT=-Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread
	INC=-I$(MKLROOT)/include
else ifeq ($(FFT),ffte)
  INC:=
  LIBFFT:=-L$(FFTE_DIR)/lib -lffte
endif

#######OPTIONS settings###########
OPT := -I$(SRCDIR) -I$(DECOMPDIR) $(FFLAGS)
LINKOPT := $(FFLAGS)

LIBIO :=
ADIOS2DIR :=
ifeq ($(IO),adios2)
  ifeq ($(ADIOS2DIR),)
    $(error Set ADIOS2DIR=/path/to/adios2/install/)
  endif
  OPT := -DADIOS2 $(OPT)
  INC := $(INC) $(shell $(ADIOS2DIR)/bin/adios2-config --fortran-flags) #$(patsubst $(shell $(ADIOS2DIR)/bin/adios2-config --fortran-libs),,$(shell $(ADIOS2DIR)/bin/adios2-config -f))
  LIBIO := $(shell $(ADIOS2DIR)/bin/adios2-config --fortran-libs)
endif

#-----------------------------------------------------------------------
# Normally no need to change anything below

all: xcompact3d

xcompact3d : $(OBJDECOMP) $(OBJ)
	$(FC) -o $@ $(LINKOPT) $(OBJDECOMP) $(OBJ) $(LIBFFT) $(LIBIO)

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
	rm -f $(DECOMPDIR)/*.o $(DECOMPDIR)/*.mod $(DECOMPDIR)/*.smod
	rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod $(SRCDIR)/*.smod
	rm -f *.o *.mod *.smod xcompact3d

.PHONY: cleanall
cleanall: clean
	rm -f *~ \#*\# out/* data/* stats/* planes/* *.xdmf *.log *.out nodefile core sauve*
