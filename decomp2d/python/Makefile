
FC = mpif90
F77 = mpif77

INSTALL_DIR = ~/.local

ifeq ($(IO),adios2)
	PIPINC = $(shell $(ADIOS2DIR)/bin/adios2-config --fortran-flags)
	PIPINC := $(PIPINC:-D%=--install-option="-D%")
	PIPINC := $(PIPINC:-I%=--install-option="-I%")
	PIPLIB = $(shell $(ADIOS2DIR)/bin/adios2-config --fortran-libs)
	PIPLIB := $(PIPLIB:%.a=-L$(dir %) -l$(notdir %))
	PIPLIB := $(PIPLIB:-llib%=-l%)
	PIPLIB := $(PIPLIB:%=--install-option="%")
	PIPDEF = -DADIOS2
	PIPDEF := $(PIPDEF:%=--install-option="%")
endif

PIP = pip3 #--verbose
PIPFLAGS = --prefix=$(INSTALL_DIR) \
	--install-option="config_fc" \
	--install-option="--verbose" \
	--install-option="--f90exec=$(FC)" \
	--install-option="--f77exec=$(F77)" \
	--install-option="build_ext" \
	$(PIPDEF) \
	$(PIPINC) \
	$(PIPLIB)
	# --install-option=build_src \
	# --install-option=--f2cmap="/Users/paulbartholomew/src/fortran/Xcompact3d/decomp2d/python/.f2py_f2cmap"

all: .f2py_f2cmap decomp2d_options.py
	$(PIP) install $(PIPFLAGS) .

clean:
	$(PIP) uninstall -y decomp-2d4py
	rm -rf build decomp_2d4py.egg-info __pycache__

# TODO: Make this settable based on presence of -DDOUBLE_PREC in a variable.
#       In principle this would then be picked up when called by a Makefile in 2DECOMP&FFT or
#       Xcompact3d, i.e.
#         with -DDOUBL_PREC:
#           mytype maps to double
#         else:
#           mytype maps to float
.PHONY: .f2py_f2cmap
.f2py_f2cmap:
	echo "{'real' : {'kind(0.0d0)':'double', 'kind(0.0)':'float'}}" > $@

.PHONY: decomp2d_options.py
decomp2d_options.py:
	echo "define_macros = [ ('DOUBLE_PREC', None) ]" > $@
