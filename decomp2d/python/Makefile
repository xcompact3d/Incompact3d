
FC = mpif90
F77 = mpif77

INSTALL_DIR = ~/.local

PIP = pip3 #--verbose
PIPFLAGS = --prefix=$(INSTALL_DIR) \
	--install-option="config_fc" \
	--install-option="--verbose" \
	--install-option="--f90exec=$(FC)" \
	--install-option="--f77exec=$(F77)" \
	# --install-option=build_src \
	# --install-option=--f2cmap="/Users/paulbartholomew/src/fortran/Xcompact3d/decomp2d/python/.f2py_f2cmap"

ifeq ($(IO),adios2)
	IO_MACROS = [ ('ADIOS2', None), ('ADIOS2_USE_MPI', None) ]
else
	IO_MACROS = []
endif

all: 2decomp4py.f90 .f2py_f2cmap decomp2d_options.py decomp2d_io_options.py
	$(PIP) install $(PIPFLAGS) .

clean:
	$(PIP) uninstall -y decomp-2d4py
	rm -f *.o *.mod
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

.PHONY: decomp2d_io_options.py
decomp2d_io_options.py:
	echo "define_macros = $(IO_MACROS)" > $@
