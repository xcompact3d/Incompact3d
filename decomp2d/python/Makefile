
FC = mpif90

INSTALL_DIR = ~/.local

PIP = pip3
PIPFLAGS = --prefix=$(INSTALL_DIR) \
	--install-option=config_fc --install-option=--f90exec=$(FC)

all: .f2py_f2cmap decomp2d_options.py
	$(PIP) install $(PIPFLAGS) .

clean:
	$(PIP) uninstall decomp-2d4py

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
