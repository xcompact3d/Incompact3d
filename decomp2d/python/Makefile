
FC = mpif90

INSTALL_DIR = ~/.local

PIP = pip3
PIPFLAGS = --prefix=$(INSTALL_DIR) \
	--install-option=config_fc --install-option=--f90exec=$(FC)

all:
	$(PIP) install $(PIPFLAGS) .

clean:
	$(PIP) uninstall decomp-2d4py
