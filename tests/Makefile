USEMPI := YES
PFUNIT_DIR := ../../unit_testing/pFUnit/build/installed
LATEST_PFUNIT_DIR := $(lastword $(shell echo $(wildcard $(PFUNIT_DIR)/PFUNIT-4.*) | xargs -n1 | sort -V))
include $(LATEST_PFUNIT_DIR)/include/PFUNIT.mk
export FFLAGS OPT DEFS DEFS2 INC

all:
	$(MAKE) -C lib all
	$(MAKE) -C verif all

clean:
	$(MAKE) -C lib clean
	$(MAKE) -C verif clean
