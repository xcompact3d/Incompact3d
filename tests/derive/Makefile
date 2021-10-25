USEMPI := YES
PFUNIT_DIR := ../../../unit_testing/pFUnit/build/installed
LATEST_PFUNIT_DIR := $(lastword $(shell echo $(wildcard $(PFUNIT_DIR)/PFUNIT-4.*) | xargs -n1 | sort -V))
include $(LATEST_PFUNIT_DIR)/include/PFUNIT.mk

all: test_x_00 test_y_00 test_y_00_str test_z_00

X3D := ../../src/x3d_tools.o ../../src/les_models.o ../../src/signal.o ../../src/genepsi3d.o ../../src/statistics.o ../../src/dynstall*.o ../../src/a*.o ../../src/turbine.o ../../src/poisson.o ../../src/implicit.o ../../src/ibm.o ../../src/schemes.o ../../src/filters.o ../../src/navier.o ../../src/derive.o ../../src/tools.o ../../src/visu.o ../../src/BC*.o ../../src/case.o ../../src/variables.o ../../src/probes.o ../../src/forces.o ../../src/module_param.o ../../src/parameters.o ../../src/thomas.o ../../decomp2d/*.o

%.o : %.F90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(DEFS2) $(INC) $(PFUNIT_EXTRA_FFLAGS) -c $<

FFLAGS += -I../../

test_x_00_TESTS := initx00.pf testx00.pf final.pf
test_x_00_OTHER_LIBRARIES := $(X3D) -L../lib -lx3d
$(eval $(call make_pfunit_test,test_x_00))

test_y_00_TESTS := inity00.pf testy00.pf final.pf
test_y_00_OTHER_LIBRARIES := $(X3D) -L../lib -lx3d
$(eval $(call make_pfunit_test,test_y_00))

test_y_00_str_TESTS := inity00_str.pf testy00_str.pf final.pf
test_y_00_str_OTHER_LIBRARIES := $(X3D) -L../lib -lx3d
$(eval $(call make_pfunit_test,test_y_00_str))

test_z_00_TESTS := initz00.pf testz00.pf final.pf
test_z_00_OTHER_LIBRARIES := $(X3D) -L../lib -lx3d
$(eval $(call make_pfunit_test,test_z_00))

clean:
	$(RM) *.o *.mod *.a  *.inc *.dat *.png
	$(RM) final.F90
	$(RM) test_x_00 initx00.F90 testx00.F90
	$(RM) test_y_00 inity00.F90 testy00.F90
	$(RM) test_y_00_str inity00_str.F90 testy00_str.F90
	$(RM) test_z_00 initz00.F90 testz00.F90
	$(RM) -rf data/ out/ probes/
	$(RM) yp.dat ypi.dat channel.dat time_evol.dat
