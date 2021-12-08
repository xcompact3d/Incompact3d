USEMPI := YES
PFUNIT_DIR := ../../../unit_testing/pFUnit/build/installed
LATEST_PFUNIT_DIR := $(lastword $(shell echo $(wildcard $(PFUNIT_DIR)/PFUNIT-4.*) | xargs -n1 | sort -V))
include $(LATEST_PFUNIT_DIR)/include/PFUNIT.mk

all: tgv2D_00 tgv2D_00_imp1 tgv2D_11 tgv2D_11_imp2

X3D := ../../src/x3d_tools.o ../../src/les_models.o ../../src/signal.o ../../src/genepsi3d.o ../../src/statistics.o ../../src/dynstall*.o ../../src/a*.o ../../src/turbine.o ../../src/poisson.o ../../src/implicit.o ../../src/ibm.o ../../src/schemes.o ../../src/filters.o ../../src/navier.o ../../src/derive.o ../../src/tools.o ../../src/visu.o ../../src/BC*.o ../../src/case.o ../../src/variables.o ../../src/probes.o ../../src/forces.o ../../src/module_param.o ../../src/parameters.o ../../src/thomas.o ../../src/transeq.o ../../src/time_integrators.o ../../decomp2d/*.o ../lib/*.o

%.o : %.F90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(DEFS2) $(INC) $(PFUNIT_EXTRA_FFLAGS) -c $<

FFLAGS += -I../../ -I../lib/

tgv2D_00_TESTS := init00.pf test00.pf final.pf
tgv2D_00_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,tgv2D_00))

tgv2D_00_imp1_TESTS := init00_imp1.pf test00_imp1.pf final.pf
tgv2D_00_imp1_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,tgv2D_00_imp1))

tgv2D_11_TESTS := init11.pf test11.pf final.pf
tgv2D_11_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,tgv2D_11))

tgv2D_11_imp2_TESTS := init11_imp2.pf test11_imp2.pf final.pf
tgv2D_11_imp2_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,tgv2D_11_imp2))

clean:
	$(RM) *.o *.mod *.a *.inc *.dat *.png *_tmp.i3d
	$(RM) final.F90
	$(RM) tgv2D_00 init00.F90 test00.F90
	$(RM) tgv2D_00_imp1 init00_imp1.F90 test00_imp1.F90
	$(RM) tgv2D_11 init11.F90 test11.F90
	$(RM) tgv2D_11_imp2 init11_imp2.F90 test11_imp2.F90
	$(RM) -rf data/ out/ probes/
	$(RM) yp.dat ypi.dat time_evol.dat
