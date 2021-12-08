USEMPI := YES
PFUNIT_DIR := ../../../unit_testing/pFUnit/build/installed
LATEST_PFUNIT_DIR := $(lastword $(shell echo $(wildcard $(PFUNIT_DIR)/PFUNIT-4.*) | xargs -n1 | sort -V))
include $(LATEST_PFUNIT_DIR)/include/PFUNIT.mk

all: verif_x_00 verif_y_00 verif_y_00_str verif_z_00 verif_x_11 verif_y_11 verif_y_11_str verif_z_11 verif_x_22 verif_y_22 verif_y_22_str verif_z_22 verif_x_12 verif_y_12 verif_y_12_str verif_z_12 verif_x_21 verif_y_21 verif_y_21_str verif_z_21

X3D := ../../src/x3d_tools.o ../../src/les_models.o ../../src/signal.o ../../src/genepsi3d.o ../../src/statistics.o ../../src/dynstall*.o ../../src/a*.o ../../src/turbine.o ../../src/poisson.o ../../src/implicit.o ../../src/ibm.o ../../src/schemes.o ../../src/filters.o ../../src/navier.o ../../src/derive.o ../../src/tools.o ../../src/visu.o ../../src/BC*.o ../../src/case.o ../../src/variables.o ../../src/probes.o ../../src/forces.o ../../src/module_param.o ../../src/parameters.o ../../src/thomas.o ../../src/transeq.o ../../src/time_integrators.o ../../decomp2d/*.o ../lib/*.o

%.o : %.F90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(DEFS2) $(INC) $(PFUNIT_EXTRA_FFLAGS) -c $<

FFLAGS += -I../../ -I../lib/

verif_x_00_TESTS := initx00.pf testx00.pf final.pf
verif_x_00_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,verif_x_00))

verif_y_00_TESTS := inity00.pf testy00.pf final.pf
verif_y_00_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,verif_y_00))

verif_y_00_str_TESTS := inity00_str.pf testy00_str.pf final.pf
verif_y_00_str_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,verif_y_00_str))

verif_z_00_TESTS := initz00.pf testz00.pf final.pf
verif_z_00_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,verif_z_00))

verif_x_11_TESTS := initx11.pf testx11.pf final.pf
verif_x_11_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,verif_x_11))

verif_y_11_TESTS := inity11.pf testy11.pf final.pf
verif_y_11_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,verif_y_11))

verif_y_11_str_TESTS := inity11_str.pf testy11_str.pf final.pf
verif_y_11_str_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,verif_y_11_str))

verif_z_11_TESTS := initz11.pf testz11.pf final.pf
verif_z_11_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,verif_z_11))

verif_x_22_TESTS := initx22.pf testx22.pf final.pf
verif_x_22_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,verif_x_22))

verif_y_22_TESTS := inity22.pf testy22.pf final.pf
verif_y_22_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,verif_y_22))

verif_y_22_str_TESTS := inity22_str.pf testy22_str.pf final.pf
verif_y_22_str_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,verif_y_22_str))

verif_z_22_TESTS := initz22.pf testz22.pf final.pf
verif_z_22_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,verif_z_22))

verif_x_12_TESTS := initx12.pf testx12.pf final.pf
verif_x_12_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,verif_x_12))

verif_y_12_TESTS := inity12.pf testy12.pf final.pf
verif_y_12_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,verif_y_12))

verif_y_12_str_TESTS := inity12_str.pf testy12_str.pf final.pf
verif_y_12_str_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,verif_y_12_str))

verif_z_12_TESTS := initz12.pf testz12.pf final.pf
verif_z_12_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,verif_z_12))

verif_x_21_TESTS := initx21.pf testx21.pf final.pf
verif_x_21_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,verif_x_21))

verif_y_21_TESTS := inity21.pf testy21.pf final.pf
verif_y_21_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,verif_y_21))

verif_y_21_str_TESTS := inity21_str.pf testy21_str.pf final.pf
verif_y_21_str_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,verif_y_21_str))

verif_z_21_TESTS := initz21.pf testz21.pf final.pf
verif_z_21_OTHER_LIBRARIES := $(X3D) -L../lib -lx3dverif
$(eval $(call make_pfunit_test,verif_z_21))

clean:
	$(RM) *.o *.mod *.a *.inc *.dat *.png *_tmp.i3d
	$(RM) final.F90
	$(RM) verif_x_00 initx00.F90 testx00.F90
	$(RM) verif_y_00 inity00.F90 testy00.F90
	$(RM) verif_y_00_str inity00_str.F90 testy00_str.F90
	$(RM) verif_z_00 initz00.F90 testz00.F90
	$(RM) verif_x_11 initx11.F90 testx11.F90
	$(RM) verif_y_11 inity11.F90 testy11.F90
	$(RM) verif_y_11_str inity11_str.F90 testy11_str.F90
	$(RM) verif_z_11 initz11.F90 testz11.F90
	$(RM) verif_x_22 initx22.F90 testx22.F90
	$(RM) verif_y_22 inity22.F90 testy22.F90
	$(RM) verif_y_22_str inity22_str.F90 testy22_str.F90
	$(RM) verif_z_22 initz22.F90 testz22.F90
	$(RM) verif_x_12 initx12.F90 testx12.F90
	$(RM) verif_y_12 inity12.F90 testy12.F90
	$(RM) verif_y_12_str inity12_str.F90 testy12_str.F90
	$(RM) verif_z_12 initz12.F90 testz12.F90
	$(RM) verif_x_21 initx21.F90 testx21.F90
	$(RM) verif_y_21 inity21.F90 testy21.F90
	$(RM) verif_y_21_str inity21_str.F90 testy21_str.F90
	$(RM) verif_z_21 initz21.F90 testz21.F90
	$(RM) -rf data/ out/ probes/
	$(RM) yp.dat ypi.dat channel.dat time_evol.dat
