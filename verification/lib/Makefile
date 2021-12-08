SRCS := x3d_unit_testing_tools.f90 verif_TGV2D.f90
OBJS := $(../decomp2d/*.o) $(../src/*.o) $(./*.o)

all: libx3dverif.a $(SRCS:%.f90=%.o)

libx3dverif.a: $(OBJS)
	$(AR) -r $@ $?

%.o : %.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(DEFS2) $(INC) -c $<

FFLAGS += -I../../

clean:
	$(RM) libx3dverif.a *.o *.mod

