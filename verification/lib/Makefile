SRCS := x3d_unit_testing_tools.f90 verif_TGV2D.f90
OBJS := $(../decomp2d/*.o) $(../src/*.o) $(./*.o)

all: libx3d.a x3d_unit_testing_tools.o verif_TGV2D.o

libx3d.a: $(OBJS)
	$(AR) -r $@ $?

%.o : %.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(DEFS2) $(INC) -c $<

FFLAGS += -I../../

clean:
	$(RM) libx3d.a *.o *.mod

