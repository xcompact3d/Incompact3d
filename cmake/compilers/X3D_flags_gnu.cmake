# Flags for GNU compiler
#
set(X3D_FFLAGS "-cpp -ffree-line-length-none -Warray-bounds -fcray-pointer -fbacktrace")
if (CMAKE_Fortran_COMPILER_VERSION GREATER_EQUAL "10")
  message(STATUS "Set New Fortran basic flags")
  set(X3D_FFLAGS "${X3D_FFLAGS} -fallow-argument-mismatch")
  set(X3D_GNU10 TRUE)
else (CMAKE_Fortran_COMPILER_VERSION GREATER_EQUAL "10")
  set(X3D_GNU10 FALSE)
endif (CMAKE_Fortran_COMPILER_VERSION GREATER_EQUAL "10")
# Flags extracted from the Makefile
set(X3D_FFLAGS_RELEASE "-O3 -funroll-loops -floop-optimize -g") #"-O3 -funroll-loops -floop-optimize -march=native")
set(X3D_FFLAGS_DEBUG "-g -Og -ffpe-trap=invalid,zero -fcheck=bounds -fimplicit-none") #"-DDEBUG -g3 -Og -ffpe-trap=invalid,zero -fcheck=all -fimplicit-none")
# Dev flag below is new and not working yet
if (FIND_MPICH AND X3D_GNU10)
  set(X3D_FFLAGS_DEV     "${X3D_FFLAGS_DEBUG} -Wall -Wno-unused-function -Wno-integer-division")
else()
  set(X3D_FFLAGS_DEV     "${X3D_FFLAGS_DEBUG} -Wall -Wpedantic -Wno-unused-function -Werror -Wno-integer-division")
endif()
if (FIND_OMPI)
  set(X3D_FFLAGS_DEV     "${X3D_FFLAGS_DEV} -Wimplicit-procedure -Wimplicit-interface")
endif()
