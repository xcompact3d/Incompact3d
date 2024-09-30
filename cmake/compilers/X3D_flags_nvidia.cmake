#Compilers Flags for NVIDIA

set(X3D_FFLAGS "-cpp -Mfree -Kieee")
set(X3D_FFLAGS_RELEASE "-O3 -fast -march=native")
set(X3D_FFLAGS_DEBUG   "-O0 -g -traceback -Mbounds -Mchkptr -Ktrap=fp")
set(X3D_FFLAGS_DEV     "${X3D_FFLAGS_DEBUG}")

