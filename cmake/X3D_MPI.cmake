# MPI CMakeLists

find_package(MPI REQUIRED)
# Stop if there is no MPI_Fortran_Compiler
if (MPI_Fortran_COMPILER)
    message(STATUS "X3D MPI_Fortran_COMPILER found: ${MPI_Fortran_COMPILER}")
    message(STATUS "X3D MPI_VERSION found: ${MPI_VERSION}")
    # Try to guess the MPI type to adapt compilation flags if necessary
    string(FIND "${MPI_Fortran_COMPILER}" "mpich" pos)
    if(pos GREATER_EQUAL "0")
      set(FIND_MPICH TRUE)
      message(STATUS "X3D MPI is MPICH type")
    endif()
    string(FIND "${MPI_Fortran_COMPILER}" "openmpi" pos)
    if(pos GREATER_EQUAL "0")
      set(FIND_OMPI TRUE)
      message(STATUS "X3D MPI is openMPI type")
    endif()
else (MPI_Fortran_COMPILER)
    message(SEND_ERROR "This application cannot compile without MPI")
endif(MPI_Fortran_COMPILER)
# Warning if Include are not found => can be fixed with more recent cmake version
if (MPI_FOUND)
  message(STATUS "X3D MPI FOUND: ${MPI_FOUND}")
  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
  message(STATUS "X3D MPI INCL ALSO FOUND: ${MPI_INCLUDE_PATH}")
  # Force the mpirun to be coherent with the mpifortran
  # We do this only at first iteration of configure after user will take charge
  # This is the same workaround used for 2DECOMP&FFT -> it is necessary to have a new
  # variable name for the PATH_TO_MPIRUN to differentiate when 2D&FFT is build autonatically via CMake
  if (NOT X3D_PATH_TO_MPIRUN)
    string(REGEX REPLACE "mpif90" "mpirun" X3D_PATH_TO_MPIRUN "${MPI_Fortran_COMPILER}")
    string(REPLACE "mpif90" "mpirun" X3D_PATH_TO_MPIRUN "${X3D_PATH_TO_MPIRUN}")
    string(REPLACE "mpiifort" "mpirun" X3D_PATH_TO_MPIRUN "${X3D_PATH_TO_MPIRUN}")
    string(REPLACE "mpiifx" "mpirun" X3D_PATH_TO_MPIRUN "${X3D_PATH_TO_MPIRUN}")
    message(STATUS "X3D Path to mpirun ${X3D_PATH_TO_MPIRUN}")
    set(MPIEXEC_EXECUTABLE "${X3D_PATH_TO_MPIRUN}" CACHE STRING
        "Force MPIRUN to be consistent with MPI_Fortran_COMPILER" FORCE)
    # Save this variable to avoid repeting the substitution
    set(X3D_PATH_TO_MPIRUN "${X3D_PATH_TO_MPIRUN}" CACHE STRING
	"Force MPIRUN to be consistent with MPI_Fortran_COMPILER" FORCE)
  endif (NOT X3D_PATH_TO_MPIRUN)
else (MPI_FOUND)
  message(STATUS "NO MPI include have been found. The executable won't be targeted with MPI include")
  message(STATUS "Code will compile but performaces can be compromised")
  message(STATUS "Alternatively use ccmake to manually set the include if available")
endif (MPI_FOUND)

