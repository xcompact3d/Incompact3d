!! SPDX-License-Identifier: BSD-3-Clause

! Constants for the 2decomp&fft library

module decomp_2d_constants

   use mpi
   use, intrinsic :: iso_fortran_env, only: real32, real64
#if defined(_GPU) && defined(_NCCL)
   use nccl
#endif

   implicit none

   !private

#ifdef DOUBLE_PREC
   integer, parameter, public :: mytype = KIND(0._real64)
   integer, parameter, public :: real_type = MPI_DOUBLE_PRECISION
   integer, parameter, public :: real2_type = MPI_2DOUBLE_PRECISION
   integer, parameter, public :: complex_type = MPI_DOUBLE_COMPLEX
#ifdef SAVE_SINGLE
   integer, parameter, public :: mytype_single = KIND(0._real32)
   integer, parameter, public :: real_type_single = MPI_REAL
#else
   integer, parameter, public :: mytype_single = KIND(0._real64)
   integer, parameter, public :: real_type_single = MPI_DOUBLE_PRECISION
#endif
#else
   integer, parameter, public :: mytype = KIND(0._real32)
   integer, parameter, public :: real_type = MPI_REAL
   integer, parameter, public :: real2_type = MPI_2REAL
   integer, parameter, public :: complex_type = MPI_COMPLEX
   integer, parameter, public :: mytype_single = KIND(0._real32)
   integer, parameter, public :: real_type_single = MPI_REAL
#endif

   !
   ! Output for the log can be changed by the external code before calling decomp_2d_init
   !
   !    0 => No log output
   !    1 => Master rank log output to stdout
   !    2 => Master rank log output to the file "decomp_2d_setup.log"
   !    3 => All ranks log output to a dedicated file
   !
   ! The default value is 2 (3 for debug builds)
   !
   integer, parameter, public :: D2D_LOG_QUIET = 0
   integer, parameter, public :: D2D_LOG_STDOUT = 1
   integer, parameter, public :: D2D_LOG_TOFILE = 2
   integer, parameter, public :: D2D_LOG_TOFILE_FULL = 3

   !
   ! Debug level can be changed by the external code before calling decomp_2d_init
   !
   ! The environment variable "DECOMP_2D_DEBUG" can be used to change the debug level
   !
   ! Debug checks are performed only when the preprocessor variable DEBUG is defined
   !
   enum, bind(c)
      enumerator :: D2D_DEBUG_LEVEL_OFF = 0
      enumerator :: D2D_DEBUG_LEVEL_CRITICAL = 1
      enumerator :: D2D_DEBUG_LEVEL_ERROR = 2
      enumerator :: D2D_DEBUG_LEVEL_WARN = 3
      enumerator :: D2D_DEBUG_LEVEL_INFO = 4
      enumerator :: D2D_DEBUG_LEVEL_DEBUG = 5
      enumerator :: D2D_DEBUG_LEVEL_TRACE = 6
   end enum

   !
   ! Profiler section
   !
   ! Integer to select the profiling tool
   !    0 => no profiling, default
   !    1 => Caliper (https://github.com/LLNL/Caliper)
   !
   integer, parameter, public :: DECOMP_PROFILER_NONE = 0
   integer, parameter, public :: DECOMP_PROFILER_CALIPER = 1

   !
   ! Supported FFT backends
   !
   integer, parameter, public :: D2D_FFT_BACKEND_GENERIC = 0
   integer, parameter, public :: D2D_FFT_BACKEND_FFTW3 = 1
   integer, parameter, public :: D2D_FFT_BACKEND_FFTW3_F03 = 2
   integer, parameter, public :: D2D_FFT_BACKEND_MKL = 3
   integer, parameter, public :: D2D_FFT_BACKEND_CUFFT = 4

   !
   ! Complex-to-complex FFT can be forward or backward
   !
   integer, parameter, public :: DECOMP_2D_FFT_FORWARD = -1
   integer, parameter, public :: DECOMP_2D_FFT_BACKWARD = 1

   !
   ! Input / output of the FFT are distributed along different pencils
   !
   integer, parameter, public :: PHYSICAL_IN_X = 1 ! Forward is input in X, output in Z
   integer, parameter, public :: PHYSICAL_IN_Z = 3 ! Forward is input in Z, output in X

   !
   ! Major and minor version number
   !
   integer, parameter :: D2D_MAJOR = 2
   integer, parameter :: D2D_MINOR = 0
   logical, parameter :: D2D_RELEASE = .true.

end module decomp_2d_constants

