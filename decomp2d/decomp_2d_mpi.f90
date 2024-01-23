!! SPDX-License-Identifier: BSD-3-Clause

! MPI module for 2decomp&fft library

module decomp_2d_mpi

   use MPI
   use decomp_2d_constants
#if defined(_GPU) && defined(_NCCL)
   use nccl
#endif

   implicit none

   integer, save, public :: mytype_bytes
   integer, save, public :: nrank = -1 ! local MPI rank
   integer, save, public :: nproc = -1 ! total number of processors
   integer, save, public :: decomp_2d_comm = MPI_COMM_NULL ! MPI communicator

   public :: decomp_2d_mpi_init, &
             decomp_2d_mpi_fin, &
             decomp_2d_mpi_comm_free, &
             decomp_2d_abort, &
             decomp_2d_warning

   interface decomp_2d_abort
      module procedure decomp_2d_abort_basic
      module procedure decomp_2d_abort_file_line
#if defined(_GPU) && defined(_NCCL)
      module procedure decomp_2d_abort_nccl_basic
      module procedure decomp_2d_abort_nccl_file_line
#endif
   end interface decomp_2d_abort

   interface decomp_2d_warning
      module procedure decomp_2d_warning_basic
      module procedure decomp_2d_warning_file_line
   end interface decomp_2d_warning

contains

   !
   ! Initialize the MPI module
   !
   subroutine decomp_2d_mpi_init(comm)

      implicit none

      integer, intent(in), optional :: comm

      integer :: ierror

      ! Use the provided MPI communicator if present
      if (present(comm)) then
         decomp_2d_comm = comm
      else
         decomp_2d_comm = MPI_COMM_WORLD
      end if

      ! If the external code has not set nrank and nproc
      if (nrank == -1) then
         call MPI_COMM_RANK(decomp_2d_comm, nrank, ierror)
         if (ierror /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_RANK")
         end if
      end if
      if (nproc == -1) then
         call MPI_COMM_SIZE(decomp_2d_comm, nproc, ierror)
         if (ierror /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_SIZE")
         end if
      end if

      if (MPI_SUCCESS /= 0) then
         call decomp_2d_abort(__FILE__, __LINE__, MPI_SUCCESS, "MPI error check is broken")
      end if

   end subroutine decomp_2d_mpi_init

   !
   ! Finalize the MPI module
   !
   subroutine decomp_2d_mpi_fin

      implicit none

      nrank = -1
      nproc = -1

   end subroutine decomp_2d_mpi_fin

   !
   ! Small wrapper to free a MPI communicator
   !
   subroutine decomp_2d_mpi_comm_free(mpi_comm)

      implicit none

      integer, intent(inout) :: mpi_comm
      integer :: ierror

      ! Return if no MPI comm to free
      if (mpi_comm == MPI_COMM_NULL) return

      ! Free the provided MPI communicator
      call MPI_COMM_FREE(mpi_comm, ierror)
      if (ierror /= 0) call decomp_2d_warning(__FILE__, __LINE__, ierror, "MPI_COMM_FREE")
      mpi_comm = MPI_COMM_NULL

   end subroutine decomp_2d_mpi_comm_free

   subroutine decomp_2d_abort_basic(errorcode, msg)

      use iso_fortran_env, only: error_unit

      implicit none

      integer, intent(IN) :: errorcode
      character(len=*), intent(IN) :: msg

      integer :: ierror

      if (nrank == 0) then
         write (*, *) '2DECOMP&FFT ERROR - errorcode: ', errorcode
         write (*, *) 'ERROR MESSAGE: '//msg
         write (error_unit, *) '2DECOMP&FFT ERROR - errorcode: ', errorcode
         write (error_unit, *) 'ERROR MESSAGE: '//msg
      end if
      call MPI_ABORT(decomp_2d_comm, errorcode, ierror)

   end subroutine decomp_2d_abort_basic

   subroutine decomp_2d_abort_file_line(file, line, errorcode, msg)

      use iso_fortran_env, only: error_unit

      implicit none

      integer, intent(IN) :: errorcode, line
      character(len=*), intent(IN) :: msg, file

      integer :: ierror

      if (nrank == 0) then
         write (*, *) '2DECOMP&FFT ERROR'
         write (*, *) '  errorcode:     ', errorcode
         write (*, *) '  error in file  '//file
         write (*, *) '           line  ', line
         write (*, *) '  error message: '//msg
         write (error_unit, *) '2DECOMP&FFT ERROR'
         write (error_unit, *) '  errorcode:     ', errorcode
         write (error_unit, *) '  error in file  '//file
         write (error_unit, *) '           line  ', line
         write (error_unit, *) '  error message: '//msg
      end if
      call MPI_ABORT(decomp_2d_comm, errorcode, ierror)

   end subroutine decomp_2d_abort_file_line

#if defined(_GPU) && defined(_NCCL)
   !
   ! This is based on the file "nccl.h" in nvhpc 22.1
   !
   function _ncclresult_to_integer(errorcode)

      implicit none

      type(ncclresult), intent(IN) :: errorcode
      integer :: _ncclresult_to_integer

      if (errorcode == ncclSuccess) then
         _ncclresult_to_integer = 0
      elseif (errorcode == ncclUnhandledCudaError) then
         _ncclresult_to_integer = 1
      elseif (errorcode == ncclSystemError) then
         _ncclresult_to_integer = 2
      elseif (errorcode == ncclInternalError) then
         _ncclresult_to_integer = 3
      elseif (errorcode == ncclInvalidArgument) then
         _ncclresult_to_integer = 4
      elseif (errorcode == ncclInvalidUsage) then
         _ncclresult_to_integer = 5
      elseif (errorcode == ncclNumResults) then
         _ncclresult_to_integer = 6
      else
         _ncclresult_to_integer = -1
         call decomp_2d_warning(__FILE__, __LINE__, _ncclresult_to_integer, &
                                "NCCL error handling needs some update")
      end if

   end function _ncclresult_to_integer
   !
   ! Small wrapper for basic NCCL errors
   !
   subroutine decomp_2d_abort_nccl_basic(errorcode, msg)

      implicit none

      type(ncclresult), intent(IN) :: errorcode
      character(len=*), intent(IN) :: msg

      call decomp_2d_abort(_ncclresult_to_integer(errorcode), &
                           msg//" "//ncclGetErrorString(errorcode))

   end subroutine decomp_2d_abort_nccl_basic

   !
   ! Small wrapper for NCCL errors
   !
   subroutine decomp_2d_abort_nccl_file_line(file, line, errorcode, msg)

      implicit none

      type(ncclresult), intent(IN) :: errorcode
      integer, intent(in) :: line
      character(len=*), intent(IN) :: msg, file

      call decomp_2d_abort(file, &
                           line, &
                           _ncclresult_to_integer(errorcode), &
                           msg//" "//ncclGetErrorString(errorcode))

   end subroutine decomp_2d_abort_nccl_file_line
#endif

   subroutine decomp_2d_warning_basic(errorcode, msg)

      use iso_fortran_env, only: error_unit

      implicit none

      integer, intent(IN) :: errorcode
      character(len=*), intent(IN) :: msg

      if (nrank == 0) then
         write (*, *) '2DECOMP&FFT WARNING - errorcode: ', errorcode
         write (*, *) 'ERROR MESSAGE: '//msg
         write (error_unit, *) '2DECOMP&FFT WARNING - errorcode: ', errorcode
         write (error_unit, *) 'ERROR MESSAGE: '//msg
      end if

   end subroutine decomp_2d_warning_basic

   subroutine decomp_2d_warning_file_line(file, line, errorcode, msg)

      use iso_fortran_env, only: error_unit

      implicit none

      integer, intent(IN) :: errorcode, line
      character(len=*), intent(IN) :: msg, file

      if (nrank == 0) then
         write (*, *) '2DECOMP&FFT WARNING'
         write (*, *) '  errorcode:     ', errorcode
         write (*, *) '  error in file  '//file
         write (*, *) '           line  ', line
         write (*, *) '  error message: '//msg
         write (error_unit, *) '2DECOMP&FFT WARNING'
         write (error_unit, *) '  errorcode:     ', errorcode
         write (error_unit, *) '  error in file  '//file
         write (error_unit, *) '           line  ', line
         write (error_unit, *) '  error message: '//msg
      end if

   end subroutine decomp_2d_warning_file_line

end module decomp_2d_mpi
