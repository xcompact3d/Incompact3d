!! SPDX-License-Identifier: BSD-3-Clause

submodule(decomp_2d) d2d_log

   use decomp_2d_constants
   use decomp_2d_mpi

   implicit none

contains

   !
   ! Get the IO unit for the log
   !
   module function d2d_listing_get_unit()

      use iso_fortran_env, only: output_unit

      implicit none

      ! Output
      integer :: d2d_listing_get_unit

      ! Local variables
      logical :: found
      integer :: ierror
      character(len=7) fname ! Sufficient for up to O(1M) ranks

      if (decomp_log == D2D_LOG_TOFILE_FULL) then
         write (fname, "(I0)") nrank ! Adapt to magnitude of nrank
         inquire (file='decomp_2d_setup_'//trim(fname)//'.log', &
                  exist=found)
         if (found) then
            open (newunit=d2d_listing_get_unit, &
                  file='decomp_2d_setup_'//trim(fname)//'.log', &
                  status="old", &
                  position="append", &
                  iostat=ierror)
         else
            open (newunit=d2d_listing_get_unit, &
                  file='decomp_2d_setup_'//trim(fname)//'.log', &
                  status="new", &
                  iostat=ierror)
         end if
      elseif (nrank == 0 .and. decomp_log == D2D_LOG_TOFILE) then
         inquire (file="decomp_2d_setup.log", &
                  exist=found)
         if (found) then
            open (newunit=d2d_listing_get_unit, &
                  file="decomp_2d_setup.log", &
                  status="old", &
                  position="append", &
                  iostat=ierror)
         else
            open (newunit=d2d_listing_get_unit, &
                  file="decomp_2d_setup.log", &
                  status="new", &
                  iostat=ierror)
         end if
      else
         d2d_listing_get_unit = output_unit
         ierror = 0
      end if
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "Could not open log file")

   end function d2d_listing_get_unit

   !
   ! Close the IO unit for the log if needed
   !
   module subroutine d2d_listing_close_unit(io_unit)

      use iso_fortran_env, only: output_unit

      implicit none

      ! Input
      integer, intent(in) :: io_unit

      ! Local variables
      integer :: ierror

      !
      ! Close the IO unit if it was not stdout
      !
      if (io_unit /= output_unit) then
         close (io_unit, iostat=ierror)
         if (ierror /= 0) call decomp_2d_abort(ierror, "Could not close log file")
      end if

   end subroutine d2d_listing_close_unit

   !
   ! Print some information about decomp_2d
   !
   module subroutine d2d_listing(given_io_unit)

      use iso_fortran_env, only: output_unit, compiler_version, compiler_options

      implicit none

      ! Argument
      integer, intent(in), optional :: given_io_unit

      ! Local variable
      integer :: io_unit
      integer :: version, subversion, ierror
#ifdef DEBUG
      character(len=512) :: fname
#endif

      ! Output log if needed
      if (decomp_log == D2D_LOG_QUIET) return
      if (decomp_log == D2D_LOG_STDOUT .and. nrank /= 0) return
      if (decomp_log == D2D_LOG_TOFILE .and. nrank /= 0) return

      ! If no IO unit provided, use stdout
      if (present(given_io_unit)) then
         io_unit = given_io_unit
      else
         io_unit = output_unit
      end if

      ! Header
      write (io_unit, *) '==========================================================='
      write (io_unit, *) '=================== Decomp2D - log ========================'
      write (io_unit, *) '==========================================================='

      ! Major and minor version number
      if (D2D_RELEASE) then
         write (io_unit, "(A,I0,A,I0)") ' Release ', D2D_MAJOR, '.', D2D_MINOR
      else
         write (io_unit, "(A,I0,A,I0,A)") ' Release ', D2D_MAJOR, '.', D2D_MINOR, '.alpha'
      end if

      ! Git hash if available
#if defined(VERSION)
      write (io_unit, *) 'Git version        : ', VERSION
#else
      write (io_unit, *) 'Git version        : unknown'
#endif

      ! Basic info
#ifdef DEBUG
      if (decomp_debug >= D2D_DEBUG_LEVEL_INFO) &
         write (io_unit, *) 'I am mpi rank ', nrank
#endif
      write (io_unit, *) 'Total ranks ', nproc
      write (io_unit, *) 'Global data size : ', nx_global, ny_global, nz_global
      write (io_unit, *) 'p_row, p_col : ', dims(1), dims(2)
      write (io_unit, *) 'Periodicity : ', periodic_x, periodic_y, periodic_z
      write (io_unit, *) 'Number of bytes / float number : ', mytype_bytes
      write (io_unit, *) '==========================================================='

      ! Show detected flags, compiler options, version of the MPI library
      write (io_unit, *) 'Compile flags detected :'
#ifdef DOUBLE_PREC
#ifdef SAVE_SINGLE
      write (io_unit, *) 'Numerical precision: Double, saving in single'
#else
      write (io_unit, *) 'Numerical precision: Double'
#endif
#else
      write (io_unit, *) 'Numerical precision: Single'
#endif
      write (io_unit, *) 'Compiled with ', compiler_version()
      write (io_unit, *) 'Compiler options : ', compiler_options()
      call MPI_Get_version(version, subversion, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_Get_version")
      write (io_unit, '(" Version of the MPI library : ",I0,".",I0)') version, subversion
#ifdef DEBUG
      write (io_unit, *) 'Compile flag DEBUG detected'
      write (io_unit, *) '   debug level : ', decomp_debug
#endif
#ifdef PROFILER
      write (io_unit, *) 'Compile flag PROFILER detected'
#endif
#ifdef EVEN
      write (io_unit, *) 'Compile flag EVEN detected'
#endif
#ifdef OVERWRITE
      write (io_unit, *) 'Compile flag OVERWRITE detected'
#endif
#ifdef HALO_DEBUG
      write (io_unit, *) 'Compile flag HALO_DEBUG detected'
#endif
#ifdef _GPU
      write (io_unit, *) 'Compile flag _GPU detected'
#endif
#ifdef _NCCL
      write (io_unit, *) 'Compile flag _NCCL detected'
#endif
      write (io_unit, *) '==========================================================='
      ! Info about each decomp_info object
      call decomp_info_print(decomp_main, io_unit, "decomp_main")
      call decomp_info_print(phG, io_unit, "phG")
      call decomp_info_print(ph1, io_unit, "ph1")
      call decomp_info_print(ph2, io_unit, "ph2")
      call decomp_info_print(ph3, io_unit, "ph3")
      call decomp_info_print(ph4, io_unit, "ph4")
      write (io_unit, *) '==========================================================='
      write (io_unit, *) '==========================================================='

#ifdef DEBUG
      !
      ! In DEBUG mode, rank 0 will also print environment variables
      !
      ! At high debug level, all ranks will print env. variables
      !
      ! The system call, if writing to a file, is not blocking if supported
      !
      if (nrank == 0 .or. decomp_debug >= D2D_DEBUG_LEVEL_INFO) then
         write (io_unit, *) '============== Environment variables ======================'
         write (io_unit, *) '==========================================================='
         write (io_unit, *) '==========================================================='
         if (io_unit == output_unit) then
            call execute_command_line("env", wait=.true.)
         else
            inquire (unit=io_unit, name=fname, iostat=ierror)
            if (ierror /= 0) then
               call decomp_2d_abort(__FILE__, __LINE__, ierror, "No name for the log file")
            end if
            ! Close the IO unit to print the environment variables
            call d2d_listing_close_unit(io_unit)
            call execute_command_line("env >> "//trim(fname), wait=.true.)
         end if
      end if
#else
      ! Close the IO unit if needed
      call d2d_listing_close_unit(io_unit)
#endif

   end subroutine d2d_listing

   !
   ! Print some information about given decomp_info object
   !
   module subroutine decomp_info_print(d2d, io_unit, d2dname)

      implicit none

      ! Arguments
      type(decomp_info), intent(in) :: d2d
      integer, intent(in) :: io_unit
      character(len=*), intent(in) :: d2dname

      ! Nothing to print if not initialized
      if (.not. allocated(d2d%x1dist)) then
         write (io_unit, *) 'Uninitialized decomp_info ', d2dname
         return
      end if

      !
      ! If DEBUG mode, print everything
      ! Otherwise, print only global size
      !
      write (io_unit, *) 'Decomp_info : ', d2dname
      write (io_unit, *) '   Global size : ', d2d%xsz(1), d2d%ysz(2), d2d%zsz(3)
#ifdef DEBUG
      write (io_unit, *) '   xsz, xst, xen : ', d2d%xsz, d2d%xst, d2d%xen
      write (io_unit, *) '   ysz, yst, yen : ', d2d%ysz, d2d%yst, d2d%yen
      write (io_unit, *) '   zsz, zst, zen : ', d2d%zsz, d2d%zst, d2d%zen
      write (io_unit, *) '   x1dist : ', d2d%x1dist
      write (io_unit, *) '   y1dist : ', d2d%y1dist
      write (io_unit, *) '   y2dist : ', d2d%y2dist
      write (io_unit, *) '   z2dist : ', d2d%z2dist
      write (io_unit, *) '   x1cnts : ', d2d%x1cnts
      write (io_unit, *) '   y1cnts : ', d2d%y1cnts
      write (io_unit, *) '   y2cnts : ', d2d%y2cnts
      write (io_unit, *) '   z2cnts : ', d2d%z2cnts
      write (io_unit, *) '   x1disp : ', d2d%x1disp
      write (io_unit, *) '   y1disp : ', d2d%y1disp
      write (io_unit, *) '   y2disp : ', d2d%y2disp
      write (io_unit, *) '   z2disp : ', d2d%z2disp
#ifdef EVEN
      write (io_unit, *) '   x1count : ', d2d%x1count
      write (io_unit, *) '   y1count : ', d2d%y1count
      write (io_unit, *) '   y2count : ', d2d%y2count
      write (io_unit, *) '   z2count : ', d2d%z2count
      write (io_unit, *) '   even : ', d2d%even
#endif
#endif

   end subroutine decomp_info_print

end submodule d2d_log
