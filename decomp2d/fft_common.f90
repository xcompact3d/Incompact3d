!! SPDX-License-Identifier: BSD-3-Clause

! This file contains common code shared by all FFT engines

integer, save :: format                 ! input X-pencil or Z-pencil

! The libary can only be initialised once
logical, save :: initialised = .false.

! Global size of the FFT
integer, save :: nx_fft, ny_fft, nz_fft

! 2D processor grid
! FIXME this is already available in the module decomp_2d
integer, save, dimension(2) :: dims

! Decomposition objects
TYPE(DECOMP_INFO), pointer, save :: ph => null()  ! physical space
TYPE(DECOMP_INFO), target, save :: sp  ! spectral space

! Workspace to store the intermediate Y-pencil data
complex(mytype), allocatable, target, dimension(:, :, :) :: wk2_c2c
complex(mytype), contiguous, pointer, dimension(:, :, :) :: wk2_r2c
! Workspace for r2c and c2r transforms
! FIXME could be removed using in-place r2c and c2r ?
complex(mytype), allocatable, dimension(:, :, :) :: wk13

public :: decomp_2d_fft_init, decomp_2d_fft_3d, &
          decomp_2d_fft_finalize, decomp_2d_fft_get_size, &
          decomp_2d_fft_get_ph, decomp_2d_fft_get_sp

! Declare generic interfaces to handle different inputs

interface decomp_2d_fft_init
   module procedure fft_init_noarg
   module procedure fft_init_arg
   module procedure fft_init_general
end interface

interface decomp_2d_fft_3d
   module procedure fft_3d_c2c
   module procedure fft_3d_r2c
   module procedure fft_3d_c2r
end interface

interface
   module subroutine decomp_2d_fft_log(backend)
      character(len=*), intent(in) :: backend
   end subroutine decomp_2d_fft_log
end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialise the FFT module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_init_noarg

   implicit none

   call fft_init_arg(PHYSICAL_IN_X)  ! default input is X-pencil data

   return
end subroutine fft_init_noarg

subroutine fft_init_arg(pencil)     ! allow to handle Z-pencil input

   implicit none

   integer, intent(IN) :: pencil

   call fft_init_general(pencil, nx_global, ny_global, nz_global)

   return
end subroutine fft_init_arg

! Initialise the FFT library to perform arbitrary size transforms
subroutine fft_init_general(pencil, nx, ny, nz)

   implicit none

   integer, intent(IN) :: pencil
   integer, intent(IN) :: nx, ny, nz

   integer :: status, errorcode

#ifdef PROFILER
   if (decomp_profiler_fft) call decomp_profiler_start("fft_init")
#endif

   ! Safety checks
   if (initialised) then
      errorcode = 4
      call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                           'FFT library should only be initialised once')
   end if
   if (nx <= 0) call decomp_2d_abort(__FILE__, __LINE__, nx, "Invalid value for nx")
   if (ny <= 0) call decomp_2d_abort(__FILE__, __LINE__, ny, "Invalid value for ny")
   if (nz <= 0) call decomp_2d_abort(__FILE__, __LINE__, nz, "Invalid value for nz")

   format = pencil
   nx_fft = nx
   ny_fft = ny
   nz_fft = nz

! determine the processor grid in use
   dims = get_decomp_dims()

! for c2r/r2c interface:
! if in physical space, a real array is of size: nx*ny*nz
! in spectral space, the complex array is of size:
!         (nx/2+1)*ny*nz, if PHYSICAL_IN_X
!      or nx*ny*(nz/2+1), if PHYSICAL_IN_Z

   if (nx_fft == nx_global .and. ny_fft == ny_global .and. nz_fft == nz_global) then
      ph => decomp_main
   else
      if (.not. associated(ph)) allocate (ph)
      call decomp_info_init(nx, ny, nz, ph)
   end if
   if (format == PHYSICAL_IN_X) then
      call decomp_info_init(nx / 2 + 1, ny, nz, sp)
   else if (format == PHYSICAL_IN_Z) then
      call decomp_info_init(nx, ny, nz / 2 + 1, sp)
   else
      call decomp_2d_abort(__FILE__, __LINE__, format, "Invalid value for format")
   end if

   !
   ! Allocate the workspace for intermediate y-pencil data
   ! The largest memory block needed is the one for c2c transforms
   !
   allocate (wk2_c2c(ph%ysz(1), ph%ysz(2), ph%ysz(3)), STAT=status)
   if (status /= 0) then
      errorcode = 3
      call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                           'Out of memory when initialising FFT')
   end if
   !
   ! A smaller memory block is needed for r2c and c2r transforms
   ! wk2_c2c and wk2_r2c start at the same memory location
   !
   !    Size of wk2_c2c : ph%ysz(1), ph%ysz(2), ph%ysz(3)
   !    Size of wk2_r2c : sp%ysz(1), sp%ysz(2), sp%ysz(3)
   !
   call c_f_pointer(c_loc(wk2_c2c), wk2_r2c, sp%ysz)
   !
   ! Allocate the workspace for r2c and c2r transforms
   !
   ! wk13 can not be easily fused with wk2_*2c due to statements such as
   ! transpose_y_to_x(wk2_r2c, wk13, sp)
   ! transpose_y_to_z(wk2_r2c, wk13, sp)
   !
   if (format == PHYSICAL_IN_X) then
      allocate (wk13(sp%xsz(1), sp%xsz(2), sp%xsz(3)), STAT=status)
   else if (format == PHYSICAL_IN_Z) then
      allocate (wk13(sp%zsz(1), sp%zsz(2), sp%zsz(3)), STAT=status)
   end if
   if (status /= 0) then
      errorcode = 3
      call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                           'Out of memory when initialising FFT')
   end if

   call init_fft_engine

   initialised = .true.

#ifdef PROFILER
   if (decomp_profiler_fft) call decomp_profiler_end("fft_init")
#endif

   return
end subroutine fft_init_general

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Final clean up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine decomp_2d_fft_finalize

   implicit none

#ifdef PROFILER
   if (decomp_profiler_fft) call decomp_profiler_start("fft_fin")
#endif

   if (nx_fft /= nx_global .or. ny_fft /= ny_global .or. nz_fft /= nz_global) then
      call decomp_info_finalize(ph)
      deallocate (ph)
   end if
   nullify (ph)
   call decomp_info_finalize(sp)

   if (allocated(wk2_c2c)) deallocate (wk2_c2c)
   if (associated(wk2_r2c)) nullify (wk2_r2c)
   if (allocated(wk13)) deallocate (wk13)

   call finalize_fft_engine

   initialised = .false.

#ifdef PROFILER
   if (decomp_profiler_fft) call decomp_profiler_end("fft_fin")
#endif

   return
end subroutine decomp_2d_fft_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Return the size, starting/ending index of the distributed array
!  whose global size is (nx/2+1)*ny*nz, for defining data structures
!  in r2c and c2r interfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine decomp_2d_fft_get_size(istart, iend, isize)

   implicit none
   integer, dimension(3), intent(OUT) :: istart, iend, isize

   if (format == PHYSICAL_IN_X) then
      istart = sp%zst
      iend = sp%zen
      isize = sp%zsz
   else if (format == PHYSICAL_IN_Z) then
      istart = sp%xst
      iend = sp%xen
      isize = sp%xsz
   end if

   return
end subroutine decomp_2d_fft_get_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Return a pointer to the decomp_info object ph
!
! The caller should not apply decomp_info_finalize on the pointer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function decomp_2d_fft_get_ph()

   implicit none

   type(decomp_info), pointer :: decomp_2d_fft_get_ph

   if (.not. associated(ph)) then
      call decomp_2d_abort(__FILE__, __LINE__, -1, 'FFT library must be initialised first')
   end if
   decomp_2d_fft_get_ph => ph

end function decomp_2d_fft_get_ph

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Return a pointer to the decomp_info object sp
!
! The caller should not apply decomp_info_finalize on the pointer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function decomp_2d_fft_get_sp()

   implicit none

   type(decomp_info), pointer :: decomp_2d_fft_get_sp

   if (.not. associated(ph)) then
      call decomp_2d_abort(__FILE__, __LINE__, -1, 'FFT library must be initialised first')
   end if
   decomp_2d_fft_get_sp => sp

end function decomp_2d_fft_get_sp
