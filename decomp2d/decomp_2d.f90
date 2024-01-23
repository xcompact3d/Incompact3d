!! SPDX-License-Identifier: BSD-3-Clause

! This is the main 2D pencil decomposition module

module decomp_2d

   use MPI
   use, intrinsic :: iso_fortran_env, only: real32, real64
   use factor
   use decomp_2d_constants
   use decomp_2d_mpi
#if defined(_GPU)
   use cudafor
   use decomp_2d_cumpi
#if defined(_NCCL)
   use nccl
   use decomp_2d_nccl
#endif
#endif

   implicit none

   private        ! Make everything private unless declared public

   ! some key global variables
   integer, save, public :: nx_global, ny_global, nz_global  ! global size

   ! parameters for 2D Cartesian topology
   integer, save, dimension(2) :: dims, coord
   integer, save, public :: DECOMP_2D_COMM_CART_X = MPI_COMM_NULL
   integer, save, public :: DECOMP_2D_COMM_CART_Y = MPI_COMM_NULL
   integer, save, public :: DECOMP_2D_COMM_CART_Z = MPI_COMM_NULL
   integer, save :: DECOMP_2D_COMM_ROW = MPI_COMM_NULL
   integer, save :: DECOMP_2D_COMM_COL = MPI_COMM_NULL

   ! define neighboring blocks (to be used in halo-cell support)
   !  first dimension 1=X-pencil, 2=Y-pencil, 3=Z-pencil
   ! second dimension 1=east, 2=west, 3=north, 4=south, 5=top, 6=bottom
   integer, save, dimension(3, 6) :: neighbour

   ! flags for periodic condition in three dimensions
   logical, save :: periodic_x, periodic_y, periodic_z

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
#ifdef DEBUG
   integer, public, save :: decomp_log = D2D_LOG_TOFILE_FULL
#else
   integer, public, save :: decomp_log = D2D_LOG_TOFILE
#endif

   !
   ! Debug level can be changed by the external code before calling decomp_2d_init
   !
   ! The environment variable "DECOMP_2D_DEBUG" can be used to change the debug level
   !
   ! Debug checks are performed only when the preprocessor variable DEBUG is defined
   !
#ifdef DEBUG
   integer(kind(D2D_DEBUG_LEVEL_OFF)), public, save :: decomp_debug = D2D_DEBUG_LEVEL_INFO
#else
   integer(kind(D2D_DEBUG_LEVEL_OFF)), public, save :: decomp_debug = D2D_DEBUG_LEVEL_OFF
#endif

   ! derived type to store decomposition info for a given global data size
   TYPE, public :: DECOMP_INFO
      ! staring/ending index and size of data held by current processor
      integer, dimension(3) :: xst, xen, xsz  ! x-pencil
      integer, dimension(3) :: yst, yen, ysz  ! y-pencil
      integer, dimension(3) :: zst, zen, zsz  ! z-pencil

      ! in addition to local information, processors also need to know
      ! some global information for global communications to work

      ! how each dimension is distributed along pencils
      integer, allocatable, dimension(:) :: &
         x1dist, y1dist, y2dist, z2dist

      ! send/receive buffer counts and displacements for MPI_ALLTOALLV
      integer, allocatable, dimension(:) :: &
         x1cnts, y1cnts, y2cnts, z2cnts
      integer, allocatable, dimension(:) :: &
         x1disp, y1disp, y2disp, z2disp

#ifdef EVEN
      ! buffer counts for MPI_ALLTOALL for padded-alltoall
      integer :: x1count, y1count, y2count, z2count
      ! evenly distributed data
      logical :: even
#endif

   END TYPE DECOMP_INFO

   ! main (default) decomposition information for global size nx*ny*nz
   TYPE(DECOMP_INFO), target, save, public :: decomp_main
   ! FIXME The extra decomp_info objects should be defined in the external code, not here
   !       Currently keeping them to avoid breaking external codes
   TYPE(DECOMP_INFO), save, public :: phG, ph1, ph2, ph3, ph4

   ! staring/ending index and size of data held by current processor
   ! duplicate 'decomp_main', needed by apps to define data structure
   integer, save, dimension(3), public :: xstart, xend, xsize  ! x-pencil
   integer, save, dimension(3), public :: ystart, yend, ysize  ! y-pencil
   integer, save, dimension(3), public :: zstart, zend, zsize  ! z-pencil

   ! These are the buffers used by MPI_ALLTOALL(V) calls
   integer, save :: decomp_buf_size = 0
   ! Shared real/complex buffers
   real(mytype), target, allocatable, dimension(:) :: work1, work2
   ! Real/complex pointers to buffers
   real(mytype), pointer, contiguous, dimension(:) :: work1_r, work2_r
   complex(mytype), pointer, contiguous, dimension(:) :: work1_c, work2_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! To define smaller arrays using every several mesh points
   integer, save, dimension(3), public :: xszS, yszS, zszS, xstS, ystS, zstS, xenS, yenS, zenS
   integer, save, dimension(3), public :: xszV, yszV, zszV, xstV, ystV, zstV, xenV, yenV, zenV
   integer, save, dimension(3), public :: xszP, yszP, zszP, xstP, ystP, zstP, xenP, yenP, zenP
   logical, save :: coarse_mesh_starts_from_1
   integer, save :: iskipS, jskipS, kskipS
   integer, save :: iskipV, jskipV, kskipV
   integer, save :: iskipP, jskipP, kskipP

   ! public user routines
   public :: decomp_2d_init, decomp_2d_finalize, &
             transpose_x_to_y, transpose_y_to_z, &
             transpose_z_to_y, transpose_y_to_x, &
             decomp_info_init, decomp_info_finalize, partition, &
             decomp_info_print, &
             init_coarser_mesh_statS, fine_to_coarseS, &
             init_coarser_mesh_statV, fine_to_coarseV, &
             init_coarser_mesh_statP, fine_to_coarseP, &
             alloc_x, alloc_y, alloc_z, &
             update_halo, &
             get_decomp_info, &
             get_decomp_dims, &
             d2d_listing_get_unit, d2d_listing_close_unit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! These are routines to perform global data transpositions
   !
   !   Four combinations are available, enough to cover all situations
   !    - transpose_x_to_y (X-pencil --> Y-pencil)
   !    - transpose_y_to_z (Y-pencil --> Z-pencil)
   !    - transpose_z_to_y (Z-pencil --> Y-pencil)
   !    - transpose_y_to_x (Y-pencil --> X-pencil)
   !
   !   Generic interface provided here to support multiple data types
   !    - real and complex types supported through generic interface
   !    - single/double precision supported through pre-processing
   !       * see 'mytype' variable at the beginning
   !    - an optional argument can be supplied to transpose data whose
   !      global size is not the default nx*ny*nz
   !       * as the case in fft r2c/c2r interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   interface decomp_2d_init
      module procedure decomp_2d_init_ref
   end interface decomp_2d_init

   interface decomp_2d_finalize
      module procedure decomp_2d_finalize_ref
   end interface decomp_2d_finalize

   interface transpose_x_to_y
      module procedure transpose_x_to_y_real_long
      module procedure transpose_x_to_y_real_short
      module procedure transpose_x_to_y_complex_long
      module procedure transpose_x_to_y_complex_short
   end interface transpose_x_to_y

   interface transpose_y_to_z
      module procedure transpose_y_to_z_real_long
      module procedure transpose_y_to_z_real_short
      module procedure transpose_y_to_z_complex_long
      module procedure transpose_y_to_z_complex_short
   end interface transpose_y_to_z

   interface transpose_z_to_y
      module procedure transpose_z_to_y_real_long
      module procedure transpose_z_to_y_real_short
      module procedure transpose_z_to_y_complex_long
      module procedure transpose_z_to_y_complex_short
   end interface transpose_z_to_y

   interface transpose_y_to_x
      module procedure transpose_y_to_x_real_long
      module procedure transpose_y_to_x_real_short
      module procedure transpose_y_to_x_complex_long
      module procedure transpose_y_to_x_complex_short
   end interface transpose_y_to_x

   interface update_halo
      module procedure update_halo_real
      module procedure update_halo_real_short
      module procedure update_halo_complex
      module procedure update_halo_complex_short
   end interface update_halo

   interface alloc_x
      module procedure alloc_x_real
      module procedure alloc_x_real_short
      module procedure alloc_x_complex
      module procedure alloc_x_complex_short
   end interface alloc_x

   interface alloc_y
      module procedure alloc_y_real
      module procedure alloc_y_real_short
      module procedure alloc_y_complex
      module procedure alloc_y_complex_short
   end interface alloc_y

   interface alloc_z
      module procedure alloc_z_real
      module procedure alloc_z_real_short
      module procedure alloc_z_complex
      module procedure alloc_z_complex_short
   end interface alloc_z

   interface

      module function d2d_listing_get_unit()
         integer :: d2d_listing_get_unit
      end function d2d_listing_get_unit

      module subroutine d2d_listing_close_unit(io_unit)
         integer, intent(in) :: io_unit
      end subroutine d2d_listing_close_unit

      module subroutine d2d_listing(given_io_unit)
         integer, intent(in), optional :: given_io_unit
      end subroutine d2d_listing

      module subroutine decomp_info_print(d2d, io_unit, d2dname)
         type(decomp_info), intent(in) :: d2d
         integer, intent(in) :: io_unit
         character(len=*), intent(in) :: d2dname
      end subroutine decomp_info_print

   end interface

contains

#include "decomp_2d_init_fin.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Return the default decomposition object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! FIXME avoid a copy and return a pointer to decomp_main
   ! TODO list the external codes using this subroutine
   subroutine get_decomp_info(decomp)

      implicit none

      ! FIXME TYPE(DECOMP_INFO), pointer :: decomp
      TYPE(DECOMP_INFO), intent(OUT) :: decomp

      ! FIXME decomp => decomp_main
      decomp = decomp_main

      return
   end subroutine get_decomp_info

   !
   ! Return the 2D processor grid
   !
   function get_decomp_dims()

      implicit none

      integer, dimension(2) :: get_decomp_dims

      get_decomp_dims = dims

   end function get_decomp_dims

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Advanced Interface allowing applications to define globle domain of
   ! any size, distribute it, and then transpose data among pencils.
   !  - generate 2D decomposition details as defined in DECOMP_INFO
   !  - the default global data size is nx*ny*nz
   !  - a different global size nx/2+1,ny,nz is used in FFT r2c/c2r
   !  - multiple global sizes can co-exist in one application, each
   !    using its own DECOMP_INFO object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine decomp_info_init(nx, ny, nz, decomp)

      use, intrinsic:: iso_c_binding, only: c_f_pointer, c_loc

      implicit none

      integer, intent(IN) :: nx, ny, nz
      TYPE(DECOMP_INFO), intent(INOUT) :: decomp

      integer :: buf_size, status, errorcode

      ! verify the global size can actually be distributed as pencils
      if (nx_global < dims(1) .or. ny_global < dims(1) .or. ny_global < dims(2) .or. nz_global < dims(2)) then
         errorcode = 6
         call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                              'Invalid 2D processor grid. '// &
                              'Make sure that min(nx,ny) >= p_row and '// &
                              'min(ny,nz) >= p_col')
      end if

      ! distribute mesh points
      allocate (decomp%x1dist(0:dims(1) - 1), decomp%y1dist(0:dims(1) - 1), &
                decomp%y2dist(0:dims(2) - 1), decomp%z2dist(0:dims(2) - 1))
      call get_dist(nx, ny, nz, decomp)

      ! generate partition information - starting/ending index etc.
      call partition(nx, ny, nz, (/1, 2, 3/), &
                     decomp%xst, decomp%xen, decomp%xsz)
      call partition(nx, ny, nz, (/2, 1, 3/), &
                     decomp%yst, decomp%yen, decomp%ysz)
      call partition(nx, ny, nz, (/2, 3, 1/), &
                     decomp%zst, decomp%zen, decomp%zsz)

      ! prepare send/receive buffer displacement and count for ALLTOALL(V)
      allocate (decomp%x1cnts(0:dims(1) - 1), decomp%y1cnts(0:dims(1) - 1), &
                decomp%y2cnts(0:dims(2) - 1), decomp%z2cnts(0:dims(2) - 1))
      allocate (decomp%x1disp(0:dims(1) - 1), decomp%y1disp(0:dims(1) - 1), &
                decomp%y2disp(0:dims(2) - 1), decomp%z2disp(0:dims(2) - 1))
      call prepare_buffer(decomp)

      ! allocate memory for the MPI_ALLTOALL(V) buffers
      ! define the buffers globally for performance reason

      buf_size = max(decomp%xsz(1) * decomp%xsz(2) * decomp%xsz(3), &
                     max(decomp%ysz(1) * decomp%ysz(2) * decomp%ysz(3), &
                         decomp%zsz(1) * decomp%zsz(2) * decomp%zsz(3)))

#ifdef EVEN
      ! padded alltoall optimisation may need larger buffer space
      buf_size = max(buf_size, &
                     max(decomp%x1count * dims(1), decomp%y2count * dims(2)))
      ! evenly distributed data ?
      if (mod(nx, dims(1)) == 0 .and. mod(ny, dims(1)) == 0 .and. &
          mod(ny, dims(2)) == 0 .and. mod(nz, dims(2)) == 0) then
         decomp%even = .true.
      else
         decomp%even = .false.
      end if
#endif

      ! check if additional memory is required
      if (buf_size > decomp_buf_size) then
         decomp_buf_size = buf_size
#if defined(_GPU)
         call decomp_2d_cumpi_init(buf_size)
#endif
         if (associated(work1_r)) nullify (work1_r)
         if (associated(work2_r)) nullify (work2_r)
         if (associated(work1_c)) nullify (work1_c)
         if (associated(work2_c)) nullify (work2_c)
         if (allocated(work1)) deallocate (work1)
         if (allocated(work2)) deallocate (work2)
         allocate (work1(2 * buf_size), STAT=status)
         if (status /= 0) then
            errorcode = 2
            call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                                 'Out of memory when allocating 2DECOMP workspace')
         end if
         allocate (work2(2 * buf_size), STAT=status)
         if (status /= 0) then
            errorcode = 2
            call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                                 'Out of memory when allocating 2DECOMP workspace')
         end if
         call c_f_pointer(c_loc(work1), work1_r, [buf_size])
         call c_f_pointer(c_loc(work2), work2_r, [buf_size])
         call c_f_pointer(c_loc(work1), work1_c, [buf_size])
         call c_f_pointer(c_loc(work2), work2_c, [buf_size])
      end if

   end subroutine decomp_info_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Release memory associated with a DECOMP_INFO object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine decomp_info_finalize(decomp)

      implicit none

      TYPE(DECOMP_INFO), intent(INOUT) :: decomp

      if (allocated(decomp%x1dist)) deallocate (decomp%x1dist)
      if (allocated(decomp%y1dist)) deallocate (decomp%y1dist)
      if (allocated(decomp%y2dist)) deallocate (decomp%y2dist)
      if (allocated(decomp%z2dist)) deallocate (decomp%z2dist)
      if (allocated(decomp%x1cnts)) deallocate (decomp%x1cnts)
      if (allocated(decomp%y1cnts)) deallocate (decomp%y1cnts)
      if (allocated(decomp%y2cnts)) deallocate (decomp%y2cnts)
      if (allocated(decomp%z2cnts)) deallocate (decomp%z2cnts)
      if (allocated(decomp%x1disp)) deallocate (decomp%x1disp)
      if (allocated(decomp%y1disp)) deallocate (decomp%y1disp)
      if (allocated(decomp%y2disp)) deallocate (decomp%y2disp)
      if (allocated(decomp%z2disp)) deallocate (decomp%z2disp)

      return
   end subroutine decomp_info_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Coarser mesh support for statistic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine init_coarser_mesh_statS(i_skip, j_skip, k_skip, from1)

      implicit none

      integer, intent(IN) :: i_skip, j_skip, k_skip
      logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
      ! .false. - save n,2n,3n...

      integer, dimension(3) :: skip
      integer :: i

      coarse_mesh_starts_from_1 = from1
      iskipS = i_skip
      jskipS = j_skip
      kskipS = k_skip

      skip(1) = iskipS
      skip(2) = jskipS
      skip(3) = kskipS

      do i = 1, 3
         if (from1) then
            xstS(i) = (xstart(i) + skip(i) - 1) / skip(i)
            if (mod(xstart(i) + skip(i) - 1, skip(i)) /= 0) xstS(i) = xstS(i) + 1
            xenS(i) = (xend(i) + skip(i) - 1) / skip(i)
         else
            xstS(i) = xstart(i) / skip(i)
            if (mod(xstart(i), skip(i)) /= 0) xstS(i) = xstS(i) + 1
            xenS(i) = xend(i) / skip(i)
         end if
         xszS(i) = xenS(i) - xstS(i) + 1
      end do

      do i = 1, 3
         if (from1) then
            ystS(i) = (ystart(i) + skip(i) - 1) / skip(i)
            if (mod(ystart(i) + skip(i) - 1, skip(i)) /= 0) ystS(i) = ystS(i) + 1
            yenS(i) = (yend(i) + skip(i) - 1) / skip(i)
         else
            ystS(i) = ystart(i) / skip(i)
            if (mod(ystart(i), skip(i)) /= 0) ystS(i) = ystS(i) + 1
            yenS(i) = yend(i) / skip(i)
         end if
         yszS(i) = yenS(i) - ystS(i) + 1
      end do

      do i = 1, 3
         if (from1) then
            zstS(i) = (zstart(i) + skip(i) - 1) / skip(i)
            if (mod(zstart(i) + skip(i) - 1, skip(i)) /= 0) zstS(i) = zstS(i) + 1
            zenS(i) = (zend(i) + skip(i) - 1) / skip(i)
         else
            zstS(i) = zstart(i) / skip(i)
            if (mod(zstart(i), skip(i)) /= 0) zstS(i) = zstS(i) + 1
            zenS(i) = zend(i) / skip(i)
         end if
         zszS(i) = zenS(i) - zstS(i) + 1
      end do

      return
   end subroutine init_coarser_mesh_statS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Coarser mesh support for visualization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine init_coarser_mesh_statV(i_skip, j_skip, k_skip, from1)

      implicit none

      integer, intent(IN) :: i_skip, j_skip, k_skip
      logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
      ! .false. - save n,2n,3n...

      integer, dimension(3) :: skip
      integer :: i

      coarse_mesh_starts_from_1 = from1
      iskipV = i_skip
      jskipV = j_skip
      kskipV = k_skip

      skip(1) = iskipV
      skip(2) = jskipV
      skip(3) = kskipV

      do i = 1, 3
         if (from1) then
            xstV(i) = (xstart(i) + skip(i) - 1) / skip(i)
            if (mod(xstart(i) + skip(i) - 1, skip(i)) /= 0) xstV(i) = xstV(i) + 1
            xenV(i) = (xend(i) + skip(i) - 1) / skip(i)
         else
            xstV(i) = xstart(i) / skip(i)
            if (mod(xstart(i), skip(i)) /= 0) xstV(i) = xstV(i) + 1
            xenV(i) = xend(i) / skip(i)
         end if
         xszV(i) = xenV(i) - xstV(i) + 1
      end do

      do i = 1, 3
         if (from1) then
            ystV(i) = (ystart(i) + skip(i) - 1) / skip(i)
            if (mod(ystart(i) + skip(i) - 1, skip(i)) /= 0) ystV(i) = ystV(i) + 1
            yenV(i) = (yend(i) + skip(i) - 1) / skip(i)
         else
            ystV(i) = ystart(i) / skip(i)
            if (mod(ystart(i), skip(i)) /= 0) ystV(i) = ystV(i) + 1
            yenV(i) = yend(i) / skip(i)
         end if
         yszV(i) = yenV(i) - ystV(i) + 1
      end do

      do i = 1, 3
         if (from1) then
            zstV(i) = (zstart(i) + skip(i) - 1) / skip(i)
            if (mod(zstart(i) + skip(i) - 1, skip(i)) /= 0) zstV(i) = zstV(i) + 1
            zenV(i) = (zend(i) + skip(i) - 1) / skip(i)
         else
            zstV(i) = zstart(i) / skip(i)
            if (mod(zstart(i), skip(i)) /= 0) zstV(i) = zstV(i) + 1
            zenV(i) = zend(i) / skip(i)
         end if
         zszV(i) = zenV(i) - zstV(i) + 1
      end do

      return
   end subroutine init_coarser_mesh_statV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Coarser mesh support for probe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine init_coarser_mesh_statP(i_skip, j_skip, k_skip, from1)

      implicit none

      integer, intent(IN) :: i_skip, j_skip, k_skip
      logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
      ! .false. - save n,2n,3n...

      integer, dimension(3) :: skip
      integer :: i

      coarse_mesh_starts_from_1 = from1
      iskipP = i_skip
      jskipP = j_skip
      kskipP = k_skip

      skip(1) = iskipP
      skip(2) = jskipP
      skip(3) = kskipP

      do i = 1, 3
         if (from1) then
            xstP(i) = (xstart(i) + skip(i) - 1) / skip(i)
            if (mod(xstart(i) + skip(i) - 1, skip(i)) /= 0) xstP(i) = xstP(i) + 1
            xenP(i) = (xend(i) + skip(i) - 1) / skip(i)
         else
            xstP(i) = xstart(i) / skip(i)
            if (mod(xstart(i), skip(i)) /= 0) xstP(i) = xstP(i) + 1
            xenP(i) = xend(i) / skip(i)
         end if
         xszP(i) = xenP(i) - xstP(i) + 1
      end do

      do i = 1, 3
         if (from1) then
            ystP(i) = (ystart(i) + skip(i) - 1) / skip(i)
            if (mod(ystart(i) + skip(i) - 1, skip(i)) /= 0) ystP(i) = ystP(i) + 1
            yenP(i) = (yend(i) + skip(i) - 1) / skip(i)
         else
            ystP(i) = ystart(i) / skip(i)
            if (mod(ystart(i), skip(i)) /= 0) ystP(i) = ystP(i) + 1
            yenP(i) = yend(i) / skip(i)
         end if
         yszP(i) = yenP(i) - ystP(i) + 1
      end do

      do i = 1, 3
         if (from1) then
            zstP(i) = (zstart(i) + skip(i) - 1) / skip(i)
            if (mod(zstart(i) + skip(i) - 1, skip(i)) /= 0) zstP(i) = zstP(i) + 1
            zenP(i) = (zend(i) + skip(i) - 1) / skip(i)
         else
            zstP(i) = zstart(i) / skip(i)
            if (mod(zstart(i), skip(i)) /= 0) zstP(i) = zstP(i) + 1
            zenP(i) = zend(i) / skip(i)
         end if
         zszP(i) = zenP(i) - zstP(i) + 1
      end do

      return
   end subroutine init_coarser_mesh_statP

   ! Copy data from a fine-resolution array to a coarse one for statistic
   subroutine fine_to_coarseS(ipencil, var_fine, var_coarse)

      implicit none

      real(mytype), dimension(:, :, :) :: var_fine
      real(mytype), dimension(:, :, :) :: var_coarse
      integer, intent(IN) :: ipencil

      real(mytype), allocatable, dimension(:, :, :) :: wk, wk2
      integer :: i, j, k

      if (ipencil == 1) then
         allocate (wk(xstS(1):xenS(1), xstS(2):xenS(2), xstS(3):xenS(3)))
         allocate (wk2(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
         wk2 = var_fine
         if (coarse_mesh_starts_from_1) then
            do k = xstS(3), xenS(3)
               do j = xstS(2), xenS(2)
                  do i = xstS(1), xenS(1)
                     wk(i, j, k) = wk2((i - 1) * iskipS + 1, (j - 1) * jskipS + 1, (k - 1) * kskipS + 1)
                  end do
               end do
            end do
         else
            do k = xstS(3), xenS(3)
               do j = xstS(2), xenS(2)
                  do i = xstS(1), xenS(1)
                     wk(i, j, k) = wk2(i * iskipS, j * jskipS, k * kskipS)
                  end do
               end do
            end do
         end if
         var_coarse = wk
      else if (ipencil == 2) then
         allocate (wk(ystS(1):yenS(1), ystS(2):yenS(2), ystS(3):yenS(3)))
         allocate (wk2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
         wk2 = var_fine
         if (coarse_mesh_starts_from_1) then
            do k = ystS(3), yenS(3)
               do j = ystS(2), yenS(2)
                  do i = ystS(1), yenS(1)
                     wk(i, j, k) = wk2((i - 1) * iskipS + 1, (j - 1) * jskipS + 1, (k - 1) * kskipS + 1)
                  end do
               end do
            end do
         else
            do k = ystS(3), yenS(3)
               do j = ystS(2), yenS(2)
                  do i = ystS(1), yenS(1)
                     wk(i, j, k) = wk2(i * iskipS, j * jskipS, k * kskipS)
                  end do
               end do
            end do
         end if
         var_coarse = wk
      else if (ipencil == 3) then
         allocate (wk(zstS(1):zenS(1), zstS(2):zenS(2), zstS(3):zenS(3)))
         allocate (wk2(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
         wk2 = var_fine
         if (coarse_mesh_starts_from_1) then
            do k = zstS(3), zenS(3)
               do j = zstS(2), zenS(2)
                  do i = zstS(1), zenS(1)
                     wk(i, j, k) = wk2((i - 1) * iskipS + 1, (j - 1) * jskipS + 1, (k - 1) * kskipS + 1)
                  end do
               end do
            end do
         else
            do k = zstS(3), zenS(3)
               do j = zstS(2), zenS(2)
                  do i = zstS(1), zenS(1)
                     wk(i, j, k) = wk2(i * iskipS, j * jskipS, k * kskipS)
                  end do
               end do
            end do
         end if
         var_coarse = wk
      end if

      deallocate (wk, wk2)

      return
   end subroutine fine_to_coarseS

   ! Copy data from a fine-resolution array to a coarse one for visualization
   subroutine fine_to_coarseV(ipencil, var_fine, var_coarse)

      implicit none

      real(mytype), dimension(:, :, :) :: var_fine
      real(mytype), dimension(:, :, :) :: var_coarse
      integer, intent(IN) :: ipencil

      real(mytype), allocatable, dimension(:, :, :) :: wk, wk2
      integer :: i, j, k

      if (ipencil == 1) then
         allocate (wk(xstV(1):xenV(1), xstV(2):xenV(2), xstV(3):xenV(3)))
         allocate (wk2(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
         wk2 = var_fine
         if (coarse_mesh_starts_from_1) then
            do k = xstV(3), xenV(3)
               do j = xstV(2), xenV(2)
                  do i = xstV(1), xenV(1)
                     wk(i, j, k) = wk2((i - 1) * iskipV + 1, (j - 1) * jskipV + 1, (k - 1) * kskipV + 1)
                  end do
               end do
            end do
         else
            do k = xstV(3), xenV(3)
               do j = xstV(2), xenV(2)
                  do i = xstV(1), xenV(1)
                     wk(i, j, k) = wk2(i * iskipV, j * jskipV, k * kskipV)
                  end do
               end do
            end do
         end if
         var_coarse = wk
      else if (ipencil == 2) then
         allocate (wk(ystV(1):yenV(1), ystV(2):yenV(2), ystV(3):yenV(3)))
         allocate (wk2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
         wk2 = var_fine
         if (coarse_mesh_starts_from_1) then
            do k = ystV(3), yenV(3)
               do j = ystV(2), yenV(2)
                  do i = ystV(1), yenV(1)
                     wk(i, j, k) = wk2((i - 1) * iskipV + 1, (j - 1) * jskipV + 1, (k - 1) * kskipV + 1)
                  end do
               end do
            end do
         else
            do k = ystV(3), yenV(3)
               do j = ystV(2), yenV(2)
                  do i = ystV(1), yenV(1)
                     wk(i, j, k) = wk2(i * iskipV, j * jskipV, k * kskipV)
                  end do
               end do
            end do
         end if
         var_coarse = wk
      else if (ipencil == 3) then
         allocate (wk(zstV(1):zenV(1), zstV(2):zenV(2), zstV(3):zenV(3)))
         allocate (wk2(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
         wk2 = var_fine
         if (coarse_mesh_starts_from_1) then
            do k = zstV(3), zenV(3)
               do j = zstV(2), zenV(2)
                  do i = zstV(1), zenV(1)
                     wk(i, j, k) = wk2((i - 1) * iskipV + 1, (j - 1) * jskipV + 1, (k - 1) * kskipV + 1)
                  end do
               end do
            end do
         else
            do k = zstV(3), zenV(3)
               do j = zstV(2), zenV(2)
                  do i = zstV(1), zenV(1)
                     wk(i, j, k) = wk2(i * iskipV, j * jskipV, k * kskipV)
                  end do
               end do
            end do
         end if
         var_coarse = wk
      end if

      deallocate (wk, wk2)

      return
   end subroutine fine_to_coarseV

   ! Copy data from a fine-resolution array to a coarse one for probe
   subroutine fine_to_coarseP(ipencil, var_fine, var_coarse)

      implicit none

      real(mytype), dimension(:, :, :) :: var_fine
      real(mytype), dimension(:, :, :) :: var_coarse
      integer, intent(IN) :: ipencil

      real(mytype), allocatable, dimension(:, :, :) :: wk, wk2
      integer :: i, j, k

      if (ipencil == 1) then
         allocate (wk(xstP(1):xenP(1), xstP(2):xenP(2), xstP(3):xenP(3)))
         allocate (wk2(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
         wk2 = var_fine
         if (coarse_mesh_starts_from_1) then
            do k = xstP(3), xenP(3)
               do j = xstP(2), xenP(2)
                  do i = xstP(1), xenP(1)
                     wk(i, j, k) = wk2((i - 1) * iskipP + 1, (j - 1) * jskipP + 1, (k - 1) * kskipP + 1)
                  end do
               end do
            end do
         else
            do k = xstP(3), xenP(3)
               do j = xstP(2), xenP(2)
                  do i = xstP(1), xenP(1)
                     wk(i, j, k) = wk2(i * iskipP, j * jskipP, k * kskipP)
                  end do
               end do
            end do
         end if
         var_coarse = wk
      else if (ipencil == 2) then
         allocate (wk(ystP(1):yenP(1), ystP(2):yenP(2), ystP(3):yenP(3)))
         allocate (wk2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
         wk2 = var_fine
         if (coarse_mesh_starts_from_1) then
            do k = ystP(3), yenP(3)
               do j = ystP(2), yenP(2)
                  do i = ystP(1), yenP(1)
                     wk(i, j, k) = wk2((i - 1) * iskipP + 1, (j - 1) * jskipP + 1, (k - 1) * kskipP + 1)
                  end do
               end do
            end do
         else
            do k = ystP(3), yenP(3)
               do j = ystP(2), yenP(2)
                  do i = ystP(1), yenP(1)
                     wk(i, j, k) = wk2(i * iskipP, j * jskipP, k * kskipP)
                  end do
               end do
            end do
         end if
         var_coarse = wk
      else if (ipencil == 3) then
         allocate (wk(zstP(1):zenP(1), zstP(2):zenP(2), zstP(3):zenP(3)))
         allocate (wk2(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
         wk2 = var_fine
         if (coarse_mesh_starts_from_1) then
            do k = zstP(3), zenP(3)
               do j = zstP(2), zenP(2)
                  do i = zstP(1), zenP(1)
                     wk(i, j, k) = wk2((i - 1) * iskipP + 1, (j - 1) * jskipP + 1, (k - 1) * kskipP + 1)
                  end do
               end do
            end do
         else
            do k = zstP(3), zenP(3)
               do j = zstP(2), zenP(2)
                  do i = zstP(1), zenP(1)
                     wk(i, j, k) = wk2(i * iskipP, j * jskipP, k * kskipP)
                  end do
               end do
            end do
         end if
         var_coarse = wk
      end if

      deallocate (wk, wk2)

      return
   end subroutine fine_to_coarseP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Find sub-domain information held by current processor
   !   INPUT:
   !     nx, ny, nz - global data dimension
   !     pdim(3)    - number of processor grid in each dimension,
   !                  valid values: 1 - distibute locally;
   !                                2 - distribute across p_row;
   !                                3 - distribute across p_col
   !   OUTPUT:
   !     lstart(3)  - starting index
   !     lend(3)    - ending index
   !     lsize(3)   - size of the sub-block (redundant)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine partition(nx, ny, nz, pdim, lstart, lend, lsize)

      implicit none

      integer, intent(IN) :: nx, ny, nz
      integer, dimension(3), intent(IN) :: pdim
      integer, dimension(3), intent(OUT) :: lstart, lend, lsize

      integer, allocatable, dimension(:) :: st, en, sz
      integer :: i, gsize

      do i = 1, 3

         if (i == 1) then
            gsize = nx
         else if (i == 2) then
            gsize = ny
         else if (i == 3) then
            gsize = nz
         end if

         if (pdim(i) == 1) then        ! all local
            lstart(i) = 1
            lend(i) = gsize
            lsize(i) = gsize
         elseif (pdim(i) == 2) then    ! distribute across dims(1)
            allocate (st(0:dims(1) - 1))
            allocate (en(0:dims(1) - 1))
            allocate (sz(0:dims(1) - 1))
            call distribute(gsize, dims(1), st, en, sz)
            lstart(i) = st(coord(1))
            lend(i) = en(coord(1))
            lsize(i) = sz(coord(1))
            deallocate (st, en, sz)
         elseif (pdim(i) == 3) then    ! distribute across dims(2)
            allocate (st(0:dims(2) - 1))
            allocate (en(0:dims(2) - 1))
            allocate (sz(0:dims(2) - 1))
            call distribute(gsize, dims(2), st, en, sz)
            lstart(i) = st(coord(2))
            lend(i) = en(coord(2))
            lsize(i) = sz(coord(2))
            deallocate (st, en, sz)
         end if

      end do
      return

   end subroutine partition

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !   - distibutes grid points in one dimension
   !   - handles uneven distribution properly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine distribute(data1, proc, st, en, sz)

      implicit none
      ! data1 -- data size in any dimension to be partitioned
      ! proc  -- number of processors in that dimension
      ! st    -- array of starting index
      ! en    -- array of ending index
      ! sz    -- array of local size  (redundent)
      integer data1, proc, st(0:proc - 1), en(0:proc - 1), sz(0:proc - 1)
      integer i, size1, nl, nu

      size1 = data1 / proc
      nu = data1 - size1 * proc
      nl = proc - nu
      st(0) = 1
      sz(0) = size1
      en(0) = size1
      do i = 1, nl - 1
         st(i) = st(i - 1) + size1
         sz(i) = size1
         en(i) = en(i - 1) + size1
      end do
      size1 = size1 + 1
      do i = nl, proc - 1
         st(i) = en(i - 1) + 1
         sz(i) = size1
         en(i) = en(i - 1) + size1
      end do
      en(proc - 1) = data1
      sz(proc - 1) = data1 - st(proc - 1) + 1

      return
   end subroutine distribute

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  Define how each dimension is distributed across processors
   !    e.g. 17 meshes across 4 processor would be distibuted as (4,4,4,5)
   !    such global information is required locally at MPI_ALLTOALLV time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_dist(nx, ny, nz, decomp)

      implicit none

      integer, intent(IN) :: nx, ny, nz
      TYPE(DECOMP_INFO), intent(INOUT) :: decomp
      integer, allocatable, dimension(:) :: st, en

      allocate (st(0:dims(1) - 1))
      allocate (en(0:dims(1) - 1))
      call distribute(nx, dims(1), st, en, decomp%x1dist)
      call distribute(ny, dims(1), st, en, decomp%y1dist)
      deallocate (st, en)

      allocate (st(0:dims(2) - 1))
      allocate (en(0:dims(2) - 1))
      call distribute(ny, dims(2), st, en, decomp%y2dist)
      call distribute(nz, dims(2), st, en, decomp%z2dist)
      deallocate (st, en)

      return
   end subroutine get_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Prepare the send / receive buffers for MPI_ALLTOALLV communications
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine prepare_buffer(decomp)

      implicit none

      TYPE(DECOMP_INFO), intent(INOUT) :: decomp

      integer :: i

      !LG : AJOUTS "bidons" pour eviter un plantage en -O3 avec gcc9.3
      !       * la fonction sortait des valeurs 'aleatoires'
      !         et le calcul plantait dans MPI_ALLTOALLV
      !       * pas de plantage en O2

      if (nrank == 0) then
         open (newunit=i, file='temp.dat', form='unformatted')
         write (i) decomp%x1dist, decomp%y1dist, decomp%y2dist, decomp%z2dist, &
            decomp%xsz, decomp%ysz, decomp%zsz
         close (i, status='delete')
      end if

      ! MPI_ALLTOALLV buffer information

      do i = 0, dims(1) - 1
         decomp%x1cnts(i) = decomp%x1dist(i) * decomp%xsz(2) * decomp%xsz(3)
         decomp%y1cnts(i) = decomp%ysz(1) * decomp%y1dist(i) * decomp%ysz(3)
         if (i == 0) then
            decomp%x1disp(i) = 0  ! displacement is 0-based index
            decomp%y1disp(i) = 0
         else
            decomp%x1disp(i) = decomp%x1disp(i - 1) + decomp%x1cnts(i - 1)
            decomp%y1disp(i) = decomp%y1disp(i - 1) + decomp%y1cnts(i - 1)
         end if
      end do

      do i = 0, dims(2) - 1
         decomp%y2cnts(i) = decomp%ysz(1) * decomp%y2dist(i) * decomp%ysz(3)
         decomp%z2cnts(i) = decomp%zsz(1) * decomp%zsz(2) * decomp%z2dist(i)
         if (i == 0) then
            decomp%y2disp(i) = 0  ! displacement is 0-based index
            decomp%z2disp(i) = 0
         else
            decomp%y2disp(i) = decomp%y2disp(i - 1) + decomp%y2cnts(i - 1)
            decomp%z2disp(i) = decomp%z2disp(i - 1) + decomp%z2cnts(i - 1)
         end if
      end do

      ! MPI_ALLTOALL buffer information
#ifdef EVEN
      ! For evenly distributed data, following is an easier implementation.
      ! But it should be covered by the more general formulation below.
      !decomp%x1count = decomp%xsz(1)*decomp%xsz(2)*decomp%xsz(3)/dims(1)
      !decomp%y1count = decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3)/dims(1)
      !decomp%y2count = decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3)/dims(2)
      !decomp%z2count = decomp%zsz(1)*decomp%zsz(2)*decomp%zsz(3)/dims(2)

      ! For unevenly distributed data, pad smaller messages. Note the
      ! last blocks along pencils always get assigned more mesh points
      ! for X <=> Y transposes
      decomp%x1count = decomp%x1dist(dims(1) - 1) * &
                       decomp%y1dist(dims(1) - 1) * decomp%xsz(3)
      decomp%y1count = decomp%x1count
      ! for Y <=> Z transposes
      decomp%y2count = decomp%y2dist(dims(2) - 1) * &
                       decomp%z2dist(dims(2) - 1) * decomp%zsz(1)
      decomp%z2count = decomp%y2count
#endif

   end subroutine prepare_buffer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Transposition routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "transpose_x_to_y.f90"
#include "transpose_y_to_z.f90"
#include "transpose_z_to_y.f90"
#include "transpose_y_to_x.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Halo cell support
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "halo.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Utility routines to help allocate 3D arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "alloc.f90"

end module decomp_2d
