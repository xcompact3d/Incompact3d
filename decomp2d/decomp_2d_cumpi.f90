!! SPDX-License-Identifier: BSD-3-Clause

! Module for the cuda aware MPI

module decomp_2d_cumpi

   use decomp_2d_constants
   use decomp_2d_mpi

   implicit none

   private        ! Make everything private unless declared public

   ! Device working arrays
   real(mytype), allocatable, dimension(:), device, public:: work1_r_d, work2_r_d
   complex(mytype), allocatable, dimension(:), device, public :: work1_c_d, work2_c_d

   public :: decomp_2d_cumpi_init, &
             decomp_2d_cumpi_fin

contains
   !
   ! init of the arrays
   !
   subroutine decomp_2d_cumpi_init(buf_size)

      implicit none

      integer, intent(in) :: buf_size
      integer :: status, errorcode

      if (allocated(work1_r_d)) deallocate (work1_r_d)
      if (allocated(work2_r_d)) deallocate (work2_r_d)
      if (allocated(work1_c_d)) deallocate (work1_c_d)
      if (allocated(work2_c_d)) deallocate (work2_c_d)
      allocate (work1_r_d(buf_size), STAT=status)
      if (status /= 0) then
         errorcode = 2
         call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                              'Out of memory when allocating 2DECOMP workspace')
      end if
      allocate (work1_c_d(buf_size), STAT=status)
      if (status /= 0) then
         errorcode = 2
         call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                              'Out of memory when allocating 2DECOMP workspace')
      end if
      allocate (work2_r_d(buf_size), STAT=status)
      if (status /= 0) then
         errorcode = 2
         call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                              'Out of memory when allocating 2DECOMP workspace')
      end if
      allocate (work2_c_d(buf_size), STAT=status)
      if (status /= 0) then
         errorcode = 2
         call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                              'Out of memory when allocating 2DECOMP workspace')
      end if

   end subroutine decomp_2d_cumpi_init
   !
   ! init of the arrays
   !
   subroutine decomp_2d_cumpi_fin

      implicit none

      deallocate (work1_r_d, work2_r_d, work1_c_d, work2_c_d)

   end subroutine decomp_2d_cumpi_fin

end module decomp_2d_cumpi

