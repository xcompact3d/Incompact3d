!! SPDX-License-Identifier: BSD-3-Clause

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Halo cell support for neighbouring pencils to exchange data
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine update_halo_real_short(in, out, level, opt_global, opt_pencil)

     implicit none

     integer, intent(IN) :: level      ! levels of halo cells required
     real(mytype), dimension(:, :, :), intent(IN) :: in
     real(mytype), allocatable, dimension(:, :, :), intent(OUT) :: out
#if defined(_GPU)
     attributes(device) :: out
#endif
     logical, optional :: opt_global
     integer, intent(in), optional :: opt_pencil

     call update_halo(in, out, level, decomp_main, opt_global, opt_pencil)

  end subroutine update_halo_real_short

  subroutine update_halo_real(in, out, level, decomp, opt_global, opt_pencil)

     implicit none

     integer, intent(IN) :: level      ! levels of halo cells required
     real(mytype), dimension(:, :, :), intent(IN) :: in
     real(mytype), allocatable, dimension(:, :, :), intent(OUT) :: out
#if defined(_GPU)
     attributes(device) :: out
#endif
     TYPE(DECOMP_INFO), intent(in) :: decomp
     logical, optional :: opt_global
     integer, intent(in), optional :: opt_pencil

     logical :: global

     ! starting/ending index of array with halo cells
     integer :: xs, ys, zs, xe, ye, ze
     ! additional start end
     integer :: ist, ien, jst, jen, kst, ken

     integer :: i, j, k, s1, s2, s3, ierror
     integer :: data_type

     integer :: icount, ilength, ijump
     integer :: halo12, halo21, halo31, halo32
     integer, dimension(4) :: requests
     integer, dimension(MPI_STATUS_SIZE, 4) :: status
     integer :: tag_e, tag_w, tag_n, tag_s, tag_t, tag_b

     integer :: ipencil
     logical, save :: first_call_x = .true., first_call_y = .true., first_call_z = .true.

     data_type = real_type

#include "halo_common.f90"

     return
  end subroutine update_halo_real

  subroutine update_halo_complex_short(in, out, level, opt_global, opt_pencil)

     implicit none

     integer, intent(IN) :: level      ! levels of halo cells required
     complex(mytype), dimension(:, :, :), intent(IN) :: in
     complex(mytype), allocatable, dimension(:, :, :), intent(OUT) :: out
#if defined(_GPU)
     attributes(device) :: out
#endif
     logical, optional :: opt_global
     integer, intent(in), optional :: opt_pencil

     call update_halo(in, out, level, decomp_main, opt_global, opt_pencil)

  end subroutine update_halo_complex_short

  subroutine update_halo_complex(in, out, level, decomp, opt_global, opt_pencil)

     implicit none

     integer, intent(IN) :: level      ! levels of halo cells required
     complex(mytype), dimension(:, :, :), intent(IN) :: in
     complex(mytype), allocatable, dimension(:, :, :), intent(OUT) :: out
#if defined(_GPU)
     attributes(device) :: out
#endif
     TYPE(DECOMP_INFO), intent(in) :: decomp
     logical, optional :: opt_global
     integer, intent(in), optional :: opt_pencil

     logical :: global

     ! starting/ending index of array with halo cells
     integer :: xs, ys, zs, xe, ye, ze
     ! additional start end
     integer :: ist, ien, jst, jen, kst, ken

     integer :: i, j, k, s1, s2, s3, ierror
     integer :: data_type

     integer :: icount, ilength, ijump
     integer :: halo12, halo21, halo31, halo32
     integer, dimension(4) :: requests
     integer, dimension(MPI_STATUS_SIZE, 4) :: status
     integer :: tag_e, tag_w, tag_n, tag_s, tag_t, tag_b

     integer :: ipencil
     logical, save :: first_call_x = .true., first_call_y = .true., first_call_z = .true.

     data_type = complex_type

#include "halo_common.f90"

     return
  end subroutine update_halo_complex
