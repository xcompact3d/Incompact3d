!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module cavity

   use decomp_2d
   use variables
   use param

   implicit none

   ! Temperature difference between both walls
   real(mytype), parameter :: deltaT = 1._mytype
   real(mytype), parameter :: temp_l = deltaT/2._mytype

   private ! All functions/subroutines private by default
   public :: init_cavity, boundary_conditions_cavity, postprocess_cavity

contains

   !############################################################################
   !!
   !! Init the cavity case
   !!
   !############################################################################
   subroutine init_cavity(ux1, uy1, uz1, ep1, phi1)

      USE decomp_2d_io
      USE MPI

      implicit none

      real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1, ep1
      real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1

      real(mytype) :: x
      integer :: i

      ! This does not apply in case of restart
      if (irestart == 0) then

         ! Velocity is zero
         ux1 = zero
         uy1 = zero
         uz1 = zero

         ! Linear temperature profile
         if (numscalar >= 1) then
            do i = 1, xsize(1)
               phi1(i, :, :, :) = temp_l - deltaT*(i - 1)/real(xsize(1) - 1, kind=mytype)
            end do
         end if

      end if

   end subroutine init_cavity

   !############################################################################
   !!
   !! Boundary conditions for the cavity case
   !!
   !############################################################################
   subroutine boundary_conditions_cavity(ux, uy, uz, phi)

      implicit none

      real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux, uy, uz
      real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi

      integer :: i, j, k

      ! Velocity
      IF (nclx1 == 2) THEN
      END IF
      IF (nclxn == 2) THEN
      END IF
      IF (ncly1 == 2) THEN
      END IF
      IF (nclyn == 2) THEN
      END IF
      IF (nclz1 == 2) THEN
      END IF
      IF (nclzn == 2) THEN
      END IF

      ! Scalar
      if (numscalar >= 1) then
         ! x=0
         if (nclxS1 == 2) then
            phi(1, :, :, :) = temp_l
         end if
         ! x=Lx
         if (nclxSn == 2) then
            phi(xsize(1), :, :, :) = temp_l - deltaT
         end if
         ! y=0
         if (nclyS1 == 2 .and. xstart(2) == 1) then
            do i = 1, xsize(1)
               phi(i, 1, :, :) = temp_l - deltaT*(i - 1)/real(xsize(1) - 1, kind=mytype)
            end do
         end if
         ! y=Ly
         if (nclySn == 2 .and. xend(2) == ny) then
            do i = 1, xsize(1)
               phi(i, xsize(2), :, :) = temp_l - deltaT*(i - 1)/real(xsize(1) - 1, kind=mytype)
            end do
         end if
      end if

      ! Clip
      if (numscalar >= 1) then
         do k = 1, xsize(3)
            do j = 1, xsize(2)
               do i = 1, xsize(1)
                  if (phi(i,j,k,1) > temp_l) phi(i,j,k,1) = temp_l
                  if (phi(i,j,k,1) < -temp_l) phi(i,j,k,1) = -temp_l
               enddo
            enddo
         enddo
      endif

   end subroutine boundary_conditions_cavity

   !############################################################################
   !!
   !! Post-processing for the cavity case
   !!
   !############################################################################
   subroutine postprocess_cavity(ux1, uy1, uz1, phi1)

      implicit none

      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1

   end subroutine postprocess_cavity

end module cavity
