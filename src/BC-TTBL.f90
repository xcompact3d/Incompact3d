!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module ttbl
   use decomp_2d
   use variables
   use param
   use MPI

   implicit none

   real(mytype), save, allocatable, dimension(:, :) :: phisxz, phipphipsxz, upphipsxz, vpphipsxz, wpphipsxz, dphidysxz
   real(mytype), save, allocatable, dimension(:, :) :: usxz, vsxz, wsxz, upupsxz
   real(mytype), save, allocatable, dimension(:, :) :: vpvpsxz, wpwpsxz, upvpsxz, vpwpsxz, upwpsxz
   real(mytype), save, allocatable, dimension(:, :) :: tkesxz, epssxz
   real(mytype), save, allocatable, dimension(:) :: ttda, ttdb, ttdc
   real(mytype), save :: thetad

   public :: init_ttbl, boundary_conditions_ttbl, momentum_forcing_ttbl, scalar_forcing_ttbl, postprocess_ttbl, visu_ttbl_init, visu_ttbl

contains
   !############################################################################
   !############################################################################
   subroutine init_ttbl(ux1, uy1, uz1, phi1)
      implicit none
      real(mytype), intent(inout), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
      real(mytype), intent(inout), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1

      real(mytype), dimension(ysize(2)) :: um
      real(mytype) :: x, y, z, aform
      integer :: i, j, k, ii, code

      ! Set initial conditions
      ux1 = zero; uy1 = zero; uz1 = zero
      if (iin > 0) then
         if (iin <= 2) then
            call system_clock(count=code)
            if (iin == 2) code = 0
            call random_seed(size=ii)
            call random_seed(put=code + 63946 * (nrank + 1) * (/(i - 1, i=1, ii)/))
            call random_number(ux1)
            call random_number(uy1)
            call random_number(uz1)
            aform = sqrt(0.10922670d0 / 2.0d0)
            do k = 1, xsize(3)
               z = real((k + xstart(3) - 1 - 1), mytype) * dz
               do j = 1, xsize(2)
                  if (istret == 0) y = real(j + xstart(2) - 1 - 1, mytype) * dy
                  if (istret /= 0) y = yp(j + xstart(2) - 1)
                  um = (one - erf(aform * y))
                  do i = 1, xsize(1)
                     x = real(i - 1, mytype) * dx
                     ux1(i, j, k) = init_noise * um(j) * (two * ux1(i, j, k) - one) + erf(y * aform) + &
                                    4 * init_noise * (y * exp(-y) / 0.3678784468818499_mytype) * &
                                    (cos(z * pi / five) * cos(x * pi / five) + &
                                     cos((x + ((one + sqrt(five)) * half)) * pi / five) * cos((z + ((one + sqrt(five)) * half)) * pi / five))
                     uy1(i, j, k) = init_noise * um(j) * (two * uy1(i, j, k) - one)
                     uz1(i, j, k) = init_noise * um(j) * (two * uz1(i, j, k) - one)
                  end do
               end do
            end do
         else
            if (nrank == 0) write (6, *) 'Invalid inflow conditions for TTBL'
            call MPI_ABORT(MPI_COMM_WORLD, -1, code)
         end if
      end if
      phi1 = zero
      if (xstart(2) == 1) then
         ux1(:, 1, :) = zero
         uy1(:, 1, :) = zero
         uz1(:, 1, :) = zero
      end if

      ! Check for valid output parameters
      if (ioutput < ilist .or. mod(ioutput, ilist) /= 0) then
         if (nrank == 0) write (6, *) 'ioutput must be exactly divisible by ilist'
         call MPI_ABORT(MPI_COMM_WORLD, -1, code)
      end if

      ! Allocate output buffers for postprocessing
      if (nrank == 0) write (6, *) 'Output buffer size for postprocessing = ', ioutput / ilist
      allocate (usxz(ysize(2), ioutput / ilist))
      allocate (vsxz(ysize(2), ioutput / ilist))
      allocate (wsxz(ysize(2), ioutput / ilist))
      allocate (upupsxz(ysize(2), ioutput / ilist))
      allocate (vpvpsxz(ysize(2), ioutput / ilist))
      allocate (wpwpsxz(ysize(2), ioutput / ilist))
      allocate (upvpsxz(ysize(2), ioutput / ilist))
      allocate (vpwpsxz(ysize(2), ioutput / ilist))
      allocate (upwpsxz(ysize(2), ioutput / ilist))
      allocate (tkesxz(ysize(2), ioutput / ilist))
      allocate (epssxz(ysize(2), ioutput / ilist))
      allocate (phisxz(ysize(2), ioutput / ilist))
      allocate (phipphipsxz(ysize(2), ioutput / ilist))
      allocate (upphipsxz(ysize(2), ioutput / ilist))
      allocate (vpphipsxz(ysize(2), ioutput / ilist))
      allocate (wpphipsxz(ysize(2), ioutput / ilist))
      allocate (dphidysxz(ysize(2), ioutput / ilist))

      ! Allocate arrays for time integration
      allocate (ttda(ysize(2)), ttdb(ysize(2)), ttdc(ysize(2)))
   end subroutine init_ttbl
   !############################################################################
   !############################################################################
   subroutine boundary_conditions_ttbl(ux, uy, uz, phi)
      implicit none
      real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux, uy, uz
      real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi

      integer :: i, k

      ! Set BCs at initialisation
      if (itime == ifirst) then

         ! Bottom boundary
         if (ncly1 == 2) then
            do k = 1, xsize(3)
               do i = 1, xsize(1)
                  byx1(i, k) = zero
                  byy1(i, k) = zero
                  byz1(i, k) = zero
               end do
            end do
            ! Semi-implicit condition
            do k = 1, ysize(3)
               do i = 1, ysize(1)
                  byx1_2(i, k) = zero
                  byy1_2(i, k) = zero
                  byz1_2(i, k) = zero
               end do
            end do
         end if

         ! Top boundary
         if (nclyn == 2) then
            do k = 1, xsize(3)
               do i = 1, xsize(1)
                  byxn(i, k) = one
                  byyn(i, k) = zero
                  byzn(i, k) = zero
               end do
            end do
            ! Semi-implicit condition
            do k = 1, ysize(3)
               do i = 1, ysize(1)
                  byxn_2(i, k) = one
                  byyn_2(i, k) = zero
                  byzn_2(i, k) = zero
               end do
            end do
         end if
         phi = zero
      end if
   end subroutine boundary_conditions_ttbl
   !############################################################################
   !############################################################################
   subroutine momentum_forcing_ttbl(dux1, duy1, duz1, ux1, uy1, uz1, phi1)
      use var, only: ta1, tb1, tc1, ux2, uy2, uz2, ta2, tb2, tc2, di2
      use ibm_param, only: ubcx, ubcy, ubcz

      implicit none
      real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1

      real(mytype), dimension(ysize(2)) :: ux2m, ydudy2m, dudy2m, du2dy22m, duxuy2pm
      real(mytype) :: y, theta, delta, tau_wall, ufric
      integer :: i, j, k, code

      ! Get velocities and derivatives
      call transpose_x_to_y(ux1, ux2)
      call transpose_x_to_y(uy1, uy2)
      call transpose_x_to_y(uz1, uz2)
      call dery(ta2, ux2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, ubcx)
      call dery(tb2, uy2, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0, ubcy)
      call dery(tc2, uz2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, ubcz)
      call transpose_y_to_x(ta2, ta1)
      call transpose_y_to_x(tb2, tb1)
      call transpose_y_to_x(tc2, tc1)
      call horizontal_avrge(ux2, ux2m)

      ! Compute thetad
      thetad = comp_thetad(thetad, ux2, uy2, ux2m)

      ! Apply forcing
      do k = 1, xsize(3)
         do j = 1, xsize(2)
            do i = 1, xsize(1)
               if (istret == 0) y = (j + xstart(2) - 1 - 1) * dy
               if (istret /= 0) y = yp(j + xstart(2) - 1)
               dux1(i, j, k, 1) = dux1(i, j, k, 1) + thetad * y * ta1(i, j, k)
               duy1(i, j, k, 1) = duy1(i, j, k, 1) + thetad * y * tb1(i, j, k)
               if (iscalar == 1 .and. ri(1) /= 0) duy1(i, j, k, 1) = duy1(i, j, k, 1) + ri(1) * phi1(i, j, k, 1)
               duz1(i, j, k, 1) = duz1(i, j, k, 1) + thetad * y * tc1(i, j, k)
            end do
         end do
      end do

      ! Write out values
      if (mod(itime, ilist) == 0) then

         ! Calculate quantities
         theta = sum((ux2m * (one - ux2m)) * ypw)
         delta = sum((one - ux2m) * ypw)
         tau_wall = sum(ta2(:, 1, :))
         call MPI_ALLREDUCE(MPI_IN_PLACE, tau_wall, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
         tau_wall = tau_wall / real(nx * nz, mytype)
         ufric = sqrt(tau_wall * xnu)

         ! Write out
         if (nrank == 0) then
            write (6, "(' thetad = ',F14.12,'    theta = ',F14.12,'    delta = ',F14.12,'    H = ',F14.12)") thetad, theta, delta, delta / theta
            write (6, "(' tau_wall = ',F14.12,'    u_tau = ',F14.12,'    cf = ',F14.12)") tau_wall, ufric, two * ufric**2
            write (6, "(' Re_theta = ',F18.12,'    Re_delta = ',F18.12)") one / (theta * xnu), one / (xnu / delta)
            write (6, "(' dx+ = ',F16.12,'    dz+ = ',F16.12,'    dy+(min,max) = ',2F16.12)") dx * ufric * re, dz * ufric * re, dyp(1) * ufric * re, dyp(ny) * ufric * re
            open (unit=67, file='out/ttbl.dat', status='unknown', form='formatted', action='write', position='append')
            if (itime == ilist) write (67, "(11A20)") 't', 'thetad', 'theta', 'delta', 'tau_wall', 'u_tau', 'cf', 'dx+', 'dz+', 'dy+_min', 'dy+_max'
            write (67, "(11E20.12)") t, thetad, theta, delta, tau_wall, ufric, two * ufric**2, dx * ufric * re, dz * ufric * re, dyp(1) * ufric * re, dyp(ny) * ufric * re
            close (67)
         end if
      end if
   end subroutine momentum_forcing_ttbl
   !############################################################################
   !############################################################################
   subroutine scalar_forcing_ttbl(uy1, dphi1, phi1)
      use var, only: ta1, ta2, di2, phi2

      implicit none
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: uy1
      real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dphi1
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1

      real(mytype), dimension(ysize(2)) :: phi2m
      real(mytype) :: omegad, omega, y
      integer :: i, j, k

      ! Compute forcing for scalars
      if (iscalar == 1 .and. ri(1) /= 0) then
         call transpose_x_to_y(phi1(:, :, :, 1), phi2(:, :, :, 1))
         call horizontal_avrge(phi2(:, :, :, 1), phi2m)
         omega = sum(phi2m * ypw)
         omegad = sc(1) * thetad
         if (nrank == 0) write (*, "(' omegad = ',E15.7,'    omega= ',F14.12)") omegad, omega
         if (omegad > 0) then
            call deryS(ta2, phi2(:, :, :, 1), di2, sy, ffypS, fsypS, fwypS, ppy, ysize(1), ysize(2), ysize(3), 1, zero)
            call transpose_y_to_x(ta2, ta1)
            do k = 1, xsize(3)
               do j = 1, xsize(2)
                  do i = 1, xsize(1)
                     if (istret == 0) y = (j + xstart(2) - 1 - 1) * dy
                     if (istret /= 0) y = yp(j + xstart(2) - 1)
                     dphi1(i, j, k, 1) = dphi1(i, j, k, 1) + omegad * y * ta1(i, j, k)
                     if (ri(1) /= 0) dphi1(i, j, k, 1) = dphi1(i, j, k, 1) + uy1(i, j, k)
                  end do
               end do
            end do
         end if
      end if
   end subroutine scalar_forcing_ttbl
   !############################################################################
   !############################################################################
   subroutine postprocess_ttbl(ux1, uy1, uz1, pp3, phi1, ep1)
      use var, only: ux2, uy2, uz2
      use var, only: ta1, tb1, tc1, td1, te1, tf1, tg1, th1, ti1, di1
      use var, only: ta2, tb2, tc2, td2, te2, tf2, di2
      use var, only: ta3, tb3, tc3, td3, te3, tf3, di3, nzmsize, phi2
      use ibm_param, only: ubcx, ubcy, ubcz

      implicit none
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1, ep1
      real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1

      real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux1p, uy1p, uz1p, phi1p, dphi1pdx
      real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: pre2, ux2p, uy2p, uz2p, phi2p, pre2p
      real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: dphi2pdx, dphi2pdy, dphi2pdz
      real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: diss1, pre1
      real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: temp2, tke3d2
      real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: dphi3pdz, phi3p

      real(mytype), dimension(ysize(2)) :: sxz1
      real(mytype), dimension(ysize(2)) :: dup2dxsxz, dvp2dxsxz, dwp2dxsxz
      real(mytype), dimension(ysize(2)) :: dup2dysxz, dvp2dysxz, dwp2dysxz
      real(mytype), dimension(ysize(2)) :: dup2dzsxz, dvp2dzsxz, dwp2dzsxz
      real(mytype), dimension(ysize(2)) :: pre_sxz, dvpppsxz, vpppsxz, ppppsxz
      real(mytype), dimension(ysize(2)) :: kurt_x, kurt_y, kurt_z, kurt_phi, skew_x, skew_y, skew_z, skew_phi
      real(mytype), dimension(ysize(2)) :: dupupvpsxz, dvpwpwpsxz, dvpvpvpsxz
      real(mytype), dimension(ysize(2)) :: dusxz, dwsxz
      real(mytype), dimension(ysize(2)) :: dtkesxzdy, d2tkesxzdy
      real(mytype), dimension(ysize(2)) :: dvpphipphipsxz
      real(mytype), dimension(ysize(2)) :: d2phipphipsxz, dphipphipsxz
      real(mytype), dimension(ysize(2)) :: dphisxz
      real(mytype), dimension(ysize(2)) :: b_eps
      real(mytype), dimension(ysize(2)) :: dusxzdy, dwsxzdy

      real(8) :: tstart
      integer :: j, b

      ! Compute quantities and store in buffer
      if (itime > ceiling(real(initstat, mytype) / real(ioutput, mytype)) * ioutput .and. mod(itime, ilist) == 0) then
         tstart = MPI_WTIME()

         ! Get index into buffer
         b = mod(itime, ioutput) / ilist
         if (mod(itime, ioutput) == 0) b = ioutput / ilist

         ! Calculate derivatives
         call pressure_x(pp3, pre1)
         call derx(ta1, ux1, di1, sx, ffx, fsx, fwx, xsize(1), xsize(2), xsize(3), 0, ubcx)
         call derx(tb1, uy1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, ubcy)
         call derx(tc1, uz1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, ubcz)
         call transpose_x_to_y(ux1, td2)
         call transpose_x_to_y(uy1, te2)
         call transpose_x_to_y(uz1, tf2)
         call dery(ta2, td2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, ubcx)
         call dery(tb2, te2, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0, ubcy)
         call dery(tc2, tf2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, ubcz)
         call transpose_y_to_z(td2, td3)
         call transpose_y_to_z(te2, te3)
         call transpose_y_to_z(tf2, tf3)
         call derz(ta3, td3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, ubcx)
         call derz(tb3, te3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, ubcy)
         call derz(tc3, tf3, di3, sz, ffz, fsz, fwz, zsize(1), zsize(2), zsize(3), 0, ubcz)
         call transpose_z_to_y(ta3, td2)
         call transpose_z_to_y(tb3, te2)
         call transpose_z_to_y(tc3, tf2)
         call transpose_y_to_x(td2, tg1)
         call transpose_y_to_x(te2, th1)
         call transpose_y_to_x(tf2, ti1)
         call transpose_y_to_x(ta2, td1)
         call transpose_y_to_x(tb2, te1)
         call transpose_y_to_x(tc2, tf1)

         ! Mean and fluctuating velocities
         call transpose_x_to_y(ux1, ux2)
         call horizontal_avrge(ux2, usxz(:, b))
         call extract_fluctuat(ux2, usxz(:, b), ux2p)
         call transpose_y_to_x(ux2p, ux1p)
         call transpose_x_to_y(uy1, uy2)
         call horizontal_avrge(uy2, vsxz(:, b))
         call extract_fluctuat(uy2, vsxz(:, b), uy2p)
         call transpose_y_to_x(uy2p, uy1p)
         call transpose_x_to_y(uz1, uz2)
         call horizontal_avrge(uz2, wsxz(:, b))
         call extract_fluctuat(uz2, wsxz(:, b), uz2p)
         call transpose_y_to_x(uz2p, uz1p)

         ! RMS of fluctuating and TKE
         tke3d2 = half * (ux2p**2 + uy2p**2 + uz2p**2)
         call horizontal_avrge(ux2p**2, upupsxz(:, b))
         call horizontal_avrge(uy2p**2, vpvpsxz(:, b))
         call horizontal_avrge(uz2p**2, wpwpsxz(:, b))
         call horizontal_avrge(ux2p * uy2p, upvpsxz(:, b))
         call horizontal_avrge(uy2p * uz2p, vpwpsxz(:, b))
         call horizontal_avrge(ux2p * uz2p, upwpsxz(:, b))
         call horizontal_avrge(tke3d2, tkesxz(:, b))

         ! Dissipation rate
         call derx(ta1, ux1p, di1, sx, ffx, fsx, fwx, xsize(1), xsize(2), xsize(3), 0, ubcx)
         call derx(tb1, uy1p, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, ubcy)
         call derx(tc1, uz1p, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, ubcz)
         call dery(ta2, ux2p, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, ubcx)
         call dery(tb2, uy2p, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0, ubcy)
         call dery(tc2, uz2p, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, ubcz)
         call transpose_y_to_z(ux2p, td3)
         call transpose_y_to_z(uy2p, te3)
         call transpose_y_to_z(uz2p, tf3)
         call derz(ta3, td3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, ubcx)
         call derz(tb3, te3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, ubcy)
         call derz(tc3, tf3, di3, sz, ffz, fsz, fwz, zsize(1), zsize(2), zsize(3), 0, ubcz)
         call transpose_z_to_y(ta3, td2)
         call transpose_z_to_y(tb3, te2)
         call transpose_z_to_y(tc3, tf2)
         call transpose_y_to_x(td2, tg1)
         call transpose_y_to_x(te2, th1)
         call transpose_y_to_x(tf2, ti1)
         call transpose_y_to_x(ta2, td1)
         call transpose_y_to_x(tb2, te1)
         call transpose_y_to_x(tc2, tf1)
         call dissipation(ta1, tb1, tc1, td1, te1, tf1, tg1, th1, ti1, diss1)
         call transpose_x_to_y(diss1, temp2)
         call horizontal_avrge(temp2, epssxz(:, b))

         ! Scalars
         if (iscalar == 1) then
            call transpose_x_to_y(phi1(:, :, :, 1), phi2(:, :, :, 1))
            call horizontal_avrge(phi2(:, :, :, 1), phisxz(:, b))
            call extract_fluctuat(phi2(:, :, :, 1), phisxz(:, b), phi2p)
            call transpose_y_to_x(phi2p, phi1p)
            call horizontal_avrge(phi2p**2, phipphipsxz(:, b))
            call horizontal_avrge(ux2p * phi2p, upphipsxz(:, b))
            call horizontal_avrge(uy2p * phi2p, vpphipsxz(:, b))
            call horizontal_avrge(uz2p * phi2p, wpphipsxz(:, b))
            call deryS(dphidysxz(:, b), phisxz(:, b), di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, zero)
         end if
         if (nrank == 0) write (*, "(' Time computing basic statistics = ',F18.12,'(s)')") MPI_WTIME() - tstart

         ! Write times to disk
         if (nrank == 0) then
            open (unit=67, file='out/times.dat', status='unknown', form='formatted', action='write', position='append')
            if (itime == ceiling(real(initstat, mytype) / real(ioutput, mytype)) * ioutput + ilist) write (67, "(A12,A20)") 'itime', 't'
            write (67, "(I12,E20.12)") itime, t; close (67)
         end if

         ! Write to disk
         if (nrank == 0 .and. mod(itime, ioutput) == 0) then
            tstart = MPI_WTIME()
            call outp(usxz, 'out/ux.dat')
            call outp(vsxz, 'out/uy.dat')
            call outp(wsxz, 'out/uz.dat')
            call outp(upupsxz, 'out/uxp2.dat')
            call outp(vpvpsxz, 'out/uyp2.dat')
            call outp(wpwpsxz, 'out/uzp2.dat')
            call outp(upvpsxz, 'out/uxp_uyp.dat')
            call outp(vpwpsxz, 'out/uyp_uzp.dat')
            call outp(upwpsxz, 'out/uxp_uzp.dat')
            call outp(tkesxz, 'out/tke.dat')
            call outp(-epssxz, 'out/diss.dat')
            if (iscalar == 1) then
               call outp(phisxz, 'out/phi.dat')
               call outp(phipphipsxz, 'out/phip2.dat')
               call outp(upphipsxz, 'out/phip_uxp.dat')
               call outp(vpphipsxz, 'out/phip_uyp.dat')
               call outp(wpphipsxz, 'out/phip_uzp.dat')
               call outp(dphidysxz, 'out/brunt_vaisala.dat')
            end if
            if (nrank == 0) write (*, "(' Time writing basic statistics = ',F18.12,'(s)')") MPI_WTIME() - tstart
         end if

         ! Calculate other quantities less often
         if (mod(itime, ioutput) == 0) then
            tstart = MPI_WTIME()

            ! Transpose and average
            call transpose_x_to_y(ta1, temp2)
            call horizontal_avrge(temp2**2, dup2dxsxz)
            call transpose_x_to_y(tb1, temp2)
            call horizontal_avrge(temp2**2, dvp2dxsxz)
            call transpose_x_to_y(tc1, temp2)
            call horizontal_avrge(temp2**2, dwp2dxsxz)
            call transpose_x_to_y(td1, temp2)
            call horizontal_avrge(temp2**2, dup2dysxz)
            call transpose_x_to_y(te1, temp2)
            call horizontal_avrge(temp2**2, dvp2dysxz)
            call transpose_x_to_y(tf1, temp2)
            call horizontal_avrge(temp2**2, dwp2dysxz)
            call transpose_x_to_y(tg1, temp2)
            call horizontal_avrge(temp2**2, dup2dzsxz)
            call transpose_x_to_y(th1, temp2)
            call horizontal_avrge(temp2**2, dvp2dzsxz)
            call transpose_x_to_y(ti1, temp2)
            call horizontal_avrge(temp2**2, dwp2dzsxz)
            call transpose_x_to_y(pre1, pre2)
            call horizontal_avrge(pre2, pre_sxz)
            call extract_fluctuat(pre2, pre_sxz, pre2p)
            call horizontal_avrge(pre2p**2, ppppsxz)

            ! Pressure transport rate
            call horizontal_avrge(uy2p * pre2p, vpppsxz)
            call dery(dvpppsxz, vpppsxz, di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, ubcy)

            ! Turbulent transport
            call horizontal_avrge(ux2p * ux2p * uy2p, sxz1)
            call dery(dupupvpsxz, sxz1, di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, ubcx * (ubcx * ubcy))
            call horizontal_avrge(uy2p * uz2p * uz2p, sxz1)
            call dery(dvpwpwpsxz, sxz1, di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, ubcy * (ubcz * ubcz))
            call horizontal_avrge(uy2p * uy2p * uy2p, sxz1)
            call dery(dvpvpvpsxz, sxz1, di2, sy, ffy, fsy, fwy, ppy, 1, ysize(2), 1, 0, ubcy * (ubcy * ubcy))

            ! Shear production rate
            call dery(dusxz, usxz(:, b), di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, ubcx)
            call dery(dwsxz, wsxz(:, b), di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, ubcz)

            ! Viscous diffusion rate
            d2tkesxzdy = zero
            dtkesxzdy = zero
            call dery(dtkesxzdy, tkesxz(:, b), di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, ubcx)
            if (iimplicit <= 0) then
               call deryy(d2tkesxzdy, tkesxz(:, b), di2, sy, sfyp, ssyp, swyp, 1, ysize(2), 1, 1, ubcx)
               if (istret /= 0) then
                  do j = 1, ysize(2)
                     d2tkesxzdy(j) = d2tkesxzdy(j) * pp2y(j) - pp4y(j) * dtkesxzdy(j)
                  end do
               end if
            else
               do j = 1, ysize(2)
                  d2tkesxzdy(j) = -pp4y(j) * dtkesxzdy(j)
               end do
            end if

            ! Derivatives and higher order moments
            call dery(dusxzdy, usxz(:, b), di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, ubcx)
            call dery(dwsxzdy, wsxz(:, b), di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, ubcz)
            call horizontal_avrge(ux2p**three, skew_x)
            call horizontal_avrge(ux2p**four, kurt_x)
            call horizontal_avrge(uy2p**three, skew_y)
            call horizontal_avrge(uy2p**four, kurt_y)
            call horizontal_avrge(uz2p**three, skew_z)
            call horizontal_avrge(uz2p**four, kurt_z)

            ! Scalars
            if (iscalar == 1) then

               ! Buoyancy production and turbulent transport
               call deryS(dphisxz, phisxz(:, b), di2, sy, ffypS, fsypS, fwypS, ppy, 1, ysize(2), 1, 1, zero)
               call horizontal_avrge(uy2p * phi2p * phi2p, sxz1)
               call deryS(dvpphipphipsxz, sxz1, di2, sy, ffypS, fsypS, fwypS, ppy, 1, ysize(2), 1, 1, zero)

               ! Buoyancy diffusion rate
               d2phipphipsxz = zero
               call deryS(dphipphipsxz, phipphipsxz(:, b), di2, sy, ffypS, fsypS, fwypS, ppy, 1, ysize(2), 1, 1, zero)
               if (iimplicit <= 0) then
                  call deryyS(d2phipphipsxz, phipphipsxz(:, b), di2, sy, sfypS, ssypS, swypS, 1, ysize(2), 1, 1, zero)
                  if (istret /= 0) then
                     do j = 1, ysize(2)
                        d2phipphipsxz(j) = d2phipphipsxz(j) * pp2y(j) - dphipphipsxz(j) * pp4y(j)
                     end do
                  end if
               else
                  do j = 1, ysize(2)
                     d2phipphipsxz(j) = -dphipphipsxz(j) * pp4y(j)
                  end do
               end if

               ! Buoyancy dissipation rate
               call derxS(dphi1pdx, phi1p, di1, sx, ffxpS, fsxpS, fwxpS, xsize(1), xsize(2), xsize(3), 1, zero)
               call transpose_x_to_y(dphi1pdx, dphi2pdx)
               call deryS(dphi2pdy, phi2p, di2, sy, ffypS, fsypS, fwypS, ppy, ysize(1), ysize(2), ysize(3), 1, zero)
               call transpose_y_to_z(phi2p, phi3p)
               call derzS(dphi3pdz, phi3p, di3, sz, ffzpS, fszpS, fwzpS, zsize(1), zsize(2), zsize(3), 1, zero)
               call transpose_z_to_y(dphi3pdz, dphi2pdz)
               call horizontal_avrge(dphi2pdx**2 + dphi2pdy**2 + dphi2pdz**2, b_eps)

               ! Higher order moments
               call horizontal_avrge(phi2p**three, skew_phi)
               call horizontal_avrge(phi2p**four, kurt_phi)
            end if

            ! Write to disk
            if (nrank == 0 .and. mod(itime, ioutput) == 0) then
               call outpd(skew_x / (upupsxz(:, b)**1.5), 'out/skew_x.dat')
               call outpd(skew_y / (vpvpsxz(:, b)**1.5), 'out/skew_y.dat')
               call outpd(skew_z / (wpwpsxz(:, b)**1.5), 'out/skew_z.dat')
               call outpd(kurt_x / (upupsxz(:, b)**2), 'out/kurt_x.dat')
               call outpd(kurt_y / (vpvpsxz(:, b)**2), 'out/kurt_y.dat')
               call outpd(kurt_z / (wpwpsxz(:, b)**2), 'out/kurt_z.dat')
               call outpd(two * dup2dxsxz / dwp2dxsxz, 'out/K1.dat')
               call outpd(two * dup2dxsxz / dvp2dxsxz, 'out/K2.dat')
               call outpd(two * dup2dxsxz / dup2dzsxz, 'out/K3.dat')
               call outpd(two * dup2dxsxz / dup2dysxz, 'out/K4.dat')
               call outpd(-two * xnu * (dup2dxsxz + dup2dysxz + dup2dzsxz), 'out/eps_x.dat')
               call outpd(-two * xnu * (dvp2dxsxz + dvp2dysxz + dvp2dzsxz), 'out/eps_y.dat')
               call outpd(-two * xnu * (dwp2dxsxz + dwp2dysxz + dwp2dzsxz), 'out/eps_z.dat')
               call outpd(dusxzdy, 'out/dudy.dat')
               call outpd(dwsxzdy, 'out/dwdy.dat')
               call outpd(pre_sxz, 'out/pre.dat')
               call outpd(ppppsxz, 'out/pre2p.dat')
               call outpd(-dvpppsxz, 'out/velp_grad.dat')
               call outpd(-half * (dupupvpsxz + dvpwpwpsxz + dvpvpvpsxz), 'out/t_diff.dat')
               call outpd(-(upvpsxz(:, b) * dusxz + vpwpsxz(:, b) * dwsxz), 'out/prod.dat')
               call outpd(xnu * d2tkesxzdy, 'out/vis_diff.dat')
               call outpd(vsxz(:, b) * dtkesxzdy, 'out/meanf_transp.dat')
               if (iscalar == 1) then
                  call outpd(skew_phi / (wpwpsxz(:, b)**1.5), 'out/skew_phi.dat')
                  call outpd(kurt_phi / (phipphipsxz(:, b)**2), 'out/kurt_phi.dat')
                  call outpd(Ri(1) * vpphipsxz(:, b), 'out/b_flux.dat')
                  call outpd(-two * Ri(1) * vpphipsxz(:, b) * dphisxz, 'out/b_prod.dat')
                  call outpd(-Ri(1) * dvpphipphipsxz, 'out/b_ttransp.dat')
                  call outpd((Ri(1) * xnu / sc(1)) * d2phipphipsxz, 'out/b_diff.dat')
                  call outpd((-two * Ri(1) * xnu / Sc(1)) * b_eps, 'out/b_diss.dat')
               end if
            end if
         end if
      end if
   end subroutine postprocess_ttbl
   !############################################################################
   !############################################################################
   subroutine visu_ttbl_init(visu_initialised)
      implicit none
      logical, intent(out) :: visu_initialised
      visu_initialised = .true.
   end subroutine visu_ttbl_init
   !############################################################################
   !############################################################################
   subroutine visu_ttbl(ux1, uy1, uz1, pp3, phi1, ep1, num)
      use var, only: nzmsize
      implicit none
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
      real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ep1
      integer, intent(in) :: num
   end subroutine visu_ttbl
   !############################################################################
   !############################################################################
   subroutine horizontal_avrge(field, profile)
      implicit none
      real(mytype), intent(in), dimension(ysize(1), ysize(2), ysize(3)) :: field
      real(mytype), intent(out), dimension(ysize(2)) :: profile
      real(mytype), dimension(ysize(2)) :: sxz
      integer :: j, code
      do j = 1, ysize(2)
         sxz(j) = sum(field(:, j, :))
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE, sxz, ny, real_type, MPI_SUM, MPI_COMM_WORLD, code)
      profile = sxz / real(nx * nz, mytype)
   end subroutine horizontal_avrge
   !############################################################################
   !############################################################################
   subroutine extract_fluctuat(field, profile, fieldp)
      implicit none
      real(mytype), intent(in), dimension(ysize(1), ysize(2), ysize(3)) :: field
      real(mytype), intent(in), dimension(ysize(2)) :: profile
      real(mytype), intent(out), dimension(ysize(1), ysize(2), ysize(3)) :: fieldp
      integer :: j
      do j = 1, ysize(2)
         fieldp(:, j, :) = field(:, j, :) - profile(j)
      end do
   end subroutine extract_fluctuat
   !############################################################################
   !############################################################################
   function comp_thetad(thetad0, ux2, uy2, ux2m) result(thetad)
      use var, only: di2
      use ibm_param, only: ubcy
      use MPI

      implicit none
      real(mytype), intent(in) :: thetad0
      real(mytype), intent(in), dimension(ysize(1), ysize(2), ysize(3)) :: ux2, uy2
      real(mytype), intent(in), dimension(ysize(2)) :: ux2m

      real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: ux2p
      real(mytype), dimension(ysize(2)) :: uxuy2pm, ux2mm
      real(mytype), dimension(ysize(2)) :: ydudy2m, dudy2m, du2dy22m, duxuy2pm
      real(mytype) :: thetad
      integer :: j, code, it
      logical :: success

      integer, parameter :: freq = 1

      ! Calculate averaged derivatives
      call extract_fluctuat(ux2, ux2m, ux2p)
      call horizontal_avrge(ux2p * uy2, uxuy2pm)
      call dery(duxuy2pm, uxuy2pm, di2, sy, ffy, fsy, fwy, ppy, 1, ysize(2), 1, 0, ubcy)
      call dery(dudy2m, ux2m, di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, ubcy)
      ydudy2m = dudy2m * yp

      ! Implicit viscous term
      if (iimplicit <= 0) then
         call deryy(du2dy22m, ux2m, di2, sy, sfyp, ssyp, swyp, 1, ysize(2), 1, 1, ubcy)
         if (istret /= 0) then
            do j = 1, ysize(2)
               du2dy22m(j) = du2dy22m(j) * pp2y(j) - pp4y(j) * dudy2m(j)
            end do
         end if
      else
         do j = 1, ysize(2)
            du2dy22m(j) = -pp4y(j) * dudy2m(j)
         end do
      end if

      ! Iteratively find optimal thetad
      thetad = thetad0
      if (itime == ifirst) then
         thetad = 0.1092267 * xnu
         if (nrank == 0) write (6, *) 'Initial thetad = ', thetad
      else if (itime >= ifirst + 2 .and. mod(itime, freq) == 0) then
         if (nrank == 0) then
            call itp(comp_theta_res, thetad0, dt, thetad, it, success)
            if (success) then
               if (mod(itime, ilist) == 0) then
                  write (6, *) 'Theta dot computed in', it, ' with thetad = ', thetad, ' theta = ', comp_theta(thetad, ydudy2m, du2dy22m, duxuy2pm, ux2m, ux2mm)
               end if
            else
               write (6, *) 'Computing theta dot failed'
               call MPI_ABORT(MPI_COMM_WORLD, -1, code)
            end if
         end if
         call MPI_BCAST(thetad, 1, real_type, 0, MPI_COMM_WORLD, code)
      end if
      ttda(:) = thetad * ydudy2m(:) + xnu * du2dy22m(:) - duxuy2pm(:)
      ttdc(:) = ttdb(:); ttdb(:) = ttda(:)

      ! Internal wrapper to pass to ITP
   contains
      function comp_theta_res(thetad) result(theta_res)
         implicit none
         real(mytype), intent(in) :: thetad
         real(mytype) :: theta_res
         theta_res = comp_theta(thetad, ydudy2m, du2dy22m, duxuy2pm, ux2m, ux2mm) - one
      end function
   end function
   !############################################################################
   !############################################################################
   function comp_theta(thetad, ydudy2m, du2dy22m, duxuy2pm, ux2m, ux2mm) result(theta)
      implicit none
      real(mytype), intent(in) :: thetad
      real(mytype), intent(in), dimension(ysize(2)) :: ydudy2m, du2dy22m, duxuy2pm, ux2m
      real(mytype), intent(out), dimension(ysize(2)) :: ux2mm
      real(mytype) :: theta
      ttda(:) = thetad * ydudy2m(:) + xnu * du2dy22m(:) - duxuy2pm(:)
      ux2mm(:) = ux2m(:) + adt(1) * ttda(:) + bdt(1) * ttdb(:) + cdt(1) * ttdc(:)
      theta = sum((ux2mm * (one - ux2mm)) * ypw)
   end function
   !############################################################################
   !############################################################################
   subroutine itp(fun, x0, dx, x, it, success)
      implicit none
      real(mytype), intent(in) :: x0, dx
      real(mytype), intent(out) :: x
      integer, intent(out) :: it
      logical, intent(out) :: success

      integer :: n12, nmax
      real(mytype) :: ax, bx, fax, fbx, axm1, bxm1, denom, xf, x12, sigma, delta, xt, r, xitp, yitp
      logical :: too_small

      integer, parameter :: maxit = 100
      real(mytype), parameter :: k1 = 0.1_mytype, k2 = 0.98_mytype * (one + (one + sqrt(five)) / two)
#ifdef DOUBLE_PREC
      real(mytype), parameter :: tol = 1.0e-12_mytype
#else
      real(mytype), parameter :: tol = 1.0e-6_mytype
#endif

      ! Interface for function
      interface
         real(mytype) function fun(x)
            use param, only: mytype
            real(mytype), intent(in) :: x
         end function
      end interface

      ! Initalise outputs
      x = huge(one)
      success = .false.

      ! Find window to search
      do it = 1, maxit
         ax = x0 - two**(it - 1) * dx
         bx = x0 + two**(it - 1) * dx
         fax = fun(ax)
         fbx = fun(bx)
         if (fax * fbx < zero) exit
         if (it == maxit) return
      end do

      ! Swap values if necessary
      if (fax > fbx) then
         call swap(ax, bx)
         call swap(fax, fbx)
      end if

      ! Initialise additional parameters
      n12 = ceiling(log((bx - ax) / (two * tol)) / log(two))
      nmax = n12 + 1
      axm1 = huge(one)
      bxm1 = huge(one)

      ! Find root
      do it = 1, maxit
         denom = fbx - fax
         too_small = abs(denom) <= tiny(one)

         ! ITP method
         if (.not. too_small) then
            ! Interpolation
            xf = (fbx * ax - fax * bx) / denom
            ! Truncation
            x12 = (ax + bx) / two
            sigma = sign(one, x12 - xf)
            delta = k1 * (bx - ax)**k2
            if (delta <= abs(x12 - xf)) then
               xt = xf + sigma * delta
            else
               xt = x12
            end if
            ! Projection
            r = max(tol * two**(nmax - (it - 1)) - (bx - ax) / two, zero)
            if (abs(xt - x12) <= r) then
               xitp = xt
            else
               xitp = x12 - sigma * r
            end if
            ! Update
            yitp = fun(xitp)
            if (abs(yitp) < tol) then
               x = xitp
               success = .true.
               return
            end if
            if (yitp > zero) then
               bx = xitp
               fbx = yitp
            elseif (yitp < zero) then
               ax = xitp
               fax = yitp
            else
               ax = xitp
               bx = xitp
            end if
         end if

         ! Bisection method
         if (too_small .or. (ax == axm1 .and. bx == bxm1)) then
            xitp = (ax + bx) / two
            yitp = fun(xitp)
            if (abs(yitp) < tol) then
               x = xitp
               success = .true.
               return
            end if
            if (fax * yitp < zero) then
               bx = xitp
               fbx = yitp
            else
               ax = xitp
               fax = yitp
            end if
         end if
         axm1 = ax
         bxm1 = bx

         ! Check convergence
         if (abs(bx - ax) < two * tol) then
            x = (ax + bx) / two
            success = .true.
            return
         end if
      end do
   end subroutine itp
   !############################################################################
   !############################################################################
   subroutine swap(x, y)
      implicit none
      real(mytype), intent(inout) :: x, y
      real(mytype) :: temp
      temp = x
      x = y
      y = temp
   end subroutine swap
   !############################################################################
   !############################################################################
   subroutine outp(field, filename)
      implicit none
      real(mytype), intent(in), dimension(ysize(2), ioutput/ilist) :: field
      character(len=*), intent(in) :: filename
      integer :: i
      character(len=300) :: fileformat
      write (fileformat, '("(",I4,"E20.12)")') ny
      open (unit=67, file=trim(filename), status='unknown', form='formatted', action='write', position='append')
      do i = 1, ioutput / ilist
         write (67, fileformat) field(:, i)
      end do
      close (67)
   end subroutine outp
   !############################################################################
   !############################################################################
   subroutine outpd(field, filename)
      implicit none
      real(mytype), intent(in), dimension(ysize(2)) :: field
      character(len=*), intent(in) :: filename
      character(len=300) :: fileformat
      write (fileformat, '("(",I4,"E20.12)")') ny
      open (unit=67, file=trim(filename), status='unknown', form='formatted', action='write', position='append')
      write (67, fileformat) field
      close (67)
   end subroutine outpd
   !############################################################################
   !############################################################################
   subroutine pressure_x(pp3, pre1)
      use var, only: pp1, di1, nxmsize
      use var, only: pp2, ppi2, dip2, ph2, nymsize
      use var, only: ppi3, dip3, ph3, nzmsize
      use var, only: npress
      use tools, only: rescale_pressure

      implicit none
      real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
      real(mytype), intent(out), dimension(xsize(1), xsize(2), xsize(3)) :: pre1

      ! Interpolate and rescale pressure
      call interzpv(ppi3, pp3(:, :, :, 1), dip3, sz, cifip6z, cisip6z, ciwip6z, cifz6, cisz6, ciwz6, &
                    (ph3%zen(1) - ph3%zst(1) + 1), (ph3%zen(2) - ph3%zst(2) + 1), nzmsize, zsize(3), 1)
      call transpose_z_to_y(ppi3, pp2, ph3)
      call interypv(ppi2, pp2, dip2, sy, cifip6y, cisip6y, ciwip6y, cify6, cisy6, ciwy6, &
                    (ph3%yen(1) - ph3%yst(1) + 1), nymsize, ysize(2), ysize(3), 1)
      call transpose_y_to_x(ppi2, pp1, ph2)
      call interxpv(pre1, pp1, di1, sx, cifip6, cisip6, ciwip6, cifx6, cisx6, ciwx6, &
                    nxmsize, xsize(1), xsize(2), xsize(3), 1)
      call rescale_pressure(pre1)
   end subroutine pressure_x
   !############################################################################
   !############################################################################
   subroutine dissipation(ta1, tb1, tc1, td1, te1, tf1, tg1, th1, ti1, diss1)
      implicit none
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ta1, tb1, tc1, td1, te1, tf1, tg1, th1, ti1
      real(mytype), intent(out), dimension(xsize(1), xsize(2), xsize(3)) :: diss1

      real(mytype), dimension(3, 3, xsize(1), xsize(2), xsize(3)) :: A
      integer :: i, j, k, l, m

      ! Instantaneous dissipation rate
      diss1 = zero
      A(:, :, :, :, :) = zero
      A(1, 1, :, :, :) = ta1(:, :, :)
      A(2, 1, :, :, :) = tb1(:, :, :)
      A(3, 1, :, :, :) = tc1(:, :, :)
      A(1, 2, :, :, :) = td1(:, :, :)
      A(2, 2, :, :, :) = te1(:, :, :)
      A(3, 2, :, :, :) = tf1(:, :, :)
      A(1, 3, :, :, :) = tg1(:, :, :)
      A(2, 3, :, :, :) = th1(:, :, :)
      A(3, 3, :, :, :) = ti1(:, :, :)
      do k = 1, xsize(3)
         do j = 1, xsize(2)
            do i = 1, xsize(1)
               do m = 1, 3
                  do l = 1, 3
                     diss1(i, j, k) = diss1(i, j, k) + two * xnu * half * half * (A(l, m, i, j, k) + A(m, l, i, j, k))**two
                  end do
               end do
            end do
         end do
      end do
   end subroutine dissipation
end module ttbl
