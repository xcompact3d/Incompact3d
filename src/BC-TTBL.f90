!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module ttbl
   use decomp_2d
   use variables
   use param
   use MPI

   implicit none

   real(mytype), save, allocatable, dimension(:, :) :: usxz, vsxz, wsxz, upupsxz, vpvpsxz, wpwpsxz, upvpsxz
   real(mytype), save, allocatable, dimension(:, :) :: presxz, prep2sxz
   real(mytype), save, allocatable, dimension(:, :) :: dudysxz, dvdysxz, dwdysxz
   real(mytype), save, allocatable, dimension(:, :) :: tkesxz, meanconvsxz, prodsxz, disssxz, viscdiffsxz, turbconvsxz, prestransxz
   real(mytype), save, allocatable, dimension(:, :) :: meanconvsxz_uu, prodsxz_uu, disssxz_uu, viscdiffsxz_uu, turbconvsxz_uu, prestransxz_uu
   real(mytype), save, allocatable, dimension(:, :) :: meanconvsxz_vv, prodsxz_vv, disssxz_vv, viscdiffsxz_vv, turbconvsxz_vv, prestransxz_vv
   real(mytype), save, allocatable, dimension(:, :) :: meanconvsxz_ww, prodsxz_ww, disssxz_ww, viscdiffsxz_ww, turbconvsxz_ww, prestransxz_ww
   real(mytype), save, allocatable, dimension(:, :) :: meanconvsxz_uv, prodsxz_uv, disssxz_uv, viscdiffsxz_uv, turbconvsxz_uv, prestransxz_uv
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

      ! Allocate arrays and read data (if restart)
      if (itime == ifirst) then
         allocate (ttda(ysize(2)), ttdb(ysize(2)), ttdc(ysize(2)))
         if (irestart == 1) then
            open (unit=67, file='checkpoint_ttbl', status='unknown', form='unformatted', action='read')
            read (67) thetad, ttda, ttdb, ttdc
            close (67)
         end if
      end if

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

      ! Write data for restart
      if (nrank == 0 .and. mod(itime, icheckpoint) == 0) then
         open (unit=67, file='checkpoint_ttbl', status='unknown', form='unformatted', action='write')
         write (67) thetad, ttda, ttdb, ttdc
         close (67)
      end if

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

      real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux1p, uy1p, uz1p, pre1, diss1
      real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: ux2p, uy2p, uz2p, pre2, pre2p
      real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: ux3p, uy3p, uz3p
      real(mytype), dimension(ysize(2)) :: tempa2, tempb2, tempc2

      real(8) :: tstart
      integer :: j, b

      ! Allocate arrays
      if (itime == ifirst) then
         allocate (usxz(ysize(2), ioutput / ilist))
         allocate (vsxz(ysize(2), ioutput / ilist))
         allocate (wsxz(ysize(2), ioutput / ilist))
         allocate (upupsxz(ysize(2), ioutput / ilist))
         allocate (vpvpsxz(ysize(2), ioutput / ilist))
         allocate (wpwpsxz(ysize(2), ioutput / ilist))
         allocate (upvpsxz(ysize(2), ioutput / ilist))
         allocate (presxz(ysize(2), ioutput / ilist))
         allocate (prep2sxz(ysize(2), ioutput / ilist))
         allocate (dudysxz(ysize(2), ioutput / ilist))
         allocate (dvdysxz(ysize(2), ioutput / ilist))
         allocate (dwdysxz(ysize(2), ioutput / ilist))
         allocate (tkesxz(ysize(2), ioutput / ilist))
         allocate (meanconvsxz(ysize(2), ioutput / ilist))
         allocate (prodsxz(ysize(2), ioutput / ilist))
         allocate (disssxz(ysize(2), ioutput / ilist))
         allocate (viscdiffsxz(ysize(2), ioutput / ilist))
         allocate (turbconvsxz(ysize(2), ioutput / ilist))
         allocate (prestransxz(ysize(2), ioutput / ilist))
         allocate (meanconvsxz_uu(ysize(2), ioutput / ilist))
         allocate (prodsxz_uu(ysize(2), ioutput / ilist))
         allocate (disssxz_uu(ysize(2), ioutput / ilist))
         allocate (viscdiffsxz_uu(ysize(2), ioutput / ilist))
         allocate (turbconvsxz_uu(ysize(2), ioutput / ilist))
         allocate (prestransxz_uu(ysize(2), ioutput / ilist))
         allocate (meanconvsxz_vv(ysize(2), ioutput / ilist))
         allocate (prodsxz_vv(ysize(2), ioutput / ilist))
         allocate (disssxz_vv(ysize(2), ioutput / ilist))
         allocate (viscdiffsxz_vv(ysize(2), ioutput / ilist))
         allocate (turbconvsxz_vv(ysize(2), ioutput / ilist))
         allocate (prestransxz_vv(ysize(2), ioutput / ilist))
         allocate (meanconvsxz_ww(ysize(2), ioutput / ilist))
         allocate (prodsxz_ww(ysize(2), ioutput / ilist))
         allocate (disssxz_ww(ysize(2), ioutput / ilist))
         allocate (viscdiffsxz_ww(ysize(2), ioutput / ilist))
         allocate (turbconvsxz_ww(ysize(2), ioutput / ilist))
         allocate (prestransxz_ww(ysize(2), ioutput / ilist))
         allocate (meanconvsxz_uv(ysize(2), ioutput / ilist))
         allocate (prodsxz_uv(ysize(2), ioutput / ilist))
         allocate (disssxz_uv(ysize(2), ioutput / ilist))
         allocate (viscdiffsxz_uv(ysize(2), ioutput / ilist))
         allocate (turbconvsxz_uv(ysize(2), ioutput / ilist))
         allocate (prestransxz_uv(ysize(2), ioutput / ilist))
      end if

      ! Compute quantities and store in buffer
      if (itime > ceiling(real(initstat, mytype) / real(ioutput, mytype)) * ioutput .and. mod(itime, ilist) == 0) then
         tstart = MPI_WTIME()

         ! Get index into buffer
         b = mod(itime, ioutput) / ilist
         if (mod(itime, ioutput) == 0) b = ioutput / ilist

         ! Mean and fluctuating quantities
         call transpose_x_to_y(ux1, ux2)
         call horizontal_avrge(ux2, usxz(:, b))
         call extract_fluctuat(ux2, usxz(:, b), ux2p)
         call transpose_x_to_y(uy1, uy2)
         call horizontal_avrge(uy2, vsxz(:, b))
         call extract_fluctuat(uy2, vsxz(:, b), uy2p)
         call transpose_x_to_y(uz1, uz2)
         call horizontal_avrge(uz2, wsxz(:, b))
         call extract_fluctuat(uz2, wsxz(:, b), uz2p)
         call horizontal_avrge(ux2p**2, upupsxz(:, b))
         call horizontal_avrge(uy2p**2, vpvpsxz(:, b))
         call horizontal_avrge(uz2p**2, wpwpsxz(:, b))
         call horizontal_avrge(ux2p * uy2p, upvpsxz(:, b))
         call pressure_x(pp3, pre1)
         call transpose_x_to_y(pre1, pre2)
         call horizontal_avrge(pre2, presxz(:, b))
         call extract_fluctuat(pre2, presxz(:, b), pre2p)
         call horizontal_avrge(pre2p**2, prep2sxz(:, b))

         ! Velocity derivatives
         call dery(dudysxz(:, b), usxz(:, b), di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, ubcx)
         call dery(dvdysxz(:, b), vsxz(:, b), di2, sy, ffy, fsy, fwy, ppy, 1, ysize(2), 1, 0, ubcy)
         call dery(dwdysxz(:, b), wsxz(:, b), di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, ubcz)

         ! TKE
         tkesxz(:, b) = half * (upupsxz(:, b) + vpvpsxz(:, b) + wpwpsxz(:, b))

         ! Mean convection
         call dery(tempb2, tkesxz(:, b), di2, sy, ffy, fsy, fwy, ppy, 1, ysize(2), 1, 0, ubcy)
         meanconvsxz(:, b) = -(vsxz(:, b) * tempb2)

         ! Production
         call horizontal_avrge(uz2p * uy2p, tempc2)
         prodsxz(:, b) = -(upvpsxz(:, b) * dudysxz(:, b) + vpvpsxz(:, b) * dvdysxz(:, b) + tempc2 * dwdysxz(:, b))

         ! Dissipation
         call transpose_y_to_x(ux2p, ux1p)
         call transpose_y_to_x(uy2p, uy1p)
         call transpose_y_to_x(uz2p, uz1p)
         call transpose_y_to_z(ux2p, ux3p)
         call transpose_y_to_z(uy2p, uy3p)
         call transpose_y_to_z(uz2p, uz3p)
         call derx(ta1, ux1p, di1, sx, ffx, fsx, fwx, xsize(1), xsize(2), xsize(3), 0, ubcx)
         call derx(tb1, uy1p, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, ubcy)
         call derx(tc1, uz1p, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, ubcz)
         call dery(ta2, ux2p, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, ubcx)
         call dery(tb2, uy2p, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0, ubcy)
         call dery(tc2, uz2p, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, ubcz)
         call derz(ta3, ux3p, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, ubcx)
         call derz(tb3, uy3p, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, ubcy)
         call derz(tc3, uz3p, di3, sz, ffz, fsz, fwz, zsize(1), zsize(2), zsize(3), 0, ubcz)
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
         call transpose_x_to_y(diss1, ta2)
         call horizontal_avrge(-ta2, disssxz(:, b))

         ! Viscous diffusion
         call deryy(tempb2, tkesxz(:, b), di2, sy, sfyp, ssyp, swyp, 1, ysize(2), 1, 1, ubcx)
         viscdiffsxz(:, b) = xnu * tempb2

         ! Turbulent convection
         call horizontal_avrge(uy2p * ux2p * ux2p, tempa2)
         call horizontal_avrge(uy2p * uy2p * uy2p, tempb2)
         call horizontal_avrge(uy2p * uz2p * uz2p, tempc2)
         call dery(turbconvsxz(:, b), -half * (tempa2 + tempb2 + tempc2), di2, sy, ffy, fsy, fwy, ppy, 1, ysize(2), 1, 0, ubcy * (ubcx * ubcx + ubcy * ubcy + ubcz * ubcz))

         ! Pressure transport
         call horizontal_avrge(uy2p * pre2p, tempa2)
         call dery(prestransxz(:, b), -tempa2, di2, sy, ffy, fsy, fwy, ppy, 1, ysize(2), 1, 0, ubcy)

         ! Stress budget: uu
         ta2 = zero
         call stress_budget(usxz(:, b), ux2p, usxz(:, b), ux2p, ta2, ta2, uy2p, meanconvsxz_uu(:, b), prodsxz_uu(:, b), disssxz_uu(:, b), viscdiffsxz_uu(:, b), turbconvsxz_uu(:, b), prestransxz_uu(:, b))

         ! Stress budget: vv
         call dery(tb2, pre2p, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0, ubcy)
         call stress_budget(vsxz(:, b), uy2p, vsxz(:, b), uy2p, tb2, tb2, uy2p, meanconvsxz_vv(:, b), prodsxz_vv(:, b), disssxz_vv(:, b), viscdiffsxz_vv(:, b), turbconvsxz_vv(:, b), prestransxz_vv(:, b))

         ! Stress budget: ww
         call stress_budget(wsxz(:, b), uz2p, wsxz(:, b), uz2p, ta2, ta2, uy2p, meanconvsxz_ww(:, b), prodsxz_ww(:, b), disssxz_ww(:, b), viscdiffsxz_ww(:, b), turbconvsxz_ww(:, b), prestransxz_ww(:, b))

         ! Stress budget: uv
         call stress_budget(usxz(:, b), ux2p, vsxz(:, b), uy2p, ta2, tb2, uy2p, meanconvsxz_uv(:, b), prodsxz_uv(:, b), disssxz_uv(:, b), viscdiffsxz_uv(:, b), turbconvsxz_uv(:, b), prestransxz_uv(:, b))

         ! Write compute time
         if (nrank == 0) write (*, "(' Time computing statistics = ',F18.12,'(s)')") MPI_WTIME() - tstart

         ! Write times to disk
         if (nrank == 0) then
            open (unit=67, file='out/times.dat', status='unknown', form='formatted', action='write', position='append')
            if (itime == ceiling(real(initstat, mytype) / real(ioutput, mytype)) * ioutput + ilist) write (67, "(A12,A20)") 'itime', 't'
            write (67, "(I12,E20.12)") itime, t; close (67)
         end if

         ! Write to disk
         if (nrank == 0 .and. mod(itime, ioutput) == 0) then
            tstart = MPI_WTIME()
            call outp(usxz, 'out/u.dat')
            call outp(vsxz, 'out/v.dat')
            call outp(wsxz, 'out/w.dat')
            call outp(upupsxz, 'out/up2.dat')
            call outp(vpvpsxz, 'out/vp2.dat')
            call outp(wpwpsxz, 'out/wp2.dat')
            call outp(upvpsxz, 'out/upvp.dat')
            call outp(presxz, 'out/pre.dat')
            call outp(prep2sxz, 'out/prep2.dat')
            call outp(dudysxz, 'out/dudy.dat')
            call outp(dvdysxz, 'out/dvdy.dat')
            call outp(dwdysxz, 'out/dwdy.dat')
            call outp(tkesxz, 'out/tke.dat')
            call outp(meanconvsxz, 'out/mean_conv.dat')
            call outp(prodsxz, 'out/prod.dat')
            call outp(disssxz, 'out/diss.dat')
            call outp(viscdiffsxz, 'out/visc_diff.dat')
            call outp(turbconvsxz, 'out/turb_conv.dat')
            call outp(prestransxz, 'out/pres_tran.dat')
            call outp(meanconvsxz_uu, 'out/mean_conv_uu.dat')
            call outp(prodsxz_uu, 'out/prod_uu.dat')
            call outp(disssxz_uu, 'out/diss_uu.dat')
            call outp(viscdiffsxz_uu, 'out/visc_diff_uu.dat')
            call outp(turbconvsxz_uu, 'out/turb_conv_uu.dat')
            call outp(prestransxz_uu, 'out/pres_tran_uu.dat')
            call outp(meanconvsxz_vv, 'out/mean_conv_vv.dat')
            call outp(prodsxz_vv, 'out/prod_vv.dat')
            call outp(disssxz_vv, 'out/diss_vv.dat')
            call outp(viscdiffsxz_vv, 'out/visc_diff_vv.dat')
            call outp(turbconvsxz_vv, 'out/turb_conv_vv.dat')
            call outp(prestransxz_vv, 'out/pres_tran_vv.dat')
            call outp(meanconvsxz_ww, 'out/mean_conv_ww.dat')
            call outp(prodsxz_ww, 'out/prod_ww.dat')
            call outp(disssxz_ww, 'out/diss_ww.dat')
            call outp(viscdiffsxz_ww, 'out/visc_diff_ww.dat')
            call outp(turbconvsxz_ww, 'out/turb_conv_ww.dat')
            call outp(prestransxz_ww, 'out/pres_tran_ww.dat')
            call outp(meanconvsxz_uv, 'out/mean_conv_uv.dat')
            call outp(prodsxz_uv, 'out/prod_uv.dat')
            call outp(disssxz_uv, 'out/diss_uv.dat')
            call outp(viscdiffsxz_uv, 'out/visc_diff_uv.dat')
            call outp(turbconvsxz_uv, 'out/turb_conv_uv.dat')
            call outp(prestransxz_uv, 'out/pres_tran_uv.dat')
            if (nrank == 0) write (*, "(' Time writing statistics = ',F18.12,'(s)')") MPI_WTIME() - tstart
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
      if (itime == ifirst .and. irestart == 0) then
         thetad = 0.1092267 * xnu
         if (nrank == 0) write (6, *) 'Initial thetad = ', thetad
      else if ((itime >= ifirst + 2 .or. irestart == 1) .and. mod(itime, freq) == 0) then
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
      integer :: l, m

      ! Instantaneous dissipation rate
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
      diss1 = zero
      do m = 1, 3
         do l = 1, 3
            diss1 = diss1 + two * xnu * half * half * (A(l, m, :, :, :) + A(m, l, :, :, :))**two
         end do
      end do
   end subroutine dissipation
   !############################################################################
   !############################################################################
   subroutine stress_budget(ui2sxz, ui2p, uj2sxz, uj2p, dprepdi2, dprepdj2, uy2p, mean_conv, prod, diss, visc_diff, turb_conv, press_tran)
      use var, only: ux2, uy2, uz2
      use var, only: ta1, tb1, tc1, td1, te1, tf1, tg1, th1, ti1, di1
      use var, only: ta2, tb2, tc2, td2, te2, tf2, di2
      use var, only: ta3, tb3, tc3, td3, te3, tf3, di3, nzmsize, phi2
      use ibm_param, only: ubcx, ubcy, ubcz

      implicit none
      real(mytype), intent(in), dimension(ysize(1), ysize(2), ysize(3)) :: ui2p, uj2p, dprepdi2, dprepdj2, uy2p
      real(mytype), intent(in), dimension(ysize(2)) :: ui2sxz, uj2sxz
      real(mytype), intent(out), dimension(ysize(2)) :: mean_conv, prod, diss, visc_diff, turb_conv, press_tran

      real(mytype), dimension(ysize(2)) :: tempa2, tempb2, tempc2, tempd2

      ! Mean convection
      call horizontal_avrge(ui2p * uj2p, tempa2)
      call dery(tempb2, tempa2, di2, sy, ffy, fsy, fwy, ppy, 1, ysize(2), 1, 0, ubcy)
      call horizontal_avrge(uy2p, tempc2)
      mean_conv = -tempc2 * tempb2

      ! Production
      call horizontal_avrge(ui2p * uy2p, tempa2)
      call dery(tempb2, uj2sxz, di2, sy, ffy, fsy, fwy, ppy, 1, ysize(2), 1, 0, ubcy)
      call horizontal_avrge(uj2p * uy2p, tempc2)
      call dery(tempd2, ui2sxz, di2, sy, ffy, fsy, fwy, ppy, 1, ysize(2), 1, 0, ubcy)
      prod = -(tempa2 * tempb2 + tempc2 * tempd2)

      ! Dissipation
      call dery(ta2, ui2p, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0, ubcy)
      call dery(tb2, uj2p, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0, ubcy)
      call horizontal_avrge(ta2 * tb2, tempa2)
      diss = -two * xnu * tempa2

      ! Viscous diffusion
      call horizontal_avrge(ui2p * uj2p, tempa2)
      call deryy(tempb2, tempa2, di2, sy, sfyp, ssyp, swyp, 1, ysize(2), 1, 1, ubcx)
      visc_diff = xnu * tempb2

      ! Turbulent convection
      call horizontal_avrge(ui2p * uj2p * uy2p, tempa2)
      call dery(turb_conv, -tempa2, di2, sy, ffy, fsy, fwy, ppy, 1, ysize(2), 1, 0, ubcy)

      ! Pressure transport
      call horizontal_avrge(-(ui2p * dprepdj2 + uj2p * dprepdi2), press_tran)
   end subroutine stress_budget
end module ttbl
