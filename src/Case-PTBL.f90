!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module ptbl
  use decomp_2d
  use decomp_2d_mpi
  use variables
  use param
  use MPI

   implicit none

   real(mytype), save, allocatable, dimension(:, :) :: usxz, vsxz, wsxz, upupsxz, vpvpsxz, wpwpsxz, upvpsxz, tkesxz
   real(mytype), save, allocatable, dimension(:, :) :: u_ins, v_ins, w_ins, p_ins
   real(mytype), save, allocatable, dimension(:, :) :: presxz, prep2sxz
   real(mytype), save, allocatable, dimension(:, :) :: dudysxz, dvdysxz, dwdysxz
   real(mytype), save, allocatable, dimension(:, :) :: oxpsxz, oypsxz, ozpsxz
   real(mytype), save, allocatable, dimension(:, :) :: meanconvsxz_uu, turbconvsxz_uu, viscdiffsxz_uu, prodsxz_uu, prestransxz_uu, disssxz_uu
   real(mytype), save, allocatable, dimension(:, :) :: meanconvsxz_vv, turbconvsxz_vv, viscdiffsxz_vv, prodsxz_vv, prestransxz_vv, disssxz_vv
   real(mytype), save, allocatable, dimension(:, :) :: meanconvsxz_ww, turbconvsxz_ww, viscdiffsxz_ww, prodsxz_ww, prestransxz_ww, disssxz_ww
   real(mytype), save, allocatable, dimension(:, :) :: meanconvsxz_uv, turbconvsxz_uv, viscdiffsxz_uv, prodsxz_uv, prestransxz_uv, disssxz_uv
   real(mytype), save, allocatable, dimension(:, :) :: meanconvsxz_k, turbconvsxz_k, viscdiffsxz_k, prodsxz_k, prestransxz_k, disssxz_k
   real(mytype), save, allocatable, dimension(:) :: ttda, ttdb, ttdc
   real(mytype), save :: theta_a, theta_b, theta_c, ZTn
   real(mytype), save :: thetad

   public :: init_ptbl, boundary_conditions_ptbl, momentum_forcing_ptbl, scalar_forcing_ptbl, postprocess_ptbl, visu_ptbl_init, visu_ptbl

contains
   !############################################################################
   !############################################################################
   subroutine init_ptbl(ux1, uy1, uz1, phi1)
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
                  !==> Pasha
                  if ((jtheta_dot == 1) .and. (jthickness == 1)) y = y * (H_12)
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
            if (nrank == 0) write (*, *) 'Invalid inflow conditions for PTBL'
            call MPI_ABORT(MPI_COMM_WORLD, -1, code)
         end if
      end if
      phi1 = zero
      if (xstart(2) == 1) then
         ux1(:, 1, :) = zero
         uy1(:, 1, :) = zero
         uz1(:, 1, :) = zero
      end if
   end subroutine init_ptbl
   !############################################################################
   !############################################################################
   subroutine boundary_conditions_ptbl(ux, uy, uz, phi)
      implicit none
      real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux, uy, uz
      real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi
      real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: ux_2, uy_2, uz_2
      integer :: i, k
      real(mytype) :: xloc

      call transpose_x_to_y(ux,ux_2)
      call transpose_x_to_y(uy,uy_2)
      call transpose_x_to_y(uz,uz_2)

!      Set BCs at initialisation
!      if (itime == ifirst) then
!      end if
         ! Bottom boundary
         if (ncly1 == 2) then
            do k = 1, xsize(3)
               do i = 1, xsize(1)
                  xloc=real(i+xstart(1)-2,mytype)*dx
                  byx1(i, k) = zero
                  byy1(i, k) = zero
                  byz1(i, k) = zero
                  if (Blowing ==1 .and.  (xloc >= Xst_Blowing) .and. (xloc <= Xen_Blowing)) then 
                     byy1(i, k) = smoothening_function(xloc,Xst_Blowing,Xen_Blowing,Range_Smooth,A_Blowing) !A_Blowing
                  end if
               end do
            end do
            ! Semi-implicit condition
            do k = 1, ysize(3)
               do i = 1, ysize(1)
                  xloc=real(i+ystart(1)-2,mytype)*dx
                  byx1_2(i, k) = zero
                  byy1_2(i, k) = zero
                  byz1_2(i, k) = zero
                  if (Blowing ==1 .and.  (xloc >= Xst_Blowing) .and. (xloc <= Xen_Blowing)) then 
                     byy1_2(i, k) = smoothening_function(xloc,Xst_Blowing,Xen_Blowing,Range_Smooth,A_Blowing) !A_Blowing
                  end if
               end do
            end do
         end if

      ! Top boundary
         if (nclyn == 2) then
            do k = 1, xsize(3)
               do i = 1, xsize(1)
                  if (FreeStream == 0) then
                     byxn(i, k) = one
                     byyn(i, k) = zero
                     byzn(i, k) = zero
                  else if (FreeStream == 1) then
                     byxn(i, k) = ux(i, xsize(2) - 1, k)
                     byyn(i, k) = zero
                     byzn(i, k) = uz(i, xsize(2) - 1, k)
                  end if   
               end do
            end do
            ! Semi-implicit condition
            do k = 1, ysize(3)
               do i = 1, ysize(1)
                  if (FreeStream == 0) then
                     byxn_2(i, k) = one
                     byyn_2(i, k) = zero
                     byzn_2(i, k) = zero
                  else if (FreeStream == 1) then
                     byxn_2(i, k) = ux_2(i, ysize(2) - 1, k)
                     byyn_2(i, k) = zero
                     byzn_2(i, k) = uz_2(i, ysize(2) - 1, k)
                  end if
               end do
            end do
         end if
         phi = zero
      !end if
   end subroutine boundary_conditions_ptbl
   function smoothening_function(x, start_transition, end_transition, transition_width, constant_value) result(output)
      real(mytype), intent(in) :: x, start_transition, end_transition, transition_width, constant_value
      real :: output
      
      ! Compute the transition regions using the tanh function
      real :: transition_start, transition_end
      transition_start = 0.5 * (1.0 + tanh((x - start_transition) / transition_width))
      transition_end = 0.5 * (1.0 + tanh((end_transition - x) / transition_width))
      
      ! Combine the transition regions and constant value
      output = transition_start * transition_end * constant_value
      
  end function smoothening_function
   !############################################################################
   !############################################################################
   subroutine momentum_forcing_ptbl(dux1, duy1, duz1, ux1, uy1, uz1, phi1)
      use var, only: ta1, tb1, tc1, ux2, uy2, uz2, ta2, tb2, tc2, di2

      implicit none
      real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1

      real(mytype), dimension(ysize(2)) :: ux2m, ydudy2m, dudy2m, du2dy22m, duxuy2pm
      real(mytype) :: y, theta, delta, tau_wall, ufric
      integer :: i, j, k, code, newunit
      real(mytype), save :: P_DpDX

      ! Allocate arrays and read data (if restart)
      if (itime == ifirst) then
         if (APG .ne. 0) P_DpDX = APG_DpDX
         write (*, "(' 1 Adverse Pressure Gradient ',F14.12)") P_DpDX
         if (jtheta_dot ==0) allocate (ttda(ysize(2)), ttdb(ysize(2)), ttdc(ysize(2))) 
         if (irestart == 1) then
            open (unit=newunit, file='checkpoint_ptbl', status='unknown', form='unformatted', action='read')
            if (jtheta_dot ==0) then
               read (newunit) thetad, ttda, ttdb, ttdc
            !==> Andy method
            else if (((jtheta_dot ==1) .and. (ilesmod /=0) ) .or. ((jtheta_dot ==1) .and. (Method_FT ==1) ) )then
               read (newunit) thetad, theta_a, theta_b, theta_c , ZTn
            end if 
            close (newunit)
            if (APG .ne. 0) then 
               open (unit=newunit, file='checkpoint_APG', status='unknown', form='unformatted', action='read')
               read (newunit) P_DpDX
               write (*, "(' Read : Adverse Pressure Gradient ',F14.12)") P_DpDX
               close(newunit)
            end if
         end if
      end if

      ! Get velocities and derivatives
      call transpose_x_to_y(ux1, ux2)
      call transpose_x_to_y(uy1, uy2)
      call transpose_x_to_y(uz1, uz2)
      call dery(ta2, ux2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, zero)
      call dery(tb2, uy2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, zero)
      call dery(tc2, uz2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, zero)
      call transpose_y_to_x(ta2, ta1)
      call transpose_y_to_x(tb2, tb1)
      call transpose_y_to_x(tc2, tc1)
      call horizontal_avrge(ux2, ux2m)

      ! Compute thetad (0: Biau, 1: Andy)
      if (jtheta_dot == 0 ) then
         thetad = comp_thetad(thetad, ux2, uy2, ux2m)
      else if (jtheta_dot == 1 ) then
         thetad = comp_thetad_II(ux2, uy2, ux2m)
      end if

      ! Write data for restart
      if (nrank == 0 .and. mod(itime, icheckpoint) == 0) then
         open (unit=67, file='checkpoint_ptbl', status='unknown', form='unformatted', action='write')
         !==> Biau method
         if (jtheta_dot ==0) then
            write (67) thetad, ttda, ttdb, ttdc
         !==> Andy method
         else if (((jtheta_dot ==1) .and. (ilesmod /=0) ) .or. ((jtheta_dot ==1) .and. (Method_FT ==1) ) )then
         !else if (jtheta_dot ==1)then
               write (67) thetad, theta_a, theta_b, theta_c, ZTn
         end if 
         close (67)
         if (APG .ne. 0) then 
            open (unit=68, file='checkpoint_APG', status='unknown', form='unformatted', action='write')
            write (68) P_DpDX
            close(68)
         end if
      end if

 ! Write out values
      if (mod(itime, ilist) == 0) then

         ! Calculate quantities
         theta = sum((ux2m * (one - ux2m)) * ypw)
         delta = sum((one - ux2m) * ypw)
         tau_wall = sum(ta2(:, 1, :))
         call MPI_ALLREDUCE(MPI_IN_PLACE, tau_wall, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
         tau_wall = tau_wall / real(nx * nz, mytype)
         ufric = sqrt(tau_wall * xnu)
         
         ! Saved for adverse pressure gradient
         if (APG == 2) then 
            P_DpDX = APG_Beta * (tau_wall/delta)
            write (*, "(' 2 Adverse Pressure Gradient ',F14.12)") P_DpDX
            write (*, "(F14.12,F14.12,F14.12)")APG_Beta,tau_wall,delta
         end if

         !Check the flow rate (what is coming in the domain is getting out)
         if (mod(itime, ilist) == 0) call tbl_flrt_Check(ux1)

         ! Write out
         if (nrank == 0) then
            write (*, "(' thetad = ',F14.12,'    theta = ',F14.12,'    delta = ',F14.12,'    H = ',F14.12)") thetad, theta, delta, delta / theta
            write (*, "(' tau_wall = ',F16.12,'    u_tau = ',F16.12,'    cf = ',F16.12)") tau_wall, ufric, two * ufric**2
            write (*, "(' Re_theta = ',F16.8,'    Re_delta = ',F16.8)") one*theta*re, one*delta*re !(xnu=one/re)
            write (*, "(' dx+ = ',F12.8,'    dz+ = ',F12.8,'    dy+(min,max) = (',F12.8,',',F12.8,')')") dx * ufric * re, dz * ufric * re, dyp(1) * ufric * re, dyp(ny) * ufric * re
            write (*, "(' Y_size1 = ',I8,' Y_size2 = ',I8,' Y_size3 = ',I8)") ysize(1), ysize(2), ysize(3)
            if (APG == 2) then 
               write (*, "(' Adverse Pressure Gradient ',F14.12)") P_DpDX
            end if   
            open (unit=newunit, file='out/ptbl.dat', status='unknown', form='formatted', action='write', position='append')
            if (itime == ilist) write (newunit, "(12A20)") 't', 'thetad', 'theta', 'delta', 'tau_wall', 'u_tau', 'cf', 'dx+', 'dz+', 'dy+_min', 'dy+_max' 
                write (newunit, "(12E20.12)") t, thetad, theta, delta, tau_wall, ufric, two * ufric**2, dx * ufric * re, dz * ufric * re, dyp(1) * ufric * re, dyp(ny) * ufric * re
            close (newunit)
         end if

      end if

      ! Apply forcing
      do k = 1, xsize(3)
         do j = 1, xsize(2)
            do i = 1, xsize(1)
               if (istret == 0) y = (j + xstart(2) - 1 - 1) * dy
               if (istret /= 0) y = yp(j + xstart(2) - 1)
               if (APG == 0) then 
                  dux1(i, j, k, 1) = dux1(i, j, k, 1) + thetad * y * ta1(i, j, k)
               else if (APG == 1) then 
                  dux1(i, j, k, 1) = dux1(i, j, k, 1) + thetad * y * ta1(i, j, k) - APG_DpDX * (one - ux2m(j+xstart(2)-1) ) 
               else if (APG == 2) then 
                  dux1(i, j, k, 1) = dux1(i, j, k, 1) + thetad * y * ta1(i, j, k) - P_DpDX * (one - ux2m(j+xstart(2)-1) )                   
               end if
               duy1(i, j, k, 1) = duy1(i, j, k, 1) + thetad * y * tb1(i, j, k)
               if (iscalar == 1 .and. ri(1) /= 0) duy1(i, j, k, 1) = duy1(i, j, k, 1) + ri(1) * phi1(i, j, k, 1)
               duz1(i, j, k, 1) = duz1(i, j, k, 1) + thetad * y * tc1(i, j, k)
            end do
         end do
      end do

   end subroutine momentum_forcing_ptbl
   !############################################################################
   !############################################################################
   subroutine scalar_forcing_ptbl(uy1, dphi1, phi1)
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
   end subroutine scalar_forcing_ptbl
   !############################################################################
   !############################################################################
   subroutine postprocess_ptbl(ux1, uy1, uz1, pp3, phi1, ep1)
      use var, only: ux2, uy2, uz2
      use var, only: ux3, uy3, uz3
      use var, only: nzmsize
      use var, only: ta1, tb1, tc1, td1, te1, tf1, di1
      use var, only: ta2, tb2, tc2, td2, te2, tf2, tg2, th2, ti2, di2
      use var, only: ta3, tb3, tc3, td3, te3, tf3, di3

      implicit none
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1, ep1
      real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1

      real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: pre1
      real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: ux2p, uy2p, uz2p, pre2, pre2p

      real(8) :: tstart
      integer :: b

      ! Allocate arrays
      if (itime == ifirst) then
         !for instantanous profiles
         if (Pro_Spectra ==1) then
            allocate (u_ins(ysize(2), ioutput / ilist))
            allocate (v_ins(ysize(2), ioutput / ilist))
            allocate (w_ins(ysize(2), ioutput / ilist))
            allocate (p_ins(ysize(2), ioutput / ilist))
         end if
         allocate (usxz(ysize(2), ioutput / ilist))
         allocate (vsxz(ysize(2), ioutput / ilist))
         allocate (wsxz(ysize(2), ioutput / ilist))
         allocate (upupsxz(ysize(2), ioutput / ilist))
         allocate (vpvpsxz(ysize(2), ioutput / ilist))
         allocate (wpwpsxz(ysize(2), ioutput / ilist))
         allocate (upvpsxz(ysize(2), ioutput / ilist))
         allocate (tkesxz(ysize(2), ioutput / ilist))
         allocate (presxz(ysize(2), ioutput / ilist))
         allocate (prep2sxz(ysize(2), ioutput / ilist))
         allocate (dudysxz(ysize(2), ioutput / ilist))
         allocate (dvdysxz(ysize(2), ioutput / ilist))
         allocate (dwdysxz(ysize(2), ioutput / ilist))
         allocate (oxpsxz(ysize(2), ioutput / ilist))
         allocate (oypsxz(ysize(2), ioutput / ilist))
         allocate (ozpsxz(ysize(2), ioutput / ilist))
         allocate (meanconvsxz_uu(ysize(2), ioutput / ilist))
         allocate (turbconvsxz_uu(ysize(2), ioutput / ilist))
         allocate (viscdiffsxz_uu(ysize(2), ioutput / ilist))
         allocate (prodsxz_uu(ysize(2), ioutput / ilist))
         allocate (prestransxz_uu(ysize(2), ioutput / ilist))
         allocate (disssxz_uu(ysize(2), ioutput / ilist))
         allocate (meanconvsxz_vv(ysize(2), ioutput / ilist))
         allocate (turbconvsxz_vv(ysize(2), ioutput / ilist))
         allocate (viscdiffsxz_vv(ysize(2), ioutput / ilist))
         allocate (prodsxz_vv(ysize(2), ioutput / ilist))
         allocate (prestransxz_vv(ysize(2), ioutput / ilist))
         allocate (disssxz_vv(ysize(2), ioutput / ilist))
         allocate (meanconvsxz_ww(ysize(2), ioutput / ilist))
         allocate (turbconvsxz_ww(ysize(2), ioutput / ilist))
         allocate (viscdiffsxz_ww(ysize(2), ioutput / ilist))
         allocate (prodsxz_ww(ysize(2), ioutput / ilist))
         allocate (prestransxz_ww(ysize(2), ioutput / ilist))
         allocate (disssxz_ww(ysize(2), ioutput / ilist))
         allocate (meanconvsxz_uv(ysize(2), ioutput / ilist))
         allocate (turbconvsxz_uv(ysize(2), ioutput / ilist))
         allocate (viscdiffsxz_uv(ysize(2), ioutput / ilist))
         allocate (prodsxz_uv(ysize(2), ioutput / ilist))
         allocate (prestransxz_uv(ysize(2), ioutput / ilist))
         allocate (disssxz_uv(ysize(2), ioutput / ilist))
         allocate (meanconvsxz_k(ysize(2), ioutput / ilist))
         allocate (turbconvsxz_k(ysize(2), ioutput / ilist))
         allocate (viscdiffsxz_k(ysize(2), ioutput / ilist))
         allocate (prodsxz_k(ysize(2), ioutput / ilist))
         allocate (prestransxz_k(ysize(2), ioutput / ilist))
         allocate (disssxz_k(ysize(2), ioutput / ilist))
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
         if (Pro_Spectra ==1) call Gathering_Probe (ux2, u_ins(:, b))

         call transpose_x_to_y(uy1, uy2)
         call horizontal_avrge(uy2, vsxz(:, b))
         call extract_fluctuat(uy2, vsxz(:, b), uy2p)
         if (Pro_Spectra ==1) call Gathering_Probe (uy2, v_ins(:, b))

         call transpose_x_to_y(uz1, uz2)
         call horizontal_avrge(uz2, wsxz(:, b))
         call extract_fluctuat(uz2, wsxz(:, b), uz2p)
         if (Pro_Spectra ==1) call Gathering_Probe (uz2, w_ins(:, b))

         call horizontal_avrge(ux2p**2, upupsxz(:, b))
         call horizontal_avrge(uy2p**2, vpvpsxz(:, b))
         call horizontal_avrge(uz2p**2, wpwpsxz(:, b))
         call horizontal_avrge(ux2p * uy2p, upvpsxz(:, b))
         tkesxz(:, b) = half * (upupsxz(:, b) + vpvpsxz(:, b) + wpwpsxz(:, b))
         call pressure_x(pp3, pre1)
         call transpose_x_to_y(pre1, pre2)
         call horizontal_avrge(pre2, presxz(:, b))
         call extract_fluctuat(pre2, presxz(:, b), pre2p)
         call horizontal_avrge(pre2p**2, prep2sxz(:, b))
         if (Pro_Spectra ==1) call Gathering_Probe (pre2, p_ins(:, b))

         ! Wall-normal derivatives of mean velocity
         call dery(dudysxz(:, b), usxz(:, b), di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, zero)
         call dery(dvdysxz(:, b), vsxz(:, b), di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, zero)
         call dery(dwdysxz(:, b), wsxz(:, b), di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, zero)

         ! Cross derivatives of fluctuating velocity
         call transpose_y_to_x(uy2p, tb1)
         call transpose_y_to_x(uz2p, tc1)
         call transpose_y_to_z(ux2p, ta3)
         call transpose_y_to_z(uy2p, tb3)
         call derx(te1, tb1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, zero)
         call derx(tf1, tc1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, zero)
         call dery(td2, ux2p, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, zero)
         call dery(tf2, uz2p, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, zero)
         call derz(td3, ta3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, zero)
         call derz(te3, tb3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, zero)
         call transpose_x_to_y(te1, tb2)
         call transpose_x_to_y(tf1, tc2)
         call transpose_z_to_y(td3, tg2)
         call transpose_z_to_y(te3, th2)

         ! Vorticity
         ta2 = tf2 - th2
         te2 = tg2 - tc2
         ti2 = tb2 - td2
         call horizontal_avrge(ta2**2, oxpsxz(:, b))
         call horizontal_avrge(te2**2, oypsxz(:, b))
         call horizontal_avrge(ti2**2, ozpsxz(:, b))

         ! Reynolds stress components
         td2 = ux2p**2
         te2 = uy2p**2
         tf2 = uz2p**2
         tg2 = ux2p * uy2p

         ! Mean convection
         call mean_convection(td2, usxz(:, b), vsxz(:, b), wsxz(:, b), meanconvsxz_uu(:, b))
         call mean_convection(te2, usxz(:, b), vsxz(:, b), wsxz(:, b), meanconvsxz_vv(:, b))
         call mean_convection(tf2, usxz(:, b), vsxz(:, b), wsxz(:, b), meanconvsxz_ww(:, b))
         call mean_convection(tg2, usxz(:, b), vsxz(:, b), wsxz(:, b), meanconvsxz_uv(:, b))
         meanconvsxz_k(:, b) = half * (meanconvsxz_uu(:, b) + meanconvsxz_vv(:, b) + meanconvsxz_ww(:, b))

         ! Turbulent convection
         call turb_convection(td2, ux2p, uy2p, uz2p, turbconvsxz_uu(:, b))
         call turb_convection(te2, ux2p, uy2p, uz2p, turbconvsxz_vv(:, b))
         call turb_convection(tf2, ux2p, uy2p, uz2p, turbconvsxz_ww(:, b))
         call turb_convection(tg2, ux2p, uy2p, uz2p, turbconvsxz_uv(:, b))
         turbconvsxz_k(:, b) = half * (turbconvsxz_uu(:, b) + turbconvsxz_vv(:, b) + turbconvsxz_ww(:, b))

         ! Viscous diffusion
         call visc_diffusion(td2, viscdiffsxz_uu(:, b))
         call visc_diffusion(te2, viscdiffsxz_vv(:, b))
         call visc_diffusion(tf2, viscdiffsxz_ww(:, b))
         call visc_diffusion(tg2, viscdiffsxz_uv(:, b))
         viscdiffsxz_k(:, b) = half * (viscdiffsxz_uu(:, b) + viscdiffsxz_vv(:, b) + viscdiffsxz_ww(:, b))

         ! Production
         call production(ux2, ux2, ux2p, ux2p, ux2p, uy2p, uz2p, prodsxz_uu(:, b))
         call production(uy2, uy2, uy2p, uy2p, ux2p, uy2p, uz2p, prodsxz_vv(:, b))
         call production(uz2, uz2, uz2p, uz2p, ux2p, uy2p, uz2p, prodsxz_ww(:, b))
         call production(ux2, uy2, ux2p, uy2p, ux2p, uy2p, uz2p, prodsxz_uv(:, b))
         prodsxz_k(:, b) = half * (prodsxz_uu(:, b) + prodsxz_vv(:, b) + prodsxz_ww(:, b))

         ! Pressure transport
         call pres_transport(1, 1, ux2p, ux2p, pre2p, prestransxz_uu(:, b))
         call pres_transport(2, 2, uy2p, uy2p, pre2p, prestransxz_vv(:, b))
         call pres_transport(3, 3, uz2p, uz2p, pre2p, prestransxz_ww(:, b))
         call pres_transport(1, 2, ux2p, uy2p, pre2p, prestransxz_uv(:, b))
         prestransxz_k(:, b) = half * (prestransxz_uu(:, b) + prestransxz_vv(:, b) + prestransxz_ww(:, b))

         ! Dissipation
         call dissipation(ux2p, ux2p, disssxz_uu(:, b))
         call dissipation(uy2p, uy2p, disssxz_vv(:, b))
         call dissipation(uz2p, uz2p, disssxz_ww(:, b))
         call dissipation(ux2p, uy2p, disssxz_uv(:, b))
         disssxz_k(:, b) = half * (disssxz_uu(:, b) + disssxz_vv(:, b) + disssxz_ww(:, b))

         ! Write compute time
         if (nrank == 0) write (*, "(' Time computing statistics = ',F18.12,'(s)')") MPI_WTIME() - tstart

         ! Write times to disk
         if (nrank == 0) then
                  open (unit=67, file='out/times.dat', status='unknown', form='formatted', action='write', position='append')
                  if (itime == ceiling(real(initstat, mytype) / real(ioutput, mytype)) * ioutput + ilist) write (67, "(A12,A15,A20)") 'itime', 't' , 'Stat File #'
                  write (67, "(I12,E20.12,I9)") itime, t, (itime-initstat)/ilist; close (67)
         end if

         ! Write to disk
         if (nrank == 0 .and. mod(itime, ioutput) == 0) then
            tstart = MPI_WTIME()
            if (Pro_Spectra ==1) then
               call outp(u_ins, 'out/u_ins.dat')
               call outp(v_ins, 'out/v_ins.dat')
               call outp(w_ins, 'out/w_ins.dat')
               call outp(p_ins, 'out/p_ins.dat')
            end if
            call outp(usxz, 'out/u.dat')
            call outp(vsxz, 'out/v.dat')
            call outp(wsxz, 'out/w.dat')
            call outp(upupsxz, 'out/up2.dat')
            call outp(vpvpsxz, 'out/vp2.dat')
            call outp(wpwpsxz, 'out/wp2.dat')
            call outp(upvpsxz, 'out/upvp.dat')
            call outp(tkesxz, 'out/tke.dat')
            call outp(presxz, 'out/pre.dat')
            call outp(prep2sxz, 'out/prep2.dat')
            call outp(dudysxz, 'out/dudy.dat')
            call outp(dvdysxz, 'out/dvdy.dat')
            call outp(dwdysxz, 'out/dwdy.dat')
            call outp(oxpsxz, 'out/oxp2.dat')
            call outp(oypsxz, 'out/oyp2.dat')
            call outp(ozpsxz, 'out/ozp2.dat')
            call outp(meanconvsxz_uu, 'out/mean_conv_uu.dat')
            call outp(turbconvsxz_uu, 'out/turb_conv_uu.dat')
            call outp(viscdiffsxz_uu, 'out/visc_diff_uu.dat')
            call outp(prodsxz_uu, 'out/prod_uu.dat')
            call outp(prestransxz_uu, 'out/pres_tran_uu.dat')
            call outp(disssxz_uu, 'out/diss_uu.dat')
            call outp(meanconvsxz_vv, 'out/mean_conv_vv.dat')
            call outp(turbconvsxz_vv, 'out/turb_conv_vv.dat')
            call outp(viscdiffsxz_vv, 'out/visc_diff_vv.dat')
            call outp(prodsxz_vv, 'out/prod_vv.dat')
            call outp(prestransxz_vv, 'out/pres_tran_vv.dat')
            call outp(disssxz_vv, 'out/diss_vv.dat')
            call outp(meanconvsxz_ww, 'out/mean_conv_ww.dat')
            call outp(turbconvsxz_ww, 'out/turb_conv_ww.dat')
            call outp(viscdiffsxz_ww, 'out/visc_diff_ww.dat')
            call outp(prodsxz_ww, 'out/prod_ww.dat')
            call outp(prestransxz_ww, 'out/pres_tran_ww.dat')
            call outp(disssxz_ww, 'out/diss_ww.dat')
            call outp(meanconvsxz_uv, 'out/mean_conv_uv.dat')
            call outp(turbconvsxz_uv, 'out/turb_conv_uv.dat')
            call outp(viscdiffsxz_uv, 'out/visc_diff_uv.dat')
            call outp(prodsxz_uv, 'out/prod_uv.dat')
            call outp(prestransxz_uv, 'out/pres_tran_uv.dat')
            call outp(disssxz_uv, 'out/diss_uv.dat')
            call outp(meanconvsxz_k, 'out/mean_conv_k.dat')
            call outp(turbconvsxz_k, 'out/turb_conv_k.dat')
            call outp(viscdiffsxz_k, 'out/visc_diff_k.dat')
            call outp(prodsxz_k, 'out/prod_k.dat')
            call outp(prestransxz_k, 'out/pres_tran_k.dat')
            call outp(disssxz_k, 'out/diss_k.dat')
            if (nrank == 0) write (*, "(' Time writing statistics = ',F18.12,'(s)')") MPI_WTIME() - tstart
         end if
      end if
   end subroutine postprocess_ptbl
   !############################################################################
   !############################################################################
   subroutine visu_ptbl_init(visu_initialised)
      implicit none
      logical, intent(out) :: visu_initialised
      visu_initialised = .true.
   end subroutine visu_ptbl_init
   !############################################################################
   !############################################################################
   subroutine Gathering_Probe(field, profile)
      implicit none
      real(mytype), intent(in), dimension(ysize(1), ysize(2), ysize(3)) :: field
      real(mytype), intent(out), dimension(ysize(2)) :: profile
      real(mytype), dimension(ysize(2)) :: sxz
      real(mytype) :: x,z
      integer :: i, j, k, code

         do i = 1, ysize(1)
            x = (i + ystart(1) - 1 - 1) * dx
            if ((x > X_Pro_Spectra - 0.6*dx) .and. (x < X_Pro_Spectra + 0.6*dx)) then
               do k = 1, ysize(3)
                  z = (k + ystart(3) - 1 - 1) * dz
                  if ((z > Z_Pro_Spectra - 0.6*dz) .and. (z < Z_Pro_Spectra + 0.6*dz)) then
                     do j = 1, ysize(2)
                        !write (6, "(' Probe location is : ',F16.8,F16.8)") x,z
                        !write (6, "(' Probe location is : ',I8,I8)") i,k
                        sxz(j) = field(i, j, k)
                     end do
                  end if
               end do
            end if
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE, sxz, ny, real_type, MPI_SUM, MPI_COMM_WORLD, code)
         profile = sxz 

   end subroutine Gathering_Probe
   !############################################################################
   !############################################################################

   subroutine visu_ptbl(ux1, uy1, uz1, pp3, phi1, ep1, num)

     use var, only : ux2, uy2, uz2, ux3, uy3, uz3
     USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
     USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
     use var, ONLY : nxmsize, nymsize, nzmsize
     use visu, only : write_field
     use ibm_param, only : ubcx,ubcy,ubcz
     
     implicit none

     real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
     real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3
     real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
     real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
     integer, intent(in) :: num


   ! Perform communications if needed
     if (sync_vel_needed) then
      call transpose_x_to_y(ux1,ux2)
      call transpose_x_to_y(uy1,uy2)
      call transpose_x_to_y(uz1,uz2)
      call transpose_y_to_z(ux2,ux3)
      call transpose_y_to_z(uy2,uy3)
      call transpose_y_to_z(uz2,uz3)
      sync_vel_needed = .false.
    endif

    !x-derivatives
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz)
    !y-derivatives
    call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
    call dery (tb2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
    call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
    !!z-derivatives
    call derz (ta3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
    call derz (tb3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcy)
    call derz (tc3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcz)
    !!all back to x-pencils
    call transpose_z_to_y(ta3,td2)
    call transpose_z_to_y(tb3,te2)
    call transpose_z_to_y(tc3,tf2)
    call transpose_y_to_x(td2,tg1)
    call transpose_y_to_x(te2,th1)
    call transpose_y_to_x(tf2,ti1)
    call transpose_y_to_x(ta2,td1)
    call transpose_y_to_x(tb2,te1)
    call transpose_y_to_x(tc2,tf1)
    !du/dx=ta1 du/dy=td1 and du/dz=tg1
    !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
    !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1
    !VORTICITY FIELD
    di1 = zero
    di1(:,:,:)=sqrt(  (tf1(:,:,:)-th1(:,:,:))**2 &
                    + (tg1(:,:,:)-tc1(:,:,:))**2 &
                    + (tb1(:,:,:)-td1(:,:,:))**2)
    call write_field(di1, ".", "vort", num, flush = .true.) ! Reusing temporary array, force flush

    !Q=-0.5*(ta1**2+te1**2+ti1**2)-td1*tb1-tg1*tc1-th1*tf1
    di1 = zero
    di1(:,:,: ) = - half*(ta1(:,:,:)**2+te1(:,:,:)**2+ti1(:,:,:)**2) &
                  - td1(:,:,:)*tb1(:,:,:) &
                  - tg1(:,:,:)*tc1(:,:,:) &
                  - th1(:,:,:)*tf1(:,:,:)
    call write_field(di1, ".", "critq", num, flush = .true.) ! Reusing temporary array, force flush

      
   end subroutine visu_ptbl
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
      call dery(duxuy2pm, uxuy2pm, di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, zero)
      call dery(dudy2m, ux2m, di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, zero)
      ydudy2m = dudy2m * yp

      ! Viscous term
      iimplicit = -iimplicit
      call deryy(du2dy22m, ux2m, di2, sy, sfyp, ssyp, swyp, 1, ysize(2), 1, 1, zero)
      iimplicit = -iimplicit
      if (istret /= 0) then
         du2dy22m = du2dy22m * pp2y - pp4y * dudy2m
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
      if (APG == 0) then 
         ttda(:) = thetad * ydudy2m(:) + xnu * du2dy22m(:) - duxuy2pm(:) 
      elseif (APG == 1) then 
         ttda(:) = thetad * ydudy2m(:) + xnu * du2dy22m(:) - duxuy2pm(:) - APG_DpDX*(one-ux2m(:))
      end if
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
      if (APG == 0) then 
         ttda(:) = thetad * ydudy2m(:) + xnu * du2dy22m(:) - duxuy2pm(:)
      elseif (APG == 1) then 
         ttda(:) = thetad * ydudy2m(:) + xnu * du2dy22m(:) - duxuy2pm(:) - APG_DpDX*(one-ux2m(:))
      end if 
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
   subroutine mean_convection(u1pu2p, um, vm, wm, mean_conv)
      use var, only: ta1, tb1, tc1, di1
      use var, only: ta2, tb2, tc2, di2
      use var, only: ta3, tb3, tc3, di3

      implicit none
      real(mytype), intent(in), dimension(ysize(1), ysize(2), ysize(3)) :: u1pu2p
      real(mytype), intent(in), dimension(ysize(2)) :: um, vm, wm
      real(mytype), intent(out), dimension(ysize(2)) :: mean_conv

      real(mytype), dimension(ysize(2)) :: tempa2, tempb2, tempc2

      call transpose_y_to_x(u1pu2p, ta1)
      call derx(tb1, ta1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, zero)
      call transpose_x_to_y(tb1, ta2)
      call horizontal_avrge(ta2, tempa2)
      call horizontal_avrge(u1pu2p, tempc2)
      call dery(tempb2, tempc2, di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, zero)
      call transpose_y_to_z(u1pu2p, ta3)
      call derz(tb3, ta3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, zero)
      call transpose_z_to_y(tb3, tc2)
      call horizontal_avrge(tc2, tempc2)
      mean_conv = -(um * tempa2 + vm * tempb2 + wm * tempc2)
   end subroutine mean_convection
   !############################################################################
   !############################################################################
   subroutine turb_convection(u1pu2p, up, vp, wp, turb_conv)
      use var, only: ta1, tb1, tc1, di1
      use var, only: ta2, tb2, tc2, di2
      use var, only: ta3, tb3, tc3, di3

      implicit none
      real(mytype), intent(in), dimension(ysize(1), ysize(2), ysize(3)) :: u1pu2p, up, vp, wp
      real(mytype), intent(out), dimension(ysize(2)) :: turb_conv

      real(mytype), dimension(ysize(2)) :: tempa2, tempb2, tempc2

      call transpose_y_to_x(u1pu2p * up, ta1)
      call derx(tb1, ta1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, zero)
      call transpose_x_to_y(tb1, ta2)
      call horizontal_avrge(ta2, tempa2)
      call horizontal_avrge(u1pu2p * vp, tempc2)
      call dery(tempb2, tempc2, di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, zero)
      call transpose_y_to_z(u1pu2p * wp, ta3)
      call derz(tb3, ta3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, zero)
      call transpose_z_to_y(tb3, tc2)
      call horizontal_avrge(tc2, tempc2)
      turb_conv = -(tempa2 + tempb2 + tempc2)
   end subroutine turb_convection
   !############################################################################
   !############################################################################
   subroutine visc_diffusion(u1pu2p, visc_diff)
      use var, only: ta1, tb1, tc1, di1
      use var, only: ta2, tb2, tc2, di2
      use var, only: ta3, tb3, tc3, di3

      implicit none
      real(mytype), intent(in), dimension(ysize(1), ysize(2), ysize(3)) :: u1pu2p
      real(mytype), intent(out), dimension(ysize(2)) :: visc_diff

      real(mytype), dimension(ysize(2)) :: tempa2, tempb2, tempc2, tempd2

      call transpose_y_to_x(u1pu2p, ta1)
      call derxx(tb1, ta1, di1, sx, sfxp, ssxp, swxp, xsize(1), xsize(2), xsize(3), 1, zero)
      call transpose_x_to_y(tb1, ta2)
      call horizontal_avrge(ta2, tempa2)
      call horizontal_avrge(u1pu2p, tempc2)
      iimplicit = -iimplicit
      call deryy(tempb2, tempc2, di2, sy, sfyp, ssyp, swyp, 1, ysize(2), 1, 1, zero)
      iimplicit = -iimplicit
      if (istret /= 0) then
         call dery(tempd2, tempc2, di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, zero)
         tempb2 = tempb2 * pp2y - pp4y * tempd2
      end if
      call transpose_y_to_z(u1pu2p, ta3)
      call derzz(tb3, ta3, di3, sz, sfzp, sszp, swzp, zsize(1), zsize(2), zsize(3), 1, zero)
      call transpose_z_to_y(tb3, tc2)
      call horizontal_avrge(tc2, tempc2)
      visc_diff = xnu * (tempa2 + tempb2 + tempc2)
   end subroutine visc_diffusion
   !############################################################################
   !############################################################################
   subroutine production(u1, u2, u1p, u2p, up, vp, wp, prod)
      use var, only: ta1, tb1, tc1, di1
      use var, only: ta2, tb2, tc2, di2
      use var, only: ta3, tb3, tc3, di3

      implicit none
      real(mytype), intent(in), dimension(ysize(1), ysize(2), ysize(3)) :: u1, u2, u1p, u2p, up, vp, wp
      real(mytype), intent(out), dimension(ysize(2)) :: prod

      real(mytype), dimension(ysize(2)) :: tempa2, tempb2, tempc2, dudx, dudy, dudz

      call transpose_y_to_x(u2, ta1)
      call derx(tb1, ta1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, zero)
      call transpose_x_to_y(tb1, ta2)
      call horizontal_avrge(ta2, dudx)
      call horizontal_avrge(u2, tempc2)
      call dery(dudy, tempc2, di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, zero)
      call transpose_y_to_z(u2, ta3)
      call derz(tb3, ta3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, zero)
      call transpose_z_to_y(tb3, tc2)
      call horizontal_avrge(tc2, dudz)
      call horizontal_avrge(u1p * up, tempa2)
      call horizontal_avrge(u1p * vp, tempb2)
      call horizontal_avrge(u1p * wp, tempc2)
      prod = -(tempa2 * dudx + tempb2 * dudy + tempc2 * dudz)
      call transpose_y_to_x(u1, ta1)
      call derx(tb1, ta1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, zero)
      call transpose_x_to_y(tb1, ta2)
      call horizontal_avrge(ta2, dudx)
      call horizontal_avrge(u1, tempc2)
      call dery(dudy, tempc2, di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, zero)
      call transpose_y_to_z(u1, ta3)
      call derz(tb3, ta3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, zero)
      call transpose_z_to_y(tb3, tc2)
      call horizontal_avrge(tc2, dudz)
      call horizontal_avrge(u2p * up, tempa2)
      call horizontal_avrge(u2p * vp, tempb2)
      call horizontal_avrge(u2p * wp, tempc2)
      prod = prod - (tempa2 * dudx + tempb2 * dudy + tempc2 * dudz)
   end subroutine production
   !############################################################################
   !############################################################################
   subroutine pres_transport(one, two, u1p, u2p, prep, pres_tran)
      use var, only: ta1, tb1, tc1, di1
      use var, only: ta2, tb2, tc2, di2
      use var, only: ta3, tb3, tc3, di3

      implicit none
      integer, intent(in) :: one, two
      real(mytype), intent(in), dimension(ysize(1), ysize(2), ysize(3)) :: u1p, u2p, prep
      real(mytype), intent(out), dimension(ysize(2)) :: pres_tran

      integer :: code
      real(mytype), dimension(ysize(2)) :: tempa2, tempb2, tempc2

      ta2 = u1p * prep
      if (two == 1) then
         call transpose_y_to_x(ta2, ta1)
         call derx(tb1, ta1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, zero)
         call transpose_x_to_y(tb1, ta2)
         call horizontal_avrge(ta2, tempa2)
      else if (two == 2) then
         call horizontal_avrge(ta2, tempc2)
         call dery(tempa2, tempc2, di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, zero)
      else if (two == 3) then
         call transpose_y_to_z(ta2, ta3)
         call derz(tb3, ta3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, zero)
         call transpose_z_to_y(tb3, tc2)
         call horizontal_avrge(tc2, tempa2)
      else
         if (nrank == 0) write (6, *) 'Invalid component index for pressure transport'
         call MPI_ABORT(MPI_COMM_WORLD, -1, code)
      end if
      ta2 = u2p * prep
      if (one == 1) then
         call transpose_y_to_x(ta2, ta1)
         call derx(tb1, ta1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, zero)
         call transpose_x_to_y(tb1, ta2)
         call horizontal_avrge(ta2, tempb2)
      else if (one == 2) then
         call horizontal_avrge(ta2, tempc2)
         call dery(tempb2, tempc2, di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, zero)
      else if (one == 3) then
         call transpose_y_to_z(ta2, ta3)
         call derz(tb3, ta3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, zero)
         call transpose_z_to_y(tb3, tc2)
         call horizontal_avrge(tc2, tempb2)
      else
         if (nrank == 0) write (6, *) 'Invalid component index for pressure transport'
         call MPI_ABORT(MPI_COMM_WORLD, -1, code)
      end if
      pres_tran = -(tempa2 + tempb2)
   end subroutine pres_transport
   !############################################################################
   !############################################################################1
   subroutine dissipation(u1p, u2p, diss)
      use var, only: ta1, tb1, tc1, di1
      use var, only: ta2, tb2, tc2, di2
      use var, only: ta3, tb3, tc3, di3

      implicit none
      real(mytype), intent(in), dimension(ysize(1), ysize(2), ysize(3)) :: u1p, u2p
      real(mytype), intent(out), dimension(ysize(2)) :: diss

      real(mytype), dimension(ysize(2)) :: tempa2, tempb2, tempc2

      call transpose_y_to_x(u1p, ta1)
      call derx(tb1, ta1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, zero)
      call transpose_x_to_y(tb1, ta2)
      call transpose_y_to_x(u2p, ta1)
      call derx(tb1, ta1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, zero)
      call transpose_x_to_y(tb1, tb2)
      call horizontal_avrge(ta2 * tb2, tempa2)
      call dery(ta2, u1p, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, zero)
      call dery(tb2, u2p, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, zero)
      call horizontal_avrge(ta2 * tb2, tempb2)
      call transpose_y_to_z(u1p, ta3)
      call derz(tb3, ta3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, zero)
      call transpose_z_to_y(tb3, ta2)
      call transpose_y_to_z(u2p, ta3)
      call derz(tb3, ta3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, zero)
      call transpose_z_to_y(tb3, tb2)
      call horizontal_avrge(ta2 * tb2, tempc2)
      diss = -two * xnu * (tempa2 + tempb2 + tempc2)
   end subroutine dissipation
   
   function comp_thetad_II(ux2, uy2, ux2m) result(thetad)
      use var, only: di2
      use MPI

      implicit none
      real(mytype), intent(in), dimension(ysize(1), ysize(2), ysize(3)) :: ux2, uy2
      real(mytype), intent(in), dimension(ysize(2)) :: ux2m

      real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: ux2p
      real(mytype), dimension(ysize(2)) :: uxuy2pm, ux2mm
      real(mytype), dimension(ysize(2)) :: ydudy2m, dudy2m
      real(mytype) :: thetad, RHS,GT,FT,int_GT,theta,ET,ZTnp1,Disp,Fy_Tot
      integer :: j, code

      integer, parameter :: freq = 1


      if (itime == ifirst .and. irestart == 0) then 
         ZTn = 0.0
         theta_a = 0.0
         theta_b = 0.0
         theta_c = 0.0
      end if

      ! Calculate averaged derivatives
      call extract_fluctuat(ux2, ux2m, ux2p)
      call horizontal_avrge(ux2p * uy2, uxuy2pm)
      call dery(dudy2m, ux2m, di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, zero)
      
      
       ! G(t) Model based on (0: Momentum Thickness)
      if (jthickness == 0) then
         ! Calculate quantities         
         theta = sum((ux2m * (one - ux2m)) * ypw)
         ET = one - theta
         int_GT = sum((dudy2m*(xnu * dudy2m - uxuy2pm))*ypw)
         GT = 2.0 * int_GT - xnu * dudy2m(1)

         ! For DNS simulation
         if (ilesmod == 0) then  
            if (Method_FT ==0) then
               FT = (-K_theta * ET + GT)/theta
            else if (Method_FT ==1) then
               theta_a = ET;
               ZTnp1 = ZTn + adt(1) * theta_a + bdt(1) * theta_b + cdt(1) * theta_c
               FT = (-2.0*K_theta * ET + GT - ((K_theta ** 2) * ZTnp1))/theta
               ! Saving for next time step
               ZTn = ZTnp1
               theta_c = theta_b; theta_b = theta_a
            end if
         else ! For LES simulation
            theta_a = ET;
            ZTnp1 = ZTn + adt(1) * theta_a + bdt(1) * theta_b + cdt(1) * theta_c
            FT = (-2.0*K_theta * ET + GT - ((K_theta ** 2) * ZTnp1))/theta
            ! Saving for next time step
            ZTn = ZTnp1
            theta_c = theta_b; theta_b = theta_a
         end if

      ! G(t) Model based on (1: Displacement Thickness)
      else if (jthickness == 1) then
         Disp  = sum((one - ux2m) * ypw)
         ET = one - Disp
         GT = xnu * dudy2m(1)

         ! For DNS simulation
         if (ilesmod == 0) then  
            if (Method_FT ==0) then
               FT = (-K_theta * ET + GT)/theta
            else if (Method_FT ==1) then
               theta_a = ET;
               ZTnp1 = ZTn + adt(1) * theta_a + bdt(1) * theta_b + cdt(1) * theta_c
               FT = (-2.0*K_theta * ET + GT - ((K_theta ** 2) * ZTnp1))/theta
               ! Saving for next time step
               ZTn = ZTnp1
               theta_c = theta_b; theta_b = theta_a
            end if
         else ! For LES simulation
            theta_a = ET;
            ZTnp1 = ZTn + adt(1) * theta_a + bdt(1) * theta_b + cdt(1) * theta_c
            FT = (-2.0*K_theta * ET + GT - ((K_theta ** 2) * ZTnp1))/theta
            ! Saving for next time step
            ZTn = ZTnp1
            theta_c = theta_b; theta_b = theta_a
         end if

      end if
      
      thetad = FT
      call MPI_BCAST(thetad, 1, real_type, 0, MPI_COMM_WORLD, code)
      
      if ((itime >= ifirst + 2 .or. irestart == 1) .and. mod(itime, freq) == 0) then
         ! total force require to keep Theta = 1
         Fy_Tot = sum(thetad * ypw* dudy2m)
         call MPI_BCAST(Fy_Tot, 1, real_type, 0, MPI_COMM_WORLD, code)
         if (nrank == 0) then
            if (mod(itime, ilist) == 0) then
               write (6, "(' F(t) =',F14.12,'  &   Total Force(y) =',F14.12)") thetad,Fy_Tot
               open (unit=67, file='out/Force.dat', status='unknown', form='formatted', action='write', position='append')
               if (itime == ilist) write (67, "(2A20)") 'F(T)','Total_Force(y)' 
                  write (67, "(2E20.12)") thetad,Fy_Tot
               close (67)
            end if
         end if
      end if
      

   end function

!############################################################################
!############################################################################
!********************************************************************
!*** Added by Pasha
   subroutine tbl_flrt_Check (ux1)
    
        USE decomp_2d
        USE decomp_2d_poisson
        USE variables
        USE param
        USE MPI
    
        implicit none
        real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1
        real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2
    
        integer :: j,i,k,code
        real(mytype) :: can,ut1,ut2
    
        call transpose_x_to_y(ux1,ux2)
        ! Flow rate at the inlet
        ut1=zero
        if (ystart(1)==1) then !! CPUs at the inlet
          do k=1,ysize(3)
            do j=1,ysize(2)-1
              ut1=ut1+(yp(j+1)-yp(j))*(ux2(1,j+1,k)-half*(ux2(1,j+1,k)-ux2(1,j,k)))
            enddo
          enddo
        endif
        call MPI_ALLREDUCE(MPI_IN_PLACE,ut1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
        ut1=ut1/real(nz,mytype) !! Volume flow rate per unit spanwise dist
        ! Flow rate at the outlet
        ut2=zero
        if (yend(1)==nx) then !! CPUs at the outlet
          do k=1,ysize(3)
            do j=1,ysize(2)-1
              ut2=ut2+(yp(j+1)-yp(j))*(ux2(ysize(1),j+1,k)-half*(ux2(ysize(1),j+1,k)-ux2(ysize(1),j,k)))
            enddo
          enddo
        endif
    
        call MPI_ALLREDUCE(MPI_IN_PLACE,ut2,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
        ut2=ut2/real(nz,mytype) !! Volume flow rate per unit spanwise dist
        if ((nrank==0).and.(mod(itime,ilist)==0)) then
          write(*,"(' Mass balance: L-BC, R-BC, diff :',2f12.6,f12.9)") ut1,ut2,abs(ut2-ut1)
        endif

      end subroutine tbl_flrt_Check

end module ptbl
