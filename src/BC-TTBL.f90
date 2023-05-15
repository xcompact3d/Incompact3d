!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module ttbl

   use decomp_2d
   use decomp_2d_io
   use variables
   use param
   use MPI

   implicit none

   integer :: FS!,rc1 ! rc1 is the point for periodicity (idea)
   character(len=100) :: fileformat
   character(len=1), parameter :: NL = char(10) !new line character
   logical :: if_irrot
   character(len=*), parameter :: io_ttbl = "io-ttbl"

   !real(mytype),save,allocatable,dimension(:,:,:) :: lm1,lm2
   real(mytype), save, allocatable, dimension(:, :) :: phisxz, phipphipsxz, upphipsxz, vpphipsxz, wpphipsxz, dphidysxz
   real(mytype), save, allocatable, dimension(:, :) :: enssxz, psxz, usxz, vsxz, wsxz, upupsxz, enspenspsxz
   real(mytype), save, allocatable, dimension(:, :) :: vpvpsxz, wpwpsxz, upvpsxz, vpwpsxz, upwpsxz
   real(mytype), save, allocatable, dimension(:, :) :: tkesxz, Isxz, epssxz
   real(mytype), save :: thetad, thetad_target

   real(mytype), save, allocatable, dimension(:, :) :: jet_mask1, jet_mask2

   real(mytype), save, allocatable, dimension(:) :: ttd1a, ttd1b, ttd1c
   !real(mytype),save, allocatable, dimension(:) :: ttd2a,ttd2b,ttd2c
   real(mytype), save, allocatable, dimension(:) :: ttd3a, ttd3b, ttd3c
   private
   public :: init_ttbl, boundary_conditions_ttbl, momentum_forcing_ttbl, scalar_forcing_ttbl, postprocess_ttbl_init, postprocess_ttbl, visu_ttbl_init, visu_ttbl

contains

   subroutine init_ttbl(ux1, uy1, uz1, ep1, phi1)
      implicit none
      real(mytype), intent(inout), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ep1
      real(mytype), intent(inout), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1

      real(mytype) :: x, y, z, aform
      real(mytype), dimension(ysize(2)) :: um
      integer :: k, j, i, ii, is, code
      integer(kind=MPI_OFFSET_KIND) :: disp
      integer, dimension(:), allocatable :: seed

      ux1 = zero; uy1 = zero; uz1 = zero
      if (iin /= 0) then
         call system_clock(count=code)
         if (iin == 2) code = 0
         call random_seed(size=ii)
         call random_seed(put=code + 63946 * (nrank + 1) * (/(i - 1, i=1, ii)/))
         !call random_seed(put = code+63946*nrank*(/ (i - 1, i = 1, ii) /))
         call random_number(ux1)
         call random_number(uy1)
         call random_number(uz1)
      end if

      if (iin == 3) then

         if (nrank == 0) write (*, *) 'READING FILES FOR INITIAL CONDITION'

         !  call decomp_2d_read_one(1,ux1,'ux.ic')
         !  call decomp_2d_read_one(1,uy1,'uy.ic')
         !  call decomp_2d_read_one(1,uz1,'uz.ic')

      else

         aform = sqrt(0.10922670d0 / 2.0d0)
         if (nrank == 0) write (*, *) 'INITIAL CONDITION WITH GOLDEN RATION', aform, init_noise
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

      end if

      phi1(:, :, :, :) = zero

      return
   end subroutine init_ttbl
   !############################################################################
   !############################################################################
   subroutine boundary_conditions_ttbl(ux, uy, uz, phi, ep)

      use var, only: ux2, uy2, uz2, ux3, uy3, uz3
      use var, only: ta1, ta2, tb1, tb2, tc1, tc2
      use var, only: td2, td3, te2, te3, tf2, tf3
      use var, only: tg1, tg2, ti1, ti2, th1, th2
      use var, only: di1, di3
      use ibm_param, only: ubcx, ubcy, ubcz
      use navier, only: tbl_flrt

      implicit none
      real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux, uy, uz, ep
      real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi
      real(mytype) :: dyp2, cx, x, z, r0, r1, r2, r3, r4, r5, r6, r7, r8, D, rdx, r, dzz, dxx
      integer :: i, j, k, rc2, xx, zz

      if (itime == ifirst) then ! static BCs just computed first time step !
         phi(:, :, :, :) = zero

         if (ncly1 == 2) then ! Bottom Boundary

            do k = 1, xsize(3)
               do i = 1, xsize(1)
                  byx1(i, k) = zero
                  byy1(i, k) = zero
                  byz1(i, k) = zero
               end do
            end do

            do k = 1, ysize(3) ! semi implicit condition
               do i = 1, ysize(1)
                  byx1_2(i, k) = zero
                  byy1_2(i, k) = zero
                  byz1_2(i, k) = zero
               end do
            end do

            if (inflow_noise == 4) then  ! steady jets in crossflow via Dirichlet
               allocate (jet_mask1(xsize(1), xsize(3)))
               allocate (jet_mask2(ysize(1), ysize(3)))
               D = 3
               rdx = 0.7_mytype ! peplinkski 2015
               R = half
               dxx = 5 * D
               dzz = 2 * D
               tg1(:, :, :) = zero
               do k = 1, xsize(3)
                  z = real((k + xstart(3) - 2), mytype) * dz
                  do i = 1, xsize(1)
                     x = real(i - 1, mytype) * dx
                     do xx = 1, (int(xlx / dxx) + 2) ! centre x
                        do zz = 1, (int(zlz / dzz) + 2) ! centre y
                           r = (two / D) * sqrt((x - ((xx - 1) * dxx))**2 + (z - ((zz - 1) * dzz))**2)
                           tg1(i, :, k) = tg1(i, :, k) + abs((two * R) * (one - r**2) * exp(-(r / rdx)**4))
                        end do
                     end do
                  end do
               end do
               do k = 1, xsize(3)
                  do i = 1, xsize(1)
                     jet_mask1(i, k) = tg1(i, 1, k)
                  end do
               end do
               call transpose_x_to_y(tg1, tg2)
               do k = 1, ysize(3)
                  do i = 1, ysize(1)
                     jet_mask2(i, k) = tg2(i, 1, k)
                  end do
               end do
            end if ! inflow_noise == 4 ! jets in crossflow via Dirichlet

         end if ! ncly1 == 2

    !! Top Boundary
         if (nclyn == 2) then
            do k = 1, xsize(3)
               do i = 1, xsize(1)
                  byxn(i, k) = one
                  byyn(i, k) = zero
                  byzn(i, k) = zero
               end do
            end do
            do k = 1, ysize(3) ! semi implicit condition
               do i = 1, ysize(1)
                  byxn_2(i, k) = one
                  byyn_2(i, k) = zero
                  byzn_2(i, k) = zero
               end do
            end do
         end if ! nclyn == 2

      end if ! itime.eq.ifirst

      if (inflow_noise == 4) then  ! steady jets in crossflow via Dirichlet
         if (ncly1 == 2) then
            do k = 1, xsize(3)
               do i = 1, xsize(1)
                  byy1(i, k) = jet_mask1(i, k)!*half+half*tanh((t-ten))
               end do
            end do
            !call transpose_x_to_y(tg1,tg2)
            do k = 1, ysize(3)
               do i = 1, ysize(1)
                  byy1_2(i, k) = jet_mask2(i, k)!*half+half*tanh((t-ten))
               end do
            end do
         end if
      end if

      !if (iscalar==1.and.nclyS1.eq.2.and.xstart(2).eq.1) then
      !   do k=1,xsize(3)
      !      do i=1,xsize(1)
      !      !phi2(i,yend(2),k)= phi2(i,yend(2)-1,k) / (1.+uset(is)*dy*sc(is)/xnu) !Robin on top BC
      !      phi(i,1,k,1) = one - (jet_mask1(i,k)/maxval(jet_mask1(:,:))) !phi1(i,j+1,k,is)! dc/dn=0
      !      enddo
      !   enddo
      !endif

      ! if (nclx1 == 2) then ! new spatial option

      ! recirculation inflow condition
      !    do k=1,xsize(3)
      !       do j=1,xsize(2)
      !          bxx1(j,k)=ux(rc1,j,k)
      !          bxy1(j,k)=uy(rc1,j,k)
      !          bxz1(j,k)=uz(rc1,j,k)
      !          if (iscalar==1) phi(1,j,k,:)=zero
      !       enddo
      !    enddo

      !subrecirculation zone to ensure aperiodicity
      !rc2 = int(rc1/1.61803398875d0)
      !do k=1,xsize(3)
      !   do j=1,xsize(2)
      !      ux(rc2,j,k)=ux(rc2,j,k)
      !      uy(rc2,j,k)=uy(rc2,j,k)
      !      uz(rc2,j,k)=uz(rc2,j,k)
      !   enddo
      !enddo

      !OUTFLOW based on a 1D convection equation
      !if (nclxn == 2) then
      !   do k=1,xsize(3)
      !      do j=1,xsize(2)
      !         cx=ux(nx,j,k)*gdt(itr)*(one/dx)
      !         if(cx.lt.0.0) cx=zero
      !         bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
      !         bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
      !         bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
      !         if (iscalar==1)phi(nx,j,k,:)=phi(nx,j,k,:)-cx*(phi(nx,j,k,:)-phi(nx-1,j,k,:))
      !      enddo
      !   enddo
      !endif

       !! Top Boundary
      !if (nclyn == 2) then
      !   do k = 1, xsize(3)
      !      do i = 1, xsize(1)
      !         byxn(i, k) = ux(i, xsize(2) - 1, k)
      !         byyn(i, k) = uy(i, xsize(2) - 1, k)
      !         byzn(i, k) = uz(i, xsize(2) - 1, k)
      !      enddo
      !   enddo
      !   do k = 1, ysize(3)
      !      do i = 1, ysize(1)
      !         byxn_2(i, k) = ux(i, xsize(2) - 1, k)
      !         byyn_2(i, k) = uy(i, xsize(2) - 1, k)
      !         byzn_2(i, k) = uz(i, xsize(2) - 1, k)
      !      enddo
      !   enddo
      !endif

      !update of the flow rate (what is coming in the domain is getting out)
      !call tbl_flrt(ux,uy,uz)

      !endif

      !  if(inflow_noise==5)then  ! irrotational top BC

      !     call derxx (ta1,ux,di1,sx,sfx ,ssx ,swx ,xsize(1),xsize(2),xsize(3),0,ubcx)
      !     call derxx (tb1,uy,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1,ubcy)
      !     call derxx (tc1,uz,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1,ubcz)

      !     call transpose_x_to_y(ta1,ta2)
      !     call transpose_x_to_y(tb1,tb2)
      !     call transpose_x_to_y(tc1,tc2)

      !     call transpose_x_to_y(ux,ux2)
      !     call transpose_x_to_y(uy,uy2)
      !     call transpose_x_to_y(uz,uz2)

      !     call transpose_y_to_z(ux2,ux3)
      !     call transpose_y_to_z(uy2,uy3)
      !     call transpose_y_to_z(uz2,uz3)

      !     call derzz (td3,ux3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1,ubcx)
      !     call derzz (te3,uy3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1,ubcy)
      !     call derzz (tf3,uz3,di3,sz,sfz ,ssz ,swz ,zsize(1),zsize(2),zsize(3),0,ubcz)

      !     call transpose_z_to_y(td3,td2)
      !     call transpose_z_to_y(te3,te2)
      !     call transpose_z_to_y(tf3,tf2)

      !     !new diverge free BC to allow vertical outflow and minimize pressure bump
      !     !must compute in Y pencil on a 3D field and transpose field to X to extract plane
      !     tg2(:,:,:)=zero
      !     th2(:,:,:)=zero
      !     ti2(:,:,:)=zero
      !     !dyp2 = dyp(ny)*dyp(ny-1); print *, dyp2
      !     dyp2 = (yp(ny)-yp(ny-1))*(yp(ny-1)-yp(ny-2))!; print *, dyp2

      !     do k=1,ysize(3)
      !        do j=1,ysize(2)
      !           do i=1,ysize(1)
      !              tg2(i,j,k)= (-ta2(i,ny-1,k) -td2(i,ny-1,k))*dyp2 +ux2(i,ny-1,k) -ux2(i,ny-2,k)
      !              th2(i,j,k)= (-tb2(i,ny-1,k) -te2(i,ny-1,k))*dyp2 +uy2(i,ny-1,k) -uy2(i,ny-2,k)
      !              ti2(i,j,k)= (-tc2(i,ny-1,k) -tf2(i,ny-1,k))*dyp2 +uz2(i,ny-1,k) -uz2(i,ny-2,k)
      !           enddo
      !        enddo
      !     enddo

      !     ! saving upper BC planes *_2 in a lapis Y variable to apply BCs at iimplicit.f90 in variable bctop
      !     do k=1,ysize(3)
      !        do i=1,ysize(1)
      !           byxn_2(i, k)=tg2(i,ny,k)
      !           byyn_2(i, k)=th2(i,ny,k)
      !           byzn_2(i, k)=ti2(i,ny,k)
      !        enddo
      !     enddo

      !     call transpose_y_to_x(tg2,tg1)
      !     call transpose_y_to_x(th2,th1)
      !     call transpose_y_to_x(ti2,ti1)

      !     do k = 1, xsize(3)
      !        do i = 1, xsize(1)
      !           byxn(i, k) = tg1(i,ny,k)
      !           byyn(i, k) = th1(i,ny,k)
      !           byzn(i, k) = ti1(i,ny,k)
      !        enddo
      !     enddo

      !     !if(nrank.eq.0)then
      !        !  write(6,*) nrank,minval(tg2(:,:,:)),maxval(tg2(:,:,:))
      !        !  write(6,*) nrank,minval(byxn_2(:,:)),maxval(byxn_2(:,:))
      !        !  write(6,*) nrank,minval(tg1(:,:,:)),maxval(tg1(:,:,:))
      !     !   write(6,*) 'ux',minval(byxn(:,:)),maxval(byxn(:,:))
      !     !   write(6,*) 'uy',minval(byyn(:,:)),maxval(byyn(:,:))
      !     !   write(6,*) 'uz',minval(byzn(:,:)),maxval(byzn(:,:))
      !     !endif

      !  endif ! inflow_noise==5 ! irrotational top BC

      !du/dy = -lambda*u !lambda <<1 Newman; lambda >> 1 Dirichlet; lambda = 1/Ly -> mixed
      !byxn(i,k)=ux(:,ny-1,:) - dy*lambda*ux(:,ny-1,:) = (1-lambda*dy)*ux(:,ny-1,:)

   end subroutine boundary_conditions_ttbl
   !############################################################################
   !############################################################################
   subroutine comp_thetad(ux2, uy2, ux2m)
      use var, only: ta1, tb1, tc1, ta2, tb2, tc2, di2
      use ibm_param, only: ubcx, ubcy, ubcz
      use MPI
      implicit none
      integer :: i, j, k, code, maxit
      real(mytype), intent(in), dimension(ysize(1), ysize(2), ysize(3)) :: ux2, uy2
      real(mytype), intent(in), dimension(ysize(2)) :: ux2m
      real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: ux2p, dudy2
      real(mytype), dimension(ysize(2)) :: uxuy2pm, ux2mm
      real(mytype), dimension(ysize(2)) :: ydudy2m, dudy2m, du2dy22m, duxuy2pm
      real(mytype) :: tol, resi, thetad1, thetad2, thetad3, theta1, theta3, theta2
      !real(8)      :: tstart
      !tstart=MPI_WTIME()
      logical reset

      call extract_fluctuat(ux2, ux2m, ux2p) ! ux'
      call horizontal_avrge(ux2p * uy2, uxuy2pm) ! ux'uy' in profile
      call dery(duxuy2pm, uxuy2pm, di2, sy, ffy, fsy, fwy, ppy, 1, ysize(2), 1, 0, ubcy) !0*1=0
      call dery(dudy2m, ux2m, di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, ubcy)
      call dery(ydudy2m, ux2m * yp, di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, ubcy)
      !ydudy2m = ydudy2m - ux2m !ydudy = (duydy+u) - u = ydudy

      if (iimplicit <= 0) then
         call deryy(du2dy22m, ux2m, di2, sy, sfyp, ssyp, swyp, 1, ysize(2), 1, 1, ubcy)
         if (istret /= 0) then
            do j = 1, ysize(2)
               du2dy22m(j) = du2dy22m(j) * pp2y(j) - pp4y(j) * dudy2m(j)
            end do
         end if
      else ! (semi)implicit Y diffusion
         do j = 1, ysize(2)
            du2dy22m(j) = -pp4y(j) * dudy2m(j)
         end do
      end if

      !DEBUG
      !call outpdt(ux2m    ,'ux2m')
      !call outpdt(dudy2m  ,'dudy2m')
      !call outpdt(du2dy22m,'du2dy22m')
      !call outpdt(duxuy2pm,'duxuy2pm')
      !ux2m(:) = ux2m_old(:)+dt*(thetad1*yp(:)*dudy2m(:)+xnu*du2dy22m(:)-duxuy2pm(:))
      !theta1 = sum((ux2m*(one-ux2m))*ypw)
      !call outpdt(ux2m    ,'ux2m_theta1')
      !ux2m(:) = ux2m_old(:)+dt*(thetad2*yp(:)*dudy2m(:)+xnu*du2dy22m(:)-duxuy2pm(:))
      !theta2 = sum((ux2m*(one-ux2m))*ypw)
      !call outpdt(ux2m    ,'ux2m_theta2')

#ifdef DOUBLE_PREC
      tol = 1.0e-12
      maxit = 50
#else
      tol = 1.0e-7
      maxit = 25
#endif

      if (itime > (ifirst + 3)) then !skipping Euler and AB2

         if (nrank == 0 .and. mod(itime, 4) == 0) then

            ! get the intervals for bisection ! do not go to wide
            thetad = thetad !+! dt*(theta-one)
            i = 0; resi = one; 
            thetad1 = thetad - dt !1e-2
            thetad2 = thetad + dt !1e-2
            reset = .false.

            do while (abs(resi) >= tol) ! error >= tolerance

               i = i + 1
               if (i > maxit) then
                  write (6, *) 'skip thetad after', i, 'thtd', thetad3, 'tht3=', theta3; exit
               end if

               ! remember here : if (itimescheme.eq.3) then ! AB3
               !  adt(1)= (twentythree/twelve)*dt
               !  bdt(1)=-(    sixteen/twelve)*dt
               !  cdt(1)= (       five/twelve)*dt

               ttd1a(:) = thetad1 * ydudy2m(:) + xnu * du2dy22m(:) - duxuy2pm(:)
               ux2mm(:) = ux2m(:) + adt(1) * ttd1a(:) + bdt(1) * ttd1b(:) + cdt(1) * ttd1c(:)
               theta1 = sum((ux2mm * (one - ux2mm)) * ypw)

               thetad3 = half * (thetad1 + thetad2)
               ttd3a(:) = thetad3 * ydudy2m(:) + xnu * du2dy22m(:) - duxuy2pm(:)
               ux2mm(:) = ux2m(:) + adt(1) * ttd3a(:) + bdt(1) * ttd3b(:) + cdt(1) * ttd3c(:)
               theta3 = sum((ux2mm * (one - ux2mm)) * ypw)

               ! flip flip flip ! key point
               if ((theta1 - one) * (theta3 - one) < zero) then
                  thetad2 = thetad3
               else
                  thetad1 = thetad3
               end if

               !   if(theta1.gt.1.1)then
               !    write(6,*) 'RESET RESET RESET'
               !    reset=.true.
               !   endif
               !   if(reset)then
               !   thetad1=1e-5
               !   thetad2=0.20d0
               !   reset=.false.
               !  endif
               resi = theta1 - one
               !if(resi.lt.zero)then ! incrase bounds
               !  thetad1=thetad1-dt
               !  thetad2=thetad2+dt
               !endif
               !write(6,*)'residual',resi,thetad3

            end do

            if (abs(resi) < tol) then
               write (6, *) 'theta dot computed in', i, 'with', thetad3, 'theta3=', theta3
               thetad_target = thetad3
               !print *,'Time computing theta dot (s)', real(MPI_WTIME()-tstart,4)
            end if

         end if !nrank.eq.0

         call MPI_BCAST(thetad_target, 1, real_type, 0, MPI_COMM_WORLD, code)

      else !itime.gt.ifirst+3

         !all the same in first time-steps
         thetad1 = thetad; thetad2 = thetad
         ttd1a(:) = thetad1 * ydudy2m(:) + xnu * du2dy22m(:) - duxuy2pm(:)
         !ttd2a(:) = thetad2*ydudy2m(:)+xnu*du2dy22m(:)-duxuy2pm(:)
         ttd3a(:) = ((thetad1 + thetad2) * 0.5) * ydudy2m(:) + xnu * du2dy22m(:) - duxuy2pm(:)

      end if

      ttd1c(:) = ttd1b(:); ttd1b(:) = ttd1a(:)
      !ttd2c(:)=ttd2b(:);ttd2b(:)=ttd2a(:)
      ttd3c(:) = ttd3b(:); ttd3b(:) = ttd3a(:)

      return
   end subroutine comp_thetad
   !############################################################################
   !############################################################################
   subroutine momentum_forcing_ttbl(dux1, duy1, duz1, ux1, uy1, uz1, phi1)
      use var, only: ta1, tb1, tc1, ux2, uy2, uz2, ta2, tb2, tc2, di2
      use ibm_param, only: ubcx, ubcy, ubcz

      implicit none
      real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
      real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: dux2, duy2, duz2
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1
      real(mytype), dimension(ysize(2)) :: ux2m
      real(mytype) :: y, ttheta, delay
      integer :: i, j, k, code

      real(mytype), dimension(ysize(2)) :: ydudy2m, dudy2m, du2dy22m, duxuy2pm
      real(mytype) :: cf_ref, delta, tau_wall, tau_wall1, ufric, temp, utau
      if (itime > ifirst) then ! we need to have start postprocess_init

         call transpose_x_to_y(ux1, ux2)
         call transpose_x_to_y(uy1, uy2)
         call horizontal_avrge(ux2, ux2m)

         ttheta = sum((ux2m * (one - ux2m)) * ypw)
         call MPI_BCAST(ttheta, 1, real_type, 0, MPI_COMM_WORLD, code)
         delta = sum((one - ux2m) * ypw)
         call MPI_BCAST(delta, 1, real_type, 0, MPI_COMM_WORLD, code)

         if (itime == (ifirst + 1)) then

            thetad = 0.1092267 * xnu ! laminar value
            if (inflow_noise > 0 .and. inflow_noise < 1) thetad = inflow_noise ! we need to provide the value for restarting
            if (nrank == 0) write (*, "(' initializing thetad= ',E14.7)") thetad

         else

            ! Option 1: Damien paper
            call comp_thetad(ux2, uy2, ux2m) ! bisection method ! updates thetad_target
            thetad = thetad_target

            ! Option 1.1: Average to filter oscillations
            ! thetad = half*(thetad+thetad_target)

            ! Option 1.2: Delayed to filter strong oscillations
            ! delay = dt
            ! thetad = thetad*(one-delay) + delay*thetad_target

            ! Option 2: do not need the bisection !
            ! thetad = thetad + dt*(ttheta-one)
            ! if(nrank.eq.0)write(6,*) 'simple thetad method =',thetad

            ! broadcast final value
            call MPI_BCAST(thetad, 1, real_type, 0, MPI_COMM_WORLD, code)

         end if
         if (mod(itime, ilist) == 0 .or. itime == (ifirst + 1)) then

            call transpose_x_to_y(uz1, uz2)
            call dery(ta2, ux2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, ubcx)
            call dery(tb2, uy2, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0, ubcy)
            call dery(tc2, uz2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, ubcz)
            call transpose_y_to_x(ta2, ta1)
            call transpose_y_to_x(tb2, tb1)
            call transpose_y_to_x(tc2, tc1)

            !temp = sum(ta2(1:rc1,1,:)) ! when adding inflow/outflow
            temp = sum(ta2(:, 1, :))
            call MPI_ALLREDUCE(temp, tau_wall1, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
            !tau_wall=tau_wall1/real(rc1*nz,mytype) ! when adding inflow/outflow
            tau_wall = tau_wall1 / real(nx * nz, mytype)
            ufric = sqrt(tau_wall * xnu)

            if (nrank == 0) then

               write (6, "(' thetad= ',E13.4,' theta= ',F14.10,' H= ',F8.6)") thetad, ttheta, delta / ttheta
               write (6, "(' tau_wall=',F15.7,' u_tau=',F15.7,' cf=',F15.7)") tau_wall, ufric, 2 * ufric**2
               write (6, "(' Re_theta=',F15.7,' Re_delta=',F15.7)") 1.0 / (ttheta * xnu), 1.0 / (xnu / delta)
               write (6, "(' dx+=',F12.6,' dz+=',F12.6,' dy+min,max=',2F12.6)") dx * ufric * re, dz * ufric * re, dyp(1) * ufric * re, dyp(ny) * ufric * re

               ! outposting to file ttbl_scales.dat
               FS = 7; write (fileformat, '( "(",I4,"(E15.7),A)" )') FS; FS = FS * 15 + 1
               open (unit=67, file='ttbl_scales.dat', status='unknown', form='formatted', recl=FS, action='write', position='append')
               write (67, fileformat) t, tau_wall, ufric, 2 * ufric**2, thetad, ttheta, delta; close (67)

               ! outposting to file ttbl_scales2.dat (optional, just for eventual plotting)
               FS = 6; write (fileformat, '( "(",I4,"(E13.4),A)" )') FS; FS = FS * 15 + 1
               open (unit=67, file='ttbl_scales2.dat', status='unknown', form='formatted', recl=FS, action='write', position='append')
               write (67, fileformat) t, dx * ufric * re, dz * ufric * re, dyp(1) * ufric * re, dyp(ny) * ufric * re, dt * ufric * re; close (67)

            end if ! nrank.eq.0

         end if

         if (thetad > 0) then
            do k = 1, xsize(3)
               do j = 1, xsize(2)
                  do i = 1, xsize(1)
                     if (istret == 0) y = (j + xstart(2) - 1 - 1) * dy
                     if (istret /= 0) y = yp(j + xstart(2) - 1)
                     dux1(i, j, k, 1) = dux1(i, j, k, 1) + thetad * y * ta1(i, j, k)!*lm1(i,j,k) !+ dpdx*(one-ux2m(j))
                     duy1(i, j, k, 1) = duy1(i, j, k, 1) + thetad * y * tb1(i, j, k)!*lm1(i,j,k)
                     ! if active scalar field on -- active coupling
                     if (iscalar == 1 .and. ri(1) /= 0) duy1(i, j, k, 1) = duy1(i, j, k, 1) + ri(1) * phi1(i, j, k, 1)
                     duz1(i, j, k, 1) = duz1(i, j, k, 1) + thetad * y * tc1(i, j, k)!*lm1(i,j,k)
                  end do
               end do
            end do
         end if
      end if

      return
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

      if (iscalar == 1 .and. ri(1) /= 0 .and. itime > ifirst) then

         call transpose_x_to_y(phi1(:, :, :, 1), phi2(:, :, :, 1))
         call horizontal_avrge(phi2(:, :, :, 1), phi2m)
         omega = sum(phi2m * ypw)
         omegad = sc(1) * thetad
         if (nrank == 0) write (*, "(' omegad= ',E15.7,' omega= ',F14.12)") omegad, omega

         if (omegad > 0) then

            call deryS(ta2, phi2(:, :, :, 1), di2, sy, ffypS, fsypS, fwypS, ppy, ysize(1), ysize(2), ysize(3), 1, zero)
            call transpose_y_to_x(ta2, ta1)

            do k = 1, xsize(3)
               do j = 1, xsize(2)
                  do i = 1, xsize(1)
                     if (istret == 0) y = (j + xstart(2) - 1 - 1) * dy
                     if (istret /= 0) y = yp(j + xstart(2) - 1)
                     dphi1(i, j, k, 1) = dphi1(i, j, k, 1) + omegad * y * ta1(i, j, k)
                     if (ri(1) /= 0) dphi1(i, j, k, 1) = dphi1(i, j, k, 1) + uy1(i, j, k)!(*dphi2/dy==1) forcing Brunt-Vaysala cte
                  end do
               end do
            end do

         end if

      end if

      return
   end subroutine scalar_forcing_ttbl
   !############################################################################
   !############################################################################
   subroutine postprocess_ttbl_init
      implicit none
      real(mytype) :: dxdydz
      integer :: i, j, k, code, ierr
      character :: a

      ! Simpson method integral operator !
      !  call alloc_x(vol1, opt_global=.true.)
      !  vol1 = zero
      !  dxdydz=dx*dy*dz
      !  do k=xstart(3),xend(3)
      !     do j=xstart(2),xend(2)
      !        do i=xstart(1),xend(1)
      !           vol1(i,j,k)=dxdydz
      !           if (i .eq. 1.or. i .eq. nx) vol1(i,j,k) = vol1(i,j,k) * (five/twelve)
      !           if (j .eq. 1.or. j .eq. ny) vol1(i,j,k) = vol1(i,j,k) * (five/twelve)
      !           if (k .eq. 1.or. k .eq. nz) vol1(i,j,k) = vol1(i,j,k) * (five/twelve)
      !           if (i .eq. 2.or. i .eq. nx-1) vol1(i,j,k) = vol1(i,j,k) * (thirteen/twelve)
      !           if (j .eq. 2.or. j .eq. ny-1) vol1(i,j,k) = vol1(i,j,k) * (thirteen/twelve)
      !           if (k .eq. 2.or. k .eq. nz-1) vol1(i,j,k) = vol1(i,j,k) * (thirteen/twelve)
      !        end do
      !     end do
      !  end do

      ! Check for valid output parameters
      if (ioutput < ilist .or. mod(ioutput, ilist) /= 0) then
         if (nrank == 0) print*,"ioutput must be exactly divisible by ilist"
         call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
      end if

      if (nrank == 0) write (6, *) 'Allocating statistics variables!', ioutput / ilist
      !allocate(enssxz(ysize(2),      int(ioutput/ilist)))
      allocate (psxz(ysize(2), int(ioutput / ilist)))
      allocate (usxz(ysize(2), int(ioutput / ilist)))
      allocate (vsxz(ysize(2), int(ioutput / ilist)))
      allocate (wsxz(ysize(2), int(ioutput / ilist)))
      allocate (upupsxz(ysize(2), int(ioutput / ilist)))
      allocate (enspenspsxz(ysize(2), int(ioutput / ilist)))
      allocate (vpvpsxz(ysize(2), int(ioutput / ilist)))
      allocate (wpwpsxz(ysize(2), int(ioutput / ilist)))
      allocate (upvpsxz(ysize(2), int(ioutput / ilist)))
      allocate (vpwpsxz(ysize(2), int(ioutput / ilist)))
      allocate (upwpsxz(ysize(2), int(ioutput / ilist)))
      allocate (tkesxz(ysize(2), int(ioutput / ilist)))
      !allocate(Isxz(ysize(2),        int(ioutput/ilist)))
      allocate (epssxz(ysize(2), int(ioutput / ilist)))
      allocate (phisxz(ysize(2), int(ioutput / ilist)))
      allocate (phipphipsxz(ysize(2), int(ioutput / ilist)))
      allocate (upphipsxz(ysize(2), int(ioutput / ilist)))
      allocate (vpphipsxz(ysize(2), int(ioutput / ilist)))
      allocate (wpphipsxz(ysize(2), int(ioutput / ilist)))
      allocate (dphidysxz(ysize(2), int(ioutput / ilist)))

      allocate (ttd1a(ysize(2)), ttd1b(ysize(2)), ttd1c(ysize(2)))
      !allocate(ttd2a(ysize(2)),ttd2b(ysize(2)),ttd2c(ysize(2))) ! we might not need this intermediary step
      allocate (ttd3a(ysize(2)), ttd3b(ysize(2)), ttd3c(ysize(2)))

      ! Atempt to do a vertically periodic domain !
      !call alloc_x(lm1); lm1(:,:,:)=one
      !call alloc_y(lm2); lm2(:,:,:)=one
      !if(nclx1==0)rc1=nx ! normal case
      !if(nclx1==2)then
      !   rc1 = nx-1
      !   rc1 = int(0.40*nx) ! define a perfecentage here
      !   do i=1,xsize(1)
      !      !x=real(i-1,mytype)*dx
      !      !lambda1(i,j,k) = (half-half*tanh((x-(0.40*xlx))*150.0))
      !      !lambda1(i,j,k) = (half-half*tanh((x-(recirc*xlx))*1000.0))
      !      if(i.gt.rc1) lm1(i,:,:) = zero
      !   end do
      !   call transpose_x_to_y(lm1,lm2)
      !else
      !   rc1=nx
      !endif
      ! NOT GOOD !

      if (nrank == 0) then
         write (6, *) 'Resolution details:'
         write (6, "(' dx  =',F12.6)") dx
         write (6, "(' dz  =',F12.6)") dz
         write (6, "(' dx/dz  =',F12.6)") dx / dz
         write (6, "(' dy  min,max =',2F12.6)") dyp(1), dyp(ny)
      end if
   end subroutine postprocess_ttbl_init
   !############################################################################
   !############################################################################
   subroutine postprocess_ttbl(ux1, uy1, uz1, pp3, phi1, ep1)
      use var, only: uvisu, one, dt, t, ux2, uy2, uz2
      use var, only: pp1, ta1, tb1, tc1, td1, te1, tf1, tg1, th1, ti1, di1, nxmsize
      use var, only: pp2, ta2, tb2, tc2, td2, te2, tf2, ppi2, di2, dip2, ph2, nymsize
      use var, only: ta3, tb3, tc3, td3, te3, tf3, di3, dip3, ph3, nzmsize, phi2
      use var, only: npress
      use ibm_param, only: ubcx, ubcy, ubcz

      implicit none

      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1, ep1
      real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1

      real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux1p, uy1p, uz1p, phi1p, dphi1pdx
      real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: pre2, ux2p, uy2p, uz2p, phi2p, pre2p, ens2p
      real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: dphi2pdx, dphi2pdy, dphi2pdz, dpre2pdy
      real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: duy2pdy
      real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: temp1, tke3d1, ens1, diss1, pre1
      real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: temp2, tke3d2, ens2
      real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: temp3, phi3, dphi3pdz, phi3p

      real(mytype), dimension(ysize(2)) :: sxz, sxz1
      real(mytype), dimension(ysize(2)) :: dup2dxsxz, dvp2dxsxz, dwp2dxsxz
      real(mytype), dimension(ysize(2)) :: dup2dysxz, dvp2dysxz, dwp2dysxz
      real(mytype), dimension(ysize(2)) :: dup2dzsxz, dvp2dzsxz, dwp2dzsxz
      real(mytype), dimension(ysize(2)) :: pre_sxz, dvpppsxz, vpppsxz, vpdppsxz, ppppsxz !pressure transport
      real(mytype), dimension(ysize(2)) :: kurt_x, kurt_y, kurt_z, kurt_phi, skew_x, skew_y, skew_z, skew_phi
      real(mytype), dimension(ysize(2)) :: dupupvpsxz, dvpwpwpsxz, dvpvpvpsxz !turbulent transport
      real(mytype), dimension(ysize(2)) :: dusxz, dwsxz !shear production
      real(mytype), dimension(ysize(2)) :: d2upupsxz, d2vpvpsxz, d2wpwpsxz, dtkesxzdy, d2tkesxzdy !visocus diffusion
      real(mytype), dimension(ysize(2)) :: dvpphipphipsxz !buoyancy turbulent transport
      real(mytype), dimension(ysize(2)) :: d2phipphipsxz, dphipphipsxz !buoyancy diffusion
      real(mytype), dimension(ysize(2)) :: dphisxz !buoyancy production
      real(mytype), dimension(ysize(2)) :: b_eps !buoyancy dissipation
      real(mytype), dimension(ysize(2)) :: dusxzdy, dwsxzdy !dudy and dwdy

      real(mytype) :: mp(numscalar)
      real(8) :: tstart
      integer :: i, j, k, is, code, b
      character(len=1), parameter :: NL = char(10)
      character(len=60) :: filename
      character(len=300) :: fileformat

      if ((itime > initstat) .and. mod(itime, ilist) == 0) then
         tstart = MPI_WTIME()

         b = int(mod(itime, ioutput) / ilist)
         if (mod(itime, ioutput) == 0) b = int(ioutput / ilist)
         !if(nrank.eq.0)write(6,*)'b=',b

         call pressure_x(pp3, pre1) !interpolate to velocity mesh and rescale

         call derx(ta1, ux1, di1, sx, ffx, fsx, fwx, xsize(1), xsize(2), xsize(3), 0, ubcx)
         call derx(tb1, uy1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, ubcy)
         call derx(tc1, uz1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, ubcz)
         !y-derivatives
         call transpose_x_to_y(ux1, td2)
         call transpose_x_to_y(uy1, te2)
         call transpose_x_to_y(uz1, tf2)

         call dery(ta2, td2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, ubcx)
         call dery(tb2, te2, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0, ubcy)
         call dery(tc2, tf2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, ubcz)
       !!z-derivatives
         call transpose_y_to_z(td2, td3)
         call transpose_y_to_z(te2, te3)
         call transpose_y_to_z(tf2, tf3)
         call derz(ta3, td3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, ubcx)
         call derz(tb3, te3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, ubcy)
         call derz(tc3, tf3, di3, sz, ffz, fsz, fwz, zsize(1), zsize(2), zsize(3), 0, ubcz)
       !!all back to x-pencils
         call transpose_z_to_y(ta3, td2)
         call transpose_z_to_y(tb3, te2)
         call transpose_z_to_y(tc3, tf2)
         call transpose_y_to_x(td2, tg1)
         call transpose_y_to_x(te2, th1)
         call transpose_y_to_x(tf2, ti1)
         call transpose_y_to_x(ta2, td1)
         call transpose_y_to_x(tb2, te1)
         call transpose_y_to_x(tc2, tf1)

         ! enstrophy computation
         !  ens1=(tf1-th1)**2+(tg1-tc1)**2+(tb1-td1)**2
         !  call transpose_x_to_y(ens1,ens2)
         !  call horizontal_avrge(ens2,enssxz(:,b)) ! enstrophy
         !  call extract_fluctuat(ens2,enssxz(:,b),ens2p) ! enstrophy'
         !  call horizontal_avrge(ens2p**2,enspenspsxz(:,b)) ! enstrophy' enstrophy'

         call transpose_x_to_y(ux1, ux2)
         call horizontal_avrge(ux2, usxz(:, b)) ! ux
         call extract_fluctuat(ux2, usxz(:, b), ux2p) ! ux'
         call transpose_y_to_x(ux2p, ux1p)

         call transpose_x_to_y(uy1, uy2)
         call horizontal_avrge(uy2, vsxz(:, b)) ! uy
         call extract_fluctuat(uy2, vsxz(:, b), uy2p) ! uy'
         !uy2p = uy2 !damian case
         call transpose_y_to_x(uy2p, uy1p)

         call transpose_x_to_y(uz1, uz2)
         call horizontal_avrge(uz2, wsxz(:, b)) ! uz
         call extract_fluctuat(uz2, wsxz(:, b), uz2p) ! uz'
         call transpose_y_to_x(uz2p, uz1p)

         call horizontal_avrge(ux2p**2, upupsxz(:, b)) ! ux' ux'
         call horizontal_avrge(uy2p**2, vpvpsxz(:, b)) ! uy' uy'
         call horizontal_avrge(uz2p**2, wpwpsxz(:, b)) ! uz' uz'

         call horizontal_avrge(ux2p * uy2p, upvpsxz(:, b)) ! ux' uy'
         call horizontal_avrge(uy2p * uz2p, vpwpsxz(:, b)) ! uy' uz'
         call horizontal_avrge(ux2p * uz2p, upwpsxz(:, b)) ! ux' uz'

         tke3d2 = half * (ux2p**2 + uy2p**2 + uz2p**2)
         call horizontal_avrge(tke3d2, tkesxz(:, b)) ! tke

         !OUTPOST TKE FIELD TO DISK start (OPTIONAL)
         !call transpose_y_to_x(tke3d2,tke3d1)
         !call tke(tke3d1*vol1)
         !if (mod(itime,ioutput).eq.0) then
         ! write(filename,"('./data/tke',I4.4)") itime/ioutput
         ! call decomp_2d_write_one(1,tke3d1uvisu,filename,2)
         !endif
         !OUTPOST TKE FIELD TO DISK end

         !compute turbulent intensity! (optional)
         !call horizontal_avrge(sqrt((two*tke3d2)/(ux2**2+uy2**2+uz2**2)),tkesxz(:,b)) ! turbulent intensity !! sqrt(w*tke)/|U|

         if (iscalar == 1) then
            call transpose_x_to_y(phi1(:, :, :, 1), phi2(:, :, :, 1))
            call horizontal_avrge(phi2(:, :, :, 1), phisxz(:, b)) ! phi
            call extract_fluctuat(phi2(:, :, :, 1), phisxz(:, b), phi2p) ! phi'
            !phi2p = phi2 !new idea - solve only for fluctuations
            call transpose_y_to_x(phi2p, phi1p)
            call horizontal_avrge(phi2p**2, phipphipsxz(:, b)) ! phi' phi'
            call horizontal_avrge(ux2p * phi2p, upphipsxz(:, b)) ! ux' phi'
            call horizontal_avrge(uy2p * phi2p, vpphipsxz(:, b)) ! uy' phi'
            call horizontal_avrge(uz2p * phi2p, wpphipsxz(:, b)) ! uz' phi'
            call deryS(dphidysxz(:, b), phisxz(:, b), di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, zero)!1
         end if

         call derx(ta1, ux1p, di1, sx, ffx, fsx, fwx, xsize(1), xsize(2), xsize(3), 0, ubcx)
         call derx(tb1, uy1p, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, ubcy)
         call derx(tc1, uz1p, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, ubcz)
         !y-derivatives
         call dery(ta2, ux2p, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, ubcx)
         call dery(tb2, uy2p, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0, ubcy)
         call dery(tc2, uz2p, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, ubcz)
       !!z-derivatives
         call transpose_y_to_z(ux2p, td3)
         call transpose_y_to_z(uy2p, te3)
         call transpose_y_to_z(uz2p, tf3)
         call derz(ta3, td3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, ubcx)
         call derz(tb3, te3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, ubcy)
         call derz(tc3, tf3, di3, sz, ffz, fsz, fwz, zsize(1), zsize(2), zsize(3), 0, ubcz)
       !!all back to x-pencils
         call transpose_z_to_y(ta3, td2)
         call transpose_z_to_y(tb3, te2)
         call transpose_z_to_y(tc3, tf2)
         call transpose_y_to_x(td2, tg1)
         call transpose_y_to_x(te2, th1)
         call transpose_y_to_x(tf2, ti1)
         call transpose_y_to_x(ta2, td1)
         call transpose_y_to_x(tb2, te1)
         call transpose_y_to_x(tc2, tf1)
         !du/dx=ta1 du/dy=td1 and du/dz=tg1
         !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
         !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1

         call dissipation(ta1, tb1, tc1, td1, te1, tf1, tg1, th1, ti1, diss1) !sij sij - viscous dissipation rate
         call transpose_x_to_y(diss1, temp2)
         call horizontal_avrge(temp2, epssxz(:, b)) ! diss
         ! Optional (write to disk, needs update to new way to write to disk)
         !  if (mod(itime,ioutput).eq.0) then
         !     write(filename,"('./data/dissp',I4.4)") itime/ioutput
         !     call decomp_2d_write_one(1,diss1,filename,2)
         !  endif

         if (nrank == 0) write (*, "(' Time computing basic statistics=',F15.7,'(s)')") MPI_WTIME() - tstart

         if (mod(itime, ioutput) == 0 .and. nrank == 0 .and. itime >= ioutput) then
            tstart = MPI_WTIME()
            !call outp(enssxz,    'out/02_ens') ! mean enstrophy - not necessary now
            call outp(usxz, 'out/03_ux') ! mean u
            call outp(vsxz, 'out/04_uy') ! mean v
            call outp(wsxz, 'out/05_uz') ! mean w
            !call outp(enspenspsxz,'out/06_ensp2') ! enstrophy fluct - not necessary now
            call outp(upupsxz, 'out/07_uxp2') ! u'
            call outp(vpvpsxz, 'out/08_uyp2') ! v'
            call outp(wpwpsxz, 'out/09_uzp2') ! w"
            call outp(upvpsxz, 'out/10_uxp_uyp') ! u'v'
            call outp(vpwpsxz, 'out/11_uyp_uzp') ! v'w'
            call outp(upwpsxz, 'out/12_uxp_uzp') ! u'w'
            call outp(tkesxz, 'out/13_tke') ! TKE
            !call outp(Isxz,      'out/14_tint') ! turbulent intensity - not necessary now
            call outp(-epssxz, 'out/15_diss') ! mean eps

            if (iscalar == 1) then

               call outp(phisxz, 'out/16_phi')
               call outp(phipphipsxz, 'out/17_phip2')
               call outp(upphipsxz, 'out/18_phip_uxp')
               call outp(vpphipsxz, 'out/19_phip_uyp')
               call outp(wpphipsxz, 'out/20_phip_uzp')
               call outp(dphidysxz, 'out/21_brunt_vaisala')

            end if

            if (nrank == 0) write (*, "(' Time writing statistics=',F15.7,'(s)')") MPI_WTIME() - tstart

         end if

         if (mod(itime, ioutput) == 0) then

         !!! HERE WE DO BUDGET FOR VELOCITY AND SCALAR
         !!! NOT COMPUTED AT EVERY TIME

            tstart = MPI_WTIME()

            call transpose_x_to_y(ta1, temp2) !ta1=du/dx
            call horizontal_avrge(temp2**2, dup2dxsxz)

            call transpose_x_to_y(tb1, temp2) !tb1=dv/dx
            call horizontal_avrge(temp2**2, dvp2dxsxz)

            call transpose_x_to_y(tc1, temp2) !tc1=dw/dx
            call horizontal_avrge(temp2**2, dwp2dxsxz)

            call transpose_x_to_y(td1, temp2) !td1=du/dy
            call horizontal_avrge(temp2**2, dup2dysxz)

            call transpose_x_to_y(te1, temp2) !te1=dv/dy
            call horizontal_avrge(temp2**2, dvp2dysxz)

            call transpose_x_to_y(tf1, temp2) !tf1=dw/dy
            call horizontal_avrge(temp2**2, dwp2dysxz)

            call transpose_x_to_y(tg1, temp2) !tg1=du/dz
            call horizontal_avrge(temp2**2, dup2dzsxz)

            call transpose_x_to_y(th1, temp2) !th1=dv/dz
            call horizontal_avrge(temp2**2, dvp2dzsxz)

            call transpose_x_to_y(ti1, temp2) !ti1=dw/dz
            call horizontal_avrge(temp2**2, dwp2dzsxz)

            !du/dx=ta1 du/dy=td1 and du/dz=tg1
            !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
            !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1

            call transpose_x_to_y(pre1, pre2) ! p
            call horizontal_avrge(pre2, pre_sxz)
            call extract_fluctuat(pre2, pre_sxz, pre2p) ! p'
            call horizontal_avrge(pre2p**2, ppppsxz) ! p' p'

            !pressure transport rate !d(w' p')/dy
            call horizontal_avrge(uy2p * pre2p, vpppsxz) ! p' v'
            call dery(dvpppsxz, vpppsxz, di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, ubcy)

            !turbulent transport
            call horizontal_avrge(ux2p * ux2p * uy2p, sxz1) ! u' u' v'
            call dery(dupupvpsxz, sxz1, di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, ubcx * (ubcx * ubcy)) !1*1*0=1
            call horizontal_avrge(uy2p * uz2p * uz2p, sxz1) ! v' w' w'
            call dery(dvpwpwpsxz, sxz1, di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, ubcy * (ubcz * ubcz)) !1*1*0=1
            call horizontal_avrge(uy2p * uy2p * uy2p, sxz1) ! v' v' v'
            call dery(dvpvpvpsxz, sxz1, di2, sy, ffy, fsy, fwy, ppy, 1, ysize(2), 1, 0, ubcy * (ubcy * ubcy)) !0*0*0=0

            !shear production rate
            call dery(dusxz, usxz(:, b), di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, ubcx) !1
            call dery(dwsxz, wsxz(:, b), di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, ubcz) !1

            !viscous diffusion rate !ux/y -> 1, uy/y -> 0, uz/y -> 1
            d2tkesxzdy = zero
            dtkesxzdy = zero !auxiliary field - ignore the name
            call dery(dtkesxzdy, tkesxz(:, b), di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, ubcx) !0*0=1
            if (iimplicit <= 0) then
               call deryy(d2tkesxzdy, tkesxz(:, b), di2, sy, sfyp, ssyp, swyp, 1, ysize(2), 1, 1, ubcx) !0*0=1
               if (istret /= 0) then
                  do j = 1, ysize(2)
                     d2tkesxzdy(j) = d2tkesxzdy(j) * pp2y(j) - pp4y(j) * dtkesxzdy(j)
                  end do
               end if
            else ! (semi)implicit Y diffusion
               do j = 1, ysize(2)
                  d2tkesxzdy(j) = -pp4y(j) * dtkesxzdy(j)
               end do
            end if

            if (iscalar == 1) then

               !buoyancy production
               call deryS(dphisxz, phisxz(:, b), di2, sy, ffypS, fsypS, fwypS, ppy, 1, ysize(2), 1, 1, zero) !1

               !buoyancy turbulent transport
               call horizontal_avrge(uy2p * phi2p * phi2p, sxz1) ! vi phi' phi'
               call deryS(dvpphipphipsxz, sxz1, di2, sy, ffypS, fsypS, fwypS, ppy, 1, ysize(2), 1, 1, zero) !1*1*0=1

               !buoyancy diffusion rate !ux/y -> 1, uy/y -> 0, uz/y -> 1
               d2phipphipsxz = zero
               call deryS(dphipphipsxz, phipphipsxz(:, b), di2, sy, ffypS, fsypS, fwypS, ppy, 1, ysize(2), 1, 1, zero) !0*0=1
               if (iimplicit <= 0) then
                  call deryyS(d2phipphipsxz, phipphipsxz(:, b), di2, sy, sfypS, ssypS, swypS, 1, ysize(2), 1, 1, zero) !0*0=1
                  if (istret /= 0) then
                     do j = 1, ysize(2)
                        d2phipphipsxz(j) = d2phipphipsxz(j) * pp2y(j) - dphipphipsxz(j) * pp4y(j)
                     end do
                  end if
               else ! (semi)implicit Y diffusion
                  do j = 1, ysize(2)
                     d2phipphipsxz(j) = -dphipphipsxz(j) * pp4y(j)
                  end do
               end if

               !buoyancy dissipation rate
               call derxS(dphi1pdx, phi1p, di1, sx, ffxpS, fsxpS, fwxpS, xsize(1), xsize(2), xsize(3), 1, zero)
               call transpose_x_to_y(dphi1pdx, dphi2pdx)
               call deryS(dphi2pdy, phi2p, di2, sy, ffypS, fsypS, fwypS, ppy, ysize(1), ysize(2), ysize(3), 1, zero)
               call transpose_y_to_z(phi2p, phi3p)
               call derzS(dphi3pdz, phi3p, di3, sz, ffzpS, fszpS, fwzpS, zsize(1), zsize(2), zsize(3), 1, zero)
               call transpose_z_to_y(dphi3pdz, dphi2pdz)
               call horizontal_avrge(dphi2pdx**2 + dphi2pdy**2 + dphi2pdz**2, b_eps)

            end if

            !dudy and dwdy
            call dery(dusxzdy, usxz(:, b), di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, ubcx)
            call dery(dwsxzdy, wsxz(:, b), di2, sy, ffyp, fsyp, fwyp, ppy, 1, ysize(2), 1, 1, ubcz)

            !higher order moments
            call horizontal_avrge(ux2p**three, skew_x)
            call horizontal_avrge(ux2p**four, kurt_x)
            call horizontal_avrge(uy2p**three, skew_y)
            call horizontal_avrge(uy2p**four, kurt_y)
            call horizontal_avrge(uz2p**three, skew_z)
            call horizontal_avrge(uz2p**four, kurt_z)

            if (iscalar == 1) then
               call horizontal_avrge(phi2p**three, skew_phi)
               call horizontal_avrge(phi2p**four, kurt_phi)
            end if

            if (mod(itime, ioutput) == 0 .and. nrank == 0 .and. itime >= ioutput) then

               call outpd(skew_x / (upupsxz(:, b)**1.5), 'out/34_skew_x')!mean(y^3)/mean(y^2)^1.5
               call outpd(skew_y / (vpvpsxz(:, b)**1.5), 'out/34_skew_y')
               call outpd(skew_z / (wpwpsxz(:, b)**1.5), 'out/34_skew_z')
               if (iscalar == 1) call outpd(skew_phi / (wpwpsxz(:, b)**1.5), 'out/34_skew_phi')

               call outpd(kurt_x / (upupsxz(:, b)**2), 'out/35_kurt_x')!mean(y^4)/mean(y^2)^2
               call outpd(kurt_y / (vpvpsxz(:, b)**2), 'out/35_kurt_y')
               call outpd(kurt_z / (wpwpsxz(:, b)**2), 'out/35_kurt_z')
               if (iscalar == 1) call outpd(kurt_phi / (phipphipsxz(:, b)**2), 'out/35_kurt_phi')

               call outpd(two * dup2dxsxz / dwp2dxsxz, 'out/20_K1')
               call outpd(two * dup2dxsxz / dvp2dxsxz, 'out/20_K2')
               call outpd(two * dup2dxsxz / dup2dzsxz, 'out/20_K3')
               call outpd(two * dup2dxsxz / dup2dysxz, 'out/20_K4')

               call outpd(-two * xnu * (dup2dxsxz + dup2dysxz + dup2dzsxz), 'out/31_eps_x')
               call outpd(-two * xnu * (dvp2dxsxz + dvp2dysxz + dvp2dzsxz), 'out/31_eps_y')
               call outpd(-two * xnu * (dwp2dxsxz + dwp2dysxz + dwp2dzsxz), 'out/31_eps_z')

               call outpd(dusxzdy, 'out/32_dudy')
               call outpd(dwsxzdy, 'out/32_dwdy')

               call outpd(pre_sxz, 'out/21_pre')
               call outpd(ppppsxz, 'out/22_pre2p')

               call outpd(-dvpppsxz, 'out/28_velp_grad')
               call outpd(-half * (dupupvpsxz + dvpwpwpsxz + dvpvpvpsxz), 'out/28_t_diff')
               call outpd(-(upvpsxz(:, b) * dusxz + vpwpsxz(:, b) * dwsxz), 'out/28_prod')
               call outpd(xnu * d2tkesxzdy, 'out/28_vis_diff')
               call outpd(vsxz(:, b) * dtkesxzdy, 'out/28_meanf_transp')

               if (iscalar == 1) then
                  call outpd(Ri(1) * vpphipsxz(:, b), 'out/30_b_flux')
                  call outpd(-two * Ri(1) * vpphipsxz(:, b) * dphisxz, 'out/30_b_prod')
                  call outpd(-Ri(1) * dvpphipphipsxz, 'out/30_b_ttransp')
                  call outpd((Ri(1) * xnu / sc(1)) * d2phipphipsxz, 'out/30_b_diff')
                  call outpd((-two * Ri(1) * xnu / Sc(1)) * b_eps, 'out/30_b_diss')
               end if

            end if!(mod(itime,ioutput).eq.0 .AND. nrank.eq.0 .AND. itime.ge.ioutput)

            !       !  mp=zero
            !       !  call budget(ux1,uy1,uz1,phi1,vol1)
            !       !  !call suspended(phi1,vol1,mp)
            !       !  if (nrank .eq. 0) then
            !       !     FS = 1+numscalar !Number of columns
            !       !     write(fileformat, '( "(",I4,"(E14.6),A)" )' ) FS
            !       !     FS = FS*14+1  !Line width
            !       !     open(67,file='./out/0_suspended',status='unknown',form='formatted',access='direct',recl=FS)
            !       !     write(67,fileformat,rec=itime/ioutput+1) t,mp,NL
            !       !     close(67)
            !       !     print *,'Time computing budgets (s)', real(MPI_WTIME()-tstart,4)
            !       !  endif

         end if !(mod(itime,ioutput).eq.0)
      end if !mod(itime,ilist).eq.0

      return
   end subroutine postprocess_ttbl
   !############################################################################
   !############################################################################
   subroutine outp(field, filename) ! Written by R Frantz
      !outpost to disk buffered quantities - avoid IO at every iteration
      implicit none
      real(mytype), intent(IN), dimension(ysize(2), ioutput) :: field
      character(len=*), intent(IN) :: filename
      integer :: i
      character(len=1), parameter :: NL = char(10)
      character(len=300) :: fileformat
      write (fileformat, '( "(",I4,"(E16.8),A)" )') ny
      open (67, file=trim(filename), status='unknown', form='formatted', access='direct', recl=(ny * 16 + 1))
      do i = 1, ioutput / ilist
         write (67, fileformat, rec=itime / ilist - ioutput / ilist + i) field(:, i), NL
      end do
      close (67)
   end subroutine outp
   !############################################################################
   !############################################################################
   subroutine outpd(field, filename) ! Written by R Frantz
      implicit none
      real(mytype), intent(IN), dimension(ysize(2)) :: field
      character(len=*), intent(IN) :: filename
      integer :: i
      character(len=1), parameter :: NL = char(10)
      character(len=300) :: fileformat
      write (fileformat, '( "(",I4,"(E16.8),A)" )') ny
      open (67, file=trim(filename), status='unknown', form='formatted', access='direct', recl=(ny * 16 + 1))
      write (67, fileformat, rec=itime / ioutput) field, NL
      close (67)
   end subroutine outpd
   !############################################################################
   !############################################################################
   subroutine horizontal_avrge(field2, profile) ! Written by R Frantz
      implicit none
      real(mytype), intent(in), dimension(ysize(1), ysize(2), ysize(3)) :: field2
      real(mytype), dimension(ysize(2)) :: sxz, sxz1
      real(mytype), intent(out), dimension(ysize(2)) :: profile
      integer :: i, j, k, code
      sxz = zero
      do k = 1, ysize(3)
         do j = 1, ysize(2)
            do i = 1, ysize(1)
               sxz(j) = sxz(j) + field2(i, j, k)
            end do
         end do
      end do
      ! Why not ?
      !  do j=1,ysize(2)
      !     sxz(j) = sum(field2(:,j,:))
      !  enddo
      call MPI_ALLREDUCE(sxz, sxz1, ny, real_type, MPI_SUM, MPI_COMM_WORLD, code)
      profile = sxz1 / real(nx * nz, mytype)
   end subroutine horizontal_avrge
   !############################################################################
   !############################################################################
   subroutine extract_fluctuat(field2, profile, field2p) ! Written by R Frantz
      implicit none
      real(mytype), intent(in), dimension(ysize(1), ysize(2), ysize(3)) :: field2
      real(mytype), intent(in), dimension(ysize(2)) :: profile
      real(mytype), intent(out), dimension(ysize(1), ysize(2), ysize(3)) :: field2p
      integer :: i, j, k
      do k = 1, ysize(3)
         do j = 1, ysize(2)
            do i = 1, ysize(1)
               field2p(i, j, k) = field2(i, j, k) - profile(j)
            end do
         end do
      end do
      ! why not ?
      !  do j=1,ysize(2)
      !     field2p(:,j,:) = field2(:,j,:) - profile(j)
      !  enddo
   end subroutine extract_fluctuat
   !############################################################################
   !############################################################################
   subroutine pressure_x(pp3, pre1)
      use var, only: pp1, ta1, di1, nxmsize
      use var, only: pp2, ppi2, dip2, ph2, nymsize
      use var, only: ppi3, dip3, ph3, nzmsize
      use var, only: npress
      use tools, only: rescale_pressure

      implicit none
      real(mytype), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(in) :: pp3
      real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(out) :: pre1

      call interzpv(ppi3, pp3(:, :, :, 1), dip3, sz, cifip6z, cisip6z, ciwip6z, cifz6, cisz6, ciwz6, &
                    (ph3%zen(1) - ph3%zst(1) + 1), (ph3%zen(2) - ph3%zst(2) + 1), nzmsize, zsize(3), 1)
      call transpose_z_to_y(ppi3, pp2, ph3) !nxm nym nz
      call interypv(ppi2, pp2, dip2, sy, cifip6y, cisip6y, ciwip6y, cify6, cisy6, ciwy6, &
                    (ph3%yen(1) - ph3%yst(1) + 1), nymsize, ysize(2), ysize(3), 1)
      call transpose_y_to_x(ppi2, pp1, ph2) !nxm ny nz
      call interxpv(pre1, pp1, di1, sx, cifip6, cisip6, ciwip6, cifx6, cisx6, ciwx6, &
                    nxmsize, xsize(1), xsize(2), xsize(3), 1)

      call rescale_pressure(pre1)

      return
   end subroutine pressure_x
   !############################################################################
   !############################################################################
   subroutine dissipation(ta1, tb1, tc1, td1, te1, tf1, tg1, th1, ti1, diss1)
      implicit none
      real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ta1, tb1, tc1, td1, te1, tf1, tg1, th1, ti1, diss1
      real(mytype), dimension(3, 3, xsize(1), xsize(2), xsize(3)) :: A
      integer :: k, j, i, m, l
      !INSTANTANEOUS DISSIPATION RATE
      diss1 = zero
      A(:, :, :, :, :) = zero
      A(1, 1, :, :, :) = ta1(:, :, :) !du/dx=ta1
      A(2, 1, :, :, :) = tb1(:, :, :) !dv/dx=tb1
      A(3, 1, :, :, :) = tc1(:, :, :) !dw/dx=tc1
      A(1, 2, :, :, :) = td1(:, :, :) !du/dy=td1
      A(2, 2, :, :, :) = te1(:, :, :) !dv/dy=te1
      A(3, 2, :, :, :) = tf1(:, :, :) !dw/dy=tf1
      A(1, 3, :, :, :) = tg1(:, :, :) !du/dz=tg1
      A(2, 3, :, :, :) = th1(:, :, :) !dv/dz=th1
      A(3, 3, :, :, :) = ti1(:, :, :) !dw/dz=ti1
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
      return
   end subroutine dissipation
   !############################################################################
   !############################################################################
   subroutine visu_ttbl_init(visu_initialised)

      use decomp_2d, only: mytype
      use decomp_2d_io, only: decomp_2d_register_variable
      use visu, only: io_name, output2D

      implicit none

      logical, intent(out) :: visu_initialised

      visu_initialised = .true.

   end subroutine visu_ttbl_init
   !############################################################################
   !############################################################################
   subroutine visu_ttbl(ux1, uy1, uz1, pp3, phi1, ep1, num)

      use var, only: ux2, uy2, uz2, ux3, uy3, uz3
      use var, only: ta1, tb1, tc1, td1, te1, tf1, tg1, th1, ti1, di1
      use var, only: ta2, tb2, tc2, td2, te2, tf2, di2, ta3, tb3, tc3, td3, te3, tf3, di3
      use var, only: nxmsize, nymsize, nzmsize
      use visu, only: write_field
      use ibm_param, only: ubcx, ubcy, ubcz

      implicit none

      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
      real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ep1
      integer, intent(in) :: num

   end subroutine visu_ttbl

end module ttbl
