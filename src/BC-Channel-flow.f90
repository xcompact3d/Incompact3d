!################################################################################
!This file is part of Xcompact3d.
!
!Xcompact3d
!Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!    Xcompact3d is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation.
!
!    Xcompact3d is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the code.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    We kindly request that you cite Xcompact3d/Incompact3d in your
!    publications and presentations. The following citations are suggested:
!
!    1-Laizet S. & Lamballais E., 2009, High-order compact schemes for
!    incompressible flows: a simple and efficient method with the quasi-spectral
!    accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
!
!    2-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence
!    problems with up to 0(10^5) computational cores, Int. J. of Numerical
!    Methods in Fluids, vol 67 (11), pp 1735-1757
!################################################################################
module channel

  use decomp_2d
  use variables
  use param

  implicit none

  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  !probes
  integer, save :: nprobes, ntimes1, ntimes2
  integer, save, allocatable, dimension(:) :: rankprobes, nxprobes, nyprobes, nzprobes

  real(mytype),save,allocatable,dimension(:) :: usum,vsum,wsum,uusum,uvsum,uwsum,vvsum,vwsum,wwsum

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_channel, boundary_conditions_channel, postprocess_channel, &
       momentum_forcing_channel, &
       geomcomplex_channel

contains
  !############################################################################
  subroutine init_channel (ux1,uy1,uz1,ep1,phi1)

    use decomp_2d
    use decomp_2d_io
    use variables
    use param
    use MPI
    use dbg_schemes, only: exp_prec

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: y,r,um,r3,x,z,h,ct
    real(mytype) :: cx0,cy0,cz0,hg,lg
    integer :: k,j,i,fh,ierror,ii,is,it,code

    integer, dimension (:), allocatable :: seed

!!! CM: Variables for "consistent" initialisation

    real(mytype), dimension(:,:,:), allocatable :: urand

    integer :: xshift, yshift, zshift

    integer ( kind = 4 ) :: seed1, seed2, seed3
    real( kind = 8 ) :: r8_random
    external r8_random

    if (iscalar==1) then
      if (nrank==0) write(*,*) 'Imposing linear temperature profile'
      do k=1,xsize(3)
         do j=1,xsize(2)
            if (istret==0) y=real(j+xstart(2)-2,mytype)*dy
            if (istret/=0) y=yp(j+xstart(2)-1)
            do i=1,xsize(1)
               phi1(i,j,k,:) = one - y/yly
            enddo
         enddo
      enddo

      phi1(:,:,:,:) = zero !change as much as you want
      if ((nclyS1 == 2).and.(xstart(2) == 1)) then
        !! Generate a hot patch on bottom boundary
        phi1(:,1,:,:) = one
      endif
      if ((nclySn == 2).and.(xend(2) == ny)) then
        phi1(:,xsize(2),:,:) = zero
      endif
    endif
!
    ux1=zero
    uy1=zero
    uz1=zero
!
    goto 111 !!! CM: Original scheme
!
    if (iin/=0) then
       call system_clock(count=code)
       if (iin.eq.2) code=0
       call random_seed(size = ii)
       call random_seed(put = code+63946*nrank*(/ (i - 1, i = 1, ii) /))

       call random_number(ux1)
       call random_number(uy1)
       call random_number(uz1)
    endif
!
 111 continue
!

!!! CM: Use a geometrical value for urand: Compute nz_global * ny_global (assuming 10,000 * 10,000 is the max)
!!! CM: and then count how many slabs we could get (split nz_global)
    allocate(urand(nx_global, ny_global, nz_global))
!
! For ux1
!
    seed1 = 12345
    seed2 = 23456
    seed3 = 34567
!
    do k = 1, nz_global
       do j = 1, ny_global
          do i = 1, nx_global
!
             urand(i, j, k) = r8_random(seed1, seed2, seed3)
!
          enddo
       enddo
    enddo
!
    if (iin == 1) then
!
       do k = xstart(3), xend(3)
          zshift = k - xstart(3) + 1
          do j = xstart(2), xend(2)
             yshift = j - xstart(2) + 1
             do i = xstart(1), xend(1)
                xshift = i - xstart(1) + 1
!
                ux1(xshift, yshift, zshift) = urand(i, j, k)
!
             enddo
          enddo
       enddo
!
    else if (iin == 3) then
!
       do k = xstart(3), xend(3)
          zshift = k - xstart(3) + 1
          do j = xstart(2), xend(2)
             yshift = j - xstart(2) + 1
             do i = xstart(1), xend(1)
                xshift = i - xstart(1) + 1
!
                uz1(xshift, yshift, zshift) = urand(i, j, k)
!
             enddo
          enddo
       enddo
!
    endif
!
! For uy1
!
    seed1 = 2 * 12345
    seed2 = 2 * 23456
    seed3 = 2 * 34567
!
    do k = 1, nz_global
       do j = 1, ny_global
          do i = 1, nx_global
!
             urand(i, j, k) = r8_random(seed1, seed2, seed3)
!
          enddo
       enddo
    enddo
!
    do k = xstart(3), xend(3)
       zshift = k - xstart(3) + 1
       do j = xstart(2), xend(2)
          yshift = j - xstart(2) + 1
          do i = xstart(1), xend(1)
             xshift = i - xstart(1) + 1
!
             uy1(xshift, yshift, zshift) = urand(i, j, k)
!
          enddo
       enddo
    enddo
!
! For uz1
!
    seed1 = 3 * 12345
    seed2 = 3 * 23456
    seed3 = 3 * 34567
 !
    do k = 1, nz_global
       do j = 1, ny_global
          do i = 1, nx_global
!
             urand(i, j, k) = r8_random(seed1, seed2, seed3)
!
          enddo
       enddo
    enddo
!
    if (iin == 1) then
!
       do k = xstart(3), xend(3)
          zshift = k - xstart(3) + 1
          do j = xstart(2), xend(2)
             yshift = j - xstart(2) + 1
             do i = xstart(1), xend(1)
                xshift = i - xstart(1) + 1
!
                uz1(xshift, yshift, zshift) = urand(i, j, k)
!
             enddo
          enddo
       enddo
!
    else if (iin == 3) then
!
       do k = xstart(3), xend(3)
          zshift = k - xstart(3) + 1
          do j = xstart(2), xend(2)
             yshift = j - xstart(2) + 1
             do i = xstart(1), xend(1)
                xshift = i - xstart(1) + 1
!
                ux1(xshift, yshift, zshift) = urand(i, j, k)
!
             enddo
          enddo
       enddo
!
    endif
!
    deallocate(urand)
!
    !modulation of the random noise + initial velocity profile
    do k=1,xsize(3)
       do j=1,xsize(2)
          if (istret==0) y=real(j+xstart(2)-1-1,mytype)*dy-yly*half
          if (istret/=0) y=yp(j+xstart(2)-1)-yly*half
          um=exp_prec(-zptwo*y*y)
          do i=1,xsize(1)
             ux1(i,j,k)=init_noise*um*(two*ux1(i,j,k)-one)+one-y*y
             uy1(i,j,k)=init_noise*um*(two*uy1(i,j,k)-one)
             uz1(i,j,k)=init_noise*um*(two*uz1(i,j,k)-one)
          enddo
       enddo
    enddo

    !INIT FOR G AND U=MEAN FLOW + NOISE
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             ux1(i,j,k)=ux1(i,j,k)+bxx1(j,k)
             uy1(i,j,k)=uy1(i,j,k)+bxy1(j,k)
             uz1(i,j,k)=uz1(i,j,k)+bxz1(j,k)
          enddo
       enddo
    enddo

#ifdef DEBG
    if (nrank .eq. 0) write(*,*) '# init end ok'
#endif

    return
  end subroutine init_channel
  !############################################################################
  !############################################################################
  subroutine boundary_conditions_channel (ux,uy,uz,phi)

    use param
    use variables
    use decomp_2d
    use tools, only : channel_cfr_x

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
!!$  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ut

    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: gx
    real(mytype) :: x, y, z
    integer :: i, j, k, is

    if (icpg /= 1) then ! if not constant pressure gradient
      if (icfr == 1) then ! constant flow rate without transposition
        call channel_cfr_x(ux,two/three)
      else if (icfr == 2) then
        call transpose_x_to_y(ux,gx)
        call channel_flrt(gx,two/three)
        call transpose_y_to_x(gx,ux)
      end if
    end if

    if (iscalar /= 0) then
       if (iimplicit <= 0) then
          if ((nclyS1 == 2).and.(xstart(2) == 1)) then
             !! Generate a hot patch on bottom boundary
             phi(:,1,:,:) = one
          endif
          if ((nclySn == 2).and.(xend(2) == ny)) THEN
             phi(:,xsize(2),:,:) = zero
          endif
       else
          !
          ! Implicit boundary conditions are usually given in input file
          ! It is possible to modify g_sc here
          ! It is not possible to modify alpha_sc and beta_sc here
          !
          ! Bottom temperature if alpha_sc(:,1)=1 and beta_sc(:,1)=0 (default)
          !if (nclyS1.eq.2) g_sc(:,1) = one
          ! Top temperature if alpha_sc(:,2)=1 and beta_sc(:,2)=0 (default)
          !if (nclySn.eq.2) g_sc(:,2) = zero
       endif
    endif

    return
  end subroutine boundary_conditions_channel
  !############################################################################
  !############################################################################
  subroutine channel_flrt (ux,constant)

    use decomp_2d
    use decomp_2d_poisson
    use variables
    use param
    use var
    use MPI

    implicit none

    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux
    real(mytype) :: constant

    integer :: j,i,k,code
    real(mytype) :: can,ut3,ut,ut4

    ut3=zero
    do k=1,ysize(3)
       do i=1,ysize(1)
          ut=zero
          do j=1,ny-1
             if (istret == 0) then
                ut=ut+dy*(ux(i,j+1,k)-half*(ux(i,j+1,k)-ux(i,j,k)))
             else
                ut=ut+(yp(j+1)-yp(j))*(ux(i,j+1,k)-half*(ux(i,j+1,k)-ux(i,j,k)))
             endif
          enddo
          ut=ut/yly
          ut3=ut3+ut
       enddo
    enddo
    ut3=ut3/(real(nx*nz,mytype))

    call MPI_ALLREDUCE(ut3,ut4,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

    can=-(constant-ut4)

    if (nrank==0) write(*,*),nrank,'correction to ensure constant flow rate',ut4,can

    do k=1,ysize(3)
       do i=1,ysize(1)
          do j=2,ny-1
             ux(i,j,k)=ux(i,j,k)-can
          enddo
       enddo
    enddo

    return
  end subroutine channel_flrt
  !############################################################################
  !############################################################################
  subroutine postprocess_channel(ux1,uy1,uz1,pp3,phi1,ep1) !By Felipe Schuch

    use MPI
    use decomp_2d
    use decomp_2d_io
    ! SR old include files (to fix)
    !use var, only : umean,vmean,wmean,pmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
    !use var, only : phimean, phiphimean
    !use var, only : ta1, pp1, di1
    !use var, only : ppi3, dip3
    !use var, only : pp2, ppi2, dip2
    ! SR New Include (to fix)
    use var, only : uvisu
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    

    use var, ONLY : nxmsize, nymsize, nzmsize
    use param, ONLY : npress

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(in) :: pp3
    character(len=30) :: filename

    integer :: is

    if ((ivisu.ne.0).and.(mod(itime, ioutput).eq.0)) then
       !! Write of Q criterion as a postprecessing example
       !x-derivatives
       call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
       call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
       call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
       !y-derivatives
       call transpose_x_to_y(ux1,td2)
       call transpose_x_to_y(uy1,te2)
       call transpose_x_to_y(uz1,tf2)
       call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
       call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
       call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
       !!z-derivatives
       call transpose_y_to_z(td2,td3)
       call transpose_y_to_z(te2,te3)
       call transpose_y_to_z(tf2,tf3)
       call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
       call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
       call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
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
       !Q=-0.5*(ta1**2+te1**2+ti1**2)-td1*tb1-tg1*tc1-th1*tf1
       di1=zero
       di1(:,:,:)=-0.5*(ta1(:,:,:)**2+te1(:,:,:)**2+ti1(:,:,:)**2)-&
                        td1(:,:,:)*tb1(:,:,:)-&
                        tg1(:,:,:)*tc1(:,:,:)-&
                        th1(:,:,:)*tf1(:,:,:)
       uvisu=zero
       call fine_to_coarseV(1,di1,uvisu)
994    format('./data/critq',I5.5)
       write(filename, 994) itime/ioutput
       call decomp_2d_write_one(1,uvisu,filename,2)
    endif

    return
  end subroutine postprocess_channel
  !############################################################################
  !############################################################################
  !!
  !!  SUBROUTINE: momentum_forcing
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Applies rotation for t < spinup_time.
  !!
  !############################################################################
  subroutine momentum_forcing_channel(dux1, duy1, ux1, uy1)

    implicit none

    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1

    if (icpg == 1) then
        !! fcpg: add constant pressure gradient in streamwise direction
        dux1(:,:,:,1) = dux1(:,:,:,1) + fcpg !* (re/re_cent)**2
    endif

    if (itime < spinup_time) then
       if (nrank==0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime == ilast)) &
          write(*,*) 'Rotating turbulent channel at speed ',wrotation
       dux1(:,:,:,1) = dux1(:,:,:,1) - wrotation*uy1(:,:,:)
       duy1(:,:,:,1) = duy1(:,:,:,1) + wrotation*ux1(:,:,:)
    endif

  end subroutine momentum_forcing_channel
  !############################################################################
  !############################################################################
  subroutine geomcomplex_channel(epsi,nxi,nxf,ny,nyi,nyf,nzi,nzf,yp,remp)

    use decomp_2d, only : mytype
    use param, only : zero, one, two, ten
    use ibm

    implicit none

    integer                    :: nxi,nxf,ny,nyi,nyf,nzi,nzf
    real(mytype),dimension(nxi:nxf,nyi:nyf,nzi:nzf) :: epsi
    real(mytype),dimension(ny) :: yp
    real(mytype)               :: remp
    integer                    :: j
    real(mytype)               :: ym
    real(mytype)               :: zeromach
    real(mytype)               :: h

    epsi(:,:,:) = zero
    h = (yly - two) / two

    zeromach=one
    do while ((one + zeromach / two) .gt. one)
       zeromach = zeromach/two
    end do
    zeromach = ten*zeromach

    do j=nyi,nyf
       ym=yp(j)
       if ((ym.le.h).or.(ym.ge.(h+two))) then
          epsi(:,j,:)=remp
       endif
    enddo

    return
  end subroutine geomcomplex_channel
  !############################################################################
end module channel
