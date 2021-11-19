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

module lockexch

  use decomp_2d, only : mytype, real_type, real2_type
  use decomp_2d, only : xsize, ysize, zsize
  use decomp_2d, only : xstart, ystart, zstart
  use decomp_2d, only : xend, yend, zend
  use decomp_2d, only : transpose_x_to_y, transpose_y_to_z, transpose_z_to_y, transpose_y_to_x

  use variables, only : numscalar

  use var, only : xnu, ri, uset, sc, Fr, prandtl
  use var, only : gravy
  use var, only : nrank
  use var, only : dens1, dens2

  use var, only : zero, half, one, two, five, twelve, thirteen
  use var, only : dx, dy, dz, nx, ny, nz

  use var, only : nrhotime

  use param, only : ilmn, ibirman_eos
  use param, only : itime, ioutput, iprocessing
  use param, only : t

  implicit none

  real(mytype), save :: pfront
  real(mytype), save, allocatable :: area2(:,:), vol1(:,:,:)

  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  logical, save :: init = .FALSE.

  character(len=*), parameter :: io_bcle = "BC-lock-exchange-io", &
       bcle_dir = "data-lock-exchange"
  
  private
  public :: init_lockexch, boundary_conditions_lockexch, postprocess_lockexch, &
       pfront, set_fluid_properties_lockexch, visu_lockexch_init

contains

  subroutine boundary_conditions_lockexch (rho1, phi1)

    USE param
    USE variables
    USE decomp_2d
    USE MPI

    implicit none

    integer  :: i,j,k,is
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1

    if (xstart(2).eq.1) then
       j = 1

       if (ncly1.eq.2) then !! Set velocity BCs
          byx1(:,:) = zero
          byy1(:,:) = zero
          byz1(:,:) = zero

          do k = 1, xsize(3)
             do i = 1, xsize(1)
                rho1(i, j, k, 1) = rho1(i, j + 1, k, 1) !! drho/dy=0
             enddo
          enddo
       endif

       if (nclyS1.eq.2) then !! Set scalar BCs
          do is=1, numscalar
             do k=1,xsize(3)
                do i=1,xsize(1)
                   !phi2(i,yend(2),k)= phi2(i,yend(2)-1,k) / (1.+uset(is)*dy*sc(is)/xnu) !Robin on top BC

                   if ((uset(is).ne.zero).and.&
                        (phi1(i,j+1,k,is).gt.phi1(i,j,k,is))) then
                      phi1(i,j,k,is) = phi1(i,j,k,is) - &
                           ((uset(is)*gdt(itr))*gravy/dy)*(phi1(i,j+1,k,is)-phi1(i,j,k,is)) !Deposit on bottom BC
                   else
                      phi1(i,j,k,is) = phi1(i,j+1,k,is)! dc/dn=0
                   endif
                enddo
             enddo
          enddo
       endif
    endif

    return
  end subroutine boundary_conditions_lockexch

  subroutine init_lockexch (rho1,ux1,uy1,uz1,ep1,phi1)

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    USE var, only : mu1

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: um,x,y,ek,ep,ekg,epg
    integer :: k,j,i,ii,is,it,code

    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             x=real(i+xstart(1)-1-1,mytype)*dx-pfront
             do is=1,numscalar
                phi1(i,j,k,is)=half * (one - tanh((sc(is) / xnu)**(half) * x)) * cp(is)
             enddo
             rho1(i,j,k,1) = half * (one - tanh((prandtl / xnu)**half * x)) &
                  * (dens1 - dens2) + dens2
          enddo
       enddo
    enddo
    do k = 1,xsize(3)
       do j = 1, xsize(2)
          do i = 1, xsize(1)
             do is=1,numscalar
                if (phi1(i,j,k,is).gt.cp(is)) then
                   phi1(i,j,k,is) = cp(is)
                elseif (phi1(i,j,k,is).lt.zero) then
                   phi1(i,j,k,is) = zero
                endif
             enddo

             if (rho1(i,j,k,1).gt.max(dens1, dens2)) then
                rho1(i,j,k,1) = max(dens1, dens2)
             elseif (rho1(i,j,k,1).lt.min(dens1, dens2)) then
                rho1(i,j,k,1) = min(dens1, dens2)
             endif
          enddo
       enddo
    enddo

    if (ilmn) then
       call set_fluid_properties_lockexch(rho1, mu1)
    end if
    
    ux1=zero; uy1=zero; uz1=zero

    if (iin.ne.0) then
       call system_clock(count=code)
       if (iin.eq.2) code=0
       call random_seed(size = ii)
       call random_seed(put = code+63946*(nrank+1)*(/ (i - 1, i = 1, ii) /))

       call random_number(ux1)
       call random_number(uy1)
       call random_number(uz1)

       !lock-exchange
       do k=1,xsize(3)
          do j=1,xsize(2)
             do i=1,xsize(1)
                x=real(i-1,mytype)*dx-pfront
                um=exp(-twentyfive*x*x)*init_noise
                ux1(i,j,k)=um*(two*ux1(i,j,k)-one)
                uy1(i,j,k)=um*(two*uy1(i,j,k)-one)
                uz1(i,j,k)=um*(two*uz1(i,j,k)-one)
             enddo
          enddo
       enddo

       ek = zero
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                ek = ek + half * rho1(i, j, k, 1) &
                     * (ux1(i, j, k)**2 + uy1(i, j, k)**2 + uz1(i, j, k)**2)
             enddo
          enddo
       enddo
       ep = zero
       do is = 1, numscalar
          do k = 1, xsize(3)
             y = (j + xstart(2) - 2) * dy
             do j = 1, xsize(2)
                do i = 1, xsize(1)
                   ep = ep - phi1(i, j, k, is) * ri(is) * (gravy * y)
                enddo
             enddo
          enddo
       enddo

       if (ilmn.and.((Fr**2).gt.zero)) then
          do k = 1, xsize(3)
             do j = 1, xsize(2)
                y = (j + xstart(2) - 2) * dy
                do i = 1, xsize(1)
                   ep = ep - (rho1(i, j, k, 1) - min(dens1, dens2)) * (gravy * y) / Fr**2
                enddo
             enddo
          enddo
       endif

       call MPI_ALLREDUCE(ek,ekg,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       call MPI_ALLREDUCE(ep,epg,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

       if ((epg.ne.zero).and.(ekg.ne.zero)) then
          um = ekg / epg
          um = init_noise / um
          um = sqrt(um)

          ux1(:,:,:) = um * ux1(:,:,:)
          uy1(:,:,:) = um * uy1(:,:,:)
          uz1(:,:,:) = um * uz1(:,:,:)

          ek = zero
          do k = 1, xsize(3)
             do j = 1, xsize(2)
                do i = 1, xsize(1)
                   ek = ek + half * rho1(i, j, k, 1) &
                        * (ux1(i, j, k)**2 + uy1(i, j, k)**2 + uz1(i, j, k)**2)
                enddo
             enddo
          enddo
          call MPI_ALLREDUCE(ek,ekg,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

          if (nrank == 0) then
             write(*,*)  "Ek / Ep: ", ekg / epg, ekg, epg
          endif
       endif
    endif

#ifdef DEBG
    if (nrank  ==  0) write(*,*) '# init end ok'
#endif

    return
  end subroutine init_lockexch

  subroutine visu_lockexch_init(visu_initialised)

    use decomp_2d, only : mytype
    use decomp_2d_io, only : decomp_2d_register_variable
    
    implicit none

    logical, intent(out) :: visu_initialised

    call decomp_2d_register_variable(io_bcle, "dissm", 3, 0, 3, mytype)
    call decomp_2d_register_variable(io_bcle, "dep", 2, 0, 2, mytype)

    visu_initialised = .true.
    
  end subroutine visu_lockexch_init

  subroutine postprocess_lockexch(rho1,ux1,uy1,uz1,phi1,ep1) !By Felipe Schuch

    use decomp_2d, only : alloc_x

    use var, only : phi2, rho2
    use var, only : phi3, rho3
    use tools, only : mean_plane_z

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: phisum1
    real(mytype),dimension(ysize(1),ysize(3),numscalar) :: dep2
    real(mytype),dimension(zsize(1),zsize(2),numscalar) :: phim3
    real(mytype),dimension(zsize(1),zsize(2)) :: rhom3
    integer :: i,j,k,is

    real(mytype) :: mp(numscalar),dms(numscalar),xf(1:2,1:3),xf2d(1:2,1:2)

    if (.not.init) then
       call alloc_x(vol1, opt_global=.true.)
       do k=xstart(3),xend(3)
          do j=xstart(2),xend(2)
             do i=xstart(1),xend(1)
                vol1(i,j,k)=dx*dy*dz
                if (i  ==  1 .or. i  ==  nx) vol1(i,j,k) = vol1(i,j,k) * half !five/twelve
                if (j  ==  1 .or. j  ==  ny) vol1(i,j,k) = vol1(i,j,k) * half !five/twelve
                if (k  ==  1 .or. k  ==  nz) vol1(i,j,k) = vol1(i,j,k) * half !five/twelve
                ! if (i  ==  2 .or. i  ==  nx-1) vol1(i,j,k) = vol1(i,j,k) * thirteen/twelve
                ! if (j  ==  2 .or. j  ==  ny-1) vol1(i,j,k) = vol1(i,j,k) * thirteen/twelve
                ! if (k  ==  2 .or. k  ==  nz-1) vol1(i,j,k) = vol1(i,j,k) * thirteen/twelve
             end do
          end do
       end do

       allocate(area2(ystart(1):yend(1),ystart(3):yend(3)))
       area2=dx*dz
       do k=ystart(3),yend(3)
          do i=ystart(1),yend(1)
             if (i  ==  1 .or. i  ==  nx) area2(i,k) = area2(i,k)/two
             if (k  ==  1 .or. k  ==  nz)  area2(i,k) = area2(i,k)/two
          end do
       end do

       init = .TRUE.
    endif

    if (mod(itime,iprocessing).ne.0) return

    mp=zero; dms=zero; xf=zero; xf2d=zero

    phisum1(:,:,:)=zero
    do is=1, numscalar
       do k=1,xsize(3)
          do j=1,xsize(2)
             do i=1,xsize(1)
                phisum1(i,j,k)=phisum1(i,j,k)+phi1(i,j,k,is)
             enddo
          end do
       end do
       call transpose_x_to_y(phi1(:,:,:,is),phi2(:,:,:,is))
       call transpose_y_to_z(phi2(:,:,:,is),phi3(:,:,:,is))
       call mean_plane_z(phi3(:,:,:,is),zsize(1),zsize(2),zsize(3),phim3(:,:,is))
    enddo

    do is = 2, numscalar
       phim3(:,:,1) = phim3(:,:,1) + phim3(:,:,is)
    enddo

    if (ilmn) then
       call transpose_x_to_y(rho1(:,:,:,1), rho2)
       call transpose_y_to_z(rho2, rho3)
       call mean_plane_z(rho3, zsize(1), zsize(2), zsize(3), rhom3)

       if (ibirman_eos) then
          phisum1(:,:,:) = phisum1(:,:,:) + (rho1(:,:,:,1) - dens2) / (dens1 - dens2)
          if (numscalar.gt.0) then
             do j = 1, zsize(2)
                do i = 1, zsize(1)
                   phim3(i,j,1) = phim3(i,j,1) + (rhom3(i,j) - dens2) / (dens1 - dens2)
                enddo
             enddo
          endif
       endif
    endif

    call budget(rho1,ux1,uy1,uz1,phi1,vol1)
    call dep(phi1,dep2)
    call suspended(phi1,vol1,mp)
    call depositrate (dep2,dms)
    call front(phisum1,xf)
    if (numscalar.gt.0) then
       call front2d(phim3(:,:,1),xf2d)
    elseif (ilmn.and.ibirman_eos) then
       call front2d(rhom3(:,:),xf2d)
    endif

    if (nrank .eq. 0) then
       FS = 1+numscalar+numscalar+3+2 !Number of columns
       write(fileformat, '( "(",I4,"(E14.6),A)" )' ) FS
       FS = FS*14+1  !Line width
       open(67,file='./out/statistics',status='unknown',form='formatted',&
            access='direct',recl=FS)
       write(67,fileformat,rec=itime/iprocessing+1) t,& !1
            mp,&                                    !numscalar
            dms,&                                   !numscalar
            xf(1,:),&                               !3
            xf2d(1,:),&                             !2
            NL                                      !+1
       close(67)
    end if

  end subroutine postprocess_lockexch

  subroutine budget(rho1,ux1,uy1,uz1,phi1,vol1)

    USE decomp_2d
    USE decomp_2d_io
    USE MPI

    use param, only : iimplicit

    use variables, only : derx, dery, derys, derz
    use variables, only : derxx, deryy, derzz
    use variables, only : derxxs, deryys, derzzs

    use var, only : ffx, ffxp, fsx, fsxp, fwx, fwxp, sfxps, ssxps, swxps, sx, &
         sfxp, ssxp, swxp
    use var, only : ffy, ffyp, ffys, fsy, fsyp, fsys, fwy, fwyp, fwys, ppy, sfyps, ssyps, swyps, sy, &
         sfyp, ssyp, swyp
    use var, only : ffz, ffzp, fsz, fszp, fwz, fwzp, sfzps, sszps, swzps, sz, &
         sfzp, sszp, swzp

    use var, only : mu1
    use var, only : rho2, ux2, uy2, uz2, phi2
    use var, only : rho3, ux3, uy3, uz3, phi3

    use var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    use var, only : ta2,tb2,tc2,td2,te2,tf2,di2
    use var, only : ta3,tb3,tc3,di3
    use ibm_param
    use param, only : zero
    use tools, only : mean_plane_z

    implicit none

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,vol1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: diss1, dphixx1
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: dphiy2, dphixx2, dphiyy2, dphizz2, ddphi2, vol2, temp2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: dphizz3, temp3

    real(mytype),dimension(3,3,xsize(1),xsize(2),xsize(3)) :: A

    real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu

    real(8) :: ek,ek1,dek,dek1,ep,ep1,dep,dep1,xvol
    integer :: ijk,i,j,k,l,m,is,code
    character(len=30) :: filename

    real(mytype) :: visc
    real(mytype) :: y

    ek=zero;ek1=zero;dek=zero;dek1=zero;ep=zero;ep1=zero;dep=zero;dep1=zero;diss1=zero

    !x-derivatives
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz)
    !y-derivatives
    call transpose_x_to_y(ux1,ux2)
    call transpose_x_to_y(uy1,uy2)
    call transpose_x_to_y(uz1,uz2)
    call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
    call dery (tb2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
    call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
    !!z-derivatives
    call transpose_y_to_z(ux2,ux3)
    call transpose_y_to_z(uy2,uy3)
    call transpose_y_to_z(uz2,uz3)
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

    A(:,:,:,:,:)=zero
    A(1,1,:,:,:)=ta1(:,:,:)
    A(2,1,:,:,:)=tb1(:,:,:)
    A(3,1,:,:,:)=tc1(:,:,:)
    A(1,2,:,:,:)=td1(:,:,:)
    A(2,2,:,:,:)=te1(:,:,:)
    A(3,2,:,:,:)=tf1(:,:,:)
    A(1,3,:,:,:)=tg1(:,:,:)
    A(2,3,:,:,:)=th1(:,:,:)
    A(3,3,:,:,:)=ti1(:,:,:)

    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             if (ilmn) then
                visc = mu1(i, j, k)
             else
                visc = one
             end if
             do m=1,3
                do l=1,3
                   diss1(i,j,k) = diss1(i,j,k) &
                        + (two * xnu * visc) &
                        * (half * (A(l,m,i,j,k) + A(m,l,i,j,k)))**2
                enddo
             enddo
          enddo
       enddo
    enddo

    do ijk=1,xsize(1)*xsize(2)*xsize(3)
       xvol=real(vol1(ijk,1,1),8)
       ek = ek + half * xvol * rho1(ijk,1,1,1) * (ux1(ijk,1,1)**2+uy1(ijk,1,1)**2+uz1(ijk,1,1)**2)
       dek = dek + xvol * diss1(ijk,1,1)
    enddo

    call transpose_x_to_y(vol1,vol2)

    ! if (ivirt==2) then
    !    ilag=0
    ! endif
    do is=1, numscalar
       if (ri(is) .eq. zero) cycle
       call derxxS (dphixx1,phi1(:,:,:,is),di1,sx,sfxpS,ssxpS,swxpS,xsize(1),xsize(2),xsize(3),1,zero)

       call transpose_x_to_y(dphixx1,dphixx2)

       call transpose_x_to_y(phi1(:,:,:,is),phi2(:,:,:,is))

       call deryS (dphiy2,phi2(:,:,:,is),di2,sy,ffyS,fsyS,fwyS,ppy,ysize(1),ysize(2),ysize(3),1,zero)

       iimplicit = -iimplicit
       call deryyS (dphiyy2,phi2(:,:,:,is),di2,sy,sfypS,ssypS,swypS,ysize(1),ysize(2),ysize(3),1,zero)
       iimplicit = -iimplicit

       call transpose_y_to_z(phi2(:,:,:,is),phi3(:,:,:,is))

       call derzzS (dphizz3,phi3(:,:,:,is),di3,sz,sfzpS,sszpS,swzpS,zsize(1),zsize(2),zsize(3),1,zero)

       call transpose_z_to_y(dphizz3,dphizz2)

       ddphi2(:,:,:)=dphixx2(:,:,:)+dphiyy2(:,:,:)+dphizz2(:,:,:)

       do k=1,ysize(3)
          do j=1,ysize(2)
             y = (j + ystart(2) - 2) * dy
             do i=1,ysize(1)
                xvol=real(vol2(i,j,k),8)
                ep = ep - xvol * ri(is) * phi2(i,j,k,is) * (gravy * y)
                dep = dep &
                     - xvol * ri(is) * (ddphi2(i,j,k)*xnu/sc(is) &
                     - uset(is) * (gravy * dphiy2(i,j,k))) &
                     * (gravy * y)
             enddo
          enddo
       enddo
    enddo

    if (ilmn.and.((Fr**2).gt.zero)) then
       call derxx(ta1, rho1(:,:,:,1), di1, sx, sfxp, ssxp, swxp, xsize(1), xsize(2), xsize(3), 1,zero)
       call transpose_x_to_y(ta1, ta2)
       call transpose_x_to_y(rho1(:,:,:,1), rho2)

       iimplicit = -iimplicit
       call deryy(tb2, rho2, di2, sy, sfyp, ssyp, swyp, ysize(1), ysize(2), ysize(3), 1,zero)
       iimplicit = -iimplicit
       call transpose_y_to_z(rho2, rho3)

       call derzz(ta3, rho3, di3, sz, sfzp, sszp, swzp, zsize(1), zsize(2), zsize(3), 1,zero)
       call transpose_z_to_y(ta3, tc2)

       do k = 1, ysize(3)
          do j = 1, ysize(2)
             y = (j + ystart(2) - 2) * dy
             do i = 1, ysize(1)
                xvol = real(vol2(i, j, k), 8)
                ep = ep - xvol * (one / Fr**2) * rho2(i, j, k) * (gravy * y)
                dep = dep - xvol * ((xnu / prandtl / (Fr**2)) &
                     * (ta2(i, j, k) + tb2(i, j, k) + tc2(i, j, k))) &
                     * (gravy * y)
             enddo
          enddo
       enddo
    endif
    ! if (ivirt==2) then
    !    ilag=1
    ! endif

    call MPI_REDUCE(ek,ek1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
    call MPI_REDUCE(dek,dek1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
    call MPI_REDUCE(ep,ep1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
    call MPI_REDUCE(dep,dep1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

    if (nrank .eq. 0) then
       open(67,file='./out/budget',status='unknown',form='formatted',&
            access='direct',recl=71) !71=5*14+1
       write(67,"(5E14.6,A)",rec=itime/iprocessing+1) t,ek1,dek1,ep1,dep1,NL
       close(67)
    end if

    if (mod(itime,ioutput).eq.0) then
       !if (save_diss.eq.1) then
       uvisu=zero
       call fine_to_coarseV(1,diss1,uvisu)
       write(filename,"('diss',I4.4)") itime/ioutput
       call decomp_2d_write_one(1,uvisu,bcle_dir,filename,2,io_bcle)
       !endif

       !if (save_dissm.eq.1) then
       call transpose_x_to_y (diss1,temp2)
       call transpose_y_to_z (temp2,temp3)
       call mean_plane_z(temp3,zsize(1),zsize(2),zsize(3),temp3(:,:,1))
       write(filename,"('dissm',I4.4)") itime/ioutput
       call decomp_2d_write_plane(3,temp3,3,1,bcle_dir,filename,io_bcle)
       !endif
    endif

  end subroutine budget

  subroutine dep(phi1,dep2)

    USE decomp_2d_io
    USE MPI

    use var, only : phi2

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),numscalar) :: phi1
    real(mytype),intent(out),dimension(ystart(1):yend(1),ystart(3):yend(3),numscalar) :: dep2

    real(mytype),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3),numscalar) :: tempdep2

    integer :: i, k, is
    character(len=30) :: filename

    dep2 = zero
    do is=1, numscalar
       if (uset(is) .eq. zero) cycle
       call transpose_x_to_y(phi1(:,:,:,is),phi2(:,:,:,is))

       tempdep2=zero

       do k=ystart(3),yend(3)
          do i=ystart(1),yend(1)
             tempdep2(i,1,k,is) = phi2(i,1,k,is)*uset(is)
             dep2(i,k,is) = tempdep2(i,1,k,is)
             !dep2(i,k,is) = dep2(i,k,is) +  phi2(i,1,k,is)*uset(is)
          end do
       end do

       write(filename,"('dep',I1.1,I4.4)") is,itime/iprocessing
       call decomp_2d_write_plane(2,tempdep2(:,:,:,is),2,1,bcle_dir,filename,io_bcle)
    enddo

  end subroutine dep

  subroutine suspended(phi1,vol1,mp1)

    USE decomp_2d_io
    USE MPI

    implicit none

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: vol1
    real(mytype),intent(out) :: mp1(1:numscalar)

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: temp1
    real(mytype) :: mp(1:numscalar)
    integer :: is,code

    mp=zero; mp1=zero

    do is=1, numscalar
       temp1 = phi1(:,:,:,is)*vol1(:,:,:)
       mp(is)= sum(temp1)
    end do

    call MPI_REDUCE(mp,mp1,numscalar,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

    return
  end subroutine suspended

  subroutine depositrate ( dep2, dms1)

    USE decomp_2d_io
    USE MPI

    real(mytype),intent(in),dimension(ystart(1):yend(1),ystart(3):yend(3),numscalar) :: dep2
    real(mytype),intent(out) :: dms1(numscalar)

    real(mytype) :: dms(numscalar)
    integer :: i,k,is,code

    dms=zero; dms1=zero
    do is=1, numscalar
       if (uset(is) .eq. zero) cycle
       do k=ystart(3),yend(3)
          do i=ystart(1),yend(1)
             dms(is)=dms(is)+dep2(i,k,is)*area2(i,k)
          end do
       end do
    enddo

    call MPI_REDUCE(dms,dms1,numscalar,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

  end subroutine depositrate

  subroutine front ( phisum1, xp )

    USE decomp_2d_io
    USE MPI

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: phisum1
    real(mytype),intent(out) :: xp(1:2,1:3)

    real(mytype) :: xp1(1:2)
    integer :: i, j ,k, code

    xp(2,:) = real(nrank,mytype)
    xp(1,:)=zero
    xp1=zero
    kloop: do k=xstart(3),xend(3)
       jloop: do j=xstart(2),xend(2)
          iloop: do i=xend(1), xstart(1), -1
             if ( phisum1(i,j,k) .ge. 0.01_mytype) then
                xp(1,1:3) = (/ real (i-1,mytype)*dx, real (j-1,mytype)*dy, real(k-1,mytype)*dz /)
                exit kloop
             end if
          end do iloop
       end do jloop
    end do kloop

    call MPI_ALLREDUCE(xp(:,1),xp1,1,real2_type,MPI_MAXLOC,MPI_COMM_WORLD,code)
    call MPI_Bcast(xp(1,:), 3,real_type, int(xp1(2)), MPI_COMM_WORLD,code)

  end subroutine front

  subroutine front2d (phim3, xp)

    USE decomp_2d_io
    USE MPI

    real(mytype),intent(in),dimension(zstart(1):zend(1),zstart(2):zend(2)) :: phim3
    real(mytype),intent(out) :: xp(1:2,1:2)
    real(mytype) :: xp1(1:2),y
    integer :: i, j, code

    xp(2,:) = real(nrank,mytype)
    xp(1,:) = zero
    xp1=zero
    jloop: do j=zstart(2),zend(2)
       y = real( j - 1, mytype )*dy
       iloop: do i=zend(1), zstart(1), -1
          if (phim3(i,j).ge.0.01_mytype) then
             xp(1,1:2) = (/ real (i-1,mytype)*dx, real (j-1,mytype)*dy /)
             exit jloop
          end if
       end do iloop
    end do jloop

    call MPI_ALLREDUCE(xp(:,1),xp1,1,real2_type,MPI_MAXLOC,MPI_COMM_WORLD,code)
    call MPI_Bcast(xp(1,:), 2,real_type, int(xp1(2)), MPI_COMM_WORLD,code)

  end subroutine front2d

  subroutine set_fluid_properties_lockexch(rho1, mu1)

    implicit none

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: mu1

    mu1(:,:,:) = rho1(:,:,:)

  endsubroutine set_fluid_properties_lockexch

end module lockexch
