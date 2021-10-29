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

module tgv

  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  real(mytype), save, allocatable, dimension(:,:,:) :: vol1,volSimps1
  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_tgv, boundary_conditions_tgv, postprocess_tgv, visu_tgv, visu_tgv_init

contains

  subroutine init_tgv (ux1,uy1,uz1,ep1,phi1)

    use decomp_2d
    use decomp_2d_io
    use variables
    use param
    use MPI
    use dbg_schemes, only: sin_prec, cos_prec

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: y,r,um,r3,x,z,h,ct
    real(mytype) :: cx0,cy0,cz0,hg,lg
    integer :: k,j,i,ierror,is,code
    integer, dimension (:), allocatable :: seed
    integer ::  isize

    if (iscalar==1) then

       phi1(:,:,:,:) = zero
    endif

    if (iin.eq.0) then !empty domain

       if (nrank==0) write(*,*) "Empty initial domain!"

       ux1=zero; uy1=zero; uz1=zero

    endif

    if (iin.eq.1) then !generation of a random noise

       if (nrank==0) write(*,*) "Filled initial domain!"

       ux1=zero; uy1=zero; uz1=zero

       do k=1,xsize(3)
          z=real((k+xstart(3)-1-1),mytype)*dz
          do j=1,xsize(2)
             y=real((j+xstart(2)-1-1),mytype)*dy
             do i=1,xsize(1)
                x=real(i-1,mytype)*dx

                if (.not.tgv_twod) then
                   ux1(i,j,k)=+sin_prec(x)*cos_prec(y)*cos_prec(z)
                   uy1(i,j,k)=-cos_prec(x)*sin_prec(y)*cos_prec(z)
                   if (iscalar == 1) then
                      phi1(i,j,k,1:numscalar)=sin_prec(x)*sin_prec(y)*cos_prec(z)
                   endif
                else
                   ux1(i,j,k)=+sin_prec(x)*cos_prec(y)
                   uy1(i,j,k)=-cos_prec(x)*sin_prec(y)
                   if (iscalar == 1) then
                      phi1(i,j,k,1:numscalar)=sin_prec(x)*sin_prec(y)
                   endif
                endif
                uz1(i,j,k)=zero
             enddo
          enddo
       enddo

       call random_seed(size=isize)
       allocate (seed(isize))
       seed(:)=67
       call random_seed(put=seed)
       !     call random_number(ux1)
       !     call random_number(uy1)
       ! call random_number(uz1)

       do k=1,xsize(3)
          do j=1,xsize(2)
             do i=1,xsize(1)
                !              ux1(i,j,k)=noise*(ux1(i,j,k)-half)
                !              uy1(i,j,k)=noise*(uy1(i,j,k)-half)
                ! uz1(i,j,k)=0.05*(uz1(i,j,k)-half)
             enddo
          enddo
       enddo

       !     !modulation of the random noise
       !     do k=1,xsize(3)
       !        do j=1,xsize(2)
       !           if (istret.eq.0) y=(j+xstart(2)-1-1)*dy-yly/two
       !           if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/two
       !           um=exp(-0.2*y*y)
       !           do i=1,xsize(1)
       !              ux1(i,j,k)=um*ux1(i,j,k)
       !              uy1(i,j,k)=um*uy1(i,j,k)
       !              uz1(i,j,k)=um*uz1(i,j,k)
       !           enddo
       !        enddo
       !     enddo

    endif

    !  bxx1(j,k)=zero
    !  bxy1(j,k)=zero
    !  bxz1(j,k)=zero

    !INIT FOR G AND U=MEAN FLOW + NOISE
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             ux1(i,j,k)=ux1(i,j,k)!+bxx1(j,k)
             uy1(i,j,k)=uy1(i,j,k)!+bxy1(j,k)
             uz1(i,j,k)=uz1(i,j,k)!+bxz1(j,k)
          enddo
       enddo
    enddo

#ifdef DEBG
    if (nrank  ==  0) write(*,*) '# init end ok'
#endif

    return
  end subroutine init_tgv
  !********************************************************************

  subroutine boundary_conditions_tgv (ux,uy,uz,phi)

    USE param
    USE variables
    USE decomp_2d

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    if (nclx1 == 2) then
    endif
    if (nclxn == 2) then
    endif

    if (ncly1 == 2) then
    endif
    if (nclyn == 2) then
    endif

    if (nclz1 == 2) then
    endif
    if (nclzn == 2) then
    endif

  end subroutine boundary_conditions_tgv

  !********************************************************************

  !############################################################################
  subroutine init_post(ep1)

    USE MPI

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ep1
    real(mytype) :: dxdydz
    integer :: i,j,k,code
    character :: a

    call alloc_x(vol1, opt_global=.true.)
    call alloc_x(volSimps1, opt_global=.true.)

    vol1 = zero
    volSimps1 = zero

    !X PENCILS !Utilizar para integral volumétrica dentro do domínio físico
    dxdydz=dx*dy*dz
    do k=xstart(3),xend(3)
       do j=xstart(2),xend(2)
          do i=xstart(1),xend(1)
             vol1(i,j,k)=dxdydz
             if (i .eq. 1   .or. i .eq. nx) vol1(i,j,k) = vol1(i,j,k)/two
             if (j .eq. 1   .or. j .eq. ny)  vol1(i,j,k) = vol1(i,j,k)/two
             if (k .eq. 1   .or. k .eq. nz)  vol1(i,j,k) = vol1(i,j,k)/two
          end do
       end do
    end do

    !X PENCILS !Utilizar para integral volumétrica dentro do domínio físico (método de Simpson)
    do k=xstart(3),xend(3)
       do j=xstart(2),xend(2)
          do i=xstart(1),xend(1)
             volSimps1(i,j,k)=dxdydz
             if (i .eq. 1   .or. i .eq. nx) volSimps1(i,j,k) = volSimps1(i,j,k) * (five/twelve)
             if (j .eq. 1   .or. j .eq. ny) volSimps1(i,j,k) = volSimps1(i,j,k) * (five/twelve)
             if (k .eq. 1   .or. k .eq. nz) volSimps1(i,j,k) = volSimps1(i,j,k) * (five/twelve)
             if (i .eq. 2   .or. i .eq. nx-1) volSimps1(i,j,k) = volSimps1(i,j,k) * (thirteen/twelve)
             if (j .eq. 2   .or. j .eq. ny-1) volSimps1(i,j,k) = volSimps1(i,j,k) * (thirteen/twelve)
             if (k .eq. 2   .or. k .eq. nz-1) volSimps1(i,j,k) = volSimps1(i,j,k) * (thirteen/twelve)
          end do
       end do
    end do

#ifdef DEBG
    if (nrank == 0) write(*,*) '# init_post ok'
#endif

  end subroutine init_post
  !############################################################################
  subroutine postprocess_tgv(ux1,uy1,uz1,phi1,ep1)

    USE decomp_2d
    USE decomp_2d_io
    USE variables, only: nx,ny,nz
    USE MPI
    USE var, only: nut1, srt_smag 
    USE var, only : uvisu, ux2, uy2, uz2, ux3, uy3, uz3
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    USE ibm_param

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1,ep1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype) :: mp(numscalar),mps(numscalar),vl,es,es1,ek,ek1,ds,ds1
    real(mytype) :: temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9

    real(mytype) :: eek, enst, eps, eps2
    integer :: nxc, nyc, nzc, xsize1, xsize2, xsize3

    integer :: i,j,k,is,code,nvect1
    character(len=30) :: filename

    if (nclx1==1.and.xend(1)==nx) then
       xsize1=xsize(1)-1
    else
       xsize1=xsize(1)
    endif
    if (ncly1==1.and.xend(2)==ny) then
       xsize2=xsize(2)-1
    else
       xsize2=xsize(2)
    endif
    if (nclz1==1.and.xend(3)==nz) then
       xsize3=xsize(3)-1
    else
       xsize3=xsize(3)
    endif
    if (nclx1==1) then
       nxc=nxm
    else
       nxc=nx
    endif
    if (ncly1==1) then
       nyc=nym
    else
       nyc=ny
    endif
    if (nclz1==1) then
       nzc=nzm
    else
       nzc=nz
    endif


    !we only collect statistics every 10 time steps to save computational time
    if (mod(itime, 10) == 0) then

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

       !! Write vorticity as an example of post processing
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

       !SPATIALLY-AVERAGED ENSTROPHY
       temp1=zero
       do k=1,xsize3
          do j=1,xsize2
             do i=1,xsize1
                temp1=temp1+zpfive*((tf1(i,j,k)-th1(i,j,k))**2+&
                                    (tg1(i,j,k)-tc1(i,j,k))**2+&
                                    (tb1(i,j,k)-td1(i,j,k))**2)
             enddo
          enddo
       enddo
       call MPI_ALLREDUCE(temp1,enst,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       enst=enst/(nxc*nyc*nzc)
       
       !SPATIALLY-AVERAGED ENERGY DISSIPATION
       temp1=zero
       do k=1,xsize3
          do j=1,xsize2
             do i=1,xsize1
                temp1=temp1+half*xnu*((two*ta1(i,j,k))**two+(two*te1(i,j,k))**two+&
                                      (two*ti1(i,j,k))**two+two*(td1(i,j,k)+tb1(i,j,k))**two+&
                                                            two*(tg1(i,j,k)+tc1(i,j,k))**two+&
                                                            two*(th1(i,j,k)+tf1(i,j,k))**two)
                if (ilesmod /= 0 .and. jles <= 3) then
                   temp1=temp1+two*nut1(i,j,k)*srt_smag(i,j,k)
                endif
             enddo
          enddo
       enddo
       call MPI_ALLREDUCE(temp1,eps,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       eps=eps/(nxc*nyc*nzc)

       !SPATIALLY-AVERAGED TKE of velocity fields
       temp1=zero
       do k=1,xsize3
          do j=1,xsize2
             do i=1,xsize1
                temp1=temp1+zpfive*((ux1(i,j,k))**2+(uy1(i,j,k))**2+(uz1(i,j,k))**2)
             enddo
          enddo
       enddo
       call MPI_ALLREDUCE(temp1,eek,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       eek=eek/(nxc*nyc*nzc)

       !SECOND DERIVATIVES
       !x-derivatives
       call derxx (ta1,ux1,di1,sx,sfx ,ssx ,swx ,xsize(1),xsize(2),xsize(3),0,ubcx)
       call derxx (tb1,uy1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1,ubcy)
       call derxx (tc1,uz1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1,ubcz)
       !y-derivatives
       call deryy (ta2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1,ubcx) 
       call deryy (tb2,uy2,di2,sy,sfy ,ssy ,swy ,ysize(1),ysize(2),ysize(3),0,ubcy) 
       call deryy (tc2,uz2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1,ubcz) 
       !!z-derivatives
       call derzz (ta3,ux3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1,ubcx)
       call derzz (tb3,uy3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1,ubcy)
       call derzz (tc3,uz3,di3,sz,sfz ,ssz ,swz ,zsize(1),zsize(2),zsize(3),0,ubcz)
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
       !d2u/dx2=ta1 d2u/dy2=td1 and d2u/dz2=tg1
       !d2v/dx2=tb1 d2v/dy2=te1 and d2v/dz2=th1
       !d2w/dx2=tc1 d2w/dy2=tf1 and d2w/dz2=ti1
       !SPATIALLY-AVERAGED ENERGY DISSIPATION WITH SECOND DERIVATIVES
       temp1=zero
       di1  =zero
       do k=1,xsize3
          do j=1,xsize2
             do i=1,xsize1
                di1(i,j,k)=(-xnu)*( ux1(i,j,k)*(ta1(i,j,k)+td1(i,j,k)+tg1(i,j,k))+ &
                                    uy1(i,j,k)*(tb1(i,j,k)+te1(i,j,k)+th1(i,j,k))+ &
                                    uz1(i,j,k)*(tc1(i,j,k)+tf1(i,j,k)+ti1(i,j,k)) )
                temp1=temp1+di1(i,j,k)
             enddo
          enddo
       enddo
       call MPI_ALLREDUCE(temp1,eps2,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       eps2=eps2/(nxc*nyc*nzc)
       
       
       if (nrank==0) then
          write(42,'(20e20.12)') (itime-1)*dt,eek,eps,eps2,enst
          call flush(42)
       endif
    endif

  end subroutine postprocess_tgv

  !############################################################################
  !!
  !!  SUBROUTINE: visu_tgv_init
  !!      AUTHOR: PB
  !! DESCRIPTION: Initialises TGV-specific visualisation
  !!
  !############################################################################
  subroutine visu_tgv_init (visu_initialised)

    use decomp_2d, only : mytype
    use decomp_2d_io, only : decomp_2d_register_variable
    use visu, only : io_name, output2D
    
    implicit none

    logical, intent(out) :: visu_initialised

    call decomp_2d_register_variable(io_name, "vort", 1, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "critq", 1, 0, output2D, mytype)

    visu_initialised = .true.
    
  end subroutine visu_tgv_init
  !############################################################################
  !!
  !!  SUBROUTINE: visu_tgv
  !!      AUTHOR: FS
  !! DESCRIPTION: Performs TGV-specific visualization
  !!
  !############################################################################
  subroutine visu_tgv(ux1, uy1, uz1, pp3, phi1, ep1, num)

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
    character(len=32), intent(in) :: num

    ! Write vorticity as an example of post processing

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
    call write_field(di1, ".", "vort", trim(num), flush=.true.) ! Reusing temporary array, force flush

    !Q=-0.5*(ta1**2+te1**2+ti1**2)-td1*tb1-tg1*tc1-th1*tf1
    di1 = zero
    di1(:,:,: ) = - 0.5*(ta1(:,:,:)**2+te1(:,:,:)**2+ti1(:,:,:)**2) &
                  - td1(:,:,:)*tb1(:,:,:) &
                  - tg1(:,:,:)*tc1(:,:,:) &
                  - th1(:,:,:)*tf1(:,:,:)
    call write_field(di1, ".", "critq", trim(num), flush=.true.) ! Reusing temporary array, force flush

  end subroutine visu_tgv
  !############################################################################
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
  !############################################################################
  subroutine dissipation (ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,diss1)

    USE param
    USE variables
    USE decomp_2d
    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,diss1
    real(mytype),dimension(3,3,xsize(1),xsize(2),xsize(3)) :: A
    integer :: k,j,i,m,l

    !INSTANTANEOUS DISSIPATION RATE
    diss1=zero
    A(:,:,:,:,:)=zero
    A(1,1,:,:,:)=ta1(:,:,:) !du/dx=ta1
    A(2,1,:,:,:)=tb1(:,:,:) !dv/dx=tb1
    A(3,1,:,:,:)=tc1(:,:,:) !dw/dx=tc1
    A(1,2,:,:,:)=td1(:,:,:) !du/dy=td1
    A(2,2,:,:,:)=te1(:,:,:) !dv/dy=te1
    A(3,2,:,:,:)=tf1(:,:,:) !dw/dy=tf1
    A(1,3,:,:,:)=tg1(:,:,:) !du/dz=tg1
    A(2,3,:,:,:)=th1(:,:,:) !dv/dz=th1
    A(3,3,:,:,:)=ti1(:,:,:) !dw/dz=ti1

    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             do m=1,3
                do l=1,3
                   diss1(i,j,k)=diss1(i,j,k)+two*xnu*half*half*(A(l,m,i,j,k)+A(m,l,i,j,k))**two
                enddo
             enddo
          enddo
       enddo
    enddo

    return

  end subroutine dissipation
  !############################################################################
  subroutine error_tgv2D(ux1,uy1,phi1)

    use decomp_2d
    use MPI
    use param, only : one, two, xnu, ifirst, itime
    use variables, only : numscalar, sc
    use dbg_schemes, only: sin_prec, cos_prec

    implicit none

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    ! Analytical error
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: errx, erry
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: errs
    ! Time-discrete error
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: errxd1, erryd1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: errsd1
    ! Space-time discrete error
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: errxd2, erryd2
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: errsd2
    ! Temporary storage array
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tmp

    integer :: i,j,k,l,it
    real(mytype) :: x, y, xl1, xl2, xlinf, yl1, yl2, ylinf
    real(mytype) :: solx0, soly0, sols0
    real(mytype) :: solxt, solyt, solst, solxd, solyd, solsd, solxdd, solydd, solsdd
    real(mytype) :: xdamping(3), ydamping(3), sdamping(3,numscalar)

    ! Compute the errors
    call compute_tgv2D_errors(xdamping, ydamping, sdamping)

    ! Compute solutions and errors
    do k = 1,xsize(3)
       do j = 1,xsize(2)
          y = real(j+xstart(2)-1-1,mytype)*dy
          do i = 1,xsize(1)
            x = real(i+xstart(1)-1-1,mytype)*dx
            ! Initial solution
            solx0 = sin_prec(x) * cos_prec(y)
            soly0 = - cos_prec(x) * sin_prec(y)
            ! Analytical solution
            solxt = solx0 * xdamping(1)
            solyt = soly0 * ydamping(1)
            ! Time discrete solution
            solxd = solx0 * xdamping(2)
            solyd = soly0 * ydamping(2)
            ! Space-time discrete solution
            solxdd = solx0 * xdamping(3)
            solydd = soly0 * ydamping(3)
            ! Errors
            errx(i,j,k) = ux1(i,j,k) - solxt
            erry(i,j,k) = uy1(i,j,k) - solyt
            errxd1(i,j,k) = ux1(i,j,k) - solxd
            erryd1(i,j,k) = uy1(i,j,k) - solyd
            errxd2(i,j,k) = ux1(i,j,k) - solxdd
            erryd2(i,j,k) = uy1(i,j,k) - solydd
          enddo
       enddo
    enddo
    if (iscalar==1) then
       do l = 1, numscalar
          do k = 1,xsize(3)
             do j = 1,xsize(2)
                y = real(j+xstart(2)-1-1,mytype)*dy
                do i = 1,xsize(1)
                  x = real(i+xstart(1)-1-1,mytype)*dx
                  ! Initial solution
                  sols0 = sin_prec(x) * sin_prec(y)
                  ! Analytical solution
                  solst = sols0 * sdamping(1,l)
                  ! Time discrete solution
                  solsd = sols0 * sdamping(2,l)
                  ! Space-time discrete solution
                  solsdd = sols0 * sdamping(3,l)
                  ! Errors
                  errs(i,j,k,l) = phi1(i,j,k,l) - solst
                  errsd1(i,j,k,l) = phi1(i,j,k,l) - solsd
                  errsd2(i,j,k,l) = phi1(i,j,k,l) - solsdd
                enddo
             enddo
          enddo
       enddo
    endif

    !
    ! Compute and print synthetic errors indicator
    !
    ! Magnitude of the solution
    call error_L1_L2_Linf(ux1, xl1, xl2, xlinf)
    call error_L1_L2_Linf(uy1, yl1, yl2, ylinf)
    if (nrank==0) then
      write(*,*) "2D Taylor-Green errors:"
      write(*,*) "  Solution amplitude, L1, sqrt(L2), Linf:"
      write(*,*) "            ", xl1, xl2, xlinf
      write(*,*) "            ", yl1, yl2, ylinf
    endif
    if (iscalar==1) then
    do l = 1, numscalar
      tmp(:,:,:) = phi1(:,:,:,l)
      call error_L1_L2_Linf(tmp, xl1, xl2, xlinf)
      if (nrank==0) then
        write(*,*) "            ", xl1, xl2, xlinf
      endif
    enddo
    endif
    ! Analytical solution
    call error_L1_L2_Linf(errx, xl1, xl2, xlinf)
    call error_L1_L2_Linf(erry, yl1, yl2, ylinf)
    if (nrank==0) then
      write(*,*) "  1: Analytical, L1, sqrt(L2), Linf:"
      write(*,*) "     XERROR ", xl1, xl2, xlinf
      write(*,*) "     YERROR ", yl1, yl2, ylinf
    endif
    if (iscalar==1) then
    do l = 1, numscalar
      tmp(:,:,:) = errs(:,:,:,l)
      call error_L1_L2_Linf(tmp, xl1, xl2, xlinf)
      if (nrank==0) then
        write(*,*) "            ", xl1, xl2, xlinf
      endif
    enddo
    endif
    ! Time discrete solution
    call error_L1_L2_Linf(errxd1, xl1, xl2, xlinf)
    call error_L1_L2_Linf(erryd1, yl1, yl2, ylinf)
    if (nrank==0) then
      write(*,*) "  2: Time-discrete, L1, sqrt(L2), Linf:"                                 
      write(*,*) "     XERROR ", xl1, xl2, xlinf
      write(*,*) "     YERROR ", yl1, yl2, ylinf
    endif
    if (iscalar==1) then
    do l = 1, numscalar
      tmp(:,:,:) = errsd1(:,:,:,l)
      call error_L1_L2_Linf(tmp, xl1, xl2, xlinf)
      if (nrank==0) then
        write(*,*) "            ", xl1, xl2, xlinf
      endif
    enddo
    endif
    ! Space-time discrete solution
    call error_L1_L2_Linf(errxd2, xl1, xl2, xlinf)
    call error_L1_L2_Linf(erryd2, yl1, yl2, ylinf)
    if (nrank==0) then
      write(*,*) "  3: Time-space-discrete, L1, sqrt(L2), Linf:"
      write(*,*) "     XERROR ", xl1, xl2, xlinf
      write(*,*) "     YERROR ", yl1, yl2, ylinf
    endif
    if (iscalar==1) then
    do l = 1, numscalar
      tmp(:,:,:) = errsd2(:,:,:,l)
      call error_L1_L2_Linf(tmp, xl1, xl2, xlinf)
      if (nrank==0) then
        write(*,*) "            ", xl1, xl2, xlinf
      endif
    enddo
    endif

  end subroutine error_tgv2D

  ! Compute the damping factors
  subroutine compute_tgv2D_errors(xdamping, ydamping, sdamping)

    use decomp_2d
    use param, only : one, two, xnu, ifirst, itime, itimescheme, iimplicit
    use variables, only : numscalar, sc
    use dbg_schemes, only: exp_prec

    implicit none

    real(mytype), intent(out) :: xdamping(3), ydamping(3), sdamping(3, numscalar)

    integer :: it, l
    real(mytype) :: ktgv, k2tgv, coef(numscalar+1)

    ! Compute modified wavenumber
    ktgv = one
    call compute_k2(ktgv, k2tgv)

    ! Init
    xdamping(:) = one
    ydamping(:) = one
    sdamping(:,:) = one

    ! Compute analytical damping
    coef(1) = exp_prec(-two*dt*xnu)
    do it = ifirst, itime
       xdamping(1) = xdamping(1) * coef(1)
       ydamping(1) = ydamping(1) * coef(1)
       do l = 1, numscalar
         sdamping(1,l) = sdamping(1,l) * exp_prec(-two*dt*xnu/sc(l))
       enddo
    enddo

    ! Compute time and space-time discrete damping
    !
    ! Explicit Euler
    if (itimescheme.eq.1) then

      ! Time discrete errors
      if (iimplicit.eq.0) then
        coef(1) = (one - two*dt*xnu)
        do l = 1, numscalar
          coef(1+l) = one - two*dt*xnu/sc(l)
        enddo
      else if (iimplicit.eq.1) then
        coef(1) = (one - dt*xnu) / (one + dt*xnu)
        do l = 1, numscalar
          coef(1+l) = (one - dt*xnu/sc(l)) / (one + dt*xnu/sc(l))
        enddo
      else if (iimplicit.eq.2) then
        coef(1) = (one - onepfive*dt*xnu) / (one + half*dt*xnu)
        do l = 1, numscalar
          coef(1+l) = (one - onepfive*dt*xnu/sc(l)) / (one + half*dt*xnu/sc(l))
        enddo
      endif
      do it = ifirst, itime
        xdamping(2) = xdamping(2) * coef(1)
        ydamping(2) = ydamping(2) * coef(1)
        do l = 1, numscalar
          sdamping(2,l) = sdamping(2,l) * coef(1+l)
        enddo
      enddo

      ! Space-time discrete errors
      if (iimplicit.eq.0) then
        coef(1) = (one - two*k2tgv*dt*xnu)
        do l = 1, numscalar
          coef(1+l) = one - two*k2tgv*dt*xnu/sc(l)
        enddo
      else if (iimplicit.eq.1) then
        coef(1) = (one - k2tgv*dt*xnu) / (one + k2tgv*dt*xnu)
        do l = 1, numscalar
          coef(1+l) = (one - k2tgv*dt*xnu/sc(l)) / (one + k2tgv*dt*xnu/sc(l))
        enddo
      else if (iimplicit.eq.2) then
        coef(1) = (one - onepfive*k2tgv*dt*xnu) / (one + half*k2tgv*dt*xnu)
        do l = 1, numscalar
          coef(1+l) = (one - onepfive*k2tgv*dt*xnu/sc(l)) / (one + half*k2tgv*dt*xnu/sc(l))
        enddo
      endif
      do it = ifirst, itime
        xdamping(3) = xdamping(3) * coef(1)
        ydamping(3) = ydamping(3) * coef(1)
        do l = 1, numscalar
          sdamping(3,l) = sdamping(3,l) * coef(1+l)
        enddo
      enddo

    else

      if (nrank==0) write(*,*)  "TGV2D: No discrete error implemented for this time scheme."
      xdamping(2:) = xdamping(1)
      ydamping(2:) = ydamping(1)
      do l = 1, numscalar
        sdamping(2:,l) = sdamping(1,l)
      enddo

    endif

  end subroutine compute_tgv2D_errors

  ! Compute the modified wavenumber for the second derivative
  ! Warning : we use the X momentum wavenumber for Y momentum and for the scalars
  subroutine compute_k2(kin, k2out)

    use decomp_2d, only : mytype
    use param
    use derivx, only : alsaix, asix, bsix, csix, dsix
    use dbg_schemes, only: cos_prec

    implicit none

    real(mytype), intent(in) :: kin
    real(mytype), intent(out) :: k2out

    if (kin.lt.zero .or. kin.gt.pi/min(dx,dy)) then
      if (nrank==0) then
        write(*,*) "TGV2D: Warning, incorrect wavenumber provided."
      endif
    endif

    k2out = asix * two * (one - cos_prec(kin*dx)) &
          + four * bsix * half * (one - cos_prec(two*kin*dx)) &
          + nine * csix * (two / nine) * (one - cos_prec(three*kin*dx)) &
          + sixteen * dsix * (one / eight) * (one - cos_prec(four*kin*dx))
    k2out = k2out / (one + two * alsaix * cos_prec(kin*dx))

  end subroutine compute_k2

  ! Compute L1, L2 and Linf norm of given 3D array
  subroutine error_L1_L2_Linf(err, l1, l2, linf)

    USE decomp_2d
    USE MPI
    
    implicit none
      
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: err
    real(mytype),intent(out) :: l1, l2, linf

    integer :: i,j,k,code
    real(mytype) :: ll1, ll2, llinf, ntot

    ll1 = zero
    ll2 = zero
    llinf = zero
    ntot = nx*ny*nz

    do k = 1,xsize(3)
      do j = 1,xsize(2)
        do i = 1,xsize(1)
          ll1 = ll1 + abs(err(i,j,k))
          ll2 = ll2 + err(i,j,k)*err(i,j,k)
          llinf = max(llinf, abs(err(i,j,k)))
        enddo
      enddo
    enddo

    ! Parallel
    call MPI_ALLREDUCE(ll1,l1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(ll2,l2,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(llinf,linf,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)

    ! Rescaling
    l1 = l1 / ntot
    l2 = sqrt(l2 / ntot)

  end subroutine error_L1_L2_Linf

end module tgv
