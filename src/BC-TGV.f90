!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!!        FILE: BC-TGV.f90
!!!      AUTHOR: ??
!!!    MODIFIED: Paul Bartholomew
!!! DESCRIPTION: This module describes the Taylor-Green Vortex flow.
!!!   CHANGELOG: [2019-02-19] Making module private by default
!!               [2019-02-19] Turning file into a module
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module tgv

  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  real(mytype), save, allocatable, dimension(:,:,:) :: vol1,volSimps1
  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  !probes !só vai funcionar se a impressão for em relação ao lapis X!
  integer :: nprobes
  integer, save, allocatable, dimension(:) :: rankprobes, nxprobes, nyprobes, nzprobes

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_tgv, boundary_conditions_tgv, postprocessing_tgv

contains

  subroutine init_tgv (ux1,uy1,uz1,ep1,phi1,dux1,duy1,duz1,dphi1)

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime,numscalar) :: dphi1

    real(mytype) :: y,r,um,r3,x,z,h,ct
    real(mytype) :: cx0,cy0,cz0,hg,lg
    integer :: k,j,i,ijk,fh,ierror,is,code
    integer (kind=MPI_OFFSET_KIND) :: disp
    integer, dimension (:), allocatable :: seed
    integer ::  isize

    if (iscalar==1) then

       phi1(:,:,:,:) = zero
       dphi1(:,:,:,1,:) = phi1(:,:,:,:)
       do is = 2,ntime
          dphi1(:,:,:,is,:) = dphi1(:,:,:,is - 1,:)
       enddo
    endif

    if (iin.eq.0) then !empty domain

       if (nrank==0) write(*,*) "Empty initial domain!"

       ux1=zero; uy1=zero; uz1=zero

    endif

    if (iin.eq.1) then !generation of a random noise

       if (nrank==0) write(*,*) "Filled initial domain!"

       ux1=zero; uy1=zero; uz1=zero

       !t=zero
       !xv=1._mytype/10zero
       !xxk1=twopi/xlx
       !xxk2=twopi/yly
       do k=1,xsize(3)
          z=real((k+xstart(3)-1-1),mytype)*dz
          do j=1,xsize(2)
             y=real((j+xstart(2)-1-1),mytype)*dy
             do i=1,xsize(1)
                x=real(i-1,mytype)*dx
                !ux1(i,j,k)=sin(2.*pi*x)*cos(2.*pi*y)*cos(2.*pi*z)
                !sin(xxk1*x)*cos(xxk2*y)*exp(-(xxk1*xxk1+xxk2*xxk2)*xv*t)
                !uy1(i,j,k)=sin(2.*pi*y)*cos(2.*pi*x)*cos(2.*pi*z)
                !-xxk1/xxk2*sin(xxk2*y)*cos(xxk1*x)*exp(-(xxk1*xxk1+xxk2*xxk2)*xv*t)
                !uz1(i,j,k)=sin(2.*pi*z)*cos(2.*pi*x)*cos(2.*pi*y)
#ifndef TWOD
                ux1(i,j,k)=+sin(x)*cos(y)*cos(z)
                uy1(i,j,k)=-cos(x)*sin(y)*cos(z)
#else
                ux1(i,j,k)=+sin(x)*cos(y)
                uy1(i,j,k)=-cos(x)*sin(y)
#endif
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
             dux1(i,j,k,1)=ux1(i,j,k)
             duy1(i,j,k,1)=uy1(i,j,k)
             duz1(i,j,k,1)=uz1(i,j,k)
             do is = 2, ntime
                dux1(i,j,k,is)=dux1(i,j,k,is - 1)
                duy1(i,j,k,is)=duy1(i,j,k,is - 1)
                duz1(i,j,k,is)=duz1(i,j,k,is - 1)
             enddo
          enddo
       enddo
    enddo

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init end ok'
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

    IF (nclx1.EQ.2) THEN
    ENDIF
    IF (nclxn.EQ.2) THEN
    ENDIF

    IF (ncly1.EQ.2) THEN
    ENDIF
    IF (nclyn.EQ.2) THEN
    ENDIF

    IF (nclz1.EQ.2) THEN
    ENDIF
    IF (nclzn.EQ.2) THEN
    ENDIF

  end subroutine boundary_conditions_tgv

  !********************************************************************

  !############################################################################
  subroutine init_post(ep1)

    USE MPI

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ep1
    real(mytype) :: dxdydz, dxdz, x, xprobes, yprobes, zprobes
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

    !probes
    !WORK X-PENCILS
    open(10,file='probes.prm',status='unknown',form='formatted')
    read (10,*) nprobes
    read (10,*) a
    if (nprobes .gt. 0) then
       allocate(nxprobes(nprobes), nyprobes(nprobes), nzprobes(nprobes), rankprobes(nprobes))
       rankprobes(:)=0
       do i=1, nprobes
          read (10,*) xprobes, yprobes, zprobes
          !x
          if (nclx) then
             nxprobes(i)=int(xprobes/dx)
          else
             nxprobes(i)=int(xprobes/dx+1)
          end if
          !y
          if (ncly) then
             nyprobes(i)=int(yprobes/dy)
          else
             nyprobes(i)=int(yprobes/dy+1)
          end if
          !z
          if (nclz) then
             nzprobes(i)=int(zprobes/dz)
          else
             nzprobes(i)=int(zprobes/dz+1)
          end if
          if     (xstart(1) .le. nxprobes(i) .and. nxprobes(i) .le. xend(1)) then
             if   (xstart(2) .le. nyprobes(i) .and. nyprobes(i) .le. xend(2)) then
                if (xstart(3) .le. nzprobes(i) .and. nzprobes(i) .le. xend(3)) then
                   rankprobes(i)=1
                endif
             endif
          endif
       enddo
    endif
    close(10)

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init_post ok'
#endif

  end subroutine init_post
  !############################################################################
  subroutine write_probes(ux1,uy1,uz1,phi1) !By Felipe Schuch

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ux1, uy1, uz1
    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),numscalar) :: phi1

    integer :: i
    character(len=30) :: filename
    FS = 1+3+numscalar !Number of columns
    write(fileformat, '( "(",I4,"(E14.6),A)" )' ) FS
    FS = FS*14+1  !Line width

    do i=1, nprobes
       if (rankprobes(i) .eq. 1) then
          write(filename,"('./out/probe',I4.4)") i
          open(67,file=trim(filename),status='unknown',form='formatted'&
               ,access='direct',recl=FS)
          write(67,fileformat,rec=itime) t,&                         !1
               ux1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !2
               uy1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !3
               uz1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !4
               phi1(nxprobes(i),nyprobes(i),nzprobes(i),:),&         !numscalar
               NL                                                    !+1
          close(67)
       endif
    enddo

  end subroutine write_probes
  !############################################################################
  subroutine postprocessing_tgv(ux1,uy1,uz1,phi1,ep1)

    USE decomp_2d
    USE MPI

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,ep1,diss1
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
    real(mytype) :: mp(numscalar),mps(numscalar),vl,es,es1,ek,ek1,ds,ds1
    real(mytype) :: temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9

    real(mytype) :: eek, enst, eps, eps2
    integer :: nxc, nyc, nzc, xsize1, xsize2, xsize3

    integer :: i,j,k,is,ijk,code,nvect1
    nvect1=xsize(1)*xsize(2)*xsize(3)

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

    !SPATIALLY-AVERAGED TKE of velocity fields (filtered or not)
    temp1=0.
    temp2=0.
    temp3=0.
    do k=1,xsize3
       do j=1,xsize2
          do i=1,xsize1
             temp1=temp1+0.5*((ux1(i,j,k))**2+(uy1(i,j,k))**2+(uz1(i,j,k))**2) 
          enddo
       enddo
    enddo
    call MPI_ALLREDUCE(temp1,eek,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    eek=eek/(nxc*nyc*nzc)

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

    !SPATIALLY-AVERAGED ENSTROPHY
    temp1=0.
    do k=1,xsize3
       do j=1,xsize2
          do i=1,xsize1
             temp1=temp1+0.5*((tf1(i,j,k)-th1(i,j,k))**2+&
                  (tg1(i,j,k)-tc1(i,j,k))**2+&
                  (tb1(i,j,k)-td1(i,j,k))**2)
          enddo
       enddo
    enddo
    call MPI_ALLREDUCE(temp1,enst,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    enst=enst/(nxc*nyc*nzc)

    !SPATIALLY-AVERAGED ENERGY DISSIPATION
    temp1=0.
    do k=1,xsize3
       do j=1,xsize2
          do i=1,xsize1
             temp1=temp1+0.5*xnu*((2.*ta1(i,j,k))**2+(2.*te1(i,j,k))**2+&
                  (2.*ti1(i,j,k))**2+2.*(td1(i,j,k)+tb1(i,j,k))**2+&
                  2.*(tg1(i,j,k)+tc1(i,j,k))**2+&
                  2.*(th1(i,j,k)+tf1(i,j,k))**2)
          enddo
       enddo
    enddo
    call MPI_ALLREDUCE(temp1,eps,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    eps=eps/(nxc*nyc*nzc)

    !du/dx=ta1 du/dy=td1 and du/dz=tg1
    !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
    !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1
    !SPATIALLY-AVERAGED SQUARE VELOCITY GRADIENTS
    temp1=0.
    temp2=0.
    temp3=0.
    temp4=0.
    temp5=0.
    temp6=0.
    temp7=0.
    temp8=0.
    temp9=0.
    do k=1,xsize3
       do j=1,xsize2
          do i=1,xsize1
             temp1=temp1+ta1(i,j,k)*ta1(i,j,k)
             temp2=temp2+td1(i,j,k)*td1(i,j,k)
             temp3=temp3+tg1(i,j,k)*tg1(i,j,k)
             temp4=temp4+tb1(i,j,k)*tb1(i,j,k)
             temp5=temp5+te1(i,j,k)*te1(i,j,k)
             temp6=temp6+th1(i,j,k)*th1(i,j,k)
             temp7=temp7+tc1(i,j,k)*tc1(i,j,k)
             temp8=temp8+tf1(i,j,k)*tf1(i,j,k)
             temp9=temp9+ti1(i,j,k)*ti1(i,j,k)
          enddo
       enddo
    enddo

    !!SECOND DERIVATIVES
    !x-derivatives
    call derxx (ta1,ux1,di1,sx,sfx ,ssx ,swx ,xsize(1),xsize(2),xsize(3),0)
    call derxx (tb1,uy1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
    call derxx (tc1,uz1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
    !y-derivatives
    call transpose_x_to_y(ux1,td2)
    call transpose_x_to_y(uy1,te2)
    call transpose_x_to_y(uz1,tf2)
    call deryy (ta2,td2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1) 
    call deryy (tb2,te2,di2,sy,sfy ,ssy ,swy ,ysize(1),ysize(2),ysize(3),0) 
    call deryy (tc2,tf2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1) 
    !!z-derivatives
    call transpose_y_to_z(td2,td3)
    call transpose_y_to_z(te2,te3)
    call transpose_y_to_z(tf2,tf3)
    call derzz(ta3,td3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
    call derzz(tb3,te3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
    call derzz(tc3,tf3,di3,sz,sfz ,ssz ,swz ,zsize(1),zsize(2),zsize(3),0)
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
    temp1=0.
    do k=1,xsize3
       do j=1,xsize2
          do i=1,xsize1
             di1(i,j,k)=(-xnu)*(ux1(i,j,k)*(ta1(i,j,k)+td1(i,j,k)+tg1(i,j,k))+ &
                  uy1(i,j,k)*(tb1(i,j,k)+te1(i,j,k)+th1(i,j,k))+ &
                  uz1(i,j,k)*(tc1(i,j,k)+tf1(i,j,k)+ti1(i,j,k)))
             temp1=temp1+di1(i,j,k)
          enddo
       enddo
    enddo
    call MPI_ALLREDUCE(temp1,eps2,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    eps2=eps2/(nxc*nyc*nzc)

    if (nrank==0) then
       write(42,'(20e20.12)') (itime-1)*dt,eek,eps,enst,eps2
       call flush(42)
    endif

  end subroutine postprocessing_tgv
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
end module tgv
