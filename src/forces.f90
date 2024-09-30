!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

!=======================================================================
! This program computes the drag and lift coefficients alongo a 
! cylinder by the control ! volume (CV) technique for 2D (pencil) 
! decomposition. 
!
! Adpated from Leandro Pinto PhD Thesis (2012) by Gabriel Narvaez Campo 
! 08-2018 Nucleo de Estudos em Transicao e Turbulencia (NETT/IPH/UFRGS)
!
!=======================================================================

module forces

  use decomp_2d_constants
  use decomp_2d_mpi
  USE decomp_2d

  implicit none

  integer,save :: nvol,iforces,i2dsim
  real(mytype),save,allocatable,dimension(:,:,:) :: ux01, uy01, ux11, uy11, ppi1
  real(mytype),save,allocatable,dimension(:) :: xld,xrd,yld,yud,zfr,zbk
  integer,save,allocatable,dimension(:) :: icvlf,icvrt,jcvlw,jcvup,kcvfr,kcvbk
  integer,save,allocatable,dimension(:) :: icvlf_lx,icvrt_lx,icvlf_ly,icvrt_ly
  integer,save,allocatable,dimension(:) :: jcvlw_lx,jcvup_lx,jcvlw_ly,jcvup_ly
  integer,save,allocatable,dimension(:) :: kcvfr_lx,kcvbk_lx,kcvfr_ly,kcvbk_ly

  character(len=*), parameter :: io_restart_forces = "restart-forces-io", &
       resfile = "restart-forces"
  
contains

  subroutine init_forces

    USE decomp_2d_io, only : decomp_2d_register_variable, decomp_2d_init_io
    USE param
    USE variables
    implicit none

    integer :: iv,stp1,stp2,h

    call alloc_x(ux01)
    call alloc_x(uy01)
    call alloc_x(ux11)
    call alloc_x(uy11)
    call alloc_x(ppi1)

    ux01 = zero
    uy01 = zero
    ux11 = zero
    uy11 = zero

    allocate(icvlf(nvol),icvrt(nvol),jcvlw(nvol),jcvup(nvol),kcvfr(nvol),kcvbk(nvol))
    allocate(icvlf_lx(nvol),icvrt_lx(nvol),icvlf_ly(nvol),icvrt_ly(nvol))
    allocate(jcvlw_lx(nvol),jcvup_lx(nvol),jcvlw_ly(nvol),jcvup_ly(nvol))
    allocate(kcvfr_lx(nvol),kcvbk_lx(nvol),kcvfr_ly(nvol),kcvbk_ly(nvol))

    !     Definition of the Control Volume
    !*****************************************************************
    !! xld,xrd,yld,yud: limits of control volume (!!don't use cex and cey anymore!!)

    do iv=1,nvol
       ! ok for istret=0 (!!to do for istret=1!!)
       icvlf(iv) = nint(xld(iv)/dx)+1
       icvrt(iv) = nint(xrd(iv)/dx)+1
       kcvfr(iv) = nint(zfr(iv)/dz)+1
       kcvbk(iv) = nint(zbk(iv)/dz)+1
       if (istret.eq.0) then 
         jcvlw(iv) = nint(yld(iv)/dy)+1
         jcvup(iv) = nint(yud(iv)/dy)+1
       else
         stp1=0
         stp2=0
         do h = 1, ny-1  
           if ((-yp(h+1)-yp(h)+two*yld(iv)).lt.(yld(iv)-yp(h)).and.(stp1.eq.0)) then 
             jcvlw(iv) = h+1
             stp1=1
           endif
           if ((-yp(h+1)-yp(h)+two*yud(iv)).lt.(yud(iv)-yp(h)).and.(stp2.eq.0)) then
             jcvup(iv) = h
             stp2=1 
           endif
         enddo
       endif
       icvlf_lx(iv) = icvlf(iv)
       icvrt_lx(iv) = icvrt(iv)
       jcvlw_lx(iv) = max(jcvlw(iv)+1-xstart(2),1)
       jcvup_lx(iv) = min(jcvup(iv)+1-xstart(2),xsize(2))
       kcvfr_lx(iv) = max(kcvfr(iv)+1-xstart(3),1)
       kcvbk_lx(iv) = min(kcvbk(iv)+1-xstart(3),xsize(3))
       icvlf_ly(iv) = max(icvlf(iv)+1-ystart(1),1)
       icvrt_ly(iv) = min(icvrt(iv)+1-ystart(1),ysize(1))
       jcvlw_ly(iv) = jcvlw(iv)
       jcvup_ly(iv) = jcvup(iv)
       kcvfr_ly(iv) = max(kcvfr(iv)+1-ystart(3),1)
       kcvbk_ly(iv) = min(kcvbk(iv)+1-ystart(3),ysize(3))
    enddo

    if (nrank==0) then
       if (i2dsim==1) then
         write(*,*) '========================Forces============================='
         write(*,*) '     (icvlf)      (icvrt) '
         write(*,*) '  (jcvup) B____________C  '
         write(*,*) '          \            \  '
         write(*,*) '          \     __     \  '
         write(*,*) '          \    \__\    \  '
         write(*,*) '          \            \  '
         write(*,*) '          \       CV   \  '
         write(*,*) '  (jcvlw) A____________D  '
         do iv=1,nvol
            write(*,"(' Control Volume     : #',I1)") iv
            write(*,"('     xld, icvlf     : (',F6.2,',',I6,')')") xld(iv), icvlf(iv)
            write(*,"('     xrd, icvrt     : (',F6.2,',',I6,')')") xrd(iv), icvrt(iv)
            write(*,"('     yld, jcvlw     : (',F6.2,',',I6,')')") yld(iv), jcvlw(iv)
            write(*,"('     yud, jcvup     : (',F6.2,',',I6,')')") yud(iv), jcvup(iv)
         enddo
         write(*,*) '==========================================================='
       elseif (i2dsim==0) then
         write(*,*) '========================Forces============================='
         write(*,*) '     (icvlf)      (icvrt)    (kcvbk)    (kcvfr)'
         write(*,*) '  (jcvup) B____________C     B`_____________B  '
         write(*,*) '          \            \  |  \              \  '
         write(*,*) '          \     __     \  |  \     ____     \  '
         write(*,*) '          \    \__\    \  |  \     \___\    \  '
         write(*,*) '          \            \  |  \              \  '
         write(*,*) '          \       CV   \  |  \    (Front)   \  '
         write(*,*) '  (jcvlw) A____________D  |  A`_____________A  '
         do iv=1,nvol
            write(*,"(' Control Volume     : #',I1)") iv
            write(*,"('     xld, icvlf     : (',F6.2,',',I6,')')") xld(iv), icvlf(iv)
            write(*,"('     xrd, icvrt     : (',F6.2,',',I6,')')") xrd(iv), icvrt(iv)
            write(*,"('     yld, jcvlw     : (',F6.2,',',I6,')')") yld(iv), jcvlw(iv)
            write(*,"('     yud, jcvup     : (',F6.2,',',I6,')')") yud(iv), jcvup(iv)
            write(*,"('     zfr, kcvfr     : (',F6.2,',',I6,')')") zfr(iv), kcvfr(iv)
            write(*,"('     zbk, kcvbk     : (',F6.2,',',I6,')')") zbk(iv), kcvbk(iv)
            enddo
         write(*,*) '==========================================================='
       endif
    endif

    call decomp_2d_init_io(io_restart_forces)
    call decomp_2d_register_variable(io_restart_forces, "ux01", 1, 0, 0, mytype)
    call decomp_2d_register_variable(io_restart_forces, "uy01", 1, 0, 0, mytype)
    call decomp_2d_register_variable(io_restart_forces, "ux11", 1, 0, 0, mytype)
    call decomp_2d_register_variable(io_restart_forces, "uy11", 1, 0, 0, mytype)
    
  end subroutine init_forces

  !
  ! Allocate 1D arrays and initialize variables before reading the forces namelist
  !
  subroutine setup_forces(iounit)

    implicit none

    ! Argument
    integer, intent(in) :: iounit

    NAMELIST /ForceCVs/ xld, xrd, yld, yud, zfr, zbk, i2dsim

    ! Safety check
    if (allocated(xld)) then
       call decomp_2d_abort(1, "Error in setup_forces")
    end if
    if (nvol < 1) then
       call decomp_2d_abort(nvol, "Invalid nvol in setup_forces")
    end if

    ! Allocate 1D arrays
    allocate(xld(nvol), xrd(nvol), yld(nvol), yud(nvol), zfr(nvol), zbk(nvol))

    ! Default values in the forces namelist
    xld = 0._mytype
    xrd = 0._mytype
    yld = 0._mytype
    yud = 0._mytype
    zfr = 0._mytype
    zbk = 0._mytype
    i2dsim = 1

    ! Read a part of the namelist and rewind
    read(iounit, nml=ForceCVs)
    rewind(iounit)

    ! Safety check
    if (i2dsim < 0 .or. i2dsim > 1) then
       call decomp_2d_abort(i2dsim, "Invalid value for the parameter i2dsim")
    end if

  end subroutine setup_forces

  subroutine restart_forces(itest1)

    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none

    ! Argument
    integer, intent(in) :: itest1

    ! Exit if writing and invalid time step
    if (itest1 == 1 .and. mod(itime, icheckpoint).ne.0) then
       return
    end if

    if (itest1 == 1) then
       call decomp_2d_open_io(io_restart_forces, resfile, decomp_2d_write_mode)
    else
       call decomp_2d_open_io(io_restart_forces, resfile, decomp_2d_read_mode)
    endif
    call decomp_2d_start_io(io_restart_forces, resfile)
    
    if (itest1==1) then
       !write
       call decomp_2d_write_one(1,ux01,resfile,"ux01",0,io_restart_forces)
       call decomp_2d_write_one(1,uy01,resfile,"uy01",0,io_restart_forces)
       call decomp_2d_write_one(1,ux11,resfile,"ux11",0,io_restart_forces)
       call decomp_2d_write_one(1,uy11,resfile,"uy11",0,io_restart_forces)
    else
       !read
       call decomp_2d_read_one(1,ux01,resfile,"ux01",io_restart_forces)
       call decomp_2d_read_one(1,uy01,resfile,"uy01",io_restart_forces)
       call decomp_2d_read_one(1,ux11,resfile,"ux11",io_restart_forces)
       call decomp_2d_read_one(1,uy11,resfile,"uy11",io_restart_forces)
    endif

    call decomp_2d_end_io(io_restart_forces, resfile)
    call decomp_2d_close_io(io_restart_forces, resfile)
    
  end subroutine restart_forces

  subroutine force(ux1,uy1,uz1,ep1)

    USE param
    USE variables
    USE MPI
    USE ibm_param

    use var, only : ta1, tb1, tc1, td1, te1, di1, tg1, tg2, tg3, th1, th2, th3, tf2, tf1
    use var, only : ux2, ux3, uy2, uy3, uz2, ta2, tb2, td2, te2, di2, di3
    use var, only : tc2

    implicit none
    character(len=30) :: filename, filename2
    integer :: nzmsize
    integer                                             :: i, iv, j, k, kk, code, jj
    integer                                             :: nvect1,nvect2,nvect3
    integer :: ierror

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1, uy1, uz1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ep1

    !real(mytype), dimension(ysize(1),ysize(2),ysize(3)) :: ppi2 ! we'll use tc2

    real(mytype), dimension(nz) :: yLift,xDrag
    real(mytype) :: yLift_mean,xDrag_mean

    real(mytype), dimension(ny-1) :: del_y

    real(mytype), dimension(nz) :: tunstx,tunsty
    real(mytype), dimension(nz) :: tconvx,tconvy
    real(mytype), dimension(nz) :: tpresx,tpresy
    real(mytype), dimension(nz) :: tdiffx,tdiffy

    real(mytype) :: uxmid,uymid,uzmid,prmid
    real(mytype) :: dudxmid,dudymid,dvdxmid,dvdymid,dudzmid,dwdxmid,dwdymid,dvdzmid
    real(mytype) :: fac,tsumx,tsumy
    real(mytype) :: fcvx,fcvy,fprx,fpry,fdix,fdiy
    real(mytype) :: xmom,ymom
    real(mytype) :: convx,convy,pressx,pressy,stressx,stressy

    nvect1=xsize(1)*xsize(2)*xsize(3)
    nvect2=ysize(1)*ysize(2)*ysize(3)
    nvect3=zsize(1)*zsize(2)*zsize(3)

    do jj = 1, ny-1
       if (istret.eq.0) then
          del_y(jj)=dy
       else
          del_y(jj)=yp(jj+1)-yp(jj) 
       endif
    enddo

    if (itime.eq.1) then
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                ux11(i,j,k)=ux1(i,j,k)
                uy11(i,j,k)=uy1(i,j,k)
             enddo
          enddo
       enddo
       return
    elseif (itime.eq.2) then
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                ux01(i,j,k)=ux1(i,j,k)
                uy01(i,j,k)=uy1(i,j,k)
             enddo
          enddo
       enddo
       return
    endif

    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)    ! dudx
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy) ! dvdx
    call transpose_x_to_y(ta1,ta2) ! dudx
    call transpose_x_to_y(tb1,tb2) ! dvdx

    call transpose_x_to_y(ux1,ux2)
    call transpose_x_to_y(uy1,uy2)
    call transpose_x_to_y(ppi1,tc2)

    call dery (td2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx) ! dudy
    call dery (te2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)    ! dvdy
    call transpose_y_to_x(td2,td1) ! dudy
    call transpose_y_to_x(te2,te1) ! dvdy

    if (i2dsim==0) then
      call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz) ! dwdx

      call transpose_y_to_z(ux2, ux3)
      call transpose_y_to_z(uy2, uy3)

      call transpose_x_to_y(uz1, uz2)
      call dery (tf2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz) ! dwdy
      call transpose_y_to_x(tf2,tf1) ! dwdy

      call derz (tg3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcx) !dudz
      call derz (th3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcy) !dvdz
      call transpose_z_to_y(tg3,tg2)
      call transpose_y_to_x(tg2,tg1)
      call transpose_z_to_y(th3,th2)
      call transpose_y_to_x(th2,th1)
    endif

    !*****************************************************************
    !      Drag and Lift coefficients
    !*****************************************************************
    if(i2dsim==1) then
      do iv=1,nvol

         !*****************************************************************
         !        Calculation of the momentum terms
         !*****************************************************************
         !
         !     Calculation of the momentum terms. First we integrate the 
         !     time rate of momentum along the CV.
         !
         !         Excluding the body internal cells. If the centroid 
         !         of the cell falls inside the body the cell is 
         !         excluded.

         tunstx=zero
         tunsty=zero
         do k=1,xsize(3)
            tsumx=zero
            tsumy=zero
            do j=jcvlw_lx(iv),jcvup_lx(iv)
               do i=icvlf_lx(iv),icvrt_lx(iv)
                  !     The velocity time rate has to be relative to the cell center, 
                  !     and not to the nodes, because, here, we have an integral 
                  !     relative to the volume, and, therefore, this has a sense 
                  !     of a "source".
                  !         fac   = (1.5*ux1(i,j,k)-2.0*ux01(i,j,k)+0.5*ux11(i,j,k))*epcv1(i,j,k)
                  fac   = (onepfive*ux1(i,j,k)-two*ux01(i,j,k)+half*ux11(i,j,k))*(one-ep1(i,j,k))
                  tsumx = tsumx+fac*dx*del_y(j+(xstart(2)-1))/dt    !tsumx+fac*dx*dy/dt
                  !sumx(k) = sumx(k)+dudt1*dx*dy

                  !         fac   = (1.5*uy1(i,j,k)-2.0*uy01(i,j,k)+0.5*uy11(i,j,k))*epcv1(i,j,k)
                  fac   = (onepfive*uy1(i,j,k)-two*uy01(i,j,k)+half*uy11(i,j,k))*(one-ep1(i,j,k))
                  tsumy = tsumy+fac*dx*del_y(j+(xstart(2)-1))/dt !tsumy+fac*dx*dy/dt
                  !sumy(k) = sumy(k)+dudt1*dx*dy
               enddo
            enddo
            tunstx(xstart(3)-1+k)=tsumx
            tunsty(xstart(3)-1+k)=tsumy
         enddo
         call MPI_ALLREDUCE(MPI_IN_PLACE,tunstx,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
         if (code /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_ALLREDUCE")
         end if
         call MPI_ALLREDUCE(MPI_IN_PLACE,tunsty,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
         if (code /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_ALLREDUCE")
         end if

         !!$!*********************************************************************************
         !!$!     Secondly, the surface momentum fluxes
         !!$!*********************************************************************************
         !!$
         !!$!       (icvlf)      (icvrt)
         !!$!(jcvup) B____________C  
         !!$!        \            \
         !!$!        \     __     \
         !!$!        \    \__\    \
         !!$!        \            \
         !!$!        \       CV   \
         !!$!(jcvlw) A____________D      

         tconvx=zero
         tconvy=zero
         tdiffx=zero
         tdiffy=zero
         tpresx=zero
         tpresy=zero
         !BC and AD : x-pencils
         !AD
         if ((jcvlw(iv).ge.xstart(2)).and.(jcvlw(iv).le.xend(2))) then
            j=jcvlw(iv)-xstart(2)+1
            do k=1,xsize(3)
               kk=xstart(3)-1+k
               fcvx=zero
               fcvy=zero
               fpry=zero
               fdix=zero
               fdiy=zero
               do i=icvlf_lx(iv),icvrt_lx(iv)-1
                  !momentum flux
                  !FIXME avoid interpolation for the non-linear term
                  uxmid = half*(ux1(i,j,k)+ux1(i+1,j,k))
                  uymid = half*(uy1(i,j,k)+uy1(i+1,j,k))
                  fcvx  = fcvx -uxmid*uymid*dx
                  fcvy  = fcvy -uymid*uymid*dx

                  !pressure
                  prmid = half*(ppi1(i,j,k)+ppi1(i+1,j,k))
                  fpry  = fpry +prmid*dx

                  !viscous term
                  dudymid = half*(td1(i,j,k)+td1(i+1,j,k))
                  dvdxmid = half*(tb1(i,j,k)+tb1(i+1,j,k))
                  dvdymid = half*(te1(i,j,k)+te1(i+1,j,k))
                  fdix    = fdix -(xnu*(dudymid+dvdxmid)*dx)
                  fdiy    = fdiy -two*xnu*dvdymid*dx

               enddo

               tconvx(kk)=tconvx(kk)+fcvx
               tconvy(kk)=tconvy(kk)+fcvy
               tpresy(kk)=tpresy(kk)+fpry
               tdiffx(kk)=tdiffx(kk)+fdix
               tdiffy(kk)=tdiffy(kk)+fdiy
            enddo
         endif
         !BC
         if ((jcvup(iv).ge.xstart(2)).and.(jcvup(iv).le.xend(2))) then
            j=jcvup(iv)-xstart(2)+1
            do k=1,xsize(3)
               kk=xstart(3)-1+k
               fcvx=zero
               fcvy=zero
               fpry=zero
               fdix=zero
               fdiy=zero
               do i=icvlf_lx(iv),icvrt_lx(iv)-1
                  !momentum flux
                  !FIXME avoid interpolation for the non-linear term
                  uxmid = half*(ux1(i,j,k)+ux1(i+1,j,k))
                  uymid = half*(uy1(i,j,k)+uy1(i+1,j,k))
                  fcvx= fcvx +uxmid*uymid*dx
                  fcvy= fcvy +uymid*uymid*dx

                  !pressure
                  prmid = half*(ppi1(i,j,k)+ppi1(i+1,j,k))
                  fpry = fpry -prmid*dx

                  !viscous term
                  dudymid = half*(td1(i,j,k)+td1(i+1,j,k))
                  dvdxmid = half*(tb1(i,j,k)+tb1(i+1,j,k))
                  dvdymid = half*(te1(i,j,k)+te1(i+1,j,k))
                  fdix = fdix +(xnu*(dudymid+dvdxmid)*dx)
                  fdiy = fdiy +two*xnu*dvdymid*dx

               enddo
               tconvx(kk)=tconvx(kk)+fcvx
               tconvy(kk)=tconvy(kk)+fcvy
               tpresy(kk)=tpresy(kk)+fpry
               tdiffx(kk)=tdiffx(kk)+fdix
               tdiffy(kk)=tdiffy(kk)+fdiy
            enddo
         endif
         !AB and DC : y-pencils
         !AB
         if ((icvlf(iv).ge.ystart(1)).and.(icvlf(iv).le.yend(1))) then
            i=icvlf(iv)-ystart(1)+1
            do k=1,ysize(3)
               kk=ystart(3)-1+k
               fcvx=zero
               fcvy=zero
               fprx=zero
               fdix=zero
               fdiy=zero
               do j=jcvlw_ly(iv),jcvup_ly(iv)-1
                  !momentum flux
                  !FIXME avoid interpolation for the non-linear term
                  uxmid = half*(ux2(i,j,k)+ux2(i,j+1,k))
                  uymid = half*(uy2(i,j,k)+uy2(i,j+1,k))
                  fcvx= fcvx -uxmid*uxmid*del_y(j)
                  fcvy= fcvy -uxmid*uymid*del_y(j)

                  !pressure
                  prmid=half*(tc2(i,j,k)+tc2(i,j+1,k))
                  fprx = fprx +prmid*del_y(j)

                  !viscous term
                  dudxmid = half*(ta2(i,j,k)+ta2(i,j+1,k))
                  dudymid = half*(td2(i,j,k)+td2(i,j+1,k))
                  dvdxmid = half*(tb2(i,j,k)+tb2(i,j+1,k))
                  fdix = fdix -two*xnu*dudxmid*del_y(j)
                  fdiy = fdiy -xnu*(dvdxmid+dudymid)*del_y(j)
               enddo
               tconvx(kk)=tconvx(kk)+fcvx
               tconvy(kk)=tconvy(kk)+fcvy
               tpresx(kk)=tpresx(kk)+fprx
               tdiffx(kk)=tdiffx(kk)+fdix
               tdiffy(kk)=tdiffy(kk)+fdiy
            enddo
         endif
         !DC
         if ((icvrt(iv).ge.ystart(1)).and.(icvrt(iv).le.yend(1))) then
            i=icvrt(iv)-ystart(1)+1
            do k=1,ysize(3)
               kk=ystart(3)-1+k
               fcvx=zero
               fcvy=zero
               fprx=zero
               fdix=zero
               fdiy=zero
               do j=jcvlw_ly(iv),jcvup_ly(iv)-1
                  !momentum flux
                  !FIXME avoid interpolation for the non-linear term
                  uxmid = half*(ux2(i,j,k)+ux2(i,j+1,k))
                  uymid = half*(uy2(i,j,k)+uy2(i,j+1,k))
                  fcvx= fcvx +uxmid*uxmid*del_y(j)
                  fcvy= fcvy +uxmid*uymid*del_y(j)

                  !pressure
                  prmid=half*(tc2(i,j,k)+tc2(i,j+1,k))
                  fprx = fprx -prmid*del_y(j)

                  !viscous term
                  dudxmid = half*(ta2(i,j,k)+ta2(i,j+1,k))
                  dudymid = half*(td2(i,j,k)+td2(i,j+1,k))
                  dvdxmid = half*(tb2(i,j,k)+tb2(i,j+1,k))
                  fdix = fdix +two*xnu*dudxmid*del_y(j)
                  fdiy = fdiy +xnu*(dvdxmid+dudymid)*del_y(j)
               enddo
               tconvx(kk)=tconvx(kk)+fcvx
               tconvy(kk)=tconvy(kk)+fcvy
               tpresx(kk)=tpresx(kk)+fprx
               tdiffx(kk)=tdiffx(kk)+fdix
               tdiffy(kk)=tdiffy(kk)+fdiy
            enddo
         endif
         call MPI_ALLREDUCE(MPI_IN_PLACE,tconvx,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
         if (code /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_ALLREDUCE!")
         end if
         call MPI_ALLREDUCE(MPI_IN_PLACE,tconvy,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
         if (code /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_ALLREDUCE!")
         end if
         call MPI_ALLREDUCE(MPI_IN_PLACE,tpresx,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
         if (code /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_ALLREDUCE!")
         end if
         call MPI_ALLREDUCE(MPI_IN_PLACE,tpresy,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
         if (code /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_ALLREDUCE!")
         end if
         call MPI_ALLREDUCE(MPI_IN_PLACE,tdiffx,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
         if (code /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_ALLREDUCE!")
         end if
         call MPI_ALLREDUCE(MPI_IN_PLACE,tdiffy,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
         if (code /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_ALLREDUCE!")
         end if

         do k=1,zsize(3)

            tpresx(k)=tpresx(k)/dt
            tpresy(k)=tpresy(k)/dt

            xmom    = tunstx(k)+tconvx(k)
            ymom    = tunsty(k)+tconvy(k)
            xDrag(k) = two*(tdiffx(k)+tpresx(k)-xmom)
            yLift(k) = two*(tdiffy(k)+tpresy(k)-ymom)

         enddo

         !Edited by F. Schuch
         xDrag_mean = sum(xDrag(:))/real(nz,mytype)
         yLift_mean = sum(yLift(:))/real(nz,mytype)

         !     if ((itime==ifirst).or.(itime==0)) then
         !        if (nrank .eq. 0) then
         !        write(filename,"('aerof',I1.1)") iv
         !        open(38+(iv-1),file=filename,status='unknown',form='formatted')
         !        endif
         !     endif
         if (nrank .eq. 0) then
            write(38,*) t,xDrag_mean,yLift_mean
            flush(38)
         endif
      enddo
    elseif(i2dsim==0) then
      do iv=1, nvol

         !*****************************************************************
         !        3D Control Volume Method (Added by Gaurav Gupta, IIST India)
         !*****************************************************************
         !
         !The following code outputs drag and lift force computed for a 3D 
         !object like sphere and coefficients of drag and lift can be 
         !calculated by dividing the forces by (0.5*rho*v*v*A). 
 
 
         tsumx=zero
         tsumy=zero
         do k=kcvfr_lx(iv), kcvbk_lx(iv)
             do j=jcvlw_lx(iv),jcvup_lx(iv)
               do i=icvlf_lx(iv),icvrt_lx(iv)
                   !     The velocity time rate has to be relative to the cell center, 
                   !     and not to the nodes, because, here, we have an integral 
                   !     relative to the volume, and, therefore, this has a sense 
                   !     of a "source".
                   !         fac   = (1.5*ux1(i,j,k)-2.0*ux01(i,j,k)+0.5*ux11(i,j,k))*epcv1(i,j,k)
                   fac   = (onepfive*ux1(i,j,k)-two*ux01(i,j,k)+half*ux11(i,j,k))*(one-ep1(i,j,k))
                   tsumx = tsumx+fac*dx*del_y(j+(xstart(2)-1))*dz/dt    !tsumx+fac*dx*dy*dz/dt
                   !sumx(k) = sumx(k)+dudt1*dx*dy
 
                   !         fac   = (1.5*uy1(i,j,k)-2.0*uy01(i,j,k)+0.5*uy11(i,j,k))*epcv1(i,j,k)
                   fac   = (onepfive*uy1(i,j,k)-two*uy01(i,j,k)+half*uy11(i,j,k))*(one-ep1(i,j,k))
                   tsumy = tsumy+fac*dx*del_y(j+(xstart(2)-1))*dz/dt !tsumy+fac*dx*dy*dz/dt
                   !sumy(k) = sumy(k)+dudt1*dx*dy
               enddo
             enddo
         enddo
 
         call MPI_ALLREDUCE(tsumx,xmom,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
         if (code /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_ALLREDUCE!")
         end if
         call MPI_ALLREDUCE(tsumy,ymom,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
         if (code /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_ALLREDUCE!")
         end if
 
         convx=zero
         convy=zero
         pressx=zero
         pressy=zero
         stressx=zero
         stressy=zero
         
         !y-pencil
         !calculating the flux through inlet wall (x-direction)
         if((icvlf(iv).ge.ystart(1)).and.(icvlf(iv).le.yend(1))) then
             i=icvlf(iv)-ystart(1)+1
             fcvx=zero
             fcvy=zero
             fprx=zero
             fdix=zero
             fdiy=zero
             do k=kcvfr_ly(iv),kcvbk_ly(iv)
               do j=jcvlw_ly(iv),jcvup_ly(iv)-1
                   !momentum flux
                   !FIXME avoid interpolation for the non-linear term
                   uxmid = half*(ux2(i,j,k) + ux2(i,j+1,k))
                   uymid = half*(uy2(i,j,k) + uy2(i,j+1,k))
                   fcvx = fcvx + uxmid*uxmid*del_y(j)*dz
                   fcvy = fcvy + uxmid*uymid*del_y(j)*dz
 
                   !pressure
                   prmid=half*(tc2(i,j,k)+tc2(i,j+1,k))
                   fprx = fprx -prmid*del_y(j)*dz
 
                   !viscous term
                   dudxmid = half*(ta2(i,j,k)+ta2(i,j+1,k))
                   dudymid = half*(td2(i,j,k)+td2(i,j+1,k))
                   dvdxmid = half*(tb2(i,j,k)+tb2(i,j+1,k))
                   fdix = fdix + two*xnu*dudxmid*del_y(j)*dz
                   fdiy = fdiy + xnu*(dvdxmid+dudymid)*del_y(j)*dz
               enddo
             enddo
             convx = convx + fcvx
             pressx = pressx + fprx
             stressx = stressx + fdix
             convy = convy + fcvy
             stressy = stressy + fdiy 
         endif
 
         !calculating the flux through outlet wall (x-direction)
         if((icvrt(iv).ge.ystart(1)).and.(icvrt(iv).le.yend(1))) then
             i=icvrt(iv)-ystart(1)+1
             fcvx=zero
             fcvy=zero
             fprx=zero
             fdix=zero
             fdiy=zero
             do k=kcvfr_ly(iv),kcvbk_ly(iv)
               do j=jcvlw_ly(iv),jcvup_ly(iv)-1
                   !momentum flux
                   !FIXME avoid interpolation for the non-linear term
                   uxmid = half*(ux2(i,j,k) + ux2(i,j+1,k))
                   uymid = half*(uy2(i,j,k) + uy2(i,j+1,k))
                   fcvx = fcvx - uxmid*uxmid*del_y(j)*dz
                   fcvy = fcvy - uxmid*uymid*del_y(j)*dz
 
                   !pressure
                   prmid=half*(tc2(i,j,k)+tc2(i,j+1,k))
                   fprx = fprx + prmid*del_y(j)*dz
 
                   !viscous term
                   dudxmid = half*(ta2(i,j,k)+ta2(i,j+1,k))
                   dudymid = half*(td2(i,j,k)+td2(i,j+1,k))
                   dvdxmid = half*(tb2(i,j,k)+tb2(i,j+1,k))
                   fdix = fdix - two*xnu*dudxmid*del_y(j)*dz
                   fdiy = fdiy -xnu*(dvdxmid+dudymid)*del_y(j)*dz
               enddo
             enddo
             convx = convx + fcvx
             pressx = pressx + fprx
             stressx = stressx + fdix
             convy = convy + fcvy
             stressy = stressy + fdiy
         endif
 
         !x-pencil
         !calculating the flux through top wall (y-direction)
         if ((jcvup(iv).ge.xstart(2)).and.(jcvup(iv).le.xend(2))) then
             j=jcvup(iv)-xstart(2)+1
             fcvx=zero
             fcvy=zero
             fpry=zero
             fdix=zero
             fdiy=zero
             do k=kcvfr_lx(iv),kcvbk_lx(iv)
               do i=icvlf_lx(iv),icvrt_lx(iv)-1
                   !momentum flux
                   !FIXME avoid interpolation for the non-linear term
                   uxmid = half*(ux1(i,j,k)+ux1(i+1,j,k))
                   uymid = half*(uy1(i,j,k)+uy1(i+1,j,k))
                   fcvx= fcvx -uxmid*uymid*dx*dz
                   fcvy= fcvy -uymid*uymid*dx*dz
 
                   !pressure
                   prmid = half*(ppi1(i,j,k)+ppi1(i+1,j,k))
                   fpry = fpry +prmid*dx*dz
 
                   !viscous term
                   dudymid = half*(td1(i,j,k)+td1(i+1,j,k))
                   dvdxmid = half*(tb1(i,j,k)+tb1(i+1,j,k))
                   dvdymid = half*(te1(i,j,k)+te1(i+1,j,k))
                   fdix = fdix -(xnu*(dudymid+dvdxmid)*dx*dz)
                   fdiy = fdiy -two*xnu*dvdymid*dx*dz
 
               enddo
             enddo
             convx = convx + fcvx
             stressx = stressx + fdix
             convy = convy + fcvy
             stressy = stressy + fdiy
             pressy = pressy + fpry
         endif
 
         !calculating the flux through bottom wall (y-direction)
         if ((jcvlw(iv).ge.xstart(2)).and.(jcvlw(iv).le.xend(2))) then
             j=jcvlw(iv)-xstart(2)+1
             fcvx=zero
             fcvy=zero
             fpry=zero
             fdix=zero
             fdiy=zero
             do k=kcvfr_lx(iv),kcvbk_lx(iv)
               do i=icvlf_lx(iv),icvrt_lx(iv)-1
                   !momentum flux
                   !FIXME avoid interpolation for the non-linear term
                   uxmid = half*(ux1(i,j,k)+ux1(i+1,j,k))
                   uymid = half*(uy1(i,j,k)+uy1(i+1,j,k))
                   fcvx= fcvx +uxmid*uymid*dx*dz
                   fcvy= fcvy +uymid*uymid*dx*dz
 
                   !pressure
                   prmid = half*(ppi1(i,j,k)+ppi1(i+1,j,k))
                   fpry = fpry -prmid*dx*dz
 
                   !viscous term
                   dudymid = half*(td1(i,j,k)+td1(i+1,j,k))
                   dvdxmid = half*(tb1(i,j,k)+tb1(i+1,j,k))
                   dvdymid = half*(te1(i,j,k)+te1(i+1,j,k))
                   fdix = fdix +(xnu*(dudymid+dvdxmid)*dx*dz)
                   fdiy = fdiy +two*xnu*dvdymid*dx*dz
               enddo
             enddo
             convx = convx + fcvx
             stressx = stressx + fdix
             convy = convy + fcvy
             stressy = stressy + fdiy
             pressy = pressy + fpry
         endif
 
         !calculating the flux through the front wall (z-direction)
         if((kcvfr(iv).ge.xstart(3).and.(kcvfr(iv).le.xend(3)))) then
             k = kcvfr(iv)-xstart(3)+1
             fcvx=zero
             fcvy=zero
             fpry=zero
             fdix=zero
             fdiy=zero
             do j=jcvlw_lx(iv), jcvup_lx(iv)
               do i=icvlf_lx(iv), icvrt_lx(iv)-1
                   !momentum flux
                   !FIXME avoid interpolation for the non-linear term
                   uxmid = half*(ux1(i,j,k)+ux1(i+1,j,k))
                   uymid = half*(uy1(i,j,k)+uy1(i+1,j,k))
                   uzmid = half*(uy1(i,j,k)+uy1(i+1,j,k))
                   fcvx= fcvx + uxmid*uzmid*dx*del_y(j)
                   fcvy= fcvy + uymid*uzmid*dx*del_y(j)
 
                   !viscous flux
                   dudzmid = half*(tg1(i,j,k)+tg1(i+1,j,k))
                   dwdxmid = half*(tc1(i,j,k)+tc1(i+1,j,k))
                   dwdymid = half*(tf1(i,j,k)+tf1(i+1,j,k))
                   dvdzmid = half*(th1(i,j,k)+th1(i+1,j,k))
                   fdix = fdix + (xnu*(dudzmid+dwdxmid)*dx*del_y(j))
                   fdiy = fdiy + (xnu*(dvdzmid+dwdymid)*dx*del_y(j))
               enddo
             enddo   
             convx = convx + fcvx
             stressx = stressx + fdix
             convy = convy + fcvy
             stressy = stressy + fdiy
         endif
 
         !calculating the flux through the front wall (z-direction)
         if((kcvbk(iv).ge.xstart(3).and.(kcvbk(iv).le.xend(3)))) then
             k = kcvbk(iv)-xstart(3)+1
             fcvx=zero
             fcvy=zero
             fpry=zero
             fdix=zero
             fdiy=zero
             do j=jcvlw_lx(iv), jcvup_lx(iv)
               do i=icvlf_lx(iv), icvrt_lx(iv)-1
                   !momentum flux
                   !FIXME avoid interpolation for the non-linear term
                   uxmid = half*(ux1(i,j,k)+ux1(i+1,j,k))
                   uymid = half*(uy1(i,j,k)+uy1(i+1,j,k))
                   uzmid = half*(uy1(i,j,k)+uy1(i+1,j,k))
                   fcvx= fcvx - uxmid*uzmid*dx*del_y(j)
                   fcvy= fcvy - uymid*uzmid*dx*del_y(j)
 
                   !viscous flux
                   dudzmid = half*(tg1(i,j,k)+tg1(i+1,j,k))
                   dwdxmid = half*(tc1(i,j,k)+tc1(i+1,j,k))
                   dwdymid = half*(tf1(i,j,k)+tf1(i+1,j,k))
                   dvdzmid = half*(th1(i,j,k)+th1(i+1,j,k))
                   fdix = fdix - (xnu*(dudzmid+dwdxmid)*dx*del_y(j))
                   fdiy = fdiy - (xnu*(dvdzmid+dwdymid)*dx*del_y(j))
               enddo
             enddo   
             convx = convx + fcvx
             stressx = stressx + fdix
             convy = convy + fcvy
             stressy = stressy + fdiy
         endif
 
         call MPI_ALLREDUCE(MPI_IN_PLACE, convx, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
         if (code /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_ALLREDUCE!")
         end if
         call MPI_ALLREDUCE(MPI_IN_PLACE, convy, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
         if (code /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_ALLREDUCE!")
         end if
         call MPI_ALLREDUCE(MPI_IN_PLACE, pressx, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
         if (code /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_ALLREDUCE!")
         end if
         call MPI_ALLREDUCE(MPI_IN_PLACE, pressy, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
         if (code /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_ALLREDUCE!")
         end if
         call MPI_ALLREDUCE(MPI_IN_PLACE, stressx, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
         if (code /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_ALLREDUCE!")
         end if
         call MPI_ALLREDUCE(MPI_IN_PLACE, stressy, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
         if (code /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_ALLREDUCE!")
         end if
   
         pressx = pressx/dt
         pressy = pressy/dt 
         
         xDrag_mean = xmom + convx - pressx -stressx
         yLift_mean = ymom + convy - pressy -stressy
         
         if (nrank .eq. 0) then
             write(38,*) t,xDrag_mean,yLift_mean 
             flush(38)
         endif
       enddo
    endif
    
    if (mod(itime, icheckpoint).eq.0) then
      if (nrank .eq. 0) then
        write(filename, '(A,I7.7,A)') 'forces', itime, '.dat'
        call execute_command_line("cp forces.dat " //filename)
      endif
    endif

    do k = 1, xsize(3)
       do j = 1, xsize(2)
          do i = 1, xsize(1)
             ux11(i,j,k)=ux01(i,j,k)
             uy11(i,j,k)=uy01(i,j,k)
             ux01(i,j,k)=ux1(i,j,k)
             uy01(i,j,k)=uy1(i,j,k)
          enddo
       enddo
    enddo

    return

  end subroutine force

end module forces
