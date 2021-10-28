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
  USE decomp_2d
  implicit none

  integer :: nvol,iforces
  real(mytype),save,allocatable,dimension(:,:,:) :: ux01, uy01, ux11, uy11, ppi1
  real(mytype),allocatable,dimension(:) :: xld,xrd,yld,yud
  integer,allocatable,dimension(:) :: icvlf,icvrt,jcvlw,jcvup
  integer,allocatable,dimension(:) :: icvlf_lx,icvrt_lx,icvlf_ly,icvrt_ly
  integer,allocatable,dimension(:) :: jcvlw_lx,jcvup_lx,jcvlw_ly,jcvup_ly

  character(len=*), parameter :: io_restart_forces = "restart-forces-io", &
       resfile = "restart-forces"
  
contains

  subroutine init_forces

    USE decomp_2d
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

    allocate(icvlf(nvol),icvrt(nvol),jcvlw(nvol),jcvup(nvol))
    allocate(icvlf_lx(nvol),icvrt_lx(nvol),icvlf_ly(nvol),icvrt_ly(nvol))
    allocate(jcvlw_lx(nvol),jcvup_lx(nvol),jcvlw_ly(nvol),jcvup_ly(nvol))

    !     Definition of the Control Volume
    !*****************************************************************
    !! xld,xrd,yld,yud: limits of control volume (!!don't use cex and cey anymore!!)

    do iv=1,nvol
       ! ok for istret=0 (!!to do for istret=1!!)
       icvlf(iv) = nint(xld(iv)/dx)+1
       icvrt(iv) = nint(xrd(iv)/dx)+1
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
       icvlf_ly(iv) = max(icvlf(iv)+1-ystart(1),1)
       icvrt_ly(iv) = min(icvrt(iv)+1-ystart(1),ysize(1))
       jcvlw_ly(iv) = jcvlw(iv)
       jcvup_ly(iv) = jcvup(iv)
    enddo

    if (nrank==0) then
       write(*,*) '========================Forces============================='
       write(*,*) '                       (icvlf)      (icvrt) '
       write(*,*) '                (jcvup) B____________C '
       write(*,*) '                        \            \ '
       write(*,*) '                        \     __     \ '
       write(*,*) '                        \    \__\    \ '
       write(*,*) '                        \            \ '
       write(*,*) '                        \       CV   \ '
       write(*,*) '                (jcvlw) A____________D '
       do iv=1,nvol
          write(*,"(' Control Volume     : #',I1)") iv
          write(*,"('     xld, icvlf     : (',F6.2,',',I6,')')") xld(iv), icvlf(iv)
          write(*,"('     xrd, icvrt     : (',F6.2,',',I6,')')") xrd(iv), icvrt(iv)
          write(*,"('     yld, jcvlw     : (',F6.2,',',I6,')')") yld(iv), jcvlw(iv)
          write(*,"('     yud, jcvup     : (',F6.2,',',I6,')')") yud(iv), jcvup(iv)
       enddo
       write(*,*) '==========================================================='
    endif

    call decomp_2d_init_io(io_restart_forces)
    call decomp_2d_register_variable(io_restart_forces, "ux01", 1, 0, 0, mytype)
    call decomp_2d_register_variable(io_restart_forces, "uy01", 1, 0, 0, mytype)
    call decomp_2d_register_variable(io_restart_forces, "ux11", 1, 0, 0, mytype)
    call decomp_2d_register_variable(io_restart_forces, "uy11", 1, 0, 0, mytype)
    
  end subroutine init_forces

  subroutine restart_forces(itest1)

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none

    integer :: ierror,code,itest1
    integer :: ierror_o=0 !error to open sauve file during restart
    character(len=30) :: filename, filestart


#ifndef ADIOS2
    write(filename, "('restart-forces',I7.7)") itime
    write(filestart,"('restart-forces',I7.7)") ifirst-1
#else
    write(filename, *) "restart-forces"
#endif
    if (itest1 == 1) then
       call decomp_2d_open_io(io_restart_forces, resfile, decomp_2d_write_mode)
    else
       call decomp_2d_open_io(io_restart_forces, resfile, decomp_2d_read_mode)
    endif
    call decomp_2d_start_io(io_restart_forces, resfile)
    
    if (itest1==1) then
       !write
       if (mod(itime, icheckpoint).ne.0) then
          return
       endif

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
    
    if (nrank.eq.0) then
       if (ierror_o .ne. 0) then !Included by Felipe Schuch
          write(*,*) '==========================================================='
          write(*,*) 'Error: Impossible to read '//trim(filestart)
          write(*,*) '==========================================================='
          call MPI_ABORT(MPI_COMM_WORLD,code,ierror)
       endif
    endif

  end subroutine restart_forces
end module forces

!***********************************************************************
subroutine force(ux1,uy1,ep1)

  !***********************************************************************

  USE forces
  USE param
  USE variables
  USE decomp_2d
  USE MPI
  USE ibm_param

  use var, only : ta1, tb1, tc1, td1, di1
  use var, only : ux2, uy2, ta2, tb2, tc2, td2, di2

  implicit none
  character(len=30) :: filename, filename2
  integer :: nzmsize
  integer                                             :: i, iv, j, k, kk, code, jj
  integer                                             :: nvect1,nvect2,nvect3

  real(mytype), dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1, uy1
  real(mytype), dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ep1

  real(mytype), dimension(ysize(1),ysize(2),ysize(3)) :: ppi2

  real(mytype), dimension(nz) :: yLift,xDrag
  real(mytype) :: yLift_mean,xDrag_mean

  real(mytype), dimension(ny-1) :: del_y

  real(mytype), dimension(nz) :: tunstxl, tunstyl
  real(mytype), dimension(nz) :: tconvxl,tconvyl
  real(mytype), dimension(nz) :: tpresxl,tpresyl
  real(mytype), dimension(nz) :: tdiffxl,tdiffyl

  real(mytype), dimension(nz) :: tunstx, tunsty
  real(mytype), dimension(nz) :: tconvx,tconvy
  real(mytype), dimension(nz) :: tpresx,tpresy
  real(mytype), dimension(nz) :: tdiffx,tdiffy

  real(mytype) :: uxmid,uymid,prmid
  real(mytype) :: dudxmid,dudymid,dvdxmid,dvdymid
  real(mytype) :: fac,tsumx,tsumy
  real(mytype) :: fcvx,fcvy,fprx,fpry,fdix,fdiy
  real(mytype) :: xmom,ymom

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
  call transpose_x_to_y(ppi1,ppi2)

  call dery (tc2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx) ! dudy
  call dery (td2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)    ! dvdy
  call transpose_y_to_x(tc2,tc1) ! dudy
  call transpose_y_to_x(td2,td1) ! dvdy

  !*****************************************************************
  !      Drag and Lift coefficients
  !*****************************************************************
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

     tunstxl=zero
     tunstyl=zero
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
        tunstxl(xstart(3)-1+k)=tsumx
        tunstyl(xstart(3)-1+k)=tsumy
     enddo
     call MPI_ALLREDUCE(tunstxl,tunstx,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(tunstyl,tunsty,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)

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

     tconvxl=zero
     tconvyl=zero
     tdiffxl=zero
     tdiffyl=zero
     tpresxl=zero
     tpresyl=zero
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
              uxmid = half*(ux1(i,j,k)+ux1(i+1,j,k))
              uymid = half*(uy1(i,j,k)+uy1(i+1,j,k))
              fcvx  = fcvx -uxmid*uymid*dx
              fcvy  = fcvy -uymid*uymid*dx

              !pressure
              prmid = half*(ppi1(i,j,k)+ppi1(i+1,j,k))
              fpry  = fpry +prmid*dx

              !viscous term
              dudymid = half*(tc1(i,j,k)+tc1(i+1,j,k))
              dvdxmid = half*(tb1(i,j,k)+tb1(i+1,j,k))
              dvdymid = half*(td1(i,j,k)+td1(i+1,j,k))
              fdix    = fdix -(xnu*(dudymid+dvdxmid)*dx)
              fdiy    = fdiy -two*xnu*dvdymid*dx

           enddo

           tconvxl(kk)=tconvxl(kk)+fcvx
           tconvyl(kk)=tconvyl(kk)+fcvy
           tpresyl(kk)=tpresyl(kk)+fpry
           tdiffxl(kk)=tdiffxl(kk)+fdix
           tdiffyl(kk)=tdiffyl(kk)+fdiy
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
              uxmid = half*(ux1(i,j,k)+ux1(i+1,j,k))
              uymid = half*(uy1(i,j,k)+uy1(i+1,j,k))
              fcvx= fcvx +uxmid*uymid*dx
              fcvy= fcvy +uymid*uymid*dx

              !pressure
              prmid = half*(ppi1(i,j,k)+ppi1(i+1,j,k))
              fpry = fpry -prmid*dx

              !viscous term
              dudymid = half*(tc1(i,j,k)+tc1(i+1,j,k))
              dvdxmid = half*(tb1(i,j,k)+tb1(i+1,j,k))
              dvdymid = half*(td1(i,j,k)+td1(i+1,j,k))
              fdix = fdix +(xnu*(dudymid+dvdxmid)*dx)
              fdiy = fdiy +two*xnu*dvdymid*dx

           enddo
           tconvxl(kk)=tconvxl(kk)+fcvx
           tconvyl(kk)=tconvyl(kk)+fcvy
           tpresyl(kk)=tpresyl(kk)+fpry
           tdiffxl(kk)=tdiffxl(kk)+fdix
           tdiffyl(kk)=tdiffyl(kk)+fdiy
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
              uxmid = half*(ux2(i,j,k)+ux2(i,j+1,k))
              uymid = half*(uy2(i,j,k)+uy2(i,j+1,k))
              fcvx= fcvx -uxmid*uxmid*del_y(j)
              fcvy= fcvy -uxmid*uymid*del_y(j)

              !pressure
              prmid=half*(ppi2(i,j,k)+ppi2(i,j+1,k))
              fprx = fprx +prmid*del_y(j)

              !viscous term
              dudxmid = half*(ta2(i,j,k)+ta2(i,j+1,k))
              dudymid = half*(tc2(i,j,k)+tc2(i,j+1,k))
              dvdxmid = half*(tb2(i,j,k)+tb2(i,j+1,k))
              fdix = fdix -two*xnu*dudxmid*del_y(j)
              fdiy = fdiy -xnu*(dvdxmid+dudymid)*del_y(j)
           enddo
           tconvxl(kk)=tconvxl(kk)+fcvx
           tconvyl(kk)=tconvyl(kk)+fcvy
           tpresxl(kk)=tpresxl(kk)+fprx
           tdiffxl(kk)=tdiffxl(kk)+fdix
           tdiffyl(kk)=tdiffyl(kk)+fdiy
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
              uxmid = half*(ux2(i,j,k)+ux2(i,j+1,k))
              uymid = half*(uy2(i,j,k)+uy2(i,j+1,k))
              fcvx= fcvx +uxmid*uxmid*del_y(j)
              fcvy= fcvy +uxmid*uymid*del_y(j)

              !pressure
              prmid=half*(ppi2(i,j,k)+ppi2(i,j+1,k))
              fprx = fprx -prmid*del_y(j)

              !viscous term
              dudxmid = half*(ta2(i,j,k)+ta2(i,j+1,k))
              dudymid = half*(tc2(i,j,k)+tc2(i,j+1,k))
              dvdxmid = half*(tb2(i,j,k)+tb2(i,j+1,k))
              fdix = fdix +two*xnu*dudxmid*del_y(j)
              fdiy = fdiy +xnu*(dvdxmid+dudymid)*del_y(j)
           enddo
           tconvxl(kk)=tconvxl(kk)+fcvx
           tconvyl(kk)=tconvyl(kk)+fcvy
           tpresxl(kk)=tpresxl(kk)+fprx
           tdiffxl(kk)=tdiffxl(kk)+fdix
           tdiffyl(kk)=tdiffyl(kk)+fdiy
        enddo
     endif
     call MPI_ALLREDUCE(tconvxl,tconvx,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(tconvyl,tconvy,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(tpresxl,tpresx,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(tpresyl,tpresy,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(tdiffxl,tdiffx,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(tdiffyl,tdiffy,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)

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
        call flush(38)
     endif
     if (mod(itime, icheckpoint).eq.0) then
        if (nrank .eq. 0) then
           write(filename,"('forces.dat',I7.7)") itime
           call system("cp forces.dat " //filename)
        endif
     endif
  enddo

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
