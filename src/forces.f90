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
  USE decomp_2d
  implicit none

  integer :: nvol,iforces
  real(mytype),save,allocatable,dimension(:,:,:) :: ux01, uy01, ux11, uy11, ppi1, uz01, uz11
  real(mytype),allocatable,dimension(:) :: xld, xrd, yld, yud, zld, zrd, xld2, xrd2, yld2, yud2, zld2, zrd2
  integer,allocatable,dimension(:) :: icvlf,icvrt,jcvlw,jcvup,zcvlf,zcvrt
  integer,allocatable,dimension(:) :: icvlf_lx, icvrt_lx, icvlf_ly, icvrt_ly, icvlf_lz, icvrt_lz
  integer,allocatable,dimension(:) :: jcvlw_lx, jcvup_lx, jcvlw_ly, jcvup_ly, jcvlw_lz, jcvup_lz
  integer,allocatable,dimension(:) :: zcvlf_lx, zcvrt_lx, zcvlf_ly, zcvrt_ly, zcvlf_lz, zcvrt_lz
  

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

   !  write(*,*) 'Inside INIT_FORCES'

    call alloc_x(ux01)
    call alloc_x(uy01)
    call alloc_x(ux11)
    call alloc_x(uy11)
    call alloc_x(ppi1)
    call alloc_x(uz01)
    call alloc_x(uz11)

    ux01 = zero
    uy01 = zero
    ux11 = zero
    uy11 = zero
    uz01 = zero
    uz11 = zero    

   !  write(*,*) 'Alloc_x called'

    allocate(icvlf(nvol), icvrt(nvol), jcvlw(nvol), jcvup(nvol), zcvlf(nvol), zcvrt(nvol))
    allocate(icvlf_lx(nvol), icvrt_lx(nvol), icvlf_ly(nvol), icvrt_ly(nvol), icvlf_lz(nvol), icvrt_lz(nvol))
    allocate(jcvlw_lx(nvol), jcvup_lx(nvol), jcvlw_ly(nvol), jcvup_ly(nvol), jcvlw_lz(nvol), jcvup_lz(nvol))
    allocate(zcvlf_lx(nvol), zcvrt_lx(nvol), zcvlf_ly(nvol), zcvrt_ly(nvol), zcvlf_lz(nvol), zcvrt_lz(nvol))
    allocate(xld2(nvol), xrd2(nvol), yld2(nvol), yud2(nvol), zld2(nvol), zrd2(nvol))

   !  write(*,*) 'allocate called'

   !  if ((iibm.ne.0).and.(t.ne.0.)) then
   !    xld2(:) = xld(:) + (t-ifirst*dt)*ubcx
   !    xrd2(:) = xrd(:) + (t-ifirst*dt)*ubcx
   !    yld2(:) = yld(:) + (t-ifirst*dt)*ubcy
   !    yud2(:) = yud(:) + (t-ifirst*dt)*ubcy
   ! else
      xld2(:) = xld(:)
      xrd2(:) = xrd(:)
      yld2(:) = yld(:)
      yud2(:) = yud(:)
      zld2(:) = zld(:)
      zrd2(:) = zrd(:)
   ! endif
  
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
       jcvlw_lz(iv) = max(jcvlw(iv)+1-zstart(2),1)
       jcvup_lz(iv) = min(jcvup(iv)+1-zstart(2),zsize(2))       
  
       icvlf_ly(iv) = max(icvlf(iv)+1-ystart(1),1)
       icvrt_ly(iv) = min(icvrt(iv)+1-ystart(1),ysize(1))
       icvlf_lz(iv) = max(icvlf(iv)+1-zstart(1),1)
       icvrt_lz(iv) = min(icvrt(iv)+1-zstart(1),zsize(1))   
       jcvlw_ly(iv) = jcvlw(iv)
       jcvup_ly(iv) = jcvup(iv)

       zcvlf(iv) = nint(zld(iv)/dz)+1
       zcvrt(iv) = nint(zrd(iv)/dz)+1
       zcvlf_lx(iv) = max(zcvlf(iv)+1-xstart(3),1)
       zcvrt_lx(iv) = min(zcvrt(iv)+1-xstart(3),xsize(3)) 
       zcvlf_ly(iv) = max(zcvlf(iv)+1-ystart(3),1)
       zcvrt_ly(iv) = min(zcvrt(iv)+1-ystart(3),ysize(3))       
       zcvlf_lz(iv) = zcvlf(iv)
       zcvrt_lz(iv) = zcvrt(iv) 
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
          write(*,"('     zld, zcvlf     : (',F6.2,',',I6,')')") zld(iv), zcvlf(iv)
          write(*,"('     zrd, zcvrt     : (',F6.2,',',I6,')')") zrd(iv), zcvrt(iv)      
       enddo
       write(*,*) '==========================================================='
    endif

    call decomp_2d_init_io(io_restart_forces)
    call decomp_2d_register_variable(io_restart_forces, "ux01", 1, 0, 0, mytype)
    call decomp_2d_register_variable(io_restart_forces, "uy01", 1, 0, 0, mytype)
    call decomp_2d_register_variable(io_restart_forces, "ux11", 1, 0, 0, mytype)
    call decomp_2d_register_variable(io_restart_forces, "uy11", 1, 0, 0, mytype)
    call decomp_2d_register_variable(io_restart_forces, "uz01", 1, 0, 0, mytype)
    call decomp_2d_register_variable(io_restart_forces, "uz11", 1, 0, 0, mytype)

  end subroutine init_forces

  subroutine update_forces

   USE decomp_2d
   USE decomp_2d_io, only : decomp_2d_register_variable, decomp_2d_init_io
   USE param
   USE variables
   implicit none

   integer :: iv,stp1,stp2,h

!   !  write(*,*) 'Inside INIT_FORCES'

!    call alloc_x(ux01)
!    call alloc_x(uy01)
!    call alloc_x(ux11)
!    call alloc_x(uy11)
!    call alloc_x(ppi1)
!    call alloc_x(uz01)
!    call alloc_x(uz11)

!    ux01 = zero
!    uy01 = zero
!    ux11 = zero
!    uy11 = zero
!    uz01 = zero
!    uz11 = zero    

!   !  write(*,*) 'Alloc_x called'

!    allocate(icvlf(nvol), icvrt(nvol), jcvlw(nvol), jcvup(nvol), zcvlf(nvol), zcvrt(nvol))
!    allocate(icvlf_lx(nvol), icvrt_lx(nvol), icvlf_ly(nvol), icvrt_ly(nvol), icvlf_lz(nvol), icvrt_lz(nvol))
!    allocate(jcvlw_lx(nvol), jcvup_lx(nvol), jcvlw_ly(nvol), jcvup_ly(nvol), jcvlw_lz(nvol), jcvup_lz(nvol))
!    allocate(zcvlf_lx(nvol), zcvrt_lx(nvol), zcvlf_ly(nvol), zcvrt_ly(nvol))
!    allocate(xld2(nvol), xrd2(nvol), yld2(nvol), yud2(nvol), zld2(nvol), zrd2(nvol))

  !  write(*,*) 'allocate called'

  !  if ((iibm.ne.0).and.(t.ne.0.)) then
  !    xld2(:) = xld(:) + (t-ifirst*dt)*ubcx
  !    xrd2(:) = xrd(:) + (t-ifirst*dt)*ubcx
  !    yld2(:) = yld(:) + (t-ifirst*dt)*ubcy
  !    yud2(:) = yud(:) + (t-ifirst*dt)*ubcy
  ! else
     xld2(:) = xld(:)
     xrd2(:) = xrd(:)
     yld2(:) = yld(:)
     yud2(:) = yud(:)
     zld2(:) = zld(:)
     zrd2(:) = zrd(:)
  ! endif
 
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
      jcvlw_lz(iv) = max(jcvlw(iv)+1-zstart(2),1)
      jcvup_lz(iv) = min(jcvup(iv)+1-zstart(2),zsize(2))       
 
      icvlf_ly(iv) = max(icvlf(iv)+1-ystart(1),1)
      icvrt_ly(iv) = min(icvrt(iv)+1-ystart(1),ysize(1))
      icvlf_lz(iv) = max(icvlf(iv)+1-zstart(1),1)
      icvrt_lz(iv) = min(icvrt(iv)+1-zstart(1),zsize(1))   
      jcvlw_ly(iv) = jcvlw(iv)
      jcvup_ly(iv) = jcvup(iv)

      zcvlf(iv) = nint(zld(iv)/dz)+1
      zcvrt(iv) = nint(zrd(iv)/dz)+1
      zcvlf_lx(iv) = max(zcvlf(iv)+1-xstart(3),1)
      zcvrt_lx(iv) = min(zcvrt(iv)+1-xstart(3),xsize(3)) 
      zcvlf_ly(iv) = max(zcvlf(iv)+1-ystart(3),1)
      zcvrt_ly(iv) = min(zcvrt(iv)+1-ystart(3),ysize(3))        
   enddo

   ! if (nrank==0) then
   !    write(*,*) '========================Forces============================='
   !    write(*,*) '                       (icvlf)      (icvrt) '
   !    write(*,*) '                (jcvup) B____________C '
   !    write(*,*) '                        \            \ '
   !    write(*,*) '                        \     __     \ '
   !    write(*,*) '                        \    \__\    \ '
   !    write(*,*) '                        \            \ '
   !    write(*,*) '                        \       CV   \ '
   !    write(*,*) '                (jcvlw) A____________D '
   !    do iv=1,nvol
   !       write(*,"(' Control Volume     : #',I1)") iv
   !       write(*,"('     xld, icvlf     : (',F6.2,',',I6,')')") xld(iv), icvlf(iv)
   !       write(*,"('     xrd, icvrt     : (',F6.2,',',I6,')')") xrd(iv), icvrt(iv)
   !       write(*,"('     yld, jcvlw     : (',F6.2,',',I6,')')") yld(iv), jcvlw(iv)
   !       write(*,"('     yud, jcvup     : (',F6.2,',',I6,')')") yud(iv), jcvup(iv)
   !       write(*,"('     zld, zcvlf     : (',F6.2,',',I6,')')") zld(iv), zcvlf(iv)
   !       write(*,"('     zrd, zcvrt     : (',F6.2,',',I6,')')") zrd(iv), zcvrt(iv)      
   !    enddo
   !    write(*,*) '==========================================================='
   ! endif

   ! call decomp_2d_init_io(io_restart_forces)
   ! call decomp_2d_register_variable(io_restart_forces, "ux01", 1, 0, 0, mytype)
   ! call decomp_2d_register_variable(io_restart_forces, "uy01", 1, 0, 0, mytype)
   ! call decomp_2d_register_variable(io_restart_forces, "ux11", 1, 0, 0, mytype)
   ! call decomp_2d_register_variable(io_restart_forces, "uy11", 1, 0, 0, mytype)
   ! call decomp_2d_register_variable(io_restart_forces, "uz01", 1, 0, 0, mytype)
   ! call decomp_2d_register_variable(io_restart_forces, "uz11", 1, 0, 0, mytype)


  end subroutine update_forces
!   if ((iibm.ne.0).and.(t.ne.0.)) then
   !    xld2(:) = xld(:) + (t-ifirst*dt)*ubcx
   !    xrd2(:) = xrd(:) + (t-ifirst*dt)*ubcx
   !    yld2(:) = yld(:) + (t-ifirst*dt)*ubcy
   !    yud2(:) = yud(:) + (t-ifirst*dt)*ubcy
   ! else
   !    xld2(:) = xld(:)
   !    xrd2(:) = xrd(:)
   !    yld2(:) = yld(:)
   !    yud2(:) = yud(:)
   ! endif
  subroutine restart_forces(itest1)

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none

    integer :: ierror,code,itest1
    integer :: ierror_o=0 !error to open save file during restart
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

  subroutine force(ux1,uy1,uz1,ep1,dra1,dra2,dra3,record_var)

    USE param
    USE variables
    USE decomp_2d
    USE MPI
    USE ibm_param
    USE ellipsoid_utils, only : CrossProduct,centrifugal_force,coriolis_force

    use var, only : ta1, tb1, tc1, td1, te1, tf1, tg1, th1, ti1, di1
    use var, only : ux2, uy2, uz2, ta2, tb2, tc2, td2, te2, tf2, tg2, th2, ti2, di2
    use var, only : ux3, uy3, uz3, ta3, tb3, tc3, td3, te3, tf3, tg3, th3, ti3, di3
      
  
    implicit none
    character(len=30) :: filename, filename2
    integer :: nzmsize
    integer                                             :: i, iv, j, k, kk, code, jj, ii
    integer                                             :: nvect1,nvect2,nvect3

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1, uy1, uz1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ep1
    integer, intent(in) ::record_var
    real(mytype), intent(out)                                       :: dra1,dra2,dra3

    real(mytype), dimension(ysize(1),ysize(2),ysize(3)) :: ppi2
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: ppi3

    real(mytype), dimension(nz) :: yLift,xDrag, zLat
    real(mytype) :: yLift_mean,xDrag_mean,zLat_mean
    real(mytype) :: xm,ym,zm,rotationalComponent(3)

    real(mytype), dimension(ny-1) :: del_y

    real(mytype), dimension(nz) :: tunstxl, tunstyl, tunstzl
    real(mytype), dimension(nz) :: tconvxl,tconvyl,tconvzl
    real(mytype), dimension(nz) :: tpresxl,tpresyl
    real(mytype), dimension(nz) :: tdiffxl,tdiffyl,tdiffzl

    real(mytype), dimension(nz) :: tunstx, tunsty, tunstz
    real(mytype), dimension(nz) :: tconvx,tconvy,tconvz
    real(mytype), dimension(nz) :: tpresx,tpresy
    real(mytype), dimension(nz) :: tdiffx,tdiffy,tdiffz

        
    
    real(mytype), dimension(ny) :: tconvxl2, tconvyl2, tconvzl2
    real(mytype), dimension(ny) :: tdiffxl2, tdiffyl2, tdiffzl2
    real(mytype), dimension(ny) :: tconvx2, tconvy2, tconvz2
    real(mytype), dimension(ny) :: tdiffx2, tdiffy2, tdiffz2
    real(mytype), dimension(ny) :: tpreszl, tpresz
    
    
    real(mytype) :: uxmid,uymid,uzmid,prmid
    real(mytype) :: dudxmid,dudymid,dudzmid,dvdxmid,dvdymid,dvdzmid
    real(mytype) :: dwdxmid,dwdymid,dwdzmid
    real(mytype) :: fac,fac1,fac2,fac3,tsumx,tsumy,tsumz,centrifugal(3),coriolis(3)
    real(mytype) :: fcvx,fcvy,fcvz,fprx,fpry,fprz,fdix,fdiy,fdiz
    real(mytype) :: xmom,ymom,zmom
    real(mytype), dimension(ny) :: ztpresx, ztpresy
    real(mytype), dimension(nz) :: zyLift, zxDrag, zzLat
    real(mytype) :: zyLift_mean, zxDrag_mean, zzLat_mean


    real(mytype), dimension(nz) :: drag1, drag2, drag11, drag22
    real(mytype), dimension(nz) :: drag3, drag4, drag33, drag44
    real(mytype) :: mom1, mom2, mom3, tp1, tp2, tp3

   !  write(*,*) 'Inside FORCE'

  

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
      do iv=1,nvol
         ! if (nrank.eq.0) then 
         !    open(12,file="Body1.dat",status='unknown',form='formatted')
         ! endif
         if ((nrank .eq. 0).and.(record_var.eq.1)) then
            write(filename,"('forces.dat',I1.1)") iv
            open(38+(iv-1),file=filename,status='unknown',form='formatted')
            ! write(*,*) 'Opened file: ', filename, 'number = ', 38+(iv-1)
         endif
      enddo
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                ux11(i,j,k)=ux1(i,j,k)
                uy11(i,j,k)=uy1(i,j,k)
                uz11(i,j,k)=uz1(i,j,k)
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
                uz01(i,j,k)=uz1(i,j,k)
             enddo
          enddo
       enddo
       return
    endif

    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,1)    ! dudx !x is 1
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,2) ! dvdx !y is 2
    call derx (te1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,3) ! dw/dx!z is 3

    call transpose_x_to_y(ta1,ta2) ! dudx
    call transpose_x_to_y(tb1,tb2) ! dvdx
    call transpose_x_to_y(te1,te2) ! dw/dx

    call transpose_x_to_y(ux1,ux2)
    call transpose_x_to_y(uy1,uy2)
    call transpose_x_to_y(uz1,uz2)
    call transpose_x_to_y(ppi1,ppi2)

    call dery (tc2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,1) ! dudy !x is 1
    call dery (td2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,2)    ! dvdy !y is 2
    call dery (tf2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,3) ! dw/dy!z is 3
    

    call transpose_y_to_z(ux2,ux3)
    call transpose_y_to_z(uy2,uy3)
    call transpose_y_to_z(uz2,uz3)

    call derz (tg3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,1)  ! du/dz
    call derz (th3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,2)  ! dv/dz
    call derz (ti3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,3)     ! dw/dz
  
    call transpose_z_to_y(tg3,tg2) ! du/dz
    call transpose_z_to_y(th3,th2) ! dv/dz
    call transpose_z_to_y(ti3,ti2) ! 

    call transpose_y_to_x(tc2,tc1) ! dudy
    call transpose_y_to_x(td2,td1) ! dvdy
    call transpose_y_to_x(th2,th1) ! dv/dz
    call transpose_y_to_x(tf2,tf1) ! dw/dy
    call transpose_y_to_x(tg2,tg1) ! 
    call transpose_y_to_x(ti2,ti1) ! 
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
       tunstzl=zero
       do k=zcvlf_lx(iv),zcvrt_lx(iv)
          tsumx=zero
          tsumy=zero
          tsumz=zero
          zm=real(xstart(3)+k-1,mytype)*dz
          do j=jcvlw_lx(iv),jcvup_lx(iv)
            ym=real(xstart(2)+j-1,mytype)*dy
             do i=icvlf_lx(iv),icvrt_lx(iv)
               xm=real(xstart(1)+i-1,mytype)*dx

               fac1   = (onepfive*ux1(i,j,k)-two*ux01(i,j,k)+half*ux11(i,j,k))*(one-ep1(i,j,k))
               fac2   = (onepfive*uy1(i,j,k)-two*uy01(i,j,k)+half*uy11(i,j,k))*(one-ep1(i,j,k))
               fac3   = (onepfive*uz1(i,j,k)-two*uz01(i,j,k)+half*uz11(i,j,k))*(one-ep1(i,j,k))

                call coriolis_force(angularVelocity(2:4),[fac1,fac2,fac3],coriolis)
                call centrifugal_force(angularVelocity(2:4), [xm,ym,zm]-position,centrifugal)
                !     The velocity time rate has to be relative to the cell center, 
                !     and not to the nodes, because, here, we have an integral 
                !     relative to the volume, and, therefore, this has a sense 
                !     of a "source".
                !         fac   = (1.5*ux1(i,j,k)-2.0*ux01(i,j,k)+0.5*ux11(i,j,k))*epcv1(i,j,k)
                !  tsumx = tsumx+(fac1-coriolis(1)-centrifugal(1))*dx*del_y(j+(xstart(2)-1))*dz/dt    !tsumx+fac*dx*dy/dt
                tsumx = tsumx+fac1*dx*del_y(min(ny-1,j+xstart(2)-1))*dz/dt
                !sumx(k) = sumx(k)+dudt1*dx*dy

                !         fac   = (1.5*uy1(i,j,k)-2.0*uy01(i,j,k)+0.5*uy11(i,j,k))*epcv1(i,j,k)
               !  tsumy = tsumy+(fac2-coriolis(2)-centrifugal(2))*dx*del_y(j+(xstart(2)-1))*dz/dt !tsumy+fac*dx*dy/dt
                tsumy = tsumy+fac2*dx*del_y(min(ny-1,j+xstart(2)-1))*dz/dt
                !sumy(k) = sumy(k)+dudt1*dx*dy

               !  tsumz = tsumz+(fac3-coriolis(3)-centrifugal(3))*dx*del_y(j+(xstart(2)-1))*dz/dt     
                tsumz = tsumz+fac3*dx*del_y(min(ny-1,j+xstart(2)-1))*dz/dt
             enddo
          enddo
          tunstxl(xstart(3)-1+k)=tsumx
          tunstyl(xstart(3)-1+k)=tsumy
          tunstzl(xstart(3)-1+k)=tsumz
       enddo
       call MPI_ALLREDUCE(tunstxl,tunstx,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       call MPI_ALLREDUCE(tunstyl,tunsty,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       call MPI_ALLREDUCE(tunstzl,tunstz,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)

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
       
  drag1(:)=0.
  drag2(:)=0.
  drag3(:)=0.
  drag4(:)=0.
  
  drag11(:)=0.
  drag22(:)=0.
  drag33(:)=0.
  drag44(:)=0.

       tconvxl=zero
       tconvyl=zero
       tconvzl=zero
       tdiffxl=zero
       tdiffyl=zero
       tdiffzl=zero
       tpresxl=zero
       tpresyl=zero
       tpreszl=zero

       tconvxl2=zero
       tconvyl2=zero
       tconvzl2=zero
       tdiffxl2=zero
       tdiffyl2=zero
       tdiffzl2=zero  
       !BC and AD : x-pencils
       !AD
       if ((jcvlw(iv).ge.xstart(2)).and.(jcvlw(iv).le.xend(2))) then
          j=jcvlw(iv)-xstart(2)+1
          jj=jcvlw(iv)
          ym=real(jj,mytype)*dy
          do k=zcvlf_lx(iv),zcvrt_lx(iv)
             kk=xstart(3)-1+k
             zm=real(kk,mytype)*dz
             fcvx=zero
             fcvy=zero
             fcvz=zero
             fpry=zero
             fdix=zero
             fdiy=zero
             fdiz=zero
             do i=icvlf_lx(iv),icvrt_lx(iv)-1
                
                ii=xstart(1)+i-1
                xm=real(ii,mytype)*dx
               !  write(*,*) 'Calculating force at upper y boundary', [xm,ym,zm]

                !momentum flux
                call crossProduct(angularVelocity(2:4),[xm,ym,zm]-position,rotationalComponent)
                uxmid = half*(ux1(i,j,k)+ux1(i+1,j,k)) - linearVelocity(1) + rotationalComponent(1)
                uymid = half*(uy1(i,j,k)+uy1(i+1,j,k)) - linearVelocity(2) + rotationalComponent(2)
                uzmid = half*(uz1(i,j,k)+uz1(i+1,j,k)) - linearVelocity(3) + rotationalComponent(3)

                fcvx  = fcvx -uxmid*uymid*dx*dz
                fcvy  = fcvy -uymid*uymid*dx*dz
                fcvz  = fcvz -uymid*uzmid*dx*dz


                !pressure
                prmid = half*(ppi1(i,j,k)+ppi1(i+1,j,k))
                fpry  = fpry +prmid*dx*dz

                !viscous term
                dudymid = half*(tc1(i,j,k)+tc1(i+1,j,k))
                dvdxmid = half*(tb1(i,j,k)+tb1(i+1,j,k))
                dvdymid = half*(td1(i,j,k)+td1(i+1,j,k))
                dwdymid = half*(tf1(i,j,k)+tf1(i+1,j,k))
                dvdzmid = half*(th1(i,j,k)+th1(i+1,j,k))

                fdix    = fdix -(xnu*(dudymid+dvdxmid)*dx*dz)
                fdiy    = fdiy -two*xnu*dvdymid*dx*dz
                fdiz    = fdiz -(xnu*(dwdymid+dvdzmid)*dx*dz)


             enddo

             tconvxl(kk)=tconvxl(kk)+fcvx
             tconvyl(kk)=tconvyl(kk)+fcvy
             tconvzl(kk)=tconvzl(kk)+fcvz
             tpresyl(kk)=tpresyl(kk)+fpry
             tdiffxl(kk)=tdiffxl(kk)+fdix
             tdiffyl(kk)=tdiffyl(kk)+fdiy
             tdiffzl(kk)=tdiffzl(kk)+fdiz
          enddo
       endif
       !BC
       if ((jcvup(iv).ge.xstart(2)).and.(jcvup(iv).le.xend(2))) then
          j=jcvup(iv)-xstart(2)+1
          jj=jcvup(iv)
          ym=real(jj,mytype)*dy
          do k=zcvlf_lx(iv),zcvrt_lx(iv)
             kk=xstart(3)-1+k   
             zm=real(kk,mytype)*dz
             fcvx=zero
             fcvy=zero
             fcvz=zero
             fpry=zero
             fdix=zero
             fdiy=zero
             fdiz=zero
             do i=icvlf_lx(iv),icvrt_lx(iv)-1
               ii=xstart(1)+i-1
               xm=real(ii,mytype)*dx
               ! write(*,*) 'xm = ', xm
               ! write(*,*) 'Calculating force at lower y boundary', [xm,ym,zm]

               !momentum flux
               call crossProduct(angularVelocity(2:4),[xm,ym,zm]-position,rotationalComponent)
               uxmid = half*(ux1(i,j,k)+ux1(i+1,j,k)) - linearVelocity(1) + rotationalComponent(1)
               uymid = half*(uy1(i,j,k)+uy1(i+1,j,k)) - linearVelocity(2) + rotationalComponent(2)
               uzmid = half*(uz1(i,j,k)+uz1(i+1,j,k)) - linearVelocity(3) + rotationalComponent(3)

                fcvx = fcvx +uxmid*uymid*dx*dz
                fcvy = fcvy +uymid*uymid*dx*dz
                fcvz = fcvz +uymid*uzmid*dx*dz


                !pressure
                prmid = half*(ppi1(i,j,k)+ppi1(i+1,j,k))
                fpry = fpry -prmid*dx*dz

                !viscous term
                dudymid = half*(tc1(i,j,k)+tc1(i+1,j,k))
                dvdxmid = half*(tb1(i,j,k)+tb1(i+1,j,k))
                dvdymid = half*(td1(i,j,k)+td1(i+1,j,k))
                dwdymid = half*(tf1(i,j,k)+tf1(i+1,j,k))
                dvdzmid = half*(th1(i,j,k)+th1(i+1,j,k))
  
                fdix = fdix + (xnu*(dudymid+dvdxmid)*dx*dz)
                fdiy = fdiy + two*xnu*dvdymid*dx*dz
                fdiz = fdiz + (xnu*(dwdymid+dvdzmid)*dx*dz)

             enddo
             tconvxl(kk)=tconvxl(kk)+fcvx
             tconvyl(kk)=tconvyl(kk)+fcvy
             tconvzl(kk)=tconvzl(kk)+fcvz

             tpresyl(kk)=tpresyl(kk)+fpry
             tdiffxl(kk)=tdiffxl(kk)+fdix
             tdiffyl(kk)=tdiffyl(kk)+fdiy
             tdiffzl(kk)=tdiffzl(kk)+fdiz

          enddo
       endif
       !AB and DC : y-pencils
       !AB
       if ((icvlf(iv).ge.ystart(1)).and.(icvlf(iv).le.yend(1))) then
          i=icvlf(iv)-ystart(1)+1
          ii=icvlf(iv)-1
          xm=real(ii,mytype)*dx
          do k=zcvlf_ly(iv),zcvrt_ly(iv)
             kk=ystart(3)+k-1
             zm=real(kk,mytype)*dz
             fcvx=zero
             fcvy=zero
             fcvz=zero
             fprx=zero
             fdix=zero
             fdiy=zero
             fdiz=zero
             do j=jcvlw_ly(iv),jcvup_ly(iv)-1

                jj=ystart(2)+j-1
                ym=real(jj,mytype)*dz
               !  write(*,*) 'Calculating force at left x boundary', [xm,ym,zm]
                !momentum flux
                call crossProduct(angularVelocity(2:4),[xm,ym,zm]-position,rotationalComponent)
                uxmid = half*(ux2(i,j,k)+ux2(i,j+1,k)) - linearVelocity(1) + rotationalComponent(1)
                uymid = half*(uy2(i,j,k)+uy2(i,j+1,k)) - linearVelocity(2) + rotationalComponent(2)
                uzmid = half*(uz2(i,j,k)+uz2(i,j+1,k)) - linearVelocity(3) + rotationalComponent(3)


                fcvx = fcvx -uxmid*uxmid*del_y(j)*dz
                fcvy = fcvy -uxmid*uymid*del_y(j)*dz
                fcvz = fcvz -uxmid*uzmid*del_y(j)*dz


                !pressure
                prmid = half*(ppi2(i,j,k)+ppi2(i,j+1,k))
                fprx = fprx +prmid*del_y(j)*dz

                !viscous term
                dudxmid = half*(ta2(i,j,k)+ta2(i,j+1,k))
                dudymid = half*(tc2(i,j,k)+tc2(i,j+1,k))
                dvdxmid = half*(tb2(i,j,k)+tb2(i,j+1,k))
                dwdxmid = half*(te2(i,j,k)+te2(i,j+1,k))
                dudzmid = half*(tg2(i,j,k)+tg2(i,j+1,k))

                fdix = fdix -two*xnu*dudxmid*del_y(j)*dz
                fdiy = fdiy -xnu*(dvdxmid+dudymid)*del_y(j)*dz
                fdiz = fdiz -xnu*(dwdxmid+dudzmid)*del_y(j)*dz
             enddo
             tconvxl(kk)=tconvxl(kk)+fcvx
             tconvyl(kk)=tconvyl(kk)+fcvy
             tconvzl(kk)=tconvzl(kk)+fcvz

             tpresxl(kk)=tpresxl(kk)+fprx
             tdiffxl(kk)=tdiffxl(kk)+fdix
             tdiffyl(kk)=tdiffyl(kk)+fdiy
             tdiffzl(kk)=tdiffzl(kk)+fdiz

          enddo
       endif
       !DC
       if ((icvrt(iv).ge.ystart(1)).and.(icvrt(iv).le.yend(1))) then
          i=icvrt(iv)-ystart(1)+1
          ii=icvrt(iv)
          xm=real(ii,mytype)*dx
          do k=zcvlf_ly(iv),zcvrt_ly(iv)
             kk=ystart(3)-1+k
             zm=real(kk,mytype)*dz
             fcvx=zero
             fcvy=zero
             fcvz=zero
             fprx=zero
             fdix=zero
             fdiy=zero
             fdiz=zero
             do j=jcvlw_ly(iv),jcvup_ly(iv)-1 !!!What's going on here?
                jj=ystart(2)+j-1
                ym=real(jj,mytype)*dy
               !  write(*,*) 'Calculating force at right x boundary', [xm,ym,zm]

                !momentum flux
                call crossProduct(angularVelocity(2:4),[xm,ym,zm]-position,rotationalComponent)
                uxmid = half*(ux2(i,j,k)+ux2(i,j+1,k)) - linearVelocity(1) + rotationalComponent(1)
                uymid = half*(uy2(i,j,k)+uy2(i,j+1,k)) - linearVelocity(2) + rotationalComponent(2)
                uzmid = half*(uz2(i,j,k)+uz2(i,j+1,k)) - linearVelocity(3) + rotationalComponent(3)


                fcvx = fcvx + uxmid*uxmid*del_y(j)*dz
                fcvy = fcvy + uxmid*uymid*del_y(j)*dz
                fcvz = fcvz + uxmid*uzmid*del_y(j)*dz


                !pressure
                prmid = half*(ppi2(i,j,k)+ppi2(i,j+1,k))
                fprx = fprx -prmid*del_y(j)*dz

                !viscous term
                dudxmid = half*(ta2(i,j,k)+ta2(i,j+1,k))
                dudymid = half*(tc2(i,j,k)+tc2(i,j+1,k))
                dvdxmid = half*(tb2(i,j,k)+tb2(i,j+1,k))
                dwdxmid = half*(te2(i,j,k)+te2(i,j+1,k))
                dudzmid = half*(tg2(i,j,k)+tg2(i,j+1,k))

                fdix = fdix + two*xnu*dudxmid*del_y(j)*dz
                fdiy = fdiy + xnu*(dvdxmid+dudymid)*del_y(j)*dz
                fdiz = fdiz + xnu*(dwdxmid+dudzmid)*del_y(j)*dz

             enddo
             tconvxl(kk)=tconvxl(kk)+fcvx
             tconvyl(kk)=tconvyl(kk)+fcvy
             tconvzl(kk)=tconvzl(kk)+fcvz

             tpresxl(kk)=tpresxl(kk)+fprx
             tdiffxl(kk)=tdiffxl(kk)+fdix
             tdiffyl(kk)=tdiffyl(kk)+fdiy
             tdiffzl(kk)=tdiffzl(kk)+fdiz

          enddo
       endif

       !Left & Right : 
       !Left
       if ((zcvlf(iv).ge.xstart(3)).and.(zcvlf(iv).le.xend(3))) then
         k=zcvlf(iv)-xstart(3)+1
         kk=zcvlf(iv)
         zm=real(kk,mytype)*dz
 
         fcvx=zero
         fcvy=zero
         fcvz=zero
         fprz=zero
         fdix=zero
         fdiy=zero
         fdiz=zero
         do j=jcvlw_lx(iv),jcvup_lx(iv)
          kk = xstart(2)-1+j
          jj = xstart(2)-1+j

          ym=real(jj,mytype)*dy
            do i=icvlf_lx(iv),icvrt_lx(iv)-1
               ii=xstart(1)+i-1
               xm=real(ii,mytype)*dx
               ! write(*,*) 'Calculating force at left z boundary', [xm,ym,zm]

               !momentum flux
               call crossProduct(angularVelocity(2:4),[xm,ym,zm]-position,rotationalComponent)
               uxmid = half*(ux1(i,j,k)+ux1(i+1,j,k)) - linearVelocity(1) + rotationalComponent(1)
               uymid = half*(uy1(i,j,k)+uy1(i+1,j,k)) - linearVelocity(2) + rotationalComponent(2)
               uzmid = half*(uz1(i,j,k)+uz1(i+1,j,k)) - linearVelocity(3) + rotationalComponent(3)
 
               fcvx= fcvx +uxmid*uzmid*dx*dy
               fcvy= fcvy +uymid*uzmid*dx*dy
               fcvz= fcvz +uzmid*uzmid*dx*dy
 
               !pressure
               prmid = half*(ppi1(i,j,k)+ppi1(i+1,j,k))
               fprz = fprz -prmid*dx*dy
 
               !viscous term
               dudzmid = half*(tg1(i,j,k)+tg1(i+1,j,k))
               dwdxmid = half*(te1(i,j,k)+te1(i+1,j,k))
               dvdzmid = half*(th1(i,j,k)+th1(i+1,j,k))
               dwdymid = half*(tf1(i,j,k)+tf1(i+1,j,k))
               dwdzmid = half*(ti1(i,j,k)+ti1(i+1,j,k))
                                                                   
               fdix = fdix +(xnu*(dudzmid+dwdxmid)*dx*dy)
               fdiy = fdiy +(xnu*(dvdzmid+dwdymid)*dx*dy)
               fdiz = fdiz +two*xnu*dwdzmid*dx*dy
            enddo
         enddo
 !print*, kk
 !        drag3(kk)=drag3(kk)+fcvx   ! Should be size ny
 !        print*, drag3(kk)
         tconvxl2(kk)=tconvxl2(kk)+fcvx
         tconvyl2(kk)=tconvyl2(kk)+fcvy
         tconvzl2(kk)=tconvzl2(kk)+fcvz
         tpreszl(kk) =tpreszl(kk) +fprz
         tdiffxl2(kk)=tdiffxl2(kk)+fdix
         tdiffyl2(kk)=tdiffyl2(kk)+fdiy
         tdiffzl2(kk)=tdiffzl2(kk)+fdiz        
      endif 
      !Right
      if ((zcvrt(iv).ge.xstart(3)).and.(zcvrt(iv).le.xend(3))) then
         k=zcvrt(iv)-xstart(3)+1
         kk=zcvrt(iv)
         zm=real(kk,mytype)*dz
 !        kk=nrank+1
 
         fcvx=zero
         fcvy=zero
         fcvz=zero
         fprz=zero
         fdix=zero
         fdiy=zero
         fdiz=zero
 !        do k=1,xsize(3)
         do j=jcvlw_lx(iv),jcvup_lx(iv)
         !  kk = xstart(2)-1+j 
          jj = xstart(2)-1+j
          ym=real(jj,mytype)*dy
           do i=icvlf_lx(iv),icvrt_lx(iv)-1
             ii=xstart(1)+i-1
             xm=real(ii,mytype)*dx
               !momentum flux
             call crossProduct(angularVelocity(2:4),[xm,ym,zm]-position,rotationalComponent)
            !  write(*,*) 'Calculating force at right z boundary', [xm,ym,zm]

             uxmid = half*(ux1(i,j,k)+ux1(i+1,j,k)) - linearVelocity(1) + rotationalComponent(1)
             uymid = half*(uy1(i,j,k)+uy1(i+1,j,k)) - linearVelocity(2) + rotationalComponent(2)
             uzmid = half*(uz1(i,j,k)+uz1(i+1,j,k)) - linearVelocity(3) + rotationalComponent(3)
 
               fcvx= fcvx -uxmid*uzmid*dx*dy
               fcvy= fcvy -uymid*uzmid*dx*dy
               fcvz= fcvz -uzmid*uzmid*dx*dy
 
               !pressure
               prmid = half*(ppi1(i,j,k)+ppi1(i+1,j,k))
               fprz = fprz +prmid*dx*dy
 
               !viscous term
               dudzmid = half*(tg1(i,j,k)+tg1(i+1,j,k))
               dwdxmid = half*(te1(i,j,k)+te1(i+1,j,k))
               dvdzmid = half*(th1(i,j,k)+th1(i+1,j,k))
               dwdymid = half*(tf1(i,j,k)+tf1(i+1,j,k))
               dwdzmid = half*(ti1(i,j,k)+ti1(i+1,j,k))
                                                   
               fdix = fdix -(xnu*(dudzmid+dwdxmid)*dx*dy)
               fdiy = fdiy -(xnu*(dvdzmid+dwdymid)*dx*dy)
               fdiz = fdiz -two*xnu*dwdzmid*dx*dy
 
            enddo
         enddo
 !        drag4(kk)=drag4(kk)+fcvx    ! Should be size ny
         tconvxl2(kk)=tconvxl2(kk)+fcvx
         tconvyl2(kk)=tconvyl2(kk)+fcvy
         tconvzl2(kk)=tconvzl2(kk)+fcvz
         tpreszl(kk) =tpreszl(kk) +fprz
         tdiffxl2(kk)=tdiffxl2(kk)+fdix
         tdiffyl2(kk)=tdiffyl2(kk)+fdiy
         tdiffzl2(kk)=tdiffzl2(kk)+fdiz
      endif     
     
       call MPI_ALLREDUCE(tconvxl,tconvx,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       call MPI_ALLREDUCE(tconvyl,tconvy,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       call MPI_ALLREDUCE(tconvzl,tconvz,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)

       call MPI_ALLREDUCE(tpresxl,tpresx,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       call MPI_ALLREDUCE(tpresyl,tpresy,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       call MPI_ALLREDUCE(tdiffxl,tdiffx,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       call MPI_ALLREDUCE(tdiffyl,tdiffy,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       call MPI_ALLREDUCE(tdiffzl,tdiffz,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)

       call MPI_ALLREDUCE(tconvxl2,tconvx2,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       call MPI_ALLREDUCE(tconvyl2,tconvy2,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       call MPI_ALLREDUCE(tconvzl2,tconvz2,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)     
       call MPI_ALLREDUCE(tpreszl, tpresz ,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       call MPI_ALLREDUCE(tdiffxl2,tdiffx2,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       call MPI_ALLREDUCE(tdiffyl2,tdiffy2,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       call MPI_ALLREDUCE(tdiffzl2,tdiffz2,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
  

       tp1 = sum(tpresx(:))/dt
       tp2 = sum(tpresy(:))/dt
       tp3 = sum(tpresz(:))/dt
    
       mom1 = sum(tunstx(:)) + sum(tconvx(:)) + (-2.0*tconv2_sign+1.0)*sum(tconvx2(:)) !if tconv2sign == 1.0, multiply by -1
       mom2 = sum(tunsty(:)) + sum(tconvy(:)) + (-2.0*tconv2_sign+1.0)*sum(tconvy2(:))
       mom3 = sum(tunstz(:)) + sum(tconvz(:)) + (-2.0*tconv2_sign+1.0)*sum(tconvz2(:))
  
       dra1 = (sum(tdiffx) + sum(tdiffx2) + tp1 - mom1)
       dra2 = (sum(tdiffy) + sum(tdiffy2) + tp2 - mom2)
       dra3 = (sum(tdiffz) + sum(tdiffz2) + tp3 - mom3)
       
       do k=1,zsize(3)

          tpresx(k)=tpresx(k)/dt
          tpresy(k)=tpresy(k)/dt
          tpresz(k)=tpresz(k)/dt


          xmom    = tunstx(k)+tconvx(k)+tconvx2(k)
          ymom    = tunsty(k)+tconvy(k)+tconvy2(k)
          zmom    = tunstz(k)+tconvz(k)+tconvz2(k)
          xDrag(k) = (tdiffx(k)+tdiffx2(k)+tpresx(k)-xmom)
          yLift(k) = (tdiffy(k)+tdiffy2(k)+tpresy(k)-ymom)
          zLat(k)  = (tdiffz(k)+tdiffz2(k)+tpresz(k)-zmom)

       enddo

       !Edited by F. Schuch
       xDrag_mean = sum(xDrag(:))/real(nz,mytype)
       yLift_mean = sum(yLift(:))/real(nz,mytype)

      !  xDrag_tot = sum(xDrag(:))
      !  yLift_tot = sum(yLift(:))
      !  zLat_tot  = sum(zLat(:))
      
      if ((itime==ifirst).or.(itime==0)) then
         
      endif
       if ((nrank .eq. 0).and.(record_var.eq.1)) then
         ! write(*,*) 'TIME STEP = ', itime
         write(38+(iv-1),*) t,dra1,dra2,dra3, sum(tdiffx), sum(tdiffx2), tp1, -mom1, -sum(tunstx(:)), -sum(tconvx(:)), -sum(tconvx2(:))
         !  write(*,*) 'written to file number', 38+(iv-1), t, dra1,dra2,dra3
          call flush(38+(iv-1))

         !  write(12 ,*) t, position(1), position(2), position(3), orientation(1), orientation(2), orientation(3), orientation(4), linearVelocity(1), linearVelocity(2), linearVelocity(3), angularVelocity(1), angularVelocity(2), angularVelocity(3), linearAcceleration(1), linearAcceleration(2), linearAcceleration(3)
       endif
      !  if (mod(itime, ioutput).eq.0) then
      !     if (nrank .eq. 0) then
      !        write(filename,"('forces.dat',I7.7)") itime
      !        call system("cp forces.dat " //filename)
      !     endif
      !  endif
    enddo

    if (torques_flag.eq.0) then !only update velocity fields if torque calculation not called.
      do k = 1, xsize(3)
         do j = 1, xsize(2)
            do i = 1, xsize(1)
               ux11(i,j,k)=ux01(i,j,k)
               uy11(i,j,k)=uy01(i,j,k)
               uz11(i,j,k)=uz01(i,j,k)
               ux01(i,j,k)=ux1(i,j,k)
               uy01(i,j,k)=uy1(i,j,k)
               uz01(i,j,k)=uz1(i,j,k)
            enddo
         enddo
      enddo
   endif

    return

  end subroutine force


  subroutine torque_calc(ux1,uy1,uz1,ep1,dra1,dra2,dra3,record_var)

   USE param
   USE variables
   USE decomp_2d
   USE MPI
   USE ibm_param
   USE ellipsoid_utils, only : CrossProduct,centrifugal_force,coriolis_force

   use var, only : ta1, tb1, tc1, td1, te1, tf1, tg1, th1, ti1, di1
   use var, only : ux2, uy2, uz2, ta2, tb2, tc2, td2, te2, tf2, tg2, th2, ti2, di2
   use var, only : ux3, uy3, uz3, ta3, tb3, tc3, td3, te3, tf3, tg3, th3, ti3, di3
     
 
   implicit none
   character(len=30) :: filename, filename2
   integer :: nzmsize
   integer                                             :: i, iv, j, k, kk, code, jj, ii
   integer                                             :: nvect1,nvect2,nvect3

   real(mytype), dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1, uy1, uz1
   real(mytype), dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ep1
   integer, intent(in) ::record_var
   real(mytype), intent(out)                                       :: dra1,dra2,dra3

   real(mytype), dimension(ysize(1),ysize(2),ysize(3)) :: ppi2
   real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: ppi3

   real(mytype), dimension(nz) :: yLift,xDrag, zLat
   real(mytype) :: yLift_mean,xDrag_mean,zLat_mean
   real(mytype) :: xm,ym,zm,rotationalComponent(3)

   real(mytype), dimension(ny-1) :: del_y

   real(mytype), dimension(nz) :: tunstxl, tunstyl, tunstzl
   real(mytype), dimension(nz) :: tconvxl,tconvyl,tconvzl
   real(mytype), dimension(nz) :: tpresxl,tpresyl
   real(mytype), dimension(nz) :: tdiffxl,tdiffyl,tdiffzl

   real(mytype), dimension(nz) :: tunstx, tunsty, tunstz
   real(mytype), dimension(nz) :: tconvx,tconvy,tconvz
   real(mytype), dimension(nz) :: tpresx,tpresy
   real(mytype), dimension(nz) :: tdiffx,tdiffy,tdiffz

       
   
   real(mytype), dimension(ny) :: tconvxl2, tconvyl2, tconvzl2
   real(mytype), dimension(ny) :: tdiffxl2, tdiffyl2, tdiffzl2
   real(mytype), dimension(ny) :: tconvx2, tconvy2, tconvz2
   real(mytype), dimension(ny) :: tdiffx2, tdiffy2, tdiffz2
   real(mytype), dimension(ny) :: tpreszl, tpresz
   
   
   real(mytype) :: uxmid,uymid,uzmid,prmid
   real(mytype) :: dudxmid,dudymid,dudzmid,dvdxmid,dvdymid,dvdzmid
   real(mytype) :: dwdxmid,dwdymid,dwdzmid
   real(mytype) :: fac,fac1,fac2,fac3,tsumx,tsumy,tsumz,centrifugal(3),coriolis(3)
   real(mytype) :: fcvx,fcvy,fcvz,fprx,fpry,fprz,fdix,fdiy,fdiz
   real(mytype) :: xmom,ymom,zmom
   real(mytype), dimension(ny) :: ztpresx, ztpresy
   real(mytype), dimension(nz) :: zyLift, zxDrag, zzLat
   real(mytype) :: zyLift_mean, zxDrag_mean, zzLat_mean


   real(mytype), dimension(nz) :: drag1, drag2, drag11, drag22
   real(mytype), dimension(nz) :: drag3, drag4, drag33, drag44
   real(mytype) :: mom1, mom2, mom3, tp1, tp2, tp3
   real(mytype) :: radial(3),angular_velocity_result(3)

  !  write(*,*) 'Inside FORCE'

 

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
     do iv=1,nvol
        if ((nrank .eq. 0).and.(record_var.eq.1)) then
           write(filename,"('torques.dat',I1.1)") iv
           open(45+(iv-1),file=filename,status='unknown',form='formatted')
           ! write(*,*) 'Opened file: ', filename, 'number = ', 38+(iv-1)
        endif
        if ((nrank.eq.0).and.(torq_debug.eq.1)) then
            open(100,file="backside.dat",status='unknown',form='formatted')
            open(101,file="frontside.dat",status='unknown',form='formatted')
            open(102,file="lowerside.dat",status='unknown',form='formatted')
            open(103,file="upperside.dat",status='unknown',form='formatted')
            open(104,file="leftside.dat",status='unknown',form='formatted')
            open(105,file="rightside.dat",status='unknown',form='formatted')
            open(106,file="ytorque_decomposition.dat",status='unknown',form='formatted')
            open(107,file="ztorque_decomposition.dat",status='unknown',form='formatted')
        endif

     enddo
      do k = 1, xsize(3)
         do j = 1, xsize(2)
            do i = 1, xsize(1)
               ux11(i,j,k)=ux1(i,j,k)
               uy11(i,j,k)=uy1(i,j,k)
               uz11(i,j,k)=uz1(i,j,k)
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
               uz01(i,j,k)=uz1(i,j,k)
            enddo
         enddo
      enddo
      return
   endif

   call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,1)    ! dudx !x is 1
   call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,2) ! dvdx !y is 2
   call derx (te1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,3) ! dw/dx!z is 3

   call transpose_x_to_y(ta1,ta2) ! dudx
   call transpose_x_to_y(tb1,tb2) ! dvdx
   call transpose_x_to_y(te1,te2) ! dw/dx

   call transpose_x_to_y(ux1,ux2)
   call transpose_x_to_y(uy1,uy2)
   call transpose_x_to_y(uz1,uz2)
   call transpose_x_to_y(ppi1,ppi2)

   call dery (tc2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,1) ! dudy !x is 1
   call dery (td2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,2)    ! dvdy !y is 2
   call dery (tf2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,3) ! dw/dy!z is 3
   

   call transpose_y_to_z(ux2,ux3)
   call transpose_y_to_z(uy2,uy3)
   call transpose_y_to_z(uz2,uz3)

   call derz (tg3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,1)  ! du/dz
   call derz (th3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,2)  ! dv/dz
   call derz (ti3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,3)     ! dw/dz
 
   call transpose_z_to_y(tg3,tg2) ! du/dz
   call transpose_z_to_y(th3,th2) ! dv/dz
   call transpose_z_to_y(ti3,ti2) ! 

   call transpose_y_to_x(tc2,tc1) ! dudy
   call transpose_y_to_x(td2,td1) ! dvdy
   call transpose_y_to_x(th2,th1) ! dv/dz
   call transpose_y_to_x(tf2,tf1) ! dw/dy
   call transpose_y_to_x(tg2,tg1) ! 
   call transpose_y_to_x(ti2,ti1) ! 
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
      tunstzl=zero
      do k=zcvlf_lx(iv),zcvrt_lx(iv)
         tsumx=zero
         tsumy=zero
         tsumz=zero
         zm=real(xstart(3)+k-1,mytype)*dz
         do j=jcvlw_lx(iv),jcvup_lx(iv)
           ym=real(xstart(2)+j-1,mytype)*dy
            do i=icvlf_lx(iv),icvrt_lx(iv)
              xm=real(xstart(1)+i-1,mytype)*dx

              fac1   = (onepfive*ux1(i,j,k)-two*ux01(i,j,k)+half*ux11(i,j,k))*(one-ep1(i,j,k))
              fac2   = (onepfive*uy1(i,j,k)-two*uy01(i,j,k)+half*uy11(i,j,k))*(one-ep1(i,j,k))
              fac3   = (onepfive*uz1(i,j,k)-two*uz01(i,j,k)+half*uz11(i,j,k))*(one-ep1(i,j,k))

               ! call coriolis_force(angularVelocity,[fac1,fac2,fac3],coriolis)
               ! call centrifugal_force(angularVelocity, [xm,ym,zm]-position,centrifugal)
               !     The velocity time rate has to be relative to the cell center, 
               !     and not to the nodes, because, here, we have an integral 
               !     relative to the volume, and, therefore, this has a sense 
               !     of a "source".
               !         fac   = (1.5*ux1(i,j,k)-2.0*ux01(i,j,k)+0.5*ux11(i,j,k))*epcv1(i,j,k)
               call crossProduct(([fac1,fac2,fac3]),[xm,ym,zm]-position, angular_velocity_result)
               tsumx = tsumx+angular_velocity_result(1)*dx*del_y(j+(xstart(2)-1))*dz/dt    !tsumx+fac*dx*dy/dt
               ! tsumx = tsumx+fac1*dx*del_y(j+xstart(2)-1)*dz/dt
               !sumx(k) = sumx(k)+dudt1*dx*dy

               !         fac   = (1.5*uy1(i,j,k)-2.0*uy01(i,j,k)+0.5*uy11(i,j,k))*epcv1(i,j,k)
               tsumy = tsumy+angular_velocity_result(2)*dx*del_y(j+(xstart(2)-1))*dz/dt !tsumy+fac*dx*dy/dt
               ! tsumy = tsumy+fac2*dx*del_y(j+xstart(2)-1)*dz/dt
               !sumy(k) = sumy(k)+dudt1*dx*dy

               tsumz = tsumz+angular_velocity_result(3)*dx*del_y(j+(xstart(2)-1))*dz/dt     
               ! tsumz = tsumz+fac3*dx*del_y(j+xstart(2)-1)*dz/dt
            enddo
         enddo
         tunstxl(xstart(3)-1+k)=tsumx
         tunstyl(xstart(3)-1+k)=tsumy
         tunstzl(xstart(3)-1+k)=tsumz
      enddo
      call MPI_ALLREDUCE(tunstxl,tunstx,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
      call MPI_ALLREDUCE(tunstyl,tunsty,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
      call MPI_ALLREDUCE(tunstzl,tunstz,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)

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
      
 drag1(:)=0.
 drag2(:)=0.
 drag3(:)=0.
 drag4(:)=0.
 
 drag11(:)=0.
 drag22(:)=0.
 drag33(:)=0.
 drag44(:)=0.

      tconvxl=zero
      tconvyl=zero
      tconvzl=zero
      tdiffxl=zero
      tdiffyl=zero
      tdiffzl=zero
      tpresxl=zero
      tpresyl=zero
      tpreszl=zero

      tconvxl2=zero
      tconvyl2=zero
      tconvzl2=zero
      tdiffxl2=zero
      tdiffyl2=zero
      tdiffzl2=zero  
      !BC and AD : x-pencils
      !AD
      if ((jcvlw(iv).ge.xstart(2)).and.(jcvlw(iv).le.xend(2))) then
         j=jcvlw(iv)-xstart(2)+1
         jj=jcvlw(iv)
         ym=real(jj,mytype)*dy
         do k=zcvlf_lx(iv),zcvrt_lx(iv)
            kk=xstart(3)-1+k
            zm=real(kk,mytype)*dz
            fcvx=zero
            fcvy=zero
            fcvz=zero
            fprx=zero
            fpry=zero
            fprz=zero
            fdix=zero
            fdiy=zero
            fdiz=zero
            do i=icvlf_lx(iv),icvrt_lx(iv)-1
               
               ii=xstart(1)+i-1
               xm=real(ii,mytype)*dx
              !  write(*,*) 'Calculating force at upper y boundary', [xm,ym,zm]
               radial = [xm,ym,zm]-position

               !momentum flux
               call crossProduct(angularVelocity(2:4),[xm,ym,zm]-position,rotationalComponent)
               uxmid = half*(ux1(i,j,k)+ux1(i+1,j,k)) - linearVelocity(1) + rotationalComponent(1)
               uymid = half*(uy1(i,j,k)+uy1(i+1,j,k)) - linearVelocity(2) + rotationalComponent(2)
               uzmid = half*(uz1(i,j,k)+uz1(i+1,j,k)) - linearVelocity(3) + rotationalComponent(3)

               fcvx = fcvx -(uymid*radial(3)-uzmid*radial(2))*uymid*dx*dz
               fcvy = fcvy -(uzmid*radial(1)-uxmid*radial(3))*uymid*dx*dz
               fcvz = fcvz -(uxmid*radial(2)-uymid*radial(1))*uymid*dx*dz


               !pressure
               prmid = half*(ppi1(i,j,k)+ppi1(i+1,j,k))
               ! fpry  = fpry +prmid*dx*dz*(radial(1)-radial(3))
               fprz = fprz +prmid*dx*dz*(radial(1))
               fprx = fprx +prmid*dx*dz*(-radial(3))

               !viscous term
               dudymid = half*(tc1(i,j,k)+tc1(i+1,j,k))
               dvdxmid = half*(tb1(i,j,k)+tb1(i+1,j,k))
               dvdymid = half*(td1(i,j,k)+td1(i+1,j,k))
               dwdymid = half*(tf1(i,j,k)+tf1(i+1,j,k))
               dvdzmid = half*(th1(i,j,k)+th1(i+1,j,k))

               fdix = fdix - (xnu*(two*dvdymid)*radial(3)-xnu*(dwdymid+dvdzmid)*radial(2))*dx*dz
               fdiy = fdiy - (xnu*(dwdymid+dvdzmid)*radial(1)-xnu*(dudymid+dvdxmid)*radial(3))*dx*dz
               fdiz = fdiz - (xnu*(dudymid+dvdxmid)*radial(2)-xnu*(two*dvdymid)*radial(1))*dx*dz

               ! fdix    = fdix -(xnu*(dudymid+dvdxmid)*dx*dz)
               ! fdiy    = fdiy -two*xnu*dvdymid*dx*dz
               ! fdiz    = fdiz -(xnu*(dwdymid+dvdzmid)*dx*dz)


            enddo

            tconvxl(1)=tconvxl(1)+fcvx
            tconvyl(1)=tconvyl(1)+fcvy
            tconvzl(1)=tconvzl(1)+fcvz

            tpresxl(1)=tpresxl(1)+fprx
            ! tpresyl(1)=tpresyl(1)+fpry
            tpreszl(1)=tpreszl(1)+fprz

            tdiffxl(1)=tdiffxl(1)+fdix
            tdiffyl(1)=tdiffyl(1)+fdiy
            tdiffzl(1)=tdiffzl(1)+fdiz
         enddo

      endif
      !BC
      if ((jcvup(iv).ge.xstart(2)).and.(jcvup(iv).le.xend(2))) then
         j=jcvup(iv)-xstart(2)+1
         jj=jcvup(iv)
         ym=real(jj,mytype)*dy
         do k=zcvlf_lx(iv),zcvrt_lx(iv)
            kk=xstart(3)-1+k   
            zm=real(kk,mytype)*dz
            fcvx=zero
            fcvy=zero
            fcvz=zero
            fprx=zero
            fpry=zero
            fprz=zero
            fdix=zero
            fdiy=zero
            fdiz=zero
            do i=icvlf_lx(iv),icvrt_lx(iv)-1
              ii=xstart(1)+i-1
              xm=real(ii,mytype)*dx
              ! write(*,*) 'xm = ', xm
              ! write(*,*) 'Calculating force at lower y boundary', [xm,ym,zm]

              !momentum flux
              call crossProduct(angularVelocity(2:4),[xm,ym,zm]-position,rotationalComponent)
              radial = [xm,ym,zm]-position
              uxmid = half*(ux1(i,j,k)+ux1(i+1,j,k)) - linearVelocity(1) + rotationalComponent(1)
              uymid = half*(uy1(i,j,k)+uy1(i+1,j,k)) - linearVelocity(2) + rotationalComponent(2)
              uzmid = half*(uz1(i,j,k)+uz1(i+1,j,k)) - linearVelocity(3) + rotationalComponent(3)

               fcvx = fcvx +(uymid*radial(3)-uzmid*radial(2))*uymid*dx*dz
               fcvy = fcvy +(uzmid*radial(1)-uxmid*radial(3))*uymid*dx*dz
               fcvz = fcvz +(uxmid*radial(2)-uymid*radial(1))*uymid*dx*dz


               !pressure
               prmid = half*(ppi1(i,j,k)+ppi1(i+1,j,k))
               ! fpry = fpry -prmid*dx*dz*(radial(1)-radial(3))
               fprz = fprz -prmid*dx*dz*(radial(1))
               fprx = fprx -prmid*dx*dz*(-radial(3))

               !viscous term
               dudymid = half*(tc1(i,j,k)+tc1(i+1,j,k))
               dvdxmid = half*(tb1(i,j,k)+tb1(i+1,j,k))
               dvdymid = half*(td1(i,j,k)+td1(i+1,j,k))
               dwdymid = half*(tf1(i,j,k)+tf1(i+1,j,k))
               dvdzmid = half*(th1(i,j,k)+th1(i+1,j,k))
 
               fdix = fdix + (xnu*(two*dvdymid)*radial(3)-xnu*(dwdymid+dvdzmid)*radial(2))*dx*dz
               fdiy = fdiy + (xnu*(dwdymid+dvdzmid)*radial(1)-xnu*(dudymid+dvdxmid)*radial(3))*dx*dz
               fdiz = fdiz + (xnu*(dudymid+dvdxmid)*radial(2)-xnu*(two*dvdymid)*radial(1))*dx*dz

               ! fdix = fdix + (xnu*(dudymid+dvdxmid)*dx*dz)
               ! fdiy = fdiy + (xnu*(two*dvdymid)*dx*dz)
               ! fdiz = fdiz + (xnu*(dwdymid+dvdzmid)*dx*dz)

            enddo
            tconvxl(2)=tconvxl(2)+fcvx
            tconvyl(2)=tconvyl(2)+fcvy
            tconvzl(2)=tconvzl(2)+fcvz

            tpresxl(2)=tpresxl(2)+fprx
            ! tpresyl(2)=tpresyl(2)+fpry
            tpreszl(2)=tpreszl(2)+fprz

            tdiffxl(2)=tdiffxl(2)+fdix
            tdiffyl(2)=tdiffyl(2)+fdiy
            tdiffzl(2)=tdiffzl(2)+fdiz

         enddo
      endif
      !AB and DC : y-pencils
      !AB
      if ((icvlf(iv).ge.ystart(1)).and.(icvlf(iv).le.yend(1))) then
         i=icvlf(iv)-ystart(1)+1
         ii=icvlf(iv)-1
         xm=real(ii,mytype)*dx
         do k=zcvlf_ly(iv),zcvrt_ly(iv)
            kk=ystart(3)+k-1
            zm=real(kk,mytype)*dz
            fcvx=zero
            fcvy=zero
            fcvz=zero
            fprx=zero
            fpry=zero
            fprz=zero
            fdix=zero
            fdiy=zero
            fdiz=zero
            do j=jcvlw_ly(iv),jcvup_ly(iv)-1

               jj=ystart(2)+j-1
               ym=real(jj,mytype)*dz
              !  write(*,*) 'Calculating force at left x boundary', [xm,ym,zm]
               !momentum flux
               radial = [xm,ym,zm]-position

               call crossProduct(angularVelocity(2:4),[xm,ym,zm]-position,rotationalComponent)
               uxmid = half*(ux2(i,j,k)+ux2(i,j+1,k)) - linearVelocity(1) + rotationalComponent(1)
               uymid = half*(uy2(i,j,k)+uy2(i,j+1,k)) - linearVelocity(2) + rotationalComponent(2)
               uzmid = half*(uz2(i,j,k)+uz2(i,j+1,k)) - linearVelocity(3) + rotationalComponent(3)

               fcvx = fcvx -(uymid*radial(3)-uzmid*radial(2))*uxmid*del_y(j)*dz
               fcvy = fcvy -(uzmid*radial(1)-uxmid*radial(3))*uxmid*del_y(j)*dz
               fcvz = fcvz -(uxmid*radial(2)-uymid*radial(1))*uxmid*del_y(j)*dz


               !pressure
               prmid = half*(ppi2(i,j,k)+ppi2(i,j+1,k))
               ! fprx = fprx +prmid*del_y(j)*dz*(radial(3)-radial(2))

               fpry = fpry + prmid*del_y(j)*dz*(radial(3))
               fprz = fprz + prmid*del_y(j)*dz*(-radial(2))

               !viscous term
               dudxmid = half*(ta2(i,j,k)+ta2(i,j+1,k))
               dudymid = half*(tc2(i,j,k)+tc2(i,j+1,k))
               dvdxmid = half*(tb2(i,j,k)+tb2(i,j+1,k))
               dwdxmid = half*(te2(i,j,k)+te2(i,j+1,k))
               dudzmid = half*(tg2(i,j,k)+tg2(i,j+1,k))

               fdix = fdix - (xnu*(dvdxmid+dudymid)*radial(3)-xnu*(dwdxmid+dudzmid)*radial(2))*del_y(j)*dz
               fdiy = fdiy - (xnu*(dwdxmid+dudzmid)*radial(1)-xnu*(two*dudxmid)*radial(3))*del_y(j)*dz
               fdiz = fdiz - (xnu*(two*dudxmid)*radial(2)-xnu*(dvdxmid+dudymid)*radial(1))*del_y(j)*dz


               ! fdix = fdix -two*xnu*dudxmid*del_y(j)*dz
               ! fdiy = fdiy -xnu*(dvdxmid+dudymid)*del_y(j)*dz
               ! fdiz = fdiz -xnu*(dwdxmid+dudzmid)*del_y(j)*dz
            enddo
            tconvxl(3)=tconvxl(3)+fcvx
            tconvyl(3)=tconvyl(3)+fcvy
            tconvzl(3)=tconvzl(3)+fcvz

            ! tpresxl(3)=tpresxl(3)+fprx
            tpresyl(3)=tpresyl(3)+fpry
            tpreszl(3)=tpreszl(3)+fprz

            tdiffxl(3)=tdiffxl(3)+fdix
            tdiffyl(3)=tdiffyl(3)+fdiy
            tdiffzl(3)=tdiffzl(3)+fdiz

         enddo
      endif
      !DC
      if ((icvrt(iv).ge.ystart(1)).and.(icvrt(iv).le.yend(1))) then
         i=icvrt(iv)-ystart(1)+1
         ii=icvrt(iv)
         xm=real(ii,mytype)*dx
         do k=zcvlf_ly(iv),zcvrt_ly(iv)
            kk=ystart(3)-1+k
            zm=real(kk,mytype)*dz
            fcvx=zero
            fcvy=zero
            fcvz=zero
            fprx=zero
            fpry=zero
            fprz=zero
            fdix=zero
            fdiy=zero
            fdiz=zero
            do j=jcvlw_ly(iv),jcvup_ly(iv)-1 
               jj=ystart(2)+j-1
               ym=real(jj,mytype)*dy
              !  write(*,*) 'Calculating force at right x boundary', [xm,ym,zm]

               radial = [xm,ym,zm]-position

               !momentum flux
               call crossProduct(angularVelocity(2:4),[xm,ym,zm]-position,rotationalComponent)
               uxmid = half*(ux2(i,j,k)+ux2(i,j+1,k)) - linearVelocity(1) + rotationalComponent(1)
               uymid = half*(uy2(i,j,k)+uy2(i,j+1,k)) - linearVelocity(2) + rotationalComponent(2)
               uzmid = half*(uz2(i,j,k)+uz2(i,j+1,k)) - linearVelocity(3) + rotationalComponent(3)


               fcvx = fcvx +(uymid*radial(3)-uzmid*radial(2))*uxmid*del_y(j)*dz
               fcvy = fcvy +(uzmid*radial(1)-uxmid*radial(3))*uxmid*del_y(j)*dz
               fcvz = fcvz +(uxmid*radial(2)-uymid*radial(1))*uxmid*del_y(j)*dz



               !pressure
               prmid = half*(ppi2(i,j,k)+ppi2(i,j+1,k))
               ! fprx = fprx -prmid*del_y(j)*dz*(radial(3)-radial(2))

               fpry = fpry - prmid*del_y(j)*dz*(radial(3))
               fprz = fprz - prmid*del_y(j)*dz*(-radial(2))  

               !viscous term
               dudxmid = half*(ta2(i,j,k)+ta2(i,j+1,k))
               dudymid = half*(tc2(i,j,k)+tc2(i,j+1,k))
               dvdxmid = half*(tb2(i,j,k)+tb2(i,j+1,k))
               dwdxmid = half*(te2(i,j,k)+te2(i,j+1,k))
               dudzmid = half*(tg2(i,j,k)+tg2(i,j+1,k))

               fdix = fdix + (xnu*(dvdxmid+dudymid)*radial(3)-xnu*(dwdxmid+dudzmid)*radial(2))*del_y(j)*dz
               fdiy = fdiy + (xnu*(dwdxmid+dudzmid)*radial(1)-xnu*(two*dudxmid)*radial(3))*del_y(j)*dz
               fdiz = fdiz + (xnu*(two*dudxmid)*radial(2)-xnu*(dvdxmid+dudymid)*radial(1))*del_y(j)*dz


               ! fdix = fdix + two*xnu*dudxmid*del_y(j)*dz
               ! fdiy = fdiy + xnu*(dvdxmid+dudymid)*del_y(j)*dz
               ! fdiz = fdiz + xnu*(dwdxmid+dudzmid)*del_y(j)*dz

            enddo
            tconvxl(4)=tconvxl(4)+fcvx
            tconvyl(4)=tconvyl(4)+fcvy
            tconvzl(4)=tconvzl(4)+fcvz

            ! tpresxl(4)=tpresxl(4)+fprx
            tpresyl(4)=tpresyl(4)+fpry
            tpreszl(4)=tpreszl(4)+fprz

            tdiffxl(4)=tdiffxl(4)+fdix
            tdiffyl(4)=tdiffyl(4)+fdiy
            tdiffzl(4)=tdiffzl(4)+fdiz

         enddo
      endif

      !Left & Right : 
      !Left
      if ((zcvlf(iv).ge.xstart(3)).and.(zcvlf(iv).le.xend(3))) then
        k=zcvlf(iv)-xstart(3)+1
        kk=zcvlf(iv)
        zm=real(kk,mytype)*dz

        fcvx=zero
        fcvy=zero
        fcvz=zero
        fprx=zero
        fpry=zero
        fprz=zero
        fdiy=zero
        fdiz=zero
        do j=jcvlw_lx(iv),jcvup_lx(iv)
         kk = xstart(2)-1+j
         jj = xstart(2)-1+j

         ym=real(jj,mytype)*dy
           do i=icvlf_lx(iv),icvrt_lx(iv)-1
              ii=xstart(1)+i-1
              xm=real(ii,mytype)*dx
              ! write(*,*) 'Calculating force at left z boundary', [xm,ym,zm]

              !momentum flux
              radial = [xm,ym,zm]-position

              call crossProduct(angularVelocity(2:4),[xm,ym,zm]-position,rotationalComponent)
              uxmid = half*(ux1(i,j,k)+ux1(i+1,j,k)) - linearVelocity(1) + rotationalComponent(1)
              uymid = half*(uy1(i,j,k)+uy1(i+1,j,k)) - linearVelocity(2) + rotationalComponent(2)
              uzmid = half*(uz1(i,j,k)+uz1(i+1,j,k)) - linearVelocity(3) + rotationalComponent(3)

              fcvx = fcvx -(uymid*radial(3)-uzmid*radial(2))*uzmid*dx*del_y(j) !!!CHANGE
              fcvy = fcvy -(uzmid*radial(1)-uxmid*radial(3))*uzmid*dx*del_y(j)
              fcvz = fcvz -(uxmid*radial(2)-uymid*radial(1))*uzmid*dx*del_y(j)

              !pressure
              prmid = half*(ppi1(i,j,k)+ppi1(i+1,j,k))
            !   fprz = fprz -prmid*dx*dy*(radial(2)-radial(1))
              fprx = fprx +prmid*dx*dy*(radial(2))
              fpry = fpry +prmid*dx*dy*(-radial(1))

              !viscous term
              dudzmid = half*(tg1(i,j,k)+tg1(i+1,j,k))
              dwdxmid = half*(te1(i,j,k)+te1(i+1,j,k))
              dvdzmid = half*(th1(i,j,k)+th1(i+1,j,k))
              dwdymid = half*(tf1(i,j,k)+tf1(i+1,j,k))
              dwdzmid = half*(ti1(i,j,k)+ti1(i+1,j,k))
              
              fdix = fdix - (xnu*(dvdzmid+dwdymid)*radial(3)-xnu*(two*dwdzmid)*radial(2))*dx*dy
              fdiy = fdiy - (xnu*(two*dwdzmid)*radial(1)-xnu*(dudzmid+dwdxmid)*radial(3))*dx*dy
              fdiz = fdiz - (xnu*(dudzmid+dwdxmid)*radial(2)-xnu*(dvdzmid+dwdymid)*radial(1))*dx*dy


            !   fdix = fdix +(xnu*(dudzmid+dwdxmid)*dx*dy)
            !   fdiy = fdiy +(xnu*(dvdzmid+dwdymid)*dx*dy)
            !   fdiz = fdiz +(xnu*(two*dwdzmid)*dx*dy)
           enddo
        enddo
!print*, kk
!        drag3(kk)=drag3(kk)+fcvx   ! Should be size ny
!        print*, drag3(kk)
        tconvxl2(5)=tconvxl2(5)+fcvx
        tconvyl2(5)=tconvyl2(5)+fcvy
        tconvzl2(5)=tconvzl2(5)+fcvz
        tpresxl(5) =tpresxl(5) +fprx
        tpresyl(5) =tpresyl(5) +fpry
      !   tpreszl(5) =tpreszl(5) +fprz
        tdiffxl2(5)=tdiffxl2(5)+fdix
        tdiffyl2(5)=tdiffyl2(5)+fdiy
        tdiffzl2(5)=tdiffzl2(5)+fdiz        
     endif 
     !Right
     if ((zcvrt(iv).ge.xstart(3)).and.(zcvrt(iv).le.xend(3))) then
        k=zcvrt(iv)-xstart(3)+1
        kk=zcvrt(iv)
        zm=real(kk,mytype)*dz
!        kk=nrank+1

        fcvx=zero
        fcvy=zero
        fcvz=zero
        fprx=zero
        fpry=zero
        fprz=zero
        fdix=zero
        fdiy=zero
        fdiz=zero
      !  do k=1,xsize(3)
        do j=jcvlw_lx(iv),jcvup_lx(iv)
        !  kk = xstart(2)-1+j 
         jj = xstart(2)-1+j
         ym=real(jj,mytype)*dy
          do i=icvlf_lx(iv),icvrt_lx(iv)-1
            ii=xstart(1)+i-1
            xm=real(ii,mytype)*dx
              !momentum flux
            radial = [xm,ym,zm]-position

            call crossProduct(angularVelocity(2:4),[xm,ym,zm]-position,rotationalComponent)
           !  write(*,*) 'Calculating force at right z boundary', [xm,ym,zm]

            uxmid = half*(ux1(i,j,k)+ux1(i+1,j,k)) - linearVelocity(1) + rotationalComponent(1)
            uymid = half*(uy1(i,j,k)+uy1(i+1,j,k)) - linearVelocity(2) + rotationalComponent(2)
            uzmid = half*(uz1(i,j,k)+uz1(i+1,j,k)) - linearVelocity(3) + rotationalComponent(3)

            fcvx = fcvx +(uymid*radial(3)-uzmid*radial(2))*uzmid*dx*del_y(j)
            fcvy = fcvy +(uzmid*radial(1)-uxmid*radial(3))*uzmid*dx*del_y(j)
            fcvz = fcvz +(uxmid*radial(2)-uymid*radial(1))*uzmid*dx*del_y(j)


              !pressure
              prmid = half*(ppi1(i,j,k)+ppi1(i+1,j,k))
            !   fprz = fprz +prmid*dx*dy*(radial(2)-radial(1))

              fprx = fprx -prmid*dx*dy*(radial(2))
              fpry = fpry -prmid*dx*dy*(-radial(1))

              !viscous term
              dudzmid = half*(tg1(i,j,k)+tg1(i+1,j,k))
              dwdxmid = half*(te1(i,j,k)+te1(i+1,j,k))
              dvdzmid = half*(th1(i,j,k)+th1(i+1,j,k))
              dwdymid = half*(tf1(i,j,k)+tf1(i+1,j,k))
              dwdzmid = half*(ti1(i,j,k)+ti1(i+1,j,k))
                            
              fdix = fdix + (xnu*(dvdzmid+dwdymid)*radial(3)-xnu*(two*dwdzmid)*radial(2))*dx*dy
              fdiy = fdiy + (xnu*(two*dwdzmid)*radial(1)-xnu*(dudzmid+dwdxmid)*radial(3))*dx*dy
              fdiz = fdiz + (xnu*(dudzmid+dwdxmid)*radial(2)-xnu*(dvdzmid+dwdymid)*radial(1))*dx*dy

            !   fdix = fdix -(xnu*(dudzmid+dwdxmid)*dx*dy)
            !   fdiy = fdiy -(xnu*(dvdzmid+dwdymid)*dx*dy)
            !   fdiz = fdiz -two*xnu*dwdzmid*dx*dy

           enddo
        enddo
!        drag4(kk)=drag4(kk)+fcvx    ! Should be size ny
        tconvxl2(6)=tconvxl2(6)+fcvx
        tconvyl2(6)=tconvyl2(6)+fcvy
        tconvzl2(6)=tconvzl2(6)+fcvz
        
        tpresxl(6) =tpresxl(6) +fprx !!!!!!!!
        tpresyl(6) =tpresyl(6) +fpry
      !   tpreszl(6) =tpreszl(6) +fprz

        tdiffxl2(6)=tdiffxl2(6)+fdix
        tdiffyl2(6)=tdiffyl2(6)+fdiy
        tdiffzl2(6)=tdiffzl2(6)+fdiz
     endif     
    
      call MPI_ALLREDUCE(tconvxl,tconvx,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
      call MPI_ALLREDUCE(tconvyl,tconvy,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
      call MPI_ALLREDUCE(tconvzl,tconvz,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)

      call MPI_ALLREDUCE(tpresxl,tpresx,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
      call MPI_ALLREDUCE(tpresyl,tpresy,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
      call MPI_ALLREDUCE(tdiffxl,tdiffx,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
      call MPI_ALLREDUCE(tdiffyl,tdiffy,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
      call MPI_ALLREDUCE(tdiffzl,tdiffz,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)

      call MPI_ALLREDUCE(tconvxl2,tconvx2,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
      call MPI_ALLREDUCE(tconvyl2,tconvy2,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
      call MPI_ALLREDUCE(tconvzl2,tconvz2,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)     
      call MPI_ALLREDUCE(tpreszl, tpresz ,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
      call MPI_ALLREDUCE(tdiffxl2,tdiffx2,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
      call MPI_ALLREDUCE(tdiffyl2,tdiffy2,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
      call MPI_ALLREDUCE(tdiffzl2,tdiffz2,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
 

      tp1 = sum(tpresx(:))/dt
      tp2 = sum(tpresy(:))/dt
      tp3 = sum(tpresz(:))/dt
   
      mom1 = sum(tunstx(:)) + sum(tconvx(:)) + (-2.0*tconv2_sign+1.0)*sum(tconvx2(:)) !if tconv2sign == 1.0, multiply by -1
      mom2 = sum(tunsty(:)) + sum(tconvy(:)) + (-2.0*tconv2_sign+1.0)*sum(tconvy2(:))
      mom3 = sum(tunstz(:)) + sum(tconvz(:)) + (-2.0*tconv2_sign+1.0)*sum(tconvz2(:))
 
      dra1 = -(sum(tdiffx) + sum(tdiffx2) + tp1 - mom1)*(1.0-2.0*torq_flip)
      dra2 = -(sum(tdiffy) + sum(tdiffy2) + tp2 - mom2)*(1.0-2.0*torq_flip)
      dra3 = -(sum(tdiffz) + sum(tdiffz2) + tp3 - mom3)*(1.0-2.0*torq_flip)
      
      ! do k=1,zsize(3)

      !    tpresx(k)=tpresx(k)/dt
      !    tpresy(k)=tpresy(k)/dt
      !    tpresz(k)=tpresz(k)/dt


      !    xmom    = tunstx(k)+tconvx(k)+tconvx2(k)
      !    ymom    = tunsty(k)+tconvy(k)+tconvy2(k)
      !    zmom    = tunstz(k)+tconvz(k)+tconvz2(k)
      !    xDrag(k) = (tdiffx(k)+tdiffx2(k)+tpresx(k)-xmom)
      !    yLift(k) = (tdiffy(k)+tdiffy2(k)+tpresy(k)-ymom)
      !    zLat(k)  = (tdiffz(k)+tdiffz2(k)+tpresz(k)-zmom)

      ! enddo

      !Edited by F. Schuch
      xDrag_mean = sum(xDrag(:))/real(nz,mytype)
      yLift_mean = sum(yLift(:))/real(nz,mytype)

     !  xDrag_tot = sum(xDrag(:))
     !  yLift_tot = sum(yLift(:))
     !  zLat_tot  = sum(zLat(:))
     
     if ((itime==ifirst).or.(itime==0)) then
        
     endif
      if ((nrank .eq. 0).and.(record_var.eq.1)) then
        ! write(*,*) 'TIME STEP = ', itime
        write(45+(iv-1),*) t,dra1,dra2,dra3, sum(tdiffx), sum(tdiffx2), tp1, -mom1, -sum(tunstx(:)), -sum(tconvx(:)), -sum(tconvx2(:))
        !  write(*,*) 'written to file number', 38+(iv-1), t, dra1,dra2,dra3
         call flush(45+(iv-1))

         if (torq_debug.eq.1) then 
            write(100,*) t,dra1,dra2,dra3, tdiffx(1), tdiffx2(1), tpresx(1)/dt,  -tunstx(1), -tconvx(1), -tconvx2(1)
            write(101,*) t,dra1,dra2,dra3, tdiffx(2), tdiffx2(2), tpresx(2)/dt,  -tunstx(2), -tconvx(2), -tconvx2(2)
            write(102,*) t,dra1,dra2,dra3, tdiffx(3), tdiffx2(3), tpresx(3)/dt,  -tunstx(3), -tconvx(3), -tconvx2(3)
            write(103,*) t,dra1,dra2,dra3, tdiffx(4), tdiffx2(4), tpresx(4)/dt,  -tunstx(4), -tconvx(4), -tconvx2(4)
            write(104,*) t,dra1,dra2,dra3, tdiffx(5), tdiffx2(5), tpresx(5)/dt,  -tunstx(5), -tconvx(5), -tconvx2(5)
            write(105,*) t,dra1,dra2,dra3, tdiffx(6), tdiffx2(6), tpresx(6)/dt,  -tunstx(6), -tconvx(6), -tconvx2(6)
   
            write(106,*) t,dra1,dra2,dra3, sum(tdiffy), sum(tdiffy2), tp2, -mom2, -sum(tunsty(:)), -sum(tconvy(:)), -sum(tconvy2(:))
            write(107,*) t,dra1,dra2,dra3, sum(tdiffz), sum(tdiffz2), tp3, -mom3, -sum(tunstz(:)), -sum(tconvz(:)), -sum(tconvz2(:))
            
            call flush(100)
            call flush(101)
            call flush(102)
            call flush(103)
            call flush(104)
            call flush(105)

            call flush(106)
            call flush(107)
            endif
      endif
     !  if (mod(itime, ioutput).eq.0) then
     !     if (nrank .eq. 0) then
     !        write(filename,"('forces.dat',I7.7)") itime
     !        call system("cp forces.dat " //filename)
     !     endif
     !  endif
   enddo

   do k = 1, xsize(3)
      do j = 1, xsize(2)
         do i = 1, xsize(1)
            ux11(i,j,k)=ux01(i,j,k)
            uy11(i,j,k)=uy01(i,j,k)
            uz11(i,j,k)=uz01(i,j,k)
            ux01(i,j,k)=ux1(i,j,k)
            uy01(i,j,k)=uy1(i,j,k)
            uz01(i,j,k)=uz1(i,j,k)
         enddo
      enddo
   enddo

   return

 end subroutine torque_calc
end module forces
