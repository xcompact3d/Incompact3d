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
    USE ibm
    implicit none
  
    integer :: nvol,iforces
    real(mytype),save,allocatable,dimension(:,:,:) :: ux01, uy01, ux11, uy11, ppi1, uz01, uz11 
    real(mytype),allocatable,dimension(:) :: xld, xrd, yld, yud, xld2, xrd2, yld2, yud2, zld, zrd
    integer,allocatable,dimension(:) :: icvlf, icvrt, jcvlw, jcvup, zcvlf, zcvrt
    integer,allocatable,dimension(:) :: icvlf_lx, icvrt_lx, icvlf_ly, icvrt_ly, icvlf_lz, icvrt_lz
    integer,allocatable,dimension(:) :: jcvlw_lx, jcvup_lx, jcvlw_ly, jcvup_ly, jcvlw_lz, jcvup_lz
    integer,allocatable,dimension(:) :: zcvlf_lx, zcvrt_lx, zcvlf_ly, zcvrt_ly
    
  
  contains
  
    subroutine init_forces
  
      USE decomp_2d
      USE param
      USE variables
      implicit none
  
      integer :: iv
  
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
  
      allocate(icvlf(nvol), icvrt(nvol), jcvlw(nvol), jcvup(nvol), zcvlf(nvol), zcvrt(nvol))
      allocate(icvlf_lx(nvol), icvrt_lx(nvol), icvlf_ly(nvol), icvrt_ly(nvol), icvlf_lz(nvol), icvrt_lz(nvol))
      allocate(jcvlw_lx(nvol), jcvup_lx(nvol), jcvlw_ly(nvol), jcvup_ly(nvol), jcvlw_lz(nvol), jcvup_lz(nvol))
      allocate(zcvlf_lx(nvol), zcvrt_lx(nvol), zcvlf_ly(nvol), zcvrt_ly(nvol))
      allocate(xld2(nvol), xrd2(nvol), yld2(nvol), yud2(nvol))
  
      ! Update Control Volume based on moving cylinder case
      if ((iibm.ne.0).and.(t.ne.0.)) then
         xld2(:) = xld(:) + (t-ifirst*dt)*ubcx
         xrd2(:) = xrd(:) + (t-ifirst*dt)*ubcx
         yld2(:) = yld(:) + (t-ifirst*dt)*ubcy
         yud2(:) = yud(:) + (t-ifirst*dt)*ubcy
      else
         xld2(:) = xld(:)
         xrd2(:) = xrd(:)
         yld2(:) = yld(:)
         yud2(:) = yud(:)
      endif
      
      !     Definition of the Control Volume
      !*****************************************************************
      !! xld,xrd,yld,yud: limits of control volume (!!don't use cex and cey anymore!!)
      do iv=1,nvol
         ! ok for istret=0 (!!to do for istret=1!!)
         icvlf(iv) = nint(xld2(iv)/dx)+1
         icvrt(iv) = nint(xrd2(iv)/dx)+1
         jcvlw(iv) = nint(yld2(iv)/dy)+1
         jcvup(iv) = nint(yud2(iv)/dy)+1
  
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
  
      if (nrank==0) then
         print *,'========================Forces============================='
         print *,'                       (icvlf)      (icvrt) '
         print *,'                (jcvup) B____________C '
         print *,'                        \            \ '
         print *,'                        \     __     \ '
         print *,'                        \    \__\    \ '
         print *,'                        \            \ '
         print *,'                        \       CV   \ '
         print *,'                (jcvlw) A____________D '
         do iv=1,nvol
            write(*,"(' Control Volume     : #',I1)") iv
            write(*,"('     xld, icvlf     : (',F6.2,',',I6,')')") xld(iv), icvlf(iv)
            write(*,"('     xrd, icvrt     : (',F6.2,',',I6,')')") xrd(iv), icvrt(iv)
            write(*,"('     yld, jcvlw     : (',F6.2,',',I6,')')") yld(iv), jcvlw(iv)
            write(*,"('     yud, jcvup     : (',F6.2,',',I6,')')") yud(iv), jcvup(iv)
         enddo
         print *,'==========================================================='
      endif
    end subroutine init_forces
  
  !***********************************************************************
  ! 
!     subroutine update_forces
!   ! 
!   !***********************************************************************
  
!       USE decomp_2d
!       USE param
!       USE variables
!       implicit none
  
!       integer :: iv
  
!       ! Update Control Volume based on moving cylinder case
!       if ((iibm.ne.0).and.(t.ne.0.)) then
!          xld2(:) = xld(:) + (t-ifirst*dt)*ubcx
!          xrd2(:) = xrd(:) + (t-ifirst*dt)*ubcx
!          yld2(:) = yld(:) + (t-ifirst*dt)*ubcy
!          yud2(:) = yud(:) + (t-ifirst*dt)*ubcy
!       else
!          xld2(:) = xld(:)
!          xrd2(:) = xrd(:)
!          yld2(:) = yld(:)
!          yud2(:) = yud(:)
!       endif
  
!       !     Definition of the Control Volume
!       !*****************************************************************
!       !! xld,xrd,yld,yud: limits of control volume (!!don't use cex and cey anymore!!)
!       do iv=1,nvol
!          ! ok for istret=0 (!!to do for istret=1!!)
!          icvlf(iv) = nint(xld2(iv)/dx)+1
!          icvrt(iv) = nint(xrd2(iv)/dx)+1
!          jcvlw(iv) = nint(yld2(iv)/dy)+1
!          jcvup(iv) = nint(yud2(iv)/dy)+1
  
!          icvlf_lx(iv) = icvlf(iv)
!          icvrt_lx(iv) = icvrt(iv)
!          jcvlw_lx(iv) = max(jcvlw(iv)+1-xstart(2),1)
!          jcvup_lx(iv) = min(jcvup(iv)+1-xstart(2),xsize(2))
!          jcvlw_lz(iv) = max(jcvlw(iv)+1-zstart(2),1)
!          jcvup_lz(iv) = min(jcvup(iv)+1-zstart(2),zsize(2))       
  
!          icvlf_ly(iv) = max(icvlf(iv)+1-ystart(1),1)
!          icvrt_ly(iv) = min(icvrt(iv)+1-ystart(1),ysize(1))
!          icvlf_lz(iv) = max(icvlf(iv)+1-zstart(1),1)
!          icvrt_lz(iv) = min(icvrt(iv)+1-zstart(1),zsize(1))       
!          jcvlw_ly(iv) = jcvlw(iv)
!          jcvup_ly(iv) = jcvup(iv)
  
!          zcvlf(iv) = nint(zld(iv)/dz)+1
!          zcvrt(iv) = nint(zrd(iv)/dz)+1
!          zcvlf_lx(iv) = max(zcvlf(iv)+1-xstart(3),1)
!          zcvrt_lx(iv) = min(zcvrt(iv)+1-xstart(3),xsize(3)) 
!          zcvlf_ly(iv) = max(zcvlf(iv)+1-ystart(3),1)
!          zcvrt_ly(iv) = min(zcvrt(iv)+1-ystart(3),ysize(3)) 
!       enddo
!     end subroutine update_forces
  
!   !***********************************************************************
  ! 
    subroutine restart_forces(itest1)
  ! 
  !***********************************************************************
  
      USE decomp_2d
      USE decomp_2d_io
      USE variables
      USE param
      USE MPI
  
      implicit none
  
      integer :: fh,ierror,code,itest1
      integer :: ierror_o=0 !error to open sauve file during restart
      character(len=30) :: filename, filestart
      integer (kind=MPI_OFFSET_KIND) :: filesize, disp
  
      write(filename, "('restart-forces',I7.7)") itime
      write(filestart,"('restart-forces',I7.7)") ifirst-1
  
      if (itest1==1) then !write
         if (mod(itime, icheckpoint).ne.0) then
            return
         endif
  
         call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,1,ux01)
         call decomp_2d_write_var(fh,disp,1,uy01)
         call decomp_2d_write_var(fh,disp,1,ux11)
         call decomp_2d_write_var(fh,disp,1,uy11)
         call decomp_2d_write_var(fh,disp,1,uz01)
         call decomp_2d_write_var(fh,disp,1,uz11)       
         call MPI_FILE_CLOSE(fh,ierror)
      else !read
         call MPI_FILE_OPEN(MPI_COMM_WORLD, filestart, &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror_o)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,1,ux01)
         call decomp_2d_read_var(fh,disp,1,uy01)
         call decomp_2d_read_var(fh,disp,1,ux11)
         call decomp_2d_read_var(fh,disp,1,uy11)
         call decomp_2d_read_var(fh,disp,1,uz01)
         call decomp_2d_read_var(fh,disp,1,uz11)       
         call MPI_FILE_CLOSE(fh,ierror_o)
      endif
  
      if (nrank.eq.0) then
         if (ierror_o .ne. 0) then !Included by Felipe Schuch
            print *,'==========================================================='
            print *,'Error: Impossible to read '//trim(filestart)
            print *,'==========================================================='
            call MPI_ABORT(MPI_COMM_WORLD,code,ierror)
         endif
      endif
  
    end subroutine restart_forces
  end module forces
  
  !***********************************************************************
  subroutine force(ux1,uy1,uz1,ep1)
  !***********************************************************************
  
    USE forces
    USE param
    USE variables
    USE decomp_2d
    USE MPI
    USE ibm
    use var, only : ta1, tb1, tc1, td1, te1, tf1, tg1, th1, ti1, di1
    use var, only : ux2, uy2, uz2, ta2, tb2, tc2, td2, te2, tf2, tg2, th2, ti2, di2
    use var, only : ux3, uy3, uz3, ta3, tb3, tc3, td3, te3, tf3, tg3, th3, ti3, di3
      
  
    implicit none
    character(len=30) :: filename, filename2
    integer :: nzmsize
    integer                                             :: i, iv, j, k, kk, code, jj
    integer                                             :: nvect1,nvect2,nvect3
  
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ep1
  
    real(mytype), dimension(ysize(1),ysize(2),ysize(3)) :: ppi2
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: ppi3
  
    real(mytype), dimension(nz) :: yLift, xDrag, zLat
    real(mytype) :: yLift_mean, xDrag_mean, zLat_mean
  
    real(mytype), dimension(nz) :: tunstxl, tunstyl, tunstzl
    real(mytype), dimension(nz) :: tconvxl, tconvyl, tconvzl
    real(mytype), dimension(nz) :: tpresxl, tpresyl
    real(mytype), dimension(nz) :: tdiffxl, tdiffyl, tdiffzl
  
    real(mytype), dimension(nz) :: tunstx, tunsty, tunstz  
    real(mytype), dimension(nz) :: tconvx, tconvy, tconvz
    real(mytype), dimension(nz) :: tpresx, tpresy 
    real(mytype), dimension(nz) :: tdiffx, tdiffy, tdiffz
    
    
    real(mytype), dimension(ny) :: tconvxl2, tconvyl2, tconvzl2
    real(mytype), dimension(ny) :: tdiffxl2, tdiffyl2, tdiffzl2
    real(mytype), dimension(ny) :: tconvx2, tconvy2, tconvz2
    real(mytype), dimension(ny) :: tdiffx2, tdiffy2, tdiffz2
    real(mytype), dimension(ny) :: tpreszl, tpresz
    
    
    
  
    real(mytype) :: uxmid, uymid, uzmid, prmid
    real(mytype) :: dudxmid, dudymid, dudzmid, dvdxmid, dvdymid, dvdzmid
    real(mytype) :: dwdxmid, dwdymid, dwdzmid
    real(mytype) :: fac,tsumx, tsumy, tsumz
    real(mytype) :: fcvx, fcvy, fcvz, fprx, fpry, fprz, fdix, fdiy, fdiz
    real(mytype) :: xmom, ymom, zmom
    real(mytype), dimension(ny) :: ztpresx, ztpresy
    real(mytype), dimension(nz) :: zyLift, zxDrag, zzLat
    real(mytype) :: zyLift_mean, zxDrag_mean, zzLat_mean
    
    
    
    real(mytype), dimension(nz) :: drag1, drag2, drag11, drag22
    real(mytype), dimension(nz) :: drag3, drag4, drag33, drag44
    real(mytype) :: mom1, mom2, mom3, tp1, tp2, tp3, dra1, dra2, dra3
  
  !  if (imove.eq.1) then
  !     ux1(:,:,:) = ux1(:,:,:) + 0.5
  !  endif
  
    nvect1=xsize(1)*xsize(2)*xsize(3)
    nvect2=ysize(1)*ysize(2)*ysize(3)
    nvect3=zsize(1)*zsize(2)*zsize(3)
  
    if (itime.eq.1) then
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
  !print*, t
  !if (nrank.eq.0) print*, ppi1
  
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)    ! du/dx
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy) ! dv/dx
    call derx (te1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz) ! dw/dx
      
    call transpose_x_to_y(ta1,ta2) ! du/dx
    call transpose_x_to_y(tb1,tb2) ! dv/dx
    call transpose_x_to_y(te1,te2) ! dw/dx
  
    call transpose_x_to_y(ux1,ux2)
    call transpose_x_to_y(uy1,uy2)
    call transpose_x_to_y(uz1,uz2)
    call transpose_x_to_y(ppi1,ppi2)
  
  
    call dery (tc2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx) ! du/dy
    call dery (td2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)    ! dv/dy
    call dery (tf2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz) ! dw/dy
  
    call transpose_y_to_z(ux2,ux3)
    call transpose_y_to_z(uy2,uy3)
    call transpose_y_to_z(uz2,uz3)
  !!!!!!!  call transpose_y_to_z(ppi2,ppi3)
    
  !!!!!!!  call transpose_y_to_z(te2,te3) ! dw/dx
  !!!!!!!  call transpose_y_to_z(tf2,tf3) ! dw/dy
    
    
    call derz (tg3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcx)  ! du/dz
    call derz (th3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcy)  ! dv/dz
    call derz (ti3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcz)     ! dw/dz
  
  
    call transpose_z_to_y(tg3,tg2) ! du/dz
    call transpose_z_to_y(th3,th2) ! dv/dz
    call transpose_z_to_y(ti3,ti2) ! 
    
    
    call transpose_y_to_x(tc2,tc1) ! du/dy
    call transpose_y_to_x(td2,td1) ! dv/dy
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
  !     do k=1,xsize(3)
       do k=zcvlf_lx(iv),zcvrt_lx(iv)
          tsumx=zero
          tsumy=zero
          tsumz=zero
          do j=jcvlw_lx(iv),jcvup_lx(iv)
             do i=icvlf_lx(iv),icvrt_lx(iv)
                !     The velocity time rate has to be relative to the cell center, 
                !     and not to the nodes, because, here, we have an integral 
                !     relative to the volume, and, therefore, this has a sense 
                !     of a "source".
                fac   = (onepfive*ux1(i,j,k)-two*ux01(i,j,k)+half*ux11(i,j,k))*(one-ep1(i,j,k))
                tsumx = tsumx+fac*dx*dy*dz/dt
  
                fac   = (onepfive*uy1(i,j,k)-two*uy01(i,j,k)+half*uy11(i,j,k))*(one-ep1(i,j,k))
                tsumy = tsumy+fac*dx*dy*dz/dt
                
                fac   = (onepfive*uz1(i,j,k)-two*uz01(i,j,k)+half*uz11(i,j,k))*(one-ep1(i,j,k))
                tsumz = tsumz+fac*dx*dy*dz/dt              
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
  !        do k=1,xsize(3)
          do k=zcvlf_lx(iv),zcvrt_lx(iv)
             kk=xstart(3)-1+k
             fcvx=zero
             fcvy=zero
             fcvz=zero
             fpry=zero
             fdix=zero
             fdiy=zero
             fdiz=zero
             do i=icvlf_lx(iv),icvrt_lx(iv)-1
                !momentum flux
                uxmid = half*(ux1(i,j,k)+ux1(i+1,j,k))
                uymid = half*(uy1(i,j,k)+uy1(i+1,j,k))
                uzmid = half*(uz1(i,j,k)+uz1(i+1,j,k))
  
                fcvx= fcvx -uxmid*uymid*dx*dz
                fcvy= fcvy -uymid*uymid*dx*dz
                fcvz= fcvz -uymid*uzmid*dx*dz
  
                !pressure
                prmid = half*(ppi1(i,j,k)+ppi1(i+1,j,k))
                fpry = fpry +prmid*dx*dz
  
                !viscous term
                dudymid = half*(tc1(i,j,k)+tc1(i+1,j,k))
                dvdxmid = half*(tb1(i,j,k)+tb1(i+1,j,k))
                dvdymid = half*(td1(i,j,k)+td1(i+1,j,k))
                dwdymid = half*(tf1(i,j,k)+tf1(i+1,j,k))
                dvdzmid = half*(th1(i,j,k)+th1(i+1,j,k))
                
                fdix = fdix -(xnu*(dudymid+dvdxmid)*dx*dz)
                fdiy = fdiy -two*xnu*dvdymid*dx*dz
                fdiz = fdiz -(xnu*(dwdymid+dvdzmid)*dx*dz)
  
             enddo
  !           drag1(kk)=drag1(kk)+fcvx
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
  !        do k=1,xsize(3)
          do k=zcvlf_lx(iv),zcvrt_lx(iv)
             kk=xstart(3)-1+k
             fcvx=zero
             fcvy=zero
             fcvz=zero
             fpry=zero
             fdix=zero
             fdiy=zero
             fdiz=zero
             do i=icvlf_lx(iv),icvrt_lx(iv)-1
                !momentum flux
                uxmid = half*(ux1(i,j,k)+ux1(i+1,j,k))
                uymid = half*(uy1(i,j,k)+uy1(i+1,j,k))
                uzmid = half*(uz1(i,j,k)+uz1(i+1,j,k))
                fcvx= fcvx +uxmid*uymid*dx*dz
                fcvy= fcvy +uymid*uymid*dx*dz
                fcvz= fcvz +uymid*uzmid*dx*dz
  
                !pressure
                prmid = half*(ppi1(i,j,k)+ppi1(i+1,j,k))
                fpry = fpry -prmid*dx*dz
  
                !viscous term
                dudymid = half*(tc1(i,j,k)+tc1(i+1,j,k))
                dvdxmid = half*(tb1(i,j,k)+tb1(i+1,j,k))
                dvdymid = half*(td1(i,j,k)+td1(i+1,j,k))
                dwdymid = half*(tf1(i,j,k)+tf1(i+1,j,k))
                dvdzmid = half*(th1(i,j,k)+th1(i+1,j,k))
  
                fdix = fdix +(xnu*(dudymid+dvdxmid)*dx*dz)
                fdiy = fdiy +two*xnu*dvdymid*dx*dz
                fdiz = fdiz +(xnu*(dwdymid+dvdzmid)*dx*dz)
  
             enddo
  !           drag2(kk)=drag2(kk)+fcvx
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
       
  !     print*,icvlf(iv) 
       if ((icvlf(iv).ge.ystart(1)).and.(icvlf(iv).le.yend(1))) then
          i=icvlf(iv)-ystart(1)+1
  !        do k=1,ysize(3)
          do k=zcvlf_ly(iv),zcvrt_ly(iv)
             kk=ystart(3)-1+k
             fcvx=zero
             fcvy=zero
             fcvz=zero
             fprx=zero
             fdix=zero
             fdiy=zero
             fdiz=zero
             do j=1,ysize(2)-1
  !           do j=jcvlw_ly(iv),jcvup_ly(iv)-1
                !momentum flux
                uxmid = half*(ux2(i,j,k)+ux2(i,j+1,k))
                uymid = half*(uy2(i,j,k)+uy2(i,j+1,k))
                uzmid = half*(uz2(i,j,k)+uz2(i,j+1,k))
                fcvx= fcvx -uxmid*uxmid*dy*dz
                fcvy= fcvy -uxmid*uymid*dy*dz
                fcvz= fcvz -uxmid*uzmid*dy*dz
  
                !pressure
                prmid=half*(ppi2(i,j,k)+ppi2(i,j+1,k))
                fprx = fprx +prmid*dy*dz
  
                !viscous term
                dudxmid = half*(ta2(i,j,k)+ta2(i,j+1,k))
                dudymid = half*(tc2(i,j,k)+tc2(i,j+1,k))
                dvdxmid = half*(tb2(i,j,k)+tb2(i,j+1,k))
                dwdxmid = half*(te2(i,j,k)+te2(i,j+1,k))
                dudzmid = half*(tg2(i,j,k)+tg2(i,j+1,k))
                
                fdix = fdix -two*xnu*dudxmid*dy*dz
                fdiy = fdiy -xnu*(dvdxmid+dudymid)*dy*dz
                fdiz = fdiz -xnu*(dwdxmid+dudzmid)*dy*dz
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
  !        do k=1,ysize(3)
          do k=zcvlf_ly(iv),zcvrt_ly(iv)
             kk=ystart(3)-1+k
             fcvx=zero
             fcvy=zero
             fcvz=zero
             fprx=zero
             fdix=zero
             fdiy=zero
             fdiz=zero
             do j=1,ysize(2)-1
  !           do j=jcvlw_ly(iv),jcvup_ly(iv)-1
                !momentum flux
                uxmid = half*(ux2(i,j,k)+ux2(i,j+1,k))
                uymid = half*(uy2(i,j,k)+uy2(i,j+1,k))
                uzmid = half*(uz2(i,j,k)+uz2(i,j+1,k))
                fcvx= fcvx +uxmid*uxmid*dy*dz
                fcvy= fcvy +uxmid*uymid*dy*dz
                fcvz= fcvz +uxmid*uzmid*dy*dz
  
                !pressure
                prmid=half*(ppi2(i,j,k)+ppi2(i,j+1,k))
                fprx = fprx -prmid*dy*dz
  
                !viscous term
                dudxmid = half*(ta2(i,j,k)+ta2(i,j+1,k))
                dudymid = half*(tc2(i,j,k)+tc2(i,j+1,k))
                dvdxmid = half*(tb2(i,j,k)+tb2(i,j+1,k))
                dwdxmid = half*(te2(i,j,k)+te2(i,j+1,k))
                dudzmid = half*(tg2(i,j,k)+tg2(i,j+1,k))
  
                fdix = fdix +two*xnu*dudxmid*dy*dz
                fdiy = fdiy +xnu*(dvdxmid+dudymid)*dy*dz
                fdiz = fdiz +xnu*(dwdxmid+dudzmid)*dy*dz
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
          
  
          fcvx=zero
          fcvy=zero
          fcvz=zero
          fprz=zero
          fdix=zero
          fdiy=zero
          fdiz=zero
          do j=jcvlw_lx(iv),jcvup_lx(iv)
           kk = xstart(2)-1+j
             do i=icvlf_lx(iv),icvrt_lx(iv)-1
                !momentum flux
                uxmid = half*(ux1(i,j,k)+ux1(i+1,j,k))
                uymid = half*(uy1(i,j,k)+uy1(i+1,j,k))
                uzmid = half*(uz1(i,j,k)+uz1(i+1,j,k))
  
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
           kk = xstart(2)-1+j 
            do i=icvlf_lx(iv),icvrt_lx(iv)-1
                !momentum flux
                uxmid = half*(ux1(i,j,k)+ux1(i+1,j,k))
                uymid = half*(uy1(i,j,k)+uy1(i+1,j,k))
                uzmid = half*(uz1(i,j,k)+uz1(i+1,j,k))
  
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
  
       
  !     call MPI_ALLREDUCE(drag1,drag11,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
  !     call MPI_ALLREDUCE(drag2,drag22,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
  !     call MPI_ALLREDUCE(drag3,drag33,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
  !     call MPI_ALLREDUCE(drag4,drag44,nz,real_type,MPI_SUM,MPI_COMM_WORLD,code)
                      
  
       tp1 = sum(tpresx(:))/dt
       tp2 = sum(tpresy(:))/dt
       tp3 = sum(tpresz(:))/dt
    
       mom1 = sum(tunstx(:) + tconvx(:) + tconvx2(:))
       mom2 = sum(tunsty(:) + tconvy(:) + tconvy2(:))
       mom3 = sum(tunstz(:) + tconvz(:) + tconvz2(:))
  
       dra1 = 2.0*(sum(tdiffx) + sum(tdiffx2) + tp1 - mom1)
       dra2 = 2.0*(sum(tdiffy) + sum(tdiffy2) + tp2 - mom2)
       dra3 = 2.0*(sum(tdiffz) + sum(tdiffz2) + tp3 - mom3)
  !print*, dra1, tp1, mom1
  !print*, dra1, mom1, tp1, sum(tpresx(:))/dt, sum(tunstx(:)), sum(tconvx(:)), sum(tdiffx)
  
  !     do k=zcvlf(iv),zcvrt(iv)
  !        tpresx(k)=tpresx(k)/dt
  !        tpresy(k)=tpresy(k)/dt
  !        tpresz(k)=tpresz(k)/dt
  
  !        xmom    = tunstx(k)+tconvx(k)
  !        ymom    = tunsty(k)+tconvy(k)
  !        zmom    = tunstz(k)+tconvz(k)
  !        xDrag(k) = two*(tdiffx(k)+tpresx(k)-xmom)
  !        yLift(k) = two*(tdiffy(k)+tpresy(k)-ymom)
  !        zLat(k) = two*(tdiffz(k)+tpresz(k)-zmom)
  !     enddo
  !     xDrag_mean = sum(xDrag(zcvlf(iv):zcvrt(iv)))/real(zcvrt(iv)-zcvlf(iv))
  !     yLift_mean = sum(yLift(zcvlf(iv):zcvrt(iv)))/real(zcvrt(iv)-zcvlf(iv))
  !     zLat_mean = sum(zLat(zcvlf(iv):zcvrt(iv)))/real(zcvrt(iv)-zcvlf(iv))
  
            
       if ((itime==ifirst).or.(itime==0)) then
          if (nrank .eq. 0) then
          write(filename,"('aerof',I1.1)") iv
          open(38+(iv-1),file=filename,status='unknown',form='formatted')
          endif
       endif
       if (nrank .eq. 0) then
          write(38+(iv-1),*) t, dra1, dra2, dra3!, sum(drag11), sum(drag22), sum(drag11)+sum(drag22), sum(drag33), sum(drag44), sum(drag33)+sum(drag44)  
  !        write(38+(iv-1),*) t, dra1, dra2, dra3
       endif
       if (itime==ilast) then
          if (nrank .eq. 0) then
             close(38+(iv-1))
             write(filename,"('aerof',I1.1)") iv
             write(filename2,"('aerof',I1.1,'-',I7.7)") iv, itime
             call system("mv " //filename //filename2)
          endif
       endif
       
       
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
  
  !  if (imove.eq.1) then
  !     ux1(:,:,:) = ux1(:,:,:) - 0.5
  !  endif
  
    return
  
  end subroutine force
  
  
  
  
  
  