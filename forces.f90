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
  USE flow_type
  implicit none

  integer :: nvol
  real(mytype),allocatable,dimension(:,:,:) :: ux01, uy01, ux11, uy11
  real(mytype),allocatable,dimension(:) :: xld,xrd,yld,yud
  integer,allocatable,dimension(:) :: icvlf,icvrt,jcvlw,jcvup
  integer,allocatable,dimension(:) :: icvlf_lx,icvrt_lx,icvlf_ly,icvrt_ly
  integer,allocatable,dimension(:) :: jcvlw_lx,jcvup_lx,jcvlw_ly,jcvup_ly

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
       jcvlw(iv) = nint(yld(iv)/dy)+1
       jcvup(iv) = nint(yud(iv)/dy)+1
       icvlf_lx(iv) = icvlf(iv)
       icvrt_lx(iv) = icvrt(iv)
       jcvlw_lx(iv) = max(jcvlw(iv)+1-xstart(2),1)
       jcvup_lx(iv) = min(jcvup(iv)+1-xstart(2),xsize(2))
       icvlf_ly(iv) = max(icvlf(iv)+1-ystart(1),1)
       icvrt_ly(iv) = min(icvrt(iv)+1-ystart(1),ysize(1))
       jcvlw_ly(iv) = jcvlw(iv)
       jcvup_ly(iv) = jcvup(iv)
    enddo
    
  end subroutine init_forces

  subroutine restart_forces(irestart)

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none

    integer :: irestart,fh,ierror,code
    integer :: ierror_o=0 !error to open sauve file during restart
    character(len=30) :: filename, filestart
    integer (kind=MPI_OFFSET_KIND) :: filesize, disp

    write(filename, "('sauve-forces',I7.7)") itime
    write(filestart,"('sauve-forces',I7.7)") ifirst-1

    if (irestart==1) then !write
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
subroutine force(ux1,uy1,ep1,pp3,nzmsize,phG,ph2,ph3)

  !***********************************************************************

  USE forces
  USE param
  USE variables
  USE decomp_2d
  USE MPI

  implicit none
  TYPE(DECOMP_INFO) :: phG,ph2,ph3
  character(len=30) :: filename
  integer :: nzmsize
  integer                                             :: i, iv, j, k, kk, ijk, code
  integer                                             :: nvect1,nvect2,nvect3

  real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1
  real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
  real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,di1
  real(mytype), dimension(ysize(1),ysize(2),ysize(3)) :: ux2, uy2
  real(mytype), dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,di2

  !Pressure
  real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzm) :: pp3  
  !Z PENCILS NXM NYM NZM-->NXM NYM NZ
  real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)) :: tta3,ddi3
  !Y PENCILS NXM NYM NZ -->NXM NY NZ
  real(mytype),dimension(ph3%yst(1):ph3%yen(1),nym,ysize(3)) :: tta2
  real(mytype),dimension(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)) :: ttb2,ddi2
  !X PENCILS NXM NY NZ  -->NX NY NZ
  real(mytype),dimension(nxm,xsize(2),xsize(3)) :: tta1
  !NX NY NZ X --> Y --> Z
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ddi1,ppi1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ppi2

  real(mytype), dimension(nz) :: yLift,xDrag
  real(mytype) :: yLift_mean,xDrag_mean

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

  if (itime.eq.1) then
     do ijk=1,nvect1
        ux11(ijk,1,1)=ux1(ijk,1,1)
        uy11(ijk,1,1)=uy1(ijk,1,1)
     enddo
     return
  elseif (itime.eq.2) then
     do ijk=1,nvect1
        ux01(ijk,1,1)=ux1(ijk,1,1)
        uy01(ijk,1,1)=uy1(ijk,1,1)
     enddo
     return
  endif

  !WORK Z-PENCILS
  call interiz6(tta3,pp3,ddi3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
       (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
  !WORK Y-PENCILS
  call transpose_z_to_y(tta3,tta2,ph3) !nxm nym nz
  call interiy6(ttb2,tta2,ddi2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
       (ph3%yen(1)-ph3%yst(1)+1),nym,ysize(2),ysize(3),1)
  !WORK X-PENCILS
  call transpose_y_to_x(ttb2,tta1,ph2) !nxm ny nz
  call interi6(ppi1,tta1,ddi1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
       nxm,xsize(1),xsize(2),xsize(3),1)

  call derx (ta1,ux1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1) ! dudx
  call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1) ! dvdx
  call transpose_x_to_y(ta1,ta2) ! dudx
  call transpose_x_to_y(tb1,tb2) ! dvdx

  call transpose_x_to_y(ux1,ux2)
  call transpose_x_to_y(uy1,uy2)
  call transpose_x_to_y(ppi1,ppi2)

  call dery (tc2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) ! dudy
  call dery (td2,uy2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) ! dvdy
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

     tunstxl(:)=0.
     tunstyl(:)=0.
     do k=1,xsize(3)
        tsumx=0.0
        tsumy=0.0
        do j=jcvlw_lx(iv),jcvup_lx(iv)
           do i=icvlf_lx(iv),icvrt_lx(iv)
              !     The velocity time rate has to be relative to the cell center, 
              !     and not to the nodes, because, here, we have an integral 
              !     relative to the volume, and, therefore, this has a sense 
              !     of a "source".
              !         fac   = (1.5*ux1(i,j,k)-2.0*ux01(i,j,k)+0.5*ux11(i,j,k))*epcv1(i,j,k)
              fac   = (1.5*ux1(i,j,k)-2.0*ux01(i,j,k)+0.5*ux11(i,j,k))*(1.-ep1(i,j,k))
              tsumx = tsumx+fac*dx*dy/dt
              !sumx(k) = sumx(k)+dudt1*dx*dy

              !         fac   = (1.5*uy1(i,j,k)-2.0*uy01(i,j,k)+0.5*uy11(i,j,k))*epcv1(i,j,k)
              fac   = (1.5*uy1(i,j,k)-2.0*uy01(i,j,k)+0.5*uy11(i,j,k))*(1.-ep1(i,j,k))
              tsumy = tsumy+fac*dx*dy/dt
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
!!$!        A____________D
!!$!(jcvlw)                    
!!$!      

     tconvxl(:)=0.
     tconvyl(:)=0.
     tdiffxl(:)=0.
     tdiffyl(:)=0.
     tpresxl(:)=0.
     tpresyl(:)=0.
     !BC and AD : x-pencils
     !AD
     if ((jcvlw(iv).ge.xstart(2)).and.(jcvlw(iv).le.xend(2))) then
        j=jcvlw(iv)-xstart(2)+1
        do k=1,xsize(3)
           kk=xstart(3)-1+k
           fcvx=0.0
           fcvy=0.0
           fpry=0.0
           fdix=0.0
           fdiy=0.0
           do i=icvlf_lx(iv),icvrt_lx(iv)-1
              !momentum flux
              uxmid = 0.5*(ux1(i,j,k)+ux1(i+1,j,k))
              uymid = 0.5*(uy1(i,j,k)+uy1(i+1,j,k))
              fcvx= fcvx -uxmid*uymid*dx
              fcvy= fcvy -uymid*uymid*dx


              !pressure
              prmid = 0.5*(ppi1(i,j,k)+ppi1(i+1,j,k))
              fpry = fpry +prmid*dx

              !viscous term
              dudymid = 0.5*(tc1(i,j,k)+tc1(i+1,j,k))
              dvdxmid = 0.5*(tb1(i,j,k)+tb1(i+1,j,k))
              dvdymid = 0.5*(td1(i,j,k)+td1(i+1,j,k))
              fdix = fdix -(xnu*(dudymid+dvdxmid)*dx)
              fdiy = fdiy -2.0*xnu*dvdymid*dx

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
           fcvx=0.0
           fcvy=0.0
           fpry=0.0
           fdix=0.0
           fdiy=0.0
           do i=icvlf_lx(iv),icvrt_lx(iv)-1
              !momentum flux
              uxmid = 0.5*(ux1(i,j,k)+ux1(i+1,j,k))
              uymid = 0.5*(uy1(i,j,k)+uy1(i+1,j,k))
              fcvx= fcvx +uxmid*uymid*dx
              fcvy= fcvy +uymid*uymid*dx

              !pressure
              prmid = 0.5*(ppi1(i,j,k)+ppi1(i+1,j,k))
              fpry = fpry -prmid*dx

              !viscous term
              dudymid = 0.5*(tc1(i,j,k)+tc1(i+1,j,k))
              dvdxmid = 0.5*(tb1(i,j,k)+tb1(i+1,j,k))
              dvdymid = 0.5*(td1(i,j,k)+td1(i+1,j,k))
              fdix = fdix +(xnu*(dudymid+dvdxmid)*dx)
              fdiy = fdiy +2.0*xnu*dvdymid*dx

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
           fcvx=0.0
           fcvy=0.0
           fprx=0.0
           fdix=0.0
           fdiy=0.0
           do j=jcvlw_ly(iv),jcvup_ly(iv)-1
              !momentum flux
              uxmid = 0.5*(ux2(i,j,k)+ux2(i,j+1,k))
              uymid = 0.5*(uy2(i,j,k)+uy2(i,j+1,k))
              fcvx= fcvx -uxmid*uxmid*dy
              fcvy= fcvy -uxmid*uymid*dy

              !pressure
              prmid=0.5*(ppi2(i,j,k)+ppi2(i,j+1,k))
              fprx = fprx +prmid*dy

              !viscous term
              dudxmid = 0.5*(ta2(i,j,k)+ta2(i,j+1,k))
              dudymid = 0.5*(tc2(i,j,k)+tc2(i,j+1,k))
              dvdxmid = 0.5*(tb2(i,j,k)+tb2(i,j+1,k))
              fdix = fdix -2.0*xnu*dudxmid*dy
              fdiy = fdiy -xnu*(dvdxmid+dudymid)*dy
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
           fcvx=0.0
           fcvy=0.0
           fprx=0.0
           fdix=0.0
           fdiy=0.0
           do j=jcvlw_ly(iv),jcvup_ly(iv)-1
              !momentum flux
              uxmid = 0.5*(ux2(i,j,k)+ux2(i,j+1,k))
              uymid = 0.5*(uy2(i,j,k)+uy2(i,j+1,k))
              fcvx= fcvx +uxmid*uxmid*dy
              fcvy= fcvy +uxmid*uymid*dy

              !pressure
              prmid=0.5*(ppi2(i,j,k)+ppi2(i,j+1,k))
              fprx = fprx -prmid*dy

              !viscous term
              dudxmid = 0.5*(ta2(i,j,k)+ta2(i,j+1,k))
              dudymid = 0.5*(tc2(i,j,k)+tc2(i,j+1,k))
              dvdxmid = 0.5*(tb2(i,j,k)+tb2(i,j+1,k))
              fdix = fdix +2.0*xnu*dudxmid*dy
              fdiy = fdiy +xnu*(dvdxmid+dudymid)*dy
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
        xDrag(k) = 2.*(tdiffx(k)+tpresx(k)-xmom)
        yLift(k) = 2.*(tdiffy(k)+tpresy(k)-ymom)

        !if (nrank==0) print *,'xDrag, yLift', xDrag(k), yLift(k)!, icvlf, icvrt, jcvlw, jcvup
     enddo
     !Edited by F. Schuch
     xDrag_mean = sum(xDrag(:))/real(nz,mytype)
     yLift_mean = sum(yLift(:))/real(nz,mytype)
     if (nrank .eq. 0) then
        write(filename,"('./out/aerof_avr',I1.1)") iv
        open(67,file=filename,status='unknown',form='formatted',access='direct',recl=43) !43 = 3*14+1
        !Using the direct access, each value for the coefficients will be written
        !in the line itime-2, eliminating any problems with possible restart
        write(67,'(3E14.6,A)',rec=itime-2) t,&                     !1
             xDrag_mean,&                                          !2
             yLift_mean,&                                          !3
             char(10) !new line character                          !+1
        close(67)
     endif
  enddo

  do ijk=1,nvect1
     ux11(ijk,1,1)=ux01(ijk,1,1)
     uy11(ijk,1,1)=uy01(ijk,1,1)
  enddo

  do ijk=1,nvect1
     ux01(ijk,1,1)=ux1(ijk,1,1)
     uy01(ijk,1,1)=uy1(ijk,1,1)
  enddo

  return

end subroutine force
