!=======================================================================
! This program computes the drag and lift coefficients alongo a 
! cylinder by the control ! volume (CV) technique for 2D (pencil) 
! decomposition. 
!
! Adpated from Leandro Pinto PhD Thesis (2012) by Gabriel Narvaez Campo 
! 03-2017 Nucleo de Estudos em Transicao e Turbulencia (NETT/IPH/UFRGS)
!
!=======================================================================
module forces
  USE decomp_2d
  USE flow_type
  implicit none

  real(mytype),allocatable,dimension(:,:,:) :: ux03, uy03, ux13, uy13, epcv3
  real(mytype) :: xld,xrd,yld,yud

contains

  subroutine init_forces(ep1)

    USE decomp_2d
    USE param
    USE variables

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    real(mytype), dimension(ysize(1),ysize(2),ysize(3)) :: ep2
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: ep3
    integer :: i,j,k,icvlf,icvrt,jcvlw,jcvup
    call alloc_z(ux03);call alloc_z(uy03);call alloc_z(ux13);call alloc_z(uy13)
    call alloc_z(epcv3)

    !     Definition of the Control Volume
    !*****************************************************************
    !     icvlf, icvrt, jcvlw and jcvup are the counters that
    !     define the borders of the CV.

    call transpose_x_to_y(ep1,ep2)
    call transpose_y_to_z(ep2,ep3)
    if (istret.ne.0) then
       icvlf = nint((cex-xld)/dx)+1
       icvrt = nint((cex+xrd)/dx)+1
       do j=1,ny
          if (yp(j).le.(cey-yld)) then
             jcvlw = j
          else if (yp(j).le.(cey+yud)) then
             jcvup = j
          endif
       enddo
    else
       icvlf = nint((cex-xld)/dx)+1
       icvrt = nint((cex+xrd)/dx)+1
       jcvlw = nint((cey-yld)/dy)+1
       jcvup = nint((cey+yud)/dy)+1
    endif

    epcv3(:,:,:) = zero
    do i=1,zsize(1)
       if (i+zstart(1)-1 < icvlf) cycle
       if (i+zstart(1)-1 > icvrt) cycle
       do j=1,zsize(2)
          if (j+zstart(2)-1 < jcvlw) cycle
          if (j+zstart(2)-1 > jcvup) cycle
          do k=1,zsize(3)
             epcv3(i,j,k) = (one - ep3(i,j,k))
          enddo
       enddo
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
       call decomp_2d_write_var(fh,disp,3,ux03)
       call decomp_2d_write_var(fh,disp,3,uy03)
       call decomp_2d_write_var(fh,disp,3,ux13)
       call decomp_2d_write_var(fh,disp,3,uy13)
       call MPI_FILE_CLOSE(fh,ierror)
    else !read
       call MPI_FILE_OPEN(MPI_COMM_WORLD, filestart, &
            MPI_MODE_RDONLY, MPI_INFO_NULL, &
            fh, ierror_o)
       disp = 0_MPI_OFFSET_KIND
       call decomp_2d_read_var(fh,disp,3,ux03)
       call decomp_2d_read_var(fh,disp,3,uy03)
       call decomp_2d_read_var(fh,disp,3,ux13)
       call decomp_2d_read_var(fh,disp,3,uy13)
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
subroutine force(ux1,uy1,uz1,ux03,ux13,uy03,uy13,ep1,epcv3,pp3,nzmsize,phG,ph2,ph3)

  USE param
  USE variables
  USE decomp_2d
  USE MPI

  implicit none

  TYPE(DECOMP_INFO) :: phG,ph2,ph3
  integer :: nzmsize
  integer                                             :: i, j, k, ijk, code
  integer                                             :: nvect1,nvect2,nvect3
  real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
  real(mytype), dimension(ysize(1),ysize(2),ysize(3)) :: ux2, uy2, uz2, ep2
  real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: ux3, uy3, uz3, ep3, epcv3
  real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: ux03, uy03, ux13, uy13
  integer :: icvlf,icvrt,jcvlw,jcvup !counters that
  !      real(mytype),dimension(ny) :: yp
  !      real(mytype) :: xcil,ycil,xcv,ycv,vxcil,vycil

  !rmomentum
  real(mytype) :: fac,dudt1,dudt2,dudt3,dudt4,duxdt,duydt,tsumx, tsumy, tsumx1, tsumy1
  real(mytype),dimension(zsize(3)) :: sumx, sumy, sumx1, sumy1

  !Surface momentum fluxes
  real(mytype) :: fxab,fyab,fxbc,fybc,fxcd,fycd,fxda,fyda
  real(mytype) :: fxab1,fyab1,fxbc1,fybc1,fxcd1,fycd1,fxda1,fyda1
  real(mytype) :: pxab,pyab,pxbc,pybc,pxcd,pycd,pxda,pyda
  real(mytype) :: xmom,ymom,uxm,uym
  real(mytype),dimension(zsize(3)) :: tvcx,tvcy

  !External forces CV
  real(mytype),dimension(zsize(3)) :: f2xx,f2yy
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,di1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,di2
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3
  real(mytype) :: foxab, foyab, foxbc, foybc, foxcd, foycd, foxda, foyda
  real(mytype) :: foxab1,foyab1,foxbc1,foybc1,foxcd1,foycd1,foxda1,foyda1
  real(mytype) :: dudxm,dudym,dvdxm,dvdym
  real(mytype) :: pbc,pyp,foypr,foxpr,rho
  real(mytype) :: pab,pcd,difp,pxp,padx
  real(mytype) :: pab1,pcd1,padx1,pbc1

  ! Hydrodynamic coefficients
  real(mytype),dimension(zsize(3)) :: yLift,xDrag
  real(mytype) :: yLift_mean,xDrag_mean

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
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ppi3

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

  !X ==> Y
  call transpose_x_to_y(ux1,ux2)
  call transpose_x_to_y(uy1,uy2)
  call transpose_x_to_y(ep1,ep2)
  call transpose_x_to_y(ppi1,ppi2)
  !Y ==> Z
  call transpose_y_to_z(ux2,ux3)
  call transpose_y_to_z(uy2,uy3)
  call transpose_y_to_z(ep2,ep3)
  call transpose_y_to_z(ppi2,ppi3)

  nvect1=xsize(1)*xsize(2)*xsize(3)
  nvect2=ysize(1)*ysize(2)*ysize(3)
  nvect3=zsize(1)*zsize(2)*zsize(3)

  if (itime.eq.1)then
     do ijk=1,nvect3
        ux13(ijk,1,1)=ux3(ijk,1,1)
        uy13(ijk,1,1)=uy3(ijk,1,1)
     enddo
     return
  elseif (itime.eq.2)then
     do ijk=1,nvect3
        ux03(ijk,1,1)=ux3(ijk,1,1)
        uy03(ijk,1,1)=uy3(ijk,1,1)
     enddo
     return
  endif

  !*****************************************************************
  !      Drag and Lift coefficients
  !*****************************************************************

!!$  !     Definition of the Control Volume
!!$  !*****************************************************************
!!$  !     icvlf, icvrt, jcvlw and jcvup are the counters that
!!$  !     define the borders of the CV.
!!$
!!$  if (istret.ne.0) then
!!$     icvlf = nint((cex-xld)/dx)+1
!!$     icvrt = nint((cex+xrd)/dx)+1
!!$     do j=1,ny
!!$        if (yp(j).le.(cey-yld)) then
!!$           jcvlw = j
!!$        else if (yp(j).le.(cey+yud)) then
!!$           jcvup = j
!!$        endif
!!$     enddo
!!$  else
!!$     icvlf = nint((cex-xld)/dx)+1
!!$     icvrt = nint((cex+xrd)/dx)+1
!!$     jcvlw = nint((cey-yld)/dy)+1
!!$     jcvup = nint((cey+yud)/dy)+1
!!$  endif
!!$
!!$  epcv3(:,:,:) = 0.0
!!$  do i=1,zsize(1)
!!$     if (i+zstart(1)-1 < icvlf) cycle
!!$     if (i+zstart(1)-1 > icvrt) cycle
!!$     do j=1,zsize(2)
!!$        if (j+zstart(2)-1 < jcvlw) cycle
!!$        if (j+zstart(2)-1 > jcvup) cycle
!!$        do k=1,zsize(3)
!!$           epcv3(i,j,k) = (1. - ep3(i,j,k))
!!$        enddo
!!$     enddo
!!$  enddo

  !*****************************************************************
  !     Definition of parallelized variables
  !*****************************************************************

  call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
  call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  !X ==> Y
  !call transpose_x_to_y(ux1,ux2)
  !call transpose_x_to_y(uy1,uy2)
  !call transpose_x_to_y(uz1,uz2)
  !call transpose_x_to_y(ppi1,ppi2)
  call transpose_x_to_y(ta1,ta2) ! dudx
  call transpose_x_to_y(tb1,tb2) ! dvdx
  call dery (tc2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) 
  call dery (td2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0) 
  !Y ==> Z
  !call transpose_y_to_z(ux2,ux3)
  !call transpose_y_to_z(uy2,uy3)
  !call transpose_y_to_z(uz2,uz3)
  !call transpose_y_to_z(ppi2,ppi3)
  call transpose_y_to_z(ta2,ta3) ! dudx
  call transpose_y_to_z(tb2,tb3) ! dvdx
  call transpose_y_to_z(tc2,tc3) ! dudy
  call transpose_y_to_z(td2,td3) ! dvdy
  !call transpose_y_to_z(ep2,ep3)
  !call transpose_y_to_z(epcv2,epcv3)

  !Z PENCILS

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

  do k=1,zsize(3)
     tsumx=zero
     tsumy=zero
     sumx(k)=zero
     sumy(k)=zero
     do j=1,zsize(2)  !zstart(2),zend(2)
#ifdef STRETCHING
        dy = ppy(j-1+zstart(2))
#endif
        do i=1,zsize(1) !zstart(1),zend(1)
           !     The velocity time rate has to be relative to the cell center, 
           !     and not to the nodes, because, here, we have an integral 
           !     relative to the volume, and, therefore, this has a sense 
           !     of a "source".
           fac   = (onepfive*ux3(i,j,k)-two*ux03(i,j,k)+half*ux13(i,j,k))*epcv3(i,j,k)
           dudt1 = fac/dt
           tsumx = tsumx+dudt1*dx*dy
           !sumx(k) = sumx(k)+dudt1*dx*dy

           fac   = (onepfive*uy3(i,j,k)-two*uy03(i,j,k)+half*uy13(i,j,k))*epcv3(i,j,k)
           dudt1 = fac/dt
           tsumy = tsumy+dudt1*dx*dy
           !sumy(k) = sumy(k)+dudt1*dx*dy
        enddo
     enddo
     call MPI_ALLREDUCE(tsumx,tsumx1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(tsumy,tsumy1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
     sumx(k) = tsumx1
     sumy(k) = tsumy1
     !*********************************************************************************
     !     Secondly, the surface momentum fluxes
     !*********************************************************************************

     !       (icvlf)      (icvrt)
     !(jcvup) B____________C  
     !        \            \
     !        \     __     \
     !        \    \__\    \
     !        \            \
     !        \       CV   \
     !        A____________D
     !(jcvlw)                    
     !      
     fxab=zero
     fyab=zero
     fxbc=zero
     fybc=zero
     fxcd=zero
     fycd=zero
     fxda=zero
     fyda=zero
     do  j=2,zsize(2)-1 !j=zstart(2)+1, zend(2)-1
#ifdef STRETCHING
        dy = ppy(j-1+zstart(2))
#endif
        do  i=2,zsize(1)-1!i=zstart(1)+1, zend(1)-1
           !         AB face
           if ((epcv3(i,j,k).eq.1).and.(epcv3(i-1,j,k).eq.0).and.(epcv3(i,j+1,k).eq.1).and.(ep3(i,j,k).eq.0)) then
              uxm = (ux3(i,j,k)+ux3(i,j+1,k))*half*epcv3(i,j,k)
              uym = (uy3(i,j,k)+uy3(i,j+1,k))*half*epcv3(i,j,k)
              pxab=-uxm*uxm*dy
              fxab= fxab+pxab
              pyab=-uxm*uym*dy
              fyab= fyab+pyab
              !         CD face
           elseif ((epcv3(i,j,k).eq.1).and.(epcv3(i+1,j,k).eq.0).and.(epcv3(i,j+1,k).eq.1).and.(ep3(i,j,k).eq.0)) then
              uxm = (ux3(i,j,k)+ux3(i,j+1,k))*half*epcv3(i,j,k)
              uym = (uy3(i,j,k)+uy3(i,j+1,k))*half*epcv3(i,j,k)
              pxcd=+uxm*uxm*dy
              fxcd= fxcd+pxcd
              pycd=+uxm*uym*dy
              fycd= fycd+pycd
              !         DA face
           elseif ((epcv3(i,j,k).eq.1).and.(epcv3(i,j-1,k).eq.0).and.(epcv3(i+1,j,k).eq.1).and.(ep3(i,j,k).eq.0)) then
              uxm = (ux3(i,j,k)+ux3(i+1,j,k))*half*epcv3(i,j,k)
              uym = (uy3(i,j,k)+uy3(i+1,j,k))*half*epcv3(i,j,k)
              pxda=-uxm*uym*dx
              fxda= fxda+pxda
              pyda=-uym*uym*dx
              fyda= fyda+pyda
              !         BC face
           elseif ((epcv3(i,j,k).eq.1).and.(epcv3(i,j+1,k).eq.0).and.(epcv3(i+1,j,k).eq.1).and.(ep3(i,j,k).eq.0)) then
              uxm = (ux3(i,j,k)+ux3(i+1,j,k))*half*epcv3(i,j,k)
              uym = (uy3(i,j,k)+uy3(i+1,j,k))*half*epcv3(i,j,k)
              pxbc=+uxm*uym*dx
              fxbc= fxbc+pxbc
              pybc=+uym*uym*dx
              fybc= fybc+pybc
           endif
        enddo
     enddo




     call MPI_ALLREDUCE(fxbc,fxbc1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(fybc,fybc1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

     call MPI_ALLREDUCE(fxda,fxda1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(fyda,fyda1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

     call MPI_ALLREDUCE(fxab,fxab1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(fyab,fyab1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

     call MPI_ALLREDUCE(fxcd,fxcd1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(fycd,fycd1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

     !     The components of the total momentum

     tvcx(k) = fxab1+fxbc1+fxcd1+fxda1
     tvcy(k) = fyab1+fybc1+fycd1+fyda1
     xmom    = sumx(k)+fxab1+fxbc1+fxcd1+fxda1
     ymom    = sumy(k)+fyab1+fybc1+fycd1+fyda1

     !********************************************************************************************

     !*************************************************************
     !       Calculation of forces on the external CV surface
     !*************************************************************


     !      call transpose_y_to_z(ta2,ta3) ! dudx
     !      call transpose_y_to_z(tb2,tb3) ! dvdx
     !      call transpose_y_to_z(tc2,tc3) ! dudy
     !      call transpose_y_to_z(td2,td3) ! dvdy

     !     VISCOUSE FORCES CONTRIBUTION.
     !
     !        THE PRESSURE CONTRIBUTION. 
     !
     !        Note: We shall pretend that pressure is given by 
     !        a vector pp(i,j,1). In this instance, we may obtain
     !        the pressure difference between two given points in
     !        the field. In reality this pressure difference has
     !        to be drawn from somewhere in the code.
     !
     !
     !       (icvlf)      (icvrt)
     !(jcvup) B____________C  
     !        \            \
     !        \     __     \
     !        \    \__\    \
     !        \            \
     !        \       CV   \
     !        A____________D
     !(jcvlw)  


     rho  = one
     foxpr= zero
     foypr= zero
     pab1  = zero
     pcd1  = zero
     padx1 = zero
     pbc1  = zero


     foxab= zero
     foyab= zero
     foxbc= zero
     foybc= zero
     foxcd= zero
     foycd= zero
     foxda= zero
     foyda= zero

     do j=1,zsize(2)-1
#ifdef STRETCHING
        dy = ppy(j-1+zstart(2))
#endif
        do i=2,zsize(1)-1
           pab  = zero
           pcd  = zero
           !         Along CV's entrance face AB
           if ((epcv3(i,j,k).eq.1).and.(epcv3(i-1,j,k).eq.0).and.(epcv3(i,j+1,k).eq.1).and.(ep3(i-1,j,k).eq.0)) then
              pab  = (ppi3(i,j,k) + ppi3(i,j+1,k))*half*epcv3(i,j,k)
              !if (k==2) print *,'pab, i, j, k, nrank', pab, i, j, k, nrank

              !         Viscouse force
              dudxm = (ta3(i,j,k)+ta3(i,j+1,k))*half
              dudym = (tc3(i,j,k)+tc3(i,j+1,k))*half
              dvdxm = (tb3(i,j,k)+tb3(i,j+1,k))*half
              pxab  =-two*xnu*dudxm*dy
              foxab = foxab+pxab
              pyab  =-xnu*(dvdxm+dudym)*dy
              foyab = foyab+pyab

              !         Along CV's exit face CD
           elseif ((epcv3(i,j,k).eq.1).and.(epcv3(i+1,j,k).eq.0).and.(epcv3(i,j+1,k).eq.1).and.(ep3(i+1,j,k).eq.0)) then
              pcd  = (ppi3(i,j,k) + ppi3(i,j+1,k))*half*epcv3(i,j,k)
              !if (k==2) print *,'pcd, i, j, k, nrank', pcd, i, j, k, nrank

              !         Viscouse force
              dudxm = (ta3(i,j,k)+ta3(i,j+1,k))*half
              dudym = (tc3(i,j,k)+tc3(i,j+1,k))*half
              dvdxm = (tb3(i,j,k)+tb3(i,j+1,k))*half
              pxcd  =+two*xnu*dudxm*dy
              foxcd = foxcd+pxcd
              pycd  =+xnu*(dvdxm+dudym)*dy
              foycd = foycd+pycd
           endif
           pab1 = pab1 + pab
           pcd1 = pcd1 + pcd
        enddo
     enddo
     !if (k==2)   print *,'pab, nrank', pab1, nrank
     !if (k==2)   print *,'pcd, nrank', pcd1,  nrank
     !        The pressure force between planes AB and CD
     call MPI_ALLREDUCE(pab1,pab,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(pcd1,pcd,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
     difp  = (pab/dt-pcd/dt)/rho
     pxp   = difp*dy !a possible problem for stretching
     foxpr = foxpr+pxp





     do j=2,zsize(2)-1 
        do i=1,zsize(1)-1
           padx = zero
           pbc  = zero
           !         Along CV's lower face DA
           if ((epcv3(i,j,k).eq.1).and.(epcv3(i,j-1,k).eq.0).and.(epcv3(i+1,j,k).eq.1).and.((1-ep3(i,j-1,k)).eq.1)) then
              padx = (ppi3(i,j,k) + ppi3(i+1,j,k))/2.0*epcv3(i,j,k)
              !if (k==2) print *,'padx, i, j, k, nrank', padx, i, j, k, nrank

              !         Viscouse force
              dudym = (tc3(i,j,k)+tc3(i+1,j,k))*half
              dvdxm = (tb3(i,j,k)+tb3(i+1,j,k))*half
              dvdym = (td3(i,j,k)+td3(i+1,j,k))*half
              pxda  =-xnu*(dudym+dvdxm)*dx
              foxda = foxda+pxda
              pyda  =-two*xnu*dvdym*dx
              foyda = foyda+pyda

              !         Along CV's upper face BC
           elseif ((epcv3(i,j,k).eq.1).and.(epcv3(i,j+1,k).eq.0).and.(epcv3(i+1,j,k).eq.1).and.((one-ep3(i,j+1,k)).eq.one)) then
              pbc  = (ppi3(i,j,k) + ppi3(i+1,j,k))*half*epcv3(i,j,k)
              !if (k==2) print *,'pbc, i, j, k, nrank', pbc, i, j, k, nrank

              !         Viscouse force
              dudym = (tc3(i,j,k)+tc3(i+1,j,k))*half
              dvdxm = (tb3(i,j,k)+tb3(i+1,j,k))*half
              dvdym = (td3(i,j,k)+td3(i+1,j,k))*half
              pxbc  =+xnu*(dudym+dvdxm)*dx
              foxbc = foxbc+pxbc
              pybc  =+two*xnu*dvdym*dx
              foybc = foybc+pybc
           endif
           padx1 = padx1 + padx
           pbc1  = pbc1  + pbc
        enddo
     enddo
     !if (k==2)   print *,'pad, nrank', padx1, nrank
     !if (k==2)   print *,'pbc, nrank', pbc1,  nrank
     !        The pressure force between planes AD and BC
     call MPI_ALLREDUCE(padx1,padx,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(pbc1,pbc,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
     difp  = (padx/dt-pbc/dt)/rho
     pyp   = difp*dx
     foypr = foypr+pyp
     !if (k==2)   print *,'difp, pyp, foypr, nrank', difp, pyp, foypr,  nrank

     !if (k==1) print *,'foypr', foypr
     call MPI_ALLREDUCE(foxab,foxab1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(foyab,foyab1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

     call MPI_ALLREDUCE(foxcd,foxcd1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(foycd,foycd1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

     call MPI_ALLREDUCE(foxbc,foxbc1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(foybc,foybc1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

     call MPI_ALLREDUCE(foxda,foxda1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(foyda,foyda1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)






     !        The resultant components along x and y
     !
     f2xx(k)=foxab1+foxbc1+foxcd1+foxda1 + foxpr
     f2yy(k)=foyab1+foybc1+foycd1+foyda1 + foypr

     !       Calculation of the aerodynamic force components
     !
     !
     !       Adjusting length scales reference. Multiply by the 
     !       airfoil maximum thickness and divide by chord. One 
     !       should observe that, in fact, calculated values   
     !       of lift and drag correspond already to the coefficients,
     !       because the program works with nondimensonalized 
     !       variables. The point here is that the original reference
     !       for lengths is the maximum thickness, and the force 
     !       coefficients are obtained considering as reference 
     !       the airfoil chord.


     xDrag(k) = two*(f2xx(k)-xmom)
     yLift(k) = two*(f2yy(k)-ymom) !Ja foi dividido por dz em rmomentum e sforce 

     !if (nrank==0) print *,'xDrag, yLift', 2*f2xx(k), 2*f2yy(k), 2*xmom, 2*ymom, xDrag(k), yLift(k)

     !**************************************************************************************************


  enddo

!!!
  ! Adições por Felipe Schuch
!!!
  xDrag_mean = sum(xDrag)/real(nz,mytype)
  yLift_mean = sum(yLift)/real(nz,mytype)
  if (nrank .eq. 0) then
     open(67,file='./out/aerof_avr',status='unknown',form='formatted',access='direct',recl=43) !43 = 3*14+1
     !Utilizando o acesso direto, cada valor para os coeficientes será escrito na linha itime-2,
     !eliminando problemas com possível restart
     write(67,'(3E14.6,A)',rec=itime-2) t,&                     !1
          xDrag_mean,&                                          !2
          yLift_mean,&                                          !3
          char(10) !new line character                          !+1
     close(67)
  endif
!!!
  !
!!!

  do ijk=1,nvect3
     ux13(ijk,1,1)=ux03(ijk,1,1)
     uy13(ijk,1,1)=uy03(ijk,1,1)
  enddo
  !
  do ijk=1,nvect3
     ux03(ijk,1,1)=ux3(ijk,1,1)
     uy03(ijk,1,1)=uy3(ijk,1,1)
  enddo

  return

end subroutine force



