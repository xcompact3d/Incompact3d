!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!!        FILE: BC-Periodic-hill.f90
!!!      AUTHOR: Sylvain Laizet
!!!    MODIFIED: Paul Bartholomew
!!! DESCRIPTION: This module describes the periodic hill flow.
!!!   CHANGELOG: [2019-02-19] Making module private by default
!!               [2019-02-19] Turning file into a module
!!               [2023-06-01] Implementing ready-to-run input file
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module hill

  USE decomp_2d_constants
  USE decomp_2d_mpi
  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_hill, boundary_conditions_hill, postprocess_hill, geomcomplex_hill, &
       visu_hill, visu_hill_init

contains

!############################################################################
  subroutine geomcomplex_hill(epsi,nxi,nxf,ny,nyi,nyf,nzi,nzf,dx,yp,dz,remp)
!############################################################################

    use param, only : zero, one, two, three, nine, fourteen, twenty, twentyeight
    use ibm

    implicit none

    integer                    :: nxi,nxf,ny,nyi,nyf,nzi,nzf
    real(mytype)               :: dx,dz
    real(mytype)               :: remp
    integer                    :: i,ic,j,k
    real(mytype)               :: xm,ym,r
    real(mytype)               :: zeromach
    real(mytype), dimension(nxi:nxf) :: dune
    real(mytype) :: y_bump
    real(mytype),dimension(nxi:nxf,nyi:nyf,nzi:nzf) :: epsi
    real(mytype),dimension(ny) :: yp

    zeromach=one
    do while ((one + zeromach / two) .gt. one)
       zeromach = zeromach/two
    end do
    zeromach = ten*zeromach
    !
    y_bump=zero
    dune=zero
    do i=nxi,nxf
       xm=real(i-1,mytype)*dx
       if (xm.gt.xlx/two) then
          xm = (xlx-xm)*twentyeight
       else
          xm = xm*twentyeight
       endif
       if ((xm >= zero).and.(xm<nine)) then
          y_bump=min(twentyeight,twentyeight+0.006775070969851_mytype*xm**two-2.124527775800E-03_mytype*xm**three)
       endif
       if ((xm >= nine).and.(xm<fourteen)) then
          y_bump=  2.507355893131E+01_mytype         +9.754803562315E-01_mytype*xm&
                  -1.016116352781E-01_mytype*xm**two +1.889794677828E-03_mytype*xm**three
       endif
       if ((xm >= fourteen).and.(xm<twenty)) then
          y_bump=  2.579601052357E+01_mytype         +8.206693007457E-01_mytype*xm &
                  -9.055370274339E-02_mytype*xm**two +1.626510569859E-03_mytype*xm**three
       endif
       if ((xm >= twenty).and.(xm<thirty)) then
          y_bump= 4.046435022819E+01_mytype         -1.379581654948E+00_mytype*xm &
                 +1.945884504128E-02_mytype*xm**two -2.070318932190E-04_mytype*xm**three
       endif
       if ((xm >= thirty).and.(xm<forty)) then
          y_bump= 1.792461334664E+01_mytype         +8.743920332081E-01_mytype*xm &
                 -5.567361123058E-02_mytype*xm**two +6.277731764683E-04_mytype*xm**three
       endif
       if ((xm >= forty).and.(xm <= fiftyfour)) then
          y_bump=max(zero,5.639011190988E+01_mytype          -2.010520359035E+00_mytype*xm &
                         +1.644919857549E-02_mytype*xm**two  +2.674976141766E-05_mytype*xm**three)
       endif
       dune(i)=y_bump/twentyeight
    enddo

    do k=nzi,nzf
       do j=nyi,nyf
          ym=yp(j)
          do i=nxi,nxf
             if (ym-dune(i).le.zeromach) then
                epsi(i,j,k)=remp
             endif
          enddo
       enddo
    enddo

    return
  end subroutine geomcomplex_hill

!############################################################################
  subroutine boundary_conditions_hill (ux,uy,uz,phi,ep1)
!############################################################################
  

    USE param
    USE variables

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: gx

    ux = ux*(one-ep1)
    call transpose_x_to_y(ux,gx)
    call hill_flrt(gx,(two/three)*(two/yly))
    call transpose_y_to_x(gx,ux)

    return
  end subroutine boundary_conditions_hill

!############################################################################
  subroutine init_hill (ux1,uy1,uz1,ep1,phi1)
!############################################################################

    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: y,r,um,r3,x,z,h,ct
    real(mytype) :: cx0,cy0,cz0,hg,lg
    integer :: k,j,i,ierror,ii,is,code

    integer, dimension (:), allocatable :: seed

    if (iscalar==1) then

       phi1(:,:,:,:) = zero !change as much as you want

    endif
    ux1=zero;uy1=zero;uz1=zero
    if (iin.ne.0) then
       call system_clock(count=code)
       if (iin.eq.2) code=0
       call random_seed(size = ii)
       call random_seed(put = code+63946*(nrank+1)*(/ (i - 1, i = 1, ii) /))

       call random_number(ux1)
       call random_number(uy1)
       call random_number(uz1)
    endif

    !modulation of the random noise
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             ux1(i,j,k)=init_noise*(two*ux1(i,j,k)-one)
             uy1(i,j,k)=init_noise*(two*uy1(i,j,k)-one)
             uz1(i,j,k)=init_noise*(two*uz1(i,j,k)-one)
          enddo
       enddo
    enddo

    !initial velocity profile
    do k=1,xsize(3)
       do j=1,xsize(2)
          if (istret.eq.0) y=real(j+xstart(2)-1-1,mytype)*dy-yly*half
          if (istret.ne.0) y=yp(j+xstart(2)-1)-yly*half
          do i=1,xsize(1)
             ux1(i,j,k)=ux1(i,j,k)+one-y*y
             uy1(i,j,k)=0.
             uz1(i,j,k)=0.
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
    if (nrank  ==  0) write(*,*) '# init end ok'
#endif

    return
  end subroutine init_hill

!############################################################################
  subroutine init_post(ep1)
!############################################################################

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ep1

  end subroutine init_post


!############################################################################
  subroutine hill_flrt (ux,constant)
!############################################################################

    USE decomp_2d_poisson
    USE variables
    USE param
    USE var
    USE MPI

    implicit none

    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux
    real(mytype) :: constant

    integer :: j,i,k,code
    real(mytype) :: can,ut3,ut

    ut3=zero
    do k=1,ysize(3)
       do i=1,ysize(1)
          ut=zero
          do j=1,ny-1
             if (istret.eq.0) then
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

    call MPI_ALLREDUCE(MPI_IN_PLACE,ut3,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    
    can=-(constant-ut3)

    do k=1,ysize(3)
       do i=1,ysize(1)
          do j=2,ny-1
             ux(i,j,k)=ux(i,j,k)-can
          enddo
       enddo
    enddo

    return
  end subroutine hill_flrt
  
!############################################################################
  subroutine postprocess_hill(ux1,uy1,uz1,pp3,phi1,ep1)
!############################################################################
    use var, ONLY : nzmsize

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3

  end subroutine postprocess_hill
  subroutine visu_hill_init(visu_initialised)

    use decomp_2d_io, only : decomp_2d_register_variable
    use visu, only : io_name, output2D
    
    implicit none

    logical, intent(out) :: visu_initialised

    call decomp_2d_register_variable(io_name, "critq", 1, 0, output2D, mytype)

    visu_initialised = .true.
    
  end subroutine visu_hill_init
  
!############################################################################
  subroutine visu_hill(ux1, uy1, uz1, pp3, phi1, ep1, num)
!############################################################################

    use var, only : ux2, uy2, uz2, ux3, uy3, uz3
    use var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    use var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    use var, ONLY : nzmsize
    use visu, only : write_field
    
    use ibm_param, only : ubcx,ubcy,ubcz

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    integer, intent(in) :: num

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

    !Q=-0.5*(ta1**2+te1**2+di1**2)-td1*tb1-tg1*tc1-th1*tf1
    di1 = zero
    di1(:,:,:) = - half*(ta1(:,:,:)**2 + te1(:,:,:)**2 + ti1(:,:,:)**2) &
                 - td1(:,:,:) * tb1(:,:,:) &
                 - tg1(:,:,:) * tc1(:,:,:) &
                 - th1(:,:,:) * tf1(:,:,:)
    call write_field(di1, ".", "critq", num, flush = .true.) ! Reusing temporary array, force flush

  end subroutine visu_hill
  
end module hill

