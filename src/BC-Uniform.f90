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

module uniform

contains

  !*******************************************************************************
  !
  subroutine init_uniform (ux1,uy1,uz1,ep1,phi1)
  !
  !*******************************************************************************

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    integer :: i,j,k,code,ii
    integer, dimension (:), allocatable :: seed
    real(mytype) :: um

    ux1=zero; uy1=zero; uz1=zero;
    if (iscalar.eq.1) then
      phi1(:,:,:,:) = zero
    endif

    if (iin.ne.0) then
      call system_clock(count=code)
      if (iin.eq.2) code=0
      call random_seed(size = ii)
      call random_seed(put = code+63946*(nrank+1)*(/ (i - 1, i = 1, ii) /))

      call random_number(ux1)
      call random_number(uy1)
      call random_number(uz1)

      um=0.5*(u1+u2)
      do k=1,xsize(3)
        do j=1,xsize(2)
          do i=1,xsize(1)
            ux1(i,j,k)=init_noise*um*(two*ux1(i,j,k)-one)
            uy1(i,j,k)=init_noise*um*(two*uy1(i,j,k)-one)
            uz1(i,j,k)=init_noise*um*(two*uz1(i,j,k)-one)
          enddo
        enddo
      enddo
    endif
            
    do k=1,xsize(3)
    do j=1,xsize(2)
    do i=1,xsize(1)
      ux1(i,j,k)=ux1(i,j,k)+um
      uy1(i,j,k)=uy1(i,j,k)+zero
      uz1(i,j,k)=uz1(i,j,k)+zero
    enddo
    enddo
    enddo

    return
  end subroutine init_uniform

  !*******************************************************************************
  !
  subroutine boundary_conditions_uniform (ux,uy,uz,phi)
  !
  !*******************************************************************************

    USE param
    USE variables
    USE decomp_2d
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    if (nclx1.eq.2) then
      if (iscalar.eq.0.or.(iscalar.eq.1.and.nclxS1.eq.2)) then
        call inflow(ux,uy,uz,phi)
      endif
    endif

    if (nclxn.eq.2) then
      if (iscalar.eq.0.or.(iscalar.eq.1.and.nclxSn.eq.2)) then
        call outflow(ux,uy,uz,phi)
      endif
    endif

    return
  end subroutine boundary_conditions_uniform

  !*******************************************************************************
  !
  subroutine inflow (ux,uy,uz,phi)
  !
  !*******************************************************************************
  
    USE param
    USE variables
    USE decomp_2d
    USE MPI
    USE var, only: ux_inflow, uy_inflow, uz_inflow
  
    implicit none
  
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    real(mytype) :: um
    integer :: i,j,k,itime_input
    
    um=0.5*(u1+u2)
    do k=1,xsize(3)
    do j=1,xsize(2)
      bxx1(j,k)=um
      bxy1(j,k)=zero
      bxz1(j,k)=zero
    enddo
    enddo
 
    if (iin.eq.1.or.iin.eq.2) then
      call random_number(bxo)
      call random_number(byo)
      call random_number(bzo)

      do k=1,xsize(3)
      do j=1,xsize(2)
        bxx1(j,k)=bxx1(j,k)+(two*bxo(j,k)-one)*inflow_noise*um
        bxy1(j,k)=bxy1(j,k)+(two*byo(j,k)-one)*inflow_noise*um
        bxz1(j,k)=bxz1(j,k)+(two*bzo(j,k)-one)*inflow_noise*um
        if (iscalar.eq.1) then
          phi(1,j,k,:)=one
        endif
      enddo
      enddo
    else if (iin.eq.3) then
      ! Reading from files (when precursor simulations exist)
      itime_input=mod(itime,ntimesteps)
      if (itime_input==0) itime_input=ntimesteps
      if (nrank==0) print *,'Reading inflow from a file, time step: ', itime_input
      do k=1,xsize(3)
      do j=1,xsize(2)
        ! Case 1: Inflow is turbulence added to mean flow profile
        !bxx1(j,k)=bxx1(j,k)+ux_inflow(itime_input,j,k)
        !bxy1(j,k)=bxy1(j,k)+uy_inflow(itime_input,j,k)
        !bxz1(j,k)=bxz1(j,k)+uz_inflow(itime_input,j,k)
        ! Case 2: Inflow is full velocity field
        bxx1(j,k)=ux_inflow(itime_input,j,k)
        bxy1(j,k)=uy_inflow(itime_input,j,k)
        bxz1(j,k)=uz_inflow(itime_input,j,k)
      enddo
      enddo
    endif

    return
  end subroutine inflow 
    
  !*******************************************************************************
  !
  subroutine outflow (ux,uy,uz,phi)
  !
  !*******************************************************************************

    USE param
    USE variables
    USE decomp_2d
    USE MPI

    implicit none

    integer :: j,k,code
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx,uxmin,uxmax,uxmin1,uxmax1

    udx=one/dx; udy=one/dy; udz=one/dz; uddx=half/dx; uddy=half/dy; uddz=half/dz

    uxmax=-1609.
    uxmin=1609.
    do k=1,xsize(3)
      do j=1,xsize(2)
        if (ux(nx-1,j,k).gt.uxmax) uxmax=ux(nx-1,j,k)
        if (ux(nx-1,j,k).lt.uxmin) uxmin=ux(nx-1,j,k)
      enddo
    enddo

    call MPI_ALLREDUCE(uxmax,uxmax1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(uxmin,uxmin1,1,real_type,MPI_MIN,MPI_COMM_WORLD,code)

    cx=0.5*(uxmax1+uxmin1)*gdt(itr)*udx
    do k=1,xsize(3)
      do j=1,xsize(2)
        bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
        bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
        bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
        if (iscalar.eq.1) then
          phi(nx,j,k,:)=phi(nx,j,k,:)-cx*(phi(nx,j,k,:)-phi(nx-1,j,k,:))
        endif
      enddo
    enddo

    !if (nrank==0) write(*,*) "Outflow velocity ux nx=n min max=",real(uxmin1,4),real(uxmax1,4)

    return
  end subroutine outflow 

  !*******************************************************************************
  !
  subroutine postprocess_uniform(ux1,uy1,uz1,ep1)
  !
  !*******************************************************************************

    USE MPI
    USE decomp_2d
    USE decomp_2d_io
    USE var, only : umean,vmean,wmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
    USE var, only : uvisu
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    character(len=30) :: filename
    
    ! Write vorticity as an example of post processing
!    if ((ivisu.ne.0).and(mod(itime, ioutput).eq.0)) then
!      !x-derivatives
!      call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
!      call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
!      call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
!      !y-derivatives
!      call transpose_x_to_y(ux1,td2)
!      call transpose_x_to_y(uy1,te2)
!      call transpose_x_to_y(uz1,tf2)
!      call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
!      call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
!      call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
!      !!z-derivatives
!      call transpose_y_to_z(td2,td3)
!      call transpose_y_to_z(te2,te3)
!      call transpose_y_to_z(tf2,tf3)
!      call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
!      call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
!      call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
!      !!all back to x-pencils
!      call transpose_z_to_y(ta3,td2)
!      call transpose_z_to_y(tb3,te2)
!      call transpose_z_to_y(tc3,tf2)
!      call transpose_y_to_x(td2,tg1)
!      call transpose_y_to_x(te2,th1)
!      call transpose_y_to_x(tf2,ti1)
!      call transpose_y_to_x(ta2,td1)
!      call transpose_y_to_x(tb2,te1)
!      call transpose_y_to_x(tc2,tf1)
!      !du/dx=ta1 du/dy=td1 and du/dz=tg1
!      !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
!      !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1
!
!      ! vorticity
!      di1(:,:,:)=sqrt((tf1(:,:,:)-th1(:,:,:))**2+(tg1(:,:,:)-tc1(:,:,:))**2+&
!           (tb1(:,:,:)-td1(:,:,:))**2)
!      if (iibm==2) then
!        di1(:,:,:) = (one - ep1(:,:,:)) * di1(:,:,:)
!      endif
!      uvisu=0.
!      call fine_to_coarseV(1,di1,uvisu)
!994   format('./data/vort',I5.5)
!      write(filename, 994) itime/ioutput
!      call decomp_2d_write_one(1,uvisu,filename,2)
!    endif

    return
  end subroutine postprocess_uniform

end module uniform
