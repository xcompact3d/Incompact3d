!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module uniform

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
  PUBLIC :: init_uniform, boundary_conditions_uniform, postprocess_uniform, &
            visu_uniform, visu_uniform_init

contains

  !*******************************************************************************
  !
  subroutine init_uniform (ux1,uy1,uz1,ep1,phi1)
  !
  !*******************************************************************************

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
    USE MPI

    implicit none

    integer :: j,k,code
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx,uxmin,uxmax

    udx=one/dx; udy=one/dy; udz=one/dz; uddx=half/dx; uddy=half/dy; uddz=half/dz

    uxmax=-1609.
    uxmin=1609.
    do k=1,xsize(3)
      do j=1,xsize(2)
        if (ux(nx-1,j,k).gt.uxmax) uxmax=ux(nx-1,j,k)
        if (ux(nx-1,j,k).lt.uxmin) uxmin=ux(nx-1,j,k)
      enddo
    enddo

    call MPI_ALLREDUCE(MPI_IN_PLACE,uxmax,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(MPI_IN_PLACE,uxmin,1,real_type,MPI_MIN,MPI_COMM_WORLD,code)

    cx=0.5*(uxmax+uxmin)*gdt(itr)*udx
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
    USE decomp_2d_io
    USE var, only : uvisu
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    USE ibm_param
    
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
 
  end subroutine postprocess_uniform
  
   subroutine visu_uniform_init (visu_initialised)

    use decomp_2d_io, only : decomp_2d_register_variable
    use visu, only : io_name, output2D
    
    implicit none

    logical, intent(out) :: visu_initialised

    call decomp_2d_register_variable(io_name, "vort", 1, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "critq", 1, 0, output2D, mytype)

    visu_initialised = .true.
    
  end subroutine visu_uniform_init

    subroutine visu_uniform(ux1, uy1, uz1, pp3, phi1, ep1, num)

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
    !VORTICITY FIELD
    di1 = zero
    di1(:,:,:)=sqrt(  (tf1(:,:,:)-th1(:,:,:))**2 &
                    + (tg1(:,:,:)-tc1(:,:,:))**2 &
                    + (tb1(:,:,:)-td1(:,:,:))**2)
    call write_field(di1, ".", "vort", num, flush = .true.) ! Reusing temporary array, force flush

    !Q=-0.5*(ta1**2+te1**2+ti1**2)-td1*tb1-tg1*tc1-th1*tf1
    di1 = zero
    di1(:,:,: ) = - half*(ta1(:,:,:)**2+te1(:,:,:)**2+ti1(:,:,:)**2) &
                  - td1(:,:,:)*tb1(:,:,:) &
                  - tg1(:,:,:)*tc1(:,:,:) &
                  - th1(:,:,:)*tf1(:,:,:)
    call write_field(di1, ".", "critq", num, flush = .true.) ! Reusing temporary array, force flush

  end subroutine visu_uniform

end module uniform
