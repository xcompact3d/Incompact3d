!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module abl

use decomp_2d_constants
use decomp_2d_mpi
use decomp_2d

contains

  !*******************************************************************************
  !
  subroutine init_abl(ux1,uy1,uz1,ep1,phi1)
  !
  !*******************************************************************************

    use decomp_2d_io
    use variables
    use param
    use MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: y, phinoise
    integer :: k,j,i,ierror,ii,code
    integer, dimension (:), allocatable :: seed

    ux1=zero
    uy1=zero
    uz1=zero
    phi1=zero
    bxx1=zero;bxy1=zero;bxz1=zero
    byx1=zero;byy1=zero;byz1=zero
    bzx1=zero;bzy1=zero;bzz1=zero

    ! ABL not yet set up for iLES, stretched grids, and non-constant explicit models. 
    if (ilesmod == 0.or.istret /= 0.or.jles.gt.1) then
        write(*,*)  'Simulation stopped: run with different options'
       call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
    endif

    ! Generation of a random noise
    if (iin /= 0) then
      call system_clock(count=code)
      call random_seed(size = ii)
      call random_seed(put = code+63946*(nrank+1)*(/ (i - 1, i = 1, ii) /)) !

      call random_number(ux1)
      call random_number(uy1)
      call random_number(uz1)

      do k=1,xsize(3)
      do j=1,xsize(2)
      do i=1,xsize(1)
        ux1(i,j,k)=init_noise*(ux1(i,j,k)*two-one)
        uy1(i,j,k)=init_noise*(uy1(i,j,k)*two-one)
        uz1(i,j,k)=init_noise*(uz1(i,j,k)*two-one)
      enddo
      enddo
      enddo
    endif

    ! Initialize with log-law or geostrophic wind
    do k=1,xsize(3)
    do j=1,xsize(2)
       if (istret == 0) y=real(j+xstart(2)-1-1,mytype)*dy
       if (istret /= 0) y=yp(j)
       if (iPressureGradient.eq.1.or.imassconserve.eq.1) then
           bxx1(j,k)=ustar/k_roughness*log((y+z_zero)/z_zero)
       else
           bxx1(j,k)=UG(1)
       endif
       bxy1(j,k)=UG(2)
       bxz1(j,k)=UG(3)
    enddo
    enddo

    ! Add boundary to noisy velocities
    do k=1,xsize(3)
    do j=1,xsize(2)
    do i=1,xsize(1)
       ux1(i,j,k)=bxx1(j,k)*(one+ux1(i,j,k))
       uy1(i,j,k)=uy1(i,j,k)
       uz1(i,j,k)=uz1(i,j,k)
    enddo
    enddo
    enddo

    ! Initialize temperature profiles
    if (iscalar == 1) then
      do j=1,xsize(2)
        if (istret == 0) y=real(j + xstart(2)-1-1,mytype)*dy
        if (istret /= 0) y=yp(j+xstart(2)-1)
        if (ibuoyancy == 1) then 
          Tstat(j,1)=T_wall + (T_top-T_wall)*y/yly
        else 
          Tstat(j,1)=zero
        endif
        ! Initialize GABLS-1 case
        if (istrat==0) then
          if (y>onehundred) then
            phi1(:,j,:,1)=T_wall-Tstat(j,1) + (y-onehundred)*one/onehundred
          else
            phi1(:,j,:,1)=T_wall-Tstat(j,1)
          endif
        ! Initialize case from Gadde et al. (2020)
        else if (istrat==1) then
          if (y>1062._mytype) then
            phi1(:,j,:,1)=T_wall-Tstat(j,1) + eight + (y-1062._mytype)*three/onethousand
          else if (y>937._mytype) then
            phi1(:,j,:,1)=T_wall-Tstat(j,1) + (y-937._mytype)*eight/125._mytype
          else 
            phi1(:,j,:,1)=T_wall-Tstat(j,1)
          endif
        endif
      enddo

      ! Add random noise
      do j=1,xsize(2)
        if (istret==0) y=real(j + xstart(2)-1-1,mytype)*dy
        if (istret/=0) y=yp(j+xstart(2)-1)
        !if (y.lt.50) then 
        !  do k=1,xsize(3)
        !  do i=1,xsize(1)
        !    call random_number(phinoise)
        !    phinoise=0.1*(phinoise*2.-1.)
        !    phi1(i,j,k,1)=phi1(i,j,k,1)+phinoise
        !  enddo
        !  enddo
        !endif
      enddo
    endif

    return
  end subroutine init_abl

  !*******************************************************************************
  !
  subroutine boundary_conditions_abl(ux,uy,uz,phi)
  !
  !*******************************************************************************

    USE param
    USE variables

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: gx
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    if (imassconserve==1) then
      call transpose_x_to_y(ux,gx)
      call forceabl(gx)
      call transpose_y_to_x(gx,ux)
    endif

    if ((iconcprec.eq.1).or.(ishiftedper.eq.1)) then
      call fringe_region(ux,uy,uz)
    endif

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
  end subroutine boundary_conditions_abl

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
      if (mod(itime,ilist)==0.and.nrank==0) print *,'Reading inflow from a file, time step: ', itime_input
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

    return
  end subroutine outflow 

  !*******************************************************************************
  !
  subroutine momentum_forcing_abl(dux1,duy1,duz1,ux1,uy1,uz1,phi1)
  !
  !*******************************************************************************

    USE param
    USE variables

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3), ntime) :: dux1, duy1, duz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3), numscalar) :: phi1

    integer :: i
   
    ! BL Forcing (Pressure gradient or geostrophic wind)
    if (iPressureGradient==1) then
       dux1(:,:,:,1)=dux1(:,:,:,1)+ustar**2./dBL
       if (iconcprec.eq.1) then
          do i=1,xsize(1)
             if (real(i-1,mytype)*dx >= pdl) then
                dux1(i,:,:,1)=dux1(i,:,:,1)-ustar**2./dBL
             endif
          enddo
       endif
    else if (iCoriolis==1 .and. iPressureGradient==0) then
       dux1(:,:,:,1)=dux1(:,:,:,1)+CoriolisFreq*(-UG(3))
       duz1(:,:,:,1)=duz1(:,:,:,1)-CoriolisFreq*(-UG(1))
    endif

    ! Coriolis terms
    if (iCoriolis==1) then
       dux1(:,:,:,1)=dux1(:,:,:,1)+CoriolisFreq*uz1(:,:,:)
       duz1(:,:,:,1)=duz1(:,:,:,1)-CoriolisFreq*ux1(:,:,:)
    endif

    ! Damping zone
    if (idamping==1) then
       call damping_zone(dux1,duy1,duz1,ux1,uy1,uz1)
    endif

    ! Buoyancy terms
    if (iscalar==1.and.ibuoyancy==1) then
       duy1(:,:,:,1)=duy1(:,:,:,1)+gravv*phi1(:,:,:,1)/Tref
    endif

    return
  end subroutine momentum_forcing_abl

  !*******************************************************************************
  !
  subroutine scalar_forcing_abl(uy1,dphi1,phi1)
  !
  !*******************************************************************************

    USE param
    USE variables

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3), ntime) :: dphi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: uy1, phi1

    ! Damping zone
    if (idamping==1) then
       call damping_zone_scalar(dphi1,phi1)
    endif

    ! Terms from decomposition
    if (ibuoyancy==1) then
       dphi1(:,:,:,1) = dphi1(:,:,:,1) + (T_wall-T_top)*uy1(:,:,:)/yly
    endif

    return
  end subroutine scalar_forcing_abl

  !*******************************************************************************
  !
  subroutine wall_sgs_slip(ux,uy,uz,phi,nut1,wallfluxx,wallfluxy,wallfluxz)
  !
  ! Outputs fluxes, only compatible with iconserv=0
  !
  !*******************************************************************************

    use MPI
    use param
    use variables
    use var, only: uxf1, uzf1, phif1, uxf3, uzf3, phif3
    use var, only: di1, di3
    use var, only: sxy1, syz1, heatflux, ta2, tb2, ta3, tb3
    use ibm_param, only : ubcx, ubcz
   
    implicit none
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux,uy,uz, nut1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar),intent(in) :: phi
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(out) :: wallfluxx,wallfluxy,wallfluxz
    real(mytype),dimension(xsize(1),xsize(3)) :: tauwallxy, tauwallzy
    real(mytype),dimension(xsize(1),xsize(3)) :: Obukhov, zeta
    integer :: i,j,k,ii,code
    integer :: nxc, nyc, nzc, xsize1, xsize2, xsize3
    real(mytype) :: delta
    real(mytype) :: ux_HAve_local, uz_HAve_local, Phi_HAve_local
    real(mytype) :: ux_HAve, uz_HAve,S_HAve,Phi_HAve,ux12,uz12,S12,Phi12,Tstat12
    real(mytype) :: PsiM_HAve_local, PsiM_HAve, PsiH_HAve_local, PsiH_HAve
    real(mytype) :: L_HAve_local, L_HAve, Q_HAve_local, Q_HAve, zL, zeta_HAve
    real(mytype) :: Lold, OL_diff

    ! Filter the velocity with twice the grid scale according to Bou-Zeid et al. (2005)

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

    call filter(zero)
    call filx(uxf1,ux,di1,fisx,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0,ubcx)
    call filx(uzf1,uz,di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1,ubcz)
    call transpose_x_to_y(uxf1,ta2)
    call transpose_x_to_y(uzf1,tb2)
    call transpose_y_to_z(ta2,ta3)
    call transpose_y_to_z(tb2,tb3)
    call filz(uxf3,ta3,di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
    call filz(uzf3,tb3,di3,fisz,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),0,ubcz)
    call transpose_z_to_y(uxf3,ta2)
    call transpose_z_to_y(uzf3,tb2)
    call transpose_y_to_x(ta2,uxf1)
    call transpose_y_to_x(tb2,uzf1)

    if (iscalar==1) then
      call filx(phif1,phi(:,:,:,1),di1,fisx,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0,zero)
      call transpose_x_to_y(phif1,ta2)
      call transpose_y_to_z(ta2,ta3)
      call filz(phif3,ta3,di3,fisz,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),0,zero)
      call transpose_z_to_y(phif3,ta2)
      call transpose_y_to_x(ta2,phif1)
    endif

    ! Reset average values
    ux_HAve_local  = zero
    uz_HAve_local  = zero
    Phi_HAve_local = zero

    ! dy to y=1/2
    if (istret/=0) delta=half*(yp(2)-yp(1))
    if (istret==0) delta=half*dy

    ! Find horizontally averaged velocities at j=1.5
    if (xstart(2)==1) then
      do k=1,xsize(3)
        do i=1,xsize(1)
           ux_HAve_local=ux_HAve_local+half*(uxf1(i,1,k)+uxf1(i,2,k))
           uz_HAve_local=uz_HAve_local+half*(uzf1(i,1,k)+uzf1(i,2,k))
           if (iscalar==1) Phi_HAve_local=Phi_HAve_local+half*(phif1(i,1,k)+phif1(i,2,k))
        enddo
      enddo
      ux_HAve_local=ux_HAve_local
      uz_HAve_local=uz_HAve_local
      Phi_HAve_local=Phi_HAve_local
    endif

    call MPI_ALLREDUCE(ux_HAve_local,ux_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(uz_HAve_local,uz_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    if (iscalar==1) call MPI_ALLREDUCE(Phi_HAve_local,Phi_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

    ux_HAve=ux_HAve/(nxc*nzc)
    uz_HAve=uz_HAve/(nxc*nzc)
    S_HAve=sqrt(ux_HAve**2.+uz_HAve**2.)
    if (iscalar==1) then 
      Phi_HAve=Phi_HAve/(nxc*nzc)
      if (ibuoyancy==1) then 
        Tstat12 =T_wall + (T_top-T_wall)*delta/yly
      else 
        Tstat12 =zero
      endif
      Phi_HAve=Phi_HAve + Tstat12
    endif

    ! Reset wall flux values
    wallfluxx=zero
    wallfluxy=zero
    wallfluxz=zero

    ! Initialize stratification variables
    if (iscalar==1.and.ibuoyancy == 1.and.xstart(2)==1) then 
      PsiM_HAve= zero
      PsiH_HAve= zero
      ii       = 0
      OL_diff  = one
      Lold     = one
      do while (OL_diff > 1.0e-14_mytype)
        if (itherm==0) then
          Q_HAve = TempFlux
        else if (itherm==1) then
          Q_HAve =-k_roughness**two*S_HAve*(Phi_HAve-(T_wall+TempRate*t))/((log(delta/z_zero)-PsiM_HAve)*(log(delta/z_zero)-PsiH_HAve))
        endif
        L_HAve=-(k_roughness*S_HAve/(log(delta/z_zero)-PsiM_HAve))**three*Phi_HAve/(k_roughness*gravv*Q_HAve)
        if (istrat==0) then
          PsiM_HAve=-4.8_mytype*delta/L_HAve
          PsiH_HAve=-7.8_mytype*delta/L_HAve
        else if (istrat==1) then
          zeta_HAve=(one-sixteen*delta/L_HAve)**zptwofive
          PsiM_HAve=two*log(half*(one+zeta_HAve))+log(zpfive*(one+zeta_HAve**two))-two*atan(zeta_HAve)+pi/two
          PsiH_HAve=two*log(half*(one+zeta_HAve**two))
        endif
        ii      = ii + 1
        OL_diff = abs((L_HAve - Lold)/Lold)
        Lold    = L_HAve
        if (ii==50) exit
      enddo
      heatflux=Q_Have
      Obukhov=L_HAve
      PsiM=PsiM_HAve
      PsiH=PsiH_HAve
      if (istrat==1) zeta=zeta_HAve
    else   
      heatflux =zero
      Obukhov  =zero
      PsiM     =zero
      PsiH     =zero
      PsiM_HAve=zero
      PsiH_HAve=zero
    endif

    ! Apply BCs locally
    if (xstart(2)==1) then
      do k=1,xsize(3)
      do i=1,xsize(1)
         ! Horizontally-averaged formulation
         if(iwallmodel==1) then
           tauwallxy(i,k)=-(k_roughness/(log(delta/z_zero)-PsiM_HAve))**two*ux_HAve*S_HAve
           tauwallzy(i,k)=-(k_roughness/(log(delta/z_zero)-PsiM_HAve))**two*uz_HAve*S_HAve
         ! Local formulation
         else
           ux12=half*(uxf1(i,1,k)+uxf1(i,2,k))
           uz12=half*(uzf1(i,1,k)+uzf1(i,2,k))
           S12=sqrt(ux12**2.+uz12**2.)
           if (iscalar==1) then
             Phi12= half*(phif1(i,1,k)+ phif1(i,2,k)) + Tstat12
             do ii=1,10
                if (itherm==1) heatflux(i,k)=-k_roughness**two*S12*(Phi12-(T_wall+TempRate*t))/((log(delta/z_zero)-PsiM(i,k))*(log(delta/z_zero)-PsiH(i,k)))
                Obukhov(i,k)=-(k_roughness*S12/(log(delta/z_zero)-PsiM(i,k)))**three*Phi12/(k_roughness*gravv*heatflux(i,k))
                if (istrat==0) then
                  PsiM(i,k)=-4.8_mytype*delta/Obukhov(i,k)
                  PsiH(i,k)=-7.8_mytype*delta/Obukhov(i,k)
                else if (istrat==1) then
                  zeta(i,k)=(one-sixteen*delta/Obukhov(i,k))**zptwofive
                  PsiM(i,k)=two*log(half*(one+zeta(i,k)))+log(zpfive*(one+zeta(i,k)**2.))-two*atan(zeta(i,k))+pi/two
                  PsiH(i,k)=two*log(half*(one+zeta(i,k)**two))
                endif
             enddo
           endif
           tauwallxy(i,k)=-(k_roughness/(log(delta/z_zero)-PsiM(i,k)))**two*ux12*S12
           tauwallzy(i,k)=-(k_roughness/(log(delta/z_zero)-PsiM(i,k)))**two*uz12*S12
         endif
         ! Apply second-order upwind scheme for the near wall
         ! Below should change for non-uniform grids, same for wall_sgs_slip_scalar
         wallfluxx(i,1,k) = -(-half*(-two*nut1(i,3,k)*sxy1(i,3,k))+&
                            two*(-two*nut1(i,2,k)*sxy1(i,2,k))-three/two*tauwallxy(i,k))/(two*delta)
         wallfluxy(i,1,k) = zero
         wallfluxz(i,1,k) = -(-half*(-two*nut1(i,3,k)*syz1(i,3,k))+&
                            two*(-two*nut1(i,2,k)*syz1(i,2,k))-three/two*tauwallzy(i,k))/(two*delta)
      enddo
      enddo
    endif

    ! Reset average values
    PsiM_HAve_local=zero
    PsiH_HAve_local=zero
    L_HAve_local   =zero
    Q_HAve_local   =zero

    ! Find horizontally averaged values
    if (iscalar==1) then
       do k=1,xsize(3)
          do i=1,xsize(1)
            PsiM_HAve_local=PsiM_HAve_local+PsiM(i,k)
            PsiH_HAve_local=PsiH_HAve_local+PsiH(i,k)
            L_HAve_local=L_HAve_local+Obukhov(i,k)
            Q_HAve_local=Q_HAve_local+heatflux(i,k)
          enddo
       enddo
       PsiM_HAve_local=PsiM_HAve_local
       PsiH_HAve_local=PsiH_HAve_local
       L_HAve_local=L_HAve_local
       Q_HAve_local=Q_HAve_local

       call MPI_ALLREDUCE(PsiM_HAve_local,PsiM_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       call MPI_ALLREDUCE(PsiH_HAve_local,PsiH_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       call MPI_ALLREDUCE(L_HAve_local,L_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       call MPI_ALLREDUCE(Q_HAve_local,Q_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

       PsiM_HAve=PsiM_HAve/(nxc*nzc)
       PsiH_HAve=PsiH_HAve/(nxc*nzc)
       L_HAve=L_HAve/(nxc*nzc)
       Q_HAve=Q_HAve/(nxc*nzc)
    endif
    
    ! Compute friction velocity u_shear and boundary layer height
    u_shear=k_roughness*S_HAve/(log(delta/z_zero)-PsiM_HAve)
    if (iheight==1) call boundary_height(ux,uy,uz,dBL)
    if (iscalar==1) zL=dBL/L_HAve

    if (mod(itime,ilist)==0.and.nrank==0) then
         write(*,*)  ' '
         write(*,*)  ' ABL:'
         write(*,*)  ' Horizontally-averaged velocity at y=1/2: ', ux_HAve,uz_Have
         write(*,*)  ' BL height: ', dBL
         write(*,*)  ' Friction velocity: ', u_shear
    
        if (iscalar==1) then
           write(*,*)  ' Temperature: ', Phi_HAve
           write(*,*)  ' PsiM: ', PsiM_HAve
           write(*,*)  ' PsiH: ', PsiH_HAve
           write(*,*)  ' Obukhov L: ', L_HAve
           write(*,*)  ' Heatflux: ', Q_HAve
           write(*,*)  ' z/L: ', zL
        endif

         write(*,*)  'Maximum wall shear stress for x and z', maxval(tauwallxy), maxval(tauwallzy)
         write(*,*)  'Minimum wall shear stress for x and z', minval(tauwallxy), minval(tauwallzy)
         write(*,*)  'Max flux x and z ', maxval(wallfluxx), maxval(wallfluxz)
         write(*,*)  'Min flux x and z ', minval(wallfluxx), minval(wallfluxz)
    endif

    return
  end subroutine wall_sgs_slip

  !*******************************************************************************
  !
  subroutine wall_sgs_slip_scalar(sgsphi1,nut1,dphidy1)
  !
  !*******************************************************************************

    use param
    use var, only: heatflux
    use variables

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: nut1, dphidy1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(inout) :: sgsphi1

    real(mytype) :: delta, Pr
    integer ::  i,k

    Pr=Sc(1)

    if (xstart(2)==1) then
       if (istret/=0) delta=(yp(2)-yp(1))/two
       if (istret==0) delta=dy/two
       do k=1,xsize(3)
          do i=1,xsize(1)
             sgsphi1(i,1,k) =-(-half*(-nut1(i,3,k)*dphidy1(i,3,k))/Pr+&
                             two*(-nut1(i,2,k)*dphidy1(i,2,k))/Pr-three/two*heatflux(i,k))/(two*delta)
          enddo
       enddo
    endif

  end subroutine wall_sgs_slip_scalar

  !*******************************************************************************
  !
  subroutine wall_sgs_noslip(ux1,uy1,uz1,nut1,wallsgsx1,wallsgsy1,wallsgsz1)
  !
  ! Outputs stresses if iconserv=1 and fluxes if iconserv=0 (wallsgsx,wallsgsy,wallsgsz)
  !
  !*******************************************************************************

    use MPI
    use param
    use variables
    use var, only: di1, di2, di3
    use var, only: sxy1, syz1, tb1, ta2, tb2
    use ibm_param, only : ubcx, ubcy, ubcz
   
    implicit none
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1,nut1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(out) :: wallsgsx1,wallsgsy1,wallsgsz1
    real(mytype),dimension(ysize(1),ysize(3)) :: tauwallxy2, tauwallzy2
    integer :: i,j,k,code,j0
    integer :: nxc, nyc, nzc, xsize1, xsize2, xsize3
    real(mytype) :: delta
    real(mytype) :: ux_HAve_local, uz_HAve_local
    real(mytype) :: ux_HAve, uz_HAve, S_HAve, ux_delta, uz_delta, S_delta
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: txy1,tyz1,dtwxydx
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: txy2,tyz2,wallsgsx2,wallsgsz2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: tyz3,dtwyzdz
    
    real(mytype)               :: y_sampling
    
    ! Reset wall flux/stresses values
    wallsgsx1 = zero
    wallsgsy1 = zero
    wallsgsz1 = zero
    
    ! Set sampling distance for the wall model
    delta=dsampling*dy
    
    if(iconserv==0) then
      ! Construct Smag SGS stress tensor 
      txy1 = -2.0*nut1*sxy1
      tyz1 = -2.0*nut1*syz1
      call transpose_x_to_y(txy1,txy2)
      call transpose_x_to_y(tyz1,tyz2)
    endif

    ! Work on Y-pencil
    call transpose_x_to_y(ux1,ta2)
    call transpose_x_to_y(uz1,tb2)
 
    ! Apply BCs locally
    do k=1,ysize(3)
    do i=1,ysize(1)
      !sampling at dsampling*dy from wall
      j0=floor(delta/dy)
      y_sampling=delta-real(j0,mytype)*dy
      ux_delta=(1-y_sampling/dy)*ta2(i,j0+1,k)+(y_sampling/dy)*ta2(i,j0+2,k)
      uz_delta=(1-y_sampling/dy)*tb2(i,j0+1,k)+(y_sampling/dy)*tb2(i,j0+2,k)
      S_delta=sqrt(ux_delta**2.+uz_delta**2.)
      tauwallxy2(i,k)=-(k_roughness/(log(delta/z_zero)))**two*ux_delta*S_delta
      tauwallzy2(i,k)=-(k_roughness/(log(delta/z_zero)))**two*uz_delta*S_delta
      txy2(i,2,k) = tauwallxy2(i,k)
      tyz2(i,2,k) = tauwallzy2(i,k)
    enddo
    enddo

    if (iconserv==0) then
      ! Derivative of wallmodel-corrected SGS stress tensor
      call dery_22(wallsgsx2,txy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),1,ubcy)
      call dery_22(wallsgsz2,tyz2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),1,ubcy)
      call transpose_y_to_x(wallsgsx2,wallsgsx1)
      call transpose_y_to_x(wallsgsz2,wallsgsz1)
      call derx(dtwxydx,txy1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)
      call transpose_y_to_z(tyz2,tyz3)
      call derz(dtwyzdz,tyz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcz)
      call transpose_z_to_y(dtwyzdz,tb2)
      call transpose_y_to_x(tb2,tb1)
      wallsgsy1 = dtwxydx + tb1
    elseif (iconserv==1) then 
      call transpose_y_to_x(txy2,wallsgsx1)
      call transpose_y_to_x(tyz2,wallsgsz1)
    endif
    

    ! Print information at y=5*dy
    ux_HAve_local  =zero
    uz_HAve_local  =zero
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
    do k=1,ysize(3)
    do i=1,ysize(1)
       ux_HAve_local=ux_HAve_local+ta2(i,6,k)
       uz_HAve_local=uz_HAve_local+tb2(i,6,k)
    enddo
    enddo
    ux_HAve_local=ux_HAve_local
    uz_HAve_local=uz_HAve_local
    call MPI_ALLREDUCE(ux_HAve_local,ux_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(uz_HAve_local,uz_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    ux_HAve=ux_HAve/(nxc*nzc)
    uz_HAve=uz_HAve/(nxc*nzc)
    S_HAve=sqrt(ux_HAve**2.+uz_HAve**2.)
    u_shear=k_roughness*S_HAve/log(5*dy/z_zero)
    
    if (mod(itime,ilist)==0.and.nrank==0) then
       ! Write u_shear in file
       write(42,'(20e20.12)') t,u_shear
       flush(42)
       ! Print in terminal
       write(*,*)  ' ABL:'
       write(*,*)  ' Horizontally-averaged velocity at 5*dy: ', ux_HAve,uz_HAve
       write(*,*)  ' Friction velocity at 5*dy: ', u_shear
    endif

    return
  end subroutine wall_sgs_noslip

  !*******************************************************************************
  !
  subroutine forceabl (ux) ! Routine to force constant flow rate
  !
  !*******************************************************************************

    use decomp_2d_poisson
    use param
    use var
    use MPI

    implicit none

    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux
    integer :: j,i,k,code
    real(mytype) :: can,ut3,ut,ut4,xloc

    ut3=zero
    do k=1,ysize(3)
      do i=1,ysize(1)
        xloc=(i+ystart(1)-1-1)*dx
        if (iconcprec.eq.1.and.xloc>=pdl) then
          continue  
        else
          ut=zero
          do j=1,ny-1
            if (istret/=0) ut=ut+(yp(j+1)-yp(j))*(ux(i,j+1,k)-half*(ux(i,j+1,k)-ux(i,j,k)))
            if (istret==0) ut=ut+(yly/real(ny-1,mytype))*(ux(i,j+1,k)-half*(ux(i,j+1,k)-ux(i,j,k)))
          enddo
          ut3=ut3+ut
        endif
      enddo
    enddo
    ut3=ut3/ysize(1)/ysize(3)

    call MPI_ALLREDUCE(ut3,ut4,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    ut4=ut4/nproc
    if (iconcprec.eq.1) ut4=ut4*(xlx/pdl)

    ! Flow rate for a logarithmic profile
    !can=-(ustar/k_roughness*yly*(log(yly/z_zero)-1.)-ut4)
    can=-(ustar/k_roughness*(yly*log(dBL/z_zero)-dBL)-ut4)

    if (nrank==0.and.mod(itime,ilist)==0)  write(*,*) '# Rank ',nrank,'correction to ensure constant flow rate',ut4,can

    do k=1,ysize(3)
      do i=1,ysize(1)
        xloc=real(i+ystart(1)-1-1,mytype)*dx
        if (iconcprec.eq.1.and.xloc>=pdl) then
          continue
        else
          do j=1,ny
            ux(i,j,k)=ux(i,j,k)-can/yly
          enddo
        endif
      enddo
    enddo

    return
  end subroutine forceabl
     
  !*******************************************************************************
  !
  subroutine fringe_region (ux,uy,uz)
  !
  !*******************************************************************************

    USE param
    USE var

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux_s,uy_s,uz_s
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ux3_s,uy3_s,uz3_s
    integer :: i, j, k, nshift, npe, nfe
    real(mytype) :: dshift, x, edl, frl, fre, frs, frd, lambda

    if (ishiftedper==1) then
      ! Shifting distance
      dshift = zlz/8 ! could be yly/4 etc..
      nshift = int(dshift/dz)

      ! Compute spanwise-shifted velocity field
      call transpose_x_to_y(ux,td2)
      call transpose_x_to_y(uy,te2)
      call transpose_x_to_y(uz,tf2)
      call transpose_y_to_z(td2,td3)
      call transpose_y_to_z(te2,te3)
      call transpose_y_to_z(tf2,tf3)

      do k=1,nz-nshift
        ux3_s(:,:,k+nshift) = td3(:,:,k)
        uy3_s(:,:,k+nshift) = te3(:,:,k)
        uz3_s(:,:,k+nshift) = tf3(:,:,k)
      enddo
      do k=1,nshift
        ux3_s(:,:,k) = td3(:,:,nz-nshift+k)
        uy3_s(:,:,k) = te3(:,:,nz-nshift+k)
        uz3_s(:,:,k) = tf3(:,:,nz-nshift+k)
      enddo 

      call transpose_z_to_y(ux3_s,td2)
      call transpose_z_to_y(uy3_s,te2)
      call transpose_z_to_y(uz3_s,tf2)
      call transpose_y_to_x(td2,ux_s)
      call transpose_y_to_x(te2,uy_s)
      call transpose_y_to_x(tf2,uz_s)
    else
      ux_s=ux
      uy_s=uy
      uz_s=uz
    endif

    ! Fringe region(s) parameters
    edl = twothird 
    frl = 1._mytype/6._mytype
    if (ishiftedper==1.and.iconcprec==0) then
      fre = xlx*edl
      npe = nx
    elseif (ishiftedper==1.and.iconcprec==1) then
      fre = pdl*edl
      npe = int(nx*pdl/xlx)
    elseif (ishiftedper==0.and.iconcprec==1) then
      fre = pdl
    endif
    frs = fre-fre*frl
    frd = fre-frs
    nfe = int(nx*fre/xlx)

    ! Apply fringe region(s)
    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,nfe
          x=real(i-1,mytype)*dx
          if (x<frs) then
            lambda=zero
          elseif ( (x>=frs) .and. (x<(fre-frd/four)) ) then
            lambda=half*(1.-cos(four*pi/three*(x-frs)/frd))
          elseif ( (x>=(fre-frd/four)) .and. (x<fre) ) then
            lambda=one
          else 
            lambda=zero
          endif
          if (ishiftedper==1) then
            ux(i+npe-nfe,j,k)=lambda*ux_s(i,j,k)+(one-lambda)*ux(i+npe-nfe,j,k)
            uy(i+npe-nfe,j,k)=lambda*uy_s(i,j,k)+(one-lambda)*uy(i+npe-nfe,j,k)
            uz(i+npe-nfe,j,k)=lambda*uz_s(i,j,k)+(one-lambda)*uz(i+npe-nfe,j,k)
          endif
          if (iconcprec==1) then
            ux(i+nx-nfe,j,k)=lambda*ux_s(i,j,k)+(one-lambda)*ux(i+nx-nfe,j,k)
            uy(i+nx-nfe,j,k)=lambda*uy_s(i,j,k)+(one-lambda)*uy(i+nx-nfe,j,k)
            uz(i+nx-nfe,j,k)=lambda*uz_s(i,j,k)+(one-lambda)*uz(i+nx-nfe,j,k)
          endif
        enddo
      enddo
    enddo

    return
  end subroutine fringe_region

  !*******************************************************************************
  !
  subroutine damping_zone (dux1,duy1,duz1,ux1,uy1,uz1) ! Damping zone for ABL
  !
  !*******************************************************************************

    use param
    use var, only: yp

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime),intent(inout) :: dux1, duy1, duz1
    integer :: i,j,k
    real(mytype) :: y, lambda, xloc
    real(mytype) :: damp_lo, coeff, wvar, dheight
   
    if (ibuoyancy.eq.1) then
      damp_lo = 300._mytype
      coeff   = 0.0016_mytype !0.5*ustar/dBL
    else
      dheight = zpone*dBL
      wvar    = fifteen !1./(1./k_roughness*(yly*log(yly/dBL)-yly+dBL)/(yly-dBL))
      coeff   = wvar*ustar/dBL
    endif

    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          if (istret==0) y=real(j+xstart(2)-1-1,mytype)*dy
          if (istret/=0) y=yp(j+xstart(2)-1)
          ! Damping for non-neutral ABL
          if (ibuoyancy==1) then
            if (y>=damp_lo) then
              lambda=sin(half*pi*(y-damp_lo)/(yly-damp_lo))**two
            else
              lambda=zero
            endif
            dux1(i,j,k,1)=dux1(i,j,k,1)-coeff*lambda*(ux1(i,j,k)-UG(1))
            duy1(i,j,k,1)=duy1(i,j,k,1)-coeff*lambda*(uy1(i,j,k)-UG(2))
            duz1(i,j,k,1)=duz1(i,j,k,1)-coeff*lambda*(uz1(i,j,k)-UG(3))
          ! Damping for neutral ABL
          else
            if (y>=(dBL+half*dheight)) then
              lambda=one
            elseif (y>=(dBL-half*dheight).and.y<(dBL+half*dheight)) then
              lambda=half*(one-cos(pi*(y-(dBL-half*dheight))/dheight))
            else
             lambda=zero
            endif
            xloc=real(i-1,mytype)*dx
            if (iconcprec.eq.1.and.xloc.ge.pdl) lambda=0.
            dux1(i,j,k,1)=dux1(i,j,k,1)-coeff*lambda*(ux1(i,j,k)-ustar/k_roughness*log(dBL/z_zero))
            duy1(i,j,k,1)=duy1(i,j,k,1)-coeff*lambda*(uy1(i,j,k)-UG(2))
            duz1(i,j,k,1)=duz1(i,j,k,1)-coeff*lambda*(uz1(i,j,k)-UG(3))
          endif
        enddo
      enddo
    enddo

  end subroutine damping_zone
  
  !*******************************************************************************
  !
  subroutine damping_zone_scalar (dphi1,phi1)
  !
  !*******************************************************************************

    use param
    use var, only: yp

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime),intent(inout) :: dphi1
    integer :: i,j,k
    real(mytype) :: y, lambda
    real(mytype) :: damp_lo, coeff

    damp_lo = 300._mytype
    coeff   = zero!0.0016_mytype

    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          if (istret==0) y=real(j+xstart(2)-1-1,mytype)*dy
          if (istret/=0) y=yp(j+xstart(2)-1)
          if (y>=damp_lo) then
            lambda=sin(half*pi*(y-damp_lo)/(yly-damp_lo))**two
          else
            lambda=zero
          endif
          dphi1(i,j,k,1)=dphi1(i,j,k,1)-coeff*lambda*(phi1(i,j,k)-zero) !SR why -0
        enddo
      enddo
    enddo

  end subroutine damping_zone_scalar

  !*******************************************************************************
  !
  subroutine boundary_height(ux,uy,uz,hBL) ! routine to find the height of the boundary layer
  !
  !*******************************************************************************

    use MPI
    use param
    use variables
    use var, only : sxz1, nut1, sxy1, syz1, phi1

    implicit none
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux, uy, uz
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: uxy, uyz, temp
    real(mytype), dimension(ny) :: u_HAve_local, u_HAve, v_HAve_local, v_HAve, w_HAve_local, w_HAve
    real(mytype), dimension(ny) :: uxy_HAve_local, uxy_HAve, uyz_HAve_local, uyz_HAve, momfl
    real(mytype), dimension(ny) :: txy_HAve_local, txy_HAve, tyz_HAve_local, tyz_HAve
    real(mytype) :: h, hBL
    integer :: i,j,k,code

    ! Convective BL: buoyancy flux (not implemented)
    !if (iscalar==1.and.istrat==1) then
    !  temp = phi1(:,:,:,1)
    !  do j=1,xsize(2)
    !  do k=1,xsize(3)
    !  do i=1,xsize(1)
    !    !temp(i,j,k) = phi1(i,j,k,1) + Tstat(j,1)
    !  enddo
    !  enddo
    !  enddo
    !  !uyt = uy*temp
    !endif    

    ! Compute BL height (see Kosovic & Curry, 2000)
    u_HAve_local  =zero
    v_HAve_local  =zero
    w_HAve_local  =zero
    uxy_HAve_local=zero
    uyz_HAve_local=zero
    txy_HAve_local=zero
    tyz_HAve_local=zero
    do j=1,xsize(2)
      do k=1,xsize(3)
        do i=1,xsize(1)
          u_HAve_local(j+xstart(2)-1) = u_HAve_local(j+xstart(2)-1) + ux(i,j,k)
          v_HAve_local(j+xstart(2)-1) = v_HAve_local(j+xstart(2)-1) + uy(i,j,k)
          w_HAve_local(j+xstart(2)-1) = w_HAve_local(j+xstart(2)-1) + uz(i,j,k)
          uxy_HAve_local(j+xstart(2)-1) = uxy_HAve_local(j+xstart(2)-1) + ux(i,j,k)*uy(i,j,k)
          uyz_HAve_local(j+xstart(2)-1) = uyz_HAve_local(j+xstart(2)-1) + uy(i,j,k)*uz(i,j,k)
          txy_HAve_local(j+xstart(2)-1) = txy_HAve_local(j+xstart(2)-1) - two*nut1(i,j,k)*sxy1(i,j,k)
          tyz_HAve_local(j+xstart(2)-1) = tyz_HAve_local(j+xstart(2)-1) - two*nut1(i,j,k)*syz1(i,j,k)
        enddo
      enddo
    enddo

    call MPI_ALLREDUCE(u_HAve_local,u_HAve,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(v_HAve_local,v_HAve,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(w_HAve_local,w_HAve,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(uxy_HAve_local,uxy_HAve,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(uyz_HAve_local,uyz_HAve,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(txy_HAve_local,txy_HAve,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(tyz_HAve_local,tyz_HAve,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    u_HAve  =  u_HAve/real(nx*nz,mytype)
    v_HAve  =  v_HAve/real(nx*nz,mytype)
    w_HAve  =  w_HAve/real(nx*nz,mytype)
    uxy_HAve=uxy_HAve/real(nx*nz,mytype)
    uyz_HAve=uyz_HAve/real(nx*nz,mytype)
    txy_HAve=txy_HAve/real(nx*nz,mytype)
    tyz_HAve=tyz_HAve/real(nx*nz,mytype)

    momfl = sqrt((uxy_HAve+txy_HAve-u_HAVE*v_HAve)**two + (uyz_HAve+tyz_HAve-v_HAVE*w_HAve)**two)
    momfl(1) = u_shear**two

    do j=2,ny
      if (momfl(j)<=0.05_mytype*u_shear**two) then
         if (istret==0) then
           h=real(j-1,mytype)*dy +  dy*(0.05_mytype*u_shear**two-momfl(j-1))/(momfl(j)-momfl(j-1)) 
         else if (istret==1) then
           h=yp(j-1) + (yp(j)-yp(j-1))*(0.05_mytype*u_shear**two-momfl(j-1))/(momfl(j)-momfl(j-1)) 
         endif
         exit
      endif
    enddo
    hBL=h/0.95_mytype

    if (nrank==0)   write(*,*)  'boundary layer height', hBL

  end subroutine boundary_height

  !*******************************************************************************
  !
  subroutine postprocess_abl(ux1, uy1, uz1, ep1)
  !
  !*******************************************************************************

    USE MPI
    USE decomp_2d_io
    USE var, only : umean,vmean,wmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
    USE var, only : uvisu
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    USE variables
    !USE param, only: ivisu, ioutput, itime
    USE param
    USE tools
    
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    character(len=30) :: filename
    integer :: out_cntr

    ! Record outflow 
    if (ioutflow.eq.1) then
      out_cntr=mod(itime,ntimesteps)
      if (out_cntr==0) out_cntr=ntimesteps
      call append_outflow(ux1,uy1,uz1,out_cntr)
      if (mod(itime,ntimesteps)==0) then
        call write_outflow(itime/ntimesteps)
      endif
    endif        

    ! Write vorticity as an example of post processing
!    if ((ivisu/=0).and.(mod(itime,ioutput)==0)) then
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
!    
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
  end subroutine postprocess_abl    

end module abl
