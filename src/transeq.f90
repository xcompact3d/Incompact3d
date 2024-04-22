!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module transeq

  use decomp_2d_constants
  use decomp_2d_mpi

  private
  public :: calculate_transeq_rhs

contains
  !############################################################################
  !!  SUBROUTINE: calculate_transeq_rhs
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Calculates the right hand sides of all transport
  !!              equations - momentum, scalar transport, etc.
  !############################################################################
  subroutine calculate_transeq_rhs(drho1,dux1,duy1,duz1,dphi1,rho1,ux1,uy1,uz1,ep1,phi1,divu3)

    use decomp_2d, only : xsize, zsize
    use variables, only : numscalar
    use param, only : ntime, ilmn, nrhotime, ilmn_solve_temp
    use mhd,   only : mhd_active,mhd_equation,calculate_mhd_transeq_rhs

    implicit none

    !! Inputs
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ux1, uy1, uz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), nrhotime), intent(in) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar), intent(in) :: phi1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ep1
    real(mytype), dimension(zsize(1), zsize(2), zsize(3)), intent(in) :: divu3

    !! Outputs
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: drho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime, numscalar) :: dphi1

    !! Momentum equations
    call momentum_rhs_eq(dux1,duy1,duz1,rho1,ux1,uy1,uz1,ep1,phi1,divu3)

    !! Scalar equations
    !! XXX Not yet LMN!!!
    call scalar(dphi1, rho1, ux1, uy1, uz1, phi1)

    if(mhd_active .and. mhd_equation) then
      call calculate_mhd_transeq_rhs(ux1,uy1,uz1)
    endif

    !! Other (LMN, ...)
    if (ilmn) THEN
       if (ilmn_solve_temp) THEN
          call temperature_rhs_eq(drho1, rho1, ux1, uy1, uz1, phi1)
       else
          call continuity_rhs_eq(drho1, rho1, ux1, divu3)
       endif
    endif

  end subroutine calculate_transeq_rhs
  !############################################################################
  !############################################################################
  !!
  !!  subroutine: momentum_rhs_eq
  !!      AUTHOR: ?
  !!    MODIFIED: Kay Schäfer
  !! DESCRIPTION: Calculation of convective and diffusion terms of momentum
  !!              equation
  !!
  !############################################################################
  !############################################################################
  subroutine momentum_rhs_eq(dux1,duy1,duz1,rho1,ux1,uy1,uz1,ep1,phi1,divu3)

    use param
    use variables
    use decomp_2d
    use var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,mu1,mu2,mu3
    use var, only : rho2,ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
    use var, only : rho3,ux3,uy3,uz3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
    use var, only : sgsx1,sgsy1,sgsz1
    use var, only : FTx, FTy, FTz, Fdiscx, Fdiscy, Fdiscz
    use ibm_param, only : ubcx,ubcy,ubcz
    use les, only : compute_SGS
    use mpi

    use case, only : momentum_forcing

    implicit none

    !! INPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    real(mytype),intent(in),dimension(zsize(1),zsize(2),zsize(3)) :: divu3

    !! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1
    
    integer :: i,j,k,is

#ifdef DEBG 
    real(mytype) :: dep, dep1
    integer :: code
#endif

    !SKEW SYMMETRIC FORM
    !WORK X-PENCILS
    if (ilmn) then
      ta1(:,:,:) = rho1(:,:,:,1) * ux1(:,:,:) * ux1(:,:,:)
      tb1(:,:,:) = rho1(:,:,:,1) * ux1(:,:,:) * uy1(:,:,:)
      tc1(:,:,:) = rho1(:,:,:,1) * ux1(:,:,:) * uz1(:,:,:)
    else
      ta1(:,:,:) = ux1(:,:,:) * ux1(:,:,:)
      tb1(:,:,:) = ux1(:,:,:) * uy1(:,:,:)
      tc1(:,:,:) = ux1(:,:,:) * uz1(:,:,:)
    endif

#ifdef DEBG
    dep=maxval(abs(ta1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## SUB momentum_rhs_eq VAR ta1 (uu) MAX ', dep1
#endif

    call derx (td1,ta1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcx*ubcx)
    call derx (te1,tb1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx*ubcy)
    call derx (tf1,tc1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx*ubcz)
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz)

#ifdef DEBG
    dep=maxval(abs(ta1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## SUB momentum_rhs_eq VAR ta1 (du) MAX ', dep1
#endif

    ! Convective terms of x-pencil are stored in tg1,th1,ti1
    if (ilmn) then
      tg1(:,:,:) = td1(:,:,:) + rho1(:,:,:,1) * ux1(:,:,:) * ta1(:,:,:)
      th1(:,:,:) = te1(:,:,:) + rho1(:,:,:,1) * ux1(:,:,:) * tb1(:,:,:)
      ti1(:,:,:) = tf1(:,:,:) + rho1(:,:,:,1) * ux1(:,:,:) * tc1(:,:,:)
    else
      tg1(:,:,:) = td1(:,:,:) + ux1(:,:,:) * ta1(:,:,:)
      th1(:,:,:) = te1(:,:,:) + ux1(:,:,:) * tb1(:,:,:)
      ti1(:,:,:) = tf1(:,:,:) + ux1(:,:,:) * tc1(:,:,:)
    endif
    ! TODO: save the x-convective terms already in dux1, duy1, duz1
#ifdef DEBG
    dep=maxval(abs(tg1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## SUB momentum_rhs_eq VAR tg1 (duu+udu) MAX ', dep1
#endif

    if (ilmn) then
       !! Quasi-skew symmetric terms
       call derx (td1,rho1(:,:,:,1),di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1, zero)
       tg1(:,:,:) = tg1(:,:,:) + ux1(:,:,:) * ux1(:,:,:) * td1(:,:,:)
       th1(:,:,:) = th1(:,:,:) + uy1(:,:,:) * ux1(:,:,:) * td1(:,:,:)
       ti1(:,:,:) = ti1(:,:,:) + uz1(:,:,:) * ux1(:,:,:) * td1(:,:,:)
    endif

    call transpose_x_to_y(ux1,ux2)
    call transpose_x_to_y(uy1,uy2)
    call transpose_x_to_y(uz1,uz2)
#ifdef DEBG
    dep=maxval(abs(ux2))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## SUB momentum_rhs_eq VAR ux2 (transpose) MAX ', dep1
    dep=maxval(abs(uy2))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## SUB momentum_rhs_eq VAR uy2 (transpose) MAX ', dep1
#endif

    if (ilmn) then
       call transpose_x_to_y(rho1(:,:,:,1),rho2)
       call transpose_x_to_y(mu1,mu2)
    else
       rho2(:,:,:) = one
    endif

    !WORK Y-PENCILS
    if (ilmn) then
      td2(:,:,:) = rho2(:,:,:) * ux2(:,:,:) * uy2(:,:,:)
      te2(:,:,:) = rho2(:,:,:) * uy2(:,:,:) * uy2(:,:,:)
      tf2(:,:,:) = rho2(:,:,:) * uz2(:,:,:) * uy2(:,:,:)
    else
      td2(:,:,:) = ux2(:,:,:) * uy2(:,:,:)
      te2(:,:,:) = uy2(:,:,:) * uy2(:,:,:)
      tf2(:,:,:) = uz2(:,:,:) * uy2(:,:,:)
    endif
#ifdef DEBG
    dep=maxval(abs(td2))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## SUB momentum_rhs_eq VAR td2 (uu) MAX ', dep1
#endif

    call dery (tg2,td2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcx*ubcy)
    call dery (th2,te2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcy*ubcy)
    call dery (ti2,tf2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcz*ubcy)
    call dery (td2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
    call dery (te2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
    call dery (tf2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)

#ifdef DEBG
    dep=maxval(abs(td2))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## SUB momentum_rhs_eq VAR td2 (du) MAX ', dep1
#endif

    ! Convective terms of y-pencil in tg2,th2,ti2
    if (ilmn) then
      tg2(:,:,:) = tg2(:,:,:) + rho2(:,:,:) * uy2(:,:,:) * td2(:,:,:)
      th2(:,:,:) = th2(:,:,:) + rho2(:,:,:) * uy2(:,:,:) * te2(:,:,:)
      ti2(:,:,:) = ti2(:,:,:) + rho2(:,:,:) * uy2(:,:,:) * tf2(:,:,:)
    else
      tg2(:,:,:) = tg2(:,:,:) + uy2(:,:,:) * td2(:,:,:)
      th2(:,:,:) = th2(:,:,:) + uy2(:,:,:) * te2(:,:,:)
      ti2(:,:,:) = ti2(:,:,:) + uy2(:,:,:) * tf2(:,:,:)
    endif
#ifdef DEBG
    dep=maxval(abs(tg2))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## SUB momentum_rhs_eq VAR tg2 (duu+udu) MAX ', dep1
#endif


    if (ilmn) then
       !! Quasi-skew symmetric terms
       call dery (te2,rho2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,zero)
       tg2(:,:,:) = tg2(:,:,:) + ux2(:,:,:) * uy2(:,:,:) * te2(:,:,:)
       th2(:,:,:) = th2(:,:,:) + uy2(:,:,:) * uy2(:,:,:) * te2(:,:,:)
       ti2(:,:,:) = ti2(:,:,:) + uz2(:,:,:) * uy2(:,:,:) * te2(:,:,:)
    endif

    call transpose_y_to_z(ux2,ux3)
    call transpose_y_to_z(uy2,uy3)
    call transpose_y_to_z(uz2,uz3)

    !WORK Z-PENCILS
    if (ilmn) then
       call transpose_y_to_z(rho2,rho3)
       call transpose_y_to_z(mu2,mu3)

       td3(:,:,:) = rho3(:,:,:) * ux3(:,:,:) * uz3(:,:,:)
       te3(:,:,:) = rho3(:,:,:) * uy3(:,:,:) * uz3(:,:,:)
       tf3(:,:,:) = rho3(:,:,:) * uz3(:,:,:) * uz3(:,:,:)
    else
       td3(:,:,:) = ux3(:,:,:) * uz3(:,:,:)
       te3(:,:,:) = uy3(:,:,:) * uz3(:,:,:)
       tf3(:,:,:) = uz3(:,:,:) * uz3(:,:,:)
    endif
#ifdef DEBG
    dep=maxval(abs(td3))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## SUB momentum_rhs_eq VAR td3 (uu) MAX ', dep1
#endif

    call derz (tg3,td3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcx*ubcz)
    call derz (th3,te3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcy*ubcz)
    call derz (ti3,tf3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcz*ubcz)
    call derz (td3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
    call derz (te3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcy)
    call derz (tf3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcz)

    ! Convective terms of z-pencil in ta3,tb3,tc3
    if (ilmn) then
      ta3(:,:,:) = tg3(:,:,:) + rho3(:,:,:) * uz3(:,:,:) * td3(:,:,:)
      tb3(:,:,:) = th3(:,:,:) + rho3(:,:,:) * uz3(:,:,:) * te3(:,:,:)
      tc3(:,:,:) = ti3(:,:,:) + rho3(:,:,:) * uz3(:,:,:) * tf3(:,:,:)
    else
      ta3(:,:,:) = tg3(:,:,:) + uz3(:,:,:) * td3(:,:,:)
      tb3(:,:,:) = th3(:,:,:) + uz3(:,:,:) * te3(:,:,:)
      tc3(:,:,:) = ti3(:,:,:) + uz3(:,:,:) * tf3(:,:,:)
    endif

    if (ilmn) then
       !! Quasi-skew symmetric terms
       call derz (tf3,rho3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,zero)
       ta3(:,:,:) = ta3(:,:,:) + ux3(:,:,:) * uz3(:,:,:) * tf3(:,:,:)
       tb3(:,:,:) = tb3(:,:,:) + uy3(:,:,:) * uz3(:,:,:) * tf3(:,:,:)
       tc3(:,:,:) = tc3(:,:,:) + uz3(:,:,:) * uz3(:,:,:) * tf3(:,:,:)

       !! Add the additional divu terms
       ta3(:,:,:) = ta3(:,:,:) + rho3(:,:,:) * ux3(:,:,:) * divu3(:,:,:)
       tb3(:,:,:) = tb3(:,:,:) + rho3(:,:,:) * uy3(:,:,:) * divu3(:,:,:)
       tc3(:,:,:) = tc3(:,:,:) + rho3(:,:,:) * uz3(:,:,:) * divu3(:,:,:)
    endif
#ifdef DEBG
    dep=maxval(abs(ta3))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## SUB momentum_rhs_eq VAR ta3 (duu+udu) MAX ', dep1
#endif

    ! Convective terms of z-pencil are in ta3 -> td3, tb3 -> te3, tc3 -> tf3
    td3(:,:,:) = ta3(:,:,:)
    te3(:,:,:) = tb3(:,:,:)
    tf3(:,:,:) = tc3(:,:,:)

    !DIFFUSIVE TERMS IN Z
    call derzz (ta3,ux3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1,ubcx)
    call derzz (tb3,uy3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1,ubcy)
    call derzz (tc3,uz3,di3,sz,sfz ,ssz ,swz ,zsize(1),zsize(2),zsize(3),0,ubcz)


    ! Add convective and diffusive terms of z-pencil (half for skew-symmetric)
    if (ilmn) then
      td3(:,:,:) = mu3(:,:,:) * xnu*ta3(:,:,:) - half * td3(:,:,:)
      te3(:,:,:) = mu3(:,:,:) * xnu*tb3(:,:,:) - half * te3(:,:,:)
      tf3(:,:,:) = mu3(:,:,:) * xnu*tc3(:,:,:) - half * tf3(:,:,:)
    else
      td3(:,:,:) = xnu*ta3(:,:,:) - half * td3(:,:,:)
      te3(:,:,:) = xnu*tb3(:,:,:) - half * te3(:,:,:)
      tf3(:,:,:) = xnu*tc3(:,:,:) - half * tf3(:,:,:)
    endif

    !WORK Y-PENCILS
    call transpose_z_to_y(td3,td2)
    call transpose_z_to_y(te3,te2)
    call transpose_z_to_y(tf3,tf2)

    ! Convective terms of y-pencil (tg2,th2,ti2) and sum of convective and diffusive terms of z-pencil (td2,te2,tf2) are now in tg2, th2, ti2 (half for skew-symmetric)
    tg2(:,:,:) = td2(:,:,:) - half * tg2(:,:,:)
    th2(:,:,:) = te2(:,:,:) - half * th2(:,:,:)
    ti2(:,:,:) = tf2(:,:,:) - half * ti2(:,:,:)
#ifdef DEBG
    dep=maxval(abs(tg2))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## SUB momentum_rhs_eq VAR tg2 (Conv+Diff)) MAX ', dep1
#endif


    !DIFFUSIVE TERMS IN Y
    if (iimplicit.le.0) then
       !-->for ux
       call deryy (td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1,ubcx)
       if (istret.ne.0) then
          call dery (te2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
          do k = 1,ysize(3)
             do j = 1,ysize(2)
                do i = 1,ysize(1)
                   td2(i,j,k) = td2(i,j,k)*pp2y(j)-pp4y(j)*te2(i,j,k)
                enddo
             enddo
          enddo
       endif

       !-->for uy
       call deryy (te2,uy2,di2,sy,sfy,ssy,swy,ysize(1),ysize(2),ysize(3),0,ubcy)
       if (istret.ne.0) then
          call dery (tf2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
          do k = 1,ysize(3)
             do j = 1,ysize(2)
                do i = 1,ysize(1)
                   te2(i,j,k) = te2(i,j,k)*pp2y(j)-pp4y(j)*tf2(i,j,k)
                enddo
             enddo
          enddo
       endif

       !-->for uz
       call deryy (tf2,uz2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1,ubcz)
       if (istret.ne.0) then
          call dery (tj2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
          do k = 1,ysize(3)
             do j = 1,ysize(2)
                do i = 1,ysize(1)
                   tf2(i,j,k) = tf2(i,j,k)*pp2y(j)-pp4y(j)*tj2(i,j,k)
                enddo
             enddo
          enddo
       endif
    else ! (semi)implicit Y diffusion
       if (istret.ne.0) then

          !-->for ux
          call dery (te2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
          do k=1,ysize(3)
             do j=1,ysize(2)
                do i=1,ysize(1)
                   td2(i,j,k)=-pp4y(j)*te2(i,j,k)
                enddo
             enddo
          enddo
          !-->for uy
          call dery (tf2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
          do k=1,ysize(3)
             do j=1,ysize(2)
                do i=1,ysize(1)
                   te2(i,j,k)=-pp4y(j)*tf2(i,j,k)
                enddo
             enddo
          enddo
          !-->for uz
          call dery (tj2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
          do k=1,ysize(3)
             do j=1,ysize(2)
                do i=1,ysize(1)
                   tf2(i,j,k)=-pp4y(j)*tj2(i,j,k)
                enddo
             enddo
          enddo

       else
       
          td2(:,:,:) = zero
          te2(:,:,:) = zero
          tf2(:,:,:) = zero
          
       endif
    endif
#ifdef DEBG
    dep=maxval(abs(td2))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## SUB momentum_rhs_eq VAR td2 (Diff Y) MAX ', dep1
    dep=maxval(abs(te2))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## SUB momentum_rhs_eq VAR te2 (Diff Y) MAX ', dep1
    dep=maxval(abs(tf2))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## SUB momentum_rhs_eq VAR tf2 (Diff Y) MAX ', dep1
#endif


    ! Add diffusive terms of y-pencil to convective and diffusive terms of y- and z-pencil
    if (ilmn) then
      ta2(:,:,:) = mu2(:,:,:) * xnu*td2(:,:,:) + tg2(:,:,:)
      tb2(:,:,:) = mu2(:,:,:) * xnu*te2(:,:,:) + th2(:,:,:)
      tc2(:,:,:) = mu2(:,:,:) * xnu*tf2(:,:,:) + ti2(:,:,:)
    else
      ta2(:,:,:) = xnu*td2(:,:,:) + tg2(:,:,:)
      tb2(:,:,:) = xnu*te2(:,:,:) + th2(:,:,:)
      tc2(:,:,:) = xnu*tf2(:,:,:) + ti2(:,:,:)
    endif

    !WORK X-PENCILS
    call transpose_y_to_x(ta2,ta1)
    call transpose_y_to_x(tb2,tb1)
    call transpose_y_to_x(tc2,tc1) !diff+conv. terms

    !DIFFUSIVE TERMS IN X
    call derxx (td1,ux1,di1,sx,sfx ,ssx ,swx ,xsize(1),xsize(2),xsize(3),0,ubcx)
    call derxx (te1,uy1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1,ubcy)
    call derxx (tf1,uz1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1,ubcz)

    if (ilmn) then
      td1(:,:,:) = mu1(:,:,:) * xnu * td1(:,:,:)
      te1(:,:,:) = mu1(:,:,:) * xnu * te1(:,:,:)
      tf1(:,:,:) = mu1(:,:,:) * xnu * tf1(:,:,:)
    else
      td1(:,:,:) = xnu * td1(:,:,:)
      te1(:,:,:) = xnu * te1(:,:,:)
      tf1(:,:,:) = xnu * tf1(:,:,:)
    endif
#ifdef DEBG
    dep=maxval(abs(td1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## SUB momentum_rhs_eq VAR td1 (Diff X) MAX ', dep1
    dep=maxval(abs(te1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## SUB momentum_rhs_eq VAR te1 (Diff X) MAX ', dep1
    dep=maxval(abs(tf1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## SUB momentum_rhs_eq VAR tf1 (Diff X) MAX ', dep1
#endif

    !FINAL SUM: DIFF TERMS + CONV TERMS
    dux1(:,:,:,1) = ta1(:,:,:) - half*tg1(:,:,:)  + td1(:,:,:)
    duy1(:,:,:,1) = tb1(:,:,:) - half*th1(:,:,:)  + te1(:,:,:)
    duz1(:,:,:,1) = tc1(:,:,:) - half*ti1(:,:,:)  + tf1(:,:,:)
#ifdef DEBG
    dep=maxval(abs(dux1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## MomRHS dux1 ', dep1
    dep=maxval(abs(duy1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## MomRHS duy1 ', dep1
    dep=maxval(abs(duz1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## MomRHS duz1 ', dep1
#endif
    if (ilmn) then
       call momentum_full_viscstress_tensor(dux1(:,:,:,1), duy1(:,:,:,1), duz1(:,:,:,1), divu3, mu1)
    endif
#ifdef DEBG
    dep=maxval(abs(dux1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## MomRHS VisTau dux1 ', dep1
    dep=maxval(abs(duy1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## MomRHS VisTau duy1 ', dep1
    dep=maxval(abs(duz1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## MomRHS VisTau duz1 ', dep1
#endif

    ! If LES modelling is enabled, add the SGS stresses
    if (ilesmod.ne.0.and.jles.le.3.and.jles.gt.0) then
       call compute_SGS(sgsx1,sgsy1,sgsz1,ux1,uy1,uz1,phi1,ep1)
       dux1(:,:,:,1) = dux1(:,:,:,1) + sgsx1(:,:,:)
       duy1(:,:,:,1) = duy1(:,:,:,1) + sgsy1(:,:,:)
       duz1(:,:,:,1) = duz1(:,:,:,1) + sgsz1(:,:,:)
    endif
#ifdef DEBG
    dep=maxval(abs(dux1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## MomRHS LES dux1 ', dep1
    dep=maxval(abs(duy1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## MomRHS LES duy1 ', dep1
    dep=maxval(abs(duz1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## MomRHS LES duz1 ', dep1
#endif

    if (ilmn) then
      !! Gravity
      if ((Fr**2).gt.zero) then
        call momentum_gravity(dux1, duy1, duz1, rho1(:,:,:,1) - one, one / Fr**2)
      endif
    endif
    if (iscalar.eq.1) then
      do is = 1, numscalar
        call momentum_gravity(dux1, duy1, duz1, phi1(:,:,:,is), ri(is))
      enddo
    endif
#ifdef DEBG
    dep=maxval(abs(dux1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## MomRHS ILMN dux1 ', dep1
    dep=maxval(abs(duy1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## MomRHS ILMN duy1 ', dep1
    dep=maxval(abs(duz1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## MomRHS ILMN duz1 ', dep1
#endif
    !! Additional forcing
    call momentum_forcing(dux1, duy1, duz1, rho1, ux1, uy1, uz1, phi1)
#ifdef DEBG
    dep=maxval(abs(dux1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## MomRHS Forc dux1 ', dep1
    dep=maxval(abs(duy1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## MomRHS Forc duy1 ', dep1
    dep=maxval(abs(duz1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## MomRHS Forc duz1 ', dep1
#endif

    !! Turbine forcing
    if (iturbine.eq.1) then
       dux1(:,:,:,1)=dux1(:,:,:,1)+FTx(:,:,:)/rho_air
       duy1(:,:,:,1)=duy1(:,:,:,1)+FTy(:,:,:)/rho_air
       duz1(:,:,:,1)=duz1(:,:,:,1)+FTz(:,:,:)/rho_air
    else if (iturbine.eq.2) then
       dux1(:,:,:,1)=dux1(:,:,:,1)+Fdiscx(:,:,:)/rho_air
       duy1(:,:,:,1)=duy1(:,:,:,1)+Fdiscy(:,:,:)/rho_air
       duz1(:,:,:,1)=duz1(:,:,:,1)+Fdiscz(:,:,:)/rho_air
    endif
#ifdef DEBG
    dep=maxval(abs(dux1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## MomRHS Turb dux1 ', dep1
    dep=maxval(abs(duy1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## MomRHS Turb duy1 ', dep1
    dep=maxval(abs(duz1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## MomRHS Turb duz1 ', dep1
#endif

    if (itrip == 1) then
       !call tripping(tb1,td1)
       call tbl_tripping(duy1(:,:,:,1),td1)
       if ((nrank==0).and.(mod(itime,ilist)==0)) write(*,*) 'TRIPPING!!'
    endif
#ifdef DEBG
    dep=maxval(abs(dux1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## MomRHS Trip dux1 ', dep1
    dep=maxval(abs(duy1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## MomRHS Trip duy1 ', dep1
    dep=maxval(abs(duz1))
    call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (nrank == 0) write(*,*)'## MomRHS Trip duz1 ', dep1
#endif

  end subroutine momentum_rhs_eq
  !############################################################################
  !############################################################################
  !!
  !!  SUBROUTINE: momentum_full_viscstress_tensor
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: In an incompressible flow the viscous stress
  !!              tensor reduces to
  !!                d2u^j / {dx^i}^2
  !!              however if div(u) != 0 we have
  !!                d2u^j / {dx^i}^2 + 1/3 d/dx^j div(u)
  !!              and further if \mu != const.
  !!                \mu (d2u^j / {dx^i}^2 + 1/3 d/dx^j div(u))
  !!                  + d\mu/dx^i (du^j/dx^i + du^i/dx^j
  !!                  - 2/3 div(u) \delta^{ij})
  !!              This subroutine computes the additional
  !!              contributions not accounted for in the
  !!              incompressible solver.
  !!
  !############################################################################
  subroutine momentum_full_viscstress_tensor(dux1, duy1, duz1, divu3, mu1)

    use param
    use variables
    use decomp_2d

    use var, only : ux1,uy1,uz1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    use var, only : ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,di2
    use var, only : ux3,uy3,uz3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
    use ibm_param
    
    implicit none

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: dux1, duy1, duz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: mu1
    real(mytype), dimension(zsize(1), zsize(2), zsize(3)), intent(in) :: divu3

    real(mytype) :: one_third

    one_third = one / three

    call derz (tc3,divu3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,zero)
    call transpose_z_to_y(tc3, tc2)
    call transpose_z_to_y(divu3, th2)

    call dery(tb2,th2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,zero)
    call transpose_y_to_x(tb2, te1)
    call transpose_y_to_x(tc2, tf1)
    call transpose_y_to_x(th2, tg1)

    call derx(td1,tg1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,zero)

    dux1(:,:,:) = dux1(:,:,:) + mu1(:,:,:) * one_third * xnu * td1(:,:,:)
    duy1(:,:,:) = duy1(:,:,:) + mu1(:,:,:) * one_third * xnu * te1(:,:,:)
    duz1(:,:,:) = duz1(:,:,:) + mu1(:,:,:) * one_third * xnu * tf1(:,:,:)

    !! Variable viscosity part
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz)
    call derx (td1,mu1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,zero)
    ta1(:,:,:) = two * ta1(:,:,:) - (two * one_third) * tg1(:,:,:)

    ta1(:,:,:) = td1(:,:,:) * ta1(:,:,:)
    tb1(:,:,:) = td1(:,:,:) * tb1(:,:,:)
    tc1(:,:,:) = td1(:,:,:) * tc1(:,:,:)

    call transpose_x_to_y(ta1, ta2)
    call transpose_x_to_y(tb1, tb2)
    call transpose_x_to_y(tc1, tc2)
    call transpose_x_to_y(td1, tg2)

    call dery (td2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
    call dery (te2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
    call dery (tf2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
    te2(:,:,:) = two * te2(:,:,:) - (two * one_third) * th2(:,:,:)

    call transpose_x_to_y(mu1, ti2)
    call dery (th2,ti2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,zero)

    ta2(:,:,:) = ta2(:,:,:) + th2(:,:,:) * td2(:,:,:)
    tb2(:,:,:) = tb2(:,:,:) + th2(:,:,:) * te2(:,:,:) + tg2(:,:,:) * td2(:,:,:)
    tc2(:,:,:) = tc2(:,:,:) + th2(:,:,:) * tf2(:,:,:)

    call transpose_y_to_z(ta2, ta3)
    call transpose_y_to_z(tb2, tb3)
    call transpose_y_to_z(tc2, tc3)
    call transpose_y_to_z(tg2, tg3) !! dmudx
    call transpose_y_to_z(th2, th3) !! dmudy
    call transpose_y_to_z(ti2, ti3) !! mu

    call derz (td3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
    call derz (te3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcy)
    call derz (tf3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcz)
    tf3(:,:,:) = two * tf3(:,:,:) - (two * one_third) * divu3(:,:,:)

    tc3(:,:,:) = tc3(:,:,:) + tg3(:,:,:) * td3(:,:,:) + th3(:,:,:) * te3(:,:,:)

    call derz (th3,ti3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,zero)

    ta3(:,:,:) = ta3(:,:,:) + th3(:,:,:) * td3(:,:,:)
    tb3(:,:,:) = tb3(:,:,:) + th3(:,:,:) * te3(:,:,:)
    tc3(:,:,:) = tc3(:,:,:) + th3(:,:,:) * tf3(:,:,:)

    call transpose_z_to_y(ta3, ta2)
    call transpose_z_to_y(tb3, tb2)
    call transpose_z_to_y(tc3, tc2)
    call transpose_z_to_y(th3, ti2) !! dmudz

    tb2(:,:,:) = tb2(:,:,:) + ti2(:,:,:) * tf2(:,:,:)

    call transpose_y_to_x(ta2, ta1)
    call transpose_y_to_x(tb2, tb1)
    call transpose_y_to_x(tc2, tc1)
    call transpose_y_to_x(th2, te1) !! dmudy
    call transpose_y_to_x(ti2, tf1) !! dmudz

    call derx (th1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
    call derx (ti1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz)
    ta1(:,:,:) = ta1(:,:,:) + te1(:,:,:) * th1(:,:,:) + tf1(:,:,:) * ti1(:,:,:)

    dux1(:,:,:) = dux1(:,:,:) + xnu * ta1(:,:,:)
    duy1(:,:,:) = duy1(:,:,:) + xnu * tb1(:,:,:)
    duz1(:,:,:) = duz1(:,:,:) + xnu * tc1(:,:,:)


  end subroutine momentum_full_viscstress_tensor
  !############################################################################
  !############################################################################
  subroutine momentum_gravity(dux1, duy1, duz1, peculiar_density1, richardson)

    use decomp_2d
    use param
    use variables

    implicit none

    !! Inputs
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: peculiar_density1
    real(mytype), intent(in) :: richardson

    !! InOut
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1

    !! Locals
    integer :: istart, jstart, kstart
    integer :: iend, jend, kend
    integer :: i, j, k

    ! Return directly if richardson is zero
    if (abs(richardson) < tiny(richardson)) return

    !! X-gravity
    if (abs(gravx) > tiny(gravx)) then
    if ((nclx1.eq.0).and.(nclxn.eq.0)) then
       istart = 1
       iend = xsize(1)
    else
       istart = 2
       iend = xsize(1) - 1
    endif
    if ((xstart(2).eq.1).and.(ncly1.eq.2)) then
       jstart = 2
    else
       jstart = 1
    endif
    if ((xend(2).eq.ny).and.(nclyn.eq.2)) then
       jend = xsize(2) - 1
    else
       jend = xsize(2)
    endif
    if ((xstart(3).eq.1).and.(nclz1.eq.2)) then
       kstart = 2
    else
       kstart = 1
    endif
    if ((xend(3).eq.nz).and.(nclzn.eq.2)) then
       kend = xsize(3) - 1
    else
       kend = xsize(3)
    endif

    do k = kstart, kend
       do j = jstart, jend
          do i = istart, iend
             dux1(i, j, k, 1) = dux1(i, j, k, 1) + peculiar_density1(i, j, k) * richardson * gravx
          enddo
       enddo
    enddo
    endif

    !! Y-gravity
    if (abs(gravy) > tiny(gravy)) then
    if (nclx1.eq.2) then
       istart = 2
    else
       istart = 1
    endif
    if (nclxn.eq.2) then
       iend = xsize(1) - 1
    else
       iend = xsize(1)
    endif
    if ((xstart(2).eq.1).and.(ncly1.ne.0)) then
       jstart = 2
    else
       jstart = 1
    endif
    if ((xend(2).eq.ny).and.(nclyn.ne.0)) then
       jend = xsize(2) - 1
    else
       jend = xsize(2)
    endif
    if ((xstart(3).eq.1).and.(nclz1.eq.2)) then
       kstart = 2
    else
       kstart = 1
    endif
    if ((xend(3).eq.nz).and.(nclzn.eq.2)) then
       kend = xsize(3) - 1
    else
       kend = xsize(3)
    endif
    do k = kstart, kend
       do j = jstart, jend
          do i = istart, iend
             duy1(i, j, k, 1) = duy1(i, j, k, 1) + peculiar_density1(i, j, k) * richardson * gravy
          enddo
       enddo
    enddo
    endif

    !! Z-gravity
    if (abs(gravz) > tiny(gravz)) then
    if (nclx1.eq.2) then
       istart = 2
    else
       istart = 1
    endif
    if (nclxn.eq.2) then
       iend = xsize(1) - 1
    else
       iend = xsize(1)
    endif
    if ((xstart(2).eq.1).and.(ncly1.eq.2)) then
       jstart = 2
    else
       jstart = 1
    endif
    if ((xend(2).eq.ny).and.(nclyn.eq.2)) then
       jend = xsize(2) - 1
    else
       jend = xsize(2)
    endif
    if ((xstart(3).eq.1).and.(nclz1.ne.0)) then
       kstart = 2
    else
       kstart = 1
    endif
    if ((xend(3).eq.nz).and.(nclzn.ne.0)) then
       kend = xsize(3) - 1
    else
       kend = xsize(3)
    endif
    do k = kstart, kend
       do j = jstart, jend
          do i = istart, iend
             duz1(i, j, k, 1) = duz1(i, j, k, 1) + peculiar_density1(i, j, k) * richardson * gravz
          enddo
       enddo
    enddo
    endif

  end subroutine momentum_gravity
  !############################################################################
  !############################################################################
  subroutine scalar_transport_eq(dphi1, rho1, ux1, uy1, uz1, phi1, schmidt, is)

    use param
    use variables
    use decomp_2d
    use case, only : scalar_forcing

    use var, only : ta1,tb1,tc1,di1
    use var, only : rho2,uy2,ta2,tb2,tc2,td2,te2,di2
    use var, only : rho3,uz3,ta3,tb3,tc3,td3,di3

    implicit none

    !! INPUTS
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: phi1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    real(mytype), intent(in) :: schmidt
    integer, optional, intent(in) :: is

    !! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dphi1

    !! LOCALS
    logical :: evensc, skewsc
    integer :: i, j, k
    real(mytype) :: xalpha

    evensc = .true.
    if (present(is)) then
      if (.not.sc_even(is)) then
        evensc = .false.
      endif
    endif
    skewsc = .false.
    if (present(is)) then
      if (sc_skew(is)) then
        skewsc = .true.
      endif
    endif

    xalpha = xnu/schmidt

    !X PENCILS
    if (skewsc) ta1(:,:,:) = ux1(:,:,:) * phi1(:,:,:)
    if (evensc) then
      call derxS (tb1,phi1(:,:,:),di1,sx,ffxpS,fsxpS,fwxpS,xsize(1),xsize(2),xsize(3),1,zero)
      if (skewsc) call derxS (tc1,ta1,di1,sx,ffxS,fsxS,fwxS,xsize(1),xsize(2),xsize(3),0,zero)
    else
      call derxS (tb1,phi1(:,:,:),di1,sx,ffxS,fsxS,fwxS,xsize(1),xsize(2),xsize(3),0,zero)
      if (skewsc) call derxS (tc1,ta1,di1,sx,ffxpS,fsxpS,fwxpS,xsize(1),xsize(2),xsize(3),1,zero) 
    endif
    if (ilmn) then
      tb1(:,:,:) = rho1(:,:,:,1) * ux1(:,:,:) * tb1(:,:,:)
      if (skewsc) tc1(:,:,:) = rho1(:,:,:,1) * tc1(:,:,:)
    else
      tb1(:,:,:) = ux1(:,:,:) * tb1(:,:,:)
    endif
    if (skewsc) then
      tb1(:,:,:) = tb1(:,:,:) + half * (tc1(:,:,:) - tb1(:,:,:))
    endif

    if (evensc) then
      call derxxS (ta1,phi1(:,:,:),di1,sx,sfxpS,ssxpS,swxpS,xsize(1),xsize(2),xsize(3),1,zero)
    else
      call derxxS (ta1,phi1(:,:,:),di1,sx,sfxS,ssxS,swxS,xsize(1),xsize(2),xsize(3),0,zero)
    endif

    ! Add convective and diffusive scalar terms of x-pencil
    ta1(:,:,:) = xalpha*ta1(:,:,:) - tb1(:,:,:)

    call transpose_x_to_y(phi1(:,:,:),td2(:,:,:))

    !Y PENCILS
    if (skewsc) tb2(:,:,:) = uy2(:,:,:) * td2(:,:,:)
    ! Explicit viscous diffusion
    if (iimplicit.le.0) then
      if (evensc) then
        call deryS (tc2,td2(:,:,:),di2,sy,ffypS,fsypS,fwypS,ppy,ysize(1),ysize(2),ysize(3),1,zero)
        if (skewsc) call deryS (te2,tb2,di2,sy,ffyS,fsyS,fwyS,ppy,ysize(1),ysize(2),ysize(3),0,zero)
        call deryyS (ta2,td2(:,:,:),di2,sy,sfypS,ssypS,swypS,ysize(1),ysize(2),ysize(3),1,zero)
      else
        call deryS (tc2,td2(:,:,:),di2,sy,ffyS,fsyS,fwyS,ppy,ysize(1),ysize(2),ysize(3),0,zero)
        if (skewsc) call deryS (te2,tb2,di2,sy,ffypS,fsypS,fwypS,ppy,ysize(1),ysize(2),ysize(3),1,zero)
        call deryyS (ta2,td2(:,:,:),di2,sy,sfyS,ssyS,swyS,ysize(1),ysize(2),ysize(3),0,zero)
      endif

      if (istret.ne.0) then
         do k = 1,ysize(3)
            do j = 1,ysize(2)
               do i = 1,ysize(1)
                  ta2(i,j,k) = ta2(i,j,k)*pp2y(j)-pp4y(j)*tc2(i,j,k)
               enddo
            enddo
         enddo
      endif

    ! (semi)implicit Y viscous diffusion
    else
      if (evensc) then
        call deryS (tc2,td2(:,:,:),di2,sy,ffypS,fsypS,fwypS,ppy,ysize(1),ysize(2),ysize(3),1,zero)
        if (skewsc) call deryS (te2,tb2,di2,sy,ffyS,fsyS,fwyS,ppy,ysize(1),ysize(2),ysize(3),0,zero)
      else      
        call deryS (tc2,td2(:,:,:),di2,sy,ffyS,fsyS,fwyS,ppy,ysize(1),ysize(2),ysize(3),0,zero)
        if (skewsc) call deryS (te2,tb2,di2,sy,ffypS,fsypS,fwypS,ppy,ysize(1),ysize(2),ysize(3),1,zero)
      endif                                                                             
       
      if (istret.ne.0) then
         do k = 1,ysize(3)                                                              
            do j = 1,ysize(2)
               do i = 1,ysize(1)
                  ta2(i,j,k) = -pp4y(j)*tc2(i,j,k)                    
               enddo
            enddo
         enddo
      else
         ta2(:,:,:) = zero
      endif     

    endif

    if (ilmn) then
       tb2(:,:,:) = rho2(:,:,:) * uy2(:,:,:) * tc2(:,:,:)
       if (skewsc) te2(:,:,:) = rho2(:,:,:) * te2(:,:,:)
    else
       tb2(:,:,:) = uy2(:,:,:) * tc2(:,:,:)
    endif
    if (skewsc) then
      tb2(:,:,:) = tb2(:,:,:) + half * (te2(:,:,:) - tb2(:,:,:))
    endif

    ! Add convective and diffusive scalar terms of y-pencil
    tc2(:,:,:) = xalpha*ta2(:,:,:) - tb2(:,:,:)

    call transpose_y_to_z(td2(:,:,:),td3(:,:,:))

    !Z PENCILS
    if (skewsc) ta3(:,:,:) = uz3(:,:,:) * td3(:,:,:)
    if (evensc) then
      call derzS (tb3,td3(:,:,:),di3,sz,ffzpS,fszpS,fwzpS,zsize(1),zsize(2),zsize(3),1,zero)
      if (skewsc) call derzS (tc3,ta3,di3,sz,ffzS,fszS,fwzS,zsize(1),zsize(2),zsize(3),0,zero)
    else
      call derzS (tb3,td3(:,:,:),di3,sz,ffzS,fszS,fwzS,zsize(1),zsize(2),zsize(3),0,zero)
      if (skewsc) call derzS (tc3,ta3,di3,sz,ffzpS,fszpS,fwzpS,zsize(1),zsize(2),zsize(3),1,zero)
    endif

    ! convective terms
    if (ilmn) then
      tb3(:,:,:) = rho3(:,:,:) * uz3(:,:,:) * tb3(:,:,:)
      if (skewsc) tc3(:,:,:) = rho3(:,:,:) * tc3(:,:,:)
    else
      tb3(:,:,:) = uz3(:,:,:) * tb3(:,:,:)
    endif
    if (skewsc) then
      tb3(:,:,:) = tb3(:,:,:) + half * (tc3(:,:,:) - tb3(:,:,:))
    endif

    ! diffusive terms
    if (evensc) then
      call derzzS (ta3,td3(:,:,:),di3,sz,sfzpS,sszpS,swzpS,zsize(1),zsize(2),zsize(3),1,zero)
    else
      call derzzS (ta3,td3(:,:,:),di3,sz,sfzS,sszS,swzS,zsize(1),zsize(2),zsize(3),0,zero)
    endif

    ! Add convective and diffusive scalar terms of z-pencil
    ta3(:,:,:) = xalpha*ta3(:,:,:) - tb3(:,:,:)

    call transpose_z_to_y(ta3,ta2)

    !Y PENCILS
    ! Add convective and diffusive scalar terms of z-pencil to y-pencil
    tc2(:,:,:) = tc2(:,:,:) + ta2(:,:,:)

    call transpose_y_to_x(tc2,tc1)

    !X PENCILS
    ! Add convective and diffusive scalar terms to final sum
    dphi1(:,:,:,1) = ta1(:,:,:) + tc1(:,:,:)

    !! Additional forcing
    call scalar_forcing(dphi1, rho1, ux1, uy1, uz1, phi1)

    !! XXX We have computed rho dphidt, want dphidt
    if (ilmn) then
      dphi1(:,:,:,1) = dphi1(:,:,:,1) / rho1(:,:,:,1)
    endif

  endsubroutine scalar_transport_eq
  !############################################################################
  !############################################################################
  subroutine scalar(dphi1, rho1, ux1, uy1, uz1, phi1)

    use param
    use variables
    use decomp_2d
    use les
    use var, only : sgsphi1, nut1
    USE param, only : zero

    implicit none

    !! INPUTS
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1

    !! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime,numscalar) :: dphi1

    !! LOCALS
    integer :: is

    !!=====================================================================
    !! XXX It is assumed that ux,uy,uz are already updated in all pencils!
    !!=====================================================================
    do is = 1, numscalar

       if (is.ne.primary_species) then
          !! For mass fractions enforce primary species Y_p = 1 - sum_s Y_s
          !! So don't solve a transport equation
          call scalar_transport_eq(dphi1(:,:,:,:,is), rho1, ux1, uy1, uz1, phi1(:,:,:,is), sc(is), is=is)
          if (uset(is).ne.zero) then
             call scalar_settling(dphi1, phi1(:,:,:,is), is)
          endif
          ! If LES modelling is enabled, add the SGS stresses
          if (ilesmod.ne.0.and.jles.le.3.and.jles.gt.0) then
             call sgs_scalar_nonconservative(sgsphi1(:,:,:,is),nut1,phi1(:,:,:,is),is)
             dphi1(:,:,:,1,is) = dphi1(:,:,:,1,is) + sgsphi1(:,:,:,is)
          endif
       endif

    end do !loop numscalar

    if (primary_species.ge.1) then
       !! Compute rate of change of primary species
       dphi1(:,:,:,1,primary_species) = zero
       do is = 1, numscalar
          if (is.ne.primary_species) then
             dphi1(:,:,:,1,primary_species) = dphi1(:,:,:,1,primary_species) - dphi1(:,:,:,1,is)
          endif
       enddo
    endif

    ! phi1, phi2 and phi3 are no longer synchronized
    sync_scal_needed = .true.

  end subroutine scalar
  !############################################################################
  !############################################################################
  subroutine scalar_settling(dphi1, phi1, is)

    use param
    use variables
    use decomp_2d

    use var, only : ta1, di1
    use var, only : ta2, tb2, di2, phi2 => tc2
    use var, only : ta3, di3, phi3 => tb3

    implicit none

    !! INPUTS
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: phi1
    integer,intent(in) :: is

    !! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime,numscalar) :: dphi1

    call transpose_x_to_y(phi1, phi2)
    call transpose_y_to_z(phi2, phi3)

    call derz (ta3, phi3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, zero)
    ta3(:,:,:) = uset(is) * gravz * ta3(:,:,:)

    call transpose_z_to_y(ta3, tb2)

    call dery (ta2, phi2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, zero)
    ta2(:,:,:) = uset(is) * gravy * ta2(:,:,:)
    ta2(:,:,:) = ta2(:,:,:) + tb2(:,:,:)

    call transpose_y_to_x(ta2, ta1)

    dphi1(:,:,:,1,is) = dphi1(:,:,:,1,is) - ta1(:,:,:)
    call derx (ta1, phi1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, zero)
    dphi1(:,:,:,1,is) = dphi1(:,:,:,1,is) - uset(is) * gravx * ta1(:,:,:)

  endsubroutine scalar_settling
  !############################################################################
  !############################################################################
  subroutine temperature_rhs_eq(drho1, rho1, ux1, uy1, uz1, phi1)

    use param
    use variables
    use decomp_2d

    use var, only : te1, tb1

    implicit none

    !! INPUTS
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    !! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: drho1

    !! Get temperature
    call calc_temp_eos(te1, rho1(:,:,:,1), phi1, tb1, xsize(1), xsize(2), xsize(3))

    !!=====================================================================
    !! XXX It is assumed that ux,uy,uz are already updated in all pencils!
    !!=====================================================================
    call scalar_transport_eq(drho1, rho1, ux1, uy1, uz1, te1, prandtl)

  end subroutine temperature_rhs_eq
  !############################################################################
  !############################################################################
  subroutine continuity_rhs_eq(drho1, rho1, ux1, divu3)

    use decomp_2d, only : xsize, ysize, zsize
    use decomp_2d, only : transpose_z_to_y, transpose_y_to_x
    use param, only : ntime, nrhotime, ibirman_eos, zero
    use param, only : xnu, prandtl
    use param, only : iimplicit
    use variables

    use var, only : ta1, tb1, di1
    use var, only : rho2, uy2, ta2, tb2, di2
    use var, only : rho3, uz3, ta3, di3

    implicit none

    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1
    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    real(mytype), intent(in), dimension(zsize(1), zsize(2), zsize(3)) :: divu3

    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: drho1

    real(mytype) :: invpe

    invpe = xnu / prandtl

    !! XXX All variables up to date - no need to transpose

    call derz (ta3, rho3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, zero)
    ta3(:,:,:) = uz3(:,:,:) * ta3(:,:,:) + rho3(:,:,:) * divu3(:,:,:)

    call transpose_z_to_y(ta3, tb2)
    call dery (ta2, rho2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, zero)
    ta2(:,:,:) = uy2(:,:,:) * ta2(:,:,:) + tb2(:,:,:)

    call transpose_y_to_x(ta2, ta1)
    call derx (drho1(:,:,:,1), rho1(:,:,:,1), &
         di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, zero)
    drho1(:,:,:,1) = -(ux1(:,:,:) * drho1(:,:,:,1) + ta1(:,:,:))

    if (ibirman_eos) THEN !! Add a diffusion term
       call derzz (ta3,rho3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1, zero)
       call transpose_z_to_y(ta3, tb2)

       iimplicit = -iimplicit
       call deryy (ta2,rho2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1, zero)
       iimplicit = -iimplicit
       ta2(:,:,:) = ta2(:,:,:) + tb2(:,:,:)
       call transpose_y_to_x(ta2, ta1)

       call derxx (tb1,rho1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1, zero)
       ta1(:,:,:) = ta1(:,:,:) + tb1(:,:,:)

       drho1(:,:,:,1) = drho1(:,:,:,1) + invpe * ta1(:,:,:)
    endif

  end subroutine continuity_rhs_eq
  !############################################################################
  !############################################################################
end module transeq
