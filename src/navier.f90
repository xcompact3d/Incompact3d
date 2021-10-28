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
module navier

  implicit none

  private

  public :: solve_poisson, divergence, calc_divu_constraint
  public :: pre_correc, cor_vel
  public :: lmn_t_to_rho_trans, momentum_to_velocity, velocity_to_momentum
  public :: gradp, tbl_flrt

contains
  !############################################################################
  !!  SUBROUTINE: solve_poisson
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Takes the intermediate momentum field as input,
  !!              computes div and solves pressure-Poisson equation.
  !############################################################################
  SUBROUTINE solve_poisson(pp3, px1, py1, pz1, rho1, ux1, uy1, uz1, ep1, drho1, divu3)

    USE decomp_2d, ONLY : mytype, xsize, zsize, ph1, nrank
    USE decomp_2d_poisson, ONLY : poisson
    USE var, ONLY : nzmsize
    USE var, ONLY : dv3
    USE param, ONLY : ntime, nrhotime, npress
    USE param, ONLY : ilmn, ivarcoeff, zero, one 

    implicit none

    !! Inputs
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ux1, uy1, uz1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ep1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime), INTENT(IN) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime), INTENT(IN) :: drho1
    REAL(mytype), DIMENSION(zsize(1), zsize(2), zsize(3)), INTENT(IN) :: divu3

    !! Outputs
    REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: px1, py1, pz1

    !! Locals
    INTEGER :: nlock, poissiter
    LOGICAL :: converged
    REAL(mytype) :: atol, rtol, rho0, divup3norm
#ifdef DEBG
    real(mytype) avg_param
#endif

    nlock = 1 !! Corresponds to computing div(u*)
    converged = .FALSE.
    poissiter = 0
    rho0 = one
#ifdef DOUBLE_PREC
    atol = 1.0e-14_mytype !! Absolute tolerance for Poisson solver
    rtol = 1.0e-14_mytype !! Relative tolerance for Poisson solver
#else
    atol = 1.0e-9_mytype !! Absolute tolerance for Poisson solver
    rtol = 1.0e-9_mytype !! Relative tolerance for Poisson solver
#endif

    IF (ilmn.AND.ivarcoeff) THEN
       !! Variable-coefficient Poisson solver works on div(u), not div(rho u)
       !! rho u -> u
       CALL momentum_to_velocity(rho1, ux1, uy1, uz1)
    ENDIF

    CALL divergence(pp3(:,:,:,1),rho1,ux1,uy1,uz1,ep1,drho1,divu3,nlock)
    IF (ilmn.AND.ivarcoeff) THEN
       dv3(:,:,:) = pp3(:,:,:,1)
    ENDIF

    do while(.not.converged)
       if (ivarcoeff) then
          !! Test convergence
          CALL test_varcoeff(converged, divup3norm, pp3, dv3, atol, rtol, poissiter)

          IF (.NOT.converged) THEN
             !! Evaluate additional RHS terms
             CALL calc_varcoeff_rhs(pp3(:,:,:,1), rho1, px1, py1, pz1, dv3, drho1, ep1, divu3, rho0, &
                  poissiter)
          ENDIF

       ENDIF

       IF (.NOT.converged) THEN
#ifdef DEBG
          avg_param = zero
          call avg3d (pp3(:,:,:,1), avg_param)
          if (nrank == 0) write(*,*)'## Solve Poisson before1 pp3', avg_param
#endif
          CALL poisson(pp3(:,:,:,1))
#ifdef DEBG
          avg_param = zero
          call avg3d (pp3(:,:,:,1), avg_param)
          if (nrank == 0) write(*,*)'## Solve Poisson after call  pp3', avg_param
#endif

          !! Need to update pressure gradient here for varcoeff
          CALL gradp(px1,py1,pz1,pp3(:,:,:,1))
#ifdef DEBG
          avg_param = zero
          call avg3d (pp3(:,:,:,1), avg_param)
          if (nrank == 0) write(*,*)'## Solve Poisson pp3', avg_param
          avg_param = zero
          call avg3d (px1, avg_param)
          if (nrank == 0) write(*,*)'## Solve Poisson px', avg_param
          avg_param = zero
          call avg3d (py1, avg_param)
          if (nrank == 0) write(*,*)'## Solve Poisson py', avg_param
          avg_param = zero
          call avg3d (pz1, avg_param)
          if (nrank == 0) write(*,*)'## Solve Poisson pz', avg_param
#endif
         

          IF ((.NOT.ilmn).OR.(.NOT.ivarcoeff)) THEN
             !! Once-through solver
             !! - Incompressible flow
             !! - LMN - constant-coefficient solver
             converged = .TRUE.
          ENDIF
       ENDIF

       poissiter = poissiter + 1
    enddo

    IF (ilmn.AND.ivarcoeff) THEN
       !! Variable-coefficient Poisson solver works on div(u), not div(rho u)
       !! u -> rho u
       CALL velocity_to_momentum(rho1, ux1, uy1, uz1)
    ENDIF

  END SUBROUTINE solve_poisson
  !############################################################################
  !!
  !!  SUBROUTINE: lmn_t_to_rho_trans
  !! DESCRIPTION: Converts the temperature transient to the density transient
  !!              term. This is achieved by application of EOS and chain rule.
  !!      INPUTS: dtemp1 - the RHS of the temperature equation.
  !!                rho1 - the density field.
  !!     OUTPUTS:  drho1 - the RHS of the density equation.
  !!      AUTHOR: Paul Bartholomew
  !!
  !############################################################################
  SUBROUTINE lmn_t_to_rho_trans(drho1, dtemp1, rho1, dphi1, phi1)

    USE decomp_2d
    USE param, ONLY : zero
    USE param, ONLY : imultispecies, massfrac, mol_weight
    USE param, ONLY : ntime
    USE var, ONLY : numscalar
    USE var, ONLY : ta1, tb1

    IMPLICIT NONE

    !! INPUTS
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3)) :: dtemp1, rho1
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), numscalar) :: phi1
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), ntime, numscalar) :: dphi1

    !! OUTPUTS
    REAL(mytype), INTENT(OUT), DIMENSION(xsize(1), xsize(2), xsize(3)) :: drho1

    !! LOCALS
    INTEGER :: is

    drho1(:,:,:) = zero

    IF (imultispecies) THEN
       DO is = 1, numscalar
          IF (massfrac(is)) THEN
             drho1(:,:,:) = drho1(:,:,:) - dphi1(:,:,:,1,is) / mol_weight(is)
          ENDIF
       ENDDO

       ta1(:,:,:) = zero !! Mean molecular weight
       DO is = 1, numscalar
          IF (massfrac(is)) THEN
             ta1(:,:,:) = ta1(:,:,:) + phi1(:,:,:,is) / mol_weight(is)
          ENDIF
       ENDDO
       drho1(:,:,:) = drho1(:,:,:) / ta1(:,:,:)  !! XXX ta1 is the inverse molecular weight
    ENDIF

    CALL calc_temp_eos(ta1, rho1, phi1, tb1, xsize(1), xsize(2), xsize(3))
    drho1(:,:,:) = drho1(:,:,:) - dtemp1(:,:,:) / ta1(:,:,:)

    drho1(:,:,:) = rho1(:,:,:) * drho1(:,:,:)

  ENDSUBROUTINE lmn_t_to_rho_trans
  !############################################################################
  !subroutine COR_VEL
  !Correction of u* by the pressure gradient to get a divergence free
  !field
  ! input : px,py,pz
  ! output : ux,uy,uz
  !written by SL 2018
  !############################################################################
  subroutine cor_vel (ux,uy,uz,px,py,pz)

    USE decomp_2d
    USE variables
    USE param

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: px,py,pz
#ifdef DEBG
    real(mytype) avg_param
#endif

#ifdef DEBG
    avg_param = zero
    call avg3d (ux, avg_param)
    if (nrank == 0) write(*,*)'## Cor Vel ux', avg_param
    avg_param = zero
    call avg3d (uy, avg_param)
    if (nrank == 0) write(*,*)'## Cor Vel uy', avg_param
    avg_param = zero
    call avg3d (uz, avg_param)
    if (nrank == 0) write(*,*)'## Cor Vel uz', avg_param
    avg_param = zero
    call avg3d (px, avg_param)
    if (nrank == 0) write(*,*)'## Cor Vel px', avg_param
    avg_param = zero
    call avg3d (py, avg_param)
    if (nrank == 0) write(*,*)'## Cor Vel py', avg_param
    avg_param = zero
    call avg3d (pz, avg_param)
    if (nrank == 0) write(*,*)'## Cor Vel pz', avg_param
#endif

    ux(:,:,:)=ux(:,:,:)-px(:,:,:)
    uy(:,:,:)=uy(:,:,:)-py(:,:,:)
    uz(:,:,:)=uz(:,:,:)-pz(:,:,:)

    sync_vel_needed = .true.

    return
  end subroutine cor_vel
  !############################################################################
  !subroutine DIVERGENCe
  !Calculation of div u* for nlock=1 and of div u^{n+1} for nlock=2
  ! input : ux1,uy1,uz1,ep1 (on velocity mesh)
  ! output : pp3 (on pressure mesh)
  !written by SL 2018
  !############################################################################
  subroutine divergence (pp3,rho1,ux1,uy1,uz1,ep1,drho1,divu3,nlock)

    USE param
    USE decomp_2d
    USE variables
    USE var, ONLY: ta1, tb1, tc1, pp1, pgy1, pgz1, di1, &
         duxdxp2, uyp2, uzp2, duydypi2, upi2, ta2, dipp2, &
         duxydxyp3, uzp3, po3, dipp3, nxmsize, nymsize, nzmsize
    USE MPI
    USE ibm_param

    implicit none

    !  TYPE(DECOMP_INFO) :: ph1,ph3,ph4

    !X PENCILS NX NY NZ  -->NXM NY NZ
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime),intent(in) :: drho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime),intent(in) :: rho1
    !Z PENCILS NXM NYM NZ  -->NXM NYM NZM
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)),intent(in) :: divu3
    real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize) :: pp3

    integer :: nvect3,i,j,k,nlock
    integer :: code
    real(mytype) :: tmax,tmoy,tmax1,tmoy1

    nvect3=(ph1%zen(1)-ph1%zst(1)+1)*(ph1%zen(2)-ph1%zst(2)+1)*nzmsize

    if (iibm.eq.0) then
       ta1(:,:,:) = ux1(:,:,:)
       tb1(:,:,:) = uy1(:,:,:)
       tc1(:,:,:) = uz1(:,:,:)
    else
       ta1(:,:,:) = (one - ep1(:,:,:)) * ux1(:,:,:) + ep1(:,:,:)*ubcx
       tb1(:,:,:) = (one - ep1(:,:,:)) * uy1(:,:,:) + ep1(:,:,:)*ubcy
       tc1(:,:,:) = (one - ep1(:,:,:)) * uz1(:,:,:) + ep1(:,:,:)*ubcz
    endif

    !WORK X-PENCILS

    call derxvp(pp1,ta1,di1,sx,cfx6,csx6,cwx6,xsize(1),nxmsize,xsize(2),xsize(3),0)

    if (ilmn.and.(nlock.gt.0)) then
       if ((nlock.eq.1).and.(.not.ivarcoeff)) then
          !! Approximate -div(rho u) using ddt(rho)
          call extrapol_drhodt(ta1, rho1, drho1)
       elseif ((nlock.eq.2).or.ivarcoeff) then
          !! Need to check our error against divu constraint
          !! Or else we are solving the variable-coefficient Poisson equation
          call transpose_z_to_y(-divu3, ta2)
          call transpose_y_to_x(ta2, ta1)
       endif
       call interxvp(pgy1,ta1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)
       pp1(:,:,:) = pp1(:,:,:) + pgy1(:,:,:)
    endif

    call interxvp(pgy1,tb1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)
    call interxvp(pgz1,tc1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)

    call transpose_x_to_y(pp1,duxdxp2,ph4)!->NXM NY NZ
    call transpose_x_to_y(pgy1,uyp2,ph4)
    call transpose_x_to_y(pgz1,uzp2,ph4)

    !WORK Y-PENCILS
    call interyvp(upi2,duxdxp2,dipp2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),1)
    call deryvp(duydypi2,uyp2,dipp2,sy,cfy6,csy6,cwy6,ppyi,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),0)

    !! Compute sum dudx + dvdy
    duydypi2(:,:,:) = duydypi2(:,:,:) + upi2(:,:,:)

    call interyvp(upi2,uzp2,dipp2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),1)

    call transpose_y_to_z(duydypi2,duxydxyp3,ph3)!->NXM NYM NZ
    call transpose_y_to_z(upi2,uzp3,ph3)

    !WORK Z-PENCILS
    call interzvp(pp3,duxydxyp3,dipp3,sz,cifzp6,ciszp6,ciwzp6,(ph1%zen(1)-ph1%zst(1)+1),&
         (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,1)
    call derzvp(po3,uzp3,dipp3,sz,cfz6,csz6,cwz6,(ph1%zen(1)-ph1%zst(1)+1),&
         (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,0)

    !! Compute sum dudx + dvdy + dwdz
    pp3(:,:,:) = pp3(:,:,:) + po3(:,:,:)

    if (nlock==2) then
       pp3(:,:,:)=pp3(:,:,:)-pp3(ph1%zst(1),ph1%zst(2),nzmsize)
    endif

    tmax=-1609._mytype
    tmoy=zero
    do k=1,nzmsize
       do j=ph1%zst(2),ph1%zen(2)
          do i=ph1%zst(1),ph1%zen(1)
             if (pp3(i,j,k).gt.tmax) tmax=pp3(i,j,k)
             tmoy=tmoy+abs(pp3(i,j,k))
          enddo
       enddo
    enddo
    tmoy=tmoy/nvect3

    call MPI_REDUCE(tmax,tmax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
    call MPI_REDUCE(tmoy,tmoy1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

    if ((nrank == 0) .and. (nlock > 0).and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime==ilast)) then
       if (nlock == 2) then
          write(*,*) 'DIV U  max mean=',real(tmax1,mytype),real(tmoy1/real(nproc),mytype)
       else
          write(*,*) 'DIV U* max mean=',real(tmax1,mytype),real(tmoy1/real(nproc),mytype)
       endif
    endif

    return
  end subroutine divergence
  !############################################################################
  !subroutine GRADP
  !Computation of the pressure gradient from the pressure mesh to the
  !velocity mesh
  !Saving pressure gradients on boundaries for correct imposition of
  !BCs on u* via the fractional step methodi (it is not possible to
  !impose BC after correction by pressure gradient otherwise lost of
  !incompressibility--> BCs are imposed on u*
  !
  ! input: pp3 - pressure field (on pressure mesh)
  ! output: px1, py1, pz1 - pressure gradients (on velocity mesh)
  !written by SL 2018
  !############################################################################
  subroutine gradp(px1,py1,pz1,pp3)

    USE param
    USE decomp_2d
    USE variables
    USE MPI
    USE var, only: pp1,pgy1,pgz1,di1,pp2,ppi2,pgy2,pgz2,pgzi2,dip2,&
         pgz3,ppi3,dip3,nxmsize,nymsize,nzmsize

    USE forces, only : iforces, ppi1

    implicit none

    integer :: i,j,k

    real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize) :: pp3
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: px1,py1,pz1

    !WORK Z-PENCILS
    call interzpv(ppi3,pp3,dip3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
         (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
    call derzpv(pgz3,pp3,dip3,sz,cfip6z,csip6z,cwip6z,cfz6,csz6,cwz6,&
         (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)

    !WORK Y-PENCILS
    call transpose_z_to_y(pgz3,pgz2,ph3) !nxm nym nz
    call transpose_z_to_y(ppi3,pp2,ph3)

    call interypv(ppi2,pp2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
         (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
    call derypv(pgy2,pp2,dip2,sy,cfip6y,csip6y,cwip6y,cfy6,csy6,cwy6,ppy,&
         (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
    call interypv(pgzi2,pgz2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
         (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)

    !WORK X-PENCILS

    call transpose_y_to_x(ppi2,pp1,ph2) !nxm ny nz
    call transpose_y_to_x(pgy2,pgy1,ph2)
    call transpose_y_to_x(pgzi2,pgz1,ph2)

    call derxpv(px1,pp1,di1,sx,cfip6,csip6,cwip6,cfx6,csx6,cwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)
    call interxpv(py1,pgy1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)
    call interxpv(pz1,pgz1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)

    if (iforces.eq.1) then
       call interxpv(ppi1,pp1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
            nxmsize,xsize(1),xsize(2),xsize(3),1)
    endif

    !we are in X pencils:
    if (nclx1.eq.2) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             dpdyx1(j,k)=py1(1,j,k)/gdt(itr)
             dpdzx1(j,k)=pz1(1,j,k)/gdt(itr)
          enddo
       enddo
    endif
    if (nclxn.eq.2) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             dpdyxn(j,k)=py1(nx,j,k)/gdt(itr)
             dpdzxn(j,k)=pz1(nx,j,k)/gdt(itr)
          enddo
       enddo
    endif

    if (ncly1.eq.2) then
       if (xstart(2)==1) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                dpdxy1(i,k)=px1(i,1,k)/gdt(itr)
                dpdzy1(i,k)=pz1(i,1,k)/gdt(itr)
             enddo
          enddo
       endif
    endif
    if (nclyn.eq.2) then
       if (xend(2)==ny) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                dpdxyn(i,k)=px1(i,xsize(2),k)/gdt(itr)
                dpdzyn(i,k)=pz1(i,xsize(2),k)/gdt(itr)
             enddo
          enddo
       endif
    endif

    if (nclz1.eq.2) then
       if (xstart(3)==1) then
          do j=1,xsize(2)
             do i=1,xsize(1)
                dpdxz1(i,j)=py1(i,j,1)/gdt(itr)
                dpdyz1(i,j)=pz1(i,j,1)/gdt(itr)
             enddo
          enddo
       endif
    endif
    if (nclzn.eq.2) then
       if (xend(3)==nz) then
          do j=1,xsize(2)
             do i=1,xsize(1)
                dpdxzn(i,j)=py1(i,j,xsize(3))/gdt(itr)
                dpdyzn(i,j)=pz1(i,j,xsize(3))/gdt(itr)
             enddo
          enddo
       endif
    endif

    return
  end subroutine gradp
  !############################################################################
  !############################################################################
  subroutine pre_correc(ux,uy,uz,ep)

    USE decomp_2d
    USE variables
    USE param
    USE var
    USE MPI
    use ibm, only : corgp_ibm, body

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,ep
    integer :: i,j,k,is
    real(mytype) :: ut,ut1,utt,ut11

    integer :: code
    integer, dimension(2) :: dims, dummy_coords
    logical, dimension(2) :: dummy_periods
#ifdef DEBG
    real(mytype) avg_param
#endif

    call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, dummy_periods, dummy_coords, code)

    !********NCLX==2*************************************
    !we are in X pencils:
    if ((itype.eq.itype_channel.or.itype.eq.itype_uniform).and.(nclx1==2.and.nclxn==2)) then

       !Computation of the flow rate Inflow/Outflow
       ut1=zero
       do k=1,xsize(3)
          do j=1,xsize(2)
             ut1=ut1+bxx1(j,k)
          enddo
       enddo
       call MPI_ALLREDUCE(ut1,ut11,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       ut11=ut11/(real(ny*nz,mytype))
       ut=zero
       do k=1,xsize(3)
          do j=1,xsize(2)
             ut=ut+bxxn(j,k)
          enddo
       enddo
       call MPI_ALLREDUCE(ut,utt,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       utt=utt/(real(ny*nz,mytype))
       if ((nrank==0).and.(mod(itime,ilist)==0)) &
          write(*,*) 'Flow rate x I/O/O-I',real(ut11,4),real(utt,4),real(utt-ut11,4)
       do k=1,xsize(3)
          do j=1,xsize(2)
             bxxn(j,k)=bxxn(j,k)-utt+ut11
          enddo
       enddo

    endif

    if (itype.eq.itype_tbl) call tbl_flrt(ux,uy,uz)

    if (nclx1==2) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             dpdyx1(j,k)=dpdyx1(j,k)*gdt(itr)
             dpdzx1(j,k)=dpdzx1(j,k)*gdt(itr)
          enddo
       enddo
       do k=1,xsize(3)
          do j=1,xsize(2)
             ux(1 ,j,k)=bxx1(j,k)
             uy(1 ,j,k)=bxy1(j,k)+dpdyx1(j,k)
             uz(1 ,j,k)=bxz1(j,k)+dpdzx1(j,k)
          enddo
       enddo
    endif
    if (nclxn==2) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             dpdyxn(j,k)=dpdyxn(j,k)*gdt(itr)
             dpdzxn(j,k)=dpdzxn(j,k)*gdt(itr)
          enddo
       enddo
       do k=1,xsize(3)
          do j=1,xsize(2)
             ux(nx,j,k)=bxxn(j,k)
             uy(nx,j,k)=bxyn(j,k)+dpdyxn(j,k)
             uz(nx,j,k)=bxzn(j,k)+dpdzxn(j,k)
          enddo
       enddo
    endif


    !********NCLX==1*************************************
    if (nclx1==1) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             ux(1 ,j,k)=zero
          enddo
       enddo
    endif
    if (nclxn==1) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             ux(nx,j,k)=zero
          enddo
       enddo
    endif

    !********NCLY==2*************************************
    if (ncly1==2) then
       if (xstart(2)==1) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                dpdxy1(i,k)=dpdxy1(i,k)*gdt(itr)
                dpdzy1(i,k)=dpdzy1(i,k)*gdt(itr)
             enddo
          enddo
          do k=1,xsize(3)
             do i=1,xsize(1)
                ux(i,1,k)=byx1(i,k)+dpdxy1(i,k)
                uy(i,1,k)=byy1(i,k)
                uz(i,1,k)=byz1(i,k)+dpdzy1(i,k)
             enddo
          enddo
       endif
    endif

    if (nclyn==2) then
       if (xend(2)==ny) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                dpdxyn(i,k)=dpdxyn(i,k)*gdt(itr)
                dpdzyn(i,k)=dpdzyn(i,k)*gdt(itr)
             enddo
          enddo
       endif
       if (dims(1)==1) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                ux(i,xsize(2),k)=byxn(i,k)+dpdxyn(i,k)
                uy(i,xsize(2),k)=byyn(i,k)
                uz(i,xsize(2),k)=byzn(i,k)+dpdzyn(i,k)
             enddo
          enddo
       elseif (ny - (nym / dims(1)) == xstart(2)) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                ux(i,xsize(2),k)=byxn(i,k)+dpdxyn(i,k)
                uy(i,xsize(2),k)=byyn(i,k)
                uz(i,xsize(2),k)=byzn(i,k)+dpdzyn(i,k)
             enddo
          enddo
       endif
    endif

    !********NCLY==1*************************************
    if (ncly1==1) then
       if (xstart(2)==1) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                uy(i,1,k)=zero
             enddo
          enddo
       endif
    endif

    if (nclyn==1) then
       if (xend(2)==ny) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                uy(i,xsize(2),k)=zero
             enddo
          enddo
       endif
    endif


    !********NCLZ==2*************************************
    if (nclz1==2) then
       if (xstart(3)==1) then
          do j=1,xsize(2)
             do i=1,xsize(1)
                dpdxz1(i,j)=dpdxz1(i,j)*gdt(itr)
                dpdyz1(i,j)=dpdyz1(i,j)*gdt(itr)
             enddo
          enddo
          do j=1,xsize(2)
             do i=1,xsize(1)
                ux(i,j,1)=bzx1(i,j)+dpdxz1(i,j)
                uy(i,j,1)=bzy1(i,j)+dpdyz1(i,j)
                uz(i,j,1)=bzz1(i,j)
             enddo
          enddo
       endif
    endif

    if (nclzn==2) then
       if (xend(3)==nz) then
          do j=1,xsize(2)
             do i=1,xsize(1)
                dpdxzn(i,j)=dpdxzn(i,j)*gdt(itr)
                dpdyzn(i,j)=dpdyzn(i,j)*gdt(itr)
             enddo
          enddo
          do j=1,xsize(2)
             do i=1,xsize(1)
                ux(i,j,xsize(3))=bzxn(i,j)+dpdxzn(i,j)
                uy(i,j,xsize(3))=bzyn(i,j)+dpdyzn(i,j)
                uz(i,j,xsize(3))=bzzn(i,j)
             enddo
          enddo
       endif
    endif
    !********NCLZ==1************************************* !just to reforce free-slip condition
    if (nclz1==1) then
       if (xstart(3)==1) then
          do j=1,xsize(2)
             do i=1,xsize(1)
                uz(i,j,1)=zero
             enddo
          enddo
       endif
    endif

    if (nclzn==1) then
       if (xend(3)==nz) then
          do j=1,xsize(2)
             do i=1,xsize(1)
                uz(i,j,xsize(3))=zero
             enddo
          enddo
       endif
    endif
#ifdef DEBG
    avg_param = zero
    call avg3d (ux, avg_param)
    if (nrank == 0) write(*,*)'## Pres corr ux ', avg_param
    avg_param = zero
    call avg3d (uy, avg_param)
    if (nrank == 0) write(*,*)'## Pres corr uy ', avg_param
    avg_param = zero
    call avg3d (uz, avg_param)
    if (nrank == 0) write(*,*)'## Pres corr uz ', avg_param
#endif

    if (iibm==1) then !solid body old school
       call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,1)
       call body(ux1,uy1,uz1,ep1)
       call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,2)
    endif

    return
  end subroutine pre_correc
  !############################################################################
  !############################################################################
  !! Convert to/from conserved/primary variables
  SUBROUTINE primary_to_conserved(rho1, var1)

    USE decomp_2d, ONLY : mytype, xsize
    USE param, ONLY : nrhotime

    IMPLICIT NONE

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: var1

    var1(:,:,:) = rho1(:,:,:,1) * var1(:,:,:)

  ENDSUBROUTINE primary_to_conserved
  !############################################################################
  !############################################################################
  SUBROUTINE velocity_to_momentum (rho1, ux1, uy1, uz1)

    USE decomp_2d, ONLY : mytype, xsize
    USE param, ONLY : nrhotime
    USE var, ONLY : ilmn

    IMPLICIT NONE

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1

    IF (.NOT.ilmn) THEN
       RETURN
    ENDIF

    CALL primary_to_conserved(rho1, ux1)
    CALL primary_to_conserved(rho1, uy1)
    CALL primary_to_conserved(rho1, uz1)

  ENDSUBROUTINE velocity_to_momentum
  !############################################################################
  !############################################################################
  SUBROUTINE conserved_to_primary(rho1, var1)

    USE decomp_2d, ONLY : mytype, xsize
    USE param, ONLY : nrhotime

    IMPLICIT NONE

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: var1

    var1(:,:,:) = var1(:,:,:) / rho1(:,:,:,1)

  ENDSUBROUTINE conserved_to_primary
  !############################################################################
  !############################################################################
  SUBROUTINE momentum_to_velocity (rho1, ux1, uy1, uz1)

    USE decomp_2d, ONLY : mytype, xsize
    USE param, ONLY : nrhotime
    USE var, ONLY : ilmn

    IMPLICIT NONE

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1

    IF (.NOT.ilmn) THEN
       RETURN
    ENDIF

    CALL conserved_to_primary(rho1, ux1)
    CALL conserved_to_primary(rho1, uy1)
    CALL conserved_to_primary(rho1, uz1)

  ENDSUBROUTINE momentum_to_velocity
  !############################################################################
  !############################################################################
  !! Calculate velocity-divergence constraint
  SUBROUTINE calc_divu_constraint(divu3, rho1, phi1)

    USE decomp_2d, ONLY : mytype, xsize, ysize, zsize
    USE decomp_2d, ONLY : transpose_x_to_y, transpose_y_to_z
    USE param, ONLY : nrhotime, zero, ilmn, pressure0, imultispecies, massfrac, mol_weight
    USE param, ONLY : ibirman_eos
    USE param, ONLY : xnu, prandtl
    USE param, ONLY : one
    USE param, ONLY : iimplicit
    USE variables

    USE var, ONLY : ta1, tb1, tc1, td1, di1
    USE var, ONLY : phi2, ta2, tb2, tc2, td2, te2, di2
    USE var, ONLY : phi3, ta3, tb3, tc3, td3, rho3, di3
    USE param, only : zero
    IMPLICIT NONE

    INTEGER :: is, tmp

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), numscalar) :: phi1
    REAL(mytype), INTENT(OUT), DIMENSION(zsize(1), zsize(2), zsize(3)) :: divu3

    IF (ilmn.and.(.not.ibirman_eos)) THEN
       !!------------------------------------------------------------------------------
       !! X-pencil

       !! We need temperature
       CALL calc_temp_eos(ta1, rho1(:,:,:,1), phi1, tb1, xsize(1), xsize(2), xsize(3))

       CALL derxx (tb1, ta1, di1, sx, sfxp, ssxp, swxp, xsize(1), xsize(2), xsize(3), 1, zero)
       IF (imultispecies) THEN
          tb1(:,:,:) = (xnu / prandtl) * tb1(:,:,:) / ta1(:,:,:)

          !! Calc mean molecular weight
          td1(:,:,:) = zero
          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                td1(:,:,:) = td1(:,:,:) + phi1(:,:,:,is) / mol_weight(is)
             ENDIF
          ENDDO
          td1(:,:,:) = one / td1(:,:,:)

          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                CALL derxx (tc1, phi1(:,:,:,is), di1, sx, sfxp, ssxp, swxp, xsize(1), xsize(2), xsize(3), 1, zero)
                tb1(:,:,:) = tb1(:,:,:) + (xnu / sc(is)) * (td1(:,:,:) / mol_weight(is)) * tc1(:,:,:)
             ENDIF
          ENDDO
       ENDIF

       CALL transpose_x_to_y(ta1, ta2)        !! Temperature
       CALL transpose_x_to_y(tb1, tb2)        !! d2Tdx2
       IF (imultispecies) THEN
          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                CALL transpose_x_to_y(phi1(:,:,:,is), phi2(:,:,:,is))
             ENDIF
          ENDDO
       ENDIF

       !!------------------------------------------------------------------------------
       !! Y-pencil
       tmp = iimplicit
       iimplicit = 0
       CALL deryy (tc2, ta2, di2, sy, sfyp, ssyp, swyp, ysize(1), ysize(2), ysize(3), 1, zero)
       iimplicit = tmp
       IF (imultispecies) THEN
          tc2(:,:,:) = (xnu / prandtl) * tc2(:,:,:) / ta2(:,:,:)

          !! Calc mean molecular weight
          te2(:,:,:) = zero
          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                te2(:,:,:) = te2(:,:,:) + phi2(:,:,:,is) / mol_weight(is)
             ENDIF
          ENDDO
          te2(:,:,:) = one / te2(:,:,:)

          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                tmp = iimplicit
                iimplicit = 0
                CALL deryy (td2, phi2(:,:,:,is), di2, sy, sfyp, ssyp, swyp, ysize(1), ysize(2), ysize(3), 1, zero)
                iimplicit = tmp
                tc2(:,:,:) = tc2(:,:,:) + (xnu / sc(is)) * (te2(:,:,:) / mol_weight(is)) * td2(:,:,:)
             ENDIF
          ENDDO
       ENDIF
       tb2(:,:,:) = tb2(:,:,:) + tc2(:,:,:)

       CALL transpose_y_to_z(ta2, ta3)        !! Temperature
       CALL transpose_y_to_z(tb2, tb3)        !! d2Tdx2 + d2Tdy2
       IF (imultispecies) THEN
          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                CALL transpose_y_to_z(phi2(:,:,:,is), phi3(:,:,:,is))
             ENDIF
          ENDDO
       ENDIF

       !!------------------------------------------------------------------------------
       !! Z-pencil
       CALL derzz (divu3, ta3, di3, sz, sfzp, sszp, swzp, zsize(1), zsize(2), zsize(3), 1, zero)
       IF (imultispecies) THEN
          divu3(:,:,:) = (xnu / prandtl) * divu3(:,:,:) / ta3(:,:,:)

          !! Calc mean molecular weight
          td3(:,:,:) = zero
          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                td3(:,:,:) = td3(:,:,:) + phi3(:,:,:,is) / mol_weight(is)
             ENDIF
          ENDDO
          td3(:,:,:) = one / td3(:,:,:)

          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                CALL derzz (tc3, phi3(:,:,:,is), di3, sz, sfzp, sszp, swzp, zsize(1), zsize(2), zsize(3), 1, zero)
                divu3(:,:,:) = divu3(:,:,:) + (xnu / sc(is)) * (td3(:,:,:) / mol_weight(is)) * tc3(:,:,:)
             ENDIF
          ENDDO
       ENDIF
       divu3(:,:,:) = divu3(:,:,:) + tb3(:,:,:)

       IF (imultispecies) THEN
          !! Thus far we have computed rho * divu, want divu
          CALL calc_rho_eos(rho3, ta3, phi3, tb3, zsize(1), zsize(2), zsize(3))
          divu3(:,:,:) = divu3(:,:,:) / rho3(:,:,:)
       ELSE
          divu3(:,:,:) = (xnu / prandtl) * divu3(:,:,:) / pressure0
       ENDIF
    ELSE
       divu3(:,:,:) = zero
    ENDIF

  ENDSUBROUTINE calc_divu_constraint

  !############################################################################
  !############################################################################
  ! Calculate extrapolation drhodt 
  SUBROUTINE extrapol_drhodt(drhodt1_next, rho1, drho1)

    USE decomp_2d, ONLY : mytype, xsize, nrank
    USE param, ONLY : ntime, nrhotime, itime, itimescheme, itr, dt, gdt, irestart
    USE param, ONLY : half, three, four
    USE param, ONLY : ibirman_eos

    IMPLICIT NONE

    INTEGER :: subitr

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), ntime) :: drho1
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), INTENT(OUT), DIMENSION(xsize(1), xsize(2), xsize(3)) :: drhodt1_next

    IF (itimescheme.EQ.1) THEN
       !! EULER
       drhodt1_next(:,:,:) = (rho1(:,:,:,1) - rho1(:,:,:,2)) / dt
    ELSEIF (itimescheme.EQ.2) THEN
       !! AB2
       IF ((itime.EQ.1).AND.(irestart.EQ.0)) THEN
          drhodt1_next(:,:,:) = (rho1(:,:,:,1) - rho1(:,:,:,2)) / dt
       ELSE
          drhodt1_next(:,:,:) = three * rho1(:,:,:,1) - four * rho1(:,:,:,2) + rho1(:,:,:,3)
          drhodt1_next(:,:,:) = half * drhodt1_next(:,:,:) / dt
       ENDIF
       ! ELSEIF (itimescheme.EQ.3) THEN
       !    !! AB3
       ! ELSEIF (itimescheme.EQ.4) THEN
       !    !! AB4
    ELSEIF (itimescheme.EQ.5) THEN
       !! RK3
       IF (itime.GT.1) THEN
          drhodt1_next(:,:,:) = rho1(:,:,:,2)
          DO subitr = 1, itr
             drhodt1_next(:,:,:) = drhodt1_next(:,:,:) + (gdt(subitr) / dt) &
                  * (rho1(:,:,:,2) - rho1(:,:,:,3))
          ENDDO
       ELSE
          drhodt1_next(:,:,:) = drho1(:,:,:,1)
       ENDIF
    ELSE
       IF (nrank.EQ.0) THEN
          PRINT *, "Extrapolating drhodt not implemented for timescheme:", itimescheme
          STOP
       ENDIF
    ENDIF

    IF (ibirman_eos) THEN
       CALL birman_drhodt_corr(drhodt1_next, rho1)
    ENDIF

  ENDSUBROUTINE extrapol_drhodt
  !############################################################################
  !!
  !! Subroutine : birman_drhodt_corr
  !! Author     :
  !! Description: Calculate extrapolation drhodt correction
  !!
  !############################################################################

  SUBROUTINE birman_drhodt_corr(drhodt1_next, rho1)

    USE decomp_2d, ONLY : mytype, xsize, ysize, zsize
    USE decomp_2d, ONLY : transpose_x_to_y, transpose_y_to_z, transpose_z_to_y, transpose_y_to_x
    USE variables, ONLY : derxx, deryy, derzz
    USE param, ONLY : nrhotime
    USE param, ONLY : xnu, prandtl
    USE param, ONLY : iimplicit

    USE var, ONLY : td1, te1, di1, sx, sfxp, ssxp, swxp
    USE var, ONLY : rho2, ta2, tb2, di2, sy, sfyp, ssyp, swyp
    USE var, ONLY : rho3, ta3, di3, sz, sfzp, sszp, swzp
    USE param, only : zero
    IMPLICIT NONE

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: drhodt1_next

    REAL(mytype) :: invpe

    invpe = xnu / prandtl

    CALL transpose_x_to_y(rho1(:,:,:,1), rho2)
    CALL transpose_y_to_z(rho2, rho3)

    !! Diffusion term
    CALL derzz (ta3,rho3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1, zero)
    CALL transpose_z_to_y(ta3, tb2)

    iimplicit = -iimplicit
    CALL deryy (ta2,rho2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1, zero)
    iimplicit = -iimplicit
    ta2(:,:,:) = ta2(:,:,:) + tb2(:,:,:)
    CALL transpose_y_to_x(ta2, te1)

    CALL derxx (td1,rho1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1, zero)
    td1(:,:,:) = td1(:,:,:) + te1(:,:,:)

    drhodt1_next(:,:,:) = drhodt1_next(:,:,:) - invpe * td1(:,:,:)

  ENDSUBROUTINE birman_drhodt_corr
  !############################################################################
  !!
  !!  SUBROUTINE: test_varcoeff
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Tests convergence of the variable-coefficient Poisson solver
  !!
  !############################################################################
  SUBROUTINE test_varcoeff(converged, divup3norm, pp3, dv3, atol, rtol, poissiter)

    USE MPI
    USE decomp_2d, ONLY: mytype, ph1, real_type, nrank
    USE var, ONLY : nzmsize
    USE param, ONLY : npress, itime
    USE variables, ONLY : nxm, nym, nzm, ilist

    IMPLICIT NONE

    !! INPUTS
    REAL(mytype), INTENT(INOUT), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
    REAL(mytype), INTENT(IN), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize) :: dv3
    REAL(mytype), INTENT(IN) :: atol, rtol
    INTEGER, INTENT(IN) :: poissiter

    !! OUTPUTS
    LOGICAL, INTENT(OUT) :: converged
    REAL(mytype) :: divup3norm

    !! LOCALS
    INTEGER :: ierr
    REAL(mytype) :: errloc, errglob

    IF (poissiter.EQ.0) THEN
       errloc = SUM(dv3**2)
       CALL MPI_ALLREDUCE(errloc,divup3norm,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
       divup3norm = SQRT(divup3norm / nxm / nym / nzm)

       if (nrank.eq.0.and.mod(itime, ilist) == 0) then
          write(*,*)  "solving variable-coefficient poisson equation:"
          write(*,*)  "+ rms div(u*) - div(u): ", divup3norm
       endif
    else
       !! Compute RMS change
       errloc = SUM((pp3(:,:,:,1) - pp3(:,:,:,2))**2)
       CALL MPI_ALLREDUCE(errloc,errglob,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
       errglob = SQRT(errglob / nxm / nym / nzm)

       if (nrank.eq.0.and.mod(itime, ilist) == 0) then
          write(*,*)  "+ RMS change in pressure: ", errglob
       endif

       if (errglob.le.atol) then
          converged = .true.
          if (nrank.eq.0.and.mod(itime, ilist) == 0) then
             write(*,*)  "- Converged: atol"
          endif
       endif

       !! Compare RMS change to size of |div(u*) - div(u)|
       if (errglob.lt.(rtol * divup3norm)) then
          converged = .true.
          if (nrank.eq.0.and.mod(itime, ilist) == 0) then
             write(*,*)  "- Converged: rtol"
          endif
       endif

       if (.not.converged) then
          pp3(:,:,:,2) = pp3(:,:,:,1)
       endif
    endif

  ENDSUBROUTINE test_varcoeff
  !############################################################################
  !!
  !!  SUBROUTINE: calc_varcoeff_rhs
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Computes RHS of the variable-coefficient Poisson solver
  !!
  !############################################################################
  SUBROUTINE calc_varcoeff_rhs(pp3, rho1, px1, py1, pz1, dv3, drho1, ep1, divu3, rho0, poissiter)

    USE MPI

    USE decomp_2d

    USE param, ONLY : nrhotime, ntime
    USE param, ONLY : one

    USE var, ONLY : ta1, tb1, tc1
    USE var, ONLY : nzmsize

    IMPLICIT NONE

    !! INPUTS
    INTEGER, INTENT(IN) :: poissiter
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3)) :: px1, py1, pz1
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), ntime) :: drho1
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ep1
    REAL(mytype), INTENT(IN), DIMENSION(zsize(1), zsize(2), zsize(3)) :: divu3
    REAL(mytype), INTENT(IN), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize) :: dv3
    real(mytype) :: rho0

    !! OUTPUTS
    REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize) :: pp3

    !! LOCALS
    INTEGER :: nlock, ierr
    REAL(mytype) :: rhomin

    IF (poissiter.EQ.0) THEN
       !! Compute rho0
       rhomin = MINVAL(rho1(:,:,:,1))

       CALL MPI_ALLREDUCE(rhomin,rho0,1,real_type,MPI_MIN,MPI_COMM_WORLD,ierr)
    ENDIF

    ta1(:,:,:) = (one - rho0 / rho1(:,:,:,1)) * px1(:,:,:)
    tb1(:,:,:) = (one - rho0 / rho1(:,:,:,1)) * py1(:,:,:)
    tc1(:,:,:) = (one - rho0 / rho1(:,:,:,1)) * pz1(:,:,:)

    nlock = -1 !! Don't do any funny business with LMN
    CALL divergence(pp3,rho1,ta1,tb1,tc1,ep1,drho1,divu3,nlock)

    !! lapl(p) = div((1 - rho0/rho) grad(p)) + rho0(div(u*) - div(u))
    !! dv3 contains div(u*) - div(u)
    pp3(:,:,:) = pp3(:,:,:) + rho0 * dv3(:,:,:)

  ENDSUBROUTINE calc_varcoeff_rhs
  !############################################################################
  !********************************************************************
  !
  subroutine tbl_flrt (ux1,uy1,uz1)
  !
  !********************************************************************

    USE decomp_2d
    USE decomp_2d_poisson
    USE variables
    USE param
    USE MPI

    implicit none
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2

    integer :: j,i,k,code
    real(mytype) :: can,ut1,ut2,ut3,ut4,utt1,utt2,utt3,utt4,udif

    ux1(1,:,:)=bxx1(:,:)
    ux1(nx,:,:)=bxxn(:,:)

    call transpose_x_to_y(ux1,ux2)
    call transpose_x_to_y(uy1,uy2)
    ! Flow rate at the inlet
    ut1=zero;utt1=zero
    if (ystart(1)==1) then !! CPUs at the inlet
      do k=1,ysize(3)
        do j=1,ysize(2)-1
          ut1=ut1+(yp(j+1)-yp(j))*(ux2(1,j+1,k)-half*(ux2(1,j+1,k)-ux2(1,j,k)))
        enddo
      enddo
      ! ut1=ut1/real(ysize(3),mytype)
    endif
    call MPI_ALLREDUCE(ut1,utt1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    utt1=utt1/real(nz,mytype) !! Volume flow rate per unit spanwise dist
    ! Flow rate at the outlet
    ut2=zero;utt2=zero
    if (yend(1)==nx) then !! CPUs at the outlet
      do k=1,ysize(3)
        do j=1,ysize(2)-1
          ut2=ut2+(yp(j+1)-yp(j))*(ux2(ysize(1),j+1,k)-half*(ux2(ysize(1),j+1,k)-ux2(ysize(1),j,k)))
        enddo
      enddo
      ! ut2=ut2/real(ysize(3),mytype)
    endif

    call MPI_ALLREDUCE(ut2,utt2,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    utt2=utt2/real(nz,mytype) !! Volume flow rate per unit spanwise dist

    ! Flow rate at the top and bottom
    ut3=zero
    ut4=zero
    do k=1,ysize(3)
      do i=1,ysize(1)
        ut3=ut3+uy2(i,1,k)
        ut4=ut4+uy2(i,ny,k)
      enddo
    enddo
    call MPI_ALLREDUCE(ut3,utt3,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(ut4,utt4,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    utt3=utt3/(real(nx*nz,mytype))*xlx  !!! Volume flow rate per unit spanwise dist
    utt4=utt4/(real(nx*nz,mytype))*xlx  !!! Volume flow rate per unit spanwise dist

    !! velocity correction
    udif=(utt1-utt2+utt3-utt4)/yly
    if ((nrank==0).and.(mod(itime,ilist)==0)) then
      write(*,"(' Mass balance: L-BC, R-BC,',2f12.6)") utt1,utt2
      write(*,"(' Mass balance: B-BC, T-BC, Crr-Vel',3f11.5)") utt3,utt4,udif
    endif
    ! do k=1,xsize(3)
    !   do j=1,xsize(2)
    !     ux1(nx,i,k)=ux1(nx,i,k)+udif
    !   enddo
    ! enddo
    do k=1,xsize(3)
      do j=1,xsize(2)
        bxxn(j,k)=bxxn(j,k)+udif
      enddo
    enddo

  end subroutine tbl_flrt
!############################################################################
!!
!!  SUBROUTINE: avg3d
!!      AUTHOR: Stefano Rolfo
!! DESCRIPTION: Compute the total sum of a a 3d field
!!
!############################################################################
subroutine avg3d (var, avg)

  use decomp_2d, only: real_type, xsize, xend
  use param
  use dbg_schemes, only: sqrt_prec
  use variables, only: nx,ny,nz,nxm,nym,nzm
  use mpi

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: var
  real(mytype), intent(out) :: avg
  real(mytype)              :: dep

  integer :: i,j,k, code
  integer :: nxc, nyc, nzc, xsize1, xsize2, xsize3

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

  dep=zero
  do k=1,xsize3
     do j=1,xsize2
        do i=1,xsize1
           !dep=dep+var(i,j,k)**2
           dep=dep+var(i,j,k)
        enddo
     enddo
  enddo
  call MPI_ALLREDUCE(dep,avg,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
  !avg=sqrt_prec(avg)/(nxc*nyc*nzc)
  avg=avg/(nxc*nyc*nzc)

  return

end subroutine avg3d

endmodule navier
