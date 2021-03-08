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
  public :: gradp

contains
  !############################################################################
  !!  SUBROUTINE: solve_poisson
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Takes the intermediate momentum field as input,
  !!              computes div and solves pressure-Poisson equation.
  !############################################################################
  SUBROUTINE solve_poisson(pp3, px1, py1, pz1, rho1, ux1, uy1, uz1, ep1, drho1, divu3)

    use decomp_2d, only : mytype, xsize, zsize, ph1
    use decomp_2d_poisson, only : poisson
    use var, only : nzmsize
    use var, only : dv3
    use param, only : ntime, nrhotime, npress
    use param, only : ilmn, ivarcoeff, one


    implicit none

    !! Inputs
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ux1, uy1, uz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ep1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), nrhotime), intent(in) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime), intent(in) :: drho1
    real(mytype), dimension(zsize(1), zsize(2), zsize(3)), intent(in) :: divu3

    !! Outputs
    real(mytype), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: px1, py1, pz1

    !! Locals
    integer :: nlock, poissiter
    logical :: converged
    real(mytype) :: atol, rtol, rho0

    nlock = 1 !! Corresponds to computing div(u*)
    converged = .false.
    poissiter = 0
    rho0 = one

    atol = 1.0e-14_mytype !! Absolute tolerance for Poisson solver
    rtol = 1.0e-14_mytype !! Relative tolerance for Poisson solver

    if (ilmn.and.ivarcoeff) then
       !! variable-coefficient poisson solver works on div(u), not div(rho u)
       !! rho u -> u
       call momentum_to_velocity(rho1, ux1, uy1, uz1)
    endif

    call divergence(pp3(:,:,:,1),rho1,ux1,uy1,uz1,ep1,drho1,divu3,nlock)
    if (ilmn.and.ivarcoeff) then
       dv3(:,:,:) = pp3(:,:,:,1)
    endif

    do while(.not.converged)
       if (ivarcoeff) then

          !! test convergence
          call test_varcoeff(converged, pp3, dv3, atol, rtol, poissiter)

          if (.not.converged) then
             !! evaluate additional rhs terms
             call calc_varcoeff_rhs(pp3(:,:,:,1), rho1, px1, py1, pz1, dv3, drho1, ep1, divu3, rho0, &
                  poissiter)
          endif
       endif

       if (.not.converged) then
          call poisson(pp3(:,:,:,1))

          !! need to update pressure gradient here for varcoeff
          call gradp(px1,py1,pz1,pp3(:,:,:,1))

          if ((.not.ilmn).or.(.not.ivarcoeff)) then
             !! once-through solver
             !! - incompressible flow
             !! - lmn - constant-coefficient solver
             converged = .true.
          endif
       endif

       poissiter = poissiter + 1
    enddo

    if (ilmn.and.ivarcoeff) then
       !! variable-coefficient poisson solver works on div(u), not div(rho u)
       !! u -> rho u
       call velocity_to_momentum(rho1, ux1, uy1, uz1)
    endif

  end subroutine solve_poisson
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

    use decomp_2d
    use param, only : zero
    use param, only : imultispecies, massfrac, mol_weight
    use param, only : ntime
    use var, only : numscalar
    use var, only : ta1, tb1

    implicit none

    !! INPUTS
    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: dtemp1, rho1
    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1
    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), ntime, numscalar) :: dphi1

    !! OUTPUTS
    real(mytype), intent(out), dimension(xsize(1), xsize(2), xsize(3)) :: drho1

    !! LOCALS
    integer :: is

    drho1(:,:,:) = zero

    if (imultispecies) then
       do is = 1, numscalar
          if (massfrac(is)) then
             drho1(:,:,:) = drho1(:,:,:) - dphi1(:,:,:,1,is) / mol_weight(is)
          endif
       enddo

       ta1(:,:,:) = zero !! Mean molecular weight
       do is = 1, numscalar
          if (massfrac(is)) then
             ta1(:,:,:) = ta1(:,:,:) + phi1(:,:,:,is) / mol_weight(is)
          endif
       enddo
       drho1(:,:,:) = drho1(:,:,:) / ta1(:,:,:)  !! XXX ta1 is the inverse molecular weight
    endif

    call calc_temp_eos(ta1, rho1, phi1, tb1, xsize(1), xsize(2), xsize(3))
    drho1(:,:,:) = drho1(:,:,:) - dtemp1(:,:,:) / ta1(:,:,:)

    drho1(:,:,:) = rho1(:,:,:) * drho1(:,:,:)

  endsubroutine lmn_t_to_rho_trans
  !############################################################################
  !subroutine COR_VEL
  !Correction of u* by the pressure gradient to get a divergence free
  !field
  ! input : px,py,pz
  ! output : ux,uy,uz
  !written by SL 2018
  !############################################################################
  subroutine cor_vel (ux,uy,uz,px,py,pz)

    use decomp_2d
    use variables
    use param

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: px,py,pz

    ux(:,:,:)=ux(:,:,:)-px(:,:,:)
    uy(:,:,:)=uy(:,:,:)-py(:,:,:)
    uz(:,:,:)=uz(:,:,:)-pz(:,:,:)

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

    use param
    use decomp_2d
    use variables
    use var, only: ta1, tb1, tc1, pp1, pgy1, pgz1, di1, &
         duxdxp2, uyp2, uzp2, duydypi2, upi2, ta2, dipp2, &
         duxydxyp3, uzp3, po3, dipp3, nxmsize, nymsize, nzmsize
    use mpi

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

    if (iibm == 0) then
       ta1(:,:,:) = ux1(:,:,:)
       tb1(:,:,:) = uy1(:,:,:)
       tc1(:,:,:) = uz1(:,:,:)
    else
       ta1(:,:,:) = (one - ep1(:,:,:)) * ux1(:,:,:)
       tb1(:,:,:) = (one - ep1(:,:,:)) * uy1(:,:,:)
       tc1(:,:,:) = (one - ep1(:,:,:)) * uz1(:,:,:)
    endif

    !WORK X-PENCILS

    call derxvp(pp1,ta1,di1,sx,cfx6,csx6,cwx6,xsize(1),nxmsize,xsize(2),xsize(3),0)

    if (ilmn.and.(nlock > 0)) then
       if ((nlock == 1).and.(.not.ivarcoeff)) then
          !! Approximate -div(rho u) using ddt(rho)
          call extrapol_drhodt(ta1, rho1, drho1)
       elseif ((nlock == 2).or.ivarcoeff) then
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

    if (nlock == 2) then
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
    USE TBL, ONLY:tbl_flrt
    use ibm, only : corgp_ibm, body

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,ep
    integer :: i,j,k,is
    real(mytype) :: ut,ut1,utt,ut11

    integer :: code
    integer, dimension(2) :: dims, dummy_coords
    logical, dimension(2) :: dummy_periods

    call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, dummy_periods, dummy_coords, code)

    !********NCLX==2*************************************
    !we are in X pencils:
    if ((itype.eq.itype_channel).and.(nclx1==2.and.nclxn==2)) then

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
       if (nrank==0) write(*,*) 'Flow rate x I/O/O-I',real(ut11,4),real(utt,4),real(utt-ut11,4)
       do k=1,xsize(3)
          do j=1,xsize(2)
             bxxn(j,k)=bxxn(j,k)-utt+ut11
          enddo
       enddo

    endif

    if (itype.eq.itype_tbl) call tbl_flrt(ux,uy,uz)

    if (nclx1 == 2) then
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
    if (nclxn == 2) then
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
    if (nclx1 == 1) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             ux(1 ,j,k)=zero
          enddo
       enddo
    endif
    if (nclxn == 1) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             ux(nx,j,k)=zero
          enddo
       enddo
    endif

    !********NCLY==2*************************************
    if (ncly1 == 2) then
       if (xstart(2) == 1) then
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

    if (nclyn == 2) then
       if (xend(2) == ny) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                dpdxyn(i,k)=dpdxyn(i,k)*gdt(itr)
                dpdzyn(i,k)=dpdzyn(i,k)*gdt(itr)
             enddo
          enddo
       endif
       if (dims(1) == 1) then
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
    if (ncly1 == 1) then
       if (xstart(2) == 1) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                uy(i,1,k)=zero
             enddo
          enddo
       endif
    endif

    if (nclyn == 1) then
       if (xend(2) == ny) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                uy(i,xsize(2),k)=zero
             enddo
          enddo
       endif
    endif


    !********NCLZ==2*************************************
    if (nclz1 == 2) then
       if (xstart(3) == 1) then
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

    if (nclzn == 2) then
       if (xend(3) == nz) then
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
    if (nclz1 == 1) then
       if (xstart(3) == 1) then
          do j=1,xsize(2)
             do i=1,xsize(1)
                uz(i,j,1)=zero
             enddo
          enddo
       endif
    endif

    if (nclzn == 1) then
       if (xend(3) == nz) then
          do j=1,xsize(2)
             do i=1,xsize(1)
                uz(i,j,xsize(3))=zero
             enddo
          enddo
       endif
    endif

    if (iibm == 1) then !solid body old school
       call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,1)
       call body(ux1,uy1,uz1,ep1)
       call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,2)
    endif

    return
  end subroutine pre_correc
  !############################################################################
  !############################################################################
  !! Convert to/from conserved/primary variables
  subroutine primary_to_conserved(rho1, var1)

    use decomp_2d, only : mytype, xsize
    use param, only : nrhotime

    implicit none

    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: var1

    var1(:,:,:) = rho1(:,:,:,1) * var1(:,:,:)

  endsubroutine primary_to_conserved
  !############################################################################
  !############################################################################
  subroutine velocity_to_momentum (rho1, ux1, uy1, uz1)

    use decomp_2d, only : mytype, xsize
    use param, only : nrhotime
    use var, only : ilmn

    implicit none

    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1

    if (.not.ilmn) then
       return
    endif

    call primary_to_conserved(rho1, ux1)
    call primary_to_conserved(rho1, uy1)
    call primary_to_conserved(rho1, uz1)

  endsubroutine velocity_to_momentum
  !############################################################################
  !############################################################################
  subroutine conserved_to_primary(rho1, var1)

    use decomp_2d, only : mytype, xsize
    use param, only : nrhotime

    implicit none

    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: var1

    var1(:,:,:) = var1(:,:,:) / rho1(:,:,:,1)

  endsubroutine conserved_to_primary
  !############################################################################
  !############################################################################
  subroutine momentum_to_velocity (rho1, ux1, uy1, uz1)

    use decomp_2d, only : mytype, xsize
    use param, only : nrhotime
    use var, only : ilmn

    implicit none

    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1

    if (.not.ilmn) then
       return
    endif

    call conserved_to_primary(rho1, ux1)
    call conserved_to_primary(rho1, uy1)
    call conserved_to_primary(rho1, uz1)

  endsubroutine momentum_to_velocity
  !############################################################################
  !############################################################################
  !! Calculate velocity-divergence constraint
  subroutine calc_divu_constraint(divu3, rho1, phi1)

    use decomp_2d, only : mytype, xsize, ysize, zsize
    use decomp_2d, only : transpose_x_to_y, transpose_y_to_z
    use param, only : nrhotime, zero, ilmn, pressure0, imultispecies, massfrac, mol_weight
    use param, only : ibirman_eos
    use param, only : xnu, prandtl
    use param, only : one
    use param, only : iimplicit
    use variables

    use var, only : ta1, tb1, tc1, td1, di1
    use var, only : phi2, ta2, tb2, tc2, td2, te2, di2
    use var, only : phi3, ta3, tb3, tc3, td3, rho3, di3

    implicit none

    integer :: is

    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1
    real(mytype), intent(out), dimension(zsize(1), zsize(2), zsize(3)) :: divu3

    if (ilmn.and.(.not.ibirman_eos)) then
       !!------------------------------------------------------------------------------
       !! X-pencil

       !! We need temperature
       call calc_temp_eos(ta1, rho1(:,:,:,1), phi1, tb1, xsize(1), xsize(2), xsize(3))

       call derxx (tb1, ta1, di1, sx, sfxp, ssxp, swxp, xsize(1), xsize(2), xsize(3), 1)
       if (imultispecies) then
          tb1(:,:,:) = (xnu / prandtl) * tb1(:,:,:) / ta1(:,:,:)

          !! Calc mean molecular weight
          td1(:,:,:) = zero
          do is = 1, numscalar
             if (massfrac(is)) then
                td1(:,:,:) = td1(:,:,:) + phi1(:,:,:,is) / mol_weight(is)
             endif
          enddo
          td1(:,:,:) = one / td1(:,:,:)

          do is = 1, numscalar
             if (massfrac(is)) then
                call derxx (tc1, phi1(:,:,:,is), di1, sx, sfxp, ssxp, swxp, xsize(1), xsize(2), xsize(3), 1)
                tb1(:,:,:) = tb1(:,:,:) + (xnu / sc(is)) * (td1(:,:,:) / mol_weight(is)) * tc1(:,:,:)
             endif
          enddo
       endif

       call transpose_x_to_y(ta1, ta2)        !! Temperature
       call transpose_x_to_y(tb1, tb2)        !! d2Tdx2
       if (imultispecies) then
          do is = 1, numscalar
             if (massfrac(is)) then
                call transpose_x_to_y(phi1(:,:,:,is), phi2(:,:,:,is))
             endif
          enddo
       endif

       !!------------------------------------------------------------------------------
       !! Y-pencil
       iimplicit = -iimplicit
       call deryy (tc2, ta2, di2, sy, sfyp, ssyp, swyp, ysize(1), ysize(2), ysize(3), 1)
       if (imultispecies) then
          tc2(:,:,:) = (xnu / prandtl) * tc2(:,:,:) / ta2(:,:,:)

          !! Calc mean molecular weight
          te2(:,:,:) = zero
          do is = 1, numscalar
             if (massfrac(is)) then
                te2(:,:,:) = te2(:,:,:) + phi2(:,:,:,is) / mol_weight(is)
             endif
          enddo
          te2(:,:,:) = one / te2(:,:,:)

          do is = 1, numscalar
             if (massfrac(is)) then
                iimplicit = -iimplicit
                call deryy (td2, phi2(:,:,:,is), di2, sy, sfyp, ssyp, swyp, ysize(1), ysize(2), ysize(3), 1)
                iimplicit = -iimplicit
                tc2(:,:,:) = tc2(:,:,:) + (xnu / sc(is)) * (te2(:,:,:) / mol_weight(is)) * td2(:,:,:)
             endif
          enddo
       endif
       tb2(:,:,:) = tb2(:,:,:) + tc2(:,:,:)

       call transpose_y_to_z(ta2, ta3)        !! Temperature
       call transpose_y_to_z(tb2, tb3)        !! d2Tdx2 + d2Tdy2
       if (imultispecies) then
          do is = 1, numscalar
             if (massfrac(is)) then
                call transpose_y_to_z(phi2(:,:,:,is), phi3(:,:,:,is))
             endif
          enddo
       endif

       !!------------------------------------------------------------------------------
       !! Z-pencil
       call derzz (divu3, ta3, di3, sz, sfzp, sszp, swzp, zsize(1), zsize(2), zsize(3), 1)
       if (imultispecies) then
          divu3(:,:,:) = (xnu / prandtl) * divu3(:,:,:) / ta3(:,:,:)

          !! Calc mean molecular weight
          td3(:,:,:) = zero
          do is = 1, numscalar
             if (massfrac(is)) then
                td3(:,:,:) = td3(:,:,:) + phi3(:,:,:,is) / mol_weight(is)
             endif
          enddo
          td3(:,:,:) = one / td3(:,:,:)

          do is = 1, numscalar
             if (massfrac(is)) then
                call derzz (tc3, phi3(:,:,:,is), di3, sz, sfzp, sszp, swzp, zsize(1), zsize(2), zsize(3), 1)
                divu3(:,:,:) = divu3(:,:,:) + (xnu / sc(is)) * (td3(:,:,:) / mol_weight(is)) * tc3(:,:,:)
             endif
          enddo
       endif
       divu3(:,:,:) = divu3(:,:,:) + tb3(:,:,:)

       if (imultispecies) then
          !! thus far we have computed rho * divu, want divu
          call calc_rho_eos(rho3, ta3, phi3, tb3, zsize(1), zsize(2), zsize(3))
          divu3(:,:,:) = divu3(:,:,:) / rho3(:,:,:)
       else
          divu3(:,:,:) = (xnu / prandtl) * divu3(:,:,:) / pressure0
       endif
    else
       divu3(:,:,:) = zero
    endif

  endsubroutine calc_divu_constraint
  !############################################################################
  !############################################################################
  ! Calculate extrapolation drhodt 
  subroutine extrapol_drhodt(drhodt1_next, rho1, drho1)

    use decomp_2d, only : mytype, xsize, nrank
    use param, only : ntime, nrhotime, itime, itimescheme, itr, dt, gdt, irestart
    use param, only : half, three, four
    use param, only : ibirman_eos

    implicit none

    integer :: subitr

    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), ntime) :: drho1
    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    real(mytype), intent(out), dimension(xsize(1), xsize(2), xsize(3)) :: drhodt1_next

    if (itimescheme.eq.1) then
       !! EULER
       drhodt1_next(:,:,:) = (rho1(:,:,:,1) - rho1(:,:,:,2)) / dt
    elseif (itimescheme.eq.2) then
       !! AB2
       if ((itime.eq.1).and.(irestart.eq.0)) then
          drhodt1_next(:,:,:) = (rho1(:,:,:,1) - rho1(:,:,:,2)) / dt
       else
          drhodt1_next(:,:,:) = three * rho1(:,:,:,1) - four * rho1(:,:,:,2) + rho1(:,:,:,3)
          drhodt1_next(:,:,:) = half * drhodt1_next(:,:,:) / dt
       endif
       ! ELSEIF (itimescheme.EQ.3) THEN
       !    !! AB3
       ! ELSEIF (itimescheme.EQ.4) THEN
       !    !! AB4
    elseif (itimescheme.eq.5) then
       !! RK3
       if (itime.gt.1) then
          drhodt1_next(:,:,:) = rho1(:,:,:,2)
          do subitr = 1, itr
             drhodt1_next(:,:,:) = drhodt1_next(:,:,:) + (gdt(subitr) / dt) &
                  * (rho1(:,:,:,2) - rho1(:,:,:,3))
          enddo
       else
          drhodt1_next(:,:,:) = drho1(:,:,:,1)
       endif
    else
       if (nrank.eq.0) then
          write(*,*)  "extrapolating drhodt not implemented for timescheme:", itimescheme
          stop
       endif
    endif

    if (ibirman_eos) then
       call birman_drhodt_corr(drhodt1_next, rho1)
    endif

  endsubroutine extrapol_drhodt
  !############################################################################
  !!
  !! Subroutine : birman_drhodt_corr
  !! Author     :
  !! Description: Calculate extrapolation drhodt correction
  !!
  !############################################################################
  subroutine birman_drhodt_corr(drhodt1_next, rho1)

    use decomp_2d, only : mytype, xsize, ysize, zsize
    use decomp_2d, only : transpose_x_to_y, transpose_y_to_z, transpose_z_to_y, transpose_y_to_x
    use variables, only : derxx, deryy, derzz
    use param, only : nrhotime
    use param, only : xnu, prandtl
    use param, only : iimplicit

    use var, only : td1, te1, di1, sx, sfxp, ssxp, swxp
    use var, only : rho2, ta2, tb2, di2, sy, sfyp, ssyp, swyp
    use var, only : rho3, ta3, di3, sz, sfzp, sszp, swzp

    implicit none

    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: drhodt1_next

    real(mytype) :: invpe

    invpe = xnu / prandtl

    call transpose_x_to_y(rho1(:,:,:,1), rho2)
    call transpose_y_to_z(rho2, rho3)

    !! Diffusion term
    call derzz (ta3,rho3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
    call transpose_z_to_y(ta3, tb2)

    iimplicit = -iimplicit
    call deryy (ta2,rho2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
    iimplicit = -iimplicit
    ta2(:,:,:) = ta2(:,:,:) + tb2(:,:,:)
    call transpose_y_to_x(ta2, te1)

    call derxx (td1,rho1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
    td1(:,:,:) = td1(:,:,:) + te1(:,:,:)

    drhodt1_next(:,:,:) = drhodt1_next(:,:,:) - invpe * td1(:,:,:)

  endsubroutine birman_drhodt_corr
  !############################################################################
  !!
  !!  SUBROUTINE: test_varcoeff
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Tests convergence of the variable-coefficient Poisson solver
  !!
  !############################################################################
  subroutine test_varcoeff(converged, pp3, dv3, atol, rtol, poissiter)

    use mpi
    use decomp_2d, only: mytype, ph1, real_type, nrank
    use var, only : nzmsize
    use param, only : npress
    use variables, only : nxm, nym, nzm
    use dbg_schemes, only: sqrt_prec

    implicit none

    !! INPUTS
    real(mytype), intent(inout), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize) :: dv3
    real(mytype), intent(in) :: atol, rtol
    integer, intent(in) :: poissiter

    !! OUTPUTS
    logical, intent(out) :: converged

    !! LOCALS
    integer :: ierr
    real(mytype) :: errloc, errglob, divup3norm

    if (poissiter.eq.0) then
       errloc = sum(dv3**2)
       CALL MPI_ALLREDUCE(errloc,divup3norm,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
       divup3norm = sqrt_prec(divup3norm / nxm / nym / nzm)

       if (nrank.eq.0) then
          write(*,*)  "solving variable-coefficient poisson equation:"
          write(*,*)  "+ rms div(u*) - div(u): ", divup3norm
       endif
    else
       !! Compute RMS change
       errloc = sum((pp3(:,:,:,1) - pp3(:,:,:,2))**2)
       call MPI_ALLREDUCE(errloc,errglob,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
       errglob = sqrt_prec(errglob / nxm / nym / nzm)

       if (nrank.eq.0) then
          write(*,*)  "+ RMS change in pressure: ", errglob
       endif

       if (errglob.le.atol) then
          converged = .true.
          if (nrank.eq.0) then
             write(*,*)  "- Converged: atol"
          endif
       endif

       !! Compare RMS change to size of |div(u*) - div(u)|
       if (errglob.lt.(rtol * divup3norm)) then
          converged = .true.
          if (nrank.eq.0) then
             write(*,*)  "- Converged: rtol"
          endif
       endif

       if (.not.converged) then
          pp3(:,:,:,2) = pp3(:,:,:,1)
       endif
    endif

  endsubroutine test_varcoeff
  !############################################################################
  !!
  !!  SUBROUTINE: calc_varcoeff_rhs
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Computes RHS of the variable-coefficient Poisson solver
  !!
  !############################################################################
  subroutine calc_varcoeff_rhs(pp3, rho1, px1, py1, pz1, dv3, drho1, ep1, divu3, rho0, poissiter)

    use MPI

    use decomp_2d

    use param, only : nrhotime, ntime
    use param, only : one

    use var, only : ta1, tb1, tc1
    use var, only : nzmsize

    implicit none

    !! INPUTS
    integer, intent(in) :: poissiter
    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: px1, py1, pz1
    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), ntime) :: drho1
    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ep1
    real(mytype), intent(in), dimension(zsize(1), zsize(2), zsize(3)) :: divu3
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize) :: dv3
    real(mytype) :: rho0

    !! OUTPUTS
    real(mytype), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize) :: pp3

    !! LOCALS
    integer :: nlock, ierr
    real(mytype) :: rhomin

    if (poissiter.eq.0) then
       !! Compute rho0
       rhomin = minval(rho1)

       call MPI_ALLREDUCE(rhomin,rho0,1,real_type,MPI_MIN,MPI_COMM_WORLD,ierr)
    endif

    ta1(:,:,:) = (one - rho0 / rho1(:,:,:,1)) * px1(:,:,:)
    tb1(:,:,:) = (one - rho0 / rho1(:,:,:,1)) * py1(:,:,:)
    tc1(:,:,:) = (one - rho0 / rho1(:,:,:,1)) * pz1(:,:,:)

    nlock = -1 !! Don't do any funny business with LMN
    call divergence(pp3,rho1,ta1,tb1,tc1,ep1,drho1,divu3,nlock)

    !! lapl(p) = div((1 - rho0/rho) grad(p)) + rho0(div(u*) - div(u))
    !! dv3 contains div(u*) - div(u)
    pp3(:,:,:) = pp3(:,:,:) + rho0 * dv3(:,:,:)

  endsubroutine calc_varcoeff_rhs
  !############################################################################
endmodule navier
