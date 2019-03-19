!************************************************************************
! Time-marching subroutine used to advance the numerical solution in time  
!
!        input: dvar1
! 
! input/output: var1
!
!       Author: Yorgos
!     Modified: Paul
! 
!************************************************************************ 
subroutine  int_time(var1,dvar1)

  USE param
  USE variables
  USE decomp_2d
  implicit none

  integer :: ijk,nxyz

  !! INPUTS
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: var1

  !! OUTPUTS
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dvar1

#ifdef DEBG
  if (nrank .eq. 0) print *,'# intt start'
#endif

  if (itimescheme.eq.1) then 
     !>>> Euler 
     var1(:,:,:)=gdt(itr)*dvar1(:,:,:,1)+var1(:,:,:) 
  elseif(itimescheme.eq.2) then
     !>>> Adam-Bashforth second order (AB2)
     
     ! Do first time step with Euler
     if(itime.eq.1.and.irestart.eq.0) then
        var1(:,:,:)=gdt(itr)*dvar1(:,:,:,1)+var1(:,:,:) 
     else
        var1(:,:,:)=adt(itr)*dvar1(:,:,:,1)+bdt(itr)*dvar1(:,:,:,2)+var1(:,:,:)
     endif
     dvar1(:,:,:,2)=dvar1(:,:,:,1) 
  elseif(itimescheme.eq.3) then
     !>>> Adams-Bashforth third order (AB3)
     
     ! Do first time step with Euler
     if(itime.eq.1.and.irestart.eq.0) then
        var1(:,:,:)=dt*dvar1(:,:,:,1)+var1(:,:,:)
     elseif(itime.eq.2.and.irestart.eq.0) then
      	! Do second time step with AB2
      	var1(:,:,:)=onepfive*dt*dvar1(:,:,:,1)-half*dt*dvar1(:,:,:,2)+var1(:,:,:)
        dvar1(:,:,:,3)=dvar1(:,:,:,2) 
     else
      	! Finally using AB3
      	var1(:,:,:)=adt(itr)*dvar1(:,:,:,1)+bdt(itr)*dvar1(:,:,:,2)+cdt(itr)*dvar1(:,:,:,3)+var1(:,:,:)
        dvar1(:,:,:,3)=dvar1(:,:,:,2) 
     endif
     dvar1(:,:,:,2)=dvar1(:,:,:,1) 
  elseif(itimescheme.eq.4) then
     !>>> Adams-Bashforth fourth order (AB4)
     
     !if (itime.eq.1.and.ilit.eq.0) then
     !var(:,:,:)=gdt(itr)*hx(:,:,:)+var(:,:,:)
     !uy(:,:,:)=gdt(itr)*hy(:,:,:)+uy(:,:,:) 
     !uz(:,:,:)=gdt(itr)*hz(:,:,:)+uz(:,:,:)
     !gx(:,:,:)=hx(:,:,:)
     !gy(:,:,:)=hy(:,:,:)
     !gz(:,:,:)=hz(:,:,:)            
     !elseif (itime.eq.2.and.ilit.eq.0) then 	   
     !var(:,:,:)=adt(itr)*hx(:,:,:)+bdt(itr)*gx(:,:,:)+var(:,:,:)
     !uy(:,:,:)=adt(itr)*hy(:,:,:)+bdt(itr)*gy(:,:,:)+uy(:,:,:)   
     !uz(:,:,:)=adt(itr)*hz(:,:,:)+bdt(itr)*gz(:,:,:)+uz(:,:,:)
     !gox(:,:,:)=gx(:,:,:)
     !goy(:,:,:)=gy(:,:,:)
     !goz(:,:,:)=gz(:,:,:)
     !gx(:,:,:)=hx(:,:,:)
     !gy(:,:,:)=hy(:,:,:)
     !gz(:,:,:)=hz(:,:,:)            
     !elseif (itime.eq.3.and.ilit.eq.0) then 
     !var(:,:,:)=adt(itr)*hx(:,:,:)+bdt(itr)*gx(:,:,:)+cdt(itr)*gox(:,:,:)+var(:,:,:)
     !uy(:,:,:)=adt(itr)*hy(:,:,:)+bdt(itr)*gy(:,:,:)+cdt(itr)*goy(:,:,:)+uy(:,:,:)   
     !uz(:,:,:)=adt(itr)*hz(:,:,:)+bdt(itr)*gz(:,:,:)+cdt(itr)*goz(:,:,:)+uz(:,:,:)
     !gox(:,:,:)=gx(:,:,:)
     !goy(:,:,:)=gy(:,:,:)
     !goz(:,:,:)=gz(:,:,:)
     !gx(:,:,:)=hx(:,:,:)
     !gy(:,:,:)=hy(:,:,:)
     !gz(:,:,:)=hz(:,:,:)            
     !else 
     !var(:,:,:)=adt(itr)*hx(:,:,:)+bdt(itr)*gx(:,:,:)+cdt(itr)*gox(:,:,:)+ddt(itr)*gax(:,:,:)+var(:,:,:)
     !uy(:,:,:)=adt(itr)*hy(:,:,:)+bdt(itr)*gy(:,:,:)+cdt(itr)*goy(:,:,:)+ddt(itr)*gay(:,:,:)+uy(:,:,:)   
     !uz(:,:,:)=adt(itr)*hz(:,:,:)+bdt(itr)*gz(:,:,:)+cdt(itr)*goz(:,:,:)+ddt(itr)*gaz(:,:,:)+uz(:,:,:)
     !gax(:,:,:)=gox(:,:,:)
     !gay(:,:,:)=goy(:,:,:)
     !gaz(:,:,:)=goz(:,:,:)		     
     !gox(:,:,:)=gx(:,:,:)
     !goy(:,:,:)=gy(:,:,:)
     !goz(:,:,:)=gz(:,:,:)
     !gx(:,:,:)=hx(:,:,:)
     !gy(:,:,:)=hy(:,:,:)
     !gz(:,:,:)=hz(:,:,:)            
     !endif	
     !>>> Runge-Kutta (low storage) RK3
  elseif(itimescheme.eq.5) then
     if(itr.eq.1) then
        var1(:,:,:)=gdt(itr)*dvar1(:,:,:,1)+var1(:,:,:) 
     else
        var1(:,:,:)=adt(itr)*dvar1(:,:,:,1)+bdt(itr)*dvar1(:,:,:,2)+var1(:,:,:)
     endif
     dvar1(:,:,:,2)=dvar1(:,:,:,1) 
     !>>> Runge-Kutta (low storage) RK4
  elseif(itimescheme.eq.6) then

  endif

#ifdef DEBG
  if (nrank .eq. 0) print *,'# intt done'
#endif

  return

end subroutine int_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  SUBROUTINE: int_time_momentum
!! DESCRIPTION: Integrates the momentum equations in time by calling time
!!              integrator.
!!      INPUTS: dux1, duy1, duz1 - the RHS(s) of the momentum equations
!!     OUTPUTS: ux1,   uy1,  uz1 - the intermediate momentum state.
!!       NOTES: This is integrating the MOMENTUM in time (!= velocity)
!!      AUTHOR: Paul Bartholomew
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine int_time_momentum(ux1, uy1, uz1, dux1, duy1, duz1)

  USE param
  USE variables
  USE decomp_2d

  implicit none

  !! INPUTS
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1

  !! OUTPUTS
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1, duy1, duz1

  call int_time(ux1, dux1)
  call int_time(uy1, duy1)
  call int_time(uz1, duz1)
  
endsubroutine int_time_momentum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  SUBROUTINE: int_time_continuity
!! DESCRIPTION: Integrates the continuity (aka density transport) equation in
!!              time
!!      INPUTS: drho1 - the RHS(s) of the continuity equation.
!!     OUTPUTS:  rho1 - the density at new time.
!!      AUTHOR: Paul Bartholomew
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine int_time_continuity(rho1, drho1)

  USE param
  USE variables
  USE decomp_2d
  
  implicit none

  integer :: it, i, j, k
  real(mytype) :: rhomin, rhomax

  !! INPUTS
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1

  !! OUTPUTS
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: drho1

  !! First, update old density / store old transients depending on scheme
  if (itimescheme.lt.5) then
     !! Euler/AB - Store old density values
     do it = nrhotime, 2, -1
        rho1(:,:,:,it) = rho1(:,:,:,it-1)
     enddo
  elseif (itimescheme.eq.5) then
     !! RK3 - Stores old transients
     if (itr.eq.1) then
        do it = nrhotime, 2, -1
           rho1(:,:,:,it) = rho1(:,:,:,it-1)
        enddo
        rho1(:,:,:,2) = drho1(:,:,:,1)
     endif
  else
     if (nrank.eq.0) then
        print *, "int_time_continuity not implemented for itimescheme", itimescheme
        stop
     endif
  endif

  !! Now we can update current density
  call int_time(rho1(:,:,:,1), drho1)

  !! Enforce boundedness on density
  rhomin = min(dens1, dens2)
  rhomax = max(dens1, dens2)
  do k = 1, xsize(3)
     do j = 1, xsize(2)
        do i = 1, xsize(1)
           rho1(i, j, k, 1) = max(rho1(i, j, k, 1), rhomin)
           rho1(i, j, k, 1) = min(rho1(i, j, k, 1), rhomax)
        enddo
     enddo
  enddo
  
endsubroutine int_time_continuity

!********************************************************************
!subroutine CORPG
!Correction of u* by the pressure gradient to get a divergence free
!field
! input : px,py,pz
! output : ux,uy,uz
!written by SL 2018
!********************************************************************    
subroutine corpg (ux,uy,uz,px,py,pz)
  
  USE decomp_2d
  USE variables
  USE param
  USE var, only: ta2
  USE MPI
  
  implicit none

  integer :: ijk,nxyz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,px,py,pz

  ux(:,:,:)=-px(:,:,:)+ux(:,:,:)
  uy(:,:,:)=-py(:,:,:)+uy(:,:,:)
  uz(:,:,:)=-pz(:,:,:)+uz(:,:,:)
  
  return
end subroutine corpg
!********************************************************************
!subroutine DIVERGENCe
!Calculation of div u* for nlock=1 and of div u^{n+1} for nlock=2
! input : ux1,uy1,uz1,ep1 (on velocity mesh)
! output : pp3 (on pressure mesh)
!written by SL 2018
!******************************************************************** 
subroutine divergence (pp3,rho1,ux1,uy1,uz1,ep1,drho1,divu3,nlock)

  USE param
  USE decomp_2d
  USE variables
  USE var, ONLY: ta1, tb1, tc1, pp1, pgy1, pgz1, di1, &
       duxdxp2, uyp2, uzp2, duydypi2, upi2, ta2, dipp2, &
       duxydxyp3, uzp3, po3, dipp3, nxmsize, nymsize, nzmsize
  USE MPI

  implicit none

!  TYPE(DECOMP_INFO) :: ph1,ph3,ph4

  !X PENCILS NX NY NZ  -->NXM NY NZ
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1,ep1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime),intent(in) :: drho1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime),intent(in) :: rho1
  !Z PENCILS NXM NYM NZ  -->NXM NYM NZM
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)),intent(in) :: divu3
  real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize) :: pp3

  integer :: ijk,nvect1,nvect2,nvect3,i,j,k,nlock
  integer :: code
  real(mytype) :: tmax,tmoy,tmax1,tmoy1

  nvect1=xsize(1)*xsize(2)*xsize(3)
  nvect2=ysize(1)*ysize(2)*ysize(3)
  nvect3=(ph1%zen(1)-ph1%zst(1)+1)*(ph1%zen(2)-ph1%zst(2)+1)*nzmsize

  if (iibm.eq.0) then
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
  if (ilmn) then
     if (nlock.eq.1) then
        !! Approximate -div(rho u) using ddt(rho)
        call extrapol_drhodt(ta1, rho1, drho1)
     elseif (nlock.eq.2) then
        !! Need to check our error against divu constraint
        call transpose_z_to_y(-divu3, ta2)
        call transpose_y_to_x(ta2, ta1)
     endif
     call interxvp(pgy1,ta1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)
  else
     pgy1(:,:,:) = zero
  endif

  pp1(:,:,:) = pp1(:,:,:) + pgy1(:,:,:) !! pgy1 is a contribution from LMN, zero otherwise
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
  tmoy=0._mytype
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

  if (nrank==0) then
    if (nlock==2) then
      print *,'DIV U  max mean=',real(tmax1,4),real(tmoy1/real(nproc),4)
    else
      print *,'DIV U* max mean=',real(tmax1,4),real(tmoy1/real(nproc),4)
    endif
  endif

  return
end subroutine divergence


!********************************************************************
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
!********************************************************************  
subroutine gradp(px1,py1,pz1,pp3)

  USE param
  USE decomp_2d
  USE variables
  USE MPI
  USE var, only: pp1,pgy1,pgz1,di1,pp2,ppi2,pgy2,pgz2,pgzi2,dip2,&
          pgz3,ppi3,dip3,nxmsize,nymsize,nzmsize
#ifdef FORCES
  USE forces, only : ppi1
#endif
  
  implicit none

  integer :: i,j,k,code
  integer, dimension(2) :: dims, dummy_coords
  logical, dimension(2) :: dummy_periods

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

#ifdef FORCES
  call interxpv(ppi1,pp1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
  nxmsize,xsize(1),xsize(2),xsize(3),1)
#endif
  
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
    if (xsize(2)==1) then
      do k=1,xsize(3)
        do i=1,xsize(1)
          dpdxy1(i,k)=px1(i,1,k)/gdt(itr)
          dpdzy1(i,k)=pz1(i,1,k)/gdt(itr)
        enddo
      enddo
    endif
  endif
  if (nclyn.eq.2) then
    if (xsize(2)==ny) then
      do k=1,xsize(3)
        do i=1,xsize(1)
          dpdxyn(i,k)=px1(i,ny,k)/gdt(itr)
          dpdzyn(i,k)=pz1(i,ny,k)/gdt(itr)
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
!*******************************************************************
subroutine pre_correc(ux,uy,uz,ep)

  USE decomp_2d
  USE variables
  USE param
  USE var
  USE MPI

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,ep
  integer :: i,j,k,code,is
  real(mytype) :: ut,ut1,utt,ut11

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

  !********NCLX==2*************************************
  !we are in X pencils:
  if (nclx1==2.and.nclxn==2) then

     !Computatation of the flow rate Inflow/Outflow
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
     if (nrank==0) print *,'Flow rate x I/O/O-I',real(ut11,4),real(utt,4),real(utt-ut11,4)
     do k=1,xsize(3)
        do j=1,xsize(2)
           bxxn(j,k)=bxxn(j,k)-utt+ut11
        enddo
     enddo

  endif

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

  return
end subroutine pre_correc
!*******************************************************************

!! Convert to/from conserved/primary variables
SUBROUTINE primary_to_conserved(rho1, var1)

  USE decomp_2d, ONLY : mytype, xsize
  USE param, ONLY : nrhotime
  
  IMPLICIT NONE

  REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: var1

  var1(:,:,:) = rho1(:,:,:,1) * var1(:,:,:)

ENDSUBROUTINE primary_to_conserved
SUBROUTINE conserved_to_primary(rho1, var1)

  USE decomp_2d, ONLY : mytype, xsize
  USE param, ONLY : nrhotime
  
  IMPLICIT NONE

  REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: var1

  var1(:,:,:) = var1(:,:,:) / rho1(:,:,:,1)

ENDSUBROUTINE conserved_to_primary

!! Calculate velocity-divergence constraint
SUBROUTINE calc_divu_constraint(divu3, rho1)

  USE decomp_2d, ONLY : mytype, xsize, ysize, zsize
  USE decomp_2d, ONLY : transpose_x_to_y, transpose_y_to_z
  USE param, ONLY : nrhotime, zero, ilmn, pressure0
  USE param, ONLY : xnu, prandtl
  USE variables

  USE var, ONLY : ta1, tb1, di1
  USE var, ONLY : ta2, tb2, tc2, td2, di2
  USE var, ONLY : ta3, tb3, di3
  
  IMPLICIT NONE

  INTEGER :: i, j, k

  REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
  REAL(mytype), INTENT(OUT), DIMENSION(zsize(1), zsize(2), zsize(3)) :: divu3

  IF (ilmn) THEN
     !!------------------------------------------------------------------------------
     !! X-pencil
     
     !! We need temperature
     ta1(:,:,:) = pressure0 / rho1(:,:,:,1)
     
     CALL derxx (tb1, ta1, di1, sx, sfxp, ssxp, swxp, xsize(1), xsize(2), xsize(3), 1)
     
     CALL transpose_x_to_y(ta1, ta2)
     CALL transpose_x_to_y(tb1, tb2)
     
     !!------------------------------------------------------------------------------
     !! Y-pencil
     CALL deryy (tc2, ta2, di2, sy, sfyp, ssyp, swyp, ysize(1), ysize(2), ysize(3), 1)

     tb2(:,:,:) = tb2(:,:,:) + tc2(:,:,:)

     CALL transpose_y_to_z(ta2, ta3)
     CALL transpose_y_to_z(tb2, tb3)
     
     !!------------------------------------------------------------------------------
     !! Z-pencil
     CALL derzz (divu3, ta3, di3, sz, sfzp, sszp, swzp, zsize(1), zsize(2), zsize(3), 1)

     divu3(:,:,:) = divu3(:,:,:) + tb3(:,:,:)
     divu3(:,:,:) = (xnu / prandtl) * divu3(:,:,:)

     divu3(:,:,:) = divu3(:,:,:) / pressure0
  ELSE
     divu3(:,:,:) = zero
  ENDIF

ENDSUBROUTINE calc_divu_constraint

SUBROUTINE extrapol_drhodt(drhodt1_next, rho1, drho1)

  USE decomp_2d, ONLY : mytype, xsize, nrank
  USE param, ONLY : ntime, nrhotime, itime, itimescheme, itr, dt, gdt, irestart
  USE param, ONLY : half, three, four
  
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
  
ENDSUBROUTINE extrapol_drhodt

