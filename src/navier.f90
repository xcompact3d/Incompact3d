!************************************************************************
! Time-marching subroutine used to advance the numerical solution in time  
!
! output: 
! 
! input: ux, uy, uz 
! 
!************************************************************************ 
subroutine  int_time_momentum(ux1,uy1,uz1,dux1,duy1,duz1)

  USE param
  USE variables
  USE decomp_2d
  implicit none

  integer :: ijk,nxyz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1

#ifdef DEBG
  if (nrank .eq. 0) print *,'# intt start'
#endif

  if (nscheme.eq.1) then 
     !>>> Euler 
     ux1(:,:,:)=gdt(itr)*dux1(:,:,:,1)+ux1(:,:,:) 
     uy1(:,:,:)=gdt(itr)*duy1(:,:,:,1)+uy1(:,:,:)
     uz1(:,:,:)=gdt(itr)*duz1(:,:,:,1)+uz1(:,:,:)
  elseif(nscheme.eq.2) then
     !>>> Adam-Bashforth second order (AB2)
     
     ! Do first time step with Euler
     if(itime.eq.1.and.ilit.eq.0) then
        ux1(:,:,:)=gdt(itr)*dux1(:,:,:,1)+ux1(:,:,:) 
        uy1(:,:,:)=gdt(itr)*duy1(:,:,:,1)+uy1(:,:,:)
        uz1(:,:,:)=gdt(itr)*duz1(:,:,:,1)+uz1(:,:,:)
     else
        ux1(:,:,:)=adt(itr)*dux1(:,:,:,1)+bdt(itr)*dux1(:,:,:,2)+ux1(:,:,:)
        uy1(:,:,:)=adt(itr)*duy1(:,:,:,1)+bdt(itr)*duy1(:,:,:,2)+uy1(:,:,:)
        uz1(:,:,:)=adt(itr)*duz1(:,:,:,1)+bdt(itr)*duz1(:,:,:,2)+uz1(:,:,:)
     endif
     dux1(:,:,:,2)=dux1(:,:,:,1) 
     duy1(:,:,:,2)=duy1(:,:,:,1) 
     duz1(:,:,:,2)=duz1(:,:,:,1) 
  elseif(nscheme.eq.3) then
     !>>> Adams-Bashforth third order (AB3)
     
     ! Do first time step with Euler
     if(itime.eq.1.and.ilit.eq.0) then
        ux1(:,:,:)=dt*dux1(:,:,:,1)+ux1(:,:,:)
        uy1(:,:,:)=dt*duy1(:,:,:,1)+uy1(:,:,:)
        uz1(:,:,:)=dt*duz1(:,:,:,1)+uz1(:,:,:)
     elseif(itime.eq.2.and.ilit.eq.0) then
      	! Do second time step with AB2
      	ux1(:,:,:)=onepfive*dt*dux1(:,:,:,1)-half*dt*dux1(:,:,:,2)+ux1(:,:,:)
      	uy1(:,:,:)=onepfive*dt*duy1(:,:,:,1)-half*dt*duy1(:,:,:,2)+uy1(:,:,:)
      	uz1(:,:,:)=onepfive*dt*duz1(:,:,:,1)-half*dt*duz1(:,:,:,2)+uz1(:,:,:)
        dux1(:,:,:,3)=dux1(:,:,:,2) 
        duy1(:,:,:,3)=duy1(:,:,:,2) 
        duz1(:,:,:,3)=duz1(:,:,:,2) 
     else
      	! Finally using AB3
      	ux1(:,:,:)=adt(itr)*dux1(:,:,:,1)+bdt(itr)*dux1(:,:,:,2)+cdt(itr)*dux1(:,:,:,3)+ux1(:,:,:)
      	uy1(:,:,:)=adt(itr)*duy1(:,:,:,1)+bdt(itr)*duy1(:,:,:,2)+cdt(itr)*duy1(:,:,:,3)+uy1(:,:,:)
      	uz1(:,:,:)=adt(itr)*duz1(:,:,:,1)+bdt(itr)*duz1(:,:,:,2)+cdt(itr)*duz1(:,:,:,3)+uz1(:,:,:)
        dux1(:,:,:,3)=dux1(:,:,:,2) 
        duy1(:,:,:,3)=duy1(:,:,:,2) 
        duz1(:,:,:,3)=duz1(:,:,:,2) 
     endif
     dux1(:,:,:,2)=dux1(:,:,:,1) 
     duy1(:,:,:,2)=duy1(:,:,:,1) 
     duz1(:,:,:,2)=duz1(:,:,:,1) 
  elseif(nscheme.eq.4) then
     !>>> Adams-Bashforth fourth order (AB4)
     
     !if (itime.eq.1.and.ilit.eq.0) then
     !ux(:,:,:)=gdt(itr)*hx(:,:,:)+ux(:,:,:)
     !uy(:,:,:)=gdt(itr)*hy(:,:,:)+uy(:,:,:) 
     !uz(:,:,:)=gdt(itr)*hz(:,:,:)+uz(:,:,:)
     !gx(:,:,:)=hx(:,:,:)
     !gy(:,:,:)=hy(:,:,:)
     !gz(:,:,:)=hz(:,:,:)            
     !elseif (itime.eq.2.and.ilit.eq.0) then 	   
     !ux(:,:,:)=adt(itr)*hx(:,:,:)+bdt(itr)*gx(:,:,:)+ux(:,:,:)
     !uy(:,:,:)=adt(itr)*hy(:,:,:)+bdt(itr)*gy(:,:,:)+uy(:,:,:)   
     !uz(:,:,:)=adt(itr)*hz(:,:,:)+bdt(itr)*gz(:,:,:)+uz(:,:,:)
     !gox(:,:,:)=gx(:,:,:)
     !goy(:,:,:)=gy(:,:,:)
     !goz(:,:,:)=gz(:,:,:)
     !gx(:,:,:)=hx(:,:,:)
     !gy(:,:,:)=hy(:,:,:)
     !gz(:,:,:)=hz(:,:,:)            
     !elseif (itime.eq.3.and.ilit.eq.0) then 
     !ux(:,:,:)=adt(itr)*hx(:,:,:)+bdt(itr)*gx(:,:,:)+cdt(itr)*gox(:,:,:)+ux(:,:,:)
     !uy(:,:,:)=adt(itr)*hy(:,:,:)+bdt(itr)*gy(:,:,:)+cdt(itr)*goy(:,:,:)+uy(:,:,:)   
     !uz(:,:,:)=adt(itr)*hz(:,:,:)+bdt(itr)*gz(:,:,:)+cdt(itr)*goz(:,:,:)+uz(:,:,:)
     !gox(:,:,:)=gx(:,:,:)
     !goy(:,:,:)=gy(:,:,:)
     !goz(:,:,:)=gz(:,:,:)
     !gx(:,:,:)=hx(:,:,:)
     !gy(:,:,:)=hy(:,:,:)
     !gz(:,:,:)=hz(:,:,:)            
     !else 
     !ux(:,:,:)=adt(itr)*hx(:,:,:)+bdt(itr)*gx(:,:,:)+cdt(itr)*gox(:,:,:)+ddt(itr)*gax(:,:,:)+ux(:,:,:)
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
  elseif(nscheme.eq.5) then
     if(itr.eq.1) then
        ux1(:,:,:)=gdt(itr)*dux1(:,:,:,1)+ux1(:,:,:) 
        uy1(:,:,:)=gdt(itr)*duy1(:,:,:,1)+uy1(:,:,:)
        uz1(:,:,:)=gdt(itr)*duz1(:,:,:,1)+uz1(:,:,:)
     else
        ux1(:,:,:)=adt(itr)*dux1(:,:,:,1)+bdt(itr)*dux1(:,:,:,2)+ux1(:,:,:)
        uy1(:,:,:)=adt(itr)*duy1(:,:,:,1)+bdt(itr)*duy1(:,:,:,2)+uy1(:,:,:)
        uz1(:,:,:)=adt(itr)*duz1(:,:,:,1)+bdt(itr)*duz1(:,:,:,2)+uz1(:,:,:)
     endif
     dux1(:,:,:,2)=dux1(:,:,:,1) 
     duy1(:,:,:,2)=duy1(:,:,:,1) 
     duz1(:,:,:,2)=duz1(:,:,:,1) 
     !>>> Runge-Kutta (low storage) RK4
  elseif(nscheme.eq.6) then

  endif

#ifdef DEBG
  if (nrank .eq. 0) print *,'# intt done'
#endif

  return

end subroutine int_time_momentum

!********************************************************************
!subroutine CORPG
!Correction of u* by the pressure gradient to get a divergence free
!field
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
subroutine divergence (ux1,uy1,uz1,ep1,pp3,&
  nxmsize,nymsize,nzmsize,nlock)

  USE param
  USE decomp_2d
  USE variables
  USE var, ONLY: ta1, tb1, tc1, pp1, pgy1, pgz1, di1, &
       duxdxp2, uyp2, uzp2, duydypi2, upi2, dipp2, &
       duxydxyp3, uzp3, po3, dipp3
  USE MPI

  implicit none

!  TYPE(DECOMP_INFO) :: ph1,ph3,ph4

  !X PENCILS NX NY NZ  -->NXM NY NZ
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1,ep1
  !Z PENCILS NXM NYM NZ  -->NXM NYM NZM
  real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize) :: pp3

  integer :: ijk,nvect1,nvect2,nvect3,i,j,k,nlock
  integer :: nxmsize,nymsize,nzmsize,code
  real(mytype) :: tmax,tmoy,tmax1,tmoy1

  nvect1=xsize(1)*xsize(2)*xsize(3)
  nvect2=ysize(1)*ysize(2)*ysize(3)
  nvect3=(ph1%zen(1)-ph1%zst(1)+1)*(ph1%zen(2)-ph1%zst(2)+1)*nzmsize

  if (ivirt.eq.0.and.ilag.eq.0) then
     ta1(:,:,:) = ux1(:,:,:)
     tb1(:,:,:) = uy1(:,:,:)
     tc1(:,:,:) = uz1(:,:,:)
  else
     ta1(:,:,:) = (one - ep1(:,:,:)) * ux1(:,:,:)
     tb1(:,:,:) = (one - ep1(:,:,:)) * uy1(:,:,:)
     tc1(:,:,:) = (one - ep1(:,:,:)) * uz1(:,:,:)
  endif

  !WORK X-PENCILS
  call decx6(pp1,ta1,di1,sx,cfx6,csx6,cwx6,xsize(1),nxmsize,xsize(2),xsize(3),0)
  call inter6(pgy1,tb1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)
  call inter6(pgz1,tc1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)

  call transpose_x_to_y(pp1,duxdxp2,ph4)!->NXM NY NZ
  call transpose_x_to_y(pgy1,uyp2,ph4)
  call transpose_x_to_y(pgz1,uzp2,ph4)

  !WORK Y-PENCILS
  call intery6(upi2,duxdxp2,dipp2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),1)
  call decy6(duydypi2,uyp2,dipp2,sy,cfy6,csy6,cwy6,ppyi,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),0)

  !! Compute sum dudx + dvdy
  duydypi2(:,:,:) = duydypi2(:,:,:) + upi2(:,:,:)

  call intery6(upi2,uzp2,dipp2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),1)
  
  call transpose_y_to_z(duydypi2,duxydxyp3,ph3)!->NXM NYM NZ
  call transpose_y_to_z(upi2,uzp3,ph3)

  !WORK Z-PENCILS
  call interz6(pp3,duxydxyp3,dipp3,sz,cifzp6,ciszp6,ciwzp6,(ph1%zen(1)-ph1%zst(1)+1),&
  (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,1)
  call decz6(po3,uzp3,dipp3,sz,cfz6,csz6,cwz6,(ph1%zen(1)-ph1%zst(1)+1),&
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
! INPUT: pp3 - pressure field (on pressure mesh)
! OUTPUT: px1, py1, pz1 - pressure gradients (on velocity mesh)
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
  call interiz6(ppi3,pp3,dip3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
  (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
  call deciz6(pgz3,pp3,dip3,sz,cfip6z,csip6z,cwip6z,cfz6,csz6,cwz6,&
  (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)

  !WORK Y-PENCILS
  call transpose_z_to_y(pgz3,pgz2,ph3) !nxm nym nz
  call transpose_z_to_y(ppi3,pp2,ph3)

  call interiy6(ppi2,pp2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
  (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
  call deciy6(pgy2,pp2,dip2,sy,cfip6y,csip6y,cwip6y,cfy6,csy6,cwy6,ppy,&
  (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
  call interiy6(pgzi2,pgz2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
  (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)

  !WORK X-PENCILS

  call transpose_y_to_x(ppi2,pp1,ph2) !nxm ny nz
  call transpose_y_to_x(pgy2,pgy1,ph2)
  call transpose_y_to_x(pgzi2,pgz1,ph2)

  call deci6(px1,pp1,di1,sx,cfip6,csip6,cwip6,cfx6,csx6,cwx6,&
  nxmsize,xsize(1),xsize(2),xsize(3),1)
  call interi6(py1,pgy1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
  nxmsize,xsize(1),xsize(2),xsize(3),1)
  call interi6(pz1,pgz1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
  nxmsize,xsize(1),xsize(2),xsize(3),1)

#ifdef FORCES
  call interi6(ppi1,pp1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
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
subroutine corgp_IBM (ux,uy,uz,px,py,pz,nlock)
  USE param
  USE decomp_2d
  USE variables
  implicit none
  integer :: ijk,nlock,nxyz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,px,py,pz
  nxyz=xsize(1)*xsize(2)*xsize(3)
  if (nlock.eq.1) then
    if (nz.gt.1) then
      do ijk=1,nxyz
        uy(ijk,1,1)=-py(ijk,1,1)+uy(ijk,1,1)
        uz(ijk,1,1)=-pz(ijk,1,1)+uz(ijk,1,1)
        ux(ijk,1,1)=-px(ijk,1,1)+ux(ijk,1,1)
      enddo
    else
      do ijk=1,nxyz
        uy(ijk,1,1)=-py(ijk,1,1)+uy(ijk,1,1)
        ux(ijk,1,1)=-px(ijk,1,1)+ux(ijk,1,1)
      enddo
    endif
  endif
  if (nlock.eq.2) then
    if (nz.gt.1) then
      do ijk=1,nxyz
        uy(ijk,1,1)=py(ijk,1,1)+uy(ijk,1,1)
        uz(ijk,1,1)=pz(ijk,1,1)+uz(ijk,1,1)
        ux(ijk,1,1)=px(ijk,1,1)+ux(ijk,1,1)
      enddo
    else
      do ijk=1,nxyz
        uy(ijk,1,1)=py(ijk,1,1)+uy(ijk,1,1)
        ux(ijk,1,1)=px(ijk,1,1)+ux(ijk,1,1)
      enddo
    endif
  endif

  return
end subroutine corgp_IBM
!*******************************************************************
#ifdef IBM
subroutine body(ux1,uy1,uz1,ep1,arg)
  USE param
  USE decomp_2d
  USE decomp_2d_io
  USE variables
  implicit none
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
  integer :: arg,ijk,nvect1

#ifdef DEBG
  if (nrank .eq. 0) print *,'# body start'
#endif


  if (arg==0) then !First execution, initt epsi
     ep1(:,:,:)=zero
     call geomcomplex(ep1,xstart(1),xend(1),ny,xstart(2),xend(2),xstart(3),xend(3),dx,yp,dz,one)
  elseif (arg==1) then  !Any other iteration
     nvect1=xsize(1)*xsize(2)*xsize(3)
     do ijk=1,nvect1
        ux1(ijk,1,1)=(one-ep1(ijk,1,1))*ux1(ijk,1,1)
        uy1(ijk,1,1)=(one-ep1(ijk,1,1))*uy1(ijk,1,1)
        uz1(ijk,1,1)=(one-ep1(ijk,1,1))*uz1(ijk,1,1)
     enddo
  endif

  !X PENCILS

#ifdef DEBG
  if (nrank .eq. 0) print *,'# body done'
#endif

  return
end subroutine body
#endif
!*******************************************************************
subroutine square(ycenter,zcenter,xthick,xlenght,ux,uy,uz,esp)
  USE param
  USE decomp_2d
  USE variables
  implicit none
  real(mytype), dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ux,uy,uz,esp
  real(mytype) :: ycenter,zcenter,xthick,xlenght,slenght,xz1,xz2,xy1,xy2
  integer :: j1,j2,z1,z2,iep,i,j,k,k1,k2,ilen

  iep=int(xthick*ny/yly)
  ilen=int(xlenght*ny/yly)

  j1=int(ycenter*ny/yly)
  k1=int(zcenter*nz/zlz)+ilen
  k2=int(zcenter*nz/zlz)-ilen

  do k=xstart(3),xend(3)
    do j=xstart(2),xend(2)
      do i=xstart(1),xend(1)
        if ((k.ge.(k2-iep)).and.(k.le.(k1+iep))) then
          if ((j.ge.(j1+ilen-iep)).and.(j.le.(j1+ilen+iep))) then
            if ((i.ge.75).and.(i.le.83)) then
              esp(i,j,k)=one
              ux(i,j,k)=zero
              uy(i,j,k)=zero
              uz(i,j,k)=zero
            endif
          endif
        endif
        if ((k.ge.(k2-iep)).and.(k.le.(k1+iep))) then
          if ((j.ge.(j1-ilen-iep)).and.(j.le.(j1-ilen+iep))) then
            if ((i.ge.75).and.(i.le.83)) then
              esp(i,j,k)=one
              ux(i,j,k)=zero
              uy(i,j,k)=zero
              uz(i,j,k)=zero
            endif
          endif
        endif
        if ((k.ge.(k1-iep)).and.(k.le.(k1+iep))) then
          if ((j.ge.(j1-ilen-iep)).and.(j.le.(j1+ilen+iep))) then
            if ((i.ge.75).and.(i.le.83)) then
              esp(i,j,k)=one
              ux(i,j,k)=zero
              uy(i,j,k)=zero
              uz(i,j,k)=zero
            endif
          endif
        endif
        if ((k.ge.(k2-iep)).and.(k.le.(k2+iep))) then
          if ((j.ge.(j1-ilen-iep)).and.(j.le.(j1+ilen+iep))) then
            if ((i.ge.75).and.(i.le.83)) then
              esp(i,j,k)=one
              ux(i,j,k)=zero
              uy(i,j,k)=zero
              uz(i,j,k)=zero
            endif
          endif
        endif
      enddo
    enddo
  enddo


  return
end subroutine square
!*******************************************************************
subroutine forcage_square(ux,uy,uz,esp)
  USE param
  USE decomp_2d
  USE variables
  implicit none
  real(mytype), dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ux,uy,uz,esp
  real(mytype),dimension(nz) :: xx1
  real(mytype),dimension(ny) :: yy1
  integer :: j, i, k, np,i1 ,ii
  real(mytype) :: ep0,ep1,ep2,ep3
  real(mytype) :: l0,l1,l2,l3,l4
  real(mytype), dimension(4) :: l5,l6
  real(mytype) :: y,z
  esp(:,:,:)=zero

  ep3=two/five !0.4  Tr=8.5 grid tmin=2.1, t2=4.2, t1=8.6, tmax=17.5
  ep2=four/five !0.8
  ep1=nine/five!
  ep0=seventeen/five!

  l0=yly/two
  l1=yly/four
  l2=yly/eight
  l3=yly/sixteen
  l4=yly/thirytwo

  !####fisrt fractal iteration##########################
  call square(l0,l0,ep0,l1,ux,uy,uz,esp)
  !####Second fractal iteration#########################
  call square(l0+l1,l0+l1,ep1,l2,ux,uy,uz,esp)
  call square(l0+l1,l1,ep1,l2,ux,uy,uz,esp)
  call square(l1,l0+l1,ep1,l2,ux,uy,uz,esp)
  call square(l1,l1,ep1,l2,ux,uy,uz,esp)
  !####third fractal iteration########################
  call square(l0-l2,l0-l2,ep2,l3,ux,uy,uz,esp)
  call square(l0-l2,l2,ep2,l3,ux,uy,uz,esp)
  call square(l2,l2,ep2,l3,ux,uy,uz,esp)
  call square(l2,l0-l2,ep2,l3,ux,uy,uz,esp)
  call square(l0+l1+l2,l0+l1+l2,ep2,l3,ux,uy,uz,esp)
  call square(l0+l1+l2,l0+l2,ep2,l3,ux,uy,uz,esp)
  call square(l0+l2,l0+l2,ep2,l3,ux,uy,uz,esp)
  call square(l0+l2,l0+l1+l2,ep2,l3,ux,uy,uz,esp)
  call square(l0+l1+l2,l2,ep2,l3,ux,uy,uz,esp)
  call square(l0+l2,l2,ep2,l3,ux,uy,uz,esp)
  call square(l0+l1+l2,l0-l2,ep2,l3,ux,uy,uz,esp)
  call square(l0+l2,l0-l2,ep2,l3,ux,uy,uz,esp)
  call square(l0-l2,l0+l2,ep2,l3,ux,uy,uz,esp)
  call square(l0-l2,l0+l1+l2,ep2,l3,ux,uy,uz,esp)
  call square(l2,l0+l2,ep2,l3,ux,uy,uz,esp)
  call square(l2,l0+l1+l2,ep2,l3,ux,uy,uz,esp)
  !###fourth fractal iteration
  l5(1)=l3;l5(2)=-l3;l5(3)=l3;l5(4)=-l3
  l6(1)=l3;l6(2)=l3;l6(3)=-l3;l6(4)=-l3
  do ii=1,4
    call square(l0-l2+l5(ii),l0-l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l0-l2+l5(ii),l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l2+l5(ii),l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l2+l5(ii),l0-l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l0+l1+l2+l5(ii),l0+l1+l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l0+l1+l2+l5(ii),l0+l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l0+l2+l5(ii),l0+l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l0+l2+l5(ii),l0+l1+l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l0+l1+l2+l5(ii),l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l0+l2+l5(ii),l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l0+l1+l2+l5(ii),l0-l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l0+l2+l5(ii),l0-l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l0-l2+l5(ii),l0+l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l0-l2+l5(ii),l0+l1+l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l2+l5(ii),l0+l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l2+l5(ii),l0+l1+l2+l6(ii),ep3,l4,ux,uy,uz,esp)
  enddo

  return
end subroutine forcage_square
