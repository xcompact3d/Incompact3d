!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module var

  use decomp_2d_constants
  use decomp_2d
  use decomp_2d_mpi
  USE variables
  USE param
  USE complex_geometry
  USE mod_stret, only : stretching

  implicit none
  
  ! define all major arrays here
  real(mytype), save, allocatable, dimension(:,:,:) :: ux1, ux2, ux3, po3, dv3
  real(mytype), save, allocatable, dimension(:,:,:,:) :: pp3
  real(mytype), save, allocatable, dimension(:,:,:) :: uy1, uy2, uy3
  real(mytype), save, allocatable, dimension(:,:,:) :: uz1, uz2, uz3
  real(mytype), save, allocatable, dimension(:,:,:,:) :: rho1, drho1
  real(mytype), save, allocatable, dimension(:,:,:) :: rho2, rho3
  real(mytype), save, allocatable, dimension(:,:,:) :: divu3
  real(mytype), save, allocatable, dimension(:,:,:,:) :: phi1, phi2, phi3
  real(mytype), save, allocatable, dimension(:,:,:) :: px1, py1, pz1
  real(mytype), save, allocatable, dimension(:,:,:) :: ep1, diss1, pre1
  real(mytype), save, allocatable, dimension(:,:,:,:) :: dux1,duy1,duz1  ! Output of convdiff
  real(mytype), save, allocatable, dimension(:,:,:,:,:) :: dphi1
  real(mytype), save, allocatable, dimension(:,:,:) :: mu1,mu2,mu3
  real(mytype), save, allocatable, dimension(:,:,:) :: uxf1, uxf2, uxf3, uyf1, uyf2, uyf3, uzf1, uzf2, uzf3, phif1, phif2, phif3

  !arrays for post processing
  real(mytype), save, allocatable, dimension(:,:,:) :: f1,fm1

  !arrays for statistic collection
  real(mytype), save, allocatable, dimension(:,:,:) :: umean,vmean,wmean,pmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
  real(mytype), save, allocatable, dimension(:,:,:,:) :: phimean,phiphimean,uphimean,vphimean,wphimean

  !arrays for visualization
  real(mytype), save, allocatable, dimension(:,:,:) :: uvisu

  ! define all work arrays here
  real(mytype), save, allocatable, dimension(:,:,:) :: ta1,tb1,tc1,td1,&
       te1,tf1,tg1,th1,ti1,di1
  real(mytype), save, allocatable, dimension(:,:,:) :: pp1,pgy1,pgz1
  real(mytype), save, allocatable, dimension(:,:,:) :: ta2,tb2,tc2,td2,&
       te2,tf2,tg2,th2,ti2,tj2,di2
  real(mytype), save, allocatable, dimension(:,:,:) :: pp2,ppi2,pgy2,pgz2,pgzi2,dip2,dipp2,duxdxp2,uyp2,uzp2,upi2,duydypi2
  real(mytype), save, allocatable, dimension(:,:,:) :: ta3,tb3,tc3,td3,&
       te3,tf3,tg3,th3,ti3,di3
  real(mytype), save, allocatable, dimension(:,:,:) :: pgz3,ppi3,dip3,dipp3,duxydxyp3,uzp3

  integer, save :: nxmsize, nymsize, nzmsize

  ! working arrays for LES
  ! These are by far too many variables and we should be able to decrease them
  real(mytype), save, allocatable, dimension(:,:,:) :: sgsx1,sgsy1,sgsz1,nut1,sxx1,syy1,szz1,sxy1,sxz1,syz1
  real(mytype), save, allocatable, dimension(:,:,:) :: sgsx2,sgsy2,sgsz2,nut2,sxx2,syy2,szz2,sxy2,sxz2,syz2
  real(mytype), save, allocatable, dimension(:,:,:) :: sgsx3,sgsy3,sgsz3,nut3,sxx3,syy3,szz3,sxy3,sxz3,syz3

  real(mytype), save, allocatable, dimension(:,:,:) :: sdxx1,sdyy1,sdzz1,sdxy1,sdxz1,sdyz1
  real(mytype), save, allocatable, dimension(:,:,:) :: sdxx2,sdyy2,sdzz2,sdxy2,sdxz2,sdyz2
  real(mytype), save, allocatable, dimension(:,:,:) :: sdxx3,sdyy3,sdzz3,sdxy3,sdxz3,sdyz3

  real(mytype), save, allocatable, dimension(:,:,:) :: gxx1,gyx1,gzx1,gxy1,gyy1,gzy1,gxz1,gyz1,gzz1
  real(mytype), save, allocatable, dimension(:,:,:) :: gxy2,gyy2,gzy2,gxz2,gyz2,gzz2
  real(mytype), save, allocatable, dimension(:,:,:) :: gxz3,gyz3,gzz3
  real(mytype), save, allocatable, dimension(:,:,:,:) :: sgsphi1,sgsphi2,sgsphi3

  real(mytype), save, allocatable, dimension(:,:,:) :: srt_smag, srt_smag2
  real(mytype), save, allocatable, dimension(:,:,:) :: srt_wale, srt_wale2, srt_wale3, srt_wale4


  ! working arrays for ABL
  real(mytype), save, allocatable, dimension(:,:) :: heatflux
  real(mytype), save, allocatable, dimension(:,:,:,:) :: abl_T
  ! tmp fix for intel
  real(mytype), save, allocatable, dimension(:,:,:) :: wallfluxx1, wallfluxy1, wallfluxz1  

  ! arrays for turbine modelling
  real(mytype), save, allocatable, dimension(:,:,:) :: FTx, FTy, FTz, Fdiscx, Fdiscy, Fdiscz
  real(mytype), save, allocatable, dimension(:,:,:,:) ::  GammaDisc

  ! arrays for inflow/outflow - precursor simulations
  real(mytype), save, allocatable, dimension(:,:,:) :: ux_inflow, uy_inflow, uz_inflow
  real(mytype), save, allocatable, dimension(:,:,:) :: ux_recoutflow, uy_recoutflow, uz_recoutflow

contains

  subroutine init_variables

    implicit none 

    TYPE(DECOMP_INFO), save :: ph! decomposition object

    integer :: i, j, k
    
#ifdef DEBG
    if (nrank == 0) write(*,*) '# Init_variables start'
#endif

    if (nrank == 0) write(*,*) '# Initializing variables...'

    if (nclx) then
       nxmsize = xsize(1)
    else
       nxmsize = xsize(1) -1
    endif
    if (ncly) then
       nymsize = ysize(2)
    else
       nymsize = ysize(2) -1
    endif
    if (nclz) then
       nzmsize = zsize(3)
    else
       nzmsize = zsize(3) -1
    endif
    call decomp_info_init(nxmsize, nymsize, nzmsize, ph)
    !xsize(i), ysize(i), zsize(i), i=1,2,3 - sizes of the sub-domains held by the current process. The first letter refers to the pencil orientation and the three 1D array elements contain the sub-domain sizes in X, Y and Z directions, respectively. In a 2D pencil decomposition, there is always one dimension which completely resides in local memory. So by definition xsize(1)==nx_global, ysize(2)==ny_global and zsize(3)==nz_global.

    !xstart(i), ystart(i), zstart(i), xend(i), yend(i), zend(i), i=1,2,3 - the starting and ending indices for each sub-domain, as in the global coordinate system. Obviously, it can be seen that xsize(i)=xend(i)-xstart(i)+1. It may be convenient for certain applications to use global coordinate (for example when extracting a 2D plane from a 3D domain, it is easier to know which process owns the plane if global index is used).


    !X PENCILS
    call alloc_x(ux1, opt_global=.true.) !global indices
    ux1 = zero
    call alloc_x(uy1, opt_global=.true.) !global indices
    uy1 = zero
    call alloc_x(uz1, opt_global=.true.) !global indices
    uz1 = zero
    call alloc_x(px1, opt_global=.true.) !global indices
    px1 = zero
    call alloc_x(py1, opt_global=.true.) !global indices
    py1 = zero
    call alloc_x(pz1, opt_global=.true.) !global indices
    pz1 = zero
    call alloc_x(diss1, opt_global=.true.) !global indices
    diss1 = zero
    call alloc_x(pre1, opt_global=.true.) !global indices
    pre1 = zero

    allocate(phi1(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),1:numscalar)) !global indices
    phi1 = zero

    call alloc_x(ta1)
    ta1 = zero
    call alloc_x(tb1)
    tb1 = zero
    call alloc_x(tc1)
    tc1 = zero
    call alloc_x(td1)
    td1 = zero
    call alloc_x(te1)
    te1 = zero
    call alloc_x(tf1)
    tf1 = zero
    call alloc_x(tg1)
    tg1 = zero
    call alloc_x(th1)
    th1 = zero
    call alloc_x(ti1)
    ti1 = zero
    call alloc_x(di1)
    di1 = zero
    call alloc_x(ep1)
    ep1 = zero
    if (ilmn) then
      call alloc_x(mu1)
      mu1 = one
    endif

    call alloc_x(uxf1);
    uxf1 = zero
    call alloc_x(uyf1);
    uyf1 = zero
    call alloc_x(uzf1);
    uzf1 = zero
    call alloc_x(phif1);
    phif1 = zero
    allocate(pp1(nxmsize,xsize(2),xsize(3)))
    pp1 = zero
    allocate(pgy1(nxmsize,xsize(2),xsize(3)))
    pgy1 = zero
    allocate(pgz1(nxmsize,xsize(2),xsize(3)))
    pgz1 = zero

    !inflow/ouflow 2d arrays
    allocate(bxx1(xsize(2),xsize(3)))
    bxx1=zero
    allocate(bxy1(xsize(2),xsize(3)))
    bxy1=zero
    allocate(bxz1(xsize(2),xsize(3)))
    bxz1=zero

    allocate(bxxn(xsize(2),xsize(3)))
    bxxn=zero
    allocate(bxyn(xsize(2),xsize(3)))
    bxyn=zero
    allocate(bxzn(xsize(2),xsize(3)))
    bxzn=zero

    allocate(byx1(xsize(1),xsize(3)))
    byx1=zero
    allocate(byy1(xsize(1),xsize(3)))
    byy1=zero
    allocate(byz1(xsize(1),xsize(3)))
    byz1=zero

    allocate(byxn(xsize(1),xsize(3)))
    byxn=zero
    allocate(byyn(xsize(1),xsize(3)))
    byyn=zero
    allocate(byzn(xsize(1),xsize(3)))
    byzn=zero

    allocate(bzx1(xsize(1),xsize(2)))
    bzx1=zero
    allocate(bzy1(xsize(1),xsize(2)))
    bzy1=zero
    allocate(bzz1(xsize(1),xsize(2)))
    bzz1=zero

    allocate(bzxn(xsize(1),xsize(2)))
    bzxn=zero
    allocate(bzyn(xsize(1),xsize(2)))
    bzyn=zero
    allocate(bzzn(xsize(1),xsize(2)))
    bzzn=zero

    allocate(bxo(xsize(2),xsize(3)))
    bxo=zero
    allocate(byo(xsize(2),xsize(3)))
    byo=zero
    allocate(bzo(xsize(2),xsize(3)))
    bzo=zero

    ! BC values in Y-pencil - Ricardo Frantz
    allocate(byx1_2(ysize(1),ysize(3)))
    byx1_2=zero
    allocate(byxn_2(ysize(1),ysize(3)))
    byxn_2=zero
    allocate(byy1_2(ysize(1),ysize(3)))
    byy1_2=zero
    allocate(byyn_2(ysize(1),ysize(3)))
    byyn_2=zero
    allocate(byz1_2(ysize(1),ysize(3)))
    byz1_2=zero
    allocate(byzn_2(ysize(1),ysize(3)))
    byzn_2=zero

    !inflow/outflow arrays (precursor simulations)
    if (iin.eq.3) then
       allocate(ux_inflow(ntimesteps,xsize(2),xsize(3)))
       ux_inflow=zero
       allocate(uy_inflow(ntimesteps,xsize(2),xsize(3)))
       uy_inflow=zero
       allocate(uz_inflow(ntimesteps,xsize(2),xsize(3)))
       uz_inflow=zero
    endif

    if (ioutflow.eq.1) then
       allocate(ux_recoutflow(ntimesteps,xsize(2),xsize(3)))
       ux_recoutflow=zero
       allocate(uy_recoutflow(ntimesteps,xsize(2),xsize(3)))
       uy_recoutflow=zero
       allocate(uz_recoutflow(ntimesteps,xsize(2),xsize(3)))
       uz_recoutflow=zero
    endif

    !pre_correc 2d array
    allocate(dpdyx1(xsize(2),xsize(3)),dpdyxn(xsize(2),xsize(3)))
    dpdyx1=zero
    dpdyxn=zero
    allocate(dpdzx1(xsize(2),xsize(3)),dpdzxn(xsize(2),xsize(3)))
    dpdzx1=zero
    dpdzxn=zero
    allocate(dpdxy1(xsize(1),xsize(3)),dpdxyn(xsize(1),xsize(3)))
    dpdxy1=zero
    dpdxyn=zero
    allocate(dpdzy1(xsize(1),xsize(3)),dpdzyn(xsize(1),xsize(3)))
    dpdzy1=zero
    dpdzyn=zero
    allocate(dpdxz1(xsize(1),xsize(2)),dpdxzn(xsize(1),xsize(2)))
    dpdxz1=zero
    dpdxzn=zero
    allocate(dpdyz1(xsize(1),xsize(2)),dpdyzn(xsize(1),xsize(2)))
    dpdyz1=zero
    dpdyzn=zero

    !arrays for visualization!pay attention to the size!
    allocate(uvisu(xstV(1):xenV(1),xstV(2):xenV(2),xstV(3):xenV(3)))
    uvisu=zero

    !arrays statistics
    allocate (umean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    umean=zero
    allocate (vmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    vmean=zero
    allocate (wmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    wmean=zero
    allocate (pmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    pmean=zero

    allocate (uumean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    uumean=zero
    allocate (vvmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    vvmean=zero
    allocate (wwmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    wwmean=zero
    allocate (uvmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    uvmean=zero
    allocate (uwmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    uwmean=zero
    allocate (vwmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    vwmean=zero

    allocate (phimean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3),numscalar))
    phimean=zero
    allocate (phiphimean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3),numscalar))
    phiphimean=zero

    allocate (tmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    tmean=zero

    !Y PENCILS
    call alloc_y(ux2)
    ux2=zero
    call alloc_y(uy2)
    uy2=zero
    call alloc_y(uz2)
    uz2=zero
    call alloc_y(ta2)
    ta2=zero
    call alloc_y(tb2)
    tb2=zero
    call alloc_y(tc2)
    tc2=zero
    call alloc_y(td2)
    td2=zero
    call alloc_y(te2)
    te2=zero
    call alloc_y(tf2)
    tf2=zero
    call alloc_y(tg2)
    tg2=zero
    call alloc_y(th2)
    th2=zero
    call alloc_y(ti2)
    ti2=zero
    call alloc_y(tj2)
    tj2=zero
    call alloc_y(di2)
    di2=zero
    call alloc_y(uxf2)
    uxf2=zero
    call alloc_y(uyf2)
    uyf2=zero
    call alloc_y(uzf2)
    uzf2=zero
    call alloc_y(phif2)
    phif2=zero
    allocate(phi2(ysize(1),ysize(2),ysize(3),1:numscalar))
    phi2=zero
    allocate(pgz2(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)))
    pgz2=zero
    allocate(pp2(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)))
    pp2=zero
    allocate(dip2(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)))
    dip2=zero
    allocate(ppi2(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)))
    ppi2=zero
    allocate(pgy2(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)))
    pgy2=zero
    allocate(pgzi2(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)))
    pgzi2=zero
    allocate(duxdxp2(ph1%yst(1):ph1%yen(1),ysize(2),ysize(3)))
    duxdxp2=zero
    allocate(uyp2(ph1%yst(1):ph1%yen(1),ysize(2),ysize(3)))
    uyp2=zero
    allocate(uzp2(ph1%yst(1):ph1%yen(1),ysize(2),ysize(3)))
    uzp2=zero
    allocate(dipp2(ph1%yst(1):ph1%yen(1),ysize(2),ysize(3)))
    dipp2=zero
    allocate(upi2(ph1%yst(1):ph1%yen(1),nymsize,ysize(3)))
    upi2=zero
    allocate(duydypi2(ph1%yst(1):ph1%yen(1),nymsize,ysize(3)))
    duydypi2=zero

    if (ilmn) then
      call alloc_y(mu2)
      mu2=one
    endif

    !Z PENCILS
    call alloc_z(ux3)
    ux3=zero
    call alloc_z(uy3)
    uy3=zero
    call alloc_z(uz3)
    uz3=zero
    call alloc_z(ta3)
    ta3=zero
    call alloc_z(tb3)
    tb3=zero
    call alloc_z(tc3)
    tc3=zero
    call alloc_z(td3)
    td3=zero
    call alloc_z(te3)
    te3=zero
    call alloc_z(tf3)
    tf3=zero
    call alloc_z(tg3)
    tg3=zero
    call alloc_z(th3)
    th3=zero
    call alloc_z(ti3)
    ti3=zero
    call alloc_z(di3)
    di3=zero
    call alloc_z(uxf3)
    uxf3=zero
    call alloc_z(uyf3)
    uyf3=zero
    call alloc_z(uzf3) 
    uzf3=zero
    call alloc_z(phif3)
    phif3=zero
    allocate(phi3(zsize(1),zsize(2),zsize(3),1:numscalar))
    phi3=zero
    allocate(pgz3(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)))
    pgz3=zero
    allocate(ppi3(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)))
    ppi3=zero
    allocate(dip3(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)))
    dip3=zero

    allocate(duxydxyp3(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),zsize(3)))
    duxydxyp3=zero
    allocate(uzp3(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),zsize(3)))
    uzp3=zero
    allocate(dipp3(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),zsize(3)))
    dipp3=zero

    if (ilmn) then
      call alloc_z(mu3)
      mu3(:,:,:) = one
    endif

    allocate(pp3(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress))
    pp3=zero

    call alloc_z(dv3,ph,.true.)
    dv3=zero
    call alloc_z(po3,ph,.true.)
    po3=zero

    if(ilesmod /= 0.and.jles > 0) then
       call alloc_x(sgsx1)
       sgsx1=zero
       call alloc_x(sgsy1)
       sgsy1=zero
       call alloc_x(sgsz1)
       sgsz1=zero
       call alloc_x(sxx1)
       sxx1=zero
       call alloc_x(syy1)
       syy1=zero
       call alloc_x(szz1)
       szz1=zero
       call alloc_x(sxy1)
       sxy1=zero
       call alloc_x(sxz1)
       sxz1=zero
       call alloc_x(syz1)
       syz1=zero
       call alloc_x(sdxx1)
       sdxx1=zero
       call alloc_x(sdyy1)
       sdyy1=zero
       call alloc_x(sdzz1)
       sdzz1=zero
       call alloc_x(sdxy1)
       sdxy1=zero
       call alloc_x(sdxz1)
       sdxz1=zero
       call alloc_x(sdyz1)
       sdyz1=zero
       call alloc_x(nut1)
       nut1=zero
       call alloc_x(srt_smag)
       srt_smag=zero
       call alloc_x(srt_wale)
       srt_wale=zero
       call alloc_x(wallfluxx1)
       wallfluxx1=zero
       call alloc_x(wallfluxy1)
       wallfluxy1=zero
       call alloc_x(wallfluxz1)
       wallfluxz1=zero
       call alloc_y(sgsx2)
       sgsx2=zero
       call alloc_y(sgsy2)
       sgsy2=zero
       call alloc_y(sgsz2)
       sgsz2=zero
       call alloc_y(sxx2)
       sxx2=zero
       call alloc_y(syy2)
       syy2=zero
       call alloc_y(szz2)
       szz2=zero
       call alloc_y(sxy2)
       sxy2=zero
       call alloc_y(sxz2)
       sxz2=zero
       call alloc_y(syz2)
       syz2=zero
       call alloc_y(sdxx2) 
       sdxx2=zero
       call alloc_y(sdyy2)
       sdyy2=zero
       call alloc_y(sdzz2)
       sdzz2=zero
       call alloc_y(sdxy2) 
       sdxy2=zero
       call alloc_y(sdxz2)
       sdxz2=zero
       call alloc_y(sdyz2)
       sdyz2=zero
       call alloc_y(nut2)
       nut2=zero
       call alloc_y(srt_smag2)
       srt_smag2=zero
       call alloc_y(srt_wale2)
       srt_wale2=zero
       call alloc_y(srt_wale3)
       srt_wale3=zero
       call alloc_y(srt_wale4)
       srt_wale4=zero
       call alloc_z(sgsx3)
       sgsx3=zero
       call alloc_z(sgsy3)
       sgsy3=zero
       call alloc_z(sgsz3)
       sgsz3=zero
       call alloc_z(sxx3)
       sxx3=zero
       call alloc_z(syy3)
       syy3=zero
       call alloc_z(szz3)
       szz3=zero
       call alloc_z(sxy3)
       sxy3=zero
       call alloc_z(sxz3)
       sxz3=zero
       call alloc_z(syz3)
       syz3=zero
       call alloc_z(sdxx3)
       sdxx3=zero
       call alloc_z(sdyy3)
       sdyy3=zero
       call alloc_z(sdzz3)
       sdzz3=zero
       call alloc_z(sdxy3)
       sdxy3=zero
       call alloc_z(sdxz3)
       sdxz3=zero
       call alloc_z(sdyz3)
       sdyz3=zero
       call alloc_z(nut3)
       nut3=zero
       call alloc_x(gxx1)
       gxx1=zero
       call alloc_x(gyx1)
       gyx1=zero
       call alloc_x(gzx1)
       gzx1=zero
       call alloc_x(gxy1)
       gxy1=zero
       call alloc_x(gyy1)
       gyy1=zero
       call alloc_x(gzy1)
       gzy1=zero
       call alloc_x(gxz1)
       gxz1=zero
       call alloc_x(gyz1)
       gyz1=zero
       call alloc_x(gzz1)
       gzz1=zero
       call alloc_y(gxy2)
       gxy2=zero
       call alloc_y(gyy2)
       gyy2=zero
       call alloc_y(gzy2)
       gzy2=zero
       call alloc_y(gxz2)
       gxz2=zero
       call alloc_y(gyz2)
       gyz2=zero
       call alloc_y(gzz2)
       gzz2=zero
       call alloc_z(gxz3)
       gxz3=zero
       call alloc_z(gyz3)
       gyz3=zero
       call alloc_z(gzz3)
       gzz3=zero

       allocate(sgsphi1(xsize(1),xsize(2),xsize(3),1:numscalar))
       sgsphi1=zero
       allocate(sgsphi2(ysize(1),ysize(2),ysize(3),1:numscalar))
       sgsphi2=zero
       allocate(sgsphi3(zsize(1),zsize(2),zsize(3),1:numscalar))
       sgsphi3=zero
    endif


    if (iibm.ne.0) then
       !complex_geometry
       nxraf=(nxm)*nraf
       if (.not.nclx) nxraf=nxraf+1
       nyraf=(nym)*nraf
       if (.not.ncly) nyraf=nyraf+1
       nzraf=(nzm)*nraf
       if (.not.nclz) nzraf=nzraf+1
       allocate(nobjx(xsize(2),xsize(3)))
       nobjx=zero
       allocate(nobjy(ysize(1),ysize(3)))
       nobjy=zero
       allocate(nobjz(zsize(1),zsize(2)))
       nobjz=zero
       allocate(nxipif(0:nobjmax,xsize(2),xsize(3)))
       nxipif=zero
       allocate(nxfpif(0:nobjmax,xsize(2),xsize(3)))
       nxfpif=zero
       allocate(nyipif(0:nobjmax,ysize(1),ysize(3)))
       nyipif=zero
       allocate(nyfpif(0:nobjmax,ysize(1),ysize(3)))
       nyfpif=zero
       allocate(nzipif(0:nobjmax,zsize(1),zsize(2)))
       nzipif=zero
       allocate(nzfpif(0:nobjmax,zsize(1),zsize(2)))
       nzfpif=zero
       allocate(xi(nobjmax,xsize(2),xsize(3)))
       xi=zero
       allocate(xf(nobjmax,xsize(2),xsize(3)))
       xf=zero
       allocate(yi(nobjmax,ysize(1),ysize(3)))
       yi=zero
       allocate(yf(nobjmax,ysize(1),ysize(3)))
       yf=zero
       allocate(zi(nobjmax,zsize(1),zsize(2)))
       zi=zero
       allocate(zf(nobjmax,zsize(1),zsize(2)))
       zf=zero
       allocate(xepsi(nxraf,xsize(2),xsize(3))) 
       xepsi=zero
       allocate(yepsi(ysize(1),nyraf,ysize(3)))
       yepsi=zero
       allocate(zepsi(zsize(1),zsize(2),nzraf)) 
       zepsi=zero

    endif

    !module filter
    allocate(fiffx(nx))
    fiffx=zero
    allocate(fisfx(nx))
    fisfx=zero
    allocate(fifsx(nx))
    fifsx=zero
    allocate(fifwx(nx))
    fifwx=zero
    allocate(fissx(nx))
    fissx=zero
    allocate(fiswx(nx))
    fiswx=zero

    allocate(fiffxp(nx))
    fiffxp=zero
    allocate(fisfxp(nx))
    fisfxp=zero
    allocate(fifsxp(nx))
    fifsxp=zero
    allocate(fifwxp(nx))
    fifwxp=zero
    allocate(fissxp(nx))
    fissxp=zero
    allocate(fiswxp(nx))
    fiswxp=zero

    allocate(fiffy(ny))
    fiffy=zero
    allocate(fisfy(ny))
    fisfy=zero
    allocate(fifsy(ny))
    fifsy=zero
    allocate(fifwy(ny))
    fifwy=zero
    allocate(fissy(ny))
    fissy=zero
    allocate(fiswy(ny))
    fiswy=zero

    allocate(fiffyp(ny))
    fiffyp=zero
    allocate(fisfyp(ny))
    fisfyp=zero
    allocate(fifsyp(ny))
    fifsyp=zero
    allocate(fifwyp(ny))
    fifwyp=zero
    allocate(fissyp(ny))
    fissyp=zero
    allocate(fiswyp(ny))
    fiswyp=zero

    allocate(fiffz(nz))
    fiffz=zero
    allocate(fisfz(nz))
    fisfz=zero
    allocate(fifsz(nz))
    fifsz=zero
    allocate(fifwz(nz))
    fifwz=zero
    allocate(fissz(nz))
    fissz=zero
    allocate(fiswz(nz))
    fiswz=zero

    allocate(fiffzp(nz))
    fiffzp=zero
    allocate(fisfzp(nz))
    fisfzp=zero
    allocate(fifszp(nz))
    fifszp=zero
    allocate(fifwzp(nz))
    fifwzp=zero
    allocate(fisszp(nz))
    fisszp=zero
    allocate(fiswzp(nz))
    fiswzp=zero

    allocate(fisx(xsize(2),xsize(3)))
    fisx=zero
    allocate(fivx(xsize(2),xsize(3)))
    fivx=zero

    allocate(fisy(ysize(1),ysize(3)))
    fisy=zero
    allocate(fivy(ysize(1),ysize(3)))
    fivy=zero

    allocate(fisz(zsize(1),zsize(2)))
    fisz=zero
    allocate(fivz(zsize(1),zsize(2)))
    fivz=zero

    !module derivative
    allocate(ffx(nx))
    ffx=zero
    allocate(sfx(nx))
    sfx=zero
    allocate(fsx(nx))
    fsx=zero
    allocate(fwx(nx))
    fwx=zero
    allocate(ssx(nx))
    ssx=zero
    allocate(swx(nx))
    swx=zero

    allocate(ffxp(nx))
    ffxp=zero
    allocate(sfxp(nx))
    sfxp=zero
    allocate(fsxp(nx))
    fsxp=zero
    allocate(fwxp(nx))
    fwxp=zero
    allocate(ssxp(nx))
    ssxp=zero
    allocate(swxp(nx))
    swxp=zero

    allocate(ffy(ny))
    ffy=zero
    allocate(sfy(ny))
    sfy=zero
    allocate(fsy(ny))
    fsy=zero
    allocate(fwy(ny))
    fwy=zero
    allocate(ssy(ny))
    ssy=zero
    allocate(swy(ny))
    swy=zero

    allocate(ffyp(ny))
    ffyp=zero
    allocate(sfyp(ny))
    sfyp=zero
    allocate(fsyp(ny))
    fsyp=zero
    allocate(fwyp(ny))
    fwyp=zero
    allocate(ssyp(ny))
    ssyp=zero
    allocate(swyp(ny))
    swyp=zero

    allocate(ffz(nz))
    ffz=zero
    allocate(sfz(nz))
    sfz=zero
    allocate(fsz(nz))
    fsz=zero
    allocate(fwz(nz))
    fwz=zero
    allocate(ssz(nz))
    ssz=zero
    allocate(swz(nz))
    swz=zero

    allocate(ffzp(nz))
    ffzp=zero
    allocate(sfzp(nz))
    sfzp=zero
    allocate(fszp(nz))
    fszp=zero
    allocate(fwzp(nz))
    fwzp=zero
    allocate(sszp(nz))
    sszp=zero
    allocate(swzp(nz))
    swzp=zero

    allocate(ffxS(nx))
    ffxS=zero
    allocate(sfxS(nx))
    sfxS=zero
    allocate(fsxS(nx))
    fsxS=zero
    allocate(fwxS(nx))
    fwxS=zero
    allocate(ssxS(nx))
    ssxS=zero
    allocate(swxS(nx))
    swxS=zero

    allocate(ffxpS(nx))
    ffxpS=zero
    allocate(sfxpS(nx))
    sfxpS=zero
    allocate(fsxpS(nx))
    fsxpS=zero
    allocate(fwxpS(nx))
    fwxpS=zero
    allocate(ssxpS(nx))
    ssxpS=zero
    allocate(swxpS(nx))
    swxpS=zero

    allocate(ffyS(ny))
    ffyS=zero
    allocate(sfyS(ny))
    sfyS=zero
    allocate(fsyS(ny))
    fsyS=zero
    allocate(fwyS(ny))
    fwyS=zero
    allocate(ssyS(ny))
    ssyS=zero
    allocate(swyS(ny))
    swyS=zero

    allocate(ffypS(ny))
    ffypS=zero
    allocate(sfypS(ny))
    sfypS=zero
    allocate(fsypS(ny))
    fsypS=zero
    allocate(fwypS(ny))
    fwypS=zero
    allocate(ssypS(ny))
    ssypS=zero
    allocate(swypS(ny))
    swypS=zero

    allocate(ffzS(nz))
    ffzS=zero
    allocate(sfzS(nz))
    sfzS=zero
    allocate(fszS(nz))
    fszS=zero
    allocate(fwzS(nz))
    fwzS=zero
    allocate(sszS(nz))
    sszS=zero
    allocate(swzS(nz))
    swzS=zero

    allocate(ffzpS(nz))
    ffzpS=zero
    allocate(sfzpS(nz))
    sfzpS=zero
    allocate(fszpS(nz))
    fszpS=zero
    allocate(fwzpS(nz))
    fwzpS=zero
    allocate(sszpS(nz))
    sszpS=zero
    allocate(swzpS(nz))
    swzpS=zero

    !MHD -- todo if active
    allocate(ffxB(nx,3))
    ffxB=zero
    allocate(sfxB(nx,3))
    sfxB=zero
    allocate(fsxB(nx,3))
    fsxB=zero
    allocate(fwxB(nx,3))
    fwxB=zero
    allocate(ssxB(nx,3))
    ssxB=zero
    allocate(swxB(nx,3))
    swxB=zero

    allocate(ffxpB(nx,3))
    ffxpB=zero
    allocate(sfxpB(nx,3))
    sfxpB=zero
    allocate(fsxpB(nx,3))
    fsxpB=zero
    allocate(fwxpB(nx,3))
    fwxpB=zero
    allocate(ssxpB(nx,3))
    ssxpB=zero
    allocate(swxpB(nx,3))
    swxpB=zero

    allocate(ffyB(ny,3))
    ffyB=zero
    allocate(sfyB(ny,3))
    sfyB=zero
    allocate(fsyB(ny,3))
    fsyB=zero
    allocate(fwyB(ny,3))
    fwyB=zero
    allocate(ssyB(ny,3))
    ssyB=zero
    allocate(swyB(ny,3))
    swyB=zero

    allocate(ffypB(ny,3))
    ffypB=zero
    allocate(sfypB(ny,3))
    sfypB=zero
    allocate(fsypB(ny,3))
    fsypB=zero
    allocate(fwypB(ny,3))
    fwypB=zero
    allocate(ssypB(ny,3))
    ssypB=zero
    allocate(swypB(ny,3))
    swypB=zero

    allocate(ffzB(nz,3))
    ffzB=zero
    allocate(sfzB(nz,3))
    sfzB=zero
    allocate(fszB(nz,3))
    fszB=zero
    allocate(fwzB(nz,3))
    fwzB=zero
    allocate(sszB(nz,3))
    sszB=zero
    allocate(swzB(nz,3))
    swzB=zero

    allocate(ffzpB(nz,3))
    ffzpB=zero
    allocate(sfzpB(nz,3))
    sfzpB=zero
    allocate(fszpB(nz,3))
    fszpB=zero
    allocate(fwzpB(nz,3))
    fwzpB=zero
    allocate(sszpB(nz,3))
    sszpB=zero
    allocate(swzpB(nz,3))
    swzpB=zero


    allocate(sx(xsize(2),xsize(3)))
    sx=zero
    allocate(vx(xsize(2),xsize(3)))
    vx=zero

    allocate(sy(ysize(1),ysize(3)))
    sy=zero
    allocate(vy(ysize(1),ysize(3)))
    vy=zero

    allocate(sz(zsize(1),zsize(2)))
    sz=zero
    allocate(vz(zsize(1),zsize(2)))
    vz=zero

    !O6SVV
    allocate(newsm(ny))
    newsm=zero
    allocate(newtm(ny))
    newtm=zero
    allocate(newsmt(ny))
    newsmt=zero
    allocate(newtmt(ny))
    newtmt=zero

    allocate(newrm(ny))
    newrm=zero
    allocate(newrmt(ny))
    newrmt=zero

    !module derpres
    allocate(cfx6(nxm))
    cfx6=zero
    allocate(ccx6(nxm))
    ccx6=zero
    allocate(cbx6(nxm))
    cbx6=zero
    allocate(cfxp6(nxm))
    cfxp6=zero
    allocate(ciwxp6(nxm))
    ciwxp6=zero
    allocate(csxp6(nxm))
    csxp6=zero
    allocate(cwxp6(nxm))
    cwxp6=zero
    allocate(csx6(nxm))
    csx6=zero
    allocate(cwx6(nxm))
    cwx6=zero
    allocate(cifx6(nxm))
    cifx6=zero
    allocate(cicx6(nxm))
    cicx6=zero
    allocate(cisx6(nxm))
    cisx6=zero

    allocate(cibx6(nxm))
    cibx6=zero
    allocate(cifxp6(nxm))
    cifxp6=zero
    allocate(cisxp6(nxm))
    cisxp6=zero
    allocate(ciwx6(nxm))
    ciwx6=zero

    allocate(cfi6(nx))
    cfi6=zero
    allocate(cci6(nx))
    cci6=zero
    allocate(cbi6(nx))
    cbi6=zero
    allocate(cfip6(nx))
    cfip6=zero
    allocate(csip6(nx))
    csip6=zero
    allocate(cwip6(nx))
    cwip6=zero
    allocate(csi6(nx))
    csi6=zero
    allocate(cwi6(nx))
    cwi6=zero
    allocate(cifi6(nx))
    cifi6=zero
    allocate(cici6(nx))
    cici6=zero
    allocate(cibi6(nx))
    cibi6=zero
    allocate(cifip6(nx))
    cifip6=zero

    allocate(cisip6(nx))
    cisip6=zero
    allocate(ciwip6(nx))
    ciwip6=zero
    allocate(cisi6(nx))
    cisi6=zero
    allocate(ciwi6(nx))
    ciwi6=zero

    allocate(cfy6(nym))
    cfy6=zero
    allocate(ccy6(nym))
    ccy6=zero
    allocate(cby6(nym))
    cby6=zero
    allocate(cfyp6(nym))
    cfyp6=zero
    allocate(csyp6(nym))
    csyp6=zero
    allocate(cwyp6(nym))
    cwyp6=zero
    allocate(csy6(nym))
    csy6=zero

    allocate(cwy6(nym))
    cwy6=zero
    allocate(cify6(nym))
    cify6=zero
    allocate(cicy6(nym))
    cicy6=zero
    allocate(ciby6(nym))
    ciby6=zero
    allocate(cifyp6(nym))
    cifyp6=zero
    allocate(cisyp6(nym))
    cisyp6=zero
    allocate(ciwyp6(nym))
    ciwyp6=zero
    allocate(cisy6(nym))
    cisy6=zero
    allocate(ciwy6(nym))
    ciwy6=zero

    allocate(cfi6y(ny))
    cfi6y=zero
    allocate(cci6y(ny))
    cci6y=zero
    allocate(cbi6y(ny))
    cbi6y=zero
    allocate(cfip6y(ny))
    cfip6y=zero
    allocate(csip6y(ny))
    csip6y=zero
    allocate(cwip6y(ny))
    cwip6y=zero
    allocate(csi6y(ny))
    csi6y=zero
    allocate(cwi6y(ny))
    cwi6y=zero
    allocate(cifi6y(ny))
    cifi6y=zero
    allocate(cici6y(ny))
    cici6y=zero

    allocate(cibi6y(ny))
    cibi6y=zero
    allocate(cifip6y(ny))
    cifip6y=zero
    allocate(cisip6y(ny))
    cisip6y=zero
    allocate(ciwip6y(ny))
    ciwip6y=zero
    allocate(cisi6y(ny))
    cisi6y=zero
    allocate(ciwi6y(ny))
    ciwi6y=zero

    allocate(cfz6(nzm))
    cfz6=zero
    allocate(ccz6(nzm))
    ccz6=zero
    allocate(cbz6(nzm))
    cbz6=zero
    allocate(cfzp6(nzm))
    cfzp6=zero
    allocate(cszp6(nzm))
    cszp6=zero
    allocate(cwzp6(nzm))
    cwzp6=zero
    allocate(csz6(nzm))
    csz6=zero

    allocate(cwz6(nzm))
    cwz6=zero
    allocate(cifz6(nzm))
    cifz6=zero
    allocate(cicz6(nzm))
    cicz6=zero
    allocate(cibz6(nzm))
    cibz6=zero
    allocate(cifzp6(nzm))
    cifzp6=zero
    allocate(ciszp6(nzm))
    ciszp6=zero
    allocate(ciwzp6(nzm))
    ciwzp6=zero
    allocate(cisz6(nzm))
    cisz6=zero
    allocate(ciwz6(nzm))
    ciwz6=zero

    allocate(cfi6z(nz))
    cfi6z=zero
    allocate(cci6z(nz))
    cci6z=zero
    allocate(cbi6z(nz))
    cbi6z=zero
    allocate(cfip6z(nz))
    cfip6z=zero
    allocate(csip6z(nz))
    csip6z=zero
    allocate(cwip6z(nz))
    cwip6z=zero
    allocate(csi6z(nz))
    csi6z=zero
    allocate(cwi6z(nz))
    cwi6z=zero
    allocate(cifi6z(nz))
    cifi6z=zero
    allocate(cici6z(nz))
    cici6z=zero

    allocate(cibi6z(nz))
    cibi6z=zero
    allocate(cifip6z(nz))
    cifip6z=zero
    allocate(cisip6z(nz))
    cisip6z=zero
    allocate(ciwip6z(nz))
    ciwip6z=zero
    allocate(cisi6z(nz))
    cisi6z=zero
    allocate(ciwi6z(nz))
    ciwi6z=zero

    !module waves
    allocate(zkz(nz/2+1))
    zkz=zero
    allocate(zk2(nz/2+1))
    zk2=zero
    allocate(ezs(nz/2+1))
    ezs=zero

    allocate(yky(ny))
    yky=zero
    allocate(yk2(ny))
    yk2=zero
    allocate(eys(ny))
    eys=zero

    allocate(xkx(nx))
    xkx=zero
    allocate(xk2(nx))
    xk2=zero
    allocate(exs(nx))
    exs=zero

    !module mesh
    allocate(ppy(ny))
    ppy=zero
    allocate(pp2y(ny))
    pp2y=zero
    allocate(pp4y(ny))
    pp4y=zero

    allocate(ppyi(ny))
    ppyi=zero
    allocate(pp2yi(ny))
    pp2yi=zero
    allocate(pp4yi(ny))
    pp4yi=zero

    allocate(xp(nx))
    xp=zero
    allocate(xpi(nx))
    xpi=zero

    allocate(yp(ny))
    yp=zero
    allocate(ypi(ny))
    ypi=zero
    allocate(del(ny))
    del=zero

    allocate(zp(nz))
    zp=zero
    allocate(zpi(nz))
    zpi=zero

    ! x-position
    do i=1,nx
      xp(i)=real(i-1,mytype)*dx
      xpi(i)=(real(i,mytype)-half)*dx
    enddo
    ! y-position
    if (istret.eq.0) then
       do j=1,ny
          yp(j)=real(j-1,mytype)*dy
          ypi(j)=(real(j,mytype)-half)*dy
          ppy(j) = one
       enddo
       if (ncly1.eq.1 .or. ncly1.eq.2) then
          ppy(1) = two
       endif
       if (nclyn.eq.1 .or. nclyn.eq.2) then
          ppy(ny) = two
       endif
    else
       call stretching(ny, &
                       yp, ypi, &
                       ppy, pp2y, pp4y, &
                       ppyi, pp2yi, pp4yi, &
                       opt_write = .true.)

       ! compute integral weights for stretched mesh - Ricardo Frantz
       !TODO: change for exact formula
       allocate(ypw(ny))
       ypw=zero
       do j=1,ny
          if    (j==1)then
             ypw(j) = (yp(j+1)-yp(j))*half
          elseif(j==ny)then
             ypw(j) = (yp(j)-yp(j-1))*half
          else
             ypw(j) = (yp(j+1)-yp(j-1))*half
          endif
       enddo
       allocate(dyp(ny))
       ! compute dy for stretched mesh - Kay
       do j=2,ny-1
          dyp(j) = half*(yp(j+1)-yp(j-1))
       enddo
       dyp(1)  = yp(2) -yp(1)
       dyp(ny) = yp(ny)-yp(ny-1)
    endif
    ! z-position
    do k=1,nz
       zp(k)=real(k-1,mytype)*dz
       zpi(k)=(real(k,mytype)-half)*dz
    enddo
    !
    adt=zero
    bdt=zero
    cdt=zero
    gdt=zero

    if (itimescheme.eq.1) then ! Euler

       iadvance_time=1

       adt(1)=one*dt
       bdt(1)=zero
       gdt(1)=adt(1)+bdt(1)
       gdt(3)=gdt(1)

       ntime = 1
       nrhotime = 2
    elseif (itimescheme.eq.2) then ! AB2
       iadvance_time=1
       adt(1)=onepfive*dt
       bdt(1)=-half*dt
       gdt(1)=adt(1)+bdt(1)
       gdt(3)=gdt(1)

       ntime = 2
       nrhotime = 3
    elseif (itimescheme.eq.3) then ! AB3
       iadvance_time=1

       adt(1)= (twentythree/twelve)*dt
       bdt(1)=-(    sixteen/twelve)*dt
       cdt(1)= (       five/twelve)*dt
       gdt(1)=adt(1)+bdt(1)+cdt(1)
       gdt(3)=gdt(1)

       ntime = 3
       nrhotime = 4
    elseif(itimescheme==4) then  ! AB4
       iadvance_time=1

       adt(1)= (  fiftyfive/twentyfour)*dt
       bdt(1)=-(  fiftynine/twentyfour)*dt
       cdt(1)= (thirtyseven/twentyfour)*dt
       ddt(1)=-(       nine/twentyfour)*dt
       gdt(1)=adt(1)+bdt(1)+cdt(1)+ddt(1)
       gdt(3)=gdt(1)

       ntime    = 4
       nrhotime = 5
    elseif(itimescheme.eq.5) then !RK3
       iadvance_time=3

       adt(1)=(eight/fifteen)*dt
       bdt(1)= zero
       gdt(1)=adt(1)
       adt(2)=(      five/twelve)*dt
       bdt(2)=(-seventeen/ sixty)*dt
       gdt(2)=adt(2)+bdt(2)
       adt(3)=( three/four)*dt
       bdt(3)=(-five/twelve)*dt
       gdt(3)=adt(3)+bdt(3)

       ntime = 2
       nrhotime = 3
    elseif(itimescheme.eq.6) then !RK4 Carpenter and Kennedy
       iadvance_time=5
       adt(1)=zero
       adt(2)=-0.4178904745_mytype
       adt(3)=-1.192151694643_mytype
       adt(4)=-1.697784692471_mytype
       adt(5)=-1.514183444257_mytype
       bdt(1)=0.1496590219993_mytype
       bdt(2)=0.3792103129999_mytype
       bdt(3)=0.8229550293869_mytype
       bdt(4)=0.6994504559488_mytype
       bdt(5)=0.1530572479681_mytype
       gdt(1)=0.1496590219993_mytype*dt
       gdt(2)=0.220741935365_mytype*dt
       gdt(3)=0.25185480577_mytype*dt
       gdt(4)=0.33602636754_mytype*dt
       gdt(5)=0.041717869325_mytype*dt

       ntime = 2
       nrhotime = 5 ! (A guess)

    endif
    allocate(dux1(xsize(1),xsize(2),xsize(3),ntime))
    dux1=zero
    allocate(duy1(xsize(1),xsize(2),xsize(3),ntime))
    duy1=zero
    allocate(duz1(xsize(1),xsize(2),xsize(3),ntime))
    duz1=zero

    !! Scalar
    allocate(dphi1(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),ntime,1:numscalar)) !global indices
    dphi1=zero
    
    !! ABL
    allocate(heatflux(xsize(1),xsize(3)))
    heatflux = zero
    allocate(PsiM(xsize(1),xsize(3)))
    PsiM = zero
    allocate(PsiH(xsize(1),xsize(3)))
    PsiH = zero
    allocate(Tstat(xsize(2),1))
    Tstat = zero
    if (itype.eq.itype_abl.and.ibuoyancy.eq.1) then
       allocate(abl_T(xsize(1),xsize(2),xsize(3),numscalar))
       abl_T = zero
    endif

    !! Turbine Modelling
    if (iturbine.eq.1) then
       allocate(FTx(xsize(1),xsize(2),xsize(3)))
       FTx = zero
       allocate(FTy(xsize(1),xsize(2),xsize(3)))
       FTy = zero
       allocate(FTz(xsize(1),xsize(2),xsize(3)))
       FTz = zero
    else if (iturbine.eq.2) then
       allocate(Fdiscx(xsize(1),xsize(2),xsize(3)))
       Fdiscx = zero
       allocate(Fdiscy(xsize(1),xsize(2),xsize(3)))
       Fdiscy = zero
       allocate(Fdiscz(xsize(1),xsize(2),xsize(3)))
       Fdiscz = zero
       allocate(Gammadisc(xsize(1),xsize(2),xsize(3),Ndiscs))
       Gammadisc = zero
    endif

    !! LMN
    if (.not.ilmn) then
       nrhotime = 1 !! Save some space
    endif
    allocate(rho1(xsize(1),xsize(2),xsize(3),nrhotime)) !Need to store old density values to extrapolate drhodt
    rho1=one
    call alloc_y(rho2)
    rho2=zero
    call alloc_z(rho3)
    rho3=zero
    allocate(drho1(xsize(1),xsize(2),xsize(3),ntime))
    drho1=zero

    call alloc_z(divu3, opt_global=.true.) !global indices
    divu3=zero

    ! !TRIPPING
    zs_param=1.7_mytype
    randomseed=4600._mytype
    zs_tr=zs_param/2.853_mytype
    z_modes=zlz/zs_tr

    allocate(h_coeff1(z_modes))
    h_coeff1=zero
    allocate(h_coeff2(z_modes))
    h_coeff2=zero
    allocate(phase1(z_modes))
    phase1=zero
    allocate(phase2(z_modes))
    phase2=zero
    allocate(h_1(xsize(3)))
    h_1=zero
    allocate(h_2(xsize(3)))
    h_2=zero

#ifdef DEBG
    if (nrank ==  0) write(*,*) '# init_variables done'
#endif
    return
  end subroutine init_variables
end module var
