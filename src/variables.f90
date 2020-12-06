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

module var

  use decomp_2d
  USE variables
  USE param
  USE complex_geometry

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
  real(mytype), save, allocatable, dimension(:,:,:) :: ep1, diss1, pre1, depo, depof, kine
  real(mytype), save, allocatable, dimension(:,:,:,:) :: dux1,duy1,duz1  ! Output of convdiff
  real(mytype), save, allocatable, dimension(:,:,:,:,:) :: dphi1
  real(mytype), save, allocatable, dimension(:,:,:) :: mu1

  !arrays for post processing
  real(mytype), save, allocatable, dimension(:,:,:) :: f1,fm1
  real(mytype), save, allocatable, dimension(:,:,:) :: uxm1, uym1, phim1, prem1, dissm1
  real(mytype), save, allocatable, dimension(:,:,:) :: uxm2, uym2, phim2, prem2, dissm2

  !arrays for statistic collection
  real(mytype), save, allocatable, dimension(:,:,:) :: umean,vmean,wmean,pmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
  real(mytype), save, allocatable, dimension(:,:,:,:) :: phimean,phiphimean,uphimean,vphimean,wphimean
  real(mytype), save, allocatable, dimension(:,:,:) :: tik1,tik2,tak1,tak2
  real(mytype), save, allocatable, dimension(:,:,:) :: u1sum_tik,u1sum_tak
  real(mytype), save, allocatable, dimension(:,:,:) :: u1sum,v1sum,w1sum,u2sum,v2sum,w2sum
  real(mytype), save, allocatable, dimension(:,:,:) :: u3sum,v3sum,w3sum,u4sum,v4sum,w4sum
  real(mytype), save, allocatable, dimension(:,:,:) :: uvsum,uwsum,vwsum,disssum,presum,tsum

  !arrays for extra statistics collection
  real(mytype), save, allocatable, dimension(:,:,:) :: dudxsum,utmapsum

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
  real(mytype), save, allocatable, dimension(:,:,:) :: sgsx1,sgsy1,sgsz1,nut1,sxx1,syy1,szz1,sxy1,sxz1,syz1
  real(mytype), save, allocatable, dimension(:,:,:) :: sgsx2,sgsy2,sgsz2,nut2,sxx2,syy2,szz2,sxy2,sxz2,syz2
  real(mytype), save, allocatable, dimension(:,:,:) :: sgsx3,sgsy3,sgsz3,nut3,sxx3,syy3,szz3,sxy3,sxz3,syz3
  real(mytype), save, allocatable, dimension(:,:,:) :: gxx1,gyx1,gzx1,gxy1,gyy1,gzy1,gxz1,gyz1,gzz1
  real(mytype), save, allocatable, dimension(:,:,:) :: gxy2,gyy2,gzy2,gxz2,gyz2,gzz2
  real(mytype), save, allocatable, dimension(:,:,:) :: gxz3,gyz3,gzz3
  real(mytype), save, allocatable, dimension(:,:,:,:) :: sgsphi1,sgsphi2,sgsphi3

  real(mytype), save, allocatable, dimension(:,:,:) :: srt_smag, srt_smag2, srt_wale, srt_wale2

contains

  subroutine init_variables

    TYPE(DECOMP_INFO), save :: ph! decomposition object

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init_variables start'
#endif

    if (nrank==0) print *,'Initializing variables...'

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
    call alloc_x(uy1, opt_global=.true.) !global indices
    call alloc_x(uz1, opt_global=.true.) !global indices
    call alloc_x(pz1, opt_global=.true.) !global indices
    call alloc_x(px1, opt_global=.true.) !global indices
    call alloc_x(py1, opt_global=.true.) !global indices
    call alloc_x(diss1, opt_global=.true.) !global indices
    call alloc_x(pre1, opt_global=.true.) !global indices

    allocate(phi1(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),1:numscalar)) !global indices

    call alloc_x(ta1);call alloc_x(tb1);call alloc_x(tc1)
    call alloc_x(td1);call alloc_x(te1);call alloc_x(tf1)
    call alloc_x(tg1);call alloc_x(th1);call alloc_x(ti1)
    call alloc_x(di1);call alloc_x(ep1)
    call alloc_x(mu1)
    mu1(:,:,:) = one

    allocate(pp1(nxmsize,xsize(2),xsize(3)))
    allocate(pgy1(nxmsize,xsize(2),xsize(3)))
    allocate(pgz1(nxmsize,xsize(2),xsize(3)))

    !inflow/ouflow 2d arrays
    allocate(bxx1(xsize(2),xsize(3)),bxy1(xsize(2),xsize(3)))
    allocate(bxz1(xsize(2),xsize(3)),bxxn(xsize(2),xsize(3)))
    allocate(bxyn(xsize(2),xsize(3)),bxzn(xsize(2),xsize(3)))
    allocate(byx1(xsize(1),xsize(3)),byy1(xsize(1),xsize(3)))
    allocate(byz1(xsize(1),xsize(3)),byxn(xsize(1),xsize(3)))
    allocate(byyn(xsize(1),xsize(3)),byzn(xsize(1),xsize(3)))
    allocate(bzx1(xsize(1),xsize(2)),bzy1(xsize(1),xsize(2)))
    allocate(bzz1(xsize(1),xsize(2)),bzxn(xsize(1),xsize(2)))
    allocate(bzyn(xsize(1),xsize(2)),bzzn(xsize(1),xsize(2)))
    allocate(bxo(xsize(2),xsize(3)))
    allocate(byo(xsize(2),xsize(3)))
    allocate(bzo(xsize(2),xsize(3)))

    bxx1=zero;bxy1=zero;bxz1=zero
    byx1=zero;byy1=zero;byz1=zero
    bzx1=zero;bzy1=zero;bzz1=zero

    bxxn=zero;bxyn=zero;bxzn=zero
    byxn=zero;byyn=zero;byzn=zero
    bzxn=zero;bzyn=zero;bzzn=zero

    !pre_correc 2d array
    allocate(dpdyx1(xsize(2),xsize(3)),dpdyxn(xsize(2),xsize(3)))
    allocate(dpdzx1(xsize(2),xsize(3)),dpdzxn(xsize(2),xsize(3)))
    allocate(dpdxy1(xsize(1),xsize(3)),dpdxyn(xsize(1),xsize(3)))
    allocate(dpdzy1(xsize(1),xsize(3)),dpdzyn(xsize(1),xsize(3)))
    allocate(dpdxz1(xsize(1),xsize(2)),dpdxzn(xsize(1),xsize(2)))
    allocate(dpdyz1(xsize(1),xsize(2)),dpdyzn(xsize(1),xsize(2)))

    dpdyx1=zero
    dpdzx1=zero
    dpdyxn=zero
    dpdzxn=zero
    !!
    dpdxy1=zero
    dpdzy1=zero
    dpdxyn=zero
    dpdzyn=zero
    !!
    dpdxz1=zero
    dpdyz1=zero
    dpdxzn=zero
    dpdyzn=zero

    !arrays for visualization!pay attention to the size!
    allocate(uvisu(xstV(1):xenV(1),xstV(2):xenV(2),xstV(3):xenV(3)))

    !arrays statistics
    allocate (umean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (vmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (wmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (pmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (uumean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (vvmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (wwmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (uvmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (uwmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (vwmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (phimean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3),numscalar))
    allocate (phiphimean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3),numscalar))
    allocate (tmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))


    !Y PENCILS
    call alloc_y(ux2);call alloc_y(uy2);call alloc_y(uz2)
    call alloc_y(ta2);call alloc_y(tb2);call alloc_y(tc2)
    call alloc_y(td2);call alloc_y(te2);call alloc_y(tf2)
    call alloc_y(tg2);call alloc_y(th2);call alloc_y(ti2)
    call alloc_y(tj2);call alloc_y(di2)
    allocate(phi2(ysize(1),ysize(2),ysize(3),1:numscalar))
    allocate(pgz2(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)))
    allocate(pp2(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)))
    allocate(dip2(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)))
    allocate(ppi2(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)))
    allocate(pgy2(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)))
    allocate(pgzi2(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)))

    allocate(duxdxp2(ph1%yst(1):ph1%yen(1),ysize(2),ysize(3)))
    allocate(uyp2(ph1%yst(1):ph1%yen(1),ysize(2),ysize(3)))
    allocate(uzp2(ph1%yst(1):ph1%yen(1),ysize(2),ysize(3)))
    allocate(dipp2(ph1%yst(1):ph1%yen(1),ysize(2),ysize(3)))
    allocate(upi2(ph1%yst(1):ph1%yen(1),nymsize,ysize(3)))
    allocate(duydypi2(ph1%yst(1):ph1%yen(1),nymsize,ysize(3)))

    !Z PENCILS
    call alloc_z(ux3);call alloc_z(uy3);call alloc_z(uz3)
    call alloc_z(ta3);call alloc_z(tb3);call alloc_z(tc3)
    call alloc_z(td3);call alloc_z(te3);call alloc_z(tf3)
    call alloc_z(tg3);call alloc_z(th3);call alloc_z(ti3)
    call alloc_z(di3)
    allocate(phi3(zsize(1),zsize(2),zsize(3),1:numscalar))
    allocate(pgz3(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)))
    allocate(ppi3(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)))
    allocate(dip3(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)))

    allocate(duxydxyp3(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),zsize(3)))
    allocate(uzp3(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),zsize(3)))
    allocate(dipp3(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),zsize(3)))

    ! if all periodic
    !   allocate (pp3(ph%zst(1):ph%zen(1),ph%zst(2):ph%zen(2),ph%zst(3):ph%zen(3)))
    !   allocate (dv3(ph%zst(1):ph%zen(1),ph%zst(2):ph%zen(2),ph%zst(3):ph%zen(3)))
    !   allocate (po3(ph%zst(1):ph%zen(1),ph%zst(2):ph%zen(2),ph%zst(3):ph%zen(3)))
    allocate(pp3(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress))
    call alloc_z(dv3,ph,.true.)
    call alloc_z(po3,ph,.true.)

    if(ilesmod.ne.0.and.jLES.gt.0) then
       call alloc_x(sgsx1);call alloc_x(sgsy1); call alloc_x(sgsz1)
       call alloc_x(sxx1);call alloc_x(syy1); call alloc_x(szz1)
       call alloc_x(sxy1);call alloc_x(sxz1); call alloc_x(syz1)
       call alloc_x(nut1);call alloc_x(srt_smag); call alloc_x(srt_wale)
       call alloc_y(sgsx2);call alloc_y(sgsy2); call alloc_y(sgsz2)
       call alloc_y(sxx2) ;call alloc_y(syy2);  call alloc_y(szz2)
       call alloc_y(sxy2) ;call alloc_y(sxz2);  call alloc_y(syz2)
       call alloc_y(nut2) ;call alloc_y(srt_smag2); call alloc_y(srt_wale2)
       call alloc_z(sgsx3);call alloc_z(sgsy3); call alloc_z(sgsz3)
       call alloc_z(sxx3) ;call alloc_z(syy3);  call alloc_z(szz3)
       call alloc_z(sxy3) ;call alloc_z(sxz3);  call alloc_z(syz3)
       call alloc_z(nut3)
       call alloc_x(gxx1); call alloc_x(gyx1); call alloc_x(gzx1); call alloc_x(gxy1)
       call alloc_x(gyy1); call alloc_x(gzy1); call alloc_x(gxz1); call alloc_x(gyz1); call alloc_x(gzz1)
       call alloc_y(gxy2); call alloc_y(gyy2); call alloc_y(gzy2); call alloc_y(gxz2); call alloc_y(gyz2); call alloc_y(gzz2)
       call alloc_z(gxz3); call alloc_z(gyz3); call alloc_z(gzz3)
       allocate(sgsphi1(xsize(1),xsize(2),xsize(3),1:numscalar))
       allocate(sgsphi2(ysize(1),ysize(2),ysize(3),1:numscalar))
       allocate(sgsphi3(zsize(1),zsize(2),zsize(3),1:numscalar))
    endif


    if (iibm.ne.0) then
       !complex_geometry
       nxraf=(nxm)*nraf+1;nyraf=(nym)*nraf+1;nzraf=(nzm)*nraf+1
       allocate(nobjx(xsize(2),xsize(3)))
       allocate(nobjy(ysize(1),ysize(3)))
       allocate(nobjz(zsize(1),zsize(2)))
       allocate(nxipif(0:nobjmax,xsize(2),xsize(3)))
       allocate(nxfpif(0:nobjmax,xsize(2),xsize(3)))
       allocate(nyipif(0:nobjmax,ysize(1),ysize(3)))
       allocate(nyfpif(0:nobjmax,ysize(1),ysize(3)))
       allocate(nzipif(0:nobjmax,zsize(1),zsize(2)))
       allocate(nzfpif(0:nobjmax,zsize(1),zsize(2)))
       allocate(xi(nobjmax,xsize(2),xsize(3)))
       allocate(xf(nobjmax,xsize(2),xsize(3)))
       allocate(yi(nobjmax,ysize(1),ysize(3)))
       allocate(yf(nobjmax,ysize(1),ysize(3)))
       allocate(zi(nobjmax,zsize(1),zsize(2)))
       allocate(zf(nobjmax,zsize(1),zsize(2)))
    endif

    !module filter
    allocate(fiffx(nx), fisfx(nx), fifsx(nx), fifwx(nx), fissx(nx), fiswx(nx))
    allocate(fiffxp(nx),fisfxp(nx),fifsxp(nx),fifwxp(nx),fissxp(nx),fiswxp(nx))
    allocate(fiffy(ny), fisfy(ny), fifsy(ny), fifwy(ny), fissy(ny), fiswy(ny))
    allocate(fiffyp(ny),fisfyp(ny),fifsyp(ny),fifwyp(ny),fissyp(ny),fiswyp(ny))
    allocate(fiffz(nz), fisfz(nz), fifsz(nz), fifwz(nz), fissz(nz), fiswz(nz))
    allocate(fiffzp(nz),fisfzp(nz),fifszp(nz),fifwzp(nz),fisszp(nz),fiswzp(nz))
    allocate(fisx(xsize(2),xsize(3)),fivx(xsize(2),xsize(3)))
    allocate(fisy(ysize(1),ysize(3)),fivy(ysize(1),ysize(3)))
    allocate(fisz(zsize(1),zsize(2)),fivz(zsize(1),zsize(2)))

    !module derivative

    allocate(ffx(nx),sfx(nx),fsx(nx),fwx(nx),ssx(nx),swx(nx))
    allocate(ffxp(nx),sfxp(nx),fsxp(nx),fwxp(nx),ssxp(nx),swxp(nx))
    allocate(ffy(ny),sfy(ny),fsy(ny),fwy(ny),ssy(ny),swy(ny))
    allocate(ffyp(ny),sfyp(ny),fsyp(ny),fwyp(ny),ssyp(ny),swyp(ny))
    allocate(ffz(nz),sfz(nz),fsz(nz),fwz(nz),ssz(nz),swz(nz))
    allocate(ffzp(nz),sfzp(nz),fszp(nz),fwzp(nz),sszp(nz),swzp(nz))

    allocate(ffxS(nx),sfxS(nx),fsxS(nx),fwxS(nx),ssxS(nx),swxS(nx))
    allocate(ffxpS(nx),sfxpS(nx),fsxpS(nx),fwxpS(nx),ssxpS(nx),swxpS(nx))
    allocate(ffyS(ny),sfyS(ny),fsyS(ny),fwyS(ny),ssyS(ny),swyS(ny))
    allocate(ffypS(ny),sfypS(ny),fsypS(ny),fwypS(ny),ssypS(ny),swypS(ny))
    allocate(ffzS(nz),sfzS(nz),fszS(nz),fwzS(nz),sszS(nz),swzS(nz))
    allocate(ffzpS(nz),sfzpS(nz),fszpS(nz),fwzpS(nz),sszpS(nz),swzpS(nz))

    allocate(sx(xsize(2),xsize(3)),vx(xsize(2),xsize(3)))
    allocate(sy(ysize(1),ysize(3)),vy(ysize(1),ysize(3)))
    allocate(sz(zsize(1),zsize(2)),vz(zsize(1),zsize(2)))

    !O6SVV
    allocate(newsm(ny),newtm(ny),newsmt(ny),newtmt(ny))
    !allocate(newrm(ny),ttm(ny),newrmt(ny),ttmt(ny))
    allocate(newrm(ny),newrmt(ny))

    !module derpres
    allocate(cfx6(nxm),ccx6(nxm),cbx6(nxm),cfxp6(nxm),ciwxp6(nxm),csxp6(nxm),&
         cwxp6(nxm),csx6(nxm),cwx6(nxm),cifx6(nxm),cicx6(nxm),cisx6(nxm))
    allocate(cibx6(nxm),cifxp6(nxm),cisxp6(nxm),ciwx6(nxm))
    allocate(cfi6(nx),cci6(nx),cbi6(nx),cfip6(nx),csip6(nx),cwip6(nx),csi6(nx),&
         cwi6(nx),cifi6(nx),cici6(nx),cibi6(nx),cifip6(nx))
    allocate(cisip6(nx),ciwip6(nx),cisi6(nx),ciwi6(nx))
    allocate(cfy6(nym),ccy6(nym),cby6(nym),cfyp6(nym),csyp6(nym),cwyp6(nym),csy6(nym))
    allocate(cwy6(nym),cify6(nym),cicy6(nym),ciby6(nym),cifyp6(nym),cisyp6(nym),&
         ciwyp6(nym),cisy6(nym),ciwy6(nym))
    allocate(cfi6y(ny),cci6y(ny),cbi6y(ny),cfip6y(ny),csip6y(ny),cwip6y(ny),&
         csi6y(ny),cwi6y(ny),cifi6y(ny),cici6y(ny))
    allocate(cibi6y(ny),cifip6y(ny),cisip6y(ny),ciwip6y(ny),cisi6y(ny),ciwi6y(ny))
    allocate(cfz6(nzm),ccz6(nzm),cbz6(nzm),cfzp6(nzm),cszp6(nzm),cwzp6(nzm),csz6(nzm))
    allocate(cwz6(nzm),cifz6(nzm),cicz6(nzm),cibz6(nzm),cifzp6(nzm),ciszp6(nzm),&
         ciwzp6(nzm),cisz6(nzm),ciwz6(nzm))
    allocate(cfi6z(nz),cci6z(nz),cbi6z(nz),cfip6z(nz),csip6z(nz),cwip6z(nz),&
         csi6z(nz),cwi6z(nz),cifi6z(nz),cici6z(nz))
    allocate(cibi6z(nz),cifip6z(nz),cisip6z(nz),ciwip6z(nz),cisi6z(nz),ciwi6z(nz))

    !module waves
    allocate(zkz(nz/2+1),zk2(nz/2+1),ezs(nz/2+1))
    allocate(yky(ny),yk2(ny),eys(ny))
    allocate(xkx(nx),xk2(nx),exs(nx))

    !module mesh
    allocate(ppy(ny),pp2y(ny),pp4y(ny))
    allocate(ppyi(ny),pp2yi(ny),pp4yi(ny))
    allocate(yp(ny),ypi(ny),del(ny))
    allocate(yeta(ny),yetai(ny))

    if (istret.eq.0) then
       do j=1,ny
          yp(j)=real(j-1,mytype)*dy
          ypi(j)=(real(j,mytype)-half)*dy
       enddo
    else
       call stretching()
    endif

    adt(:)=zero ; bdt(:)=zero ; cdt(:)=zero ; gdt(:)=zero
    if (itimescheme.eq.1) then ! Euler
       iadvance_time=1
       adt(1)=1.0_mytype*dt
       bdt(1)=0.0_mytype*dt
       gdt(1)=adt(1)+bdt(1)
       gdt(3)=gdt(1)

       ntime = 1
       nrhotime = 2
    elseif (itimescheme.eq.2) then ! AB2
       iadvance_time=1
       adt(1)=1.5_mytype*dt
       bdt(1)=-0.5_mytype*dt
       gdt(1)=adt(1)+bdt(1)
       gdt(3)=gdt(1)

       ntime = 2
       nrhotime = 3
    elseif (itimescheme.eq.3) then ! AB3
       iadvance_time=1
       adt(1)= (23._mytype/12._mytype)*dt
       bdt(1)=-(16._mytype/12._mytype)*dt
       cdt(1)= ( 5._mytype/12._mytype)*dt
       gdt(1)=adt(1)+bdt(1)+cdt(1)
       gdt(3)=gdt(1)

       ntime = 3
       nrhotime = 4
    elseif(itimescheme==4) then  ! AB4
       iadvance_time=1
       adt(1)=(55.0_mytype/24.0_mytype)*dt
       bdt(1)=-(59.0_mytype/24.0_mytype)*dt
       cdt(1)=(37.0_mytype/24.0_mytype)*dt
       ddt(1)=-(9.0_mytype/24.0_mytype)*dt
       gdt(1)=adt(1)+bdt(1)+cdt(1)+ddt(1)
       gdt(3)=gdt(1)

       ntime = 4
       nrhotime = 5
    elseif(itimescheme.eq.5) then !RK3
       iadvance_time=3
       adt(1)=(8._mytype/15._mytype)*dt
       bdt(1)=0._mytype
       gdt(1)=adt(1)
       adt(2)=(5._mytype/12._mytype)*dt
       bdt(2)=(-17._mytype/60._mytype)*dt
       gdt(2)=adt(2)+bdt(2)
       adt(3)=(3._mytype/4._mytype)*dt
       bdt(3)=(-5._mytype/12._mytype)*dt
       gdt(3)=adt(3)+bdt(3)

       ntime = 2
       nrhotime = 3
    elseif(itimescheme.eq.6) then !RK4 Carpenter and Kennedy
       iadvance_time=5
       adt(1)=0.0_mytype
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

    elseif(itimescheme.eq.7) then !Semi-implicit
       iadvance_time=1
       adt(1)= (23./12.)*dt
       bdt(1)=-(16./12.)*dt
       cdt(1)= ( 5./12.)*dt
       gdt(1)=adt(1)+bdt(1)+cdt(1)
       gdt(3)=gdt(1)

       ntime = 3
       nrhotime = 4
    endif
    allocate(dux1(xsize(1),xsize(2),xsize(3),ntime))
    allocate(duy1(xsize(1),xsize(2),xsize(3),ntime))
    allocate(duz1(xsize(1),xsize(2),xsize(3),ntime))

    !! Scalar
    allocate(dphi1(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),ntime,1:numscalar)) !global indices

    !! LMN
    if (.not.ilmn) then
       nrhotime = 1 !! Save some space
    endif
    allocate(rho1(xsize(1),xsize(2),xsize(3),nrhotime)) !Need to store old density values to extrapolate drhodt
    call alloc_y(rho2)
    call alloc_z(rho3)
    allocate(drho1(xsize(1),xsize(2),xsize(3),ntime))
    rho1(:,:,:,:) = one
    call alloc_z(divu3, opt_global=.true.) !global indices

    ! !TRIPPING
    !A_trip=1.
    zs_param=1.7
    randomseed=4600
    zs_tr=zs_param/2.853
    !z_modes=int(zlz /zs_tr)
    z_modes=zlz /zs_tr
    allocate(h_coeff1(z_modes),h_coeff2(z_modes))
    allocate(phase1(z_modes),phase2(z_modes))
    allocate(h_1(xsize(3)), h_2(xsize(3)))

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init_variables done'
#endif
    return
  end subroutine init_variables
end module var
