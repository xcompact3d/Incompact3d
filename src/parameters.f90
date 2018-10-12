subroutine parameter()

  USE param
  USE variables
  USE complex_geometry
  USE decomp_2d

  implicit none

  real(mytype) :: theta, cfl,cf2
  integer :: longueur ,impi,j, is, total


#ifdef DEBG
  if (nrank .eq. 0) print *,'# parameter start'
#endif

  if (nrank==0) then
     print *,'==========================================================='
     print *,'======================Incompact3d=========================='
     print *,'===Copyright (c) 2018 Eric Lamballais and Sylvain Laizet==='
     print *,'===Modified by Felipe Schuch and Ricardo Frantz============'
     print *,'==========================================================='
#if defined(VERSION)
     write(*,*)'Git version        : ', VERSION
#else
     write(*,*)'Git version        : unknown'
#endif
  endif

  allocate(nsc(nphi),uset(nphi),cp(nphi),ri(nphi),group(nphi))

  ro = 99999999._mytype
  ri = zero
  nsc = one
  uset = zero
  cp = 1
  angle = zero
  u1 = 2
  u2 = 1
  noise = zero
  noise1 = zero
  iin = 0
  nscheme = 4
  istret = 0
  beta = 0
  iscalar = 0
  cont_phi = 0
  filepath = './data/'
  ilit = 0
  datapath = './data/'
  fpi2 = 4.
  nraf = 10
  nobjmax = 1
  itrip = 0
  wrotation = zero
  irotation = 0
  itest=1
  
  
  save_ux = 0
  save_uy = 0
  save_uz = 0
  save_phi = 0
  save_uxm = 0
  save_uym = 0
  save_uzm = 0
  save_phim = 0
  save_w = 0
  save_w1 = 0
  save_w2 = 0
  save_w3 = 0
  save_qc = 0
  save_pc = 0
  save_V = 0
  save_dudx = 0
  save_dudy = 0
  save_dudz = 0
  save_dvdx = 0
  save_dvdy = 0
  save_dvdz = 0
  save_dwdx = 0
  save_dwdy = 0
  save_dwdz = 0
  save_dphidx = 0
  save_dphidy = 0
  save_dphidz = 0
  save_pre = 0
  save_prem = 0
  save_dmap = 0
  save_utmap = 0
  save_ibm = 0

  ivirt = 0
  ilag = 0
  npif = 2
  izap = 1

  call ft_parameter(.false.)

  if (nclx1.eq.0.and.nclxn.eq.0) then
     nclx=.true.
     nxm=nx
  else
     nclx=.false.
     nxm=nx-1
  endif
  if (ncly1.eq.0.and.nclyn.eq.0) then
     ncly=.true.
     nym=ny
  else
     ncly=.false.
     nym=ny-1
  endif
  if (nclz1.eq.0.and.nclzn.eq.0) then
     nclz=.true.
     nzm=nz
  else
     nclz=.false.
     nzm=nz-1
  endif

  dx=xlx/real(nxm,mytype)
  dy=yly/real(nym,mytype)
  dz=zlz/real(nzm,mytype)

  if (nrank==0) call system('mkdir data out probes 2> /dev/null')

#ifdef DEBG
  if (nrank .eq. 0) print *,'# parameter incompact3d.prm done'
#endif
  if (nrank==0) then
     print *,''
     print *,'(lx,ly,lz)=',xlx,yly,zlz
     print *,'(nx,ny,nz)=',nx,ny,nz
     print *,'(dx,dy,dz)=',dx,dy,dz
     print *,'(nx*ny*nz)=',nx*ny*nz
     print *,'(p_row,p_col)=',p_row,p_col
     print *,''
#ifdef DOUBLE_PREC
#ifdef SAVE_SINGLE
     print *,'Numerical precision: Double, saving in single'
#else
     print *,'Numerical precision: Double'
#endif
#else
     print *,'Numerical precision: Single'
#endif
     write(*,"(' Boundary condition : (nclx1 ,nclxn )=(',I1,',',I1,')')") nclx1,nclxn
     write(*,"('                      (ncly1 ,nclyn )=(',I1,',',I1,')')") ncly1,nclyn
     write(*,"('                      (nclz1 ,nclzn )=(',I1,',',I1,')')") nclz1,nclzn
     write(*,"(' High and low speed : u1=',F6.2,' and u2=',F6.2)") u1,u2
     write(*,"(' Reynolds number Re : ',F15.8)") re
     write(*,"(' Time step dt       : ',F15.8)") dt
     write (*,"(' Spatial scheme     : ',F15.8)") fpi2
     if (jLES.eq.0) print *,'                   : DNS'
     if (jLES.eq.1) print *,'                   : iLES'
     if (jLES.eq.2) print *,'                   : Explicit Simple Smagorinsky'
     if (jLES.eq.3) print *,'                   : Explicit Wall-Adaptive LES'
     if (jLES.eq.4) print *,'                   : Explicit Dynamic Smagorinsky LES'
     if (nscheme.eq.1) print *,'Temporal scheme    : Adams-bashforth 2'
     if (nscheme.eq.2) print *,'Temporal scheme    : Adams-bashforth 3'
     if (nscheme.eq.3) print *,'Temporal scheme    : Runge-Kutta 3'
     if (iscalar.eq.0) print *,'Scalar             : off'
     if (iscalar.eq.1) then
        print *,'Scalar             : on'
        write(*,"(' Boundary condition : (nclxS1,nclxSn)=(',I1,',',I1,')')") nclxS1,nclxSn
        write(*,"('                      (nclyS1,nclySn)=(',I1,',',I1,')')") nclyS1,nclySn
        write(*,"('                      (nclzS1,nclzSn)=(',I1,',',I1,')')") nclzS1,nclzSn
        do is=1, nphi
           write (*,"(' Particle fraction  : #',I1)") is
           write (*,"(' Concentration      : ',F15.8)") cp(is)
           write (*,"(' Richardson number  : ',F15.8)") ri(is)
           write (*,"(' Settling velocity  : ',F15.8)") uset(is)
           write (*,"(' Schmidt number     : ',F15.8)") nsc(is)
        end do
     endif
     if (ivirt.eq.0) print *,'Immersed boundary  : off'
     if (ivirt.eq.1) print *,'Immersed boundary  : old school'
     if (ivirt.eq.2) print *,'Immersed boundary  : on with Lagrangian Poly'
     if (angle.ne.0.) write(*,"(' Solid rotation     : ',F6.2)") angle
     print *,''
  endif

  xnu=one/re

#ifdef DOUBLE_PREC
  anglex = dsin(pi*angle/180._mytype)
  angley = dcos(pi*angle/180._mytype)
#else
  anglex = sin(pi*angle/180._mytype)
  angley = cos(pi*angle/180._mytype)
#endif

  dx2=dx*dx
  dy2=dy*dy
  dz2=dz*dz

#ifdef DEBG
  if (nrank .eq. 0) print *,'# parameter done'
#endif

  return
end subroutine parameter
