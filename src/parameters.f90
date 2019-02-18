!################################################################################
!This file is part of Incompact3d.
!
!Incompact3d
!Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!    Incompact3d is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation.
!
!    Incompact3d is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the code.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    We kindly request that you cite Incompact3d in your publications and 
!    presentations. The following citations are suggested:
!
!    1-Laizet S. & Lamballais E., 2009, High-order compact schemes for 
!    incompressible flows: a simple and efficient method with the quasi-spectral 
!    accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
!
!    2-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence 
!    problems with up to 0(10^5) computational cores, Int. J. of Numerical 
!    Methods in Fluids, vol 67 (11), pp 1735-1757
!################################################################################

!================================================================================
!
!  SUBROUTINE: parameter
! DESCRIPTION: Reads the input.i3d file and sets the parameters of the
!              simulation.
!      AUTHOR: Paul Bartholomew <paul.bartholomew08@imperial.ac.uk>
!
!================================================================================
subroutine parameter(input_i3d)

  USE param
  USE variables
  USE complex_geometry
  USE decomp_2d

  implicit none

  character(len=80), intent(in) :: input_i3d
  real(mytype) :: theta, cfl,cf2
  integer :: longueur ,impi,j, is, total

  NAMELIST /BasicParam/ p_row, p_col, nx, ny, nz, istret, beta, xlx, yly, zlz, &
       itype, iin, re, u1, u2, init_noise, inflow_noise, &
       dt, ifirst, ilast, &
       iturbmod, iscalar, iibm, &
       nclx1, nclxn, ncly1, nclyn, nclz1, nclzn, &
       ivisu, ipost
  NAMELIST /NumOptions/ ifirstder, isecondder, itimescheme, rxxnu, cnu, fpi2, ipinter
  NAMELIST /InOutParam/ irestart, icheckpoint, ioutput, nvisu
  NAMELIST /Statistics/ wrotation,spinup_time, nstat, initstat
  NAMELIST /ScalarParam/ numscalar, sc
  NAMELIST /TurbulenceModel/ iles, smagcst, walecst, iwall
  NAMELIST /TurbulenceWallModel/ smagwalldamp

#ifdef DEBG
  if (nrank .eq. 0) print *,'# parameter start'
#endif

  if (nrank==0) then
     print *,'==========================================================='
     print *,'======================Xcompact3D==========================='
     print *,'===Copyright (c) 2018 Eric Lamballais and Sylvain Laizet==='
     print *,'===Modified by Felipe Schuch and Ricardo Frantz============'
     print *,'===Modified by Paul Bartholomew, Yorgos Deskos and========='
     print *,'===Sylvain Laizet -- 2018- ================================'
     print *,'==========================================================='
#if defined(VERSION)
     write(*,*)'Git version        : ', VERSION
#else
     write(*,*)'Git version        : unknown'
#endif
  endif

  call parameter_defaults()

  !! Read parameters
  open(10, file=input_i3d)

  !! These are the 'essential' parameters
  read(10, nml=BasicParam)
  read(10, nml=NumOptions)
  read(10, nml=InOutParam)
  read(10, nml=Statistics)

  ! !! These are the 'optional'/model parameters
  ! read(10, nml=ScalarParam)
  ! read(10, nml=TurbulenceModel)
  ! read(10, nml=TurbulenceWallModel)

  close(10)

! !!! BAD
!   xlx = pi
!   yly = pi
!   zlz = pi

  jLES = iles

  allocate(sc(numscalar),cp(numscalar),ri(numscalar),group(numscalar))

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
     if (iturbmod.ne.0) then
        print *,'                   : DNS'
     else
        if (jLES.eq.1) then
           print *,'                   : Phys Smag'
        else if (jLES.eq.2) then
           print *,'                   : Phys WALE'
        else if (jLES.eq.3) then
           print *,'                   : Phys dyn. Smag'
        else if (jLES.eq.4) then
           print *,'                   : iSVV'
        else
        endif
     endif
     if (itimescheme.eq.1) then
        print *,'Temporal scheme    : Forwards Euler'
     elseif (itimescheme.eq.2) then
        print *,'Temporal scheme    : Adams-bashforth 2'
     elseif (itimescheme.eq.3) then
        print *,'Temporal scheme    : Adams-bashforth 3'
     elseif (itimescheme.eq.4) then
        print *,'Temporal scheme    : Adams-bashforth 4'
        print *,'Error: Adams-bashforth 4 not implemented!'
        stop
     elseif (itimescheme.eq.5) then
        print *,'Temporal scheme    : Runge-kutta 3'
     elseif (itimescheme.eq.6) then
        print *,'Temporal scheme    : Runge-kutta 4'
        print *,'Error: Runge-kutta 4 not implemented!'
        stop
     else
        print *,'Error: itimescheme must be specified as 1-6'
        stop
     endif
     if (iscalar.eq.0) then
        print *,'Scalar             : off'
     else
        print *,'Scalar             : on'
        write(*,"(' Boundary condition : (nclxS1,nclxSn)=(',I1,',',I1,')')") nclxS1,nclxSn
        write(*,"('                      (nclyS1,nclySn)=(',I1,',',I1,')')") nclyS1,nclySn
        write(*,"('                      (nclzS1,nclzSn)=(',I1,',',I1,')')") nclzS1,nclzSn
        do is=1, numscalar
           write (*,"(' Particle fraction  : #',I1)") is
           write (*,"(' Concentration      : ',F15.8)") cp(is)
           write (*,"(' Richardson number  : ',F15.8)") ri(is)
           write (*,"(' Settling velocity  : ',F15.8)") uset(is)
           write (*,"(' Schmidt number     : ',F15.8)") sc(is)
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

  dx2 = dx * dx
  dy2 = dy * dy
  dz2 = dz * dz

#ifdef DEBG
  if (nrank .eq. 0) print *,'# parameter done'
#endif

  return
end subroutine parameter

!================================================================================
!
!  SUBROUTINE: parameter_defaults
! DESCRIPTION: Sets the default simulation parameters.
!      AUTHOR: Paul Bartholomew <paul.bartholomew08@imperial.ac.uk>
!
!================================================================================
subroutine parameter_defaults()

  USE param
  USE variables
  USE decomp_2d

  IMPLICIT NONE

  ro = 99999999._mytype
  angle = zero
  u1 = 2
  u2 = 1
  init_noise = zero
  inflow_noise = zero
  iin = 0
  itimescheme = 4
  istret = 0
  ipinter=2
  beta = 0
  iscalar = 0
  cont_phi = 0
  filepath = './data/'
  irestart = 0
  datapath = './data/'
  fpi2 = 4.
  ! nraf = 10
  ! nobjmax = 1
  itrip = 0
  wrotation = zero
  irotation = 0
  itest=1

  !! IO
  ivisu = 1
  ipost = 0
  
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

  ipost = 0

  ivirt = 0
  ilag = 0
  npif = 2
  izap = 1

  imodulo2 = 1

end subroutine parameter_defaults
