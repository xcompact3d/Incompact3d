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

!================================================================================
!
!  SUBROUTINE: parameter
! DESCRIPTION: Reads the input.i3d file and sets the parameters of the
!              simulation.
!      AUTHOR: Paul Bartholomew <paul.bartholomew08@imperial.ac.uk>
!
!================================================================================
subroutine parameter(input_i3d)

  USE iso_fortran_env

  USE param
  USE variables
  USE complex_geometry
  USE decomp_2d
  USE ibm

  USE var, ONLY : dphi1

  USE lockexch, ONLY : pfront

  USE forces, ONLY : nvol, xld, xrd, yld, yud

  implicit none

  character(len=80), intent(in) :: input_i3d
  real(mytype) :: theta, cfl,cf2
  integer :: longueur ,impi,j, is, total

  NAMELIST /BasicParam/ p_row, p_col, nx, ny, nz, istret, beta, xlx, yly, zlz, &
       itype, iin, re, u1, u2, init_noise, inflow_noise, &
       dt, ifirst, ilast, &
       numscalar, iibm, ilmn, &
       ilesmod, iscalar, &
       nclx1, nclxn, ncly1, nclyn, nclz1, nclzn, &
       ivisu, ipost, &
       gravx, gravy, gravz
  NAMELIST /NumOptions/ ifirstder, isecondder, itimescheme, nu0nu, cnu, fpi2, ipinter
  NAMELIST /InOutParam/ irestart, icheckpoint, ioutput, nvisu, iprocessing
  NAMELIST /Statistics/ wrotation,spinup_time, nstat, initstat
  NAMELIST /ScalarParam/ sc, ri, uset, cp, &
       nclxS1, nclxSn, nclyS1, nclySn, nclzS1, nclzSn, &
       scalar_lbound, scalar_ubound
  NAMELIST /LESModel/ jles, smagcst, walecst, maxdsmagcst, iwall
  NAMELIST /WallModel/ smagwalldamp

  NAMELIST /ibmstuff/ cex,cey,ra,nobjmax,nraf,nvol
  NAMELIST /ForceCVs/ xld, xrd, yld, yud
  NAMELIST /LMN/ dens1, dens2, prandtl, ilmn_bound, ivarcoeff, ilmn_solve_temp, &
       massfrac, mol_weight, imultispecies, primary_species, &
       Fr, ibirman_eos
  NAMELIST /CASE/ tgv_twod, pfront
#ifdef DEBG
  if (nrank .eq. 0) print *,'# parameter start'
#endif

  if (nrank==0) then
     print *,'==========================================================='
     print *,'======================Xcompact3D==========================='
     print *,'===Copyright (c) 2018 Eric Lamballais and Sylvain Laizet==='
     print *,'===Modified by Felipe Schuch and Ricardo Frantz============'
     print *,'===Modified by Paul Bartholomew, Georgios Deskos and======='
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
  if (iibm.ne.0) then
     read(10, nml=ibmstuff)
     if (nvol.gt.0) then
        iforces = .TRUE.
     endif
  endif

  if (iforces) then
     allocate(xld(nvol), xrd(nvol), yld(nvol), yud(nvol))
     read(10, nml=ForceCVs)
  endif

  if (numscalar.ne.0) then
     iscalar = 1
     !! Set Scalar BCs same as fluid (may be overridden)
     nclxS1 = nclx1; nclxSn = nclxn
     nclyS1 = ncly1; nclySn = nclyn
     nclzS1 = nclz1; nclzSn = nclzn

     !! Allocate scalar arrays and set sensible defaults
     allocate(massfrac(numscalar))
     allocate(mol_weight(numscalar))
     massfrac(:) = .FALSE.
     mol_weight(:) = one
     allocate(sc(numscalar), ri(numscalar), uset(numscalar), cp(numscalar))
     ri(:) = zero
     uset(:) = zero
     cp(:) = zero

     allocate(scalar_lbound(numscalar), scalar_ubound(numscalar))
     scalar_lbound(:) = -huge(one)
     scalar_ubound(:) = huge(one)
  endif

  if (ilmn) then
     read(10, nml=LMN)

     do is = 1, numscalar
        if (massfrac(is)) then
           imultispecies = .TRUE.
        endif
     enddo

     if (imultispecies) then
        if (primary_species.lt.1) then
           if (nrank.eq.0) then
              print *, "Error: you must set a primary species for multispecies flow"
              print *, "       solver will enforce Y_p = 1 - sum_s Y_s, s != p."
              stop
           endif
        else if (.not.massfrac(primary_species)) then
           if (nrank.eq.0) then
              print *, "Error: primary species must be a massfraction!"
           endif
        endif
     endif
  endif
  if (numscalar.ne.0) then
     read(10, nml=ScalarParam)
  endif
  ! !! These are the 'optional'/model parameters
  ! read(10, nml=ScalarParam)
  if(ilesmod==0) then
    nu0nu=four
    cnu=0.44_mytype
  endif
  if(ilesmod.ne.0) read(10, nml=LESModel)
  ! read(10, nml=TurbulenceWallModel)
  read(10, nml=CASE) !! Read case-specific variables
  close(10)

  ! allocate(sc(numscalar),cp(numscalar),ri(numscalar),group(numscalar))

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
  if (nrank .eq. 0) print *,'# parameter input.i3d done'
#endif
  if (nrank==0) then
     print *,''
     if (itype.eq.itype_user) then
        print *,'User-defined simulation'
     elseif (itype.eq.itype_lockexch) then
        print *,'Simulating lock-exchange'
     elseif (itype.eq.itype_tgv) then
        print *,'Simulating TGV'
     elseif (itype.eq.itype_channel) then
        print *,'Simulating channel'
     elseif (itype.eq.itype_hill) then
        print *,'Simulating periodic hill'
     elseif (itype.eq.itype_cyl) then
        print *,'Simulating cylinder'
     elseif (itype.eq.itype_dbg) then
        print *,'Debug schemes'
     elseif (itype.eq.itype_mixlayer) then
        print *,'Mixing layer'
     elseif (itype.eq.itype_jet) then
        print *,'Jet'
     elseif (itype.eq.itype_tbl) then
        print *,'Turbulent boundary layer'
     else
        print *,'Unknown itype: ', itype
        stop
     endif
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
     write(*,"(' Gravity vector     : (gx, gy, gz)=(',F15.8,',',F15.8,',',F15.8,')')") gravx, gravy, gravz
     if (ilmn) then
        print *, "LMN                : Enabled"
        if (ivarcoeff) then
           print *, "LMN-Poisson solver : Variable-coefficient"
        else
           print *, "LMN-Poisson solver : Constant-coefficient"
        endif
        if (ilmn_bound) then
           print *, "LMN boundedness    : Enforced"
        else
           print *, "LMN boundedness    : Not enforced"
        endif
        write(*,"(' dens1 and dens2    : ',F6.2' ',F6.2)") dens1, dens2
        write(*,"(' Prandtl number Re  : ',F15.8)") prandtl
     endif
     write(*,"(' Time step dt       : ',F15.8)") dt
     write (*,"(' Spatial scheme     : ',F15.8)") fpi2
     if (ilesmod.ne.0) then
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
     elseif (itimescheme.eq.7) then
        print *,'Temporal scheme    : Semi-implicit'
     else
        print *,'Error: itimescheme must be specified as 1-7'
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
           write (*,"(' Scalar             : #',I1)") is
           write (*,"(' Schmidt number     : ',F15.8)") sc(is)
           write (*,"(' Richardson number  : ',F15.8)") ri(is)
        end do
     endif
     if (iibm.eq.0) then
        print *,'Immersed boundary  : off'
     elseif (iibm.eq.1) then
        print *,'Immersed boundary  : old school'
     elseif (iibm.eq.2) then
        print *,'Immersed boundary  : on with Lagrangian Poly'
     endif
     if (angle.ne.0.) write(*,"(' Solid rotation     : ',F6.2)") angle
     print *,''

     !! Print case-specific information
     if (itype==itype_lockexch) then
        print *, "Initial front location: ", pfront
     elseif (itype==itype_tgv) then
        print *, "TGV 2D: ", tgv_twod
     endif
  endif

  xnu=one/re

  if (ilmn) then
     if (ivarcoeff) then
        npress = 2 !! Need current pressure and previous iterate
     else
        npress = 1
     endif
  endif

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
  USE complex_geometry

  USE forces, ONLY : nvol

  IMPLICIT NONE

  integer :: i

  ro = 99999999._mytype
  angle = zero
  u1 = 2
  u2 = 1
  init_noise = zero
  inflow_noise = zero
  iin = 0
  itimescheme = 4
  istret = 0
  ipinter=3
  beta = 0
  iscalar = 0
  cont_phi = 0
  filepath = './data/'
  irestart = 0
  datapath = './data/'
  fpi2 = (48._mytype / seven) / (PI**2)

  !! IBM stuff
  nraf = 0
  nobjmax = 0

  nvol = 0
  iforces = .FALSE.

  itrip = 1
  wrotation = zero
  irotation = 0
  itest=1

  !! Gravity field
  gravx = zero
  gravy = zero
  gravz = zero

  !! LMN stuff
  ilmn = .FALSE.
  ilmn_bound = .TRUE.
  pressure0 = one
  prandtl = one
  dens1 = one
  dens2 = one
  ivarcoeff = .FALSE.
  npress = 1 !! By default people only need one pressure field
  ilmn_solve_temp = .FALSE.
  imultispecies = .FALSE.
  Fr = zero
  ibirman_eos = .FALSE.

  primary_species = -1

  !! IO
  ivisu = 1
  ipost = 0
  iprocessing = huge(i)
  initstat = huge(i)

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
  iibm=0
  npif = 2
  izap = 1

  imodulo2 = 1

  !! CASE specific variables
  tgv_twod = .FALSE.

end subroutine parameter_defaults
