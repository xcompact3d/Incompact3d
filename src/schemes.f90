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

!********************************************************************
!
subroutine schemes()
  !
  !********************************************************************

  USE param
  USE derivX
  USE derivY
  USE derivZ
  USE variables
  USE var
  USE ydiff_implicit, only : init_implicit, implicit_schemes

  implicit none

  integer :: is

#ifdef DEBG
  if (nrank  ==  0) write(*,*)'# schemes start'
#endif

  !Velocity
  ! First derivative
  if (nclx1.eq.0.and.nclxn.eq.0) derx => derx_00
  if (nclx1.eq.1.and.nclxn.eq.1) derx => derx_11
  if (nclx1.eq.1.and.nclxn.eq.2) derx => derx_12
  if (nclx1.eq.2.and.nclxn.eq.1) derx => derx_21
  if (nclx1.eq.2.and.nclxn.eq.2) derx => derx_22
  !
  if (ncly1.eq.0.and.nclyn.eq.0) dery => dery_00
  if (ncly1.eq.1.and.nclyn.eq.1) dery => dery_11
  if (ncly1.eq.1.and.nclyn.eq.2) dery => dery_12
  if (ncly1.eq.2.and.nclyn.eq.1) dery => dery_21
  if (ncly1.eq.2.and.nclyn.eq.2) dery => dery_22
  !
  if (nclz1.eq.0.and.nclzn.eq.0) derz => derz_00
  if (nclz1.eq.1.and.nclzn.eq.1) derz => derz_11
  if (nclz1.eq.1.and.nclzn.eq.2) derz => derz_12
  if (nclz1.eq.2.and.nclzn.eq.1) derz => derz_21
  if (nclz1.eq.2.and.nclzn.eq.2) derz => derz_22
  ! Second derivative
  !x
  if (nclx1.eq.0.and.nclxn.eq.0) derxx => derxx_00
  if (nclx1.eq.1.and.nclxn.eq.1) derxx => derxx_11
  if (nclx1.eq.1.and.nclxn.eq.2) derxx => derxx_12
  if (nclx1.eq.2.and.nclxn.eq.1) derxx => derxx_21
  if (nclx1.eq.2.and.nclxn.eq.2) derxx => derxx_22
  !y
  if (ncly1.eq.0.and.nclyn.eq.0) deryy => deryy_00
  if (ncly1.eq.1.and.nclyn.eq.1) deryy => deryy_11
  if (ncly1.eq.1.and.nclyn.eq.2) deryy => deryy_12
  if (ncly1.eq.2.and.nclyn.eq.1) deryy => deryy_21
  if (ncly1.eq.2.and.nclyn.eq.2) deryy => deryy_22
  !z
  if (nclz1.eq.0.and.nclzn.eq.0) derzz => derzz_00
  if (nclz1.eq.1.and.nclzn.eq.1) derzz => derzz_11
  if (nclz1.eq.1.and.nclzn.eq.2) derzz => derzz_12
  if (nclz1.eq.2.and.nclzn.eq.1) derzz => derzz_21
  if (nclz1.eq.2.and.nclzn.eq.2) derzz => derzz_22

  call first_derivative(alfa1x,af1x,bf1x,cf1x,df1x,alfa2x,af2x,alfanx,afnx,bfnx,&
       cfnx,dfnx,alfamx,afmx,alfaix,afix,bfix,&
       ffx,fsx,fwx,ffxp,fsxp,fwxp,dx,nx,nclx1,nclxn)
  call first_derivative(alfa1y,af1y,bf1y,cf1y,df1y,alfa2y,af2y,alfany,afny,bfny,&
       cfny,dfny,alfamy,afmy,alfajy,afjy,bfjy,&
       ffy,fsy,fwy,ffyp,fsyp,fwyp,dy,ny,ncly1,nclyn)
  call first_derivative(alfa1z,af1z,bf1z,cf1z,df1z,alfa2z,af2z,alfanz,afnz,bfnz,&
       cfnz,dfnz,alfamz,afmz,alfakz,afkz,bfkz,&
       ffz,fsz,fwz,ffzp,fszp,fwzp,dz,nz,nclz1,nclzn)
  call second_derivative(alsa1x,as1x,bs1x,&
       cs1x,ds1x,alsa2x,as2x,alsanx,asnx,bsnx,csnx,dsnx,alsamx,&
       asmx,alsa3x,as3x,bs3x,alsatx,astx,bstx,&
       alsa4x,as4x,bs4x,cs4x,&
       alsattx,asttx,bsttx,csttx,&
       alsaix,asix,bsix,csix,dsix,&
       sfx,ssx,swx,sfxp,ssxp,swxp,dx2,nx,nclx1,nclxn)
  call second_derivative(alsa1y,as1y,bs1y,&
       cs1y,ds1y,alsa2y,as2y,alsany,asny,bsny,csny,dsny,alsamy,&
       asmy,alsa3y,as3y,bs3y,alsaty,asty,bsty,&
       alsa4y,as4y,bs4y,cs4y,&
       alsatty,astty,bstty,cstty,&
       alsajy,asjy,bsjy,csjy,dsjy,&
       sfy,ssy,swy,sfyp,ssyp,swyp,dy2,ny,ncly1,nclyn)
  call second_derivative(alsa1z,as1z,bs1z,&
       cs1z,ds1z,alsa2z,as2z,alsanz,asnz,bsnz,csnz,dsnz,alsamz,&
       asmz,alsa3z,as3z,bs3z,alsatz,astz,bstz,&
       alsa4z,as4z,bs4z,cs4z,&
       alsattz,asttz,bsttz,csttz,&
       alsakz,askz,bskz,cskz,dskz,&
       sfz,ssz,swz,sfzp,sszp,swzp,dz2,nz,nclz1,nclzn)

  if (iscalar.ne.0 .or. (ilmn)) then
     !Scalar
     ! First derivative
     if (nclxS1.eq.0.and.nclxSn.eq.0) derxS => derx_00
     if (nclxS1.eq.1.and.nclxSn.eq.1) derxS => derx_11
     if (nclxS1.eq.1.and.nclxSn.eq.2) derxS => derx_12
     if (nclxS1.eq.2.and.nclxSn.eq.1) derxS => derx_21
     if (nclxS1.eq.2.and.nclxSn.eq.2) derxS => derx_22
     !
     if (nclyS1.eq.0.and.nclySn.eq.0) deryS => dery_00
     if (nclyS1.eq.1.and.nclySn.eq.1) deryS => dery_11
     if (nclyS1.eq.1.and.nclySn.eq.2) deryS => dery_12
     if (nclyS1.eq.2.and.nclySn.eq.1) deryS => dery_21
     if (nclyS1.eq.2.and.nclySn.eq.2) deryS => dery_22
     !
     if (nclzS1.eq.0.and.nclzSn.eq.0) derzS => derz_00
     if (nclzS1.eq.1.and.nclzSn.eq.1) derzS => derz_11
     if (nclzS1.eq.1.and.nclzSn.eq.2) derzS => derz_12
     if (nclzS1.eq.2.and.nclzSn.eq.1) derzS => derz_21
     if (nclzS1.eq.2.and.nclzSn.eq.2) derzS => derz_22
     ! Second derivative
     if (nclxS1.eq.0.and.nclxSn.eq.0) derxxS => derxx_00
     if (nclxS1.eq.1.and.nclxSn.eq.1) derxxS => derxx_11
     if (nclxS1.eq.1.and.nclxSn.eq.2) derxxS => derxx_12
     if (nclxS1.eq.2.and.nclxSn.eq.1) derxxS => derxx_21
     if (nclxS1.eq.2.and.nclxSn.eq.2) derxxS => derxx_22
     !y
     if (nclyS1.eq.0.and.nclySn.eq.0) deryyS => deryy_00
     if (nclyS1.eq.1.and.nclySn.eq.1) deryyS => deryy_11
     if (nclyS1.eq.1.and.nclySn.eq.2) deryyS => deryy_12
     if (nclyS1.eq.2.and.nclySn.eq.1) deryyS => deryy_21
     if (nclyS1.eq.2.and.nclySn.eq.2) deryyS => deryy_22
     !z
     if (nclzS1.eq.0.and.nclzSn.eq.0) derzzS => derzz_00
     if (nclzS1.eq.1.and.nclzSn.eq.1) derzzS => derzz_11
     if (nclzS1.eq.1.and.nclzSn.eq.2) derzzS => derzz_12
     if (nclzS1.eq.2.and.nclzSn.eq.1) derzzS => derzz_21
     if (nclzS1.eq.2.and.nclzSn.eq.2) derzzS => derzz_22
     call first_derivative(alfa1x,af1x,bf1x,cf1x,df1x,alfa2x,af2x,alfanx,afnx,bfnx,&
          cfnx,dfnx,alfamx,afmx,alfaix,afix,bfix,&
          ffxS,fsxS,fwxS,ffxpS,fsxpS,fwxpS,dx,nx,nclxS1,nclxSn)
     call first_derivative(alfa1y,af1y,bf1y,cf1y,df1y,alfa2y,af2y,alfany,afny,bfny,&
          cfny,dfny,alfamy,afmy,alfajy,afjy,bfjy,&
          ffyS,fsyS,fwyS,ffypS,fsypS,fwypS,dy,ny,nclyS1,nclySn)
     call first_derivative(alfa1z,af1z,bf1z,cf1z,df1z,alfa2z,af2z,alfanz,afnz,bfnz,&
          cfnz,dfnz,alfamz,afmz,alfakz,afkz,bfkz,&
          ffzS,fszS,fwzS,ffzpS,fszpS,fwzpS,dz,nz,nclzS1,nclzSn)
     call second_derivative(alsa1x,as1x,bs1x,&
          cs1x,ds1x,alsa2x,as2x,alsanx,asnx,bsnx,csnx,dsnx,alsamx,&
          asmx,alsa3x,as3x,bs3x,alsatx,astx,bstx,&
          alsa4x,as4x,bs4x,cs4x,&
          alsattx,asttx,bsttx,csttx,&
          alsaix,asix,bsix,csix,dsix,&
          sfxS,ssxS,swxS,sfxpS,ssxpS,swxpS,dx2,nx,nclxS1,nclxSn)
     call second_derivative(alsa1y,as1y,bs1y,&
          cs1y,ds1y,alsa2y,as2y,alsany,asny,bsny,csny,dsny,alsamy,&
          asmy,alsa3y,as3y,bs3y,alsaty,asty,bsty,&
          alsa4y,as4y,bs4y,cs4y,&
          alsatty,astty,bstty,cstty,&
          alsajy,asjy,bsjy,csjy,dsjy,&
          sfyS,ssyS,swyS,sfypS,ssypS,swypS,dy2,ny,nclyS1,nclySn)
     call second_derivative(alsa1z,as1z,bs1z,&
          cs1z,ds1z,alsa2z,as2z,alsanz,asnz,bsnz,csnz,dsnz,alsamz,&
          asmz,alsa3z,as3z,bs3z,alsatz,astz,bstz,&
          alsa4z,as4z,bs4z,cs4z,&
          alsattz,asttz,bsttz,csttz,&
          alsakz,askz,bskz,cskz,dskz,&
          sfzS,sszS,swzS,sfzpS,sszpS,swzpS,dz2,nz,nclzS1,nclzSn)
  endif
  call interpolation(dx,nxm,nx,nclx1,nclxn,&
       alcaix6,acix6,bcix6,&
       ailcaix6,aicix6,bicix6,cicix6,dicix6,&
       cfx6,ccx6,cbx6,cfxp6,ciwxp6,csxp6,&
       cwxp6,csx6,cwx6,cifx6,cicx6,cisx6,&
       cibx6,cifxp6,cisxp6,ciwx6,&
       cfi6,cci6,cbi6,cfip6,csip6,cwip6,csi6,&
       cwi6,cifi6,cici6,cibi6,cifip6,&
       cisip6,ciwip6,cisi6,ciwi6)
  call interpolation(dy,nym,ny,ncly1,nclyn,&
       alcaiy6,aciy6,bciy6,&
       ailcaiy6,aiciy6,biciy6,ciciy6,diciy6,&
       cfy6,ccy6,cby6,cfyp6,ciwyp6,csyp6,&
       cwyp6,csy6,cwy6,cify6,cicy6,cisy6,&
       ciby6,cifyp6,cisyp6,ciwy6,&
       cfi6y,cci6y,cbi6y,cfip6y,csip6y,cwip6y,csi6y,&
       cwi6y,cifi6y,cici6y,cibi6y,cifip6y,&
       cisip6y,ciwip6y,cisi6y,ciwi6y)
  call interpolation(dz,nzm,nz,nclz1,nclzn,&
       alcaiz6,aciz6,bciz6,&
       ailcaiz6,aiciz6,biciz6,ciciz6,diciz6,&
       cfz6,ccz6,cbz6,cfzp6,ciwzp6,cszp6,&
       cwzp6,csz6,cwz6,cifz6,cicz6,cisz6,&
       cibz6,cifzp6,ciszp6,ciwz6,&
       cfi6z,cci6z,cbi6z,cfip6z,csip6z,cwip6z,csi6z,&
       cwi6z,cifi6z,cici6z,cibi6z,cifip6z,&
       cisip6z,ciwip6z,cisi6z,ciwi6z)

  if (iimplicit.ne.0) then
     call init_implicit()
     call implicit_schemes()
  endif

#ifdef DEBG
  if (nrank  ==  0) write(*,*)'# schemes end'
#endif

  return
end subroutine schemes

!*******************************************************************
!
subroutine prepare (b,c,f,s,w,n)
  !
  !*******************************************************************

  use decomp_2d, only : mytype
  use param, only : one

  implicit none

  integer, intent(in) :: n
  real(mytype), dimension(n), intent(in) :: b,c,f
  real(mytype), dimension(n), intent(out) :: s,w
  integer :: i

  do i=1,n
     w(i)=c(i)
  enddo
  do i=2,n
     s(i)=b(i-1)/w(i-1)
     w(i)=w(i)-f(i-1)*s(i)
  enddo
  do i=1,n
     w(i)=one/w(i)
  enddo

  return
end subroutine prepare

!*******************************************************************
!
subroutine first_derivative(alfa1,af1,bf1,cf1,df1,alfa2,af2,alfan,afn,bfn,&
     cfn,dfn,alfam,afm,alfai,afi,bfi,&
     ff,fs,fw,ffp,fsp,fwp,d,n,ncl1,ncln)
  !
  !*******************************************************************

  use decomp_2d, only : mytype, nrank
  use param
  use MPI

  implicit none

  real(mytype),intent(in) :: d
  integer,intent(in) :: n,ncl1,ncln
  real(mytype),dimension(n),intent(out) :: ff,fs,fw,ffp,fsp,fwp
  real(mytype),intent(out) :: alfa1,af1,bf1,cf1,df1,alfa2,af2,alfan,afn,bfn,&
       cfn,dfn,alfam,afm,alfai,afi,bfi
  integer :: i,code,ierror
  real(mytype),dimension(n) :: fb,fc

  ff=zero;fs=zero;fw=zero;ffp=zero;fsp=zero;fwp=zero
  fb=zero;fc=zero

  if (ifirstder==1) then    ! Second-order central
     alfai= zero
     afi  = one/(two*d)
     bfi  = zero
  elseif(ifirstder==2) then ! Fourth-order central
     if (nrank==0) write(*,*)'Set of coefficients not ready yet'
     call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
  elseif(ifirstder==3) then ! Fourth-order compact
     if (nrank==0) write(*,*)'Set of coefficients not ready yet'
     call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
  elseif(ifirstder==4) then ! Sixth-order compact
     alfai= one/three
     afi  = (seven/nine)/d
     bfi  = (one/thirtysix)/d
  else
     if (nrank==0) then
        write(*,*) 'This is not an option. Please use ifirstder=1,2,3,4'
     endif
     call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
  endif

  if (ifirstder==1) then
     alfa1 = zero
     af1   = zero
     bf1   = zero
     cf1   = zero
     df1   = zero
     alfa2 = zero
     af2   = zero

     alfam = zero
     afm   = zero
     alfan = zero
     afn   = zero
     bfn   = zero
     cfn   = zero
     dfn   = zero
  else
     alfa1= two
     af1  =-(five/two)/d
     bf1  = (two)/d
     cf1  = (half)/d
     df1  = zero

     alfa2= one/four
     af2  = (three/four)/d
     alfan= two
     afn  =-(five/two)/d
     bfn  = (two)/d
     cfn  = (half)/d
     dfn  = zero
     alfam= one/four
     afm  = (three/four)/d
  endif

  if     (ncl1.eq.0) then !Periodic
     ff(1)   =alfai
     ff(2)   =alfai
     fc(1)   =two
     fc(2)   =one
     fb(1)   =alfai
     fb(2)   =alfai
  elseif (ncl1.eq.1) then !Free-slip
     ff(1)   =alfai+alfai
     ff(2)   =alfai
     fc(1)   =one
     fc(2)   =one
     fb(1)   =alfai
     fb(2)   =alfai
  elseif (ncl1.eq.2) then !Dirichlet
     ff(1)   =alfa1
     ff(2)   =alfa2
     fc(1)   =one
     fc(2)   =one
     fb(1)   =alfa2
     fb(2)   =alfai
  endif
  if     (ncln.eq.0) then !Periodic
     ff(n-2)=alfai
     ff(n-1)=alfai
     ff(n)  =zero
     fc(n-2)=one
     fc(n-1)=one
     fc(n  )=one+alfai*alfai
     fb(n-2)=alfai
     fb(n-1)=alfai
     fb(n  )=zero
  elseif (ncln.eq.1) then !Free-slip
     ff(n-2)=alfai
     ff(n-1)=alfai
     ff(n)  =zero
     fc(n-2)=one
     fc(n-1)=one
     fc(n  )=one
     fb(n-2)=alfai
     fb(n-1)=alfai+alfai
     fb(n  )=zero
  elseif (ncln.eq.2) then !Dirichlet
     ff(n-2)=alfai
     ff(n-1)=alfam
     ff(n)  =zero
     fc(n-2)=one
     fc(n-1)=one
     fc(n  )=one
     fb(n-2)=alfam
     fb(n-1)=alfan
     fb(n  )=zero
  endif
  do i=3,n-3
     ff(i)=alfai
     fc(i)=one
     fb(i)=alfai
  enddo

  do i=1,n
     ffp(i)=ff(i)
  enddo

  call prepare (fb,fc,ff ,fs ,fw ,n)

  if (ncl1.eq.1) then
     ffp(1)=zero
  endif
  if (ncln.eq.1) then
     fb(n-1)=zero
  endif

  call prepare (fb,fc,ffp,fsp,fwp,n)

  return
end subroutine first_derivative

!*******************************************************************
subroutine second_derivative(alsa1,as1,bs1,&
     cs1,ds1,alsa2,as2,alsan,asn,bsn,csn,dsn,alsam,&
     asm,alsa3,as3,bs3,alsat,ast,bst,&
     alsa4,as4,bs4,cs4,&
     alsatt,astt,bstt,cstt,&
     alsai,asi,bsi,csi,dsi,&
     sf,ss,sw,sfp,ssp,swp,d2,n,ncl1,ncln)
  !*******************************************************************

  use decomp_2d, only : mytype, nrank
  use param
  use MPI
  use variables, only : nu0nu,cnu

  implicit none

  real(mytype),intent(in) :: d2
  integer,intent(in) :: n,ncl1,ncln
  integer :: code,ierror
  real(mytype),dimension(n),intent(out) :: sf,ss,sw,sfp,ssp,swp
  real(mytype),intent(out) :: alsa1,as1,bs1,&
       cs1,ds1,alsa2,as2,alsan,asn,bsn,csn,dsn,alsam,&
       asm,alsa3,as3,bs3,alsat,ast,bst,&
       alsa4,as4,bs4,cs4,&
       alsatt,astt,bstt,cstt,&
       alsai,asi,bsi,csi,dsi
  integer :: i
  real(mytype),dimension(n) :: sb,sc
  real(mytype) :: xxnu,dpis3,kppkc,kppkm,xnpi2,xmpi2,den

  sf=zero;ss=zero;sw=zero;sfp=zero;ssp=zero;swp=zero

  ! Define coefficients based on the desired formal accuracy of the numerical schemes
  if (isecondder==1) then    ! Second-order central
     alsai=zero
     asi  =one/d2  !((six-nine*alsai)/four)/d2
     bsi  =zero !((-three+twentyfour*alsai)/five)/(four*d2)
     csi  =zero !((two-eleven*alsai)/twenty)/(nine*d2)
     dsi  =zero

     alsa4= alsai
     as4  = asi
     bs4  = bsi
     cs4  = csi

     alsatt = alsai
     astt = asi
     bstt = bsi
     cstt = csi
  elseif(isecondder==2) then ! Fourth-order central
     if (nrank==0) write(*,*)'Set of coefficients not ready yet'
     call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
  elseif(isecondder==3) then ! Fourth-order compact
     if (nrank==0) write(*,*)'Set of coefficients not ready yet'
     call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
  elseif(isecondder==4) then ! Sixth-order compact Lele style (no extra dissipation)
     alsai= two / eleven
     asi  = (twelve / eleven)/d2
     bsi  = (three / 44.0_mytype )/d2
     csi  = zero
     dsi = zero

     alsa4= alsai
     as4  = asi
     bs4  = bsi
     cs4  = csi

     alsatt = alsai
     astt = asi
     bstt = bsi
     cstt = csi
  elseif(isecondder==5) then ! Sixth-order Hyperviscous operator
     if(nrank==0) write(*,*) 'Using the hyperviscous operator with (nu_0/nu,c_nu) = ', '(', nu0nu,',', cnu,')'
     dpis3=two*pi/three
     kppkc=pi*pi*(one+nu0nu)
     kppkm=dpis3*dpis3*(one+cnu*nu0nu) !exp(-((pi-dpis3)/(zpthree*pi-dpis3))**two)/xxnu+dpis3*dpis3
     xnpi2=kppkc
     xmpi2=kppkm

     den = 405._mytype * xnpi2 - 640._mytype * xmpi2 + 144._mytype

     alsai = half - (320._mytype * xmpi2 - 1296._mytype) / den
     asi = -(4329._mytype * xnpi2 / eight - 32._mytype * xmpi2 - 140._mytype * xnpi2 * xmpi2 + 286._mytype) / den / d2
     bsi = (2115._mytype * xnpi2 - 1792._mytype * xmpi2 - 280._mytype * xnpi2 * xmpi2 + 1328._mytype) / den / (four * d2)
     csi = -(7695._mytype * xnpi2 / eight + 288._mytype * xmpi2 - 180._mytype * xnpi2 * xmpi2 - 2574._mytype) / den / (nine * d2)
     dsi = (198._mytype * xnpi2 + 128._mytype * xmpi2 - 40._mytype * xnpi2 * xmpi2 - 736._mytype) / den / (four**2 * d2)
  else
     if (nrank==0) then
        write(*,*) 'This is not an option.'
     endif
  endif

  ! Defined for the bounadies when dirichlet conditions are used
  alsa1= eleven
  as1  = (thirteen)/d2
  bs1  =-(twentyseven)/d2
  cs1  = (fifteen)/d2
  ds1  =-(one)/d2

  if (isecondder==1) then
     alsa2 = zero
     as2   = one / d2
  else
     alsa2= zpone
     as2  = (six/five)/d2
  endif

  alsa3= two/eleven
  as3  = (twelve/eleven)/d2
  bs3  = (three/fortyfour)/d2

  alsa4= two/eleven
  as4  = (twelve/eleven)/d2
  bs4  = (three/fortyfour)/d2
  cs4  = zero

  alsan= eleven
  asn  = (thirteen)/d2
  bsn  =-(twentyseven)/d2
  csn  = (fifteen)/d2
  dsn  =-(one)/d2

  if (isecondder==1) then
     alsam = zero
     asm   = one / d2
  else
     alsam= zpone
     asm  = (six/five)/d2
  endif

  alsat= two/eleven
  ast  = (twelve/eleven)/d2
  bst  = (three/fortyfour)/d2

  alsatt = two/eleven
  astt = (twelve/eleven)/d2
  bstt = (three/fortyfour)/d2
  cstt = zero

  if     (ncl1.eq.0) then !Periodic
     sf(1)   =alsai
     sf(2)   =alsai
     sf(3)   =alsai
     sf(4)   =alsai
     sc(1)   =two
     sc(2)   =one
     sc(3)   =one
     sc(4)   =one
     sb(1)   =alsai
     sb(2)   =alsai
     sb(3)   =alsai
     sb(4)   =alsai
  elseif (ncl1.eq.1) then !Free-slip
     sf(1)   =alsai+alsai
     sf(2)   =alsai
     sf(3)   =alsai
     sf(4)   =alsai
     sc(1)   =one
     sc(2)   =one
     sc(3)   =one
     sc(4)   =one
     sb(1)   =alsai
     sb(2)   =alsai
     sb(3)   =alsai
     sb(4)   =alsai
  elseif (ncl1.eq.2) then !Dirichlet
     sf(1)   =alsa1
     sf(2)   =alsa2
     sf(3)   =alsa3
     sf(4)   =alsa4
     sc(1)   =one
     sc(2)   =one
     sc(3)   =one
     sc(4)   =one
     sb(1)   =alsa2
     sb(2)   =alsa3
     sb(3)   =alsa4
     sb(4)   =alsai
  endif
  if     (ncln.eq.0) then !Periodic
     sf(n-4)=alsai
     sf(n-3)=alsai
     sf(n-2)=alsai
     sf(n-1)=alsai
     sf(n)  =zero
     sc(n-4)=one
     sc(n-3)=one
     sc(n-2)=one
     sc(n-1)=one
     sc(n  )=one+alsai*alsai
     sb(n-4)=alsai
     sb(n-3)=alsai
     sb(n-2)=alsai
     sb(n-1)=alsai
     sb(n  )=zero
  elseif (ncln.eq.1) then !Free-slip
     sf(n-4)=alsai
     sf(n-3)=alsai
     sf(n-2)=alsai
     sf(n-1)=alsai
     sf(n)  =zero
     sc(n-4)=one
     sc(n-3)=one
     sc(n-2)=one
     sc(n-1)=one
     sc(n  )=one
     sb(n-4)=alsai
     sb(n-3)=alsai
     sb(n-2)=alsai
     sb(n-1)=alsai+alsai
     sb(n  )=zero
  elseif (ncln.eq.2) then !Dirichlet
     sf(n-4)=alsai
     sf(n-3)=alsatt
     sf(n-2)=alsat
     sf(n-1)=alsam
     sf(n)  =zero
     sc(n-4)=one
     sc(n-3)=one
     sc(n-2)=one
     sc(n-1)=one
     sc(n  )=one
     sb(n-4)=alsatt
     sb(n-3)=alsat
     sb(n-2)=alsam
     sb(n-1)=alsan
     sb(n  )=zero
  endif
  do i=5,n-5
     sf(i)=alsai
     sc(i)=one
     sb(i)=alsai
  enddo

  do i=1,n
     sfp(i)=sf(i)
  enddo

  if (ncl1.eq.1) then
     sf (1)=zero
  endif

  call prepare (sb,sc,sf ,ss ,sw ,n)
  call prepare (sb,sc,sfp,ssp,swp,n)

  if (ncln.eq.1) then
     sb(n-1)=zero
     call prepare (sb,sc,sf ,ss ,sw ,n)
  endif

  return
end subroutine second_derivative

!*******************************************************************
!
subroutine interpolation(dx,nxm,nx,nclx1,nclxn,&
     alcaix6,acix6,bcix6,&
     ailcaix6,aicix6,bicix6,cicix6,dicix6,&
     cfx6,ccx6,cbx6,cfxp6,ciwxp6,csxp6,&
     cwxp6,csx6,cwx6,cifx6,cicx6,cisx6,&
     cibx6,cifxp6,cisxp6,ciwx6,&
     cfi6,cci6,cbi6,cfip6,csip6,cwip6,csi6,&
     cwi6,cifi6,cici6,cibi6,cifip6,&
     cisip6,ciwip6,cisi6,ciwi6)
  !
  !*******************************************************************

  use decomp_2d, only : mytype
  use param, only : zero, half, one, two, three, four, nine, ten
  use param, only : ipinter, ifirstder

  implicit none

  real(mytype),intent(in) :: dx
  integer,intent(in) :: nxm,nx,nclx1,nclxn
  real(mytype) :: alcaix6,acix6,bcix6
  real(mytype) :: ailcaix6,aicix6,bicix6,cicix6,dicix6
  real(mytype),dimension(nxm) :: cfx6,ccx6,cbx6,cfxp6,ciwxp6,csxp6,&
       cwxp6,csx6,cwx6,cifx6,cicx6,cisx6
  real(mytype),dimension(nxm) :: cibx6,cifxp6,cisxp6,ciwx6
  real(mytype),dimension(nx) :: cfi6,cci6,cbi6,cfip6,csip6,cwip6,csi6,&
       cwi6,cifi6,cici6,cibi6,cifip6
  real(mytype),dimension(nx) :: cisip6,ciwip6,cisi6,ciwi6

  integer :: i

  if (ifirstder==1) then
     alcaix6 = zero
     acix6   = one / dx
     bcix6   = zero
  else
     alcaix6=nine/62._mytype
     acix6=(63._mytype/62._mytype)/dx
     bcix6=(17._mytype/62._mytype)/three/dx
  endif

  cfx6(1)=alcaix6
  cfx6(2)=alcaix6
  cfx6(nxm-2)=alcaix6
  cfx6(nxm-1)=alcaix6
  cfx6(nxm)=zero
  if (nclx1==0) ccx6(1)=two
  if (nclx1==1) ccx6(1)=one + alcaix6
  if (nclx1==2) ccx6(1)=one + alcaix6
  ccx6(2)=one
  ccx6(nxm-2)=one
  ccx6(nxm-1)=one
  if (nclxn==0) ccx6(nxm)=one + alcaix6*alcaix6
  if (nclxn==1) ccx6(nxm)=one + alcaix6
  if (nclxn==2) ccx6(nxm)=one + alcaix6
  cbx6(1)=alcaix6
  cbx6(2)=alcaix6
  cbx6(nxm-2)=alcaix6
  cbx6(nxm-1)=alcaix6
  cbx6(nxm)=zero
  do i=3,nxm-3
     cfx6(i)=alcaix6
     ccx6(i)=one
     cbx6(i)=alcaix6
  enddo

  cfi6(1)=alcaix6 + alcaix6
  cfi6(2)=alcaix6
  cfi6(nx-2)=alcaix6
  cfi6(nx-1)=alcaix6
  cfi6(nx)=zero
  cci6(1)=one
  cci6(2)=one
  cci6(nx-2)=one
  cci6(nx-1)=one
  cci6(nx)=one
  cbi6(1)=alcaix6
  cbi6(2)=alcaix6
  cbi6(nx-2)=alcaix6
  cbi6(nx-1)=alcaix6 + alcaix6
  cbi6(nx)=zero
  do i=3,nx-3
     cfi6(i)=alcaix6
     cci6(i)=one
     cbi6(i)=alcaix6
  enddo

  if (ifirstder == 1) then
     ailcaix6 = zero
     aicix6   = half
     bicix6   = zero
     cicix6   = zero
     dicix6   = zero
  else if (ipinter.eq.1) then
     ailcaix6=three/ten
     aicix6=three/four
     bicix6=one/(two*ten)
     cicix6=zero
     dicix6=zero
  else if (ipinter.eq.2) then
     ailcaix6=0.461658_mytype

     dicix6=0.00293016_mytype
     aicix6=one/64._mytype *(75._mytype +70._mytype *ailcaix6-320._mytype *dicix6)
     bicix6=one/128._mytype *(126._mytype *ailcaix6-25._mytype +1152._mytype *dicix6)
     cicix6=one/128._mytype *(-ten*ailcaix6+three-640._mytype *dicix6)

     aicix6=aicix6/two
     bicix6=bicix6/two
     cicix6=cicix6/two
     dicix6=dicix6/two
  else if (ipinter.eq.3) then
     ailcaix6=0.49_mytype
     aicix6=one/128._mytype *(75._mytype +70._mytype*ailcaix6)
     bicix6=one/256._mytype *(126._mytype*ailcaix6-25._mytype)
     cicix6=one/256._mytype *(-ten*ailcaix6+three)
     dicix6=zero
  endif

  cifx6(1)=ailcaix6
  cifx6(2)=ailcaix6
  cifx6(nxm-2)=ailcaix6
  cifx6(nxm-1)=ailcaix6
  cifx6(nxm)=zero
  if (nclx1==0) cicx6(1)=two
  if (nclx1==1) cicx6(1)=one + ailcaix6
  if (nclx1==2) cicx6(1)=one + ailcaix6
  cicx6(2)=one
  cicx6(nxm-2)=one
  cicx6(nxm-1)=one
  if (nclxn==0) cicx6(nxm)=one + ailcaix6*ailcaix6
  if (nclxn==1) cicx6(nxm)=one + ailcaix6
  if (nclxn==2) cicx6(nxm)=one + ailcaix6
  cibx6(1)=ailcaix6
  cibx6(2)=ailcaix6
  cibx6(nxm-2)=ailcaix6
  cibx6(nxm-1)=ailcaix6
  cibx6(nxm)=zero
  do i=3,nxm-3
     cifx6(i)=ailcaix6
     cicx6(i)=one
     cibx6(i)=ailcaix6
  enddo
  cifi6(1)=ailcaix6 + ailcaix6
  cifi6(2)=ailcaix6
  cifi6(nx-2)=ailcaix6
  cifi6(nx-1)=ailcaix6
  cifi6(nx)=zero
  cici6(1)=one
  cici6(2)=one
  cici6(nx-2)=one
  cici6(nx-1)=one
  cici6(nx)=one
  cibi6(1)=ailcaix6
  cibi6(2)=ailcaix6
  cibi6(nx-2)=ailcaix6
  cibi6(nx-1)=ailcaix6 + ailcaix6
  cibi6(nx)=zero
  do i=3,nx-3
     cifi6(i)=ailcaix6
     cici6(i)=one
     cibi6(i)=ailcaix6
  enddo

  do i=1,nxm
     cfxp6(i)=cfx6(i)
     cifxp6(i)=cifx6(i)
  enddo
  do i=1,nx
     cifip6(i)=cifi6(i)
     cfip6(i)=cfi6(i)
  enddo
  cfxp6(1)=zero
  cfip6(1)=zero
  call prepare (cbx6,ccx6,cfx6 ,csx6 ,cwx6 ,nxm)
  call prepare (cbx6,ccx6,cfxp6,csxp6,cwxp6,nxm)
  call prepare (cibx6,cicx6,cifx6 ,cisx6 ,ciwx6 ,nxm)
  call prepare (cibx6,cicx6,cifxp6,cisxp6,ciwxp6,nxm)
  call prepare (cbi6,cci6,cfi6 ,csi6 ,cwi6 ,nx)
  call prepare (cbi6,cci6,cfip6,csip6,cwip6,nx)
  call prepare (cibi6,cici6,cifi6 ,cisi6 ,ciwi6 ,nx)
  call prepare (cibi6,cici6,cifip6,cisip6,ciwip6,nx)
  if (nclxn.eq.1) then
     cbx6(nxm-1)=zero
     cibx6(nxm)=zero
     cbi6(nx-1)=zero
     cibi6(nx)=zero
     call prepare (cbx6,ccx6,cfxp6,csxp6,cwxp6,nxm)
     call prepare (cibx6,cicx6,cifxp6,cisxp6,ciwxp6,nxm)
     call prepare (cbi6,cci6,cfip6,csip6,cwip6,nx)
     call prepare (cibi6,cici6,cifip6,cisip6,ciwip6,nx)
  endif
  if (nclxn.eq.2) then
     cbx6(nxm-1)=zero
     cibx6(nxm)=zero
     cbi6(nx-1)=zero
     cibi6(nx)=zero
     call prepare (cbx6,ccx6,cfxp6,csxp6,cwxp6,nxm)
     call prepare (cibx6,cicx6,cifxp6,cisxp6,ciwxp6,nxm)
     call prepare (cbi6,cci6,cfip6,csip6,cwip6,nx)
     call prepare (cibi6,cici6,cifip6,cisip6,ciwip6,nx)
  endif

  return
end subroutine interpolation
