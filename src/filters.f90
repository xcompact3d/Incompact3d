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

subroutine filter(af)
  USE param
  USE parfiX
  USE parfiY
  USE parfiZ
  USE variables
  USE var
  !=================================================
  ! Discrete low-pass filter according to
  !=================================================
  implicit none
  real(mytype),intent(in) :: af

#ifdef DEBG
  if (nrank  ==  0) print *,'# filter calculation start'
#endif

  ! Filter functions
  if (nclx1 == 0.and.nclxn == 0) filx => filx_00
  if (nclx1 == 1.and.nclxn == 1) filx => filx_11
  if (nclx1 == 1.and.nclxn == 2) filx => filx_12
  if (nclx1 == 2.and.nclxn == 1) filx => filx_21
  if (nclx1 == 2.and.nclxn == 2) filx => filx_22
  !
  if (ncly1 == 0.and.nclyn == 0) fily => fily_00
  if (ncly1 == 1.and.nclyn == 1) fily => fily_11
  if (ncly1 == 1.and.nclyn == 2) fily => fily_12
  if (ncly1 == 2.and.nclyn == 1) fily => fily_21
  if (ncly1 == 2.and.nclyn == 2) fily => fily_22
  !
  if (nclz1 == 0.and.nclzn == 0) filz => filz_00
  if (nclz1 == 1.and.nclzn == 1) filz => filz_11
  if (nclz1 == 1.and.nclzn == 2) filz => filz_12
  if (nclz1 == 2.and.nclzn == 1) filz => filz_21
  if (nclz1 == 2.and.nclzn == 2) filz => filz_22

  ! Set coefficients for x-direction filter
  call set_filter_coefficients(af,fial1x,fia1x,fib1x,fic1x,fid1x,fial2x,fia2x,fib2x,fic2x,fid2x,fial3x,fia3x,fib3x,fic3x,fid3x,fie3x,fif3x,&
       fialnx,fianx,fibnx,ficnx,fidnx,fialmx,fiamx,fibmx,ficmx,fidmx,fialpx,fiapx,fibpx,ficpx,fidpx,fiepx,fifpx,&
       fialix,fiaix,fibix,ficix,fidix,fiffx,fifsx,fifwx,fiffxp,fifsxp,fifwxp,nx,nclx1,nclxn)
  ! Set coefficients for y-direction filter
  call set_filter_coefficients(af,fial1y,fia1y,fib1y,fic1y,fid1y,fial2y,fia2y,fib2y,fic2y,fid2y,fial3y,fia3y,fib3y,fic3y,fid3y,fie3y,fif3y,&
       fialny,fiany,fibny,ficny,fidny,fialmy,fiamy,fibmy,ficmy,fidmy,fialpy,fiapy,fibpy,ficpy,fidpy,fiepy,fifpy,&
       fialjy,fiajy,fibjy,ficjy,fidjy,fiffy,fifsy,fifwy,fiffyp,fifsyp,fifwyp,ny,ncly1,nclyn)
  ! Set coefficients for z-direction filter
  call set_filter_coefficients(af,fial1z,fia1z,fib1z,fic1z,fid1z,fial2z,fia2z,fib2z,fic2z,fid2z,fial3z,fia3z,fib3z,fic3z,fid3z,fie3z,fif3z,&
       fialnz,fianz,fibnz,ficnz,fidnz,fialmz,fiamz,fibmz,ficmz,fidmz,fialpz,fiapz,fibpz,ficpz,fidpz,fiepz,fifpz,&
       fialkz,fiakz,fibkz,fickz,fidkz,fiffz,fifsz,fifwz,fiffzp,fifszp,fifwzp,nz,nclz1,nclzn)
#ifdef DEBG
  if (nrank  ==  0) print *,'# filter calculation end'
#endif

  return

end subroutine filter


subroutine set_filter_coefficients(af,alfa1,a1,b1,c1,d1,alfa2,a2,b2,c2,d2,alfa3,a3,b3,c3,d3,e3,f3,&
     alfan,an,bn,cn,dn,alfam,am,bm,cm,dm,alfap,ap,bp,cp,dp,ep,fp,&
     alfai,ai,bi,ci,di,ff,fs,fw,ffp,fsp,fwp,n,ncl1,ncln)

  use decomp_2d, only : mytype
  use param

  implicit none

  real(mytype),intent(in) :: af
  integer,intent(in) :: n,ncl1,ncln
  real(mytype),dimension(n),intent(out) :: ff,fs,fw,ffp,fsp,fwp
  real(mytype),intent(out) :: alfa1,a1,b1,c1,d1,alfa2,a2,b2,c2,d2,alfa3,a3,b3,c3,d3,e3,f3,&
       alfan,an,bn,cn,dn,alfam,am,bm,cm,dm,alfap,ap,bp,cp,dp,ep,fp,&
       alfai,ai,bi,ci,di
  integer :: i
  real(mytype),dimension(n) :: fb,fc

  ! Set the coefficient for the discrete filter following
  ! the tridiagonal filtering of Motheau and Abraham, JCP 2016
  ! Filter should be -0.5<filax<0.5

  ! General Case (entire points)
  ! alpha*fhat(i-1)+fhat(i)+alpha*fhat(i+1)=af(i)+b/2*[f(i+1)+f(i-1)] + ...

  ! Coefficients are calculated according to the report of Gaitonde & Visbal, 1998,
  ! "High-order schemes for Navier-Stokes equations: Algorithm and implementation into FDL3DI"


  alfai=af                                       ! alpha_f
  !Interior points
  ai=(eleven + ten*af)/sixteen                   ! a
  bi=half*(fifteen +thirtyfour*af)/thirtytwo     ! b/2
  ci=half*(-three + six*af)/sixteen              ! c/2
  di=half*(one - two*af)/thirtytwo               ! d/2
  ! Explicit third/fifth-order filters near the boundaries!
  !Boundary point 1 (no-filtering)
  alfa1=zero
  a1=one                           ! a1=7./8.+af/8.! a1/2
  b1=zero                          ! b1=3./8.+5.*af/8.
  c1=zero                          ! c1=-3./8.+3./8.*af
  d1=zero                          ! d1=1./8.-1./8.*af
  !Boundary point 2 (Third order)
  alfa2=af
  a2=one/eight+three/four*af            ! a2
  b2=five/eight+three/four*af           ! b2
  c2=three/eight+af/four                ! c2
  d2=-one/eight+af/four                 ! d2
  !Boundary point 3 (Fifth order)
  alfa3=af
  a3= -one/thirtytwo+af/sixteen         ! a3
  b3= five/thirtytwo+eleven/sixteen*af  ! b3
  c3= eleven/sixteen+five*af/eight      ! c3
  d3= five/sixteen+three*af/eight       ! d3
  e3=-five/thirtytwo+five*af/sixteen    ! e3
  f3= one/thirtytwo-af/sixteen          ! f3
  !Boundary point n (no-filtering)
  alfan=zero
  an=one                                !an = 7./8.+af/8.! a1/2
  bn=zero                               !bn = 3./8.+5.*af/8.
  cn=zero                               !cn =-3./8.+3./8.*af
  dn=zero                               !dn = 1./8.-1./8.*af
  !Boundary point 2 (Third order)
  alfam=af
  am=one/eight+three/four*af            ! am
  bm=five/eight+three/four*af           ! bm
  cm=three/eight+af/four                ! cm
  dm=-one/eight+af/four                 ! dm
  !Boundary point 3 (Fifth order)
  alfap=af
  ap=-one/thirtytwo+af/sixteen          ! ap
  bp= five/thirtytwo+eleven/sixteen*af  ! bp
  cp= eleven/sixteen+five*af/eight      ! cp
  dp= five/sixteen+three*af/eight       ! dp
  ep=-five/thirtytwo+five*af/sixteen    ! ep
  fp= one/thirtytwo-af/sixteen          ! fp

  ff=zero;fs=zero;fw=zero;ffp=zero;fsp=zero;fwp=zero
  fb=zero;fc=zero

  if     (ncl1 == 0) then !Periodic
     ff(1)   =alfai
     ff(2)   =alfai
     fc(1)   =two
     fc(2)   =one
     fb(1)   =alfai
     fb(2)   =alfai
  elseif (ncl1 == 1) then !Free-slip
     ff(1)   =alfai+alfai
     ff(2)   =alfai
     fc(1)   =one
     fc(2)   =one
     fb(1)   =alfai
     fb(2)   =alfai
  elseif (ncl1 == 2) then !Dirichlet
     ff(1)   =alfa1
     ff(2)   =alfa2
     fc(1)   =one
     fc(2)   =one
     fb(1)   =alfa2
     fb(2)   =alfai
  endif
  if (ncln == 0) then !Periodic
     ff(n-2)=alfai
     ff(n-1)=alfai
     ff(n)  =zero
     fc(n-2)=one
     fc(n-1)=one
     fc(n  )=one+alfai*alfai
     fb(n-2)=alfai
     fb(n-1)=alfai
     fb(n  )=zero
  elseif (ncln == 1) then !Free-slip
     ff(n-2)=alfai
     ff(n-1)=alfai
     ff(n)  =zero
     fc(n-2)=one
     fc(n-1)=one
     fc(n  )=one
     fb(n-2)=alfai
     fb(n-1)=alfai+alfai
     fb(n  )=zero
  elseif (ncln == 2) then !Dirichlet
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

  call prepare (fb,fc,ffp ,fsp ,fwp ,n)

  if (ncl1 == 1) then
     ff(1)=zero
  endif
  if (ncln == 1) then
     fb(n-1)=zero
  endif

  call prepare (fb,fc,ff,fs,fw,n)

  return

end subroutine set_filter_coefficients

subroutine filx_00(tx,ux,rx,fisx,fiffx,fifsx,fifwx,nx,ny,nz,npaire,lind) 

  USE param
  use parfiX
  use ibm, only : lagpolx, cubsplx

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nz, npaire
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx, rx
  real(mytype), intent(inout), dimension(nx,ny,nz) :: ux
  real(mytype), intent(out), dimension(ny,nz) :: fisx
  real(mytype), intent(in), dimension(nx) :: fiffx, fifsx, fifwx
  real(mytype), intent(in) :: lind

  ! Local variables
  integer :: i, j, k

  if (iibm == 2) call lagpolx(ux)
  if (iibm == 3) call cubsplx(ux,lind)

  do concurrent (k=1:nz)
     do concurrent (j=1:ny)

        ! Compute r.h.s.
        tx(1,j,k) = fiaix*ux(1,j,k) + fibix*(ux(2,j,k)+ux(nx,j,k)) &
                  + ficix*(ux(3,j,k)+ux(nx-1,j,k)) &
                  + fidix*(ux(4,j,k)+ux(nx-2,j,k))
        tx(2,j,k) = fiaix*ux(2,j,k) + fibix*(ux(3,j,k)+ux(1,j,k)) &
                  + ficix*(ux(4,j,k)+ux(nx,j,k)) &
                  + fidix*(ux(5,j,k)+ux(nx-1,j,k))
        tx(3,j,k) = fiaix*ux(3,j,k) + fibix*(ux(4,j,k)+ux(2,j,k)) &
                  + ficix*(ux(5,j,k)+ux(1,j,k)) &
                  + fidix*(ux(6,j,k)+ux(nx,j,k))
        do concurrent (i=4:nx-3)
           tx(i,j,k) = fiaix*ux(i,j,k) + fibix*(ux(i+1,j,k)+ux(i-1,j,k)) &
                     + ficix*(ux(i+2,j,k)+ux(i-2,j,k)) &
                     + fidix*(ux(i+3,j,k)+ux(i-3,j,k))
        enddo
        tx(nx-2,j,k) = fiaix*ux(nx-2,j,k) + fibix*(ux(nx-3,j,k)+ux(nx-1,j,k)) &
                     + ficix*(ux(nx-4,j,k)+ux(nx,j,k)) &
                     + fidix*(ux(nx-5,j,k)+ux(1,j,k))
        tx(nx-1,j,k) = fiaix*ux(nx-1,j,k) + fibix*(ux(nx-2,j,k)+ux(nx,j,k)) &
                     + ficix*(ux(nx-3,j,k)+ux(1,j,k))&
                     + fidix*(ux(nx-4,j,k)+ux(2,j,k))
        tx(nx  ,j,k) = fiaix*ux(nx,j,k) + fibix*(ux(nx-1,j,k)+ux(1,j,k)) &
                     + ficix*(ux(nx-2,j,k)+ux(2,j,k)) &
                     + fidix*(ux(nx-3,j,k)+ux(3,j,k))
        rx(1,j,k) = -one
        do concurrent (i=2:nx-1)
           rx(i,j,k) = zero
        enddo
        rx(nx,j,k) = fialix

        ! Solve tri-diagonal system
        do i = 2, nx
           tx(i,j,k) = tx(i,j,k) - tx(i-1,j,k)*fifsx(i)
           rx(i,j,k) = rx(i,j,k) - rx(i-1,j,k)*fifsx(i)
        enddo
        tx(nx,j,k) = tx(nx,j,k) * fifwx(nx)
        rx(nx,j,k) = rx(nx,j,k) * fifwx(nx)
        do i=nx-1,1,-1
           tx(i,j,k) = (tx(i,j,k)-fiffx(i)*tx(i+1,j,k)) * fifwx(i)
           rx(i,j,k) = (rx(i,j,k)-fiffx(i)*rx(i+1,j,k)) * fifwx(i)
        enddo
        fisx(j,k) = (    tx(1,j,k)-fialix*tx(nx,j,k)) &
                  / (one+rx(1,j,k)-fialix*rx(nx,j,k))
        do concurrent (i=1:nx)
           tx(i,j,k) = tx(i,j,k) - fisx(j,k)*rx(i,j,k)
        enddo
     enddo
  enddo

end subroutine filx_00

subroutine filx_ij(tx,ux,fisx,fiffx,fifsx,fifwx,nx,ny,nz,npaire,lind,ncl1,ncln)

  USE param
  use parfiX
  use ibm, only : lagpolx, cubsplx

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nz, npaire, ncl1, ncln
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(inout), dimension(nx,ny,nz) :: ux
  real(mytype), intent(out), dimension(ny,nz) :: fisx
  real(mytype), intent(in), dimension(nx) :: fiffx,fifsx,fifwx
  real(mytype), intent(in) :: lind

  ! Local variables
  integer :: i, j, k

  if (iibm == 2) call lagpolx(ux)
  if (iibm == 3) call cubsplx(ux,lind)

  do concurrent (k=1:nz)
     do concurrent (j=1:ny)

        ! Compute r.h.s.
        if (ncl1 == 1) then
           if (npaire == 1) then
              tx(1,j,k) = fiaix*ux(1,j,k) + fibix*(ux(2,j,k)+ux(2,j,k)) &
                        + ficix*(ux(3,j,k)+ux(3,j,k)) &
                        + fidix*(ux(4,j,k)+ux(4,j,k))
              tx(2,j,k) = fiaix*ux(2,j,k) + fibix*(ux(3,j,k)+ux(1,j,k)) &
                        + ficix*(ux(4,j,k)+ux(2,j,k)) &
                        + fidix*(ux(5,j,k)+ux(3,j,k))
              tx(3,j,k) = fiaix*ux(3,j,k) + fibix*(ux(4,j,k)+ux(2,j,k)) &
                        + ficix*(ux(5,j,k)+ux(1,j,k)) &
                        + fidix*(ux(6,j,k)+ux(2,j,k))
           else
              tx(1,j,k) = zero
              tx(2,j,k) = fiaix*ux(2,j,k) + fibix*(ux(3,j,k)+ux(1,j,k)) &
                        + ficix*(ux(4,j,k)-ux(2,j,k)) &
                        + fidix*(ux(5,j,k)-ux(3,j,k))
              tx(3,j,k) = fiaix*ux(3,j,k) + fibix*(ux(4,j,k)+ux(2,j,k)) &
                        + ficix*(ux(5,j,k)+ux(1,j,k)) &
                        + fidix*(ux(6,j,k)-ux(2,j,k))
           endif
        else
           tx(1,j,k) = ux(1,j,k)
           tx(2,j,k) = fia2x*ux(1,j,k) + fib2x*ux(2,j,k)+fic2x*ux(3,j,k) &
                     + fid2x*ux(4,j,k)
           tx(3,j,k) = fia3x*ux(1,j,k) + fib3x*ux(2,j,k)+fic3x*ux(3,j,k) &
                     + fid3x*ux(4,j,k) + fie3x*ux(5,j,k)+fif3x*ux(6,j,k)
        endif
        do concurrent (i=4:nx-3)
           tx(i,j,k) = fiaix*ux(i,j,k) + fibix*(ux(i+1,j,k)+ux(i-1,j,k)) &
                     + ficix*(ux(i+2,j,k)+ux(i-2,j,k)) &
                     + fidix*(ux(i+3,j,k)+ux(i-3,j,k))
        enddo
        if (ncln == 1) then
           if (npaire == 1) then
              tx(nx-2,j,k) = fiaix*ux(nx-2,j,k) + fibix*(ux(nx-1,j,k)+ux(nx-3,j,k)) &
                           + ficix*(ux(  nx,j,k)+ux(nx-4,j,k)) &
                           + fidix*(ux(nx-1,j,k)+ux(nx-5,j,k))
              tx(nx-1,j,k) = fiaix*ux(nx-1,j,k) + fibix*(ux(  nx,j,k)+ux(nx-2,j,k)) &
                           + ficix*(ux(nx-1,j,k)+ux(nx-3,j,k)) &
                           + fidix*(ux(nx-2,j,k)+ux(nx-4,j,k))
              tx(nx,j,k)   = fiaix*ux(nx,j,k) + fibix*(ux(nx-1,j,k)+ux(nx-1,j,k)) &
                           + ficix*(ux(nx-2,j,k)+ux(nx-2,j,k)) &
                           + fidix*(ux(nx-3,j,k)+ux(nx-3,j,k))
           else
              tx(nx-2,j,k) = fiaix*ux(nx-2,j,k) + fibix*( ux(nx-1,j,k)+ux(nx-3,j,k)) &
                           + ficix*( ux(nx  ,j,k)+ux(nx-4,j,k)) &
                           + fidix*(-ux(nx-1,j,k)+ux(nx-5,j,k))
              tx(nx-1,j,k) = fiaix*ux(nx-1,j,k) + fibix*( ux(nx  ,j,k)+ux(nx-2,j,k)) &
                           + ficix*(-ux(nx-1,j,k)+ux(nx-3,j,k)) &
                           + fidix*(-ux(nx-2,j,k)+ux(nx-4,j,k))
              tx(nx,  j,k) = zero
           endif
        else
           tx(nx-2,j,k) = fiapx*ux(nx,j,k) + fibpx*ux(nx-1,j,k) + ficpx*ux(nx-2,j,k) &
                        + fidpx*ux(nx-3,j,k) + fiepx*ux(nx-4,j,k) + fifpx*ux(nx-5,j,k)
           tx(nx-1,j,k) = fiamx*ux(nx,j,k) + fibmx*ux(nx-1,j,k) + ficmx*ux(nx-2,j,k) &
                        + fidmx*ux(nx-3,j,k)
           tx(nx,j,k) = ux(nx,j,k)
        endif

        ! Solve tri-diagonal system
        do i=2,nx
           tx(i,j,k) = tx(i,j,k) - tx(i-1,j,k)*fifsx(i)
        enddo
        tx(nx,j,k) = tx(nx,j,k) * fifwx(nx)
        do i=nx-1,1,-1
           tx(i,j,k) = (tx(i,j,k)-fiffx(i)*tx(i+1,j,k)) * fifwx(i)
        enddo
     enddo
  enddo

end subroutine filx_ij

subroutine filx_11(tx,ux,rx,fisx,fiffx,fifsx,fifwx,nx,ny,nz,npaire,lind) 

  USE param

  implicit none

  integer, intent(in) :: nx, ny, nz, npaire
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx, rx
  real(mytype), intent(inout), dimension(nx,ny,nz) :: ux
  real(mytype), intent(out), dimension(ny,nz) :: fisx
  real(mytype), intent(in), dimension(nx) :: fiffx,fifsx,fifwx
  real(mytype), intent(in) :: lind

  call filx_ij(tx,ux,fisx,fiffx,fifsx,fifwx,nx,ny,nz,npaire,lind,1,1)

end subroutine filx_11

subroutine filx_12(tx,ux,rx,fisx,fiffx,fifsx,fifwx,nx,ny,nz,npaire,lind) 

  USE param

  implicit none

  integer, intent(in) :: nx, ny, nz, npaire
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx, rx
  real(mytype), intent(inout), dimension(nx,ny,nz) :: ux
  real(mytype), intent(out), dimension(ny,nz) :: fisx
  real(mytype), intent(in), dimension(nx) :: fiffx,fifsx,fifwx
  real(mytype), intent(in) :: lind

  call filx_ij(tx,ux,fisx,fiffx,fifsx,fifwx,nx,ny,nz,npaire,lind,1,2)

end subroutine filx_12

subroutine filx_21(tx,ux,rx,fisx,fiffx,fifsx,fifwx,nx,ny,nz,npaire,lind) 

  USE param

  implicit none

  integer, intent(in) :: nx, ny, nz, npaire
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx, rx
  real(mytype), intent(inout), dimension(nx,ny,nz) :: ux
  real(mytype), intent(out), dimension(ny,nz) :: fisx
  real(mytype), intent(in), dimension(nx) :: fiffx,fifsx,fifwx
  real(mytype), intent(in) :: lind

  call filx_ij(tx,ux,fisx,fiffx,fifsx,fifwx,nx,ny,nz,npaire,lind,2,1)

end subroutine filx_21


subroutine filx_22(tx,ux,rx,fisx,fiffx,fifsx,fifwx,nx,ny,nz,npaire,lind) 

  USE param

  implicit none

  integer, intent(in) :: nx, ny, nz, npaire
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx, rx
  real(mytype), intent(inout), dimension(nx,ny,nz) :: ux
  real(mytype), intent(out), dimension(ny,nz) :: fisx
  real(mytype), intent(in), dimension(nx) :: fiffx,fifsx,fifwx
  real(mytype), intent(in) :: lind

  call filx_ij(tx,ux,fisx,fiffx,fifsx,fifwx,nx,ny,nz,npaire,lind,2,2)

end subroutine filx_22

subroutine fily_00(ty,uy,ry,fisy,fiffy,fifsy,fifwy,nx,ny,nz,npaire,lind) 

  USE param
  use parfiY
  use ibm, only : lagpoly, cubsply

  implicit none

  integer, intent(in) :: nx, ny, nz, npaire
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty, ry
  real(mytype), intent(inout), dimension(nx,ny,nz) :: uy
  real(mytype), intent(out), dimension(nx,nz)  :: fisy
  real(mytype), intent(in), dimension(ny) :: fiffy,fifsy,fifwy
  real(mytype), intent(in) :: lind

  integer :: i, j, k

  if (iibm == 2) call lagpoly(uy)
  if (iibm == 3) call cubsply(uy,lind)

  do concurrent (k=1:nz)

     ! Solve r.h.s.
     do concurrent (i=1:nx)
        ty(i,1,k) = fiajy*uy(i,1,k) + fibjy*(uy(i,2,k)+uy(i,ny,k)) &
                  + ficjy*(uy(i,3,k)+uy(i,ny-1,k)) &
                  + fidjy*(uy(i,4,k)+uy(i,ny-2,k))
     enddo
     do concurrent (i=1:nx)
        ty(i,2,k) = fiajy*uy(i,2,k) + fibjy*(uy(i,3,k)+uy(i,1,k)) &
                  + ficjy*(uy(i,4,k)+uy(i,ny,k)) &
                  + fidjy*(uy(i,5,k)+uy(i,ny-1,k))
     enddo
     do concurrent (i=1:nx)
        ty(i,3,k) = fiajy*uy(i,3,k) + fibjy*(uy(i,4,k)+uy(i,2,k)) &
                  + ficjy*(uy(i,5,k)+uy(i,1,k)) &
                  + fidjy*(uy(i,6,k)+uy(i,ny,k))
     enddo
     do concurrent (j=4:ny-3)
        do concurrent (i=1:nx)
           ty(i,j,k) = fiajy*uy(i,j,k) + fibjy*(uy(i,j+1,k)+uy(i,j-1,k)) &
                     + ficjy*(uy(i,j+2,k)+uy(i,j-2,k)) &
                     + fidjy*(uy(i,j+3,k)+uy(i,j-3,k))
        enddo
     enddo
     do concurrent (i=1:nx)
        ty(i,ny-2,k) = fiajy*uy(i,ny-2,k) + fibjy*(uy(i,ny-3,k)+uy(i,ny-1,k)) &
                     + ficjy*(uy(i,ny-4,k)+uy(i,ny,k)) &
                     + fidjy*(uy(i,ny-5,k)+uy(i,1,k))
     enddo
     do concurrent (i=1:nx)
        ty(i,ny-1,k) = fiajy*uy(i,ny-1,k) + fibjy*(uy(i,ny-2,k)+uy(i,ny,k)) &
                     + ficjy*(uy(i,ny-3,k)+uy(i,1,k)) &
                     + fidjy*(uy(i,ny-4,k)+uy(i,2,k))
     enddo
     do concurrent (i=1:nx)
        ty(i,ny,k) = fiajy*uy(i,ny,k) + fibjy*(uy(i,ny-1,k)+uy(i,1,k)) &
                   + ficjy*(uy(i,ny-2,k)+uy(i,2,k)) &
                   + fidjy*(uy(i,ny-3,k)+uy(i,3,k))
     enddo
     do concurrent (i=1:nx)
        ry(i,1,k) = -one
     enddo
     do concurrent (j=2:ny-1)
        do concurrent (i=1:nx)
           ry(i,j,k) = zero
        enddo
     enddo
     do concurrent (i=1:nx)
        ry(i,ny,k) = fialjy
     enddo

     ! Solve tri-diagonal system
     do j=2,ny
        do concurrent (i=1:nx)
           ty(i,j,k) = ty(i,j,k) - ty(i,j-1,k)*fifsy(j)
           ry(i,j,k) = ry(i,j,k) - ry(i,j-1,k)*fifsy(j)
        enddo
     enddo
     do concurrent (i=1:nx)
        ty(i,ny,k) = ty(i,ny,k) * fifwy(ny)
     enddo
     do concurrent (i=1:nx)
        ry(i,ny,k) = ry(i,ny,k) * fifwy(ny)
     enddo
     do j=ny-1,1,-1
        do concurrent (i=1:nx)
           ty(i,j,k) = (ty(i,j,k)-fiffy(j)*ty(i,j+1,k)) * fifwy(j)
           ry(i,j,k) = (ry(i,j,k)-fiffy(j)*ry(i,j+1,k)) * fifwy(j)
        enddo
     enddo
     fisy(i,k) = (    ty(i,1,k)-fialjy*ty(i,ny,k)) &
               / (one+ry(i,1,k)-fialjy*ry(i,ny,k))
     do concurrent (j=1:ny)
        do concurrent (i=1:nx)
           ty(i,j,k) = ty(i,j,k) - fisy(i,k)*ry(i,j,k)
        enddo
     enddo

  enddo

end subroutine fily_00

subroutine fily_ij(ty,uy,fisy,fiffy,fifsy,fifwy,nx,ny,nz,npaire,lind,ncl1,ncln)
  !
  !********************************************************************

  USE param
  use parfiY
  use ibm, only : lagpoly, cubsply

  implicit none

  integer, intent(in) :: nx, ny, nz, npaire, ncl1, ncln
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(inout), dimension(nx,ny,nz) :: uy
  real(mytype), intent(out), dimension(nx,nz)  :: fisy
  real(mytype), intent(in), dimension(ny) :: fiffy,fifsy,fifwy
  real(mytype), intent(in) :: lind

  integer :: i, j, k

  if (iibm == 2) call lagpoly(uy)
  if (iibm == 3) call cubsply(uy,lind)

  do concurrent (k=1:nz)

     ! Solve r.h.s.
     if (ncl1 == 1) then
        if (npaire==1) then
           do concurrent (i=1:nx)
              ty(i,1,k) = fiajy*uy(i,1,k) + fibjy*(uy(i,2,k)+uy(i,2,k)) &
                        + ficjy*(uy(i,3,k)+uy(i,3,k)) &
                        + fidjy*(uy(i,4,k)+uy(i,4,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,2,k) = fiajy*uy(i,2,k) + fibjy*(uy(i,3,k)+uy(i,1,k)) &
                        + ficjy*(uy(i,4,k)+uy(i,2,k)) &
                        + fidjy*(uy(i,5,k)+uy(i,3,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,3,k) = fiajy*uy(i,3,k) + fibjy*(uy(i,4,k)+uy(i,2,k)) &
                        + ficjy*(uy(i,5,k)+uy(i,1,k)) &
                        + fidjy*(uy(i,6,k)+uy(i,2,k))
           enddo
        else
           do concurrent (i=1:nx)
              ty(i,1,k) = zero !fiajy*uy(i,1,k)
           enddo
           do concurrent (i=1:nx)
              ty(i,2,k) = fiajy*uy(i,2,k) + fibjy*(uy(i,3,k)+uy(i,1,k)) &
                        + ficjy*(uy(i,4,k)-uy(i,2,k)) &
                        + fidjy*(uy(i,5,k)-uy(i,3,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,3,k) = fiajy*uy(i,3,k) + fibjy*(uy(i,4,k)+uy(i,2,k)) &
                        + ficjy*(uy(i,5,k)+uy(i,1,k)) &
                        + fidjy*(uy(i,6,k)-uy(i,2,k))
           enddo
        endif
     else
     endif
     do concurrent (j=4:ny-3)
        do concurrent (i=1:nx)
           ty(i,j,k) = fiajy*uy(i,j,k) + fibjy*(uy(i,j+1,k)+uy(i,j-1,k)) &
                     + ficjy*(uy(i,j+2,k)+uy(i,j-2,k)) &
                     + fidjy*(uy(i,j+3,k)+uy(i,j-3,k))
        enddo
     enddo
     if (ncln == 1) then
        if (npaire == 1) then
           do concurrent (i=1:nx)
              ty(i,ny-2,k) = fiajy*uy(i,ny-2,k)+fibjy*(uy(i,ny-1,k)+uy(i,ny-3,k)) &
                           + ficjy*(uy(i,ny,k)+uy(i,ny-4,k)) &
                           + fidjy*(uy(i,ny-1,k)+uy(i,ny-5,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,ny-1,k) = fiajy*uy(i,ny-1,k)+fibjy*(uy(i,ny,k)  +uy(i,ny-2,k)) &
                           + ficjy*(uy(i,ny-1,k)+uy(i,ny-3,k)) &
                           + fidjy*(uy(i,ny-2,k)+uy(i,ny-4,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,ny  ,k) = fiajy*uy(i,ny,k) + fibjy*(uy(i,ny-1,k)+uy(i,ny-1,k)) &
                           + ficjy*(uy(i,ny-2,k)+uy(i,ny-2,k)) &
                           + fidjy*(uy(i,ny-3,k)+uy(i,ny-3,k))
           enddo
        else
           do concurrent (i=1:nx)
              ty(i,ny-2,k) = fiajy*uy(i,ny-2,k) + fibjy*(uy(i,ny-1,k)+uy(i,ny-3,k)) &
                           + ficjy*(uy(i,ny,k)+uy(i,ny-4,k)) &
                           + fidjy*(-uy(i,ny-1,k)+uy(i,ny-5,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,ny-1,k) = fiajy*uy(i,ny-1,k) + fibjy*(uy(i,ny,k)+uy(i,ny-2,k)) &
                           + ficjy*(-uy(i,ny-1,k)+uy(i,ny-3,k)) &
                           + fidjy*(-uy(i,ny-2,k)+uy(i,ny-4,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,ny,k) = zero !fiajy*uy(i,ny,k)
           enddo
        endif
     else
     endif

     ! Solve tri-diagonal system
     do j=2,ny
        do concurrent (i=1:nx)     
           ty(i,j,k) = ty(i,j,k) - ty(i,j-1,k)*fifsy(j)
        enddo
     enddo
     do concurrent (i=1:nx)
        ty(i,ny,k) = ty(i,ny,k) * fifwy(ny)
     enddo
     do j=ny-1,1,-1
        do concurrent (i=1:nx)
           ty(i,j,k) = (ty(i,j,k)-fiffy(j)*ty(i,j+1,k)) * fifwy(j)
        enddo
     enddo

  enddo

end subroutine fily_ij
!********************************************************************
!
subroutine fily_11(ty,uy,ry,fisy,fiffy,fifsy,fifwy,nx,ny,nz,npaire,lind) 
  !
  !********************************************************************

  USE param

  implicit none

  integer, intent(in) :: nx, ny, nz, npaire
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty, ry
  real(mytype), intent(inout), dimension(nx,ny,nz) :: uy
  real(mytype), intent(out), dimension(nx,nz)  :: fisy
  real(mytype), intent(in), dimension(ny) :: fiffy,fifsy,fifwy
  real(mytype), intent(in) :: lind

  call fily_ij(ty,uy,fisy,fiffy,fifsy,fifwy,nx,ny,nz,npaire,lind,1,1)

end subroutine fily_11


subroutine fily_12(ty,uy,ry,fisy,fiffy,fifsy,fifwy,nx,ny,nz,npaire,lind) 

  USE param

  implicit none

  integer, intent(in) :: nx, ny, nz, npaire
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty, ry
  real(mytype), intent(inout), dimension(nx,ny,nz) :: uy
  real(mytype), intent(out), dimension(nx,nz)  :: fisy
  real(mytype), intent(in), dimension(ny) :: fiffy,fifsy,fifwy
  real(mytype), intent(in) :: lind

  call fily_ij(ty,uy,fisy,fiffy,fifsy,fifwy,nx,ny,nz,npaire,lind,1,2)

end subroutine fily_12


subroutine fily_21(ty,uy,ry,fisy,fiffy,fifsy,fifwy,nx,ny,nz,npaire,lind) 

  USE param

  implicit none

  integer, intent(in) :: nx, ny, nz, npaire
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty, ry
  real(mytype), intent(inout), dimension(nx,ny,nz) :: uy
  real(mytype), intent(out), dimension(nx,nz)  :: fisy
  real(mytype), intent(in), dimension(ny) :: fiffy,fifsy,fifwy
  real(mytype), intent(in) :: lind

  call fily_ij(ty,uy,fisy,fiffy,fifsy,fifwy,nx,ny,nz,npaire,lind,2,1)

end subroutine fily_21

subroutine fily_22(ty,uy,ry,fisy,fiffy,fifsy,fifwy,nx,ny,nz,npaire,lind) 

  USE param

  implicit none

  integer, intent(in) :: nx, ny, nz, npaire
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty, ry
  real(mytype), intent(inout), dimension(nx,ny,nz) :: uy
  real(mytype), intent(out), dimension(nx,nz)  :: fisy
  real(mytype), intent(in), dimension(ny) :: fiffy,fifsy,fifwy
  real(mytype), intent(in) :: lind

  call fily_ij(ty,uy,fisy,fiffy,fifsy,fifwy,nx,ny,nz,npaire,lind,2,2)

end subroutine fily_22

subroutine filz_00(tz,uz,rz,fisz,fiffz,fifsz,fifwz,nx,ny,nz,npaire,lind) 

  USE param
  use parfiZ
  use ibm, only : lagpolz, cubsplz

  implicit none

  integer, intent(in) :: nx, ny, nz, npaire
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz, rz
  real(mytype), intent(inout), dimension(nx,ny,nz) :: uz
  real(mytype), intent(out), dimension(nx,ny) :: fisz
  real(mytype), intent(in), dimension(nz) :: fiffz,fifsz,fifwz
  real(mytype), intent(in):: lind

  integer :: i, j, k

  if (iibm == 2) call lagpolz(uz)
  if (iibm == 3) call cubsplz(uz,lind)

  ! Compute r.h.s.
  do concurrent (j=1:ny)
     do concurrent (i=1:nx)
        tz(i,j,1) = fiakz*uz(i,j,1) + fibkz*(uz(i,j,2)+uz(i,j,nz)) &
                  + fickz*(uz(i,j,3)+uz(i,j,nz-1)) &
                  + fidkz*(uz(i,j,4)+uz(i,j,nz-2))
     enddo
  enddo
  do concurrent (j=1:ny)
     do concurrent (i=1:nx)
        tz(i,j,2) = fiakz*uz(i,j,2) + fibkz*(uz(i,j,3)+uz(i,j,1)) &
                  + fickz*(uz(i,j,4)+uz(i,j,nz)) &
                  + fidkz*(uz(i,j,5)+uz(i,j,nz-1))
     enddo
  enddo
  do concurrent (j=1:ny)
     do concurrent (i=1:nx)
        tz(i,j,3) = fiakz*uz(i,j,3) + fibkz*(uz(i,j,4)+uz(i,j,2)) &
                  + fickz*(uz(i,j,5)+uz(i,j,1)) &
                  + fidkz*(uz(i,j,6)+uz(i,j,nz))
     enddo
  enddo
  do concurrent (k=4:nz-3)
     do concurrent (j=1:ny)
        do concurrent (i=1:nx)
           tz(i,j,k) = fiakz*uz(i,j,k) + fibkz*(uz(i,j,k+1)+uz(i,j,k-1)) &
                     + fickz*(uz(i,j,k+2)+uz(i,j,k-2)) &
                     + fidkz*(uz(i,j,k+3)+uz(i,j,k-3))
        enddo
     enddo
  enddo
  do concurrent (j=1:ny)
     do concurrent (i=1:nx)
        tz(i,j,nz-2) = fiakz*uz(i,j,nz-2) + fibkz*(uz(i,j,nz-3)+uz(i,j,nz-1)) &
                     + fickz*(uz(i,j,nz-4)+uz(i,j,nz)) &
                     + fidkz*(uz(i,j,nz-5)+uz(i,j,1))
     enddo
  enddo
  do concurrent (j=1:ny)
     do concurrent (i=1:nx)
        tz(i,j,nz-1) = fiakz*uz(i,j,nz-1) + fibkz*(uz(i,j,nz-2)+uz(i,j,nz)) &
                     + fickz*(uz(i,j,nz-3)+uz(i,j,1)) &
                     + fidkz*(uz(i,j,nz-4)+uz(i,j,2))
     enddo
  enddo
  do concurrent (j=1:ny)
     do concurrent (i=1:nx)
        tz(i,j,nz  ) = fiakz*uz(i,j,nz) + fibkz*(uz(i,j,nz-1)+uz(i,j,1)) &
                     + fickz*(uz(i,j,nz-2)+uz(i,j,2)) &
                     + fidkz*(uz(i,j,nz-3)+uz(i,j,3))
     enddo
  enddo
  do concurrent (j=1:ny)
     do concurrent (i=1:nx)
        rz(i,j,1) = -one
     enddo
  enddo
  do concurrent (k=2:nz-1)
     do concurrent (j=1:ny)
        do concurrent (i=1:nx)
           rz(i,j,k) = zero
        enddo
     enddo
  enddo
  do concurrent (j=1:ny)
     do concurrent (i=1:nx)
        rz(i,j,nz) = fialkz
     enddo
  enddo

  ! Solve tri-diagonal system
  do k=2,nz
     do concurrent (j=1:ny)
        do concurrent (i=1:nx)
           tz(i,j,k) = tz(i,j,k) - tz(i,j,k-1)*fifsz(k)
           rz(i,j,k) = rz(i,j,k) - rz(i,j,k-1)*fifsz(k)
        enddo
     enddo
  enddo
  do concurrent (j=1:ny)
     do concurrent (i=1:nx)
        tz(i,j,nz) = tz(i,j,nz) * fifwz(nz)
        rz(i,j,nz) = rz(i,j,nz) * fifwz(nz)
     enddo
  enddo
  do k=nz-1,1,-1
     do concurrent (j=1:ny)
        do concurrent (i=1:nx)
           tz(i,j,k) = (tz(i,j,k)-fiffz(k)*tz(i,j,k+1)) * fifwz(k)
           rz(i,j,k) = (rz(i,j,k)-fiffz(k)*rz(i,j,k+1)) * fifwz(k)
        enddo
     enddo
  enddo
  do concurrent (j=1:ny)
     do concurrent (i=1:nx)
        fisz(i,j) = (    tz(i,j,1)-fialkz*tz(i,j,nz)) &
                  / (one+rz(i,j,1)-fialkz*rz(i,j,nz))
     enddo
  enddo
  do concurrent (k=1:nz)
     do concurrent (j=1:ny)
        do concurrent (i=1:nx)
           tz(i,j,k) = tz(i,j,k) - fisz(i,j)*rz(i,j,k)
        enddo
     enddo
  enddo

  return
end subroutine filz_00

subroutine filz_ij(tz,uz,fisz,fiffz,fifsz,fifwz,nx,ny,nz,npaire,lind,ncl1,ncln)

  USE param
  use parfiZ
  use ibm, only : lagpolz, cubsplz

  implicit none

  integer, intent(in) :: nx, ny, nz, npaire, ncl1, ncln
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(inout), dimension(nx,ny,nz) :: uz
  real(mytype), intent(out), dimension(nx,ny) :: fisz
  real(mytype), intent(in), dimension(nz) :: fiffz,fifsz,fifwz
  real(mytype), intent(in) :: lind

  integer :: i, j, k

  if (iibm == 2) call lagpolz(uz)
  if (iibm == 3) call cubsplz(uz,lind)

  ! Compute r.h.s.
  if (ncl1 == 1) then
     if (npaire==1) then
        do concurrent (j=1:ny)
           do concurrent (i=1:nx)
              tz(i,j,1) = fiakz*uz(i,j,1) + fibkz*(uz(i,j,2)+uz(i,j,2)) &
                        + fickz*(uz(i,j,3)+uz(i,j,3)) &
                        + fidkz*(uz(i,j,4)+uz(i,j,4))
           enddo
        enddo
        do concurrent (j=1:ny)
           do concurrent (i=1:nx)
              tz(i,j,2) = fiakz*uz(i,j,2) + fibkz*(uz(i,j,3)+uz(i,j,1)) &
                        + fickz*(uz(i,j,4)+uz(i,j,2)) &
                        + fidkz*(uz(i,j,5)+uz(i,j,3))
           enddo
        enddo
        do concurrent (j=1:ny)
           do concurrent (i=1:nx)
              tz(i,j,3) = fiakz*uz(i,j,3) + fibkz*(uz(i,j,4)+uz(i,j,2)) &
                        + fickz*(uz(i,j,5)+uz(i,j,1)) &
                        + fidkz*(uz(i,j,6)+uz(i,j,2))
           enddo
        enddo
     else
        do concurrent (j=1:ny)
           do concurrent (i=1:nx)
              tz(i,j,1) = zero
           enddo
        enddo
        do concurrent (j=1:ny)
           do concurrent (i=1:nx)
              tz(i,j,2) = fiakz*uz(i,j,2) + fibkz*(uz(i,j,3)+uz(i,j,1)) &
                        + fickz*(uz(i,j,4)-uz(i,j,2)) &
                        + fidkz*(uz(i,j,5)-uz(i,j,3))
           enddo
        enddo
        do concurrent (j=1:ny)
           do concurrent (i=1:nx)
              tz(i,j,3) = fiakz*uz(i,j,3) + fibkz*(uz(i,j,4)+uz(i,j,2)) &
                        + fickz*(uz(i,j,5)+uz(i,j,1)) &
                        + fidkz*(uz(i,j,6)-uz(i,j,2))
           enddo
        enddo
     endif
  else
     do concurrent (j=1:ny)
        do concurrent (i=1:nx)
           tz(i,j,1) = uz(i,j,1)
        enddo
     enddo
     do concurrent (j=1:ny)
        do concurrent (i=1:nx)
           tz(i,j,2) = fia2z*uz(i,j,1) + fib2z*uz(i,j,2) + fic2z*uz(i,j,3) &
                     + fid2z*uz(i,j,4)
        enddo
     enddo
     do concurrent (j=1:ny)
        do concurrent (i=1:nx)
           tz(i,j,3) = fia3z*uz(i,j,1) + fib3z*uz(i,j,2) + fic3z*uz(i,j,3) &
                     + fid3z*uz(i,j,4) + fie3z*uz(i,j,5) + fif3z*uz(i,j,6)
        enddo
     enddo
  endif
  do concurrent (k=4:nz-3)
     do concurrent (j=1:ny)
        do concurrent (i=1:nx)
           tz(i,j,k) = fiakz*uz(i,j,k) + fibkz*(uz(i,j,k+1)+uz(i,j,k-1)) &
                     + fickz*(uz(i,j,k+2)+uz(i,j,k-2)) &
                     + fidkz*(uz(i,j,k+3)+uz(i,j,k-3))
        enddo
     enddo
  enddo
  if (ncln == 1) then
     if (npaire == 1) then
        do concurrent (j=1:ny)
           do concurrent (i=1:nx)
              tz(i,j,nz-2) = fiakz*uz(i,j,nz-2) + fibkz*(uz(i,j,nz-1)+uz(i,j,nz-3)) &
                           + fickz*(uz(i,j,nz  )+uz(i,j,nz-4)) &
                           + fidkz*(uz(i,j,nz-1)+uz(i,j,nz-5))
           enddo
        enddo
        do concurrent (j=1:ny)
           do concurrent (i=1:nx)
              tz(i,j,nz-1) = fiakz*uz(i,j,nz-1) + fibkz*(uz(i,j,nz  )+uz(i,j,nz-2)) &
                           + fickz*(uz(i,j,nz-1)+uz(i,j,nz-3)) &
                           + fidkz*(uz(i,j,nz-2)+uz(i,j,nz-4))
           enddo
        enddo
        do concurrent (j=1:ny)
           do concurrent (i=1:nx)
              tz(i,j,nz) = fiakz*uz(i,j,nz) + fibkz*(uz(i,j,nz-1)+uz(i,j,nz-1)) &
                         + fickz*(uz(i,j,nz-2)+uz(i,j,nz-2)) &
                         + fidkz*(uz(i,j,nz-3)+uz(i,j,nz-3))
           enddo
        enddo
     else
        do concurrent (j=1:ny)
           do concurrent (i=1:nx)
              tz(i,j,nz-2) = fiakz*uz(i,j,nz-2) + fibkz*( uz(i,j,nz-1)+uz(i,j,nz-3)) &
                           + fickz*( uz(i,j,nz  )+uz(i,j,nz-4)) &
                           + fidkz*(-uz(i,j,nz-1)+uz(i,j,nz-5))
           enddo
        enddo
        do concurrent (j=1:ny)
           do concurrent (i=1:nx)
              tz(i,j,nz-1) = fiakz*uz(i,j,nz-1) + fibkz*( uz(i,j,nz  )+uz(i,j,nz-2)) &
                           + fickz*(-uz(i,j,nz-1)+uz(i,j,nz-3)) &
                           + fidkz*(-uz(i,j,nz-2)+uz(i,j,nz-4))
           enddo
        enddo
        do concurrent (j=1:ny)
           do concurrent (i=1:nx)
              tz(i,j,nz  ) = zero
           enddo
        enddo
     endif
  else
     do concurrent (j=1:ny)
        do concurrent (i=1:nx)
           tz(i,j,nz-2 ) = fiapz*uz(i,j,nz  ) + fibpz*uz(i,j,nz-1) + ficpz*uz(i,j,nz-2) &
                         + fidpz*uz(i,j,nz-3) + fiepz*uz(i,j,nz-4) + fifpz*uz(i,j,nz-5)
        enddo
     enddo
     do concurrent (j=1:ny)
        do concurrent (i=1:nx)
           tz(i,j,nz-1 ) = fiamz*uz(i,j,nz  ) + fibmz*uz(i,j,nz-1) + ficmz*uz(i,j,nz-2) &
                         + fidmz*uz(i,j,nz-3)
        enddo
     enddo
     do concurrent (j=1:ny)
        do concurrent (i=1:nx)
           tz(i,j,nz   ) = uz(i,j,nz  )
        enddo
     enddo
  endif

  ! Solve tri-diagonal system
  do k=2,nz
     do concurrent (j=1:ny)
        do concurrent (i=1:nx)
           tz(i,j,k) = tz(i,j,k) - tz(i,j,k-1)*fifsz(k)
        enddo
     enddo
  enddo
  do concurrent (j=1:ny)
     do concurrent (i=1:nx)
        tz(i,j,nz) = tz(i,j,nz) * fifwz(nz)
     enddo
  enddo
  do k=nz-1,1,-1
     do concurrent (j=1:ny)
        do concurrent (i=1:nx)
           tz(i,j,k) = (tz(i,j,k)-fiffz(k)*tz(i,j,k+1)) * fifwz(k)
        enddo
     enddo
  enddo

end subroutine filz_ij

subroutine filz_11(tz,uz,rz,fisz,fiffz,fifsz,fifwz,nx,ny,nz,npaire,lind) 

  USE param

  implicit none

  integer, intent(in) :: nx, ny, nz, npaire
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz, rz
  real(mytype), intent(inout), dimension(nx,ny,nz) :: uz
  real(mytype), intent(out), dimension(nx,ny) :: fisz
  real(mytype), intent(in), dimension(nz) :: fiffz,fifsz,fifwz
  real(mytype), intent(in) :: lind

  call filz_ij(tz,uz,fisz,fiffz,fifsz,fifwz,nx,ny,nz,npaire,lind,1,1)

end subroutine filz_11

subroutine filz_12(tz,uz,rz,fisz,fiffz,fifsz,fifwz,nx,ny,nz,npaire,lind) 

  USE param

  implicit none

  integer, intent(in) :: nx, ny, nz, npaire
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz, rz
  real(mytype), intent(inout), dimension(nx,ny,nz) :: uz
  real(mytype), intent(out), dimension(nx,ny) :: fisz
  real(mytype), intent(in), dimension(nz) :: fiffz,fifsz,fifwz
  real(mytype), intent(in) :: lind

  call filz_ij(tz,uz,fisz,fiffz,fifsz,fifwz,nx,ny,nz,npaire,lind,1,2)

end subroutine filz_12

subroutine filz_21(tz,uz,rz,fisz,fiffz,fifsz,fifwz,nx,ny,nz,npaire,lind) 

  USE param

  implicit none

  integer, intent(in) :: nx, ny, nz, npaire
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz, rz
  real(mytype), intent(inout), dimension(nx,ny,nz) :: uz
  real(mytype), intent(out), dimension(nx,ny) :: fisz
  real(mytype), intent(in), dimension(nz) :: fiffz,fifsz,fifwz
  real(mytype), intent(in) :: lind

  call filz_ij(tz,uz,fisz,fiffz,fifsz,fifwz,nx,ny,nz,npaire,lind,2,1)

end subroutine filz_21


subroutine filz_22(tz,uz,rz,fisz,fiffz,fifsz,fifwz,nx,ny,nz,npaire,lind) 

  USE param

  implicit none

  integer, intent(in) :: nx, ny, nz, npaire
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz, rz
  real(mytype), intent(inout), dimension(nx,ny,nz) :: uz
  real(mytype), intent(out), dimension(nx,ny) :: fisz
  real(mytype), intent(in), dimension(nz) :: fiffz,fifsz,fifwz
  real(mytype), intent(in) :: lind

  call filz_ij(tz,uz,fisz,fiffz,fifsz,fifwz,nx,ny,nz,npaire,lind,2,2)

end subroutine filz_22
