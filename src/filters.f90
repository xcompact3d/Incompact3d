subroutine filter(af)

USE decomp_2d
USE param, only: nclx, ncly, ncly 
USE parfiX 
USE parfiY 
USE parfiZ 
USE variables
!=================================================
! Discrete low-pass filter according to 
!=================================================	
implicit none
real(mytype),intent(in) :: af 
real(mytype) :: af1,afn
integer ::  i,j,k
! Set the coefficient for the discrete filter following 
! the tridiagonal filtering of Motheau and Abraham, JCP 2016 
! Filter should be -0.5<filax<0.5

! General Case (entire points)
! alpha*fhat(i-1)+fhat(i)+alpha*fhat(i+1)=af(i)+b/2*[f(i+1)+f(i-1)] + ...

! Coefficients are calculated according to the report of Gaitonde & Visbal, 1998,
! "High-order schemes for Navier-Stokes equations: Algorithm and implementation into FDL3DI"

!========================================
! Define filter coefficients for X-pencil
!========================================
fialx=af                         ! alpha_f
!Interior points
fiaix=(11. + 10.*af)/16.         ! a
fibix=0.5*(15. + 34.*af)/32.     ! b/2 
ficix=0.5*(-3. + 6.*af)/16.      ! c/2
fidix=0.5*(1. - 2.*af)/32.       ! d/2
! Explicit third-order filters near the boundaries!
!Boundary point 1
af1=0.
fia1x=7./8.+af/8.               ! a1/2
fib1x=3./8.+5.*af/8.            ! b1/2
fic1x=-3./8.+3./8.*af           ! c1/2
fid1x=1./8.-1./8.*af            ! d1/2
!Boundary point 2
fia2x=1./8.+3./4.*af    ! a2
fib2x=5./8.+3./4.*af    ! b2/2
fic2x=3./8.+af/4.    ! c2/2
fid2x=-1./8.+af/4.   ! d2/2
!Boundary point n
afn=0.
fianx=7./8.+af/8.               ! a1/2
fibnx=3./8.+5.*af/8.            ! b1/2
ficnx=-3./8.+3./8.*af           ! c1/2
fidnx=1./8.-1./8.*af            ! d1/2
!Boundary point m=n-1
fiamx=1./8.+3./4.*af    ! a2
fibmx=5./8.+3./4.*af    ! b2/2
ficmx=3./8.+af/4.    ! c2/2
fidmx=-1./8.+af/4.   ! d2/2
! Set the coefficients for the matrix A
! Periodic case
if (nclx.eq.0) then
   fiffx(1)   =af
   fiffx(2)   =af
   fiffx(nx-2)=af
   fiffx(nx-1)=af
   fiffx(nx)  =0.
   fifcx(1)   =2.
   fifcx(2)   =1.
   fifcx(nx-2)=1.
   fifcx(nx-1)=1.
   fifcx(nx  )=1.+af*af
   fifbx(1)   =af
   fifbx(2)   =af
   fifbx(nx-2)=af
   fifbx(nx-1)=af
   fifbx(nx  )=0.
   do i=3,nx-3
      fiffx(i)=af
      fifcx(i)=1.
      fifbx(i)=af
   enddo   
endif
if (nclx.eq.1) then
   fiffx(1)   =af1
   fiffx(2)   =af
   fiffx(ny-2)=af
   fiffx(ny-1)=af
   fiffy(ny)  =0.
   fifcx(1)   =1.
   fifcx(2)   =1.
   fifcx(nx-2)=1.
   fifcx(nx-1)=1.
   fifcx(nx  )=1.
   fifbx(1)   =af 
   fifbx(2)   =af
   fifbx(nx-2)=af
   fifbx(nx-1)=afn
   fifbx(nx  )=0.
   do i=3,nx-3
      fiffy(i)=af
      fifcy(i)=1.
      fifby(i)=af
   enddo
endif
if (nclx.eq.2) then
   fiffx(1)   =af1 
   fiffx(2)   =af
   fiffx(nx-2)=af
   fiffx(nx-1)=af
   fiffx(nx)  =0.
   fifcx(1)   =1.
   fifcx(2)   =1.
   fifcx(nx-2)=1.
   fifcx(nx-1)=1.
   fifcx(nx  )=1.
   fifbx(1)   =af 
   fifbx(2)   =af
   fifbx(nx-2)=af
   fifbx(nx-1)=afn
   fifbx(nx  )=0.
   do i=3,nx-3
      fiffx(i)=af
      fifcx(i)=1.
      fifbx(i)=af
   enddo
endif
! Prepare coefficients to be used in the Thomas Algorithm
do i=1,nx
   fiffxp(i)=fiffx(i)
enddo
call prepare (fifbx,fifcx,fiffx,fifsx,fifwx,nx)
call prepare (fifbx,fifcx,fiffxp,fifsxp,fifwxp,nx)
if (nclx.eq.1) then
   fiffx(1)=0.0_mytype
   fifbx(nx-1)=0.0_mytype
call prepare (fifbx,fifcx,fiffx,fifsx,fifwx,nx)
endif

!========================================
! Define filter coefficients for Y-pencil
!========================================
fialy=af                         ! alpha_f
!Interior points
fiajy=(11. + 10.*af)/16.         ! a
fibjy=0.5*(15. + 34.*af)/32.     ! b/2 
ficjy=0.5*(-3. + 6.*af)/16.      ! c/2
fidjy=0.5*(1. - 2.*af)/32.       ! d/2
!Boundary point 1
fia1y=7./8.+af/8.               ! a1/2
fib1y=3./8.+5.*af/8.            ! b1/2
fic1y=-3./8.+3./8.*af           ! c1/2
fid1y=1./8.-1./8.*af            ! d1/2
!Boundary point 2
fia2y=1./8.+3./4.*af    ! a2
fib2y=5./8.+3./4.*af    ! b2/2
fic2y=3./8.+af/4.    ! c2/2
fid2y=-1./8.+af/4.   ! d2/2
!Boundary point n
fiany=7./8.+af/8.               ! a1/2
fibny=3./8.+5.*af/8.            ! b1/2
ficny=-3./8.+3./8.*af           ! c1/2
fidny=1./8.-1./8.*af            ! d1/2
!Boundary point m=n-1
fiamy=1./8.+3./4.*af    ! a2
fibmy=5./8.+3./4.*af    ! b2/2
ficmy=3./8.+af/4.    ! c2/2
fidmy=-1./8.+af/4.   ! d2/2
! Define coefficients
if (ncly.eq.0) then
   fiffy(1)   =af
   fiffy(2)   =af
   fiffy(ny-2)=af
   fiffy(ny-1)=af
   fiffy(ny)  =0.
   fifcy(1)   =2.
   fifcy(2)   =1.
   fifcy(ny-2)=1.
   fifcy(ny-1)=1.
   fifcy(ny  )=1.+af*af
   fifby(1)   =af
   fifby(2)   =af
   fifby(ny-2)=af
   fifby(ny-1)=af
   fifby(ny  )=0.
   do j=3,ny-3
      fiffy(j)=af
      fifcy(j)=1.
      fifby(j)=af
   enddo   
endif
if (ncly.eq.1) then
   fiffy(1)   =af+af
   fiffy(2)   =af
   fiffy(ny-2)=af
   fiffy(ny-1)=af
   fiffy(ny)  =0.0_mytype  
   fifcy(1)   =1.0_mytype
   fifcy(2)   =1.0_mytype
   fifcy(ny-2)=1.0_mytype
   fifcy(ny-1)=1.0_mytype
   fifcy(ny  )=1.0_mytype
   fifby(1)   =af 
   fifby(2)   =af
   fifby(ny-2)=af
   fifby(ny-1)=af+af
   fifby(ny  )=0.0_mytype
   do j=3,ny-3
      fiffy(j)=af
      fifcy(j)=1.0_mytype
      fifby(j)=af
   enddo
endif
if (ncly.eq.2) then
   fiffy(1)   =af1 
   fiffy(2)   =af
   fiffy(ny-2)=af
   fiffy(ny-1)=af
   fiffy(ny)  =0.0_mytype 
   fifcy(1)   =1.0_mytype
   fifcy(2)   =1.0_mytype
   fifcy(ny-2)=1.0_mytype
   fifcy(ny-1)=1.0_mytype
   fifcy(ny  )=1.0_mytype
   fifby(1)   =af 
   fifby(2)   =af
   fifby(ny-2)=af
   fifby(ny-1)=afn
   fifby(ny  )=0.0_mytype
   do j=3,ny-3
      fiffy(j)=af
      fifcy(j)=1.0_mytype
      fifby(j)=af
   enddo
endif
do j=1,ny
   fiffyp(j)=fiffy(j)
enddo
call prepare (fifby,fifcy,fiffy,fifsy,fifwy,ny)
call prepare (fifby,fifcy,fiffyp,fifsyp,fifwyp,ny)
if (ncly.eq.1) then
   fiffy(1)=0.0_mytype
   fifby(ny-1)=0.0_mytype
call prepare (fifby,fifcy,fiffy,fifsy,fifwy,ny)
endif
!========================================
! Define filter coefficients for Z-pencil
!========================================
fialz=af                         ! alpha_f
!Interior points
fiakz=(11. + 10.*af)/16.         ! a
fibkz=0.5*(15. + 34.*af)/32.     ! b/2 
fickz=0.5*(-3. + 6.*af)/16.      ! c/2
fidkz=0.5*(1. - 2.*af)/32.       ! d/2
!Boundary point 1
fia1z=7./8.+af/8.               ! a1/2
fib1z=3./8.+5.*af/8.            ! b1/2
fic1z=-3./8.+3./8.*af           ! c1/2
fid1z=1./8.-1./8.*af            ! d1/2
!Boundary point 2
fia2z=1./8.+3./4.*af    ! a2
fib2z=5./8.+3./4.*af    ! b2/2
fic2z=3./8.+af/4.    ! c2/2
fid2z=-1./8.+af/4.   ! d2/2
!Boundary point n
fianz=7./8.+af/8.               ! a1/2
fibnz=3./8.+5.*af/8.            ! b1/2
ficnz=-3./8.+3./8.*af           ! c1/2
fidnz=1./8.-1./8.*af            ! d1/2
!Boundary point m=n-1
fiamz=1./8.+3./4.*af    ! a2
fibmz=5./8.+3./4.*af    ! b2/2
ficmz=3./8.+af/4.    ! c2/2
fidmz=-1./8.+af/4.   ! d2/2
if (nclz.eq.0) then
      fiffz(1)   =af
      fiffz(2)   =af
      fiffz(nz-2)=af
      fiffz(nz-1)=af
      fiffz(nz)  =0.
      fifcz(1)   =2.
      fifcz(2)   =1.
      fifcz(nz-2)=1.
      fifcz(nz-1)=1.
      fifcz(nz  )=1.+af*af
      fifbz(1)   =af
      fifbz(2)   =af
      fifbz(nz-2)=af
      fifbz(nz-1)=af
      fifbz(nz  )=0.
      do k=3,nz-3
         fiffz(k)=af
         fifcz(k)=1.
         fifbz(k)=af
      enddo
endif
if (nclz.eq.1) then
   fiffz(1)   =af+af
   fiffz(2)   =af
   fiffz(nz-2)=af
   fiffz(nz-1)=af
   fiffz(nz)  =0.
   fifcz(1)   =1.
   fifcz(2)   =1.
   fifcz(nz-2)=1.
   fifcz(nz-1)=1.
   fifcz(nz  )=1.
   fifbz(1)   =af 
   fifbz(2)   =af
   fifbz(nz-2)=af
   fifbz(nz-1)=af+af
   fifbz(nz  )=0.
   do k=3,nz-3
      fiffz(k)=af
      fifcz(k)=1.
      fifbz(k)=af
   enddo
endif
if (nclz.eq.2) then
   fiffz(1)   =af1
   fiffz(2)   =af
   fiffz(nz-2)=af
   fiffz(nz-1)=af
   fiffz(nz)  =0.
   fifcz(1)   =1.
   fifcz(2)   =1.
   fifcz(nz-2)=1.
   fifcz(nz-1)=1.
   fifcz(nz  )=1.
   fifbz(1)   =af 
   fifbz(2)   =af
   fifbz(nz-2)=af
   fifbz(nz-1)=afn
   fifbz(nz  )=0.
   do k=3,nz-3
      fiffz(k)=af
      fifcz(k)=1.
      fifbz(k)=af
   enddo
endif
do k=1,nz
   fiffzp(k)=fiffz(k)
enddo
call prepare (fifbz,fifcz,fiffz,fifsz,fifwz,nz)
call prepare (fifbz,fifcz,fiffzp,fifszp,fifwzp,nz)
if (nclz.eq.1) then
   fiffz(1)=0.0_mytype
   fifbz(nz-1)=0.0_mytype
call prepare (fifbz,fifcz,fiffz,fifsz,fifwz,nz)
endif

return 

end subroutine filter


subroutine filx(tx,ux,rx,fisx,fiffx,fifsx,fifwx,nx,ny,nz,npaire) 
  
USE param  
USE parfiX 

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: tx,ux,rx 
real(mytype), dimension(ny,nz) :: fisx
real(mytype), dimension(nx) :: fiffx,fifsx,fifwx

! This is for periodic boundary conditions
if (nclx==0) then 
   do k=1,nz 
   do j=1,ny 
      tx(1,j,k)=fiaix*ux(1,j,k)+fibix*(ux(2,j,k)+ux(nx,j,k))& 
                               +ficix*(ux(3,j,k)+ux(nx-1,j,k))&
                               +fidix*(ux(4,j,k)+ux(nx-2,j,k)) 
      rx(1,j,k)=-1.
      tx(2,j,k)=fiaix*ux(2,j,k)+fibix*(ux(3,j,k)+ux(1,j,k))&
                               +ficix*(ux(4,j,k)+ux(nx,j,k))& 
                               +fidix*(ux(5,j,k)+ux(nx-1,j,k)) 
      rx(2,j,k)=0. 
      tx(3,j,k)=fiaix*ux(3,j,k)+fibix*(ux(4,j,k)+ux(2,j,k))&
          
                               +ficix*(ux(5,j,k)+ux(1,j,k))& 
                               +fidix*(ux(6,j,k)+ux(nx,j,k)) 
      rx(3,j,k)=0. 
      do i=4,nx-3
         tx(i,j,k)=fiaix*ux(i,j,k)+fibix*(ux(i+1,j,k)+ux(i-1,j,k))& 
                                  +ficix*(ux(i+2,j,k)+ux(i-2,j,k))&
                                  +fidix*(ux(i+3,j,k)+ux(i-3,j,k)) 
         rx(i,j,k)=0. 
      enddo
      tx(nx-2,j,k)=fiaix*ux(nx-2,j,k)+fibix*(ux(nx-3,j,k)+ux(nx-1,j,k))&
                                     +ficix*(ux(nx-4,j,k)+ux(nx,j,k))& 
                                     +fidix*(ux(nx-5,j,k)+ux(1,j,k)) 
      rx(nx-2,j,k)=0. 
      tx(nx-1,j,k)=fiaix*ux(nx-1,j,k)+fibix*(ux(nx-2,j,k)+ux(nx,j,k))&
                                     +ficix*(ux(nx-3,j,k)+ux(1,j,k))& 
                                     +fidix*(ux(nx-4,j,k)+ux(2,j,k)) 
      rx(nx-1,j,k)=0. 
      tx(nx,j,k)=fiaix*ux(nx,j,k)+fibix*(ux(nx-1,j,k)+ux(1,j,k))&
                                 +ficix*(ux(nx-2,j,k)+ux(2,j,k))& 
                                 +fidix*(ux(nx-3,j,k)+ux(3,j,k)) 
      rx(nx,j,k)=fialx           
      do i=2, nx
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fifsx(i) 
         rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*fifsx(i) 
      enddo
      tx(nx,j,k)=tx(nx,j,k)*fifwx(nx) 
      rx(nx,j,k)=rx(nx,j,k)*fifwx(nx) 
      do i=nx-1,1,-1
         tx(i,j,k)=(tx(i,j,k)-fiffx(i)*tx(i+1,j,k))*fifwx(i) 
         rx(i,j,k)=(rx(i,j,k)-fiffx(i)*rx(i+1,j,k))*fifwx(i) 
      enddo
        fisx(j,k)=(tx(1,j,k)-fialx*tx(nx,j,k))&
           /(1.+rx(1,j,k)-fialx*rx(nx,j,k)) 
      do i=1,nx 
         tx(i,j,k)=tx(i,j,k)-fisx(j,k)*rx(i,j,k) 
      enddo
   enddo
   enddo
endif

if (nclx==2) then 
   do k=1,nz
   do j=1,ny 
      tx(1,j,k)=ux(1,j,k)
      tx(2,j,k)=fia2x*ux(1,j,k)+fib2x*ux(2,j,k)+fic2x*ux(3,j,k)+&
                fid2x*ux(4,j,k)
      do i=3,nx-2
         tx(i,j,k)=fiaix*ux(i,j,k)+fibix*(ux(i+1,j,k)+ux(i-1,j,k))& 
                                  +ficix*(ux(i+2,j,k)+ux(i-2,j,k))&
                                  +fidix*(ux(i+3,j,k)+ux(i-3,j,k)) 
      enddo
      tx(nx,j,k)=ux(nx,j,k)
      tx(nx-1,j,k)=fiamx*ux(nx,j,k)+fibmx*ux(nx-1,j,k)+ficmx*ux(nx-2,j,k)+&
                            fidmx*ux(nx-3,j,k)
      do i=2,nx 
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fifsx(i) 
      enddo
      tx(nx,j,k)=tx(nx,j,k)*fifwx(nx) 
      do i=nx-1,1,-1
         tx(i,j,k)=(tx(i,j,k)-fiffx(i)*tx(i+1,j,k))*fifwx(i) 
      enddo
   enddo 
   enddo 
endif


return  
end subroutine filx


!********************************************************************
!
subroutine fily(ty,uy,ry,fisy,fiffy,fifsy,fifwy,ppy,nx,ny,nz,npaire) 
!
!********************************************************************
  
USE param  
USE parfiY 

implicit none

integer :: nx,ny,nz,i,j,k,npaire
real(mytype), dimension(nx,ny,nz) :: ty,uy 
real(mytype), dimension(nx,ny,nz) :: ry
real(mytype), dimension(nx,nz)  :: fisy
real(mytype), dimension(ny) :: fiffy,fifsy,fifwy,ppy

if (ncly==0) then 
   do k=1,nz 
   do i=1,nx 
      ty(i,1,k)=fiajy*uy(i,1,k)+fibjy*(uy(i,2,k)+uy(i,nx,k))& 
                               +ficjy*(uy(i,3,k)+uy(i,nx-1,k))&
                               +fidjy*(uy(i,4,k)+uy(i,nx-2,k)) 
      ry(i,1,k)=-1.
      ty(i,2,k)=fiajy*uy(i,2,k)+fibjy*(uy(i,3,k)+uy(i,1,k))&
                               +ficjy*(uy(i,4,k)+uy(i,nx,k))& 
                               +fidjy*(uy(i,5,k)+uy(i,nx-1,k)) 
      ry(i,2,k)=0. 
      ty(i,3,k)=fiajy*uy(i,3,k)+fibjy*(uy(i,4,k)+uy(i,2,k))&
                               +ficjy*(uy(i,5,k)+uy(i,1,k))& 
                               +fidjy*(uy(i,6,k)+uy(i,nx,k)) 
      ry(i,3,k)=0. 
      do j=4,ny-3
         ty(i,j,k)=fiajy*uy(i,j,k)+fibjy*(uy(i,j+1,k)+uy(i,j-1,k))& 
                                  +ficjy*(uy(i,j+2,k)+uy(i,j-2,k))&
                                  +fidjy*(uy(i,j+3,k)+uy(i,j-3,k)) 
         ry(i,j,k)=0. 
      enddo
      ty(i,ny-2,k)=fiajy*uy(i,ny-2,k)+fibjy*(uy(i,ny-3,k)+uy(i,ny-1,k))&
                                     +ficjy*(uy(i,ny-4,k)+uy(i,ny,k))& 
                                     +fidjy*(uy(i,ny-5,k)+uy(i,1,k)) 
      ry(i,ny-2,k)=0. 
      ty(i,ny-1,k)=fiajy*uy(i,ny-1,k)+fibjy*(uy(i,ny-2,k)+uy(i,ny,k))&
                                     +ficjy*(uy(i,ny-3,k)+uy(i,1,k))& 
                                     +fidjy*(uy(i,ny-4,k)+uy(i,2,k)) 
      ry(i,ny-1,k)=0. 
      ty(i,ny,k)=fiajy*uy(i,ny,k)+fibjy*(uy(i,ny-1,k)+uy(i,1,k))&
                                 +ficjy*(uy(i,ny-2,k)+uy(i,2,k))& 
                                 +fidjy*(uy(i,ny-3,k)+uy(i,3,k)) 
      ry(i,ny,k)=fialy           
      do j=2, ny
         ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*fifsy(j) 
         ry(i,j,k)=ry(i,j,k)-ry(i,j-1,k)*fifsy(j) 
      enddo
      ty(i,ny,k)=ty(i,ny,k)*fifwy(ny) 
      ry(i,ny,k)=ry(i,ny,k)*fifwy(ny) 
      do j=ny-1,1,-1
         ty(i,j,k)=(ty(i,j,k)-fiffy(j)*ty(i,j+1,k))*fifwy(j) 
         ry(i,j,k)=(ry(i,j,k)-fiffy(j)*ry(i,j+1,k))*fifwy(j) 
      enddo
        fisy(i,k)=(ty(i,1,k)-fialy*ty(i,ny,k))&
           /(1.+ry(i,1,k)-fialy*ry(i,ny,k)) 
      do j=1,ny 
         ty(i,j,k)=ty(i,j,k)-fisy(i,k)*ry(i,j,k) 
      enddo
   enddo
   enddo
endif


if (ncly==1) then 
    if (npaire==1) then 
    do k=1,nz 
    do i=1,nx 
         ty(i,1,k)=fiajy*uy(i,1,k)+fibjy*(uy(i,2,k)+uy(i,2,k))&
                                  +ficjy*(uy(i,3,k)+uy(i,3,k))&
                                  +fidjy*(uy(i,4,k)+uy(i,4,k))
         ty(i,2,k)=fiajy*uy(i,2,k)+fibjy*(uy(i,3,k)+uy(i,1,k))& 
                                  +ficjy*(uy(i,4,k)+uy(i,2,k))&
                                  +fidjy*(uy(i,5,k)+uy(i,3,k)) 
         ty(i,3,k)=fiajy*uy(i,3,k)+fibjy*(uy(i,4,k)+uy(i,2,k))& 
                                  +ficjy*(uy(i,5,k)+uy(i,1,k))&
                                  +fidjy*(uy(i,6,k)+uy(i,2,k)) 
    enddo
    enddo 
    do k=1,nz 
    do j=4,ny-3 
    do i=1,nx 
       ty(i,j,k)=fiajy*uy(i,j,k)+fibjy*(uy(i,j+1,k)+uy(i,j-1,k))& 
                                +ficjy*(uy(i,j+2,k)+uy(i,j-2,k))&
                                +fidjy*(uy(i,j+3,k)+uy(i,j-3,k)) 
    enddo
    enddo 
    enddo 
    do k=1,nz 
    do i=1,nx 
       ty(i,ny,k)=fiajy*uy(i,ny,k)+fibjy*(uy(i,ny-1,k)+uy(i,ny-1,k))&
                                  +ficjy*(uy(i,ny-2,k)+uy(i,ny-2,k))&
                                  +fidjy*(uy(i,ny-3,k)+uy(i,ny-3,k))
       ty(i,ny-1,k)=fiajy*uy(i,ny-1,k)+fibjy*(uy(i,ny,k)  +uy(i,ny-2,k))& 
                                      +ficjy*(uy(i,ny-1,k)+uy(i,ny-3,k))&
                                      +fidjy*(uy(i,ny-2,k)+uy(i,ny-4,k)) 
       ty(i,ny-2,k)=fiajy*uy(i,ny-2,k)+fibjy*(uy(i,ny-1,k)+uy(i,ny-3,k))& 
                                      +ficjy*(uy(i,ny,k)+uy(i,ny-4,k))&
                                      +fidjy*(uy(i,ny-1,k)+uy(i,ny-5,k)) 
    enddo
    enddo 
    do k=1,nz
    do j=2,ny  
    do i=1,nx 
       ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*fifsy(j) 
    enddo
    enddo
    enddo 
    do k=1,nz 
    do i=1,nx 
       ty(i,ny,k)=ty(i,ny,k)*fifwy(ny) 
    enddo 
    enddo 
    do k=1,nz
    do j=ny-1,1,-1  
    do i=1,nx 
    ty(i,j,k)=(ty(i,j,k)-fiffy(j)*ty(i,j+1,k))*fifwy(j) 
    enddo 
    enddo 
    enddo 
   endif
   if (npaire==0) then 
      do k=1,nz 
      do i=1,nx 
         ty(i,1,k)=fiajy*uy(i,1,k)
         ty(i,2,k)=fiajy*uy(i,2,k)+fibjy*(uy(i,3,k)+uy(i,1,k))& 
                                  +ficjy*(uy(i,4,k)-uy(i,2,k))&
                                  +fidjy*(uy(i,5,k)-uy(i,3,k)) 
         ty(i,3,k)=fiajy*uy(i,3,k)+fibjy*(uy(i,4,k)+uy(i,2,k))& 
                                  +ficjy*(uy(i,5,k)+uy(i,1,k))&
                                  +fidjy*(uy(i,6,k)-uy(i,2,k)) 
      enddo
      enddo 
      do k=1,nz 
      do j=4,ny-3 
      do i=1,nx 
         ty(i,j,k)=fiajy*uy(i,j,k)+fibjy*(uy(i,j+1,k)+uy(i,j-1,k))& 
                                  +ficjy*(uy(i,j+2,k)+uy(i,j-2,k))&
                                  +fidjy*(uy(i,j+3,k)+uy(i,j-3,k)) 
      enddo
      enddo 
      enddo 
      do k=1,nz 
      do i=1,nx 
         ty(i,ny,k)=fiajy*uy(i,ny,k)
         ty(i,ny-1,k)=fiajy*uy(i,ny-1,k)+fibjy*(uy(i,ny,k)+uy(i,ny-2,k))& 
                                        +ficjy*(-uy(i,ny-1,k)+uy(i,ny-3,k))&
                                        +fidjy*(-uy(i,ny-2,k)+uy(i,ny-4,k)) 
         ty(i,ny-2,k)=fiajy*uy(i,ny-2,k)+fibjy*(uy(i,ny-1,k)+uy(i,ny-3,k))& 
                                        +ficjy*(uy(i,ny,k)+uy(i,ny-4,k))&
                                        +fidjy*(-uy(i,ny-1,k)+uy(i,ny-5,k)) 
      enddo
      enddo 
      do k=1,nz
      do j=2,ny  
      do i=1,nx 
         ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*fifsy(j) 
      enddo
      enddo
      enddo 
      do k=1,nz 
      do i=1,nx 
         ty(i,ny,k)=ty(i,ny,k)*fifwy(ny) 
      enddo 
      enddo 
      do k=1,nz
      do j=ny-1,1,-1  
      do i=1,nx 
      ty(i,j,k)=(ty(i,j,k)-fiffy(j)*ty(i,j+1,k))*fifwy(j) 
      enddo 
      enddo 
      enddo 
   endif
endif

if (ncly==2) then 
   do k=1,nz
   do i=1,nx 
      ty(i,1,k)=uy(i,1,k)
      ty(i,2,k)=fia2y*uy(i,1,k)+fib2y*uy(i,2,k)+fic2y*uy(i,3,k)+&
                fid2y*uy(i,4,k)
      do j=3,ny-2
         ty(i,j,k)=fiajy*uy(i,j,k)+fibjy*(uy(i,j+1,k)+uy(i,j-1,k))& 
                                  +ficjy*(uy(i,j+2,k)+uy(i,j-2,k))&
                                  +fidjy*(uy(i,j+3,k)+uy(i,j-3,k)) 
      enddo
      ty(i,ny,k)=uy(i,ny,k)
      ty(i,ny-1,k)=fiamy*uy(i,ny,k)+fibmy*uy(i,ny-1,k)+ficmy*uy(i,ny-2,k)+&
                                    fidmy*uy(i,ny-3,k)
      do j=2,ny 
         ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*fifsy(j) 
      enddo
      ty(i,ny,k)=ty(i,ny,k)*fifwy(ny) 
      do j=ny-1,1,-1
         ty(i,j,k)=(ty(i,j,k)-fiffy(j)*ty(i,j+1,k))*fifwy(j) 
      enddo
   enddo 
   enddo 
endif

end subroutine fily


subroutine filz(tz,uz,rz,fisz,fiffz,fifsz,fifwz,nx,ny,nz,npaire) 
  
USE param  
USE parfiZ 

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: tz,uz,rz
real(mytype), dimension(nx,ny) :: fisz
real(mytype), dimension(nz) :: fiffz,fifsz,fifwz

! This is for periodic boundary conditions
if (nclz==0) then 
   do j=1,ny 
   do i=1,nx 
      tz(i,j,1)=fiakz*uz(i,j,1)+fibkz*(uz(i,j,2)+uz(i,j,nz))& 
                               +fickz*(uz(i,j,3)+uz(i,j,nz-1))&
                               +fidkz*(uz(i,j,4)+uz(i,j,nz-2)) 
      rz(i,j,1)=-1.
      tz(i,j,2)=fiakz*uz(i,j,2)+fibkz*(uz(i,j,3)+uz(i,j,1))&
                               +fickz*(uz(i,j,4)+uz(i,j,nz))& 
                               +fidkz*(uz(i,j,5)+uz(i,j,nz-1)) 
      rz(i,j,2)=0. 
      tz(i,j,3)=fiakz*uz(i,j,3)+fibkz*(uz(i,j,4)+uz(i,j,2))&
                               +fickz*(uz(i,j,5)+uz(i,j,1))& 
                               +fidkz*(uz(i,j,6)+uz(i,j,nz)) 
      rz(i,j,3)=0.
   enddo
   enddo 
   do k=4,nz-3
   do j=1,ny
   do i=1,nx
         tz(i,j,k)=fiakz*uz(i,j,k)+fibkz*(uz(i,j,k+1)+uz(i,j,k-1))& 
                                  +fickz*(uz(i,j,k+2)+uz(i,j,k-2))&
                                  +fidkz*(uz(i,j,k+3)+uz(i,j,k-3)) 
         rz(i,j,k)=0. 
   enddo
   enddo
   enddo 
   do j=1,ny 
   do i=1,nx 
      tz(i,j,nz-2)=fiakz*uz(i,j,nz-2)+fibkz*(uz(i,j,nz-3)+uz(i,j,nz-1))&
                                     +fickz*(uz(i,j,nz-4)+uz(i,j,nz))& 
                                     +fidkz*(uz(i,j,nz-5)+uz(i,j,1)) 
      rz(i,j,nz-2)=0. 
      tz(i,j,nz-1)=fiakz*uz(i,j,nz-1)+fibkz*(uz(i,j,nz-2)+uz(i,j,nz))&
                                     +fickz*(uz(i,j,nz-3)+uz(i,j,1))& 
                                     +fidkz*(uz(i,j,nz-4)+uz(i,j,2)) 
      rz(i,j,nz-1)=0. 
      tz(i,j,nz)=fiakz*uz(i,j,nz)+fibkz*(uz(i,j,nz-1)+uz(i,j,1))&
                                 +fickz*(uz(i,j,nz-2)+uz(i,j,2))& 
                                 +fidkz*(uz(i,j,nz-3)+uz(i,j,3)) 
      rz(i,j,nz)=fialz           
   enddo
   enddo
   do k=2,nz
   do j=1,ny
   do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*fifsz(k) 
         rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*fifsz(k) 
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx 
      tz(i,j,nz)=tz(i,j,nz)*fifwz(nz) 
      rz(i,j,nz)=rz(i,j,nz)*fifwz(nz) 
   enddo
   enddo
   do k=nz-1,1,-1
   do j=1,ny
   do i=1,nx    
         tz(i,j,k)=(tz(i,j,k)-fiffz(k)*tz(i,j,k+1))*fifwz(k) 
         rz(i,j,k)=(rz(i,j,k)-fiffz(k)*rz(i,j,k+1))*fifwz(k) 
   enddo
   enddo
   enddo   
   do j=1,ny
   do i=1,nx      
    fisz(i,j)=(tz(i,j,1)-fialz*tz(i,j,nz))&
           /(1.+rz(i,j,1)-fialz*rz(i,j,nz)) 
   enddo
   enddo    
   do k=1,nz 
   do j=1,ny
   do i=1,nx      
         tz(i,j,k)=tz(i,j,k)-fisz(i,j)*rz(i,j,k) 
   enddo
   enddo
   enddo
endif

return  
end subroutine filz

