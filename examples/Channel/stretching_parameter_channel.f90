program xcompact3d_stretching_channel

  implicit none
  
  integer, parameter :: n=128
  real(8), dimension(n) ::  yeta, phi
  real(8) :: re, ret, alpha, beta, betamin, eps, pi, betaold
  real(8) :: den, den1, den2, den3, den4, xnum1, y, f, cst
  real(8) :: dy, dymax, dymin, xcx, xl2, xnum, yly, yold, yp, yinf
  integer :: i, npvis
  
  re=4200.  !#TO CHANGE#!input Reynolds number in Xcompact3d
                        !based on centreline velocity Poisseuille profile
  ret=0.123*re**0.875   !corresponding Reynolds number based on u_tau
    
  yp=1.      !#TO CHANGE#!targeted resolution at the wall in wall units
  y=1000

  beta=1.        !initiation stretching parameter
  betamin=0.
  eps=1.e-5

  !INITIALISATION of stretching function as defined in Xcompact3d
  do while (abs(y-yp).gt.eps)
     if (y-yp.gt.0.) then
        betaold=beta
        beta=beta-0.5*beta
     else
        beta=beta+0.5*abs(beta-betaold)
     endif
     
     pi=acos(-1.)
   
     yly=2.        !domain size
     yinf=-yly/2.
     den=2.*beta*yinf
     xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf)
     alpha=abs(xnum/den)
     xcx=1./beta/alpha
     phi(1)=yinf*ret
     do i=2,n
        yeta(i)=(i-1.)*(1./(n-1.))-0.5
        den1=sqrt(alpha*beta+1.)
        xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)     
        den=2.*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
        den3=((sin(pi*yeta(i)))*(sin(pi*yeta(i)))/beta/pi)+alpha/pi
        den4=2.*alpha*beta-cos(2.*pi*yeta(i))+1.
        xnum1=(atan(xnum*tan(pi*yeta(i))))*den4/den1/den3/den
        cst=sqrt(beta)*pi/(2.*sqrt(alpha)*sqrt(alpha*beta+1.))
        if (yeta(i).lt.0.5) phi(i)=xnum1-cst+yly/2.
        if (yeta(i).eq.0.5) phi(i)=0.+yly/2.
        if (yeta(i).gt.0.5) phi(i)=xnum1+cst+yly/2.
        phi(i)=phi(i)*ret
     enddo
     y=(phi(n)-phi(n-1))   
  enddo

  !calculation of the number of mesh nodes in the viscous sublayer
  npvis=1
  xl2=0.
  do i=2,n/2+1
     if (xl2+(phi(i)-phi(i-1)).le.5.) npvis=npvis+1
     xl2=xl2+(phi(i)-phi(i-1))
  enddo

  dymin=(phi(n)-phi(n-1))
  dymax=(phi(n/2)-phi(n/2-1))
  
  print *,'number of mesh node in wall normal direction : ',n
  print *,'Reynolds number based on u_tau: ',ret
  print *
  print *,'beta parameter = ',beta
  print *
  print *,'Min mesh size in wall normal direction : y+ = ',dymin
  print *,'Max mesh size in wall normal direction : y+ = ',dymax
  print *,'Number of mesh nodes in viscous sublayer : ', npvis
  print *

  !comparison with a Tchebichev stretching
  y=-1.
  pi=acos(-1.)
  dymin= 1.e10
  dymax=-1.e10
  npvis=1
  do i=2,n/2+1
     yold=y
     f=(i-1.)/(n-1.)
     y=-cos(f*pi)
     if ((y+1)*ret.le.5.) npvis=npvis+1
     dy=y-yold
     if (dy.lt.dymin) dymin=dy
     if (dy.gt.dymax) dymax=dy
  enddo
  
  dymin=dymin*ret
  dymax=dymax*ret
  
  print *
  print *,'Tchebichev :' 
  print *,'Min mesh size in wall normal direction : y+ = ',dymin
  print *,'Max mesh size in wall normal direction : y+ = ',dymax
  print *,'Number of mesh nodes in viscous sublayer : ', npvis

end program xcompact3d_stretching_channel
