module ibm

  public

contains
  !############################################################################
  subroutine corgp_IBM (ux,uy,uz,px,py,pz,nlock)
    USE param
    USE decomp_2d
    USE variables
    implicit none
    integer :: i,j,k,nlock
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,px,py,pz
    if (nz == 1) then
       write(*,*) "2D currently unsupported - see ibm.f90"
       stop
    endif
    if (nlock == 1) then
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                ux(i,j,k)=-px(i,j,k)+ux(i,j,k)
                uy(i,j,k)=-py(i,j,k)+uy(i,j,k)
                uz(i,j,k)=-pz(i,j,k)+uz(i,j,k)
             enddo
          enddo
       enddo
    endif
    if (nlock == 2) then
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                ux(i,j,k)=px(i,j,k)+ux(i,j,k)
                uy(i,j,k)=py(i,j,k)+uy(i,j,k)
                uz(i,j,k)=pz(i,j,k)+uz(i,j,k)
             enddo
          enddo
       enddo
    endif

    return
  end subroutine corgp_IBM
  !############################################################################
  !############################################################################
  subroutine body(ux1,uy1,uz1,ep1)
    use param, only : zero, one, dx, dz
    use decomp_2d, only : xstart, xend, xsize, mytype, nrank
    !use decomp_2d_io
    use variables, only : ny
    implicit none
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    integer :: i,j,k

#ifdef DEBG
    if (nrank  ==  0) write(*,*) '# body start'
#endif

    do k = 1, xsize(3)
      do j = 1, xsize(2)
         do i = 1, xsize(1)
            ux1(i,j,k)=(one-ep1(i,j,k))*ux1(i,j,k)
            uy1(i,j,k)=(one-ep1(i,j,k))*uy1(i,j,k)
            uz1(i,j,k)=(one-ep1(i,j,k))*uz1(i,j,k)
         enddo
      enddo
    enddo

#ifdef DEBG
    if (nrank  ==  0) write(*,*) '# body done'
#endif

    return
  end subroutine body
  !############################################################################
  !############################################################################
  subroutine lagpolx(u)
    !
    USE param
    USE complex_geometry
    USE decomp_2d
    USE variables
    !
    implicit none
    !
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: u
    integer                                            :: i,j,k
    real(mytype)                                       :: x,y,z
    integer                                            :: ix              !=   position du point "zappé"
    integer                                            :: ipif,ipol,nxpif
    integer                                            :: ipoli,ipolf     !=   positions Initiales et Finales du POLynôme considéré
    real(mytype)                                       :: xpol,ypol,dypol   !|variables concernant les polynômes
    real(mytype),dimension(10)                         :: xa,ya           !|de   Lagrange. A mettre impérativement en
    integer                                            :: ia,na             !|double précision
    !
    do k=1,xsize(3)
       do j=1,xsize(2)
          if(nobjx(j,k).ne.0)then
             ia=0
             do i=1,nobjx(j,k)          !boucle sur le nombre d'objets par (j,k)
                !1ère frontière
                nxpif=npif
                ia=ia+1
                xa(ia)=xi(i,j,k)
                ya(ia)=zero
                if(xi(i,j,k) > zero)then!objet immergé
                   ix=xi(i,j,k)/dx+1
                   ipoli=ix+1
                   if(nxipif(i,j,k).lt.npif)nxpif=nxipif(i,j,k)
                   do ipif=1,nxpif
                      ia=ia+1
                      if(izap.eq.1)then!zapping
                         xa(ia)=(ix-1)*dx-ipif*dx
                         ya(ia)=u(ix-ipif,j,k)
                      else             !no zapping
                         xa(ia)=(ix-1)*dx-(ipif-1)*dx
                         ya(ia)=u(ix-ipif+1,j,k)
                      endif
                   enddo
                else                   !objet semi-immergé
                   ipoli=1
                endif
                !2ème frontière
                nxpif=npif
                ia=ia+1
                xa(ia)=xf(i,j,k)
                ya(ia)=zero
                if(xf(i,j,k) < xlx)then!objet immergé
                   ix=(xf(i,j,k)+dx)/dx+1
                   ipolf=ix-1
                   if(nxfpif(i,j,k).lt.npif)nxpif=nxfpif(i,j,k)
                   do ipif=1,nxpif
                      ia=ia+1
                      if(izap.eq.1)then!zapping
                         xa(ia)=(ix-1)*dx+ipif*dx
                         ya(ia)=u(ix+ipif,j,k)
                      else             !no zapping
                         xa(ia)=(ix-1)*dx+(ipif-1)*dx
                         ya(ia)=u(ix+ipif-1,j,k)
                      endif
                   enddo
                else                   !objet semi-immergé
                   ipolf=nx
                endif
                !calcul du polynôme
                na=ia
                do ipol=ipoli,ipolf
                   xpol=dx*(ipol-1)
                   call polint(xa,ya,na,xpol,ypol,dypol)
                   u(ipol,j,k)=ypol
                enddo
                ia=0
             enddo
          endif
       enddo
    enddo
    !
    return
  end subroutine lagpolx
  !############################################################################
  !############################################################################
  subroutine lagpoly(u)
    !
    USE param
    USE complex_geometry
    USE decomp_2d
    USE variables
    !
    implicit none
    !
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: u
    integer                                            :: i,j,k
    real(mytype)                                       :: x,y,z
    integer                                            :: jy              !=   position du point "zappé"
    integer                                            :: jpif,jpol,nypif
    integer                                            :: jpoli,jpolf     !=   positions Initiales et Finales du POLynôme considéré
    real(mytype)                                       :: xpol,ypol,dypol   !|variables concernant les polynômes
    real(mytype),dimension(10)                         :: xa,ya           !|de   Lagrange. A mettre impérativement en
    integer                                            :: ia,na             !|double précision
    !
    do k=1,ysize(3)
       do i=1,ysize(1)
          if(nobjy(i,k).ne.0)then
             ia=0
             do j=1,nobjy(i,k)          !boucle sur le nombre d'objets par (j,k)
                !1ère frontière
                nypif=npif
                ia=ia+1
                xa(ia)=yi(j,i,k)
                ya(ia)=zero
                if(yi(j,i,k) > zero)then!objet immergé
                   jy=1!jy=yi(j,i,k)/dy+1
                   do while(yp(jy) < yi(j,i,k))
                      jy=jy+1
                   enddo
                   jy=jy-1
                   jpoli=jy+1
                   if(nyipif(j,i,k).lt.npif)nypif=nyipif(j,i,k)
                   do jpif=1,nypif
                      ia=ia+1
                      if(izap.eq.1)then!zapping
                         xa(ia)=yp(jy-jpif)!(jy-1)*dy-jpif*dy
                         ya(ia)=u(i,jy-jpif,k)
                      else             !no zapping
                         xa(ia)=yp(jy-jpif+1)!(jy-1)*dy-(jpif-1)*dy
                         ya(ia)=u(i,jy-jpif+1,k)
                      endif
                   enddo
                else                   !objet semi-immergé
                   jpoli=1
                endif
                !2ème frontière
                nypif=npif
                ia=ia+1
                xa(ia)=yf(j,i,k)
                ya(ia)=zero
                if(yf(j,i,k) < yly)then!objet immergé
                   jy=1!jy=(yf(j,i,k)+dy)/dy+1
                   do while(yp(jy).lt.yf(j,i,k))  !there was a bug here yi<-->yf
                      jy=jy+1
                   enddo
                   jpolf=jy-1
                   if(nyfpif(j,i,k).lt.npif)nypif=nyfpif(j,i,k)
                   do jpif=1,nypif
                      ia=ia+1
                      if(izap.eq.1)then!zapping
                         xa(ia)=yp(jy+jpif)!(jy-1)*dy+jpif*dy
                         ya(ia)=u(i,jy+jpif,k)
                      else             !no zapping
                         xa(ia)=yp(jy+jpif-1)!(jy-1)*dy+(jpif-1)*dy
                         ya(ia)=u(i,jy+jpif-1,k)
                      endif
                   enddo
                else                   !objet semi-immergé
                   jpolf=ny
                endif
                !calcul du polynôme
                na=ia
                do jpol=jpoli,jpolf
                   xpol=yp(jpol)!dy*(jpol-1)
                   call polint(xa,ya,na,xpol,ypol,dypol)
                   u(i,jpol,k)=ypol
                enddo
                ia=0
             enddo
          endif
       enddo
    enddo
    !
    return
  end subroutine lagpoly
  !############################################################################
  !############################################################################
  subroutine lagpolz(u)
    !
    USE param
    USE complex_geometry
    USE decomp_2d
    USE variables
    !
    implicit none
    !
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: u
    integer                                            :: i,j,k
    real(mytype)                                       :: x,y,z
    integer                                            :: kz              !=   position du point "zappé"
    integer                                            :: kpif,kpol,nzpif
    integer                                            :: kpoli,kpolf     !=   positions Initiales et Finales du POLynôme considéré
    real(mytype)                                       :: xpol,ypol,dypol   !|variables concernant les polynômes
    real(mytype),dimension(10)                         :: xa,ya           !|de   Lagrange. A mettre imérativement en
    integer                                            :: ia,na             !|double précision
    !
    do j=1,zsize(2)
       do i=1,zsize(1)
          if(nobjz(i,j).ne.0)then
             ia=0
             do k=1,nobjz(i,j)          !boucle sur le nombre d'objets par couple   (i,j)
                !1ère frontière
                nzpif=npif
                ia=ia+1
                xa(ia)=zi(k,i,j)
                ya(ia)=zero
                if(zi(k,i,j) > zero)then!objet immergé
                   kz=zi(k,i,j)/dz+1
                   kpoli=kz+1
                   if(nzipif(k,i,j).lt.npif)nzpif=nzipif(k,i,j)
                   do kpif=1,nzpif
                      ia=ia+1
                      if(izap.eq.1)then!zapping
                         xa(ia)=(kz-1)*dz-kpif*dz
                         ya(ia)=u(i,j,kz-kpif)
                      else             !no zapping
                         xa(ia)=(kz-1)*dz-(kpif-1)*dz
                         ya(ia)=u(i,j,kz-kpif+1)
                      endif
                   enddo
                else                   !objet semi-immergé
                   kpoli=1
                endif
                !2ème frontière
                nzpif=npif
                ia=ia+1
                xa(ia)=zf(k,i,j)
                ya(ia)=zero
                if(zf(k,i,j) < zlz)then!objet immergé
                   kz=(zf(k,i,j)+dz)/dz+1
                   kpolf=kz-1
                   if(nzfpif(k,i,j).lt.npif)nzpif=nzfpif(k,i,j)
                   do kpif=1,nzpif
                      ia=ia+1
                      if(izap.eq.1)then!zapping
                         xa(ia)=(kz-1)*dz+kpif*dz
                         ya(ia)=u(i,j,kz+kpif)
                      else             !no zapping
                         xa(ia)=(kz-1)*dz+(kpif-1)*dz
                         ya(ia)=u(i,j,kz+kpif-1)
                      endif
                   enddo
                else                   !objet semi-immergé
                   kpolf=nz
                endif
                !calcul du polynôme
                na=ia
                do kpol=kpoli,kpolf
                   xpol=dz*(kpol-1)
                   call polint(xa,ya,na,xpol,ypol,dypol)
                   u(i,j,kpol)=ypol
                enddo
                ia=0
             enddo
          endif
       enddo
    enddo
    !
    return
  end subroutine lagpolz
  !############################################################################
  !############################################################################
  subroutine polint(xa,ya,n,x,y,dy)
    !
    use decomp_2d
    use dbg_schemes, only: abs_prec
    !
    implicit none
    !
    integer,parameter            :: nmax=30
    integer                      :: n,i,m,ns
    real(mytype)                 :: dy,x,y,den,dif,dift,ho,hp,w
    real(mytype),dimension(nmax) :: c,d
    real(mytype),dimension(n)    :: xa,ya
    ns=1
    dif=abs_prec(x-xa(1))
    do i=1,n
       dift=abs(x-xa(i))
       if(dift.lt.dif)then
          ns=i
          dif=dift
       endif
       c(i)=ya(i)
       d(i)=ya(i)
    enddo
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       do i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          !         if(den.eq.0)read(*,*)
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
       enddo
       if (2*ns.lt.n-m)then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       endif
       y=y+dy
    enddo
    return
  end subroutine polint
!*--------------------------------------------------------------------------*!
!*                                                                          *!
!*                      Cubic Spline Reconstruction                         *!
!*                                                                          *!
!*--------------------------------------------------------------------------*!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine cubsplx(u,lind)
  !
  USE param
  USE complex_geometry
  USE decomp_2d
  USE variables
  USE ibm_param
  !
  implicit none
  !
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: u
  integer                                            :: i,j,k
  real(mytype)                                       :: x,y,z
  integer                                            :: ix              ! Counter for Skipping Points
  integer                                            :: ipif,ipol,nxpif 
  integer                                            :: ipoli,ipolf     ! Starting and Ending Points for the Reconstruction
  real(mytype)                                       :: xpol,ypol       ! Position and Value of the Reconstructed Solution 
  real(mytype),dimension(10)                         :: xa,ya           ! Position and Value of the Input Data Function 
  integer                                            :: ia,na           
  real(mytype)                                       :: lind            ! Identifying which BC to Impose
  real(mytype)                                       :: bcimp           ! Imposed BC 
  integer                                            :: inxi,inxf
  real(mytype)                                       :: ana_resi,ana_resf         ! Position of Boundary (Analytically)
  !
  ! Initialise Arrays
  xa(:)=0.
  ya(:)=0.
  !
  ! Impose the Correct BC
  bcimp=lind  
  !
  do k=1,xsize(3)
     do j=1,xsize(2)
        if(nobjx(j,k).ne.0)then
           ia=0
           do i=1,nobjx(j,k)          
              !  1st Boundary
              nxpif=npif
              ia=ia+1
              if (ianal.eq.0) then
                 xa(ia)=xi(i,j,k)
                 ana_resi=xi(i,j,k)
              else
                 call analitic_x(j,xi(i,j,k),ana_resi,k) ! Calculate the position of BC analytically
                 xa(ia)=ana_resi
              endif
              ya(ia)=bcimp
              if(xi(i,j,k).gt.0.)then ! Immersed Object
                 inxi=0
                 ix=xi(i,j,k)/dx+1
                 ipoli=ix+1
                 if(nxipif(i,j,k).lt.npif)nxpif=nxipif(i,j,k)
                 do ipif=1,nxpif
                    ia=ia+1
                    if(izap.eq.1)then ! Skip First Points
                       xa(ia)=(ix-1)*dx-ipif*dx
                       ya(ia)=u(ix-ipif,j,k)
                    else              ! Don't Skip any Points
                       xa(ia)=(ix-1)*dx-(ipif-1)*dx
                       ya(ia)=u(ix-ipif+1,j,k)
                    endif
                 enddo
              else                    ! Boundary Coincides with Physical Boundary (Inlet)
                 inxi=1
                 ipoli=1
                 ix=xi(i,j,k)/dx
                 ipoli=ix+1
                 if(nxipif(i,j,k).lt.npif)nxpif=nxipif(i,j,k)
                 do ipif=1,nxpif
                    ia=ia+1
                    if(izap.eq.1)then ! Skip First Points
                       xa(ia)=(ix-1)*dx-ipif*dx
                       ya(ia)=bcimp                   
                    else              ! Don't Skip any Points
                       xa(ia)=(ix-1)*dx-(ipif-1)*dx
                       ya(ia)=bcimp      
                    endif
                 enddo
              endif
              !
              !  2nd Boundary
              nxpif=npif
              ia=ia+1
              if (ianal.eq.0) then
                 xa(ia)=xf(i,j,k)
                 ana_resf=xf(i,j,k)
              else
                 call analitic_x(j,xf(i,j,k),ana_resf,k) ! Calculate the position of BC analytically
                 xa(ia)=ana_resf
              endif              
              ya(ia)=bcimp              
              if(xf(i,j,k).lt.xlx)then ! Immersed Object
                 inxf=0
                 ix=(xf(i,j,k)+dx)/dx+1
                 ipolf=ix-1
                 if(nxfpif(i,j,k).lt.npif)nxpif=nxfpif(i,j,k)
                 do ipif=1,nxpif
                    ia=ia+1
                    if(izap.eq.1)then  ! Skip First Points
                       xa(ia)=(ix-1)*dx+ipif*dx
                       ya(ia)=u(ix+ipif,j,k)
                    else               ! Don't Skip any Points
                       xa(ia)=(ix-1)*dx+(ipif-1)*dx
                       ya(ia)=u(ix+ipif-1,j,k)
                    endif
                 enddo
              else                     ! Boundary Coincides with Physical Boundary (Outlet)
                 inxf=1
                 ipolf=nx
                 ix=(xf(i,j,k)+dx)/dx+1
                 ipolf=ix-1
                 if(nxfpif(i,j,k).lt.npif)nxpif=nxfpif(i,j,k)
                 do ipif=1,nxpif
                    ia=ia+1
                    if(izap.eq.1)then  ! Skip First Points
                       xa(ia)=(ix-1)*dx+ipif*dx
                       ya(ia)=bcimp                                   
                    else               ! Don't Skip any Points
                       xa(ia)=(ix-1)*dx+(ipif-1)*dx
                       ya(ia)=bcimp                                   
                    endif
                 enddo
              endif
              ! Special Case
              if (xi(i,j,k).eq.xf(i,j,k)) then
                  u(ipol,j,k)=bcimp                                   
              else
              ! Cubic Spline Reconstruction
		  na=ia
		  do ipol=ipoli,ipolf
		     if ((inxf.eq.1).and.(inxi.eq.1)) then ! If the Body Extends from the Inlet to the Outlet (Special Case)
                 u(ipol,j,k)=bcimp                            
             else
		         xpol=dx*(ipol-1)
		         if (xpol.eq.ana_resi) then
		            u(ipol,j,k)=bcimp
		         elseif (xpol.eq.ana_resf) then
		            u(ipol,j,k)=bcimp
		         else   
		            call cubic_spline(xa,ya,na,xpol,ypol)
		            u(ipol,j,k)=ypol
		         endif
		     endif
		  enddo
		  ia=0
	      endif    
           enddo
        endif
     enddo
  enddo
  !
  return
end subroutine cubsplx
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine cubsply(u,lind)
  !
  USE param
  USE complex_geometry
  USE decomp_2d
  USE variables
  USE ibm_param
  !
  implicit none
  !
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: u
  integer                                            :: i,j,k
  real(mytype)                                       :: x,y,z
  integer                                            :: jy              ! Counter for Skipping Points
  integer                                            :: jpif,jpol,nypif
  integer                                            :: jpoli,jpolf     ! Starting and Ending Points for the Reconstruction
  real(mytype)                                       :: xpol,ypol,dypol ! Position and Value of the Reconstructed Solution 
  real(mytype),dimension(10)                         :: xa,ya           ! Position and Value of the Input Data Function 
  integer                                            :: ia,na           
  real(mytype)                                       :: lind            ! Identifying which BC to Impose
  real(mytype)                                       :: bcimp           ! Imposed BC 
  integer                                            :: inxi,inxf  
  real(mytype)                                       :: ana_resi,ana_resf
  !
  ! Initialise Arrays
  xa(:)=0.
  ya(:)=0.
  !
  ! Impose the Correct BC
  bcimp=lind  
  !
  do k=1,ysize(3)
     do i=1,ysize(1)
        if(nobjy(i,k).ne.0)then
           ia=0
           do j=1,nobjy(i,k)          
              !  1st Boundary
              nypif=npif
              ia=ia+1
              if (ianal.eq.0) then
                 xa(ia)=yi(j,i,k)
                 ana_resi=yi(j,i,k)
              else
                 call analitic_y(i,yi(j,i,k),ana_resi,k) ! Calculate the position of BC analytically
                 xa(ia)=ana_resi
              endif  
              ya(ia)=bcimp
              if(yi(j,i,k).gt.0.)then ! Immersed Object
                 jy=1
                 do while(yp(jy).lt.yi(j,i,k))
                    jy=jy+1
                 enddo
                 jy=jy-1
                 jpoli=jy+1
                 if(nyipif(j,i,k).lt.npif)nypif=nyipif(j,i,k)
                 do jpif=1,nypif
                    ia=ia+1
                    if(izap.eq.1)then ! Skip First Points
                       xa(ia)=yp(jy-jpif)
                       ya(ia)=u(i,jy-jpif,k)
                    else              ! Don't Skip any Points
                       xa(ia)=yp(jy-jpif+1)
                       ya(ia)=u(i,jy-jpif+1,k)
                    endif
                 enddo
              else                   ! Boundary Coincides with Physical Boundary (Bottom)
                 jy=1
                 jpoli=1
                 do while(yp(jy).lt.yi(j,i,k))
                    jy=jy+1
                 enddo
                 jy=jy-1
                 jpoli=jy+1
                 if(nyipif(j,i,k).lt.npif)nypif=nyipif(j,i,k)
                 do jpif=1,nypif
                    ia=ia+1
                    if(izap.eq.1)then ! Skip First Points
                       xa(ia)=yp(1)-(jpif+1)*dy
                       ya(ia)=bcimp
                    else              ! Don't Skip any Points
                       xa(ia)=yp(1)-(jpif*dy)
                       ya(ia)=bcimp
                    endif
                 enddo
              endif
              ! 2nd Boundary
              nypif=npif
              ia=ia+1
              if (ianal.eq.0) then
                 xa(ia)=yf(j,i,k)
                 ana_resf=yf(j,i,k)
              else
                 call analitic_y(i,yf(j,i,k),ana_resf,k) ! Calculate the position of BC analytically
                 xa(ia)=ana_resf
              endif  
              ya(ia)=bcimp
              if(yf(j,i,k).lt.yly)then ! Immersed Object
                 jy=1
                 do while(yp(jy).lt.yf(j,i,k))  
                    jy=jy+1
                 enddo
                 jpolf=jy-1
                 if(nyfpif(j,i,k).lt.npif)nypif=nyfpif(j,i,k)
                 do jpif=1,nypif
                    ia=ia+1
                    if(izap.eq.1)then ! Skip First Points
                       xa(ia)=yp(jy+jpif)
                       ya(ia)=u(i,jy+jpif,k)
                    else              ! Don't Skip any Points
                       xa(ia)=yp(jy+jpif-1)
                       ya(ia)=u(i,jy+jpif-1,k)
                    endif
                 enddo
              else                   ! Boundary Coincides with Physical Boundary (Top)
                 jy=1
                 jpolf=ny
                 do while(yp(jy).lt.yf(j,i,k))  
                    jy=jy+1
                 enddo
                 jpolf=jy-1
                 if(nyfpif(j,i,k).lt.npif)nypif=nyfpif(j,i,k)
                 do jpif=1,nypif
                    ia=ia+1
                    if(izap.eq.1)then ! Skip First Points
                       xa(ia)=yp(ny)+(jpif+1)*dy
                       ya(ia)=bcimp
                    else              ! Don't Skip any Points
                       xa(ia)=yp(ny)+(jpif)*dy
                       ya(ia)=bcimp
                    endif
                 enddo
              endif
              ! Special Case
              if (yi(j,i,k).eq.yf(j,i,k)) then
                  u(i,jpol,k)=bcimp                                   
              else
		  !calcul du polynôme
		   na=ia
		   do jpol=jpoli,jpolf
		         xpol=yp(jpol)
		         if (xpol.eq.ana_resi) then
		            u(i,jpol,k)=bcimp
		         elseif (xpol.eq.ana_resf) then
		            u(i,jpol,k)=bcimp
		         else   
		            call cubic_spline(xa,ya,na,xpol,ypol)
		            u(i,jpol,k)=ypol
		         endif
		   enddo
		   ia=0
	      endif    
           enddo
        endif
     enddo
  enddo
  !
  return
end subroutine cubsply
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine cubsplz(u,lind)
  !
  USE param
  USE complex_geometry
  USE decomp_2d
  USE variables
  USE ibm_param
  !
  implicit none
  !
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: u
  integer                                            :: i,j,k
  real(mytype)                                       :: x,y,z
  integer                                            :: kz              != position du point "zappé"
  integer                                            :: kpif,kpol,nzpif
  integer                                            :: kpoli,kpolf     != positions Initiales et Finales du POLynôme considéré
  real(mytype)                                       :: xpol,ypol,dypol !|variables concernant les polynômes
  real(mytype),dimension(10)                         :: xa,ya           !|de Lagrange. A mettre imérativement en 
  integer                                            :: ia,na           !|double précision
  real(mytype)                                       :: lind            ! Identifying which BC to Impose
  real(mytype)                                       :: bcimp           ! Imposed BC 
  integer                                            :: inxi,inxf  
  real(mytype)                                       :: ana_resi,ana_resf
  !
  ! Initialise Arrays
  xa(:)=zero
  ya(:)=zero
  !
  ! Impose the Correct BC
  bcimp=lind  
  !
  do j=1,zsize(2)
     do i=1,zsize(1)
        if(nobjz(i,j).ne.0)then
           ia=0
           do k=1,nobjz(i,j)          
              !  1st Boundary
              nzpif=npif
              ia=ia+1
              if (ianal.eq.0) then
                 xa(ia)=zi(k,i,j)
                 ana_resi=zi(k,i,j)
              else
!                 call analitic_z(i,zi(k,i,j),ana_resi,j) ! Calculate the position of BC analytically
                 xa(ia)=ana_resi
              endif  
              ya(ia)=bcimp
              if(zi(k,i,j).gt.0.)then ! Immersed Object
                 inxi=0
                 kz=zi(k,i,j)/dz+1
                 kpoli=kz+1
                 if(nzipif(k,i,j).lt.npif)nzpif=nzipif(k,i,j)
                 do kpif=1,nzpif
                    ia=ia+1
                    if(izap.eq.1)then ! Skip First Points
                       xa(ia)=(kz-1)*dz-kpif*dz
                       ya(ia)=u(i,j,kz-kpif)
                    else              ! Don't Skip any Points
                       xa(ia)=(kz-1)*dz-(kpif-1)*dz
                       ya(ia)=u(i,j,kz-kpif+1)
                    endif
                 enddo
              else                   ! Boundary Coincides with Physical Boundary (Front) 
                 inxi=1
                 kz=zi(k,i,j)/dz
                 kpoli=1
                 if(nzipif(k,i,j).lt.npif)nzpif=nzipif(k,i,j)
                 do kpif=1,nzpif
                    ia=ia+1
                    if(izap.eq.1)then ! Skip First Points
                       xa(ia)=(kz-1)*dz-kpif*dz
                       ya(ia)=bcimp
                    else              ! Don't Skip any Points
                       xa(ia)=(kz-1)*dz-(kpif-1)*dz
                       ya(ia)=bcimp
                    endif
                 enddo
              endif
              !  2nd Boundary
              nzpif=npif
              ia=ia+1
              if (ianal.eq.0) then
                 xa(ia)=zf(k,i,j)
                 ana_resf=zf(k,i,j)
              else
                 !call analitic_z(i,zf(k,i,j),ana_resf,j) ! Calculate the position of BC analytically
                 xa(ia)=ana_resf
              endif
              ya(ia)=bcimp
              if(zf(k,i,j).lt.zlz)then  ! Immersed Object
                 inxf=0
                 kz=(zf(k,i,j)+dz)/dz+1
                 kpolf=kz-1
                 if(nzfpif(k,i,j).lt.npif)nzpif=nzfpif(k,i,j)
                 do kpif=1,nzpif
                    ia=ia+1
                    if(izap.eq.1)then  ! Skip First Points
                       xa(ia)=(kz-1)*dz+kpif*dz
                       ya(ia)=u(i,j,kz+kpif)
                    else               ! Don't Skip any Points
                       xa(ia)=(kz-1)*dz+(kpif-1)*dz
                       ya(ia)=u(i,j,kz+kpif-1)
                    endif
                 enddo
              else                     ! Boundary Coincides with Physical Boundary (Back)
                 inxf=1
                 kz=(zf(k,i,j)+dz)/dz+1
                 kpolf=nz
                 if(nzfpif(k,i,j).lt.npif)nzpif=nzfpif(k,i,j)
                 do kpif=1,nzpif
                    ia=ia+1
                    if(izap.eq.1)then  ! Skip First Points
                       xa(ia)=(kz-1)*dz+kpif*dz
                       ya(ia)=bcimp
                    else               ! Don't Skip any Points
                       xa(ia)=(kz-1)*dz+(kpif-1)*dz
                       ya(ia)=bcimp
                    endif
                 enddo
              endif
         !     ! Special Case
         !     if (zi(k,i,j).eq.zf(k,i,j)) then
         !         u(i,j,kpol)=bcimp                                   
         !     else              
    	      ! Cubic Spline Reconstruction
	            na=ia
	            do kpol=kpoli,kpolf 
                          ! Special Case
                          if (zi(k,i,j).eq.zf(k,i,j)) then
                                u(i,j,kpol)=bcimp
                          else
	                    if ((inxf.eq.1).and.(inxi.eq.1)) then ! If the Body Extends from the Front to the Back (Special Case)
	                          u(i,j,kpol)=bcimp                            
	                     else              
	                          xpol=dz*(kpol-1)
	                          call cubic_spline(xa,ya,na,xpol,ypol)
	                          u(i,j,kpol)=ypol
	                     endif
                           endif
	            enddo
	            ia=0
         !     endif
           enddo
        endif
     enddo
  enddo
  !
  return
end subroutine cubsplz
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine cubic_spline(xa,ya,n,x,y)
   !
   use decomp_2d
   use param, only : zero, two, three, zpfive
   !
   implicit none
   !
   integer                      :: n,i,j,nc,nk
   real(mytype)                 :: x,y,xcc
   real(mytype),dimension(n)    :: xa,ya
   real(mytype),dimension(10)   :: xaa,yaa
   real(mytype)                 :: ypri,yprf
   real(mytype),dimension(n-2)  :: xx,alpha,cc,zz,ll,aa,yy
   real(mytype),dimension(n-3)  :: hh,dd,bb,mm
     !
   ! Initialise Arrays
   xaa(:)=zero
   yaa(:)=zero
   ! Arrange Points in Correct Order (based on x-coor)
   j=n/2
   do i=1,n
      if (i <= n/2) then
         xaa(i)=xa(j)
         yaa(i)=ya(j)
         j=j-1
      else
         xaa(i)=xa(i)
         yaa(i)=ya(i)
      endif
   enddo
   !
   ypri=(yaa(3)-yaa(1))/(xaa(3)-xaa(1))
   yprf=(yaa(n)-yaa(n-2))/(xaa(n)-xaa(n-2))
   !
   nk=n-1
   !
   do i=2,nk
      yy(i-1)=yaa(i)
   enddo
   !
   do i=2,nk
      xx(i-1)=xaa(i)
   enddo
   !
   nc=nk-1
   !
   do i=1,nc
      aa(i)=yy(i)
   enddo
   !
   do i=1,nc-1
      hh(i)=xx(i+1)-xx(i)
   enddo
   !
   alpha(1)=(three*(aa(2)-aa(1)))/hh(1) - three*ypri
   alpha(nc)= three*yprf - three*(aa(nc)-aa(nc-1))/hh(nc-1)
   !
   do i=2,nc-1
      alpha(i)=(three/hh(i))*(aa(i+1)-aa(i))-(three/hh(i-1))*(aa(i)-aa(i-1))
   enddo
   ll(1)=two*hh(1)
   mm(1)=zpfive
   zz(1)=alpha(1)/ll(1)
   !
   do i=2,nc-1
      ll(i)=two*(xx(i+1)-xx(i-1))-hh(i-1)*mm(i-1);
      mm(i)=hh(i)/ll(i);
      zz(i)=(alpha(i)-hh(i-1)*zz(i-1))/ll(i);
   enddo
   !
   ll(nc)=hh(nc-1)*(two-mm(nc-1));
   zz(nc)=(alpha(nc)-hh(nc-1)*zz(nc-1))/ll(nc);
   cc(nc)=zz(nc);
   !
   do j=nc-1,1,-1
      cc(j)=zz(j)-mm(j)*cc(j+1);
      bb(j)=(aa(j+1)-aa(j))/hh(j)-(hh(j)/three)*(cc(j+1)+two*cc(j));
      dd(j)=(cc(j+1)-cc(j))/(three*hh(j));
   enddo
   !
   do j=2,nc
      xcc=x;
      if (xcc <= xx(j) .and. xcc >= xx(j-1)) then
         y= aa(j-1) + bb(j-1)*(xcc-xx(j-1)) + cc(j-1)*(xcc-xx(j-1))**2 + dd(j-1)*(xcc-xx(j-1))**3;
      endif
   enddo
   !
   return
end subroutine cubic_spline
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine ana_y_cyl(i,y_pos,ana_res)
  !
  USE param
  USE complex_geometry
  USE decomp_2d
  USE variables
  USE ibm_param
  !
  implicit none
  !
  integer                                            :: i
  real(mytype)                                       :: y_pos,ana_res 
  real(mytype)                                       :: cexx,ceyy
  !
  if (t.ne.0.) then
     cexx = cex + ubcx*(t-ifirst*dt)
     ceyy = cey + ubcy*(t-ifirst*dt)
  else
     cexx = cex
     ceyy = cey
  endif
  if (y_pos.gt.ceyy) then     ! Impose analytical BC
      ana_res=ceyy + sqrt(ra**2.0-((i+ystart(1)-1-1)*dx-cexx)**2.0)
  else
      ana_res=ceyy - sqrt(ra**2.0-((i+ystart(1)-1-1)*dx-cexx)**2.0)
  endif     
  !
  return
end subroutine ana_y_cyl
!***************************************************************************
!
subroutine ana_x_cyl(j,x_pos,ana_res)
  !
  USE param
  USE complex_geometry
  USE decomp_2d
  USE variables
  USE ibm_param
  !
  implicit none
  !
  integer                                            :: j
  real(mytype)                                       :: x_pos,ana_res 
  real(mytype)                                       :: cexx,ceyy
  !
  if (t.ne.0.) then
     cexx = cex + ubcx*(t-ifirst*dt)
     ceyy = cey + ubcy*(t-ifirst*dt)
  else
     cexx = cex
     ceyy = cey
  endif
  if (x_pos.gt.cexx) then     ! Impose analytical BC
      ana_res = cexx + sqrt(ra**2.0-(yp(j+xstart(2)-1)-ceyy)**2.0)
  else
      ana_res = cexx - sqrt(ra**2.0-(yp(j+xstart(2)-1)-ceyy)**2.0)
  endif     
  !
  return
end subroutine ana_x_cyl
!*******************************************************************
SUBROUTINE analitic_x(j,x_pos,ana_res,k)

  USE param, ONLY : itype, itype_cyl
  USE decomp_2d, ONLY : mytype
!  USE cyl, ONLY : geomcomplex_cyl

  IMPLICIT NONE

  integer                                            :: j,k
  real(mytype)                                       :: x_pos,ana_res 

  IF (itype.EQ.itype_cyl) THEN

     CALL ana_x_cyl(j,x_pos,ana_res)

  ENDIF

END SUBROUTINE analitic_x
!*******************************************************************
!*******************************************************************
SUBROUTINE analitic_y(i,y_pos,ana_res,k)

  USE param, ONLY : itype, itype_cyl
  USE decomp_2d, ONLY : mytype
!  USE cyl, ONLY : geomcomplex_cyl

  IMPLICIT NONE

  integer                                            :: i,k
  real(mytype)                                       :: y_pos,ana_res 

  IF (itype.EQ.itype_cyl) THEN

     CALL ana_y_cyl(i,y_pos,ana_res)

  ENDIF

END SUBROUTINE analitic_y
!*******************************************************************
!*******************************************************************
  
  
end module ibm
