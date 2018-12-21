subroutine stabiltemp() !from Erik, adapted by Leonardo Romero Monteiro

  use param
  use variables
  use var

  implicit none

  complex(mytype) :: z,eit,ei2t,ei3t,eimt,eim2t,eim3t
  real(mytype) :: theta, dtheta, cc, fourier, cfl
  real(mytype) :: xkm, xk, xkp, xks, xkf, x, y
  real(mytype) :: am1, a0, a1, a2, a3
  real(mytype) :: bm1, b0, b1, b2, b3
  real(mytype) :: alpha1, c1, c11
  real(mytype) :: alpha2, c2
  real(mytype) :: alpha3, beta3, c3, d3
  integer :: i,ntheta,order

  ntheta=360
  dtheta=twopi/(ntheta-1.)
  xk=(fpi2+1.)*pi*pi
  order = 6   ! ordem da hiperviscosidade 0 = sem hiperviscosidade; 4 = 4a ordem com 2 formados;  6 = 6a ordem com 1 formado

  print *,'Writing stability data!'

  if (itimescheme==0) then !Euler (not implemented)
     am1=0; a0=1.; a1=0.; a2=0.
  endif

  if (itimescheme.eq.1) then !AB2
     am1=0; a0=1.5; a1=-0.5; a2=0.; a3=0.; bm1=1.; b0=-1.; b1=0.; b2=0.; b3=0.
  endif

  if (itimescheme.eq.3) then !RK3
     if (nrank==0) write(*,*) "Non implemented for RK3"
  endif

  if (itimescheme.eq.2) then !AB3
     am1=0.; a0=23./12.; a1=-16./12.; a2=5./12; a0=3./2+a2; a1=-1./2-2*a2; a3=0.; bm1=1.; b0=-1.; b1=0.; b2=0.; b3=0.
  endif

  open(10,file='stabiltemp_1.dat',form='formatted')
  do i=1,ntheta
     theta=(i-1)*dtheta

     eit=exp(cmplx(0.,1.)*theta)
     ei2t=eit*eit
     ei3t=eit*eit*eit
     eimt=1./eit
     eim2t=1./ei2t
     eim3t=1./ei3t
     !z=(eit-1.)/a0
     !z=(eit*(eit-1.))/(a0*eit+a1)
     !z=(ei3t-ei2t)/(a0*ei2t+a1*eit+a2)
     z=(bm1*eit+b0+b1*eimt+b2*eim2t+b3*eim3t)/(a0+a1*eimt+a2*eim2t+a3*eim3t)
     !z=(eit-1.)/(am1*eit+a0+a1*eimt)
     !z=(eit-1.)/(am1*eit+a0+a1*eimt+a2*eim2t)

     write(10,*) real(z),imag(z)
  enddo
  close(10)


  alpha1=1./3.
  a1=(alpha1+9.)/6.
  b1=(32.*alpha1-9.)/15.
  c1=(-3.*alpha1+1.)/10.

  if (order.eq.0) then

     alpha2=2./11
     a2=12./11
     b2=3./11
     c2=0.

  elseif (order.eq.4) then

     c11=exp(-((pi-2.*pi/3.)/(0.3*pi-2.*pi/3.))**2 )
     xkm=(c11*fpi2+1.)*(4./9.)*pi*pi

     alpha2=(64.*xkm-27.*xk-96.)/(64.*xkm-54.*xk+48.)
     a2 = (54.*xk-15.*xkm*xk+12.)/(64.*xkm-54.*xk+48.)
     b2 = (192.*xkm-216.*xk+24.*xkm*xk-48.)/(64.*xkm-54.*xk+48.)
     c2 = 3.*(18.*xk -3.*xkm*xk-36.)/(64.*xkm-54.*xk+48.)

  elseif(order.eq.6) then

     alpha2=(45.*xk-272.)/(2*(45.*xk-208.))
     c2=(2.-11.*alpha2)/20.
     a2=(6.-9.*alpha2)/4.
     b2=(-3.+24.*alpha2)/5.

  endif

  !alpha3=0.45
  !beta3=(3.-2.*alpha3)/10.
  !a3=(2.+3.*alpha3)/4.
  !b3=(6.+7*alpha3)/8.
  !c3=(6.+alpha3)/20.
  !d3=(2-3.*alpha3)/40.

  cc=4.
  fourier=xnu*dt/(dx*dx)
  cfl=cc*dt/dx

  open(10,file='stabiltemp_2.dat',form='formatted')
  do i=1,ntheta
     theta=(i-1)*dtheta

     xkp=(a1*sin(theta)+(b1/2)*sin(2*theta) +(c1/3)*sin(3*theta))/(1+2*alpha1*cos(theta))
     xks=(2*a2*(1-cos(theta))+(b2/2)*(1-cos(2*theta)) +(2*c2/9)*(1-cos(3*theta)))/(1+2*alpha2*cos(theta))
     !xkf=(a3+b3*cos(theta)+c3*cos(2*theta)+d3*cos(3*theta)) /(1+2*alpha3*cos(theta)+2*beta3*cos(2*theta))
     x=-fourier*xks
     y=-cfl*xkp!*xkf

     write(10,*) x,y
  enddo
  close(10)

end subroutine stabiltemp
!*******************************************************************
!===================================================
! Subroutine for computing the local and global CFL
! number, according to Lele 1992.
!===================================================
subroutine cfl_compute(uxmax,uymax,uzmax)

use param
use variables
use var
    
implicit none

real(mytype),intent(in) :: uxmax,uymax,uzmax
real(mytype) :: cfl_x_adv,cfl_x_diff,cfl_y_adv,cfl_y_diff,cfl_z_adv,cfl_z_diff
real(mytype) :: cfl_conv_lim, cfl_diff_lim
real(mytype) :: sigma_conv(3), sigma_diff(3)
real(mytype) :: visc

! Set the constants (this is true for periodic boundaries)
sigma_conv=[0.0, sqrt(3.0), 2.85]
sigma_diff=[2.0, 2.5, 2.9]

if(jles==0) then
visc=xnu
elseif (jles==1) then
visc=20*fpi2*xnu
endif

! This is considering 1D peridic boundaries
! Do x-direction
cfl_x_adv=abs(uxmax)*dt/dx
cfl_x_diff=visc*dt/dx**2.0
! Do y-direction
cfl_y_adv=abs(uymax)*dt/dy
cfl_y_diff=visc*dt/dy**2.0
! Do z-direction
cfl_z_adv=abs(uzmax)*dt/dz
cfl_z_diff=visc*dt/dz**2.0

! So far we will focus on uniform grids
if(nrank==0) then
    write(*,*) ' '
    write(*,1002) cfl_x_adv, cfl_x_diff
1002  format('CFL x-direction (Adv and Diff) =',F9.4,',',F9.4)
    write(*,1003) cfl_y_adv, cfl_y_diff
1003  format('CFL y-direction (Adv and Diff) =',F9.4,',',F9.4)
    write(*,1004) cfl_z_adv, cfl_z_diff
1004  format('CFL z-direction (Adv and Diff) =',F9.4,',',F9.4)
    cfl_conv_lim=sigma_conv(itimescheme)/sqrt(3.0)
    cfl_diff_lim=sigma_diff(itimescheme)/6.0
    write(*,1005) cfl_conv_lim, cfl_diff_lim
    write(*,*) ' '
1005  format('CFL limits (Adv and Diff) : ',F9.4,',',F9.4)
endif

end subroutine
!*******************************************************************
subroutine test_scalar_min_max(phi)

USE decomp_2d
USE decomp_2d_poisson
USE variables
USE param
USE var
USE MPI

implicit none

integer :: code,ierror,ijk,is
real(mytype) :: phimax,phimin,phimax1,phimin1
real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

do is=1, numscalar
 phimax=-1609.;phimin=1609.
 do ijk=1,xsize(1)*xsize(2)*xsize(3)
   if (phi(ijk,1,1,is).gt.phimax) phimax=phi(ijk,1,1,is)
   if (phi(ijk,1,1,is).lt.phimin) phimin=phi(ijk,1,1,is)
 enddo

 call MPI_REDUCE(phimax,phimax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
 call MPI_REDUCE(phimin,phimin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)

 if (nrank.eq.0) then

  print *,'Phi'//char(48+is)//' min max=', real(phimin1,4), real(phimax1,4)

  if (abs(phimax1).ge.10.) then !if phi control turned off
   stop 'Scalar diverged! FATALITY!'
  endif
 endif

enddo

return
end subroutine test_scalar_min_max
!*******************************************************************
subroutine test_speed_min_max(ux,uy,uz)

USE decomp_2d
USE decomp_2d_poisson
USE variables
USE param
USE var
USE MPI

implicit none

integer :: code,ierror,ijk
real(mytype) :: uxmax,uymax,uzmax,uxmin,uymin,uzmin
real(mytype) :: uxmax1,uymax1,uzmax1,uxmin1,uymin1,uzmin1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz

uxmax=-1609.;uymax=-1609.;uzmax=-1609.;uxmin=1609.;uymin=1609.;uzmin=1609.
do ijk=1,xsize(1)*xsize(2)*xsize(3)
#ifdef IBM
    if (ep1(ijk,1,1).eq.0.) then
#endif
   if (ux(ijk,1,1).gt.uxmax) uxmax=ux(ijk,1,1)
   if (uy(ijk,1,1).gt.uymax) uymax=uy(ijk,1,1)
   if (uz(ijk,1,1).gt.uzmax) uzmax=uz(ijk,1,1)
   if (ux(ijk,1,1).lt.uxmin) uxmin=ux(ijk,1,1)
   if (uy(ijk,1,1).lt.uymin) uymin=uy(ijk,1,1)
   if (uz(ijk,1,1).lt.uzmin) uzmin=uz(ijk,1,1)
#ifdef IBM
   endif
#endif
enddo

call MPI_REDUCE(uxmax,uxmax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(uymax,uymax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(uzmax,uzmax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(uxmin,uxmin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(uymin,uymin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(uzmin,uzmin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)

if (nrank.eq.0) then

  print *,'U,V,W min=',real(uxmin1,4),real(uymin1,4),real(uzmin1,4)
  print *,'U,V,W max=',real(uxmax1,4),real(uymax1,4),real(uzmax1,4)
 !print *,'CFL=',real(abs(max(uxmax1,uymax1,uzmax1)*dt)/min(dx,dy,dz),4)

  if((abs(uxmax1).ge.10.).OR.(abs(uymax1).ge.10.).OR.(abs(uzmax1).ge.10.)) then
   stop 'Velocity diverged! FATALITY!'
  endif

endif

 !if (mod(itime,imodulo).eq.0) call cfl_compute(uxmax,uymax,uzmax)

return
end subroutine test_speed_min_max
!*******************************************************************
subroutine restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3,phi1,px1,py1,pz1,iresflg)

  USE decomp_2d
  USE decomp_2d_io
  USE variables
  USE param
  USE MPI

  implicit none

  integer :: i,j,k,iresflg,nzmsize,fh,ierror,is,code
  integer :: ierror_o=0 !error to open sauve file during restart
  real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
  real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: px1,py1,pz1
  real(mytype), dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1
  real(mytype), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
  real(mytype), dimension(phG%zst(1):phG%zen(1),phG%zst(2):phG%zen(2),phG%zst(3):phG%zen(3)) :: pp3
  integer (kind=MPI_OFFSET_KIND) :: filesize, disp
  real(mytype) :: xdt
  integer, dimension(2) :: dims, dummy_coords
  logical, dimension(2) :: dummy_periods
  character(len=30) :: filename, filestart

  if (iresflg .eq. 1 ) then !Writing restart
     if (nrank==0) print *,'===========================================================<<<<<'
     if (nrank==0) print *,'Writing restart point',itime/icheckpoint
     !if (nrank==0) print *,'File size',real((s3df*16.)*1e-9,4),'GB'
  end if

  write(filename,"('restart',I7.7)") itime
  write(filestart,"('restart',I7.7)") ifirst-1

  if (iresflg==1) then !write
     call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
          MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
          fh, ierror)
     filesize = 0_MPI_OFFSET_KIND
     call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
     disp = 0_MPI_OFFSET_KIND
     call decomp_2d_write_var(fh,disp,1,ux1)
     call decomp_2d_write_var(fh,disp,1,uy1)
     call decomp_2d_write_var(fh,disp,1,uz1)
     call decomp_2d_write_var(fh,disp,1,ep1)
     do is=1, ntime
        call decomp_2d_write_var(fh,disp,1,dux1(:,:,:,is))
        call decomp_2d_write_var(fh,disp,1,duy1(:,:,:,is))
        call decomp_2d_write_var(fh,disp,1,duz1(:,:,:,is))
     end do
     call decomp_2d_write_var(fh,disp,1,px1)
     call decomp_2d_write_var(fh,disp,1,py1)
     call decomp_2d_write_var(fh,disp,1,pz1)
     call decomp_2d_write_var(fh,disp,1,pp3,phG)
     if (iscalar==1) then
        do is=1, numscalar
           call decomp_2d_write_var(fh,disp,1,phi1(:,:,:,is))
        end do
     endif
     call MPI_FILE_CLOSE(fh,ierror)
  else
     if (nrank==0) print *,'RESTART from file:', filestart
     call MPI_FILE_OPEN(MPI_COMM_WORLD, filestart, &
          MPI_MODE_RDONLY, MPI_INFO_NULL, &
          fh, ierror_o)
     disp = 0_MPI_OFFSET_KIND
     call decomp_2d_read_var(fh,disp,1,ux1)
     call decomp_2d_read_var(fh,disp,1,uy1)
     call decomp_2d_read_var(fh,disp,1,uz1)
     call decomp_2d_read_var(fh,disp,1,ep1)     
     do is=1, ntime
        call decomp_2d_read_var(fh,disp,1,dux1(:,:,:,is))
        call decomp_2d_read_var(fh,disp,1,duy1(:,:,:,is))
        call decomp_2d_read_var(fh,disp,1,duz1(:,:,:,is))
     end do
     call decomp_2d_read_var(fh,disp,1,px1)
     call decomp_2d_read_var(fh,disp,1,py1)
     call decomp_2d_read_var(fh,disp,1,pz1)
     call decomp_2d_read_var(fh,disp,1,pp3,phG)
     if (iscalar==1) then
        do is=1, numscalar
           call decomp_2d_read_var(fh,disp,1,phi1(:,:,:,is))
        end do
     endif
     call MPI_FILE_CLOSE(fh,ierror_o)
  endif

  if (nrank.eq.0) then
     if (ierror_o .ne. 0) then !Included by Felipe Schuch
        print *,'==========================================================='
        print *,'Error: Impossible to read '//trim(filestart)
        print *,'==========================================================='
        call MPI_ABORT(MPI_COMM_WORLD,code,ierror)
     endif
  endif

  if (iresflg==0) then
     ! reconstruction of the dp/dx, dp/dy and dp/dz from px1,py1 and pz1
     ! Temporal scheme (1:EULER, 2:AB2, 3: AB3, 4:AB4, 5:RK3, 6:RK4)
     if (itimescheme.eq.1) then
        xdt=gdt(1)
     elseif (itimescheme.eq.2) then
        xdt=gdt(1)
     elseif (itimescheme.eq.3) then
        xdt = gdt(1)
     elseif (itimescheme.eq.5) then
        xdt=gdt(3)
     else
        if (nrank.eq.0) then
           print *, "Timescheme not implemented!"
           stop
        endif
     endif

     do k=1,xsize(3)
        do j=1,xsize(2)
           dpdyx1(j,k)=py1(1,j,k)/xdt
           dpdzx1(j,k)=pz1(1,j,k)/xdt
           dpdyxn(j,k)=py1(nx,j,k)/xdt
           dpdzxn(j,k)=pz1(nx,j,k)/xdt
        enddo
     enddo

     if (xsize(3)==1) then
        do j=1,xsize(2)
           do i=1,xsize(1)
              dpdxz1(i,j)=px1(i,j,1)/xdt
              dpdyz1(i,j)=py1(i,j,1)/xdt
           enddo
        enddo
     endif
     if (xsize(3)==nz) then
        do j=1,xsize(2)
           do i=1,xsize(1)
              dpdxzn(i,j)=px1(i,j,nz)/xdt
              dpdyzn(i,j)=py1(i,j,nz)/xdt
           enddo
        enddo
     endif

     ! determine the processor grid in use
     call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
          dims, dummy_periods, dummy_coords, code)

     if (dims(1)==1) then
        do k=1,xsize(3)
           do i=1,xsize(1)
              dpdxy1(i,k)=px1(i,1,k)/xdt
              dpdzy1(i,k)=pz1(i,1,k)/xdt
           enddo
        enddo
        do k=1,xsize(3)
           do i=1,xsize(1)
              dpdxyn(i,k)=px1(i,xsize(2),k)/xdt
              dpdzyn(i,k)=pz1(i,xsize(2),k)/xdt
           enddo
        enddo
     else
        !find j=1 and j=ny
        if (xstart(2)==1) then
           do k=1,xsize(3)
              do i=1,xsize(1)
                 dpdxy1(i,k)=px1(i,1,k)/xdt
                 dpdzy1(i,k)=pz1(i,1,k)/xdt
              enddo
           enddo
        endif
        !      print *,nrank,xstart(2),ny-(nym/p_row)
        if (ny-(nym/dims(1))==xstart(2)) then
           do k=1,xsize(3)
              do i=1,xsize(1)
                 dpdxyn(i,k)=px1(i,xsize(2),k)/xdt
                 dpdzyn(i,k)=pz1(i,xsize(2),k)/xdt
              enddo
           enddo
        endif

     endif

     if (nrank==0) print *,'reconstruction pressure gradients done!'
  endif

  if (iresflg .eq. 1 ) then !Writing restart
     if (nrank==0)  print *,'Restart point',itime/icheckpoint,'saved successfully!'
     !if (nrank==0) print *,'Elapsed time (s)',real(trestart,4)
     !if (nrank==0) print *,'Aproximated writing speed (MB/s)',real(((s3df*16.)*1e-6)/trestart,4)
     if (nrank==0) print *,'If necesseary restart from:',itime+1
  end if

end subroutine restart
!
!*******************************************************************
subroutine stretching()
  !
  !*******************************************************************
  !
  USE decomp_2d
  !USE decomp_2d_poisson
  USE variables
  USE param
  USE var
  USE MPI

  implicit none

  real(mytype) :: yinf,den,xnum,xcx,den1,den2,den3,den4,xnum1,cst
  integer :: j

  yinf=-yly/two
  den=two*beta*yinf
  xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf)
  alpha=abs(xnum/den)
  xcx=one/beta/alpha
  if (alpha.ne.0.) then
     if (istret.eq.1) yp(1)=zero
     if (istret.eq.2) yp(1)=zero
     if (istret.eq.1) yeta(1)=zero
     if (istret.eq.2) yeta(1)=-half
     if (istret.eq.3) yp(1)=zero
     if (istret.eq.3) yeta(1)=-half
     do j=2,ny
        if (istret==1) yeta(j)=real(j-1,mytype)*(one/nym)
        if (istret==2) yeta(j)=real(j-1,mytype)*(one/nym)-half
        if (istret==3) yeta(j)=real(j-1,mytype)*(half/nym)-half
      den1=sqrt(alpha*beta+one)
        xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
        den=two*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
        den3=((sin(pi*yeta(j)))*(sin(pi*yeta(j)))/beta/pi)+alpha/pi
        den4=two*alpha*beta-cos(two*pi*yeta(j))+one
        xnum1=(atan(xnum*tan(pi*yeta(j))))*den4/den1/den3/den
        cst=sqrt(beta)*pi/(two*sqrt(alpha)*sqrt(alpha*beta+one))
        if (istret==1) then
           if (yeta(j).lt.half) yp(j)=xnum1-cst-yinf
           if (yeta(j).eq.half) yp(j)=zero-yinf
           if (yeta(j).gt.half) yp(j)=xnum1+cst-yinf
        endif
        if (istret==2) then
           if (yeta(j).lt.half) yp(j)=xnum1-cst+yly
           if (yeta(j).eq.half) yp(j)=zero+yly
           if (yeta(j).gt.half) yp(j)=xnum1+cst+yly
        endif
        if (istret==3) then
           if (yeta(j).lt.half) yp(j)=(xnum1-cst+yly)*two
           if (yeta(j).eq.half) yp(j)=(zero+yly)*two
           if (yeta(j).gt.half) yp(j)=(xnum1+cst+yly)*two
        endif
     enddo
  endif
  if (alpha.eq.0.) then
     yp(1)=-1.e10
     do j=2,ny
        yeta(j)=real(j-1,mytype)*(one/ny)
        yp(j)=-beta*cos(pi*yeta(j))/sin(yeta(j)*pi)
     enddo
  endif
  if (alpha.ne.0.) then
     do j=1,ny
        if (istret==1) yetai(j)=(real(j,mytype)-half)*(one/nym)
        if (istret==2) yetai(j)=(real(j,mytype)-half)*(one/nym)-half
        if (istret==3) yetai(j)=(real(j,mytype)-half)*(half/nym)-half
        den1=sqrt(alpha*beta+one)
        xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
        den=2.*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
        den3=((sin(pi*yetai(j)))*(sin(pi*yetai(j)))/beta/pi)+alpha/pi
        den4=two*alpha*beta-cos(two*pi*yetai(j))+one
        xnum1=(atan(xnum*tan(pi*yetai(j))))*den4/den1/den3/den
        cst=sqrt(beta)*pi/(two*sqrt(alpha)*sqrt(alpha*beta+one))
        if (istret==1) then
           if (yetai(j).lt.half) ypi(j)=xnum1-cst-yinf
           if (yetai(j).eq.half) ypi(j)=zero-yinf
           if (yetai(j).gt.half) ypi(j)=xnum1+cst-yinf
        endif
        if (istret==2) then
           if (yetai(j).lt.half) ypi(j)=xnum1-cst+yly
           if (yetai(j).eq.half) ypi(j)=zero+yly
           if (yetai(j).gt.half) ypi(j)=xnum1+cst+yly
        endif
        if (istret==3) then
           if (yetai(j).lt.half) ypi(j)=(xnum1-cst+yly)*two
           if (yetai(j).eq.half) ypi(j)=(zero+yly)*two
           if (yetai(j).gt.half) ypi(j)=(xnum1+cst+yly)*two
        endif
     enddo
  endif
  if (alpha.eq.0.) then
     ypi(1)=-1.e10
     do j=2,ny
        yetai(j)=real(j-1,mytype)*(one/ny)
        ypi(j)=-beta*cos(pi*yetai(j))/sin(yetai(j)*pi)
     enddo
  endif

!Mapping!!, metric terms
if (istret .ne. 3) then
  do j=1,ny
     ppy(j)=yly*(alpha/pi+(one/pi/beta)*sin(pi*yeta(j))*sin(pi*yeta(j)))
     pp2y(j)=ppy(j)*ppy(j)
     pp4y(j)=(-two/beta*cos(pi*yeta(j))*sin(pi*yeta(j)))
  enddo
  do j=1,ny
     ppyi(j)=yly*(alpha/pi+(one/pi/beta)*sin(pi*yetai(j))*sin(pi*yetai(j)))
     pp2yi(j)=ppyi(j)*ppyi(j)
     pp4yi(j)=(-two/beta*cos(pi*yetai(j))*sin(pi*yetai(j)))
  enddo
endif

if (istret .eq. 3) then
do j=1,ny
   ppy(j)=yly*(alpha/pi+(one/pi/beta)*sin(pi*yeta(j))*sin(pi*yeta(j)))
   pp2y(j)=ppy(j)*ppy(j)
   pp4y(j)=(-two/beta*cos(pi*yeta(j))*sin(pi*yeta(j)))/two
enddo
do j=1,ny
   ppyi(j)=yly*(alpha/pi+(one/pi/beta)*sin(pi*yetai(j))*sin(pi*yetai(j)))
   pp2yi(j)=ppyi(j)*ppyi(j)
   pp4yi(j)=(-two/beta*cos(pi*yetai(j))*sin(pi*yetai(j)))/two
enddo
endif

!   yp(1) = 0.0
!   yp(2) = 0.01
!   coeff0= 1.1
!   blender1 = 0.0
!   blender2 = 0.0
!   do j=3,ny
!!      yeta(j)=(j-1.)*(1./ny)
!!      yp(j)=-beta*cos(pi*yeta(j))/sin(yeta(j)*pi)
!
!     if (yp(j-1).LE.3.5*1.0) then
!       dy_plus_target = 8.0
!       !Calculate re_tau guess somewhere
!      dy_plus_current= (yp(j-1)-yp(j-2))*85.0
!       !dy_plus_coeff is from 1 to 0
!       dy_plus_coeff = (dy_plus_target-dy_plus_current)/dy_plus_target
!       coeff = coeff0**dy_plus_coeff
!
!       dy_plus_coeff_old1 = dy_plus_coeff   !will be required for blenders
!     else if (yp(j-1).GE.39.0*1.0) then
!       dy_plus_target = 10.0
!       !Calculate re_tau guess somewhere
!       dy_plus_current= (yp(j-1)-yp(j-2))*85.0
!       !dy_plus_coeff is from 1 to 0
!       dy_plus_coeff = (dy_plus_target-dy_plus_current)/dy_plus_target
!
!       if (blender2.LT.1.0) blender2 = blender2 + 0.1   !carry the coeff smoothly
!       coeff = coeff0**((1.0-blender2)*dy_plus_coeff_old2+blender2*dy_plus_coeff)
!     else
!       dy_plus_target = 80.0
!       !Calculate re_tau guess somewhere
!       dy_plus_current= (yp(j-1)-yp(j-2))*85.0
!       !dy_plus_coeff is from 1 to 0
!       dy_plus_coeff = (dy_plus_target-dy_plus_current)/dy_plus_target
!
!       if (blender1.LT.1.0) blender1 = blender1 + 0.1   !carry the coeff smoothly
!       coeff = coeff0**((1.0-blender1)*dy_plus_coeff_old1+blender1*dy_plus_coeff)
!
!       dy_plus_coeff_old2 = dy_plus_coeff   !will be required for blenders
!     endif
!     yp(j) = yp(j-1)+(yp(j-1)-yp(j-2))*coeff
!   enddo
!
!   !Normalize to yly
!   ypmax = yp(ny)
!   yp = yp/ypmax*yly

if (nrank==0) then
open(10,file='yp.dat', form='formatted')
do j=1,ny
write(10,*)yp(j)
enddo
close(10)
endif

end subroutine stretching

!*****************************************************************
!
subroutine inversion5_v1(aaa,eee,spI)
!
!*****************************************************************

USE decomp_2d
!USE decomp_2d_poisson
USE variables
USE param
USE var
USE MPI

implicit none

! decomposition object for spectral space
TYPE(DECOMP_INFO) :: spI

#ifdef DOUBLE_PREC
  real(mytype), parameter :: epsilon = 1.e-16
#else
  real(mytype), parameter :: epsilon = 1.e-8
#endif

complex(mytype),dimension(spI%yst(1):spI%yen(1),ny/2,spI%yst(3):spI%yen(3),5) :: aaa
complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(2):spI%yen(2),spI%yst(3):spI%yen(3)) :: eee
integer :: i,j,k,m,mi,jc
integer,dimension(2) :: ja,jb
complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: sr
complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: a1,b1

real(mytype) :: tmp1,tmp2,tmp3,tmp4

do i=1,2
   ja(i)=4-i
   jb(i)=5-i
enddo
do m=1,ny/2-2
   do i=1,2
      mi=m+i
      do k=spI%yst(3),spI%yen(3)
      do j=spI%yst(1),spI%yen(1)
         if (real(aaa(j,m,k,3), kind=mytype).ne.zero) tmp1=real(aaa(j,mi,k,3-i), kind=mytype)/real(aaa(j,m,k,3), kind=mytype)
         if (aimag(aaa(j,m,k,3)).ne.zero)tmp2=aimag(aaa(j,mi,k,3-i))/aimag(aaa(j,m,k,3))
         sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
         eee(j,mi,k)=cmplx(real(eee(j,mi,k), kind=mytype)-tmp1*real(eee(j,m,k), kind=mytype),&
              aimag(eee(j,mi,k))-tmp2*aimag(eee(j,m,k)), kind=mytype)
      enddo
      enddo
      do jc=ja(i),jb(i)
      do k=spI%yst(3),spI%yen(3)
      do j=spI%yst(1),spI%yen(1)
         aaa(j,mi,k,jc)=cmplx(real(aaa(j,mi,k,jc), kind=mytype)-real(sr(j,k), kind=mytype)*real(aaa(j,m,k,jc+i), kind=mytype),&
              aimag(aaa(j,mi,k,jc))-aimag(sr(j,k))*aimag(aaa(j,m,k,jc+i)), kind=mytype)
      enddo
      enddo
      enddo
   enddo
enddo


do k=spI%yst(3),spI%yen(3)
do j=spI%yst(1),spI%yen(1)
   if (abs(real(aaa(j,ny/2-1,k,3), kind=mytype)).gt.epsilon) then
      tmp1=real(aaa(j,ny/2,k,2), kind=mytype)/real(aaa(j,ny/2-1,k,3), kind=mytype)
   else
      tmp1=zero
   endif
   if (abs(aimag(aaa(j,ny/2-1,k,3))).gt.epsilon) then
      tmp2=aimag(aaa(j,ny/2,k,2))/aimag(aaa(j,ny/2-1,k,3))
   else
      tmp2=zero
   endif
   sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   b1(j,k)=cmplx(real(aaa(j,ny/2,k,3), kind=mytype)-tmp1*real(aaa(j,ny/2-1,k,4), kind=mytype),&
        aimag(aaa(j,ny/2,k,3))-tmp2*aimag(aaa(j,ny/2-1,k,4)), kind=mytype)

   if (abs(real(b1(j,k), kind=mytype)).gt.epsilon) then
      tmp1=real(sr(j,k), kind=mytype)/real(b1(j,k), kind=mytype)
      tmp3=real(eee(j,ny/2,k), kind=mytype)/real(b1(j,k), kind=mytype)-tmp1*real(eee(j,ny/2-1,k), kind=mytype)
   else
      tmp1=zero
      tmp3=zero
   endif
   if (abs(aimag(b1(j,k))).gt.epsilon) then
      tmp2=aimag(sr(j,k))/aimag(b1(j,k))
      tmp4=aimag(eee(j,ny/2,k))/aimag(b1(j,k))-tmp2*aimag(eee(j,ny/2-1,k))
   else
      tmp2=zero
      tmp4=zero
   endif
   a1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   eee(j,ny/2,k)=cmplx(tmp3,tmp4, kind=mytype)

   if (abs(real(aaa(j,ny/2-1,k,3), kind=mytype)).gt.epsilon) then
      tmp1=one/real(aaa(j,ny/2-1,k,3), kind=mytype)
   else
      tmp1=zero
   endif
   if (abs(aimag(aaa(j,ny/2-1,k,3))).gt.epsilon) then
      tmp2=one/aimag(aaa(j,ny/2-1,k,3))
   else
      tmp2=zero
   endif
   b1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   a1(j,k)=cmplx(real(aaa(j,ny/2-1,k,4), kind=mytype)*real(b1(j,k), kind=mytype),&
        aimag(aaa(j,ny/2-1,k,4))*aimag(b1(j,k)), kind=mytype)
   eee(j,ny/2-1,k)=cmplx(real(eee(j,ny/2-1,k))*real(b1(j,k))-real(a1(j,k))*real(eee(j,ny/2,k)),&
        aimag(eee(j,ny/2-1,k))*aimag(b1(j,k))-aimag(a1(j,k))*aimag(eee(j,ny/2,k)), kind=mytype)
enddo
enddo

do i=ny/2-2,1,-1
do k=spI%yst(3),spI%yen(3)
do j=spI%yst(1),spI%yen(1)
   if (abs(real(aaa(j,i,k,3), kind=mytype)).gt.epsilon) then
      tmp1=one/real(aaa(j,i,k,3), kind=mytype)
   else
      tmp1=zero
   endif
   if (abs(aimag(aaa(j,i,k,3))).gt.epsilon) then
      tmp2=one/aimag(aaa(j,i,k,3))
   else
      tmp2=zero
   endif
   sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   a1(j,k)=cmplx(real(aaa(j,i,k,4), kind=mytype)*real(sr(j,k), kind=mytype),&
        aimag(aaa(j,i,k,4))*aimag(sr(j,k)), kind=mytype)
   b1(j,k)=cmplx(real(aaa(j,i,k,5), kind=mytype)*real(sr(j,k), kind=mytype),&
        aimag(aaa(j,i,k,5))*aimag(sr(j,k)), kind=mytype)
   eee(j,i,k)=cmplx(real(eee(j,i,k), kind=mytype)*real(sr(j,k), kind=mytype)-&
        real(a1(j,k), kind=mytype)*real(eee(j,i+1,k), kind=mytype)-&
        real(b1(j,k), kind=mytype)*real(eee(j,i+2,k), kind=mytype),&
        aimag(eee(j,i,k))*aimag(sr(j,k))-&
        aimag(a1(j,k))*aimag(eee(j,i+1,k))-aimag(b1(j,k))*aimag(eee(j,i+2,k)), kind=mytype)
enddo
enddo
enddo

return

end subroutine inversion5_v1

!*****************************************************************
!
subroutine inversion5_v2(aaa,eee,spI)
!
!*****************************************************************

USE decomp_2d
!USE decomp_2d_poisson
USE variables
USE param
USE var
USE MPI

implicit none

! decomposition object for spectral space
TYPE(DECOMP_INFO) :: spI

#ifdef DOUBLE_PREC
  real(mytype), parameter :: epsilon = 1.e-16
#else
  real(mytype), parameter :: epsilon = 1.e-8
#endif

complex(mytype),dimension(spI%yst(1):spI%yen(1),nym,spI%yst(3):spI%yen(3),5) :: aaa
complex(mytype),dimension(spI%yst(1):spI%yen(1),nym,spI%yst(3):spI%yen(3)) :: eee
integer :: i,j,k,m,mi,jc
integer,dimension(2) :: ja,jb
complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: sr
complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: a1,b1

real(mytype) :: tmp1,tmp2,tmp3,tmp4

do i=1,2
   ja(i)=4-i
   jb(i)=5-i
enddo
do m=1,nym-2
   do i=1,2
   mi=m+i
   do k=spI%yst(3),spI%yen(3)
   do j=spI%yst(1),spI%yen(1)
      if (real(aaa(j,m,k,3), kind=mytype).ne.zero) tmp1=real(aaa(j,mi,k,3-i), kind=mytype)/real(aaa(j,m,k,3), kind=mytype)
      if (aimag(aaa(j,m,k,3)).ne.zero)tmp2=aimag(aaa(j,mi,k,3-i))/aimag(aaa(j,m,k,3))
      sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
      eee(j,mi,k)=cmplx(real(eee(j,mi,k), kind=mytype)-tmp1*real(eee(j,m,k), kind=mytype),&
           aimag(eee(j,mi,k))-tmp2*aimag(eee(j,m,k)), kind=mytype)
   enddo
   enddo
   do jc=ja(i),jb(i)
   do k=spI%yst(3),spI%yen(3)
   do j=spI%yst(1),spI%yen(1)
      aaa(j,mi,k,jc)=cmplx(real(aaa(j,mi,k,jc), kind=mytype)-real(sr(j,k), kind=mytype)*real(aaa(j,m,k,jc+i), kind=mytype),&
           aimag(aaa(j,mi,k,jc))-aimag(sr(j,k))*aimag(aaa(j,m,k,jc+i)), kind=mytype)
   enddo
   enddo
   enddo
enddo
enddo
do k=spI%yst(3),spI%yen(3)
do j=spI%yst(1),spI%yen(1)
   if (abs(real(aaa(j,nym-1,k,3), kind=mytype)).gt.epsilon) then
      tmp1=real(aaa(j,nym,k,2), kind=mytype)/real(aaa(j,nym-1,k,3), kind=mytype)
   else
      tmp1=zero
   endif
   if (abs(aimag(aaa(j,nym-1,k,3))).gt.epsilon) then
      tmp2=aimag(aaa(j,nym,k,2))/aimag(aaa(j,nym-1,k,3))
   else
      tmp2=zero
   endif
   sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   b1(j,k)=cmplx(real(aaa(j,nym,k,3), kind=mytype)-tmp1*real(aaa(j,nym-1,k,4), kind=mytype),&
        aimag(aaa(j,nym,k,3))-tmp2*aimag(aaa(j,nym-1,k,4)), kind=mytype)
   if (abs(real(b1(j,k), kind=mytype)).gt.epsilon) then
      tmp1=real(sr(j,k), kind=mytype)/real(b1(j,k), kind=mytype)
      tmp3=real(eee(j,nym,k), kind=mytype)/real(b1(j,k), kind=mytype)-tmp1*real(eee(j,nym-1,k), kind=mytype)
   else
      tmp1=zero
      tmp3=zero
   endif
   if (abs(aimag(b1(j,k))).gt.epsilon) then
      tmp2=aimag(sr(j,k))/aimag(b1(j,k))
      tmp4=aimag(eee(j,nym,k))/aimag(b1(j,k))-tmp2*aimag(eee(j,nym-1,k))
   else
      tmp2=zero
      tmp4=zero
   endif
   a1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   eee(j,nym,k)=cmplx(tmp3,tmp4, kind=mytype)

   if (abs(real(aaa(j,nym-1,k,3), kind=mytype)).gt.epsilon) then
      tmp1=one/real(aaa(j,nym-1,k,3), kind=mytype)
   else
      tmp1=zero
   endif
   if (abs(aimag(aaa(j,nym-1,k,3))).gt.epsilon) then
      tmp2=one/aimag(aaa(j,nym-1,k,3))
   else
      tmp2=zero
   endif
   b1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   a1(j,k)=cmplx(real(aaa(j,nym-1,k,4), kind=mytype)*real(b1(j,k), kind=mytype),&
        aimag(aaa(j,nym-1,k,4))*aimag(b1(j,k)), kind=mytype)
   eee(j,nym-1,k)=cmplx(real(eee(j,nym-1,k), kind=mytype)*real(b1(j,k), kind=mytype)-&
        real(a1(j,k), kind=mytype)*real(eee(j,nym,k), kind=mytype),&
        aimag(eee(j,nym-1,k))*aimag(b1(j,k))-aimag(a1(j,k))*aimag(eee(j,nym,k)), kind=mytype)
enddo
enddo

do i=nym-2,1,-1
do k=spI%yst(3),spI%yen(3)
do j=spI%yst(1),spI%yen(1)
   if (abs(real(aaa(j,i,k,3), kind=mytype)).gt.epsilon) then
      tmp1=one/real(aaa(j,i,k,3), kind=mytype)
   else
      tmp1=zero
   endif
   if (abs(aimag(aaa(j,i,k,3))).gt.epsilon) then
      tmp2=one/aimag(aaa(j,i,k,3))
   else
      tmp2=zero
   endif
   sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   a1(j,k)=cmplx(real(aaa(j,i,k,4), kind=mytype)*real(sr(j,k), kind=mytype),&
        aimag(aaa(j,i,k,4))*aimag(sr(j,k)), kind=mytype)
   b1(j,k)=cmplx(real(aaa(j,i,k,5), kind=mytype)*real(sr(j,k), kind=mytype),&
        aimag(aaa(j,i,k,5))*aimag(sr(j,k)), kind=mytype)
   eee(j,i,k)=cmplx(real(eee(j,i,k), kind=mytype)*real(sr(j,k), kind=mytype)-&
        real(a1(j,k), kind=mytype)*real(eee(j,i+1,k), kind=mytype)-&
        real(b1(j,k), kind=mytype)*real(eee(j,i+2,k), kind=mytype),&
        aimag(eee(j,i,k))*aimag(sr(j,k))-&
        aimag(a1(j,k))*aimag(eee(j,i+1,k))-aimag(b1(j,k))*aimag(eee(j,i+2,k)), kind=mytype)
enddo
enddo
enddo

return

end subroutine inversion5_v2

!********************************************************************
!
subroutine channel (ux,constant)
  !
  !********************************************************************

  USE decomp_2d
  USE decomp_2d_poisson
  USE variables
  USE param
  USE var
  USE MPI

  implicit none

  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux
  real(mytype) :: constant

  integer :: j,i,k,code
  real(mytype) :: can,ut3,ut,ut4

  ut3=zero
  do k=1,ysize(3)
     do i=1,ysize(1)
        ut=zero
        do j=1,ny-1
#ifdef STRETCHING
           ut=ut+(yp(j+1)-yp(j))*(ux(i,j+1,k)-half*(ux(i,j+1,k)-ux(i,j,k)))
#else
           ut=ut+dy*(ux(i,j+1,k)-half*(ux(i,j+1,k)-ux(i,j,k)))
#endif
        enddo
        ut=ut/yly
        ut3=ut3+ut
     enddo
  enddo
  ut3=ut3/(real(nx*nz,mytype))

  call MPI_ALLREDUCE(ut3,ut4,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

  can=-(constant-ut4)

  if (nrank==0) print *,nrank,'UT',ut4,can

  do k=1,ysize(3)
     do i=1,ysize(1)
        do j=2,ny-1
           ux(i,j,k)=ux(i,j,k)-can
        enddo
     enddo
  enddo

  return
end subroutine channel

subroutine read_geomcomplex()
!
USE complex_geometry
USE decomp_2d
USE MPI
!
implicit none
!
integer :: i,j,k
integer :: code
!
if(nrank.eq.0)then
   open(11,file='nobjx.dat'  ,form='formatted', status='old')
   do k=1,nz
   do j=1,ny
      read(11,*)nobjx(j,k)
   enddo
   enddo
   close(11)
endif
call MPI_BCAST(nobjx,ny*nz,MPI_INTEGER,0,MPI_COMM_WORLD,code)
if(nrank.eq.0)then
   open(12,file='nobjy.dat'  ,form='formatted', status='old')
   do k=1,nz
   do i=1,nx
      read(12,*)nobjy(i,k)
   enddo
   enddo
   close(12)
endif
call MPI_BCAST(nobjy,nx*nz,MPI_INTEGER,0,MPI_COMM_WORLD,code)
if(nrank.eq.0)then
   open(13,file='nobjz.dat'  ,form='formatted', status='old')
   do j=1,ny
   do i=1,nx
      read(13,*)nobjz(i,j)
   enddo
   enddo
   close(13)
endif
call MPI_BCAST(nobjz,nx*ny,MPI_INTEGER,0,MPI_COMM_WORLD,code)
if(nrank.eq.0)then
   open(21,file='nxifpif.dat',form='formatted', status='old')
   do k=1,nz
   do j=1,ny
   do i=0,nobjmax
      read(21,*)nxipif(i,j,k),nxfpif(i,j,k)
   enddo
   enddo
   enddo
   close(21)
endif
call MPI_BCAST(nxipif,ny*nz*(nobjmax+1),MPI_INTEGER,0,MPI_COMM_WORLD,code)
call MPI_BCAST(nxfpif,ny*nz*(nobjmax+1),MPI_INTEGER,0,MPI_COMM_WORLD,code)
if(nrank.eq.0)then
   open(22,file='nyifpif.dat',form='formatted', status='old')
   do k=1,nz
   do i=1,nx
   do j=0,nobjmax
      read(22,*)nyipif(j,i,k),nyfpif(j,i,k)
   enddo
   enddo
   enddo
   close(22)
endif
call MPI_BCAST(nyipif,nx*nz*(nobjmax+1),MPI_INTEGER,0,MPI_COMM_WORLD,code)
call MPI_BCAST(nyfpif,nx*nz*(nobjmax+1),MPI_INTEGER,0,MPI_COMM_WORLD,code)
if(nrank.eq.0)then
   open(23,file='nzifpif.dat',form='formatted', status='old')
   do j=1,ny
   do i=1,nx
   do k=0,nobjmax
      read(23,*)nzipif(k,i,j),nzfpif(k,i,j)
   enddo
   enddo
   enddo
   close(23)
endif
call MPI_BCAST(nzipif,nx*ny*(nobjmax+1),MPI_INTEGER,0,MPI_COMM_WORLD,code)
call MPI_BCAST(nzfpif,nx*ny*(nobjmax+1),MPI_INTEGER,0,MPI_COMM_WORLD,code)
if(nrank.eq.0)then
   open(31,file='xixf.dat'   ,form='formatted', status='old')
   do k=1,nz
   do j=1,ny
   do i=1,nobjmax
      read(31,*)xi(i,j,k),xf(i,j,k)
   enddo
   enddo
   enddo
   close(31)
endif
call MPI_BCAST(xi,ny*nz*nobjmax,MPI_REAL,0,MPI_COMM_WORLD,code)
call MPI_BCAST(xf,ny*nz*nobjmax,MPI_REAL,0,MPI_COMM_WORLD,code)
if(nrank.eq.0)then
   open(32,file='yiyf.dat'   ,form='formatted', status='old')
   do k=1,nz
   do i=1,nx
   do j=1,nobjmax
      read(32,*)yi(j,i,k),yf(j,i,k)
   enddo
   enddo
   enddo
   close(32)
endif
call MPI_BCAST(yi,nx*nz*nobjmax,MPI_REAL,0,MPI_COMM_WORLD,code)
call MPI_BCAST(yf,nx*nz*nobjmax,MPI_REAL,0,MPI_COMM_WORLD,code)
if(nrank.eq.0)then
   open(33,file='zizf.dat'   ,form='formatted', status='old')
   do j=1,ny
   do i=1,nx
   do k=1,nobjmax
      read(33,*)zi(k,i,j),zf(k,i,j)
   enddo
   enddo
   enddo
   close(33)
endif
call MPI_BCAST(zi,nx*ny*nobjmax,MPI_REAL,0,MPI_COMM_WORLD,code)
call MPI_BCAST(zf,nx*ny*nobjmax,MPI_REAL,0,MPI_COMM_WORLD,code)
!
return
end subroutine read_geomcomplex
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
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
integer                                            :: ix              != position du point "zappé"
integer                                            :: ipif,ipol,nxpif
integer                                            :: ipoli,ipolf     != positions Initiales et Finales du POLynôme considéré
real(mytype)                                       :: xpol,ypol,dypol !|variables concernant les polynômes
real(mytype),dimension(10)                         :: xa,ya           !|de Lagrange. A mettre impérativement en 
integer                                            :: ia,na           !|double précision
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
         ya(ia)=0.
         if(xi(i,j,k).gt.0.)then!objet immergé
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
         ya(ia)=0.
         if(xf(i,j,k).lt.xlx)then!objet immergé
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
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
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
integer                                            :: jy              != position du point "zappé"
integer                                            :: jpif,jpol,nypif
integer                                            :: jpoli,jpolf     != positions Initiales et Finales du POLynôme considéré
real(mytype)                                       :: xpol,ypol,dypol !|variables concernant les polynômes
real(mytype),dimension(10)                         :: xa,ya           !|de Lagrange. A mettre impérativement en 
integer                                            :: ia,na           !|double précision
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
         ya(ia)=0.
         if(yi(j,i,k).gt.0.)then!objet immergé
            jy=1!jy=yi(j,i,k)/dy+1
            do while(yp(jy).lt.yi(j,i,k))
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
         ya(ia)=0.
         if(yf(j,i,k).lt.yly)then!objet immergé
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
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
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
integer                                            :: kz              != position du point "zappé"
integer                                            :: kpif,kpol,nzpif
integer                                            :: kpoli,kpolf     != positions Initiales et Finales du POLynôme considéré
real(mytype)                                       :: xpol,ypol,dypol !|variables concernant les polynômes
real(mytype),dimension(10)                         :: xa,ya           !|de Lagrange. A mettre imérativement en 
integer                                            :: ia,na           !|double précision
!
do j=1,zsize(2)
do i=1,zsize(1)
   if(nobjz(i,j).ne.0)then
      ia=0
      do k=1,nobjz(i,j)          !boucle sur le nombre d'objets par couple (i,j)
         !1ère frontière
         nzpif=npif
         ia=ia+1
         xa(ia)=zi(k,i,j)
         ya(ia)=0.
         if(zi(k,i,j).gt.0.)then!objet immergé
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
         ya(ia)=0.
         if(zf(k,i,j).lt.zlz)then!objet immergé
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
!
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine polint(xa,ya,n,x,y,dy)
!
USE decomp_2d
!
implicit none
!
   integer,parameter            :: nmax=30
   integer                      :: n,i,m,ns
   real(mytype)                 :: dy,x,y,den,dif,dift,ho,hp,w
   real(mytype),dimension(nmax) :: c,d
   real(mytype),dimension(n)    :: xa,ya
   ns=1
   dif=abs(x-xa(1))
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
!********************************************************************
subroutine tripping(tb,ta) !TRIPPING SUBROUTINE FOR TURBULENT BOUNDARY LAYERS

USE param
USE variables
USE decomp_2d
USE MPI

implicit none

integer :: i,j,k
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tb, ta
integer :: seed0, ii, code
real(mytype) :: z_pos, randx, p_tr, b_tr, x_pos, y_pos, A_tr

!Done in X-Pencils
seed0=randomseed !Seed for random number
!A_tr=A_trip*min(1.0,0.8+real(itime)/200.0)
!xs_tr=4.0/2.853
!ys_tr=2.0/2.853
!ts_tr=4.0/2.853
!x0_tr=40.0/2.853
A_tr = 0.1*dt

if ((itime.eq.ifirst).and.(nrank.eq.0)) then
call random_seed(SIZE=ii)
call random_seed(PUT=seed0*(/ (1, i = 1, ii) /))

!DEBUG:
!call random_number(randx)
!call MPI_BCAST(randx,1,real_type,0,MPI_COMM_WORLD,code)
!write(*,*) 'RANDOM:', nrank, randx, ii
!First random generation of h_nxt


  do j=1,z_modes

    call random_number(randx)
    h_coeff(j)=1.0*(randx-0.5)
  enddo
    h_coeff=h_coeff/sqrt(DBLE(z_modes)) 
endif

!Initialization h_nxt  (always bounded by xsize(3)^2 operations)
if (itime.eq.ifirst) then
  call MPI_BCAST(h_coeff,z_modes,real_type,0,MPI_COMM_WORLD,code)
  nxt_itr=0
  do k=1,xsize(3)
     h_nxt(k)=0.0
     z_pos=-zlz/2.0+(xstart(3)+(k-1)-1)*dz
   do j=1,z_modes
      h_nxt(k)= h_nxt(k)+h_coeff(j)*sin(2.0*pi*j*z_pos/zlz)
   enddo
  enddo
end if



!Time-loop
 i=int(t/ts_tr)
    if (i.ge.nxt_itr) then  !Nxt_itr is a global variable
        nxt_itr=i+1
        
        !First random generation of h
        h_i(:)=h_nxt(:)
        if (nrank .eq. 0) then
         do j=1,z_modes
            call random_number(randx)
            h_coeff(j)=1.0*(randx-0.5)
         enddo
        h_coeff=h_coeff/sqrt(DBLE(z_modes)) !Non-dimensionalization
        end if
        
        call MPI_BCAST(h_coeff,z_modes,real_type,0,MPI_COMM_WORLD,code)
          
        
        !Initialization h_nxt  (always bounded by z_steps^2 operations)
       do k=1,xsize(3)
          h_nxt(k)=0.0
          z_pos=-zlz/2.0+(xstart(3)+(k-1)-1)*dz
          do j=1,z_modes
            h_nxt(k)= h_nxt(k)+h_coeff(j)*sin(2.0*pi*j*z_pos/zlz)
          enddo
        enddo
     endif


 !Time coefficient
  p_tr=t/ts_tr-i
  b_tr=3.0*p_tr**2-2.0*p_tr**3
  
  !Creation of tripping velocity
  do i=1,xsize(1)
    x_pos=(xstart(1)+(i-1)-1)*dx
    do j=1,xsize(2) 
      !y_pos=(xstart(2)+(j-1)-1)*dy   
      y_pos=yp(xstart(2)+(j-1))
      do k=1,xsize(3)
       !g(z)*EXP_F(X,Y)
       ta(i,j,k)=((1.0-b_tr)*h_i(k)+b_tr*h_nxt(k))
       !ta(i,j,k)=A_tr*exp(-((x_pos-x0_tr)/xs_tr)**2-(y_pos/ys_tr)**2)*ta(i,j,k)
       ta(i,j,k)=A_tr*exp(-((x_pos-x0_tr)/xs_tr)**2-((y_pos-0.5)/ys_tr)**2)*ta(i,j,k)  
       tb(i,j,k)=tb(i,j,k)+ta(i,j,k)
       
       z_pos=-zlz/2.0+(xstart(3)+(k-1)-1)*dz
      ! if ((((x_pos-x0_tr)**2).le.9.0e-3).and.(y_pos.le.0.0001).and.((z_pos).le.0.03))then
      !       open(442,file='tripping.dat',form='formatted',position='APPEND')
      !  write(442,*) t,ta(i,j,k)
      !  close(442)   
      ! end if

      enddo
    enddo
  enddo
    
return   
end subroutine tripping
!********************************************************************
function rl(complexnumber)

USE param

implicit none

real(mytype) :: rl
complex(mytype) :: complexnumber

rl = real(complexnumber, kind=mytype)

end function rl
!********************************************************************
function iy(complexnumber)

USE param

implicit none

real(mytype) :: iy
complex(mytype) :: complexnumber

iy = aimag(complexnumber)

end function iy
!********************************************************************
function cx(realpart,imaginarypart)

USE param

implicit none

complex(mytype) :: cx
real(mytype) :: realpart, imaginarypart

cx = cmplx(realpart, imaginarypart, kind=mytype)

end function cx
!********************************************************************

