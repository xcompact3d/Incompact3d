! ***********************************************************
!
subroutine waves ()
!
!***********************************************************

   USE derivX 	
   USE derivY 	
   USE derivZ 
   USE decomp_2d
   USE decomp_2d_fft
   USE variables
   USE param	
  
  implicit none
  
  integer :: i,j,k	
!  complex(mytype), dimension(nz/2+1) :: zkz,zk2,ezs
!  complex(mytype), dimension(ny) :: yky,yk2,eys	
!  complex(mytype), dimension(nx) :: xkx,xk2,exs	
  real(mytype) :: w, w1	

   if (nclx==0) then  
!WAVE NUMBER IN X
  do i=1,nx/2+1
     w=2.*pi*(i-1)/nx
     xkx(i)=cmplx(nx*w/xlx,nx*w/yly)
     exs(i)=cmplx(nx*w/xlx,nx*w/yly)
     xk2(i)=cmplx((nx*w/xlx)**2,(nx*w/xlx)**2)
  enddo
  do i=nx/2+2,nx
     xkx(i)=xkx(nx-i+2)
     exs(i)=exs(nx-i+2)
     xk2(i)=xk2(nx-i+2)
  enddo
else
 do i=1,nx
     w=2.*pi*(i-1)*0.5/nx
     w1=2.*pi*(nx-i+1)*0.5/nx
     xkx(i)=cmplx(nx*w/xlx,nx*w/xlx)
     exs(i)=cmplx(nx*w/xlx,nx*w/xlx)
     xk2(i)=cmplx((nx*w/xlx)**2,(nx*w/xlx)**2)
  enddo
endif

!WAVE NUMBER IN Y
  do j=1,ny/2+1
     w=2.*pi*(j-1)/ny
     yky(j)=cmplx(ny*w/yly,ny*w/yly)
     eys(j)=cmplx(ny*w/yly,ny*w/yly)
     yk2(j)=cmplx((ny*w/yly)**2,(ny*w/yly)**2)
  enddo
  do j=ny/2+2,ny
     yky(j)=yky(ny-j+2)
     eys(j)=eys(ny-j+2)
     yk2(j)=yk2(ny-j+2)
  enddo

!WAVE NUMBER IN Z
  do k=1,nz/2+1
     w=2.*pi*(k-1)/nz
     zkz(k)=cmplx(nz*w/zlz,nz*w/zlz)
     ezs(k)=cmplx(nz*w/zlz,nz*w/zlz)
     zk2(k)=cmplx((nz*w/zlz)**2,(nz*w/zlz)**2)
  enddo

return
end subroutine waves
