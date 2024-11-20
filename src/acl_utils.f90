!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module actuator_line_model_utils

    use decomp_2d_constants, only: mytype
    use param, only: zero, one, two
    use decomp_2d_mpi, only : decomp_2d_abort
    
    implicit none

    interface trilinear_interpolation
      module procedure trilinear_interpolation_scalar
      module procedure trilinear_interpolation_vec3
    end interface trilinear_interpolation 

    public QuatRot, cross, IsoKernel, AnIsoKernel, int2str, trilinear_interpolation

contains

    !*******************************************************************************
    !
    real(mytype) function trilinear_interpolation_scalar(x0,y0,z0, &
                                                         x1,y1,z1, &
                                                         x,y,z,    &
                                                         u000,u100,u001,u101, &
                                                         u010,u110,u011,u111)
    !
    !*******************************************************************************
      use, intrinsic :: ieee_arithmetic
      
      implicit none
      real(mytype),intent(in) :: x0,y0,z0,x1,y1,z1,x,y,z,u000,u100,u001,u101,u010,u110,u011,u111
      real(mytype) :: c00,c01,c10,c11,c0,c1,xd,yd,zd

      if (x1/=x0) then
         xd=(x-x0)/(x1-x0)
      else
         xd=zero
      endif

      if (y1/=y0) then
         yd=(y-y0)/(y1-y0)
      else
         yd=zero
      endif

      if (z1/=z0) then
         zd=(z-z0)/(z1-z0)
      else
         zd=zero
      endif

      ! Interpolate along X
      c00=u000*(one-xd)+u100*xd
      c01=u001*(one-xd)+u101*xd
      c10=u010*(one-xd)+u110*xd
      c11=u011*(one-xd)+u111*xd

      ! Interpolate along Y
      c0 = c00*(one-yd)+c10*yd
      c1 = c01*(one-yd)+c11*yd

      ! Interpolate along Z
      trilinear_interpolation_scalar=c0*(one-zd)+c1*zd

      if(ieee_is_nan(trilinear_interpolation_scalar)) then
        print*,'        x0:',x0,y0,z0
        print*,'        x1:',x1,y1,z1
        print*,'         x:',x,y,z
        print*,'         u:',u000,u100,u001,u101,u010,u110,u011,u111

        call decomp_2d_abort(1,"error in trilinear_interpolation_scalar")

      endif

      return

    end function trilinear_interpolation_scalar

    function trilinear_interpolation_vec3(x0,y0,z0,x1,y1,z1,x,y,z, &
                                          u000,u100,u001,u101,u010,u110,u011,u111) result(u)
    !
    !*******************************************************************************
      use, intrinsic :: ieee_arithmetic
      
      implicit none
      real(mytype),intent(in) :: x0,y0,z0,x1,y1,z1,x,y,z
      real(mytype),intent(in) :: u000(3),u100(3),u001(3),u101(3),u010(3),u110(3),u011(3),u111(3)
      real(mytype) :: u(3)
      real(mytype) :: c00(3),c01(3),c10(3),c11(3),c0(3),c1(3),xd,yd,zd

      if (x1 .ne. x0) then
         xd=(x-x0)/(x1-x0)
      else
         xd=0._mytype
      endif

      if (y1 .ne. y0) then
         yd=(y-y0)/(y1-y0)
      else
         yd=0._mytype
      endif

      if (z1 .ne. z0) then
         zd=(z-z0)/(z1-z0)
      else
         zd=0._mytype
      endif

      ! Interpolate along X
      c00=u000*(1._mytype-xd)+u100*xd
      c01=u001*(1._mytype-xd)+u101*xd
      c10=u010*(1._mytype-xd)+u110*xd
      c11=u011*(1._mytype-xd)+u111*xd

      ! Interpolate along Y
      c0 = c00*(1._mytype-yd)+c10*yd
      c1 = c01*(1._mytype-yd)+c11*yd

      ! Interpolate along Z
      u=c0*(1._mytype-zd)+c1*zd

      if(ieee_is_nan(u(1)) .or. ieee_is_nan(u(2)) .or. ieee_is_nan(u(3))) then
        print*,'        x0:',x0,y0,z0
        print*,'        x1:',x1,y1,z1
        print*,'         x:',x,y,z
        print*,'         u:',u000,u100,u001,u101,u010,u110,u011,u111

        call decomp_2d_abort(1,"error in trilinear_interpolation_vec3")

      endif

      return

    end function trilinear_interpolation_vec3

    !*******************************************************************************
    !
    subroutine cross(ax,ay,az,bx,by,bz,cx,cy,cz)
    !
    !*******************************************************************************

      real(mytype) :: ax,ay,az,bx,by,bz,cx,cy,cz

      cx = ay*bz - az*by
      cy = az*bx - ax*bz
      cz = ax*by - ay*bx

    end subroutine cross

    !*******************************************************************************
    !
    subroutine QuatRot(vx,vy,vz,Theta,Rx,Ry,Rz,Ox,Oy,Oz,vRx,vRy,vRz)
    !
    !*******************************************************************************

      implicit none
      real(mytype), intent(in) :: vx,vy,vz,Theta,Rx,Ry,Rz,Ox,Oy,Oz
      real(mytype),intent(inout) :: vRx,vRy,vRz
      real(mytype) :: nRx,nRy,nRz
      real(mytype) :: p(4,1), pR(4,1), q(4), qbar(4), RMag, vOx, vOy, vOz
      real(mytype) :: QL(4,4), QbarR(4,4)

      ! Perform rotation of vector v around normal vector nR using the quaternion machinery.
      ! v: input vector
      ! Theta: rotation angle (rad)
      ! nR: normal vector around which to rotate
      ! Origin: origin point of rotation
      ! vR: Rotated vector

      ! Force normalize nR
      RMag=sqrt(Rx**2.0+Ry**2.0+Rz**2.0)
      nRx=Rx/RMag
      nRy=Ry/RMag
      nRz=Rz/RMag

      ! Quaternion form of v
      vOx=vx-Ox
      vOy=vy-Oy
      vOz=vz-Oz
      p=reshape([zero,vOx,vOy,vOz],[4,1])

      ! Rotation quaternion and conjugate
      q=(/cos(Theta/2),nRx*sin(Theta/2),nRy*sin(Theta/2),nRz*sin(Theta/2)/)
      qbar=(/q(1),-q(2),-q(3),-q(4)/)

      QL=transpose(reshape((/q(1), -q(2), -q(3), -q(4), &
                             q(2),  q(1), -q(4),  q(3), &
                             q(3),  q(4),  q(1), -q(2), &
                             q(4), -q(3),  q(2),  q(1)/),(/4,4/)))

      QbarR=transpose(reshape((/qbar(1), -qbar(2), -qbar(3), -qbar(4), &
                                qbar(2),  qbar(1),  qbar(4), -qbar(3), &
                                qbar(3), -qbar(4),  qbar(1),  qbar(2), &
                                qbar(4),  qbar(3), -qbar(2),  qbar(1)/),(/4,4/)))

      ! Rotate p
      pR=matmul(matmul(QbarR,QL),p)
      vRx=pR(2,1)+Ox
      vRy=pR(3,1)+Oy
      vRz=pR(4,1)+Oz

    end subroutine QuatRot

    !*******************************************************************************
    !
    subroutine IDW(Ncol,Xcol,Ycol,Zcol,Fxcol,Fycol,Fzcol,p,Xmesh,Ymesh,Zmesh,Fxmesh,Fymesh,Fzmesh)
    !
    !*******************************************************************************

      implicit none
      integer, intent(in) :: Ncol
      real(mytype), dimension(Ncol), intent(in) :: Xcol,Ycol,Zcol,Fxcol,Fycol,Fzcol
      real(mytype), intent(in) :: Xmesh,Ymesh,Zmesh
      integer,intent(in) :: p
      real(mytype), intent(inout) :: Fxmesh,Fymesh,Fzmesh
      real(mytype), dimension(Ncol) :: d(Ncol), w(Ncol)
      real(mytype) ::  wsum
      integer :: i,imin

      wsum=zero
      do i=1,Ncol
         d(i)=sqrt((Xcol(i)-Xmesh)**2+(Ycol(i)-Ymesh)**2+(Zcol(i)-Zmesh)**2)
         w(i)=one/d(i)**p
         wsum=wsum+w(i)
      enddo

      if (minval(d)<0.001_mytype) then
         imin=minloc(d,1)
         Fxmesh=Fxcol(imin)
         Fymesh=Fycol(imin)
         Fzmesh=Fzcol(imin)
      else
         Fxmesh=zero
         Fymesh=zero
         Fzmesh=zero
         do i=1,Ncol
            Fxmesh=Fxmesh+w(i)*Fxcol(i)/wsum
            Fymesh=Fymesh+w(i)*Fycol(i)/wsum
            Fzmesh=Fzmesh+w(i)*Fzcol(i)/wsum
         enddo
      endif

    end subroutine IDW

    !*******************************************************************************
    !
    real(mytype) function IsoKernel(dr,epsilon_par,dim)
    !
    !*******************************************************************************

      use constants

      implicit none
      integer, intent(in) :: dim
      real(mytype), intent(in) :: dr, epsilon_par

      if (dim==2) then
         IsoKernel = one/(epsilon_par**2*pi)*exp(-(dr/epsilon_par)**2.0)
      else if (dim==3) then
         IsoKernel = one/(epsilon_par**3.0*pi**1.5)*exp(-(dr/epsilon_par)**2.0)
      else
         write(*,*) "1D source not implemented"
         stop
      endif

    end function IsoKernel

    !*******************************************************************************
    !
    real(mytype) function AnIsoKernel(dx,dy,dz,nx,ny,nz,tx,ty,tz,sx,sy,sz,ec,et,es)
    !
    !*******************************************************************************

      use constants

      implicit none
      real(mytype), intent(in) :: dx,dy,dz,nx,ny,nz,tx,ty,tz,sx,sy,sz,ec,et,es
      real(mytype) :: n,t,s

      n=dx*nx+dy*ny+dz*nz ! normal projection
      t=dx*tx+dy*ty+dz*tz ! Chordwise projection
      s=dx*sx+dy*sy+dz*sz ! Spanwise projection

      if (abs(s)<=es) then
         AnIsoKernel = exp(-((n/et)**2.0+(t/ec)**2.0))/(ec*et*pi)
      else
         AnIsoKernel = zero
      endif

    end function AnIsoKernel

    !*******************************************************************************
    !
    integer function FindMinimum(x,Start,End)
    !
    !*******************************************************************************

      implicit none
      integer, dimension(1:), intent(in) :: x
      integer, intent(in) :: Start, End
      integer :: Minimum
      integer :: Location
      integer :: i

      minimum = x(start)
      Location = Start
      do i=start+1,End
         if (x(i) < Minimum) then
            Minimum = x(i)
            Location = i
         endif
      enddo
      FindMinimum = Location

    end function FindMinimum

    !*******************************************************************************
    !
    subroutine swap(a,b)
    !
    !*******************************************************************************

      implicit none
      integer, intent(inout) :: a,b
      integer :: Temp

      Temp = a
      a = b
      b = Temp

    end subroutine swap

    !*******************************************************************************
    !
    subroutine sort(x,size)
    !
    !*******************************************************************************

      implicit none
      integer, dimension(1:), intent(inout) :: x
      integer, intent(in) :: size
      integer :: i
      integer :: Location

      do i=1,Size-1
         location=FindMinimum(x,i,size)
         call swap(x(i),x(Location))
      enddo

    end subroutine sort

    !*******************************************************************************
    !
    function dirname(number)
    !
    !*******************************************************************************

      integer, intent(in) :: number
      character(len=6) :: dirname

      ! Cast the (rounded) number to string using 6 digits and leading zeros
      write (dirname, '(I6.1)')  number
      ! This is the same w/o leading zeros
      !write (dirname, '(I6)')  nint(number)
      ! This is for one digit (no rounding)
      !write (dirname, '(F4.1)')  number

    end function dirname

    !*******************************************************************************
    !
    function outdirname(number)
    !
    !*******************************************************************************

      integer, intent(in) :: number
      character(len=6) :: outdirname

      ! Cast the (rounded) number to string using 6 digits and leading zeros
      write (outdirname, '(I6.1)')  number
      ! This is the same w/o leading zeros
      !write (dirname, '(I6)')  nint(number)
      ! This is for one digit (no rounding)
      !write (dirname, '(F4.1)')  number
 
    end function outdirname

    !*******************************************************************************
    !
    character(20) function int2str(num)
    !
    !*******************************************************************************

      integer, intent(in) ::num
      character(20) :: str

      ! convert integer to string using formatted write
      write(str, '(i20)') num
      int2str = adjustl(str)

    end function int2str

end module actuator_line_model_utils
