!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module ellipsoid_utils

    use decomp_2d, only: mytype
    use param, only: zero, one, two
    use dbg_schemes, only: sqrt_prec, cos_prec, exp_prec, sin_prec
    
    implicit none
    ! public QuatRot, cross, IsoKernel, AnIsoKernel, int2str

contains

    ! !*******************************************************************************
    ! !
    ! real(mytype) function trilinear_interpolation(x0,y0,z0, &
    !                                               x1,y1,z1, &
    !                                               x,y,z, &
    !                                               u000,u100,u001,u101, &
    !                                               u010,u110,u011,u111)
    ! !
    ! !*******************************************************************************

    !   implicit none
    !   real(mytype),intent(in) :: x0,y0,z0,x1,y1,z1,x,y,z,u000,u100,u001,u101,u010,u110,u011,u111
    !   real(mytype) :: c00,c01,c10,c11,c0,c1,xd,yd,zd

    !   if (x1/=x0) then
    !      xd=(x-x0)/(x1-x0)
    !   else
    !      xd=zero
    !   endif

    !   if (y1/=y0) then
    !      yd=(y-y0)/(y1-y0)
    !   else
    !      yd=zero
    !   endif

    !   if (z1/=z0) then
    !      zd=(z-z0)/(z1-z0)
    !   else
    !      zd=zero
    !   endif

    !   ! Interpolate along X
    !   c00=u000*(one-xd)+u100*xd
    !   c01=u001*(one-xd)+u101*xd
    !   c10=u010*(one-xd)+u110*xd
    !   c11=u011*(one-xd)+u111*xd

    !   ! Interpolate along Y
    !   c0 = c00*(one-yd)+c10*yd
    !   c1 = c01*(one-yd)+c11*yd

    !   ! Interpolate along Z
    !   trilinear_interpolation=c0*(one-zd)+c1*zd

    !   return

    ! end function trilinear_interpolation

    ! !*******************************************************************************
    ! !
    ! subroutine cross(ax,ay,az,bx,by,bz,cx,cy,cz)
    ! !
    ! !*******************************************************************************

    !   real(mytype) :: ax,ay,az,bx,by,bz,cx,cy,cz

    !   cx = ay*bz - az*by
    !   cy = az*bx - ax*bz
    !   cz = ax*by - ay*bx

    ! end subroutine cross

    ! !*******************************************************************************
    ! !
    ! subroutine QuatRot(vx,vy,vz,Theta,Rx,Ry,Rz,Ox,Oy,Oz,vRx,vRy,vRz)
    ! !
    ! !*******************************************************************************

    !   implicit none
    !   real(mytype), intent(in) :: vx,vy,vz,Theta,Rx,Ry,Rz,Ox,Oy,Oz
    !   real(mytype),intent(inout) :: vRx,vRy,vRz
    !   real(mytype) :: nRx,nRy,nRz
    !   real(mytype) :: p(4,1), pR(4,1), q(4), qbar(4), RMag, vOx, vOy, vOz
    !   real(mytype) :: QL(4,4), QbarR(4,4)

    !   ! Perform rotation of vector v around normal vector nR using the quaternion machinery.
    !   ! v: input vector
    !   ! Theta: rotation angle (rad)
    !   ! nR: normal vector around which to rotate
    !   ! Origin: origin point of rotation
    !   ! vR: Rotated vector

    !   ! Force normalize nR
    !   RMag=sqrt_prec(Rx**2.0+Ry**2.0+Rz**2.0)
    !   nRx=Rx/RMag
    !   nRy=Ry/RMag
    !   nRz=Rz/RMag

    !   ! Quaternion form of v
    !   vOx=vx-Ox
    !   vOy=vy-Oy
    !   vOz=vz-Oz
    !   p=reshape([zero,vOx,vOy,vOz],[4,1])

    !   ! Rotation quaternion and conjugate
    !   q=(/cos_prec(Theta/2),nRx*sin_prec(Theta/2),nRy*sin_prec(Theta/2),nRz*sin_prec(Theta/2)/)
    !   qbar=(/q(1),-q(2),-q(3),-q(4)/)

    !   QL=transpose(reshape((/q(1), -q(2), -q(3), -q(4), &
    !                          q(2),  q(1), -q(4),  q(3), &
    !                          q(3),  q(4),  q(1), -q(2), &
    !                          q(4), -q(3),  q(2),  q(1)/),(/4,4/)))

    !   QbarR=transpose(reshape((/qbar(1), -qbar(2), -qbar(3), -qbar(4), &
    !                             qbar(2),  qbar(1),  qbar(4), -qbar(3), &
    !                             qbar(3), -qbar(4),  qbar(1),  qbar(2), &
    !                             qbar(4),  qbar(3), -qbar(2),  qbar(1)/),(/4,4/)))

    !   ! Rotate p
    !   pR=matmul(matmul(QbarR,QL),p)
    !   vRx=pR(2,1)+Ox
    !   vRy=pR(3,1)+Oy
    !   vRz=pR(4,1)+Oz

    ! end subroutine QuatRot



    subroutine NormalizeQuaternion(quaternion) 
      real(mytype), intent(inout) :: quaternion(4)
      real(mytype) :: normalizedQuaternion(4)
    
      ! Compute the magnitude of the quaternion
      real(mytype) :: magnitude
      magnitude = sqrt(quaternion(1)**2 + quaternion(2)**2 + quaternion(3)**2 + quaternion(4)**2)
      if (magnitude < 0.0001) then
        magnitude = one
        write(*,*) "Tried to normalize a zero quaternion"
      endif
      ! Normalize the quaternion
      quaternion = quaternion / magnitude
    
    end subroutine NormalizeQuaternion

    subroutine QuaternionNorm(q,norm)
      real(mytype),intent(in) :: q(4)
      real(mytype),intent(out):: norm

      norm = sqrt_prec(q(1)**2+q(2)**2+q(3)**2+q(4)**2)

    end subroutine QuaternionNorm

    subroutine QuaternionConjugate(q, q_c)
      real(mytype), intent(in) :: q(4)
      real(mytype), intent(out) :: q_c(4)

      q_c = [q(1), -q(2), -q(3), -q(4)]
    end subroutine


    subroutine QuaternionMultiply(q1, q2, result)
      real(mytype), intent(in) :: q1(4), q2(4)
      real(mytype), intent(out) :: result(4)
      
      result(1) = q1(1) * q2(1) - q1(2) * q2(2) - q1(3) * q2(3) - q1(4) * q2(4)
      result(2) = q1(1) * q2(2) + q1(2) * q2(1) + q1(3) * q2(4) - q1(4) * q2(3)
      result(3) = q1(1) * q2(3) - q1(2) * q2(4) + q1(3) * q2(1) + q1(4) * q2(2)
      result(4) = q1(1) * q2(4) + q1(2) * q2(3) - q1(3) * q2(2) + q1(4) * q2(1)
    end subroutine QuaternionMultiply


    subroutine RotatePoint(point, quaternion, rotatedPoint)
      real(mytype), intent(in) :: point(3), quaternion(4)
      real(mytype), intent(out) :: rotatedPoint(3)
      real(mytype) :: conjugateQuaternion(4)
      real(mytype) :: resultQuaternion(4)
      real(mytype) :: rotatedPointQuaternion(4)

      ! Convert the point to a quaternion
      real(mytype) :: pointQuaternion(4)
      pointQuaternion(1) = 0.0D0
      pointQuaternion(2:4) = point(:)
    
      ! Perform the rotation
      
      conjugateQuaternion = [quaternion(1), -quaternion(2), -quaternion(3), -quaternion(4)]
    
      call QuaternionMultiply(quaternion, pointQuaternion, resultQuaternion)
      call QuaternionMultiply(resultQuaternion, conjugateQuaternion, rotatedPointQuaternion)
    
      ! Convert the rotated quaternion back to a 3D point
      rotatedPoint = rotatedPointQuaternion(2:4)
    end subroutine RotatePoint

    subroutine EllipsoidalRadius(point, centre, orientation, shape, radius)
      real(mytype), intent(in) :: point(3), centre(3), orientation(4), shape(3)
      real(mytype), intent(out) :: radius
      real(mytype) :: trans_point(3),rotated_point(3),scaled_point(3), orientation_c(4)
      integer :: i

      !translate point to body frame
      trans_point = point-centre

      call QuaternionConjugate(orientation, orientation_c)

      !rotate point into body frame (using inverse(conjugate) of orientation)
      call RotatePoint(trans_point, orientation, rotated_point)

      do i = 1,3
         scaled_point(i)=rotated_point(i)/shape(i)
      end do 

      radius=sqrt_prec(scaled_point(1)**2+scaled_point(2)**2+scaled_point(3)**2)

      if (radius /= radius) then
        write(*,*) "Got an error in grid check!"
        write(*,*) "Radius = ", radius, "point = ", point
        write(*,*) "Rotated point = ", rotated_point
        write(*,*) "Scaled Point = ", scaled_point
      endif

   end subroutine    

  subroutine CrossProduct(a, b, result)
    real(mytype), intent(in) :: a(3), b(3)
    real(mytype), intent(inout) :: result(3)
  
    result(1) = a(2) * b(3) - a(3) * b(2)
    result(2) = a(3) * b(1) - a(1) * b(3)
    result(3) = a(1) * b(2) - a(2) * b(1)
  end subroutine CrossProduct
  
  subroutine CalculatePointVelocity(point, center, angularVelocity, linearVelocity, pointVelocity)
    real(mytype), intent(in) :: point(3), center(3), linearVelocity(3), angularVelocity(3)
    real(mytype), intent(out) :: pointVelocity(3)
    real(mytype) :: crossed(3)
    ! Compute the distance vector from the center to the point
    real(mytype) :: distance(3)
    distance = point - center
  
    ! Compute the cross product of angular velocity and distance vector
    
    call CrossProduct(angularVelocity, distance, crossed)
  
    ! Calculate the velocity at the point
    pointVelocity = crossed + linearVelocity
  end subroutine CalculatePointVelocity

  subroutine is_inside_ellipsoid(point, centre, orientation, shape, ra, zeromach, is_inside)
    real(mytype), intent(in) :: point(3), centre(3), orientation(4), shape(3), ra, zeromach
    logical,intent(out)      :: is_inside
    real(mytype)             :: r

    call EllipsoidalRadius(point,centre,orientation,shape,r)
    
    is_inside = ((r-ra).lt.zeromach)

  end subroutine is_inside_ellipsoid


  subroutine navierFieldGen(center, linearVelocity, angularVelocity, ep1, ep1_x, ep1_y, ep1_z)
    use param
    use decomp_2d
    real(mytype), intent(in) :: center(3), linearVelocity(3), angularVelocity(3)
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(out) :: ep1_x, ep1_y, ep1_z
    real(mytype) :: xm, ym, zm, point(3), x_pv, y_pv, z_pv, pointVelocity(3)
    integer :: i,j,k

    do k = 1,xsize(3)
      zm=real(k+xstart(3)-2, mytype)*dz
      do j = 1,xsize(2)
        ym=real(j+xstart(2)-2, mytype)*dy
        do i = 1,xsize(1)
          xm=real(i+xstart(1)-2, mytype)*dx
          point=[xm,ym,zm]
          if (ep1(i,j,k).eq.1) then 
            call CalculatePointVelocity(point, center, angularVelocity, linearVelocity, pointVelocity)
            x_pv=pointVelocity(1)
            y_pv=pointVelocity(2)
            z_pv=pointVelocity(3)
          else
            x_pv=0
            y_pv=0
            z_pv=0
          endif
          ep1_x(i,j,k)=x_pv
          ep1_y(i,j,k)=y_pv
          ep1_z(i,j,k)=z_pv
        end do
      end do
    end do

    end subroutine navierFieldGen

  
    subroutine body_to_lab(p_body, q, p_lab)
      real(mytype),intent(in) :: p_body(4), q(4)
      real(mytype),intent(out):: p_lab(4)
      real(mytype)            :: q_inv(4),q_m(4)

      call QuaternionConjugate(q, q_inv)

      call QuaternionMultiply(p_body,q_inv,q_m)
      call QuaternionMultiply(q,q_m, p_lab)

    end subroutine body_to_lab

    subroutine lab_to_body(p_lab, q, p_body)
      real(mytype),intent(in) :: p_lab(4), q(4)
      real(mytype),intent(out):: p_body(4)
      real(mytype)            :: q_inv(4),q_m(4)

      call QuaternionConjugate(q, q_inv)

      call QuaternionMultiply(p_lab,q,q_m)
      call QuaternionMultiply(q_inv,q_m, p_body)

    end subroutine lab_to_body

    subroutine omega_stepper(omega_n, ang_accel, time_step, omega_n1)
      real(mytype),intent(in) :: omega_n(4), ang_accel(4), time_step
      real(mytype),intent(out):: omega_n1(4)
      
      omega_n1 = omega_n + ang_accel * time_step

    end subroutine omega_stepper

    subroutine orientation_stepper(q1, omega_q, time_step, q_n1)
      use param
      real(mytype),intent(in) :: q1(4),omega_q(4),time_step
      real(mytype),intent(out):: q_n1(4)
      real(mytype)            :: mag, re_part, im_sc, im_part(3), omega_n1(4)

      call QuaternionNorm(omega_q, mag)
      re_part=cos_prec(mag*time_step*half)
      if (mag.gt.zero) then
        im_sc = sin_prec(mag*time_step*half)/mag
      else 
        im_sc = zero
      endif
      im_part=im_sc * omega_q(2:4)
      omega_n1=[re_part,im_part(1),im_part(2),im_part(3)]
      
      call QuaternionMultiply(omega_n1,q1,q_n1)

    end subroutine orientation_stepper

    SUBROUTINE ConvertToMovingRotatingFrame(vI, positionI, originO, vO, Omega, vR)
      IMPLICIT NONE
    
      ! Arguments:
      ! vI       : Velocity in the inertial frame (3-element array)
      ! positionI: Position in the inertial frame (3-element array)
      ! originO  : Position of the origin of the rotating frame in the inertial frame (3-element array)
      ! vO       : Linear velocity of the origin of the rotating frame
      ! Omega    : Angular velocity of the rotating frame (3-element array)
      ! vR       : Velocity in the moving and rotating frame (Output, 3-element array)
    
      real(mytype), INTENT(IN)  :: vI(3), positionI(3), originO(3), vO(3), Omega(3)
      real(mytype), INTENT(OUT) :: vR(3)
      real(mytype) :: r(3), crossProduct_v(3)
    
      ! Compute r = positionI - originO
      r(1) = positionI(1) - originO(1)
      r(2) = positionI(2) - originO(2)
      r(3) = positionI(3) - originO(3)
    
      ! Compute Omega x r (cross product)
      crossProduct_v(1) = Omega(2)*r(3) - Omega(3)*r(2)
      crossProduct_v(2) = Omega(3)*r(1) - Omega(1)*r(3)
      crossProduct_v(3) = Omega(1)*r(2) - Omega(2)*r(1)
    
      ! Compute vR = vI - vO - Omega x r
      vR(1) = vI(1) - vO(1) - crossProduct_v(1)
      vR(2) = vI(2) - vO(2) - crossProduct_v(2)
      vR(3) = vI(3) - vO(3) - crossProduct_v(3)
    
    END SUBROUTINE ConvertToMovingRotatingFrame

    subroutine coriolis_force(omega, vr, fcoriolis)
      implicit none
        
      ! Arguments:
      ! omega      : Angular velocity of the rotating frame (3-element array)
      ! vr         : Velocity in the rotating frame (3-element array)
      ! fcoriolis  : Coriolis force (Output, 3-element array)
    
      real(mytype), intent(in)  :: omega(3), vr(3)
      real(mytype), intent(out) :: fcoriolis(3)
    
      ! Compute 2 * omega x vr (cross product)
      fcoriolis(1) = 2.0_mytype * (omega(2)*vr(3) - omega(3)*vr(2))
      fcoriolis(2) = 2.0_mytype * (omega(3)*vr(1) - omega(1)*vr(3))
      fcoriolis(3) = 2.0_mytype * (omega(1)*vr(2) - omega(2)*vr(1))
    
    end subroutine coriolis_force

    subroutine centrifugal_force(omega, r, fcentrifugal)
      implicit none
    
      ! Parameters:
      integer, parameter :: mytype = selected_real_kind(p=15) ! Assuming double precision
    
      ! Arguments:
      ! omega        : Angular velocity of the rotating frame (3-element array)
      ! r            : Position vector in the rotating frame (3-element array)
      ! fcentrifugal : Centrifugal force (Output, 3-element array)
    
      real(mytype), intent(in)  :: omega(3), r(3)
      real(mytype), intent(out) :: fcentrifugal(3)
      real(mytype) :: cross_product_omega_r(3)
    
      ! Compute omega x r (cross product)
      cross_product_omega_r(1) = omega(2)*r(3) - omega(3)*r(2)
      cross_product_omega_r(2) = omega(3)*r(1) - omega(1)*r(3)
      cross_product_omega_r(3) = omega(1)*r(2) - omega(2)*r(1)
    
      ! Compute fcentrifugal = -omega x (omega x r)
      fcentrifugal(1) = -(omega(2)*cross_product_omega_r(3) - omega(3)*cross_product_omega_r(2))
      fcentrifugal(2) = -(omega(3)*cross_product_omega_r(1) - omega(1)*cross_product_omega_r(3))
      fcentrifugal(3) = -(omega(1)*cross_product_omega_r(2) - omega(2)*cross_product_omega_r(1))
    
    end subroutine centrifugal_force
    
    
    
    subroutine invert_3x3_matrix(matrix, inverse)
      real(mytype), intent(in) :: matrix(3, 3)
      real(mytype), intent(out) :: inverse(3, 3)
      real(mytype) :: det
    
      ! Calculate the determinant of the 3x3 matrix
      det = matrix(1, 1) * (matrix(2, 2) * matrix(3, 3) - matrix(3, 2) * matrix(2, 3)) &
          - matrix(1, 2) * (matrix(2, 1) * matrix(3, 3) - matrix(3, 1) * matrix(2, 3)) &
          + matrix(1, 3) * (matrix(2, 1) * matrix(3, 2) - matrix(3, 1) * matrix(2, 2))
    
      ! Check if the determinant is zero (singular matrix)
      if (abs(det) < 1e-10) then
        write(*, *) "Matrix is singular. Inverse does not exist."
        return
      end if
    
      ! Calculate the elements of the inverse matrix using Cramer's rule
      inverse(1, 1) = (matrix(2, 2) * matrix(3, 3) - matrix(3, 2) * matrix(2, 3)) / det
      inverse(1, 2) = (matrix(1, 3) * matrix(3, 2) - matrix(3, 3) * matrix(1, 2)) / det
      inverse(1, 3) = (matrix(1, 2) * matrix(2, 3) - matrix(2, 2) * matrix(1, 3)) / det
      inverse(2, 1) = (matrix(2, 3) * matrix(3, 1) - matrix(3, 3) * matrix(2, 1)) / det
      inverse(2, 2) = (matrix(1, 1) * matrix(3, 3) - matrix(3, 1) * matrix(1, 3)) / det
      inverse(2, 3) = (matrix(1, 3) * matrix(2, 1) - matrix(2, 3) * matrix(1, 1)) / det
      inverse(3, 1) = (matrix(2, 1) * matrix(3, 2) - matrix(3, 1) * matrix(2, 2)) / det
      inverse(3, 2) = (matrix(1, 2) * matrix(3, 1) - matrix(3, 2) * matrix(1, 1)) / det
      inverse(3, 3) = (matrix(1, 1) * matrix(2, 2) - matrix(2, 1) * matrix(1, 2)) / det
    end subroutine invert_3x3_matrix

    subroutine matrix_vector_multiply(matrix, vector, result)
      real(mytype), intent(in) :: matrix(3, 3)
      real(mytype), intent(in) :: vector(3)
      real(mytype), intent(out) :: result(3)
      integer :: i, j
    
      do i = 1, 3
        result(i) = zero
        do j = 1, 3
          result(i) = result(i) + matrix(i, j) * vector(j)
        end do
      end do
    end subroutine matrix_vector_multiply
    
    

    subroutine accel_get(omega, inertia, torque_b, ang_accel)
      real(mytype),intent(in)  :: omega(4),inertia(3,3),torque_b(4) 
      real(mytype),intent(out) :: ang_accel(4)
      real(mytype)             :: inertia_inv(3,3),omega_v(3),torque_v(3),test(3),crossed(3),ang_accel_v(3)

      ! write(*,*) 'inverting ', inertia
      call invert_3x3_matrix(inertia,inertia_inv)
      omega_v=omega(2:4)
      torque_v=torque_b(2:4)
      call matrix_vector_multiply(inertia,omega_v,test)
      call CrossProduct(omega_v,test,crossed)
      call matrix_vector_multiply(inertia_inv,(torque_v-crossed),ang_accel_v)
      ang_accel(:)=0_mytype
      ang_accel(2:4)=ang_accel_v(1:3)

    end subroutine accel_get



    subroutine ang_half_step(q, omega_q, torque_vec, q_new, o_new)
      use param
      real(mytype),intent(in)  :: q(4),omega_q(4),torque_vec(3)
      real(mytype),intent(out) :: q_new(4),o_new(4)
      real(mytype)             :: inertia(3,3)
      real(mytype)             :: omega_b(4),torque_q(4),ang_accel_b(4),torque_b(4)
      real(mytype)             :: omega_n_quarter_b(4),omega_n_quarter(4),omega_n_half_b(4),omega_n_half(4)
      real(mytype)             :: q_half_predict(4)

      call lab_to_body(omega_q,q,omega_b)
      torque_q(1)=zero
      torque_q(2:4)=torque_vec(:)

      call lab_to_body(torque_q,q,torque_b)

      call accel_get(omega_b, inertia, torque_b, ang_accel_b)

      call omega_stepper(omega_b,ang_accel_b,dt*0.25,omega_n_quarter_b)
      call omega_stepper(omega_b,ang_accel_b,dt*half,omega_n_half_b)

      call body_to_lab(omega_n_quarter_b,q,omega_n_quarter)
      call orientation_stepper(q,omega_n_quarter,dt*half,q_half_predict)
      
      call body_to_lab(omega_n_half_b,q,omega_n_half)

      q_new=q_half_predict
      o_new=omega_n_half

    end subroutine ang_half_step

    subroutine ang_full_step(q,omega_q,q_half,omega_n_half,torque_vec,q_full,omega_full)
      use param
      real(mytype),intent(in)  :: q(4),omega_q(4),q_half(4),omega_n_half(4),torque_vec(3)
      real(mytype),intent(out) :: q_full(4),omega_full(4)
      real(mytype)             :: inertia(3,3)
      real(mytype)             :: omega_b(4),omega_n_half_b(4),omega_full_b(4)
      real(mytype)             :: torque_q(4),torque_b(4)
      real(mytype)             :: ang_accel_half_b(4),omega_n_half2(4)

      call lab_to_body(omega_q, q, omega_b)
      call lab_to_body(omega_n_half,q_half,omega_n_half_b)

      torque_q(1)=zero
      torque_q(2:4)=torque_vec(:)
      call lab_to_body(torque_q,q_half,torque_b)

      call accel_get(omega_n_half_b,inertia,torque_b,ang_accel_half_b)
      call body_to_lab(omega_n_half_b,q_half,omega_n_half2)

      call orientation_stepper(q,omega_n_half2,dt,q_full)

      call omega_stepper(omega_b,ang_accel_half_b,dt,omega_full_b)
      call body_to_lab(omega_full_b,q_full,omega_full)

    end subroutine ang_full_step

    subroutine ang_step(q,omega_q,torque_vec,time_step,q1,omega1)
      use param
      use ibm_param, only: inertia
      real(mytype),intent(in) :: q(4),omega_q(4),torque_vec(3),time_step
      real(mytype),intent(out):: q1(4),omega1(4)
      real(mytype)            :: torque_q(4),omega_b(4),torque_b(4),omega_half_b(4),omega1_b(4)
      real(mytype)            :: ang_accel_b(4), omega_half(4)

      ! write(*,*) 'ang_vel_lab = ', omega_q
      call lab_to_body(omega_q, q, omega_b) !convert to body frame
      ! write(*,*) 'ang_vel_b =   ', omega_b
      torque_q(1)=zero
      torque_q(2:4)=torque_vec(:)
      ! write(*,*) 'torque =      ', torque_q
      call lab_to_body(torque_q,q, torque_b)
      ! write(*,*) 'torque_b =    ', torque_b

      call accel_get(omega_b,inertia,torque_b,ang_accel_b) !calculate acceleration
      ! write(*,*) 'acceleration =', ang_accel_b

      call omega_stepper(omega_b, ang_accel_b,time_step*half,omega_half_b)
      ! write(*,*) 'omega_half_b =', omega_half_b
      call omega_stepper(omega_b, ang_accel_b,time_step,omega1_b) !calculate omega at half and full timestep
      ! write(*,*) 'omega_full_b =', omega1_b
      call body_to_lab(omega_half_b,q,omega_half) !convert back to lab
      ! write(*,*) 'omega_half_lab', omega_half
      call orientation_stepper(q,omega_half,time_step,q1) !step forward orientation
      ! write(*,*) 'time_step    =', time_step
      ! write(*,*) 'orientation1 =', q1
      call body_to_lab(omega1_b,q1,omega1)
      ! write(*,*) 'omega_full   =', omega1

    end subroutine ang_step


    subroutine lin_step(position,linearVelocity,linearForce,time_step,position_1,linearVelocity_1)
      use ibm_param, only: ellip_m
      real(mytype),intent(in)   :: position(3),linearVelocity(3),linearForce(3),time_step
      real(mytype),intent(out)  :: position_1(3),linearVelocity_1(3)
      real(mytype)              :: linearAcceleration(3)

      linearAcceleration(:) = linearForce(:) / ellip_m
      position_1(:) = position(:) + time_step*linearVelocity(:)
      linearVelocity_1 = linearVelocity(:) + time_step*linearAcceleration(:)

    end subroutine lin_step

    subroutine ellipMassCalculate(shape,rho_s,mass)
      use constants, only: pi
      real(mytype),intent(in)  :: shape(3),rho_s
      real(mytype),intent(out) :: mass
      real(mytype)             :: a,b,c,vol

      a=shape(1)
      b=shape(2)
      c=shape(3)
      vol=(4_mytype/3_mytype)*pi*a*b*c
      mass=vol*rho_s

    end subroutine ellipMassCalculate

    subroutine ellipInertiaCalculate(shape,rho_s,inertia)
      real(mytype),intent(in)  :: shape(3),rho_s
      real(mytype),intent(out) :: inertia(3,3)
      real(mytype)             :: a,b,c,i1,i2,i3,mass

      call ellipMassCalculate(shape,rho_s,mass)

      a=shape(1)
      b=shape(2)
      c=shape(3)

      i1=mass*(b**2+c**2)*0.2
      i2=mass*(a**2+c**2)*0.2
      i3=mass*(a**2+b**2)*0.2

      inertia(:,:)=0_mytype
      inertia(1,1)=i1
      inertia(2,2)=i2
      inertia(3,3)=i3

    end subroutine ellipInertiaCalculate
      
end module ellipsoid_utils
