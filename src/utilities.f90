!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module utilities

  use decomp_2d_constants, only : mytype

  implicit none

  private

  public :: rl, iy, cx

contains

   !##################################################################
   function rl(complexnumber)

  implicit none

  real(mytype) :: rl
  complex(mytype) :: complexnumber

  rl = real(complexnumber, kind=mytype)

end function rl
!##################################################################
!##################################################################
function iy(complexnumber)

  implicit none

  real(mytype) :: iy
  complex(mytype) :: complexnumber

  iy = aimag(complexnumber)

end function iy
!##################################################################
!##################################################################
function cx(realpart,imaginarypart)

  implicit none

  complex(mytype) :: cx
  real(mytype) :: realpart, imaginarypart

  cx = cmplx(realpart, imaginarypart, kind=mytype)

end function cx
!##################################################################

end module utilities
