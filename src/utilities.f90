!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module utilities

  use decomp_2d_constants, only : mytype

  implicit none

  private

  public :: rl, iy, cx,gen_snapshotname,gen_h5path,int_to_str,gen_filename

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
 function gen_snapshotname(pathname, varname, num, ext)
    character(len=*), intent(in) :: pathname, varname, ext
    integer, intent(in) :: num
    character(len=:), allocatable :: gen_snapshotname
#ifndef ADIOS2
    gen_snapshotname = gen_filename(pathname, varname, num, ext)
#else
    gen_snapshotname = varname//'-'//int_to_str(num)//'.'//ext
#endif
  end function gen_snapshotname
  
  function gen_filename(pathname, varname, num, ext)

    character(len=*), intent(in) :: pathname, varname, ext
    integer, intent(in) :: num
#ifndef ADIOS2
    character(len=:), allocatable :: gen_filename
    gen_filename = pathname//'/'//varname//'-'//int_to_str(num)//'.'//ext
#else
    character(len=len(varname)) :: gen_filename
    write(gen_filename, "(A)") varname
#endif
    
  end function gen_filename

  function gen_h5path(filename, num)

    character(len=*), intent(in) :: filename
    integer, intent(in) :: num
#ifndef ADIOS2
    character(len=*), parameter :: path_to_h5file = "./"
    character(len=(len(path_to_h5file) + len(filename))) :: gen_h5path
    write(gen_h5path, "(A)") path_to_h5file//filename
#else
    character(len=*), parameter :: path_to_h5file = "../data.hdf5:/Step"
    character(len=:), allocatable :: gen_h5path
    gen_h5path = path_to_h5file//int_to_str(num)//"/"//filename
#endif
    
  end function gen_h5path

  ! Function converting an integer to a string.
  function int_to_str(i)

    integer, intent(in) :: i ! Integer input.

    ! String return value. The string must be long enough to contain all the characters required to
    ! represent the integer, i.e 1 + log_10(i). To protect against calling with integer 0 the value
    ! passed to log_10 must be >= 1.
    character(len=(1 + int(log10(real(max(i, 1)))))) :: int_to_str

    write(int_to_str, "(I0)") i
  end function int_to_str

end module utilities
