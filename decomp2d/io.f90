!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2013 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This module provides parallel IO facilities for applications based on
! 2D decomposition.

module decomp_2d_io

  use decomp_2d
  use MPI
#ifdef T3PIO
  use t3pio
#endif

  implicit none

  private        ! Make everything private unless declared public

  public :: decomp_2d_write_one, decomp_2d_read_one, &
       decomp_2d_write_var, decomp_2d_read_var, &
       decomp_2d_write_scalar, decomp_2d_read_scalar, &
       decomp_2d_write_plane, decomp_2d_write_every, &
       decomp_2d_write_subdomain

  ! Generic interface to handle multiple data types

  interface decomp_2d_write_one
     module procedure write_one_real
     module procedure write_one_complex
     module procedure mpiio_write_real_coarse
     module procedure mpiio_write_real_probe
  end interface decomp_2d_write_one

  interface decomp_2d_read_one
     module procedure read_one_real
     module procedure read_one_complex
  end interface decomp_2d_read_one

  interface decomp_2d_write_var
     module procedure write_var_real
     module procedure write_var_complex
  end interface decomp_2d_write_var

  interface decomp_2d_read_var
     module procedure read_var_real
     module procedure read_var_complex
  end interface decomp_2d_read_var

  interface decomp_2d_write_scalar
     module procedure write_scalar_real
     module procedure write_scalar_complex
     module procedure write_scalar_integer
     module procedure write_scalar_logical
  end interface decomp_2d_write_scalar

  interface decomp_2d_read_scalar
     module procedure read_scalar_real
     module procedure read_scalar_complex
     module procedure read_scalar_integer
     module procedure read_scalar_logical
  end interface decomp_2d_read_scalar

  interface decomp_2d_write_plane
     module procedure write_plane_3d_real
     module procedure write_plane_3d_complex
     !     module procedure write_plane_2d
  end interface decomp_2d_write_plane

  interface decomp_2d_write_every
     module procedure write_every_real
     module procedure write_every_complex
  end interface decomp_2d_write_every

  interface decomp_2d_write_subdomain
     module procedure write_subdomain
  end interface decomp_2d_write_subdomain

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Using MPI-IO library to write a single 3D array to a file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_one_real(ipencil,var,filename,opt_decomp)

    implicit none

    integer, intent(IN) :: ipencil
    real(mytype), dimension(:,:,:), intent(IN) :: var
    character(len=*), intent(IN) :: filename
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, fh, data_type, info, gs

    data_type = real_type

#include "io_write_one.inc"

    return
  end subroutine write_one_real


  subroutine write_one_complex(ipencil,var,filename,opt_decomp)

    implicit none

    integer, intent(IN) :: ipencil
    complex(mytype), dimension(:,:,:), intent(IN) :: var
    character(len=*), intent(IN) :: filename
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, fh, data_type, info, gs

    data_type = complex_type

#include "io_write_one.inc"

    return
  end subroutine write_one_complex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Using MPI-IO library to read from a file a single 3D array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_one_real(ipencil,var,filename,opt_decomp)

    implicit none

    integer, intent(IN) :: ipencil
    real(mytype), dimension(:,:,:), intent(INOUT) :: var
    character(len=*), intent(IN) :: filename
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, fh, data_type
    real(mytype_single), allocatable, dimension(:,:,:) :: varsingle

    data_type = real_type_single
    allocate (varsingle(xstV(1):xenV(1),xstV(2):xenV(2),xstV(3):xenV(3)))

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       call get_decomp_info(decomp)
    end if

    ! determine subarray parameters
    sizes(1) = decomp%xsz(1)
    sizes(2) = decomp%ysz(2)
    sizes(3) = decomp%zsz(3)

    if (ipencil == 1) then
       subsizes(1) = decomp%xsz(1)
       subsizes(2) = decomp%xsz(2)
       subsizes(3) = decomp%xsz(3)
       starts(1) = decomp%xst(1)-1  ! 0-based index
       starts(2) = decomp%xst(2)-1
       starts(3) = decomp%xst(3)-1
    else if (ipencil == 2) then
       subsizes(1) = decomp%ysz(1)
       subsizes(2) = decomp%ysz(2)
       subsizes(3) = decomp%ysz(3)
       starts(1) = decomp%yst(1)-1
       starts(2) = decomp%yst(2)-1
       starts(3) = decomp%yst(3)-1
    else if (ipencil == 3) then
       subsizes(1) = decomp%zsz(1)
       subsizes(2) = decomp%zsz(2)
       subsizes(3) = decomp%zsz(3)
       starts(1) = decomp%zst(1)-1
       starts(2) = decomp%zst(2)-1
       starts(3) = decomp%zst(3)-1
    endif

    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, data_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, &
         fh, ierror)
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,data_type, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, varsingle, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         data_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    var = real(varsingle,mytype)
    deallocate(varsingle)

    return
  end subroutine read_one_real


  subroutine read_one_complex(ipencil,var,filename,opt_decomp)

    implicit none

    integer, intent(IN) :: ipencil
    complex(mytype), dimension(:,:,:), intent(INOUT) :: var
    character(len=*), intent(IN) :: filename
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, fh, data_type

    data_type = complex_type

#include "io_read_one.inc"

    return
  end subroutine read_one_complex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a 3D array as part of a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated after the writing
  !  operation to prepare the writing of next chunk of data.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_var_real(fh,disp,ipencil,var,opt_decomp)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: ipencil
    real(mytype), dimension(:,:,:), intent(IN) :: var
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, data_type

    data_type = real_type

#include "io_write_var.inc"

    return
  end subroutine write_var_real


  subroutine write_var_complex(fh,disp,ipencil,var,opt_decomp)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: ipencil
    complex(mytype), dimension(:,:,:), intent(IN) :: var
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, data_type

    data_type = complex_type

#include "io_write_var.inc"

    return
  end subroutine write_var_complex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read a 3D array as part of a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated after the reading
  !  operation to prepare the reading of next chunk of data.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_var_real(fh,disp,ipencil,var,opt_decomp)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: ipencil
    real(mytype), dimension(:,:,:), intent(INOUT) :: var
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, data_type

    data_type = real_type

#include "io_read_var.inc"

    return
  end subroutine read_var_real


  subroutine read_var_complex(fh,disp,ipencil,var,opt_decomp)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: ipencil
    complex(mytype), dimension(:,:,:), intent(INOUT) :: var
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, data_type

    data_type = complex_type

#include "io_read_var.inc"

    return
  end subroutine read_var_complex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write scalar variables as part of a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated after the reading
  !  operation to prepare the reading of next chunk of data.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_scalar_real(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh             ! file handle
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement
    integer, intent(IN) :: n              ! number of scalars
    real(mytype), dimension(n), &
         intent(IN) :: var                ! array of scalars

    integer :: m, ierror

    call MPI_FILE_SET_VIEW(fh,disp,real_type, &
         real_type,'native',MPI_INFO_NULL,ierror)
    if (nrank==0) then
       m = n ! only one rank needs to write
    else
       m = 0
    end if
    call MPI_FILE_WRITE_ALL(fh, var, m, real_type, &
         MPI_STATUS_IGNORE, ierror)
    disp = disp + n*mytype_bytes

    return
  end subroutine write_scalar_real


  subroutine write_scalar_complex(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: n
    complex(mytype), dimension(n), intent(IN) :: var

    integer :: m, ierror

    call MPI_FILE_SET_VIEW(fh,disp,complex_type, &
         complex_type,'native',MPI_INFO_NULL,ierror)
    if (nrank==0) then
       m = n
    else
       m = 0
    end if
    call MPI_FILE_WRITE_ALL(fh, var, m, complex_type, &
         MPI_STATUS_IGNORE, ierror)
    disp = disp + n*mytype_bytes*2

    return
  end subroutine write_scalar_complex


  subroutine write_scalar_integer(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: n
    integer, dimension(n), intent(IN) :: var

    integer :: m, ierror

    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         MPI_INTEGER,'native',MPI_INFO_NULL,ierror)
    if (nrank==0) then
       m = n
    else
       m = 0
    end if
    call MPI_FILE_WRITE_ALL(fh, var, m, MPI_INTEGER, &
         MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_SIZE(MPI_INTEGER,m,ierror)
    disp = disp + n*m

    return
  end subroutine write_scalar_integer


  subroutine write_scalar_logical(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: n
    logical, dimension(n), intent(IN) :: var

    integer :: m, ierror

    call MPI_FILE_SET_VIEW(fh,disp,MPI_LOGICAL, &
         MPI_LOGICAL,'native',MPI_INFO_NULL,ierror)
    if (nrank==0) then
       m = n
    else
       m = 0
    end if
    call MPI_FILE_WRITE_ALL(fh, var, m, MPI_LOGICAL, &
         MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_SIZE(MPI_LOGICAL,m,ierror)
    disp = disp + n*m

    return
  end subroutine write_scalar_logical


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read scalar variables as part of a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated after the reading
  !  operation to prepare the reading of next chunk of data.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_scalar_real(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh             ! file handle
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement
    integer, intent(IN) :: n              ! number of scalars
    real(mytype), dimension(n), &
         intent(INOUT) :: var             ! array of scalars

    integer :: ierror

    call MPI_FILE_SET_VIEW(fh,disp,real_type, &
         real_type,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, n, real_type, &
         MPI_STATUS_IGNORE, ierror)
    disp = disp + n*mytype_bytes

    return
  end subroutine read_scalar_real


  subroutine read_scalar_complex(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: n
    complex(mytype), dimension(n), intent(INOUT) :: var

    integer :: ierror

    call MPI_FILE_SET_VIEW(fh,disp,complex_type, &
         complex_type,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, n, complex_type, &
         MPI_STATUS_IGNORE, ierror)
    disp = disp + n*mytype_bytes*2

    return
  end subroutine read_scalar_complex


  subroutine read_scalar_integer(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: n
    integer, dimension(n), intent(INOUT) :: var

    integer :: m, ierror

    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         MPI_INTEGER,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, n, MPI_INTEGER, &
         MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_SIZE(MPI_INTEGER,m,ierror)
    disp = disp + n*m

    return
  end subroutine read_scalar_integer


  subroutine read_scalar_logical(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: n
    logical, dimension(n), intent(INOUT) :: var

    integer :: m, ierror

    call MPI_FILE_SET_VIEW(fh,disp,MPI_LOGICAL, &
         MPI_LOGICAL,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, n, MPI_LOGICAL, &
         MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_SIZE(MPI_LOGICAL,m,ierror)
    disp = disp + n*m

    return
  end subroutine read_scalar_logical


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a 2D slice of the 3D data to a file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_plane_3d_real(ipencil,var,iplane,n,filename, &
       opt_decomp)

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    real(mytype), dimension(:,:,:), intent(IN) :: var
    integer, intent(IN) :: iplane !(x-plane=1; y-plane=2; z-plane=3)
    integer, intent(IN) :: n ! which plane to write (global coordinate)
    character(len=*), intent(IN) :: filename
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    real(mytype), allocatable, dimension(:,:,:) :: wk2d
    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh, data_type

    data_type = real_type

#include "io_write_plane.inc"

    return
  end subroutine write_plane_3d_real


  subroutine write_plane_3d_complex(ipencil,var,iplane,n, &
       filename,opt_decomp)

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    complex(mytype), dimension(:,:,:), intent(IN) :: var
    integer, intent(IN) :: iplane !(x-plane=1; y-plane=2; z-plane=3)
    integer, intent(IN) :: n ! which plane to write (global coordinate)
    character(len=*), intent(IN) :: filename
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    complex(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    complex(mytype), allocatable, dimension(:,:,:) :: wk2d
    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh, data_type

    data_type = complex_type

#include "io_write_plane.inc"

    return
  end subroutine write_plane_3d_complex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a 2D array to a file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !************** TO DO ***************
  !* Consider handling distributed 2D data set
  !  subroutine write_plane_2d(ipencil,var,filename)
  !    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
  !    real(mytype), dimension(:,:), intent(IN) :: var ! 2D array
  !    character(len=*), intent(IN) :: filename
  !
  !    if (ipencil==1) then
  !       ! var should be defined as var(xsize(2)
  !
  !    else if (ipencil==2) then
  !
  !    else if (ipencil==3) then
  !
  !    end if
  !
  !    return
  !  end subroutine write_plane_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write 3D array data for every specified mesh point
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_every_real(ipencil,var,iskip,jskip,kskip, &
       filename, from1)

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    real(mytype), dimension(:,:,:), intent(IN) :: var
    integer, intent(IN) :: iskip,jskip,kskip 
    character(len=*), intent(IN) :: filename
    logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
    ! .false. - save n,2n,3n...

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh, key,color,newcomm, data_type
    integer, dimension(3) :: xsz,ysz,zsz,xst,yst,zst,xen,yen,zen,skip

    data_type = real_type

#include "io_write_every.inc"

    return
  end subroutine write_every_real


  subroutine write_every_complex(ipencil,var,iskip,jskip,kskip, &
       filename, from1)

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    complex(mytype), dimension(:,:,:), intent(IN) :: var
    integer, intent(IN) :: iskip,jskip,kskip 
    character(len=*), intent(IN) :: filename
    logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
    ! .false. - save n,2n,3n...

    complex(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh, key,color,newcomm, data_type
    integer, dimension(3) :: xsz,ysz,zsz,xst,yst,zst,xen,yen,zen,skip

    data_type = complex_type

#include "io_write_every.inc"

    return
  end subroutine write_every_complex


  subroutine mpiio_write_real_coarse(ipencil,var,filename,icoarse)

    ! USE param
    ! USE variables

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    integer, intent(IN) :: icoarse !(nstat=1; nvisu=2)
    real(mytype), dimension(:,:,:), intent(IN) :: var
    real(mytype_single), allocatable, dimension(:,:,:) :: varsingle
    character(len=*) :: filename

    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh

    if (icoarse==1) then
       sizes(1) = xszS(1)
       sizes(2) = yszS(2)
       sizes(3) = zszS(3)

       if (ipencil == 1) then
          subsizes(1) = xszS(1)
          subsizes(2) = xszS(2)
          subsizes(3) = xszS(3)
          starts(1) = xstS(1)-1  ! 0-based index
          starts(2) = xstS(2)-1
          starts(3) = xstS(3)-1
       else if (ipencil == 2) then
          subsizes(1) = yszS(1)
          subsizes(2) = yszS(2)
          subsizes(3) = yszS(3)
          starts(1) = ystS(1)-1
          starts(2) = ystS(2)-1
          starts(3) = ystS(3)-1
       else if (ipencil == 3) then
          subsizes(1) = zszS(1)
          subsizes(2) = zszS(2)
          subsizes(3) = zszS(3)
          starts(1) = zstS(1)-1
          starts(2) = zstS(2)-1
          starts(3) = zstS(3)-1
       endif
    endif

    if (icoarse==2) then
       sizes(1) = xszV(1)
       sizes(2) = yszV(2)
       sizes(3) = zszV(3)

       if (ipencil == 1) then
          subsizes(1) = xszV(1)
          subsizes(2) = xszV(2)
          subsizes(3) = xszV(3)
          starts(1) = xstV(1)-1  ! 0-based index
          starts(2) = xstV(2)-1
          starts(3) = xstV(3)-1
       else if (ipencil == 2) then
          subsizes(1) = yszV(1)
          subsizes(2) = yszV(2)
          subsizes(3) = yszV(3)
          starts(1) = ystV(1)-1
          starts(2) = ystV(2)-1
          starts(3) = ystV(3)-1
       else if (ipencil == 3) then
          subsizes(1) = zszV(1)
          subsizes(2) = zszV(2)
          subsizes(3) = zszV(3)
          starts(1) = zstV(1)-1
          starts(2) = zstV(2)-1
          starts(3) = zstV(3)-1
       endif
    endif

    allocate (varsingle(xstV(1):xenV(1),xstV(2):xenV(2),xstV(3):xenV(3)))
    varsingle=var
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, real_type_single, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,real_type_single, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, varsingle, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         real_type_single, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    deallocate(varsingle)

    return
  end subroutine mpiio_write_real_coarse

  subroutine mpiio_write_real_probe(ipencil,var,filename,nlength)

    ! USE param
    ! USE variables

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    integer, intent(in) :: nlength
    real(mytype), dimension(:,:,:,:), intent(IN) :: var

    character(len=*) :: filename

    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(4) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh


    sizes(1) = xszP(1)
    sizes(2) = yszP(2)
    sizes(3) = zszP(3)
    sizes(4) = nlength
    if (ipencil == 1) then
       subsizes(1) = xszP(1)
       subsizes(2) = xszP(2)
       subsizes(3) = xszP(3)
       subsizes(4) = nlength
       starts(1) = xstP(1)-1  ! 0-based index
       starts(2) = xstP(2)-1
       starts(3) = xstP(3)-1
       starts(4) = 0
    else if (ipencil == 2) then
       subsizes(1) = yszP(1)
       subsizes(2) = yszP(2)
       subsizes(3) = yszP(3)
       starts(1) = ystP(1)-1
       starts(2) = ystP(2)-1
       starts(3) = ystP(3)-1
    else if (ipencil == 3) then
       subsizes(1) = zszP(1)
       subsizes(2) = zszP(2)
       subsizes(3) = zszP(3)
       starts(1) = zstP(1)-1
       starts(2) = zstP(2)-1
       starts(3) = zstP(3)-1
    endif
    !   print *,nrank,starts(1),starts(2),starts(3),starts(4)
    call MPI_TYPE_CREATE_SUBARRAY(4, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,real_type, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3)*subsizes(4), &
         real_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)


    return
  end subroutine mpiio_write_real_probe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a 3D data set covering a smaller sub-domain only
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_subdomain(ipencil,var,is,ie,js,je,ks,ke,filename)

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    real(mytype), dimension(:,:,:), intent(IN) :: var
    integer, intent(IN) :: is, ie, js, je, ks, ke
    character(len=*), intent(IN) :: filename

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: color, key, errorcode, newcomm, ierror
    integer :: newtype, fh, data_type, i, j, k
    integer :: i1, i2, j1, j2, k1, k2

    data_type = real_type

    ! validate the input paramters
    if (is<1 .OR. ie>nx_global .OR. js<1 .OR. je>ny_global .OR. &
         ks<1 .OR. ke>nz_global) then
       errorcode = 10
       call decomp_2d_abort(errorcode, &
            'Invalid subdomain specified in I/O')
    end if

    ! create a communicator for all those MPI ranks containing the subdomain
    color = 1
    key = 1
    if (ipencil==1) then
       if (xstart(1)>ie .OR. xend(1)<is .OR. xstart(2)>je .OR. xend(2)<js &
            .OR. xstart(3)>ke .OR. xend(3)<ks) then
          color = 2
       end if
    else if (ipencil==2) then
       if (ystart(1)>ie .OR. yend(1)<is .OR. ystart(2)>je .OR. yend(2)<js &
            .OR. ystart(3)>ke .OR. yend(3)<ks) then
          color = 2
       end if
    else if (ipencil==3) then
       if (zstart(1)>ie .OR. zend(1)<is .OR. zstart(2)>je .OR. zend(2)<js &
            .OR. zstart(3)>ke .OR. zend(3)<ks) then
          color = 2
       end if
    end if
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,newcomm,ierror)

    if (color==1) then ! only ranks in this group do IO collectively

       ! generate MPI-IO subarray information

       ! global size of the sub-domain to write
       sizes(1) = ie - is + 1
       sizes(2) = je - js + 1
       sizes(3) = ke - ks + 1

       ! 'subsizes' & 'starts' as required by MPI_TYPE_CREATE_SUBARRAY
       ! note the special code whe subdomain only occupy part of the pencil
       if (ipencil==1) then

          subsizes(1) = xsize(1)
          starts(1) = xstart(1) - is
          if (xend(1)>ie .AND. xstart(1)<is) then
             subsizes(1) = ie - is + 1
             starts(1) = 0
          else if (xstart(1)<is) then
             subsizes(1) = xend(1) - is + 1
             starts(1) = 0
          else if (xend(1)>ie) then
             subsizes(1) = ie - xstart(1) + 1
          end if
          subsizes(2) = xsize(2)
          starts(2) = xstart(2) - js
          if (xend(2)>je .AND. xstart(2)<js) then
             subsizes(2) = je - js + 1
             starts(2) = 0
          else if (xstart(2)<js) then
             subsizes(2) = xend(2) - js + 1
             starts(2) = 0
          else if (xend(2)>je) then
             subsizes(2) = je - xstart(2) + 1
          end if
          subsizes(3) = xsize(3)
          starts(3) = xstart(3) - ks
          if (xend(3)>ke .AND. xstart(3)<ks) then
             subsizes(3) = ke - ks + 1
             starts(3) = 0
          else if (xstart(3)<ks) then
             subsizes(3) = xend(3) - ks + 1
             starts(3) = 0
          else if (xend(3)>ke) then
             subsizes(3) = ke - xstart(3) + 1
          end if

       else if (ipencil==2) then

          ! TODO 

       else if (ipencil==3) then

          ! TODO

       end if


       ! copy data from orginal to a temp array
       ! pay attention to blocks only partially cover the sub-domain
       if (ipencil==1) then

          if (xend(1)>ie .AND. xstart(1)<is) then
             i1 = is
             i2 = ie
          else if (xend(1)>ie) then
             i1 = xstart(1)
             i2 = ie
          else if (xstart(1)<is) then
             i1 = is
             i2 = xend(1)
          else
             i1 = xstart(1)
             i2 = xend(1)
          end if

          if (xend(2)>je .AND. xstart(2)<js) then
             j1 = js
             j2 = je
          else if (xend(2)>je) then
             j1 = xstart(2)
             j2 = je
          else if (xstart(2)<js) then
             j1 = js
             j2 = xend(2)
          else
             j1 = xstart(2)
             j2 = xend(2)
          end if

          if (xend(3)>ke .AND. xstart(3)<ks) then
             k1 = ks
             k2 = ke
          else if (xend(3)>ke) then
             k1 = xstart(3)
             k2 = ke
          else if (xstart(3)<ks) then
             k1 = ks
             k2 = xend(3)
          else
             k1 = xstart(3)
             k2 = xend(3)
          end if

          allocate(wk(i1:i2, j1:j2, k1:k2))
          allocate(wk2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
          wk2 = var
          do k=k1,k2
             do j=j1,j2
                do i=i1,i2
                   wk(i,j,k) = wk2(i,j,k)
                end do
             end do
          end do

       else if (ipencil==2) then

          ! TODO

       else if (ipencil==3) then

          ! TODO

       end if

       deallocate(wk2)

       ! MPI-IO
       call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
            MPI_ORDER_FORTRAN, data_type, newtype, ierror)
       call MPI_TYPE_COMMIT(newtype,ierror)
       call MPI_FILE_OPEN(newcomm, filename, &
            MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
            fh, ierror)
       filesize = 0_MPI_OFFSET_KIND
       call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
       disp = 0_MPI_OFFSET_KIND
       call MPI_FILE_SET_VIEW(fh,disp,data_type, &
            newtype,'native',MPI_INFO_NULL,ierror)
       call MPI_FILE_WRITE_ALL(fh, wk, &
            subsizes(1)*subsizes(2)*subsizes(3), &
            data_type, MPI_STATUS_IGNORE, ierror)
       call MPI_FILE_CLOSE(fh,ierror)
       call MPI_TYPE_FREE(newtype,ierror)

       deallocate(wk)

    end if

    return
  end subroutine write_subdomain


end module decomp_2d_io
