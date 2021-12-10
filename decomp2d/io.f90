!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2013 Ning Li, the Numerical Algorithms Group (NAG)
! Copyright (C) 2021               the University of Edinburgh (UoE)
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

#ifdef ADIOS2
  use adios2
#endif

  implicit none

  integer, parameter, public :: decomp_2d_write_mode = 1, decomp_2d_read_mode = 2, &
       decomp_2d_append_mode = 3
  integer, parameter :: MAX_IOH = 10 ! How many live IO things should we handle?
  character(len=*), parameter :: io_sep = "::"
  integer, save :: nreg_io = 0
#ifndef ADIOS2
  integer, dimension(MAX_IOH), save :: fh_registry
  logical, dimension(MAX_IOH), target, save :: fh_live
  character(len=80), dimension(MAX_IOH), target, save :: fh_names
  integer(kind=MPI_OFFSET_KIND), dimension(MAX_IOH), save :: fh_disp
#else
  type(adios2_adios) :: adios
  character(len=80), dimension(MAX_IOH), target, save :: engine_names
  logical, dimension(MAX_IOH), target, save :: engine_live
  type(adios2_engine), dimension(MAX_IOH), save :: engine_registry
#endif
  
  private        ! Make everything private unless declared public

  public :: decomp_2d_write_one, decomp_2d_read_one, &
       decomp_2d_write_var, decomp_2d_read_var, &
       decomp_2d_write_scalar, decomp_2d_read_scalar, &
       decomp_2d_write_plane, decomp_2d_write_every, &
       decomp_2d_write_subdomain, &
       decomp_2d_write_outflow, decomp_2d_read_inflow, &
       decomp_2d_io_init, decomp_2d_io_finalise, & ! XXX: initialise/finalise 2decomp&fft IO module
       decomp_2d_init_io, & ! XXX: initialise an io process - awful naming
       decomp_2d_register_variable, &
       decomp_2d_open_io, decomp_2d_close_io, &
       decomp_2d_start_io, decomp_2d_end_io, &
       gen_iodir_name

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

  interface decomp_2d_write_outflow
     module procedure write_outflow
  end interface decomp_2d_write_outflow

  interface decomp_2d_read_inflow
     module procedure read_inflow
  end interface decomp_2d_read_inflow

contains

  subroutine decomp_2d_io_init()

#ifdef ADIOS2
  integer :: ierror
  logical :: adios2_debug_mode
  character(len=80) :: config_file="adios2_config.xml"
#endif

#ifdef ADIOS2
  !! TODO: make this a runtime-option
  adios2_debug_mode = .true.

  call adios2_init(adios, trim(config_file), MPI_COMM_WORLD, adios2_debug_mode, ierror)
  if (ierror.ne.0) then
     print *, "Error initialising ADIOS2 - is adios2_config.xml present and valid?"
     call MPI_ABORT(MPI_COMM_WORLD, -1, ierror)
  endif

  engine_live(:) = .false.
#endif
  
  end subroutine decomp_2d_io_init
  subroutine decomp_2d_io_finalise()

#ifdef ADIOS2
    use adios2
#endif

    implicit none

#ifdef ADIOS2
    integer :: ierror
#endif
    
#ifdef ADIOS2
    call adios2_finalize(adios, ierror)
#endif

  end subroutine decomp_2d_io_finalise
  
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
  subroutine read_one_real(ipencil,var,dirname,varname,io_name,opt_decomp,reduce_prec)

    implicit none

    integer, intent(IN) :: ipencil
    real(mytype), dimension(:,:,:), intent(INOUT) :: var
    character(len=*), intent(IN) :: varname, dirname, io_name
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    logical, intent(in), optional :: reduce_prec

    logical :: read_reduce_prec
    
    TYPE(DECOMP_INFO) :: decomp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, data_type
    real(mytype_single), allocatable, dimension(:,:,:) :: varsingle
    integer :: idx
    integer :: disp_bytes

    read_reduce_prec = .true.
    
    idx = get_io_idx(io_name, dirname)
#ifndef ADIOS2
    if (present(reduce_prec)) then
       if (.not. reduce_prec) then
          read_reduce_prec = .false.
       end if
    end if
    if (read_reduce_prec) then
       data_type = real_type_single
    else
       data_type = real_type
    end if
    call MPI_TYPE_SIZE(data_type,disp_bytes,ierror)

    !! Use MPIIO
    if (read_reduce_prec) then
       allocate (varsingle(xstV(1):xenV(1),xstV(2):xenV(2),xstV(3):xenV(3)))
    end if
    
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

    associate(fh => fh_registry(idx), &
         disp => fh_disp(idx))
      call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
           MPI_ORDER_FORTRAN, data_type, newtype, ierror)
      call MPI_TYPE_COMMIT(newtype,ierror)
      call MPI_FILE_SET_VIEW(fh,disp,data_type, &
           newtype,'native',MPI_INFO_NULL,ierror)
      if (read_reduce_prec) then
         call MPI_FILE_READ_ALL(fh, varsingle, &
              subsizes(1)*subsizes(2)*subsizes(3), &
              data_type, MPI_STATUS_IGNORE, ierror)
         var = real(varsingle,mytype)
         deallocate(varsingle)
      else
         call MPI_FILE_READ_ALL(fh, var, &
              subsizes(1)*subsizes(2)*subsizes(3), &
              data_type, MPI_STATUS_IGNORE, ierror)
      endif
      call MPI_TYPE_FREE(newtype,ierror)

      disp = disp + sizes(1) * sizes(2) * sizes(3) * disp_bytes
    end associate
#else
    call adios2_read_one_real(ipencil, var, dirname, varname, io_name)
#endif
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

#ifdef ADIOS2
  subroutine adios2_read_one_real(ipencil,var,engine_name,varname,io_name)

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    character(len=*), intent(in) :: engine_name
    character(len=*), intent(in) :: io_name
    character*(*), intent(in) :: varname
    real(mytype), dimension(:,:,:), intent(out) :: var
    
    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
    integer :: i,j,k, ierror, newtype, fh
    type(adios2_io) :: io_handle
    type(adios2_variable) :: var_handle
    integer :: idx

    call adios2_at_io(io_handle, adios, io_name, ierror)
    call adios2_inquire_variable(var_handle, io_handle, varname, ierror)
    if (.not.var_handle % valid) then
       print *, "ERROR: trying to read variable without registering first! ", varname
       stop
    endif

    idx = get_io_idx(io_name, engine_name)
    call adios2_get(engine_registry(idx), var_handle, var, adios2_mode_deferred, ierror)

    return
    
  end subroutine adios2_read_one_real
#endif
  
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


  subroutine write_outflow(dirname,varname,ntimesteps,var,io_name,opt_decomp)

    implicit none

    character(len=*), intent(in) :: dirname, varname, io_name
    integer, intent(IN) :: ntimesteps 
    real(mytype), dimension(:,:,:), intent(IN) :: var
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, data_type
    integer :: idx
#ifdef ADIOS2
    type(adios2_io) :: io_handle
    type(adios2_variable) :: var_handle
#endif

    data_type = real_type

#include "io_write_outflow.f90"

    return
  end subroutine write_outflow


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


  subroutine read_inflow(dirname,varname,ntimesteps,var,io_name,opt_decomp)
 
    implicit none

    character(len=*), intent(in) :: dirname, varname, io_name
    integer, intent(IN) :: ntimesteps
    real(mytype), dimension(:,:,:), intent(INOUT) :: var
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, data_type
    integer :: idx
#ifdef ADIOS2
    type(adios2_io) :: io_handle
    type(adios2_variable) :: var_handle
#endif
    
    data_type = real_type

#include "io_read_inflow.f90"

    return
  end subroutine read_inflow

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
  subroutine plane_extents (sizes, subsizes, starts, iplane, opt_decomp, opt_nplanes)

    integer, intent(in) :: iplane
    type(decomp_info), intent(in), optional :: opt_decomp
    integer, intent(in), optional :: opt_nplanes
    
    integer, dimension(3), intent(out) :: sizes, subsizes, starts

    integer :: nplanes
    type(decomp_info) :: decomp

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       call get_decomp_info(decomp)
    end if

    if (present(opt_nplanes)) then
       nplanes = opt_nplanes
    else
       nplanes = 1
    end if
    
    if (iplane == 1) then
       sizes(1) = nplanes
       sizes(2) = decomp%ysz(2)
       sizes(3) = decomp%zsz(3)
       subsizes(1) = nplanes
       subsizes(2) = decomp%xsz(2)
       subsizes(3) = decomp%xsz(3)
       starts(1) = 0
       starts(2) = decomp%xst(2)-1
       starts(3) = decomp%xst(3)-1
    else if (iplane == 2) then
       sizes(1) = decomp%xsz(1)
       sizes(2) = nplanes
       sizes(3) = decomp%zsz(3)
       subsizes(1) = decomp%ysz(1)
       subsizes(2) = nplanes
       subsizes(3) = decomp%ysz(3)
       starts(1) = decomp%yst(1)-1
       starts(2) = 0
       starts(3) = decomp%yst(3)-1
    else if (iplane == 3) then
       sizes(1) = decomp%xsz(1)
       sizes(2) = decomp%ysz(2)
       sizes(3) = nplanes
       subsizes(1) = decomp%zsz(1)
       subsizes(2) = decomp%zsz(2)
       subsizes(3) = nplanes
       starts(1) = decomp%zst(1)-1
       starts(2) = decomp%zst(2)-1
       starts(3) = 0
    else
       print *, "Can't work with plane ", iplane
       stop
    endif
    
  end subroutine plane_extents

  subroutine write_plane_3d_real(ipencil,var,iplane,n,dirname,varname,io_name, &
       opt_decomp)

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    real(mytype), dimension(:,:,:), intent(IN) :: var
    integer, intent(IN) :: iplane !(x-plane=1; y-plane=2; z-plane=3)
    integer, intent(IN) :: n ! which plane to write (global coordinate)
    character(len=*), intent(IN) :: dirname,varname,io_name
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    real(mytype), allocatable, dimension(:,:,:) :: wk2d
    TYPE(DECOMP_INFO) :: decomp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, data_type

    logical :: opened_new, dir_exists
    character(len=:), allocatable :: full_io_name
    integer :: idx
#ifdef ADIOS2
    type(adios2_io) :: io_handle
    type(adios2_variable) :: var_handle
#endif
    
    data_type = real_type

#include "io_write_plane.inc"

    return
  end subroutine write_plane_3d_real


  subroutine write_plane_3d_complex(ipencil,var,iplane,n, &
       dirname,varname,io_name,opt_decomp)

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    complex(mytype), dimension(:,:,:), intent(IN) :: var
    integer, intent(IN) :: iplane !(x-plane=1; y-plane=2; z-plane=3)
    integer, intent(IN) :: n ! which plane to write (global coordinate)
    character(len=*), intent(IN) :: dirname,varname,io_name
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    complex(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    complex(mytype), allocatable, dimension(:,:,:) :: wk2d
    TYPE(DECOMP_INFO) :: decomp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, data_type
    logical :: opened_new, dir_exists
    character(len=:), allocatable :: full_io_name
    integer :: idx
#ifdef ADIOS2
    type(adios2_io) :: io_handle
    type(adios2_variable) :: var_handle
#endif

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

  subroutine coarse_extents(ipencil, icoarse, sizes, subsizes, starts, opt_decomp)

    implicit none
    
    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    integer, intent(IN) :: icoarse !(nstat=1; nvisu=2)
    type(decomp_info), intent(in), optional :: opt_decomp

    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror
    type(decomp_info) :: decomp

    if ((icoarse.lt.0).or.(icoarse.gt.2)) then
       print *, "Error invalid value of icoarse: ", icoarse
       call MPI_ABORT(MPI_COMM_WORLD, -1, ierror)
    endif
    if ((ipencil.lt.1).or.(ipencil.gt.3)) then
       print *, "Error invalid value of ipencil: ", ipencil
       call MPI_ABORT(MPI_COMM_WORLD, -1, ierror)
    endif

    if (icoarse==0) then
       !! Use full fields
       
       if (present(opt_decomp)) then
          decomp = opt_decomp
       else
          call get_decomp_info(decomp)
       endif

       sizes(1) = decomp%xsz(1)
       sizes(2) = decomp%ysz(2)
       sizes(3) = decomp%zsz(3)

       if (ipencil == 1) then
          subsizes(1:3) = decomp%xsz(1:3)
          starts(1:3) = decomp%xst(1:3) - 1
       elseif (ipencil == 2) then
          subsizes(1:3) = decomp%ysz(1:3)
          starts(1:3) = decomp%yst(1:3) - 1
       elseif (ipencil == 3) then
          subsizes(1:3) = decomp%zsz(1:3)
          starts(1:3) = decomp%zst(1:3) - 1
       endif
    elseif (icoarse==1) then
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
    elseif (icoarse==2) then
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
    
  end subroutine coarse_extents

  subroutine mpiio_write_real_coarse(ipencil,var,dirname,varname,icoarse,io_name,opt_decomp,reduce_prec,opt_deferred_writes)

    ! USE param
    ! USE variables

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    integer, intent(IN) :: icoarse !(nstat=1; nvisu=2)
    real(mytype), dimension(:,:,:), intent(IN) :: var
    character(len=*), intent(in) :: dirname, varname, io_name
    type(decomp_info), intent(in), optional :: opt_decomp
    logical, intent(in), optional :: reduce_prec
    logical, intent(in), optional :: opt_deferred_writes

    real(mytype_single), allocatable, dimension(:,:,:) :: varsingle
    real(mytype), allocatable, dimension(:,:,:) :: varfull
    logical :: write_reduce_prec
    logical :: deferred_writes
    
    integer (kind=MPI_OFFSET_KIND) :: filesize
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype
    integer :: idx
    logical :: opened_new
    integer :: disp_bytes
    logical :: dir_exists
    character(len=:), allocatable :: full_io_name
#ifdef ADIOS2
    type(adios2_io) :: io_handle
    type(adios2_variable) :: var_handle
    integer :: write_mode
#endif

    !! Set defaults
    write_reduce_prec = .true.
    deferred_writes = .true.
    
    opened_new = .false.
    idx = get_io_idx(io_name, dirname)
#ifndef ADIOS2
    if (present(reduce_prec)) then
       if (.not. reduce_prec) then
          write_reduce_prec = .false.
       end if
    end if
    if (write_reduce_prec) then
       call MPI_TYPE_SIZE(real_type_single,disp_bytes,ierror)
    else
       call MPI_TYPE_SIZE(real_type,disp_bytes,ierror)
    end if
    
    !! Use original MPIIO writers
    if (present(opt_decomp)) then
       call coarse_extents(ipencil, icoarse, sizes, subsizes, starts, opt_decomp)
    else
       call coarse_extents(ipencil, icoarse, sizes, subsizes, starts)
    end if
    if (write_reduce_prec) then
       allocate (varsingle(xstV(1):xenV(1),xstV(2):xenV(2),xstV(3):xenV(3)))
       varsingle=real(var, mytype_single)
    end if

    if (write_reduce_prec) then
       call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
            MPI_ORDER_FORTRAN, real_type_single, newtype, ierror)
    else
       call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
            MPI_ORDER_FORTRAN, real_type, newtype, ierror)
    end if
    call MPI_TYPE_COMMIT(newtype,ierror)

    if (idx .lt. 1) then
       ! Create folder if needed
       if (nrank==0) then
          inquire(file=dirname, exist=dir_exists)
          if (.not.dir_exists) then
             call system("mkdir "//dirname//" 2> /dev/null")
          end if
       end if
       allocate(character(len(trim(dirname)) + 1 + len(trim(varname))) :: full_io_name)
       full_io_name = dirname//"/"//varname
       call decomp_2d_open_io(io_name, full_io_name, decomp_2d_write_mode)
       idx = get_io_idx(io_name, full_io_name)
       opened_new = .true.
    end if

    if (write_reduce_prec) then
       call MPI_FILE_SET_VIEW(fh_registry(idx),fh_disp(idx),real_type_single, &
            newtype,'native',MPI_INFO_NULL,ierror)
       call MPI_FILE_WRITE_ALL(fh_registry(idx), varsingle, &
            subsizes(1)*subsizes(2)*subsizes(3), &
            real_type_single, MPI_STATUS_IGNORE, ierror)
    else
       call MPI_FILE_SET_VIEW(fh_registry(idx),fh_disp(idx),real_type, &
            newtype,'native',MPI_INFO_NULL,ierror)
       call MPI_FILE_WRITE_ALL(fh_registry(idx), var, &
            subsizes(1)*subsizes(2)*subsizes(3), &
            real_type, MPI_STATUS_IGNORE, ierror)
    end if
    
    fh_disp(idx) = fh_disp(idx) + sizes(1) * sizes(2) * sizes(3) * disp_bytes
    
    if (opened_new) then
       call decomp_2d_close_io(io_name, full_io_name)
       deallocate(full_io_name)
    end if

    call MPI_TYPE_FREE(newtype,ierror)
    if (write_reduce_prec) then
       deallocate(varsingle)
    end if
#else
    if (.not. engine_live(idx)) then
       print *, "ERROR: Engine is not live!"
       stop
    end if
    
    call adios2_at_io(io_handle, adios, io_name, ierror)
    call adios2_inquire_variable(var_handle, io_handle, varname, ierror)
    if (.not.var_handle % valid) then
       print *, "ERROR: trying to write variable before registering!", varname
       stop
    endif

    if (idx .lt. 1) then
       print *, "You haven't opened ", io_name, ":", dirname
       stop
    end if

    if (present(opt_deferred_writes)) then
       deferred_writes = opt_deferred_writes
    end if

    if (deferred_writes) then
       write_mode = adios2_mode_deferred
    else
       write_mode = adios2_mode_sync
    end if

    if (engine_registry(idx)%valid) then
       call adios2_put(engine_registry(idx), var_handle, var, write_mode, ierror)
       if (ierror /= 0) then
          print *, "ERROR: something went wrong in adios2_put"
          stop
       end if
    else
       print *, "ERROR: decomp2d thinks engine is live, but adios2 engine object is not valid"
       stop
    end if
#endif

    return
  end subroutine mpiio_write_real_coarse

  subroutine decomp_2d_register_variable(io_name, varname, ipencil, icoarse, iplane, type, opt_decomp, opt_nplanes)

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    integer, intent(IN) :: icoarse !(nstat=1; nvisu=2)
    character(len=*), intent(in) :: io_name
    integer, intent(in) :: type
    integer, intent(in) :: iplane
    type(decomp_info), intent(in), optional :: opt_decomp
    integer, intent(in), optional :: opt_nplanes

    integer :: nplanes
    character*(*), intent(in) :: varname
#ifdef ADIOS2
    integer, dimension(3) :: sizes, subsizes, starts
    type(adios2_io) :: io_handle
    type(adios2_variable) :: var_handle
    integer, parameter :: ndims = 3
    logical, parameter :: adios2_constant_dims = .true.
    integer :: data_type
    integer :: ierror

    if (iplane .eq. 0) then
       if (present(opt_decomp)) then
          call coarse_extents(ipencil, icoarse, sizes, subsizes, starts, opt_decomp)
       else
          call coarse_extents(ipencil, icoarse, sizes, subsizes, starts)
       endif
    else
       if (present(opt_nplanes)) then
          nplanes = opt_nplanes
       else
          nplanes = 1
       end if
       if (present(opt_decomp)) then
          call plane_extents(sizes, subsizes, starts, iplane, opt_decomp, opt_nplanes=nplanes)
       else
          call plane_extents(sizes, subsizes, starts, iplane, opt_nplanes=nplanes)
       endif
    end if
    
    ! Check if variable already exists, if not create it
    call adios2_at_io(io_handle, adios, io_name, ierror)
    if (io_handle%valid) then
       call adios2_inquire_variable(var_handle, io_handle, varname, ierror)
       if (.not.var_handle % valid) then
          !! New variable
          if (nrank .eq. 0) then
             print *, "Registering variable for IO: ", varname
          endif

          ! Need to set the ADIOS2 data type
          if (type.eq.kind(0.0d0)) then
             !! Double
             data_type = adios2_type_dp
          else if (type.eq.kind(0.0)) then
             !! Single
             data_type = adios2_type_real
          else
             print *, "Trying to write unknown data type!"
             call MPI_ABORT(MPI_COMM_WORLD, -1, ierror)
          endif

          call adios2_define_variable(var_handle, io_handle, varname, data_type, &
               ndims, int(sizes, kind=8), int(starts, kind=8), int(subsizes, kind=8), &
               adios2_constant_dims, ierror)
          if (ierror /= 0) then
             print *, "ERROR registering variable"
             stop
          end if
       endif
    else
       print *, "ERROR trying to register variable with invalid IO!"
       stop
    end if
#endif
    
  end subroutine decomp_2d_register_variable
  
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

  subroutine decomp_2d_init_io(io_name)

    implicit none

    character(len=*), intent(in) :: io_name
    integer :: ierror
#ifdef ADIOS2
    type(adios2_io) :: io
#endif
    
    if (nrank .eq. 0) then
       print *, "Initialising IO for ", io_name
    end if

#ifdef ADIOS2
    if (adios%valid) then
       call adios2_declare_io(io, adios, io_name, ierror)
    else
       print *, "ERROR: couldn't declare IO - adios object not valid"
       stop
    end if
    
    if (ierror /= 0) then
       print *, "ERROR: Initialising IO"
       stop
    end if
#endif
    
  end subroutine decomp_2d_init_io
  
  subroutine decomp_2d_open_io(io_name, io_dir, mode)

    implicit none

    logical :: dir_exists
    character(len=*), intent(in) :: io_name, io_dir
    integer, intent(in) :: mode

    logical, dimension(:), pointer :: live_ptrh
    character(len=80), dimension(:), pointer :: names_ptr
    character(len=(len(io_name)+len(io_sep)+len(io_dir))) :: full_name
    
    integer :: idx, ierror
    integer :: access_mode
#ifndef ADIOS2
    integer(MPI_OFFSET_KIND) :: filesize
#else
    type(adios2_io) :: io
#endif

#ifndef ADIOS2
    live_ptrh => fh_live
    names_ptr => fh_names
#else
    live_ptrh => engine_live
    names_ptr => engine_names
#endif
    
    idx = get_io_idx(io_name, io_dir)
    if (idx .lt. 1) then
       !! New io destination
       if (nreg_io .lt. MAX_IOH) then
          nreg_io = nreg_io + 1
          do idx = 1, MAX_IOH
             if (.not. live_ptrh(idx)) then
                live_ptrh(idx) = .true.
                exit
             end if
          end do
          
          full_name = io_name//io_sep//io_dir
          names_ptr(idx) = full_name

          if (mode .eq. decomp_2d_write_mode) then
             !! Setup writers
#ifndef ADIOS2
             filesize = 0_MPI_OFFSET_KIND
             fh_disp(idx) = 0_MPI_OFFSET_KIND
             access_mode = MPI_MODE_CREATE + MPI_MODE_WRONLY
#else
             access_mode = adios2_mode_write
#endif
          else if (mode .eq. decomp_2d_read_mode) then
             !! Setup readers
#ifndef ADIOS2
             fh_disp(idx) = 0_MPI_OFFSET_KIND
             access_mode = MPI_MODE_RDONLY
#else
             access_mode = adios2_mode_read
#endif
          else if (mode .eq. decomp_2d_append_mode) then
#ifndef ADIOS2
             filesize = 0_MPI_OFFSET_KIND
             fh_disp(idx) = 0_MPI_OFFSET_KIND
             access_mode = MPI_MODE_CREATE + MPI_MODE_WRONLY
#else
             access_mode = adios2_mode_append
#endif
          else
             print *, "ERROR: Unknown mode!"
             stop
          endif

          !! Open IO
#ifndef ADIOS2
          call MPI_FILE_OPEN(MPI_COMM_WORLD, io_dir, &
               access_mode, MPI_INFO_NULL, &
               fh_registry(idx), ierror)
          if (mode .eq. decomp_2d_write_mode) then
             !! Guarantee overwriting
             call MPI_FILE_SET_SIZE(fh_registry(idx), filesize, ierror)
          end if
#else
          call adios2_at_io(io, adios, io_name, ierror)
          if (io%valid) then
             call adios2_open(engine_registry(idx), io, trim(gen_iodir_name(io_dir, io_name)), access_mode, ierror)
             if (ierror /= 0) then
                print *, "ERROR opening engine!"
                stop
             end if
          else
             print *, "ERROR: Couldn't find IO handle"
             stop
          end if
#endif
       end if
    end if
    
  end subroutine decomp_2d_open_io

  subroutine decomp_2d_close_io(io_name, io_dir)

    implicit none

    character(len=*), intent(in) :: io_name, io_dir
    
    character(len=80), dimension(:), pointer :: names_ptr
    logical, dimension(:), pointer :: live_ptrh
    integer :: idx, ierror

    idx = get_io_idx(io_name, io_dir)
#ifndef ADIOS2
    names_ptr => fh_names
    live_ptrh => fh_live
    call MPI_FILE_CLOSE(fh_registry(idx), ierror)
#else
    names_ptr => engine_names
    live_ptrh => engine_live
    call adios2_close(engine_registry(idx), ierror)
    if (ierror /= 0) then
       print *, "ERROR closing IO"
    end if
#endif
    names_ptr(idx) = ""
    live_ptrh(idx) = .false.
    nreg_io = nreg_io - 1

  end subroutine decomp_2d_close_io

  subroutine decomp_2d_start_io(io_name, io_dir)

    implicit none

    character(len=*), intent(in) :: io_name, io_dir
#ifdef ADIOS2
    integer :: idx, ierror

    idx = get_io_idx(io_name, io_dir)
    associate(engine => engine_registry(idx))
      if (engine%valid) then
         call adios2_begin_step(engine, ierror)
         if (ierror /= 0) then
            print *, "ERROR beginning step"
            stop
         end if
      else
         print *, "ERROR trying to begin step with invalid engine"
         stop
      end if
    end associate
#endif
    
  end subroutine decomp_2d_start_io

  subroutine decomp_2d_end_io(io_name, io_dir)

    implicit none

    character(len=*), intent(in) :: io_name, io_dir
#ifdef ADIOS2
    integer :: idx, ierror

    idx = get_io_idx(io_name, io_dir)
    associate(engine => engine_registry(idx))
      if (engine%valid) then
         call adios2_end_step(engine, ierror)
         if (ierror /= 0) then
            print *, "ERROR ending step"
            stop
         end if
      else
         print *, "ERROR trying to end step with invalid engine"
         stop
      end if
    end associate
#endif

  end subroutine decomp_2d_end_io
  
  integer function get_io_idx(io_name, engine_name)

    implicit none

    character(len=*), intent(in) :: io_name
    character(len=*), intent(in) :: engine_name

    character(len=(len(io_name)+len(io_sep)+len(engine_name))) :: full_name
    integer :: idx
    logical :: found

    character(len=80), dimension(:), pointer :: names_ptr

#ifndef ADIOS2
    names_ptr => fh_names
#else
    names_ptr => engine_names
#endif

    full_name = io_name//io_sep//engine_name
    
    found = .false.
    do idx = 1, MAX_IOH
       if (names_ptr(idx) .eq. full_name) then
          found = .true.
          exit
       end if
    end do

    if (.not. found) then
       idx = -1
    end if

    get_io_idx = idx
    
  end function get_io_idx

  function gen_iodir_name(io_dir, io_name)

    character(len=*), intent(in) :: io_dir, io_name
    character(len=(len(io_dir) + 5)) :: gen_iodir_name
#ifdef ADIOS2
    integer :: ierror
    type(adios2_io) :: io
    character(len=5) :: ext
#endif

#ifndef ADIOS2
    write(gen_iodir_name, "(A)") io_dir
#else
    call adios2_at_io(io, adios, io_name, ierror)
    if (io%engine_type .eq. "BP4") then
       ext = ".bp4"
    else if (io%engine_type .eq. "HDF5") then
       ext = ".hdf5"
    else
       print *, "ERROR: Unkown engine type! ", io%engine_type
       print *, "-  IO: ", io_name
       print *, "- DIR:", io_dir
       stop
    endif
    write(gen_iodir_name, "(A,A)") io_dir, trim(ext)
#endif
    
  end function gen_iodir_name

end module decomp_2d_io
