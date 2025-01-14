!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module x3d_io

  use decomp_2d_constants
  use decomp_2d_io
  use decomp_2d_io_adios
  use decomp_2d_io_object_adios
  use decomp_2d_io_object_mpi
  use decomp_2d_mpi
  use decomp_2d

  implicit none

  !
  ! Derived type for each active IO
  !
  type, public :: x3d_live_io
     character(:), allocatable :: label    ! Label of the active IO
     logical :: mpi, adios
     type(d2d_io_mpi) :: mpi_io
     type(d2d_io_adios) :: adios_io
     type(d2d_io_family) :: family
  contains
     procedure :: init => x3d_live_io_init ! Init the object
     procedure :: fin => x3d_live_io_fin   ! Clear the object
  end type x3d_live_io

  ! This is used when the external code does not provide a live IO object
  type(x3d_live_io), target, save :: x3d_live_io_default

  private
  public :: x3d_io_init, &
            x3d_io_fin, &
            x3d_io_open, &
            x3d_io_close, &
            x3d_io_declare, &
            x3d_io_register_var, &
            x3d_io_register_plane, &
            x3d_io_write, &
            x3d_io_read, &
            x3d_live_io_get_default, &
            gen_iodir_name

  interface x3d_io_write
    module procedure x3d_io_write_real
  end interface x3d_io_write

  interface x3d_io_read
    module procedure x3d_io_read_real
  end interface x3d_io_read

contains

  !
  ! Initialize the x3d IO module
  !
  subroutine x3d_io_init()

    implicit none

    ! Initialize decomp2d IO
    call decomp_2d_io_init()

    ! Initialize the default live IO object
    call x3d_live_io_default%init("x3d_default")

  end subroutine x3d_io_init

  !
  ! Finalize the x3d IO module
  !
  subroutine x3d_io_fin()

    implicit none

    ! Finalize the default live IO object
    call x3d_live_io_default%fin()

    ! Finalize decomp2d IO
    call decomp_2d_io_fin()

  end subroutine x3d_io_fin

  !
  ! MPI IO : open the file
  ! ADIOS2 IO : open + start the reader / writer
  !
  subroutine x3d_io_open(io, fname, mode)

    implicit none

    type(x3d_live_io), intent(inout) :: io
    character(len=*), intent(in) :: fname
    integer, intent(in) :: mode

    if (io%adios) then
      call io%adios_io%open_start(mode, opt_family=io%family)
    else
      call io%mpi_io%open(fname, mode)
    end if

  end subroutine x3d_io_open

  !
  ! MPI IO : close the file
  ! ADIOS2 IO : stop + end the reader / writer
  !
  subroutine x3d_io_close(io)

    implicit none

    type(x3d_live_io), intent(inout) :: io

    if (io%adios) then
      call io%adios_io%end_close()
    else
      call io%mpi_io%close()
    end if

  end subroutine x3d_io_close

  !
  ! MPI IO : nothing to do
  ! ADIOS2 IO : declare a new IO family
  !
  subroutine x3d_io_declare(io, name)

    implicit none

    type(x3d_live_io), intent(inout) :: io
    character(len=*), intent(in) :: name

    ! Init if needed
    if (.not.allocated(io%label)) call io%init(name)

    ! Declare a new IO family
    if (io%adios) call io%family%init(name)

  end subroutine x3d_io_declare

  !
  ! MPI IO : nothing to do
  ! ADIOS2 IO : register the 3D variable (in the given IO family)
  !
  subroutine x3d_io_register_var(io, name, ipencil, type, &
                                 opt_reduce_prec, &
                                 opt_decomp)

    implicit none

    type(x3d_live_io), intent(inout) :: io
    character(len=*), intent(in) :: name
    integer, intent(in) :: ipencil, type
    logical, intent(in), optional :: opt_reduce_prec
    type(decomp_info), intent(in), optional :: opt_decomp

    if (io%adios) call io%family%register_var(name, &
                                              ipencil, &
                                              type, &
                                              opt_reduce_prec, &
                                              opt_decomp)

  end subroutine x3d_io_register_var

  !
  ! MPI IO : nothing to do
  ! ADIOS2 IO : register the planes (in the given IO family)
  !
  subroutine x3d_io_register_plane(io, name, ipencil, type, &
                                   opt_reduce_prec, &
                                   opt_decomp, &
                                   opt_nplanes)

    implicit none

    type(x3d_live_io), intent(inout) :: io
    character(len=*), intent(in) :: name
    integer, intent(in) :: ipencil, type
    logical, intent(in), optional :: opt_reduce_prec
    type(decomp_info), intent(in), optional :: opt_decomp
    integer, intent(in), optional :: opt_nplanes

    if (io%adios) call io%family%register_plane(name, &
                                                ipencil, &
                                                type, &
                                                opt_reduce_prec, &
                                                opt_decomp, &
                                                opt_nplanes)

  end subroutine x3d_io_register_plane

  !
  ! Write a 3D real array
  !
  subroutine x3d_io_write_real(io, ipencil, var, dirname, varname, &
                               opt_reduce_prec, &
                               opt_decomp)

    implicit none

    ! Arguments
    type(x3d_live_io), intent(inout) :: io
    integer, intent(in) :: ipencil
    real(mytype), contiguous, dimension(:, :, :), intent(in) :: var
    character(len=*), intent(in) :: dirname, varname
    logical, intent(in), optional :: opt_reduce_prec
    type(decomp_info), target, intent(in), optional :: opt_decomp

    ! Local variable
    logical :: open_close

    ! Check if open / close is needed
    if (io%adios) then
      open_close = .not.io%adios_io%is_open
    else
      open_close = .not.io%mpi_io%is_open
    end if

    ! Open/start if needed
    if (open_close) then
      if (io%adios) then
        call io%adios_io%open_start(decomp_2d_write_mode, opt_family=io%family)
      else
        call io%mpi_io%open(trim(dirname)//"/"//trim(varname), decomp_2d_write_mode)
      end if
    end if

    ! Write
    if (io%adios) then
      call decomp_2d_adios_write_var(io%adios_io, var, varname, &
                                     opt_mode=decomp_2d_io_sync, &
                                     opt_family=io%family, &
                                     opt_reduce_prec=opt_reduce_prec)
    else
      call decomp_2d_write_var(io%mpi_io, ipencil, var, &
                               opt_reduce_prec=opt_reduce_prec, &
                               opt_decomp=opt_decomp)
    end if

    ! End/close if needed
    if (open_close) then
      if (io%adios) then
        call io%adios_io%end_close()
      else
        call io%mpi_io%close()
      end if
    end if

  end subroutine x3d_io_write_real

  !
  ! Read a 3D real array
  !
  subroutine x3d_io_read_real(io, ipencil, var, dirname, varname, &
                              opt_reduce_prec, &
                              opt_decomp)

    implicit none

    ! Arguments
    type(x3d_live_io), intent(inout) :: io
    integer, intent(in) :: ipencil
    real(mytype), contiguous, dimension(:, :, :), intent(out) :: var
    character(len=*), intent(in) :: dirname, varname
    logical, intent(in), optional :: opt_reduce_prec
    type(decomp_info), target, intent(in), optional :: opt_decomp

    ! Local variable
    logical :: open_close

    ! Check if open / close is needed
    if (io%adios) then
      open_close = .not.io%adios_io%is_open
    else
      open_close = .not.io%mpi_io%is_open
    end if

    ! Open/start if needed
    if (open_close) then
      if (io%adios) then
        call io%adios_io%open_start(decomp_2d_read_mode, opt_family=io%family)
      else
        call io%mpi_io%open(trim(dirname)//"/"//trim(varname), decomp_2d_read_mode)
      end if
    end if

    ! Read
    if (io%adios) then
      call decomp_2d_adios_read_var(io%adios_io, var, varname, &
                                    opt_family=io%family, &
                                    opt_reduce_prec=opt_reduce_prec)
    else
      call decomp_2d_read_var(io%mpi_io, ipencil, var, &
                              opt_reduce_prec=opt_reduce_prec, &
                              opt_decomp=opt_decomp)
    end if

    ! End/close if needed
    if (open_close) then
      if (io%adios) then
        call io%adios_io%end_close()
      else
        call io%mpi_io%close()
      end if
    end if

  end subroutine x3d_io_read_real

  !
  ! Get the default live IO object.
  !
  ! The caller should not finalize it.
  !
  function x3d_live_io_get_default() result(output)

    implicit none

    type(x3d_live_io), pointer :: output

    output => x3d_live_io_default

  end function x3d_live_io_get_default

  !
  ! Initialize the IO object
  !
  subroutine x3d_live_io_init(io, label, mpi, adios)

    implicit none

    ! Arguments
    class(x3d_live_io), intent(inout) :: io
    character(len=*), intent(in) :: label
    logical, intent(in), optional :: mpi, adios

    ! Store the name of the IO object
    if (allocated(io%label)) deallocate(io%label)
    io%label = label

    ! Default : ADIOS2 IO if possible
    if (present(adios)) then
      io%adios = adios
    else
      io%adios = .true.
    end if
    io%adios = io%adios .and. decomp_2d_with_adios2

    ! Default : MPI if no ADIOS2
    if (present(mpi)) then
      io%mpi = mpi
    else
      io%mpi = .not.io%adios
    end if

  end subroutine x3d_live_io_init

  !
  ! Finalize the IO object
  !
  subroutine x3d_live_io_fin(io)

    implicit none

    class(x3d_live_io), intent(inout) :: io

    if (allocated(io%label)) deallocate(io%label)

  end subroutine x3d_live_io_fin

  function gen_iodir_name(io_dir, io_name)

    character(len=*), intent(in) :: io_dir, io_name
    character(len=(len(io_dir) + 5)) :: gen_iodir_name
#ifdef ADIOS2
    integer :: ierror
    type(adios2_io) :: io
    character(len=5) :: ext
#endif

#ifndef ADIOS2
    write (gen_iodir_name, "(A)") io_dir
#else
    call adios2_at_io(io, adios, io_name, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "adios2_at_io "//trim(io_name))
    if (io%engine_type == "BP4") then
      ext = ".bp4"
    else if (io%engine_type == "HDF5") then
      ext = ".hdf5"
    else if (io%engine_type == "SST") then
      ext = ""
    else
      print *, "ERROR: Unkown engine type! ", io%engine_type
      print *, "-  IO: ", io_name
      print *, "- DIR:", io_dir
      stop
    end if
    write (gen_iodir_name, "(A,A)") io_dir, trim(ext)
#endif

  end function gen_iodir_name

end module x3d_io
