!!! 2decomp4py.f90 provides a light-weight interface to
!!!                2decomp&fft/Xcompact3d for wrapping with f2py.
!!! Copyright (C) 2021-      University of Edinburgh
!!! Copyright (C) 2011-2021, The Numerical Algorithms Group (NAG) All rights reserved.
!!! 
!!! Redistribution and use in source and binary forms, with or without modification, are permitted
!!! provided that the following conditions are met:
!!! 
!!!     Redistributions of source code must retain the above copyright notice, this list of
!!!     conditions and the following disclaimer.
!!!     Redistributions in binary form must reproduce the above copyright notice, this list of
!!!     conditions and the following disclaimer in the documentation and/or other materials provided
!!!     with the distribution.
!!!     Neither the name of the copyright owner nor the names of its contributors may be used to
!!!     endorse or promote products derived from this software without specific prior written
!!!     permission.
!!! 
!!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
!!! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
!!! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER BE LIABLE
!!! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
!!! OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
!!! STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!!! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

! Due to how f2py handles kind statements (see .f2py_f2cmap), we need a kind specifier that expands
! into *text* constants, not variables/parameters as used in 2decomp&fft.
#ifdef DOUBLE_PREC
#define D2DREAL real(kind(0.0d0))
#else
#define D2DREAL real(kind(0.0))
#endif

module decomp4py
  ! A module to wrap the 2decomp&fft library such that f2py can build a Python extension.
  
  use mpi
  use decomp_2d
  use decomp_2d_io
#ifdef ADIOS2
  use adios2
#endif
  
  implicit none

#ifdef ADIOS2
  logical, public, parameter :: have_adios2 = .true.
#else
  logical, public, parameter :: have_adios2 = .false.
#endif
  
  private
  public :: init_decomp4py
  public :: get_grid_size
  public :: get_grid_start
  public :: transpose
  public :: read_field
  public :: register
  
contains

  subroutine init_decomp4py(nx, ny, nz, p_row, p_col)
    ! Initialises the 2decomp&fft library for a domain discretised by :math:`n_x \times n_y \times
    ! n_z` mesh points with a :math:`p_{row} \times p_{col}` pencil decomposition.
    
    integer, intent(in) :: nx, ny, nz
    integer, intent(in) :: p_row, p_col
    !f2py intent(in) nx
    !f2py intent(in) ny
    !f2py intent(in) nz
    !f2py intent(in) p_row
    !f2py intent(in) col
    integer :: ierr
    integer :: nvisu

    ! XXX: Need to initialise the nproc and nrank variables from 2decomp&fft.
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
    
    call decomp_2d_init(nx, ny, nz, p_row, p_col)
    call decomp_2d_io_init()

    nvisu = 1
    call init_coarser_mesh_statV(nvisu, nvisu, nvisu, .true.)
    
  end subroutine init_decomp4py

  function get_grid_size(ax)

    character(len=*), intent(in) :: ax
    integer, dimension(3) :: get_grid_size
    !f2py intent(in) ax
    !f2py intent(out) get_grid_size
    
    if (ax == "x") then
       get_grid_size = xsize
    else if (ax == "y") then
       get_grid_size = ysize
    else if (ax == "z") then
       get_grid_size = zsize
    else
       print *, "ERROR: Only the x/y/z grid size can be queried."
    end if
    
  end function get_grid_size

  function get_grid_start(ax)

    character(len=*), intent(in) :: ax
    integer, dimension(3) :: get_grid_start
    !f2py intent(in) ax
    !f2py intent(out) get_grid_start
    type(decomp_info) :: decomp

    call get_decomp_info(decomp)

    if (ax == "x") then
       get_grid_start = decomp%xst
    else if (ax == "y") then
       get_grid_start = decomp%yst
    else if (ax == "z") then
       get_grid_start = decomp%zst
    else
       print *, "ERROR: Only the x/y/z grid size can be queried."
    end if

  end function get_grid_start
  
  subroutine transpose(arr_in, op, arr_tp)
    !

    real(kind(0.0d0)), dimension(:,:,:), intent(in) :: arr_in
    character(len=*), intent(in) :: op
    real(kind(0.0d0)), dimension(:,:,:), intent(inout) :: arr_tp
    !f2py intent(in) arr_in
    !f2py intent(in) op
    !f2py intent(inout) arr_tp
    
    if (op == "xy") then
       call transpose_x_to_y(arr_in, arr_tp)
    else if (op == "yz") then
       call transpose_y_to_z(arr_in, arr_tp)
    else if (op == "zy") then
       call transpose_z_to_y(arr_in, arr_tp)
    else if (op == "yx") then
       call transpose_y_to_x(arr_in, arr_tp)
    else
       print *, "ERROR: Only xy, yz, zy and yx transposes are supported in 2decomp&fft!"
    end if
       
  end subroutine transpose

  subroutine read_field(ipencil, var, dirname, varname, io_name, nx, ny, nz)

    integer, intent(in) :: ipencil
    real(kind(0.0d0)), dimension(nx, ny, nz), intent(inout) :: var
    character(len=*), intent(in) :: varname, dirname, io_name
    integer, intent(in) :: nx, ny, nz
    !f2py intent(in) ipencil
    !f2py intent(inout) var
    !f2py intent(in) varname
    !f2py intent(in) dirname
    !f2py intent(in) io_name
    
    call decomp_2d_read_one(ipencil,var,dirname,varname,io_name)
    
  end subroutine read_field

  subroutine register(field_name, io_name)

    use decomp_2d_io, only : decomp_2d_register_variable
    
    character(len=*), intent(in) :: field_name
    character(len=*), intent(in) :: io_name

    integer, parameter :: output2D = 0
    
    call decomp_2d_register_variable(io_name, field_name, 1, 0, output2D, kind(0.0d0))
    
  end subroutine register
  
end module decomp4py
