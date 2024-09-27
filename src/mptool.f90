!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause
module mptool
  !
  use mpi
  use decomp_2d_constants, only : mytype, real_type
  use decomp_2d_mpi, only: nrank, nproc
  !
  implicit none
  !
  interface psum
    module procedure psum_mytype_ary
    module procedure psum_integer
    module procedure psum_mytype
  end interface
  !
  interface pmax
    module procedure pmax_int
    module procedure pmax_mytype
  end interface
  !
  interface pwrite
    module procedure pwrite_1darray
  end interface
  !
  interface pread
    module procedure pread_1darray
  end interface
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to compute sum in MPI_COMM_WORLD
  !|    It can not use a dedicated MPI communicator
  !+-------------------------------------------------------------------+
  !
  function psum_mytype_ary(var) result(varsum)
    !
    ! arguments
    real(mytype),intent(in) :: var(:)
    real(mytype) :: varsum(size(var))
    !
    ! local data
    integer :: ierr
    !
    call mpi_allreduce(var,varsum,size(var),real_type,mpi_sum,             &
                                                    mpi_comm_world,ierr)

  end function psum_mytype_ary
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to compute sum in parallel
  !|    It can use a dedicated MPI communicator if provided
  !+-------------------------------------------------------------------+
  !
  function psum_integer(var,comm) result(varsum)
    !
    ! arguments
    integer,intent(in) :: var
    integer,optional,intent(in) :: comm
    integer :: varsum
    !
    ! local data
    integer :: ierr,comm2use
    !
    if(present(comm)) then
        comm2use=comm
    else
        comm2use=mpi_comm_world
    endif
    !
    !
    call mpi_allreduce(var,varsum,1,mpi_integer,mpi_sum,           &
                                                    comm2use,ierr)
    !
    return
    !
  end function psum_integer
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to compute sum in parallel
  !|    It can use a dedicated MPI communicator if provided
  !+-------------------------------------------------------------------+
  !
  function psum_mytype(var,comm) result(varsum)
    !
    ! arguments
    real(mytype),intent(in) :: var
    integer,optional,intent(in) :: comm
    real(mytype) :: varsum
    !
    ! local data
    integer :: ierr,comm2use
    !
    if(present(comm)) then
        comm2use=comm
    else
        comm2use=mpi_comm_world
    endif
    !
    !
    call mpi_allreduce(var,varsum,1,real_type,mpi_sum,           &
                                                    comm2use,ierr)
    !
  end function psum_mytype
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to compute the max in MPI_COMM_WORLD
  !|    It can not use a dedicated MPI communicator
  !+-------------------------------------------------------------------+
  !
  integer function  pmax_int(var)
    !
    ! arguments
    integer,intent(in) :: var
    !
    ! local data
    integer :: ierr
    !
    call mpi_allreduce(var,pmax_int,1,mpi_integer,mpi_max,             &
                                                    mpi_comm_world,ierr)
    !
  end function pmax_int
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to compute the max in MPI_COMM_WORLD
  !|    It can not use a dedicated MPI communicator
  !+-------------------------------------------------------------------+
  !
  real(mytype) function  pmax_mytype(var)
    !
    ! arguments
    real(mytype),intent(in) :: var
    !
    ! local data
    integer :: ierr
    !
    call mpi_allreduce(var,pmax_mytype,1,real_type,mpi_max,             &
                                                    mpi_comm_world,ierr)
    !
  end function pmax_mytype
  !
  pure function cross_product(a,b)
    !
    real(mytype) :: cross_product(3)
    real(mytype),dimension(3),intent(in) :: a(3), b(3)
  
    cross_product(1) = a(2) * b(3) - a(3) * b(2)
    cross_product(2) = a(3) * b(1) - a(1) * b(3)
    cross_product(3) = a(1) * b(2) - a(2) * b(1)
    !
  end function cross_product

  subroutine pwrite_1darray(filename,data2write)

    ! arguments
    character(len=*),intent(in) :: filename
    real(mytype),intent(in) :: data2write(:)

    ! local data
    integer :: local_size
    integer :: ierr, fh, datatype, status(mpi_status_size)
    integer(kind=mpi_offset_kind) :: offset

    ! get the size of the array
    local_size=size(data2write)

    ! calculate the offset for each process
    offset = prelay(local_size)*8_8

    ! Open the file in write mode with MPI-IO
    call mpi_file_open(mpi_comm_world, filename, mpi_mode_wronly + mpi_mode_create, mpi_info_null, fh, ierr)

    ! set the file view for each process
    call mpi_file_set_view(fh, offset, real_type, real_type, 'native', mpi_info_null, ierr)

    ! write the local array to the file
    call mpi_file_write(fh, data2write, local_size, real_type, status, ierr)

    ! close the file
    call mpi_file_close(fh, ierr)

    if(nrank==0) print*,'<< ',filename


  end subroutine pwrite_1darray


  subroutine pread_1darray(filename,data2read,size_data2read)

    ! arguments
    character(len=*),intent(in) :: filename
    real(mytype),intent(out),allocatable :: data2read(:)
    integer,intent(in) :: size_data2read

    ! local data
    integer :: ierr, fh, datatype, status(mpi_status_size)
    integer(kind=mpi_offset_kind) :: offset

    ! allocate local array for each process
    allocate(data2read(size_data2read))

    ! calculate the offset for each process
    offset = prelay(size_data2read)*8_8

    ! Open the file in read mode using MPI-IO
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
  
    ! Set the file view for each process
    call MPI_FILE_SET_VIEW(fh, offset, real_type, real_type, 'native', MPI_INFO_NULL, ierr)
  
    ! Read the local portion of the array from the file
    call MPI_FILE_READ(fh, data2read, size_data2read, real_type, status, ierr)
  
    ! Close the file
    call MPI_FILE_CLOSE(fh, ierr)

    if(nrank==0) print*,'>> ',filename


  end subroutine pread_1darray

  ! this function return the number add from all processors before it
  integer function prelay(number)

    ! arguments
    integer,intent(in) :: number

    ! local data
    integer :: table(nproc)
    integer :: ierr,i

    call mpi_allgather(number,1,mpi_integer,                              &
                       table,1,mpi_integer,mpi_comm_world,ierr)

    ! print*,nrank,'=',table
    prelay=0
    do i=1,nrank
      prelay=prelay+table(i)
    enddo
    ! print*,prelay

    return

  end function prelay
  !
end module mptool
