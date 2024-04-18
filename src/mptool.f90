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
  !
end module mptool
