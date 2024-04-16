!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause
module mptool
  !
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
  interface ptabupd
    module procedure ptable_update_int_arr
    module procedure updatable_int
  end interface ptabupd
  !
  !
  interface pgather
    module procedure pgather_int
  end interface
  !
  ! integer(hid_t) :: h5file_id
  !
  integer :: mpi_comm_particle,mpi_rank_part,mpi_size_part
  !
  character(len=4) :: rankname
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to finalise mpi and stop the program.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 19-July-2019: Created by J. Fang @ STFC Daresbury Laboratory      |
  !+-------------------------------------------------------------------+
  subroutine mpistop
    !
    use mpi
    !
    integer :: ierr
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    call mpi_finalize(ierr)
    !
    if(nrank==0) print*,' ** The job is done!'
    !
    stop
    !
  end subroutine mpistop
  !+-------------------------------------------------------------------+
  !| The end of the subroutine mpistop.                                |
  !+-------------------------------------------------------------------+
  !!
  function psum_mytype_ary(var) result(varsum)
    !
    use mpi
    !
    ! arguments
    real(mytype),intent(in) :: var(:)
    real(mytype),allocatable :: varsum(:)
    !
    ! local data
    integer :: ierr,nsize
    !
    nsize=size(var)
    !
    allocate(varsum(nsize))
    !
    call mpi_allreduce(var,varsum,nsize,real_type,mpi_sum,             &
                                                    mpi_comm_world,ierr)
    !
    return
    !
  end function psum_mytype_ary
  !
  function psum_integer(var,comm) result(varsum)
    !
    use mpi
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
  function psum_mytype(var,comm) result(varsum)
    !
    use mpi
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
    return
    !
  end function psum_mytype
  !!
  !+-------------------------------------------------------------------+
  !| this function is to update table based on alltoall mpi            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-Jun-2022  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  function ptable_update_int_arr(vain) result(vout)
    !
    use mpi
    !
    integer,intent(in) :: vain(:)
    integer :: vout(size(vain))
    !
    ! local variables
    integer :: nvar,ierr
    !
    nvar=size(vain)
    !
    call mpi_alltoall(vain,1,mpi_integer,                   &
                      vout,1,mpi_integer,mpi_comm_world,ierr)
    !
    return
    !
  end function ptable_update_int_arr
  !
  function updatable_int(var,offset,debug,comm,comm_size) result(table)
    !
    use mpi
    !
    ! arguments
    integer,allocatable :: table(:)
    integer,intent(in) :: var
    integer,optional,intent(out) :: offset
    logical,intent(in),optional :: debug
    integer,intent(in),optional :: comm,comm_size
    !
    ! local data
    integer :: comm2use,comm2size
    integer :: ierr,i
    integer,allocatable :: vta(:)
    logical :: ldebug
    !
    if(present(debug)) then
      ldebug=debug
    else
      ldebug=.false.
    endif
    !
    if(present(comm)) then
        comm2use=comm
    else
        comm2use=mpi_comm_world
    endif
    !
    if(present(comm_size)) then
        comm2size=comm_size
    else
        comm2size=nproc
    endif
    !
    allocate(table(0:comm2size-1),vta(0:comm2size-1))
    !
    call mpi_allgather(var,1,mpi_integer,                              &
                       vta,1,mpi_integer,comm2use,ierr)
    !
    table=vta
    !
    if(present(offset)) then
      !
      if(nrank==0) then
        offset=0
      else
        !
        offset=0
        do i=0,nrank-1
          offset=offset+vta(i)
        enddo
        !
      endif
      !
    endif
    !
  end function updatable_int
  !
  !+-------------------------------------------------------------------+
  !| The end of the subroutine ptable_update_int_arr.                  |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| The wraper of MPI_Wtime                                           |
  !+-------------------------f------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-November-2019: Created by J. Fang @ STFC Daresbury Laboratory  |
  !+-------------------------------------------------------------------+
  real(mytype) function ptime()
    !
    use mpi
    !
    ptime=MPI_Wtime()
    !
    return
    !
  end function ptime
  !
  subroutine pgather_int(var,data,mode)
    !
    use mpi
    !
    ! arguments
    integer,intent(in) :: var
    integer,intent(out),allocatable :: data(:)
    character(len=*),intent(in),optional :: mode
    !
    !
    ! local data
    integer :: counts(0:nproc-1)
    integer :: ierr,jrank,ncou
    !
    call mpi_allgather(var, 1, mpi_integer, counts, 1, mpi_integer,  &
                       mpi_comm_world, ierr)
    !
    if(present(mode) .and. mode=='noneg') then
      ! only pick >=0 values
      ncou=0
      do jrank=0,nproc-1
        if(counts(jrank)>=0) then
          ncou=ncou+1
        endif
      enddo
      !
      allocate(data(ncou))
      ncou=0
      do jrank=0,nproc-1
        if(counts(jrank)>=0) then
          ncou=ncou+1
          data(ncou)=counts(jrank)
        endif
      enddo
      !
    else
      allocate(data(0:nproc-1))
      data=counts
    endif
    !
  end subroutine pgather_int
  !
  integer function  pmax_int(var)
    !
    use mpi
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
  real(mytype) function  pmax_mytype(var)
    !
    use mpi
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
    real(8) :: cross_product(3)
    real(8),dimension(3),intent(in) :: a(3), b(3)
  
    cross_product(1) = a(2) * b(3) - a(3) * b(2)
    cross_product(2) = a(3) * b(1) - a(1) * b(3)
    cross_product(3) = a(1) * b(2) - a(2) * b(1)
    !
    return
    !
  end function cross_product
  !
end module mptool
