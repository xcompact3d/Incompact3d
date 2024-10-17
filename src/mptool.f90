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
  interface ptabupd
    module procedure ptable_update_int_arr
    module procedure updatable_int
  end interface ptabupd
  !
  interface pwrite
    module procedure pwrite_1darray
    module procedure pwrite_2darray
  end interface
  !
  interface pread
    module procedure pread_1darray
    module procedure pread_2darray
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
  !| The end of the subroutine ptabupd.                                |
  !+-------------------------------------------------------------------+

  pure function cross_product(a,b)
    !
    real(mytype) :: cross_product(3)
    real(mytype),dimension(3),intent(in) :: a(3), b(3)
  
    cross_product(1) = a(2) * b(3) - a(3) * b(2)
    cross_product(2) = a(3) * b(1) - a(1) * b(3)
    cross_product(3) = a(1) * b(2) - a(2) * b(1)
    !
  end function cross_product

  subroutine pwrite_1darray(filename,data_1,data_2,data_3,data_4,data_5,data_6)

    ! arguments
    character(len=*),intent(in) :: filename
    real(mytype),intent(in) :: data_1(:)
    real(mytype),intent(in),optional :: data_2(:),data_3(:),data_4(:),data_5(:),data_6(:)

    ! local data
    integer :: local_size,total_size
    integer :: ierr, fh, datatype, status(mpi_status_size)
    integer(kind=mpi_offset_kind) :: offset

    ! get the size of the array
    local_size=size(data_1)

    ! calculate the offset for each process
    offset = prelay(local_size)*8_8

    ! Open the file in write mode with MPI-IO
    call mpi_file_open(mpi_comm_world, filename, mpi_mode_wronly + mpi_mode_create, mpi_info_null, fh, ierr)

    ! set the file view for each process
    call mpi_file_set_view(fh, offset, real_type, real_type, 'native', mpi_info_null, ierr)

    ! write the local array to the file
    call mpi_file_write(fh, data_1, local_size, real_type, status, ierr)

    if(present(data_2)) then

      ! Barrier to synchronize processes
      call mpi_barrier(mpi_comm_world, ierr)

      total_size=psum(local_size)*8_8

      offset = offset + total_size

      ! write the data_2 to the file
      call mpi_file_set_view(fh, offset, real_type, real_type, 'native', mpi_info_null, ierr)
      call mpi_file_write(fh, data_2, local_size, real_type, status, ierr)

    endif

    if(present(data_3)) then

      ! Barrier to synchronize processes
      call mpi_barrier(mpi_comm_world, ierr)

      offset = offset + total_size

      ! write the data_3 to the file
      call mpi_file_set_view(fh, offset, real_type, real_type, 'native', mpi_info_null, ierr)
      call mpi_file_write(fh, data_3, local_size, real_type, status, ierr)

    endif

    if(present(data_4)) then

      ! Barrier to synchronize processes
      call mpi_barrier(mpi_comm_world, ierr)

      offset = offset + total_size

      ! write the data_4 to the file
      call mpi_file_set_view(fh, offset, real_type, real_type, 'native', mpi_info_null, ierr)
      call mpi_file_write(fh, data_4, local_size, real_type, status, ierr)

    endif

    if(present(data_5)) then

      ! Barrier to synchronize processes
      call mpi_barrier(mpi_comm_world, ierr)

      offset = offset + total_size

      ! write the data_5 to the file
      call mpi_file_set_view(fh, offset, real_type, real_type, 'native', mpi_info_null, ierr)
      call mpi_file_write(fh, data_5, local_size, real_type, status, ierr)

    endif

    if(present(data_6)) then

      ! Barrier to synchronize processes
      call mpi_barrier(mpi_comm_world, ierr)

      offset = offset + total_size

      ! write the data_6 to the file
      call mpi_file_set_view(fh, offset, real_type, real_type, 'native', mpi_info_null, ierr)
      call mpi_file_write(fh, data_6, local_size, real_type, status, ierr)

    endif

    ! close the file
    call mpi_file_close(fh, ierr)

    if(nrank==0) print*,'<< ',filename

  end subroutine pwrite_1darray

  subroutine pwrite_2darray(filename,data_1,data_2,data_3,data_4,data_5,data_6)

    ! arguments
    character(len=*),intent(in) :: filename
    real(mytype),intent(in) :: data_1(:,:)
    real(mytype),intent(in),optional :: data_2(:,:),data_3(:,:),data_4(:,:), &
                                        data_5(:,:),data_6(:,:)

    ! local data
    integer :: local_size,total_size
    integer :: ierr, fh, datatype, status(mpi_status_size)
    integer(kind=mpi_offset_kind) :: offset

    ! get the size of the array
    local_size=size(data_1)

    ! calculate the offset for each process
    offset = prelay(local_size)*8_8

    ! Open the file in write mode with MPI-IO
    call mpi_file_open(mpi_comm_world, filename, mpi_mode_wronly + mpi_mode_create, mpi_info_null, fh, ierr)

    ! set the file view for each process
    call mpi_file_set_view(fh, offset, real_type, real_type, 'native', mpi_info_null, ierr)

    ! write the local array to the file
    call mpi_file_write(fh, data_1, local_size, real_type, status, ierr)

    if(present(data_2)) then

      ! Barrier to synchronize processes
      call mpi_barrier(mpi_comm_world, ierr)

      total_size=psum(local_size)*8_8

      offset = offset + total_size

      ! write the data_2 to the file
      call mpi_file_set_view(fh, offset, real_type, real_type, 'native', mpi_info_null, ierr)
      call mpi_file_write(fh, data_2, local_size, real_type, status, ierr)

    endif

    if(present(data_3)) then

      ! Barrier to synchronize processes
      call mpi_barrier(mpi_comm_world, ierr)

      offset = offset + total_size

      ! write the data_3 to the file
      call mpi_file_set_view(fh, offset, real_type, real_type, 'native', mpi_info_null, ierr)
      call mpi_file_write(fh, data_3, local_size, real_type, status, ierr)

    endif

    if(present(data_4)) then

      ! Barrier to synchronize processes
      call mpi_barrier(mpi_comm_world, ierr)

      offset = offset + total_size

      ! write the data_4 to the file
      call mpi_file_set_view(fh, offset, real_type, real_type, 'native', mpi_info_null, ierr)
      call mpi_file_write(fh, data_4, local_size, real_type, status, ierr)

    endif

    if(present(data_5)) then

      ! Barrier to synchronize processes
      call mpi_barrier(mpi_comm_world, ierr)

      offset = offset + total_size

      ! write the data_5 to the file
      call mpi_file_set_view(fh, offset, real_type, real_type, 'native', mpi_info_null, ierr)
      call mpi_file_write(fh, data_5, local_size, real_type, status, ierr)

    endif

    if(present(data_6)) then

      ! Barrier to synchronize processes
      call mpi_barrier(mpi_comm_world, ierr)

      offset = offset + total_size

      ! write the data_6 to the file
      call mpi_file_set_view(fh, offset, real_type, real_type, 'native', mpi_info_null, ierr)
      call mpi_file_write(fh, data_6, local_size, real_type, status, ierr)

    endif

    ! close the file
    call mpi_file_close(fh, ierr)

    if(nrank==0) print*,'<< ',filename

  end subroutine pwrite_2darray

  subroutine pread_1darray(filename,data2read)

    ! arguments
    character(len=*),intent(in) :: filename
    real(mytype),intent(inout),allocatable :: data2read(:)

    ! local data
    integer :: local_size, ierr, fh, datatype, status(mpi_status_size)
    integer(kind=mpi_offset_kind) :: offset

    ! calcualte local size
    local_size=size(data2read)

    ! calculate the offset for each process
    offset = prelay(local_size)*8_8

    ! Open the file in read mode using MPI-IO
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
  
    ! Set the file view for each process
    call MPI_FILE_SET_VIEW(fh, offset, real_type, real_type, 'native', MPI_INFO_NULL, ierr)
  
    ! Read the local portion of the array from the file
    call MPI_FILE_READ(fh, data2read, local_size, real_type, status, ierr)
  
    ! Close the file
    call MPI_FILE_CLOSE(fh, ierr)

    if(nrank==0) print*,'>> ',filename


  end subroutine pread_1darray


  subroutine pread_2darray(filename,data2read)

    ! arguments
    character(len=*),intent(in) :: filename
    real(mytype),intent(inout),allocatable :: data2read(:,:)

    ! local data
    integer :: local_size, ierr, fh, datatype, status(mpi_status_size)
    integer(kind=mpi_offset_kind) :: offset


    ! calcualte local size
    local_size=size(data2read)

    ! calculate the offset for each process
    offset = prelay(local_size)*8_8

    ! Open the file in read mode using MPI-IO
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
  
    ! Set the file view for each process
    call MPI_FILE_SET_VIEW(fh, offset, real_type, real_type, 'native', MPI_INFO_NULL, ierr)
  
    ! Read the local portion of the array from the file
    call MPI_FILE_READ(fh, data2read, local_size, real_type, status, ierr)
  
    ! Close the file
    call MPI_FILE_CLOSE(fh, ierr)

    if(nrank==0) print*,'>> ',filename


  end subroutine pread_2darray

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
  function linintp(xx1,xx2,yy1,yy2,xx) result(yy)
    !
    real(8),intent(in) :: xx1,xx2,xx
    real(8),intent(in) ::  yy1(:,:,:),yy2(:,:,:)
    real(8) :: yy(1:size(yy1,1),1:size(yy1,2),1:size(yy1,3))
    !
    real(8) :: var1
    !
    if(abs(xx-xx1)<1.d-16) then
      yy=yy1
    elseif(abs(xx-xx2)<1.d-16) then
      yy=yy2
    else
      var1=(xx-xx1)/(xx2-xx1)
      yy=(yy2-yy1)*var1+yy1
    endif
    !
    return
    !
  end function linintp

  ! this function distritbuion evently in processors (as much as possible..)
  pure function numdist(number)
    
    integer,intent(in) :: number
    integer :: numdist,remainder

    numdist = number/nproc

    remainder=mod(number,nproc)

    if(nrank<remainder) then
        numdist=numdist+1
    endif

    return

  end function numdist

end module mptool
