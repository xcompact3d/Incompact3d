!+---------------------------------------------------------------------+
!| this module contains some tool subroutines used by mhd and particle |
!| modules                                                             |
!+---------------------------------------------------------------------+
!| change record                                                       |
!| -------------                                                       |
!| 30-Mar-2023  | Created by J. Fang STFC Daresbury Laboratory         |
!+---------------------------------------------------------------------+
module mptool
  !
  use decomp_2d, only : mytype,nrank,nproc
  ! use hdf5
  ! use h5lt
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
  Interface h5write
    !
    module procedure h5wa_r8
    module procedure h5w_real8
    module procedure h5w_int4
    !
  end Interface h5write
  !
  Interface h5sread
    !
    module procedure h5_readarray1d
    module procedure h5_read1rl8
    module procedure h5_read1int
    !
  end Interface h5sread
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
    use decomp_2d, only : real_type
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
    use decomp_2d, only : real_type
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
    use decomp_2d, only : real_type
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
  !+-------------------------------------------------------------------+
  !| The end of the function ptime.                                    |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to open the h5file interface and assign   |
  !| h5file_id. For write each new file, this will be called first, but|
  !| once it is called, the file will be overwriten.                   |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 03-Jun-2020 | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine h5io_init(filename,mode,comm)
    !
    use mpi, only: mpi_comm_world,mpi_info_null
    !
    ! arguments
    character(len=*),intent(in) :: filename
    character(len=*),intent(in) :: mode
    integer,intent(in),optional :: comm
    ! h5file_id is returned
    !
    ! ! local data
    ! integer :: h5error,comm2use
    ! integer(hid_t) :: plist_id
    ! !
    ! if(present(comm)) then
    !     comm2use=comm
    ! else
    !     comm2use=mpi_comm_world
    ! endif
    ! !
    ! call h5open_f(h5error)
    ! if(h5error.ne.0)  stop ' !! error in h5io_init call h5open_f'
    ! !
    ! ! create access property list and set mpi i/o
    ! call h5pcreate_f(h5p_file_access_f,plist_id,h5error)
    ! if(h5error.ne.0)  stop ' !! error in h5io_init call h5pcreate_f'
    ! !
    ! call h5pset_fapl_mpio_f(plist_id,comm2use,mpi_info_null,     &
    !                                                             h5error)
    ! if(h5error.ne.0)  stop ' !! error in h5io_init call h5pset_fapl_mpio_f'
    ! !
    ! if(mode=='writ') then
    !   call h5fcreate_f(filename,h5f_acc_trunc_f,h5file_id,             &
    !                                         h5error,access_prp=plist_id)
    !   if(h5error.ne.0)  stop ' !! error in h5io_init call h5fcreate_f'
    ! elseif(mode=='read') then
    !   call h5fopen_f(filename,h5f_acc_rdwr_f,h5file_id,                &
    !                                         h5error,access_prp=plist_id)
    !   if(h5error.ne.0)  stop ' !! error in h5io_init call h5fopen_f'
    ! else
    !     stop ' !! mode not defined @ h5io_init'
    ! endif
    ! !
    ! call h5pclose_f(plist_id,h5error)
    ! if(h5error.ne.0)  stop ' !! error in h5io_init call h5pclose_f'
    ! !
    ! if(nrank==0) print*,' ** open h5 file: ',filename
    !
  end subroutine h5io_init
  !+-------------------------------------------------------------------+
  !| This end of the subroutine h5io_init.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to close hdf5 interface after finish      |
  !| input/output a hdf5 file.                                         |
  !| the only data needed is h5file_id                                 |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 03-Jun-2020 | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine h5io_end
    !
    ! local data
    integer :: h5error
    !
    ! call h5fclose_f(h5file_id,h5error)
    ! if(h5error.ne.0)  stop ' !! error in h5io_end call h5fclose_f'
    ! !
    ! call h5close_f(h5error)
    ! if(h5error.ne.0)  stop ' !! error in h5io_end call h5close_f'
    !
  end subroutine h5io_end
  !+-------------------------------------------------------------------+
  !| This end of the subroutine h5io_end.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to write a 1D array with hdf5 interface.  |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 02-Jun-2020 | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine h5wa_r8(varname,var,total_size,comm,comm_size,comm_rank)
    !
    use decomp_2d, only : mytype
    use mpi, only: mpi_comm_world,mpi_info_null
    ! use parallel,only: nrank,mpistop,psum,ptabupd,nrankmax
    !
    ! arguments
    character(LEN=*),intent(in) :: varname
    real(mytype),intent(in),allocatable :: var(:)
    integer,intent(in) :: total_size
    integer,intent(in),optional :: comm,comm_size,comm_rank
    !
    ! ! local data
    ! integer :: jrk
    ! integer :: dim,dima,rank2use,comm2use,comm2size
    ! integer,allocatable :: dim_table(:)
    ! integer(hsize_t), dimension(1) :: offset
    ! integer :: h5error
    ! !
    ! integer(hid_t) :: dset_id,filespace,memspace,plist_id
    ! integer(hsize_t) :: dimt(1),dimat(1)
    ! !
    ! if(allocated(var)) then
    !   dim=size(var)
    ! else
    !   dim=0
    ! endif
    ! !
    ! if(present(comm)) then
    !     comm2use=comm
    ! else
    !     comm2use=mpi_comm_world
    ! endif
    ! !
    ! if(present(comm_size)) then
    !     comm2size=comm_size
    ! else
    !     comm2size=nproc
    ! endif
    ! !
    ! if(present(comm_rank)) then
    !     rank2use=comm_rank
    ! else
    !     rank2use=nrank
    ! endif
    ! !
    ! allocate(dim_table(0:comm2size-1))
    ! ! dima=psum(dim,comm=comm2use)
    ! !
    ! dimt=(/dim/)
    ! dimat=(/total_size/)
    ! !
    ! dim_table=ptabupd(dim,comm=comm2use,comm_size=comm2size)
    ! !
    ! offset=0
    ! do jrk=0,rank2use-1
    !   offset=offset+dim_table(jrk)
    ! enddo
    ! !
    ! ! print*,mpi_rank_part,nrank,'|',dim,offset
    ! !
    ! ! writing the data
    ! !
    ! call h5screate_simple_f(1,dimat,filespace,h5error)
    ! if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5screate_simple_f'
    ! call h5dcreate_f(h5file_id,varname,h5t_native_double,filespace,    &
    !                                                    dset_id,h5error)

    ! if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5dcreate_f'
    ! call h5screate_simple_f(1,dimt,memspace,h5error)
    ! if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5screate_simple_f'
    ! call h5sclose_f(filespace,h5error)
    ! if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5sclose_f'
    ! call h5dget_space_f(dset_id,filespace,h5error)
    ! if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5dget_space_f'
    ! call h5sselect_hyperslab_f(filespace,h5s_select_set_f,offset,      &
    !                                                       dimt,h5error)
    ! if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5sselect_hyperslab_f'
    ! call h5pcreate_f(h5p_dataset_xfer_f,plist_id,h5error)
    ! if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5pcreate_f'
    ! call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_collective_f,h5error)
    ! if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5pset_dxpl_mpio_f'
    ! call h5dwrite_f(dset_id,h5t_native_double,var,dimt,h5error,        &
    !                 file_space_id=filespace,mem_space_id=memspace,     &
    !                                                  xfer_prp=plist_id)
    ! if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5dwrite_f'
    ! call h5sclose_f(filespace,h5error)
    ! if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5sclose_f'
    ! call h5sclose_f(memspace,h5error)
    ! if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5sclose_f'
    ! call h5dclose_f(dset_id,h5error)
    ! if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5dclose_f'
    ! call h5pclose_f(plist_id,h5error)
    ! if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5pclose_f'
    ! !
    ! if(mpi_rank_part==0) print*,' << ',varname
    !
  end subroutine h5wa_r8
  !
  subroutine h5w_int4(varname,var)
    !
    use mpi, only: mpi_info_null
    !
    ! arguments
    character(LEN=*),intent(in) :: varname
    integer,intent(in) :: var
    !
    ! local data
    ! integer :: nvar(1)
    ! integer :: h5error
    ! integer(hsize_t) :: dimt(1)=(/1/)
    ! !
    ! ! writing the data
    ! !
    ! nvar=var
    ! call h5ltmake_dataset_f(h5file_id,varname,1,dimt,                  &
    !                                     h5t_native_integer,nvar,h5error)
    ! if(h5error.ne.0)  stop ' !! error in h5w_int4 call h5ltmake_dataset_f'
    ! !
    ! if(nrank==0) print*,' << ',varname
    !
  end subroutine h5w_int4
  !
  subroutine h5w_real8(varname,var)
    !
    use mpi, only: mpi_info_null
    !
    ! arguments
    character(LEN=*),intent(in) :: varname
    real(mytype),intent(in) :: var
    !
    ! ! local data
    ! real(mytype) :: rvar(1)
    ! integer :: h5error
    ! integer(hsize_t) :: dimt(1)=(/1/)
    ! !
    ! rvar=var
    ! call h5ltmake_dataset_f(h5file_id,varname,1,dimt,                  &
    !                                     h5t_native_double,rvar,h5error)
    ! if(h5error.ne.0)  stop ' !! error in h5w_real8 call h5ltmake_dataset_f'
    ! !
    ! if(nrank==0) print*,' << ',varname
    !
  end subroutine h5w_real8
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to create a sub-communicator from nranks.    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-08-2022: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine subcomm_group(rank,communicator,newrank,newsize)
    !
    use mpi
    ! arguments
    integer,intent(in) :: rank
    integer,intent(out) :: communicator,newrank,newsize
    !
    ! local data
    integer :: group_mpi,mpi_group_world
    integer :: ierr,ncout,jrank
    integer,allocatable :: rank_use(:),ranktemp(:)
    !
    allocate(ranktemp(0:nproc-1))
    !
    call pgather(rank,ranktemp)
    !
    ncout=0
    do jrank=0,nproc-1
      !
      if(ranktemp(jrank)>=0) then
        ncout=ncout+1
      endif
      !
    enddo
    !
    allocate(rank_use(1:ncout))
    !
    ncout=0
    do jrank=0,nproc-1
      !
      if(ranktemp(jrank)>=0) then
        ncout=ncout+1
        !
        rank_use(ncout)=ranktemp(jrank)
        !
      endif
      !
    enddo
    !
    call mpi_comm_group(mpi_comm_world,mpi_group_world,ierr)
    call mpi_group_incl(mpi_group_world,size(rank_use),rank_use,group_mpi,ierr)
    call mpi_comm_create(mpi_comm_world,group_mpi,communicator,ierr)
    !
    if(any(rank_use==nrank)) then
      call mpi_comm_size(communicator,newsize,ierr)
      call mpi_comm_rank(communicator,newrank,ierr)
      if(newrank==0) print*,' ** new subcomm created, size: ',newsize
      ! print*,' ** local rank:',newrank,', gloable rank:',nrank
    else
      newrank=-1
      newsize=0
    endif
    !
  end subroutine subcomm_group
  !+-------------------------------------------------------------------+
  !| The end of the subroutine subcomm_group.                          |
  !+-------------------------------------------------------------------+
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
    use decomp_2d, only : real_type
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
  !| This subroutine is used to read 1-D array via hdf5 interface.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 31-03-2022  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine h5_readarray1d(varname,var,dim,filename,explicit)
    !
    !
    real(mytype),intent(out) :: var(:)
    integer,intent(in) :: dim
    character(len=*),intent(in) :: varname,filename
    logical,intent(in), optional:: explicit
    !
    ! logical :: lexplicit
    ! !
    ! integer(hid_t) :: file_id
    ! ! file identifier
    ! integer(hid_t) :: dset_id1
    ! ! dataset identifier
    ! integer :: h5error ! error flag
    ! integer(hsize_t) :: dimt(1)
    ! !
    ! if (present(explicit)) then
    !    lexplicit = explicit
    ! else
    !    lexplicit = .true.
    ! end if
    ! !
    ! call h5open_f(h5error)
    ! !
    ! call h5fopen_f(filename,h5f_acc_rdwr_f,file_id,h5error)

    ! ! open an existing dataset.
    ! call h5dopen_f(file_id,varname,dset_id1,h5error)
    ! !
    ! dimt=(/dim/)
    ! !
    ! ! read the dataset.
    ! call h5dread_f(dset_id1,h5t_native_double,var,dimt,h5error)

    ! if(h5error.ne.0)  stop ' !! error in h5_readarray1d 1'
    ! !
    ! ! close the dataset
    ! call h5dclose_f(dset_id1, h5error)
    ! if(h5error.ne.0)  stop ' !! error in h5_readarray1d 2'
    ! ! close the file.
    ! call h5fclose_f(file_id, h5error)
    ! if(h5error.ne.0)  stop ' !! error in h5_readarray1d 3'
    ! !
    ! ! close fortran interface.
    ! call h5close_f(h5error)
    ! if(h5error.ne.0)  stop ' !! error in h5_readarray1d 4'
    ! !
    ! if(lexplicit)  print*,' >> ',varname,' from ',filename,' ... done'
    !
  end subroutine h5_readarray1d
  !
  subroutine h5_read1int(var,varname,filename,explicit)
    !
    integer,intent(out) :: var
    character(len=*),intent(in) :: varname,filename
    logical,intent(in), optional:: explicit

    ! logical :: lexplicit
    ! !
    ! integer(hid_t) :: file_id
    ! ! file identifier
    ! integer(hid_t) :: dset_id1
    ! ! dataset identifier
    ! integer :: v(1)
    ! integer :: h5error ! error flag
    ! integer(hsize_t) :: dimt(1)
    ! !
    ! if (present(explicit)) then
    !    lexplicit = explicit
    ! else
    !    lexplicit = .true.
    ! end if
    ! !
    ! dimt=(/1/)
    ! !
    ! call h5open_f(h5error)
    ! print*,' ** open hdf5 interface'
    ! !
    ! call h5fopen_f(filename,h5f_acc_rdwr_f,file_id,h5error)
    ! !
    ! call h5ltread_dataset_f(file_id,varname,h5t_native_integer,v,dimt,h5error)
    ! !
    ! call h5fclose_f(file_id,h5error)
    ! !
    ! if(h5error.ne.0)  stop ' !! error in h5_readarray1dint 1'
    ! !
    ! ! close fortran interface.
    ! call h5close_f(h5error)
    ! !
    ! var=v(1)
    ! if(h5error.ne.0)  stop ' !! error in h5_readarray1dint 2'
    ! !
    ! if(lexplicit)  print*,' >> ',varname,' from ',filename,' ... done'
    !
  end subroutine h5_read1int
  !
  subroutine h5_read1rl8(var,varname,filename,explicit)
    !
    real(mytype),intent(out) :: var
    character(len=*),intent(in) :: varname,filename
    logical,intent(in), optional:: explicit
    
    ! logical :: lexplicit
    ! !
    ! integer(hid_t) :: file_id
    ! ! file identifier
    ! integer(hid_t) :: dset_id1
    ! ! dataset identifier
    ! real(mytype) :: v(1)
    ! integer :: h5error ! error flag
    ! integer(hsize_t) :: dimt(1)
    ! !
    ! if (present(explicit)) then
    !    lexplicit = explicit
    ! else
    !    lexplicit = .true.
    ! end if
    ! !
    ! dimt=(/1/)
    ! !
    ! call h5open_f(h5error)
    ! if(lexplicit)  print*,' ** open hdf5 interface'
    ! !
    ! call h5fopen_f(filename,h5f_acc_rdwr_f,file_id,h5error)
    ! !
    ! call h5ltread_dataset_f(file_id,varname,h5t_native_double,v,dimt,h5error)
    ! !
    ! call h5fclose_f(file_id,h5error)
    ! !
    ! if(h5error.ne.0)  stop ' !! error in h5_readarray1dint 1'
    ! !
    ! ! close fortran interface.
    ! call h5close_f(h5error)
    ! !
    ! var=v(1)
    ! if(h5error.ne.0)  stop ' !! error in h5_readarray1dint 2'
    ! !
    ! if(lexplicit)  print*,' >> ',varname,' from ',filename,' ... done'
    !
    !
  end subroutine h5_read1rl8
  !+-------------------------------------------------------------------+
  !| This end of the subroutine h5_readarray1d.                        |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This function is used to get the dimension of the hdf5 array.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 30-JuL-2020 | Coped from ASTR Post by J. Fang STFC Daresbury Lab. |
  !+-------------------------------------------------------------------+
  function h5getdim3d(varname,filenma) result(dims)
    !
    character(len=*),intent(in) :: varname,filenma
    integer :: dims
    !
    ! ! local data
    ! integer(hid_t)  :: file, space, dset
    ! integer(hsize_t) :: ndims(1)
    ! integer         :: h5error ! error flag
    ! integer(hsize_t) :: dims_h5(1)
    ! !
    ! !
    ! call h5open_f(h5error)
    ! call h5fopen_f(filenma,h5f_acc_rdonly_f,file,h5error)
    ! call h5dopen_f (file,varname,dset, h5error)
    ! call h5dget_space_f(dset, space, h5error)
    ! call h5sget_simple_extent_dims_f(space,dims_h5,ndims,h5error)
    ! !
    ! dims=dims_h5(1)
    ! !
    ! call h5dclose_f(dset , h5error)
    ! call h5sclose_f(space, h5error)
    ! call h5fclose_f(file , h5error)
    ! call h5close_f(h5error)
    !
  end function h5getdim3d
  !+-------------------------------------------------------------------+
  !| This end of the function h5getdim3d.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function is to do the cross product for a 3-D vector.        | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
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
  !+-------------------------------------------------------------------+
  !| The end of the function cross_product.                            |
  !+-------------------------------------------------------------------+
  !
end module mptool
!+---------------------------------------------------------------------+
!| The end of the module mptool.                                       |
!+---------------------------------------------------------------------+