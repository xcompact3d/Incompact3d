!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2011 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This file contain common code to be included by subroutines 
! 'write_var_...' in io.f90

  ! Using MPI-IO to write a distributed 3D variable to a file. File 
  ! operations (open/close) need to be done in calling application. This
  ! allows multiple variables to be written to a single file. Together 
  ! with the corresponding read operation, this is the perfect solution
  ! for applications to perform restart/checkpointing.

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       call get_decomp_info(decomp)
    end if

    ! Create file type and set file view
    sizes(1) = ntimesteps 
    sizes(2) = decomp%ysz(2)
    sizes(3) = decomp%zsz(3)
    subsizes(1) = ntimesteps 
    subsizes(2) = decomp%xsz(2)
    subsizes(3) = decomp%xsz(3)
    starts(1) = 0  ! 0-based index
    starts(2) = decomp%xst(2)-1
    starts(3) = decomp%xst(3)-1
   
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, data_type, newtype, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_CREATE_SUBARRAY")
    call MPI_TYPE_COMMIT(newtype,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_COMMIT")
    call MPI_FILE_SET_VIEW(fh,disp,data_type, &
         newtype,'native',MPI_INFO_NULL,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_SET_VIEW")
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         data_type, MPI_STATUS_IGNORE, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_WRITE_ALL")
    call MPI_TYPE_FREE(newtype,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_FREE")

    ! update displacement for the next write operation
    disp = disp + sizes(1)*sizes(2)*sizes(3)*mytype_bytes
    if (data_type == complex_type) then
       disp = disp + sizes(1)*sizes(2)*sizes(3)*mytype_bytes
    end if
