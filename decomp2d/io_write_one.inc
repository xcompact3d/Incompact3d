!                          -*- mode: f90 -*-                            
!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG) 
! Copyright (C) 2021               the University of Edinburgh (UoE)
!
!=======================================================================

! This file contain common code to be included by subroutines 
! 'mpiio_write_one_...' in io.f90

! Using MPI-IO to write a distributed 3D array into a file

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

#ifndef ADIOS2
!! Original writers

#ifdef T3PIO
call MPI_INFO_CREATE(info, ierror)
gs = ceiling(real(sizes(1),mytype)*real(sizes(2),mytype)* &
     real(sizes(3),mytype)/1024./1024.)
call t3pio_set_info(MPI_COMM_WORLD, info, "./", ierror, &
     GLOBAL_SIZE=gs, factor=1)
#endif

call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
     MPI_ORDER_FORTRAN, data_type, newtype, ierror)
call MPI_TYPE_COMMIT(newtype,ierror)
#ifdef T3PIO
call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
     MPI_MODE_CREATE+MPI_MODE_WRONLY, info, fh, ierror)
#else
call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
     MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
     fh, ierror)
#endif
filesize = 0_MPI_OFFSET_KIND
call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
disp = 0_MPI_OFFSET_KIND
call MPI_FILE_SET_VIEW(fh,disp,data_type, &
     newtype,'native',MPI_INFO_NULL,ierror)
call MPI_FILE_WRITE_ALL(fh, var, &
     subsizes(1)*subsizes(2)*subsizes(3), &
     data_type, MPI_STATUS_IGNORE, ierror)
call MPI_FILE_CLOSE(fh,ierror)
call MPI_TYPE_FREE(newtype,ierror)
#ifdef T3PIO
call MPI_INFO_FREE(info,ierror)
#endif

#else
!! Using ADIOS2

print *, "Trying to use ADIOS2 - interface not completed (WIP), please try again later."
print *, "For time being, please build without ADIOS2 to use default MPI-IO interface."
call MPI_ABORT(MPI_COMM_WORLD, -1, ierror)

! Check if variable already exists, if not create it
call adios2_inquire_variable(var_handle, io_write_one, filename, ierror)
if (.not.var_handle%valid) then
   !! New variable
   
   ! Need to set the ADIOS2 data type
   if (data_type.eq.MPI_DOUBLE_PRECISION) then
      !! Double
      data_type = adios2_type_dp
   else if (data_type.eq.MPI_REAL) then
      !! Single
      data_type = adios2_type_real
   else if (data_type.eq.MPI_DOUBLE_COMPLEX) then
      !! Complex double
      data_type = adios2_type_complex_dp
   else if (data_type.eq.MPI_COMPLEX) then
      !! Complex single
      data_type = adios2_type_complex
   else
      print *, "Trying to write unknown data type!"
      call MPI_ABORT(MPI_COMM_WORLD, -1, ierror)
   endif

   call adios2_define_variable(var_handle, io_write_one, filename, data_type, &
        ndims, int(sizes, kind=8), int(starts, kind=8), int(subsizes, kind=8), &
        adios2_constant_dims, ierror)
endif

call adios2_put(engine_write_one, var_handle, var, ierror)
#endif
