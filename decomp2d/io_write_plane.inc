!                          -*- mode: f90 -*-                            
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
! 'mpiio_write_plane_3d_...' in io.f90

! It is much easier to implement if all mpi ranks participate I/O.
! Transpose the 3D data if necessary.

if (present(opt_decomp)) then
   decomp = opt_decomp
else
   call get_decomp_info(decomp)
end if

if (iplane==1) then
   allocate(wk(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)))
   if (ipencil==1) then
      wk = var
   else if (ipencil==2) then
      call transpose_y_to_x(var,wk,decomp)
   else if (ipencil==3) then
      allocate(wk2(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
      call transpose_z_to_y(var,wk2,decomp)
      call transpose_y_to_x(wk2,wk,decomp)
      deallocate(wk2)
   end if
   allocate(wk2d(1,decomp%xsz(2),decomp%xsz(3)))
   if (n.ge.1) then
      do k=1,decomp%xsz(3)
         do j=1,decomp%xsz(2)
            wk2d(1,j,k)=wk(n,j,k)
         end do
      end do
   else
      do k=1,decomp%xsz(3)
         do j=1,decomp%xsz(2)
            wk2d(1,j,k)=sum(wk(:,j,k))/real(decomp%xsz(1),kind=mytype)
         end do
      end do
   endif
else if (iplane==2) then
   allocate(wk(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
   if (ipencil==1) then
      call transpose_x_to_y(var,wk,decomp)
   else if (ipencil==2) then
      wk = var
   else if (ipencil==3) then
      call transpose_z_to_y(var,wk,decomp)
   end if
   allocate(wk2d(decomp%ysz(1),1,decomp%ysz(3)))
   if (n.ge.1) then
      do k=1,decomp%ysz(3)
         do i=1,decomp%ysz(1)
            wk2d(i,1,k)=wk(i,n,k)
         end do
      end do
   else
      do k=1,decomp%ysz(3)
         do i=1,decomp%ysz(1)
            wk2d(i,1,k)=sum(wk(i,:,k))/real(decomp%ysz(2),kind=mytype)
         end do
      end do
   endif
else if (iplane==3) then
   allocate(wk(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)))
   if (ipencil==1) then
      allocate(wk2(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
      call transpose_x_to_y(var,wk2,decomp)
      call transpose_y_to_z(wk2,wk,decomp)
      deallocate(wk2)
   else if (ipencil==2) then
      call transpose_y_to_z(var,wk,decomp)
   else if (ipencil==3) then
      wk = var
   end if
   allocate(wk2d(decomp%zsz(1),decomp%zsz(2),1))
   if (n.ge.1) then
      do j=1,decomp%zsz(2)
         do i=1,decomp%zsz(1) 
            wk2d(i,j,1)=wk(i,j,n)
         end do
      end do
   else
      do j=1,decomp%zsz(2)
         do i=1,decomp%zsz(1)
            wk2d(i,j,1)=sum(wk(i,j,:))/real(decomp%zsz(3),kind=mytype)
         end do
      end do
   endif
end if

#ifndef ADIOS2
!! Use default MPIIO writers
call plane_extents(sizes, subsizes, starts, iplane, decomp)
call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
     MPI_ORDER_FORTRAN, data_type, newtype, ierror)
call MPI_TYPE_COMMIT(newtype,ierror)
call MPI_FILE_OPEN(MPI_COMM_WORLD, varname, &
     MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
     fh, ierror)
filesize = 0_MPI_OFFSET_KIND
call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
disp = 0_MPI_OFFSET_KIND
call MPI_FILE_SET_VIEW(fh,disp,data_type, &
     newtype,'native',MPI_INFO_NULL,ierror)
call MPI_FILE_WRITE_ALL(fh, wk2d, &
     subsizes(1)*subsizes(2)*subsizes(3), &
     data_type, MPI_STATUS_IGNORE, ierror)
call MPI_FILE_CLOSE(fh,ierror)
call MPI_TYPE_FREE(newtype,ierror)
#else
!! Write using ADIOS2
call adios2_at_io(io_handle, adios, io_name, ierror)
call adios2_inquire_variable(var_handle, io_handle, varname, ierror)
if (.not.var_handle % valid) then
   print *, "ERROR: trying to write variable before registering!", varname
   stop
endif

idx = get_engine_idx(io_name, dirname)

!! Note - need to use sync mode as the array for the output plane gets reused.
call adios2_put(engine_registry(idx), var_handle, wk2d, adios2_mode_sync, ierror)
#endif

deallocate(wk,wk2d)
