!                          -*- mode: f90 -*-
!! SPDX-License-Identifier: BSD-3-Clause

! This file contain common code to be included by subroutines
! 'mpiio_write_plane_3d_...' in io.f90

! It is much easier to implement if all mpi ranks participate I/O.
! Transpose the 3D data if necessary.

if (present(opt_decomp)) then
   decomp = opt_decomp
else
   call get_decomp_info(decomp)
end if

! The current GPU build will not perform transpose operations when writing planes
#if defined(_GPU)
if (ipencil /= iplane) then
   call decomp_2d_warning(__FILE__, &
                          __LINE__, &
                          ipencil, &
                          "GPU build does not allow transpose operations when writing planes")
   return
endif
#endif

opened_new = .false.

if (iplane == 1) then
   allocate (wk(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)))
   if (ipencil == 1) then
      wk = var
   else if (ipencil == 2) then
      call transpose_y_to_x(var, wk, decomp)
   else if (ipencil == 3) then
      allocate (wk2(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)))
      call transpose_z_to_y(var, wk2, decomp)
      call transpose_y_to_x(wk2, wk, decomp)
      deallocate (wk2)
   end if
   allocate (wk2d(1, decomp%xsz(2), decomp%xsz(3)))
   if (n >= 1) then
      do k = 1, decomp%xsz(3)
         do j = 1, decomp%xsz(2)
            wk2d(1, j, k) = wk(n, j, k)
         end do
      end do
   else
      do k = 1, decomp%xsz(3)
         do j = 1, decomp%xsz(2)
            wk2d(1, j, k) = sum(wk(:, j, k))/real(decomp%xsz(1), kind=mytype)
         end do
      end do
   end if
else if (iplane == 2) then
   allocate (wk(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)))
   if (ipencil == 1) then
      call transpose_x_to_y(var, wk, decomp)
   else if (ipencil == 2) then
      wk = var
   else if (ipencil == 3) then
      call transpose_z_to_y(var, wk, decomp)
   end if
   allocate (wk2d(decomp%ysz(1), 1, decomp%ysz(3)))
   if (n >= 1) then
      do k = 1, decomp%ysz(3)
         do i = 1, decomp%ysz(1)
            wk2d(i, 1, k) = wk(i, n, k)
         end do
      end do
   else
      do k = 1, decomp%ysz(3)
         do i = 1, decomp%ysz(1)
            wk2d(i, 1, k) = sum(wk(i, :, k))/real(decomp%ysz(2), kind=mytype)
         end do
      end do
   end if
else if (iplane == 3) then
   allocate (wk(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)))
   if (ipencil == 1) then
      allocate (wk2(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)))
      call transpose_x_to_y(var, wk2, decomp)
      call transpose_y_to_z(wk2, wk, decomp)
      deallocate (wk2)
   else if (ipencil == 2) then
      call transpose_y_to_z(var, wk, decomp)
   else if (ipencil == 3) then
      wk = var
   end if
   allocate (wk2d(decomp%zsz(1), decomp%zsz(2), 1))
   if (n >= 1) then
      do j = 1, decomp%zsz(2)
         do i = 1, decomp%zsz(1)
            wk2d(i, j, 1) = wk(i, j, n)
         end do
      end do
   else
      do j = 1, decomp%zsz(2)
         do i = 1, decomp%zsz(1)
            wk2d(i, j, 1) = sum(wk(i, j, :))/real(decomp%zsz(3), kind=mytype)
         end do
      end do
   end if
end if

idx = get_io_idx(io_name, dirname)

#ifndef ADIOS2
!! Use default MPIIO writers

if (idx < 1) then
   ! Create folder if needed
   if (nrank == 0) then
      inquire (file=dirname, exist=dir_exists)
      if (.not. dir_exists) then
         call execute_command_line("mkdir "//dirname//" 2> /dev/null", wait=.true.)
      end if
   end if
   allocate (character(len(trim(dirname)) + 1 + len(trim(varname))) :: full_io_name)
   full_io_name = dirname//"/"//varname
   call decomp_2d_open_io(io_name, full_io_name, decomp_2d_write_mode)
   idx = get_io_idx(io_name, full_io_name)
   opened_new = .true.
end if

call plane_extents(sizes, subsizes, starts, iplane, decomp)
associate (fh => fh_registry(idx), &
           disp => fh_disp(idx))
   call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts, &
                                 MPI_ORDER_FORTRAN, data_type, newtype, ierror)
   if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_CREATE_SUBARRAY")
   call MPI_TYPE_COMMIT(newtype, ierror)
   if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_COMMIT")
   call MPI_FILE_SET_VIEW(fh, disp, data_type, &
                          newtype, 'native', MPI_INFO_NULL, ierror)
   if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_SET_VIEW")
   call MPI_FILE_WRITE_ALL(fh, wk2d, &
                           subsizes(1)*subsizes(2)*subsizes(3), &
                           data_type, MPI_STATUS_IGNORE, ierror)
   if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_WRITE_ALL")
   call MPI_TYPE_FREE(newtype, ierror)
   if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_FREE")
end associate

if (opened_new) then
   call decomp_2d_close_io(io_name, full_io_name)
   deallocate (full_io_name)
end if
#else
!! Write using ADIOS2
call adios2_at_io(io_handle, adios, io_name, ierror)
if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "adios2_at_io "//trim(io_name))
call adios2_inquire_variable(var_handle, io_handle, varname, ierror)
if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "adios2_inquire_variable "//trim(varname))
if (.not. var_handle%valid) then
   call decomp_2d_abort(__FILE__, __LINE__, -1, &
        "trying to write variable before registering!  "//trim(varname))
end if

!! Note - need to use sync mode as the array for the output plane gets reused.
call adios2_put(engine_registry(idx), var_handle, wk2d, adios2_mode_sync, ierror)
if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "adios2_put")
#endif

deallocate (wk, wk2d)
