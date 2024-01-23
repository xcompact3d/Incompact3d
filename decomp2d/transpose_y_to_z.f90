!! SPDX-License-Identifier: BSD-3-Clause

! This file contains the routines that transpose data from Y to Z pencil

  subroutine transpose_y_to_z_real_short(src, dst)

     implicit none

     real(mytype), dimension(:, :, :), intent(IN) :: src
     real(mytype), dimension(:, :, :), intent(OUT) :: dst

     call transpose_y_to_z(src, dst, decomp_main)

  end subroutine transpose_y_to_z_real_short

  subroutine transpose_y_to_z_real_long(src, dst, decomp)

     implicit none

     real(mytype), dimension(:, :, :), intent(IN) :: src
     real(mytype), dimension(:, :, :), intent(OUT) :: dst
     TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
     integer :: istat, nsize
#endif

     if (dims(2) == 1) then
#if defined(_GPU)
        nsize = product(decomp%ysz)
        !$acc host_data use_device(src,dst)
        istat = cudaMemcpy(dst, src, nsize, cudaMemcpyDeviceToDevice)
        !$acc end host_data
        if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy")
#else
        dst = src
#endif
     else
        call transpose_y_to_z_real(src, dst, decomp)
     end if

  end subroutine transpose_y_to_z_real_long

  subroutine transpose_y_to_z_real(src, dst, decomp)

     implicit none

     real(mytype), dimension(:, :, :), intent(IN) :: src
     real(mytype), dimension(:, :, :), intent(OUT) :: dst
     TYPE(DECOMP_INFO), intent(IN) :: decomp

#if defined(_GPU)
     integer :: istat
#endif

     integer :: s1, s2, s3, d1, d2, d3
     integer :: ierror

#ifdef PROFILER
     if (decomp_profiler_transpose) call decomp_profiler_start("transp_y_z_r")
#endif

     s1 = SIZE(src, 1)
     s2 = SIZE(src, 2)
     s3 = SIZE(src, 3)
     d1 = SIZE(dst, 1)
     d2 = SIZE(dst, 2)
     d3 = SIZE(dst, 3)

     ! rearrange source array as send buffer
#if defined(_GPU)
     call mem_split_yz_real(src, s1, s2, s3, work1_r_d, dims(2), &
                            decomp%y2dist, decomp)
#else
     call mem_split_yz_real(src, s1, s2, s3, work1_r, dims(2), &
                            decomp%y2dist, decomp)
#endif

     ! define receive buffer
#ifdef EVEN
     if (decomp%even) then
        call MPI_ALLTOALL(work1_r, decomp%y2count, &
                          real_type, dst, decomp%z2count, &
                          real_type, DECOMP_2D_COMM_ROW, ierror)
        if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
     else
        call MPI_ALLTOALL(work1_r, decomp%y2count, &
                          real_type, work2_r, decomp%z2count, &
                          real_type, DECOMP_2D_COMM_ROW, ierror)
        if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALL")
     end if
#else

#if defined(_GPU)
#if defined(_NCCL)
     call decomp_2d_nccl_send_recv_row(work2_r_d, &
                                       work1_r_d, &
                                       decomp%y2disp, &
                                       decomp%y2cnts, &
                                       decomp%z2disp, &
                                       decomp%z2cnts, &
                                       dims(2))
#else
     call MPI_ALLTOALLV(work1_r_d, decomp%y2cnts, decomp%y2disp, &
                        real_type, work2_r_d, decomp%z2cnts, decomp%z2disp, &
                        real_type, DECOMP_2D_COMM_ROW, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif
#else
     call MPI_ALLTOALLV(work1_r, decomp%y2cnts, decomp%y2disp, &
                        real_type, dst, decomp%z2cnts, decomp%z2disp, &
                        real_type, DECOMP_2D_COMM_ROW, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif

#endif

     ! rearrange receive buffer
#ifdef EVEN
     if (.not. decomp%even) then
        call mem_merge_yz_real(work2_r, d1, d2, d3, dst, dims(2), &
                               decomp%z2dist, decomp)
     end if
#else
     ! note the receive buffer is already in natural (i,j,k) order
     ! so no merge operation needed

#if defined(_GPU)
     !If one of the array in cuda call is not device we need to add acc host_data
     !$acc host_data use_device(dst)
     istat = cudaMemcpy(dst, work2_r_d, d1 * d2 * d3, cudaMemcpyDeviceToDevice)
     !$acc end host_data
     if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy2D")
#endif

#endif

#ifdef PROFILER
     if (decomp_profiler_transpose) call decomp_profiler_end("transp_y_z_r")
#endif

     return
  end subroutine transpose_y_to_z_real

  subroutine transpose_y_to_z_complex_short(src, dst)

     implicit none

     complex(mytype), dimension(:, :, :), intent(IN) :: src
     complex(mytype), dimension(:, :, :), intent(OUT) :: dst

     call transpose_y_to_z(src, dst, decomp_main)

  end subroutine transpose_y_to_z_complex_short

  subroutine transpose_y_to_z_complex_long(src, dst, decomp)

     implicit none

     complex(mytype), dimension(:, :, :), intent(IN) :: src
     complex(mytype), dimension(:, :, :), intent(OUT) :: dst
     TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
     integer :: istat, nsize
#endif

     if (dims(2) == 1) then
#if defined(_GPU)
        nsize = product(decomp%ysz)
        !$acc host_data use_device(src,dst)
        istat = cudaMemcpy(dst, src, nsize, cudaMemcpyDeviceToDevice)
        !$acc end host_data
        if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy")
#else
        dst = src
#endif
     else
        call transpose_y_to_z_complex(src, dst, decomp)
     end if

  end subroutine transpose_y_to_z_complex_long

  subroutine transpose_y_to_z_complex(src, dst, decomp)

     implicit none

     complex(mytype), dimension(:, :, :), intent(IN) :: src
     complex(mytype), dimension(:, :, :), intent(OUT) :: dst
     TYPE(DECOMP_INFO), intent(IN) :: decomp

#if defined(_GPU)
     integer :: istat
#endif

     integer :: s1, s2, s3, d1, d2, d3
     integer :: ierror

#ifdef PROFILER
     if (decomp_profiler_transpose) call decomp_profiler_start("transp_y_z_c")
#endif

     s1 = SIZE(src, 1)
     s2 = SIZE(src, 2)
     s3 = SIZE(src, 3)
     d1 = SIZE(dst, 1)
     d2 = SIZE(dst, 2)
     d3 = SIZE(dst, 3)

     ! rearrange source array as send buffer
#if defined(_GPU)
     call mem_split_yz_complex(src, s1, s2, s3, work1_c_d, dims(2), &
                               decomp%y2dist, decomp)
#else
     call mem_split_yz_complex(src, s1, s2, s3, work1_c, dims(2), &
                               decomp%y2dist, decomp)
#endif

     ! define receive buffer
#ifdef EVEN
     if (decomp%even) then
        call MPI_ALLTOALL(work1_c, decomp%y2count, &
                          complex_type, dst, decomp%z2count, &
                          complex_type, DECOMP_2D_COMM_ROW, ierror)
     else
        call MPI_ALLTOALL(work1_c, decomp%y2count, &
                          complex_type, work2_c, decomp%z2count, &
                          complex_type, DECOMP_2D_COMM_ROW, ierror)
     end if
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALL")
#else

#if defined(_GPU)
#if defined(_NCCL)
     call decomp_2d_nccl_send_recv_row(work2_c_d, &
                                       work1_c_d, &
                                       decomp%y2disp, &
                                       decomp%y2cnts, &
                                       decomp%z2disp, &
                                       decomp%z2cnts, &
                                       dims(2), &
                                       decomp_buf_size)
#else
     call MPI_ALLTOALLV(work1_c_d, decomp%y2cnts, decomp%y2disp, &
                        complex_type, work2_c_d, decomp%z2cnts, decomp%z2disp, &
                        complex_type, DECOMP_2D_COMM_ROW, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif
#else
     call MPI_ALLTOALLV(work1_c, decomp%y2cnts, decomp%y2disp, &
                        complex_type, dst, decomp%z2cnts, decomp%z2disp, &
                        complex_type, DECOMP_2D_COMM_ROW, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif

#endif

     ! rearrange receive buffer
#ifdef EVEN
     if (.not. decomp%even) then
        call mem_merge_yz_complex(work2_c, d1, d2, d3, dst, dims(2), &
                                  decomp%z2dist, decomp)
     end if
#else
     ! note the receive buffer is already in natural (i,j,k) order
     ! so no merge operation needed

#if defined(_GPU)
     !$acc host_data use_device(dst)
     istat = cudaMemcpy(dst, work2_c_d, d1 * d2 * d3, cudaMemcpyDeviceToDevice)
     !$acc end host_data
     if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy2D")
#endif

#endif

#ifdef PROFILER
     if (decomp_profiler_transpose) call decomp_profiler_end("transp_y_z_c")
#endif

     return
  end subroutine transpose_y_to_z_complex

  subroutine mem_split_yz_real(in, n1, n2, n3, out, iproc, dist, decomp)

     implicit none

     integer, intent(IN) :: n1, n2, n3
     real(mytype), dimension(n1, n2, n3), intent(IN) :: in
     real(mytype), dimension(*), intent(OUT) :: out
     integer, intent(IN) :: iproc
     integer, dimension(0:iproc - 1), intent(IN) :: dist
     TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
     attributes(device) :: out
     integer :: istat
#endif

     integer :: i, j, k, m, i1, i2, pos

     do m = 0, iproc - 1
        if (m == 0) then
           i1 = 1
           i2 = dist(0)
        else
           i1 = i2 + 1
           i2 = i1 + dist(m) - 1
        end if

#ifdef EVEN
        pos = m * decomp%y2count + 1
#else
        pos = decomp%y2disp(m) + 1
#endif

#if defined(_GPU)
        !$acc host_data use_device(in)
        istat = cudaMemcpy2D(out(pos), n1 * (i2 - i1 + 1), in(1, i1, 1), n1 * n2, n1 * (i2 - i1 + 1), n3, cudaMemcpyDeviceToDevice)
        !$acc end host_data
        if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy2D")
#else
        do k = 1, n3
           do j = i1, i2
              do i = 1, n1
                 out(pos) = in(i, j, k)
                 pos = pos + 1
              end do
           end do
        end do
#endif
     end do

     return
  end subroutine mem_split_yz_real

  subroutine mem_split_yz_complex(in, n1, n2, n3, out, iproc, dist, decomp)

     implicit none

     integer, intent(IN) :: n1, n2, n3
     complex(mytype), dimension(n1, n2, n3), intent(IN) :: in
     complex(mytype), dimension(*), intent(OUT) :: out
     integer, intent(IN) :: iproc
     integer, dimension(0:iproc - 1), intent(IN) :: dist
     TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
     attributes(device) :: out
     integer :: istat
#endif

     integer :: i, j, k, m, i1, i2, pos

     do m = 0, iproc - 1
        if (m == 0) then
           i1 = 1
           i2 = dist(0)
        else
           i1 = i2 + 1
           i2 = i1 + dist(m) - 1
        end if

#ifdef EVEN
        pos = m * decomp%y2count + 1
#else
        pos = decomp%y2disp(m) + 1
#endif

#if defined(_GPU)
        !$acc host_data use_device(in)
        istat = cudaMemcpy2D(out(pos), n1 * (i2 - i1 + 1), in(1, i1, 1), n1 * n2, n1 * (i2 - i1 + 1), n3, cudaMemcpyDeviceToDevice)
        !$acc end host_data
        if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy2D")
#else
        do k = 1, n3
           do j = i1, i2
              do i = 1, n1
                 out(pos) = in(i, j, k)
                 pos = pos + 1
              end do
           end do
        end do
#endif
     end do

     return
  end subroutine mem_split_yz_complex

  subroutine mem_merge_yz_real(in, n1, n2, n3, out, iproc, dist, decomp)

     implicit none

     integer, intent(IN) :: n1, n2, n3
     real(mytype), dimension(*), intent(IN) :: in
     real(mytype), dimension(n1, n2, n3), intent(OUT) :: out
     integer, intent(IN) :: iproc
     integer, dimension(0:iproc - 1), intent(IN) :: dist
     TYPE(DECOMP_INFO), intent(IN) :: decomp

     integer :: i, j, k, m, i1, i2, pos

     do m = 0, iproc - 1
        if (m == 0) then
           i1 = 1
           i2 = dist(0)
        else
           i1 = i2 + 1
           i2 = i1 + dist(m) - 1
        end if

#ifdef EVEN
        pos = m * decomp%z2count + 1
#else
        pos = decomp%z2disp(m) + 1
#endif

        do k = i1, i2
           do j = 1, n2
              do i = 1, n1
                 out(i, j, k) = in(pos)
                 pos = pos + 1
              end do
           end do
        end do
     end do

     return
  end subroutine mem_merge_yz_real

  subroutine mem_merge_yz_complex(in, n1, n2, n3, out, iproc, dist, decomp)

     implicit none

     integer, intent(IN) :: n1, n2, n3
     complex(mytype), dimension(*), intent(IN) :: in
     complex(mytype), dimension(n1, n2, n3), intent(OUT) :: out
     integer, intent(IN) :: iproc
     integer, dimension(0:iproc - 1), intent(IN) :: dist
     TYPE(DECOMP_INFO), intent(IN) :: decomp

     integer :: i, j, k, m, i1, i2, pos

     do m = 0, iproc - 1
        if (m == 0) then
           i1 = 1
           i2 = dist(0)
        else
           i1 = i2 + 1
           i2 = i1 + dist(m) - 1
        end if

#ifdef EVEN
        pos = m * decomp%z2count + 1
#else
        pos = decomp%z2disp(m) + 1
#endif

        do k = i1, i2
           do j = 1, n2
              do i = 1, n1
                 out(i, j, k) = in(pos)
                 pos = pos + 1
              end do
           end do
        end do
     end do

     return
  end subroutine mem_merge_yz_complex
