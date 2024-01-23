!! SPDX-License-Identifier: BSD-3-Clause

! This file contains the routines that transpose data from X to Y pencil

  subroutine transpose_x_to_y_real_short(src, dst)

     implicit none

     real(mytype), dimension(:, :, :), intent(IN) :: src
     real(mytype), dimension(:, :, :), intent(OUT) :: dst

     call transpose_x_to_y(src, dst, decomp_main)

  end subroutine transpose_x_to_y_real_short

  subroutine transpose_x_to_y_real_long(src, dst, decomp)

     implicit none

     real(mytype), dimension(:, :, :), intent(IN) :: src
     real(mytype), dimension(:, :, :), intent(OUT) :: dst
     TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
     integer :: istat, nsize
#endif

     if (dims(1) == 1) then
#if defined(_GPU)
        nsize = product(decomp%xsz)
        !$acc host_data use_device(src,dst)
        istat = cudaMemcpy(dst, src, nsize, cudaMemcpyDeviceToDevice)
        !$acc end host_data
        if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy")
#else
        dst = src
#endif
     else
        call transpose_x_to_y_real(src, dst, decomp)
     end if

  end subroutine transpose_x_to_y_real_long

  subroutine transpose_x_to_y_real(src, dst, decomp)

     implicit none

     Real(mytype), dimension(:, :, :), intent(IN) :: src
     real(mytype), dimension(:, :, :), intent(OUT) :: dst
     TYPE(DECOMP_INFO), intent(IN) :: decomp

     integer :: s1, s2, s3, d1, d2, d3
     integer :: ierror

#ifdef PROFILER
     if (decomp_profiler_transpose) call decomp_profiler_start("transp_x_y_r")
#endif

     s1 = SIZE(src, 1)
     s2 = SIZE(src, 2)
     s3 = SIZE(src, 3)
     d1 = SIZE(dst, 1)
     d2 = SIZE(dst, 2)
     d3 = SIZE(dst, 3)

     ! rearrange source array as send buffer
#if defined(_GPU)
     call mem_split_xy_real(src, s1, s2, s3, work1_r_d, dims(1), &
                            decomp%x1dist, decomp)
#else
     call mem_split_xy_real(src, s1, s2, s3, work1_r, dims(1), &
                            decomp%x1dist, decomp)
#endif

     ! define receive buffer
     ! transpose using MPI_ALLTOALL(V)
#ifdef EVEN
     call MPI_ALLTOALL(work1_r, decomp%x1count, &
                       real_type, work2_r, decomp%y1count, &
                       real_type, DECOMP_2D_COMM_COL, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALL")
#else

#if defined(_GPU)
#if defined(_NCCL)
     call decomp_2d_nccl_send_recv_col(work2_r_d, &
                                       work1_r_d, &
                                       decomp%x1disp, &
                                       decomp%x1cnts, &
                                       decomp%y1disp, &
                                       decomp%y1cnts, &
                                       dims(1))
#else
     call MPI_ALLTOALLV(work1_r_d, decomp%x1cnts, decomp%x1disp, &
                        real_type, work2_r_d, decomp%y1cnts, decomp%y1disp, &
                        real_type, DECOMP_2D_COMM_COL, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif
#else
     call MPI_ALLTOALLV(work1_r, decomp%x1cnts, decomp%x1disp, &
                        real_type, work2_r, decomp%y1cnts, decomp%y1disp, &
                        real_type, DECOMP_2D_COMM_COL, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif

#endif

     ! rearrange receive buffer
#if defined(_GPU)
     call mem_merge_xy_real(work2_r_d, d1, d2, d3, dst, dims(1), &
                            decomp%y1dist, decomp)
#else
     call mem_merge_xy_real(work2_r, d1, d2, d3, dst, dims(1), &
                            decomp%y1dist, decomp)
#endif

#ifdef PROFILER
     if (decomp_profiler_transpose) call decomp_profiler_end("transp_x_y_r")
#endif

     return
  end subroutine transpose_x_to_y_real

  subroutine transpose_x_to_y_complex_short(src, dst)

     implicit none

     complex(mytype), dimension(:, :, :), intent(IN) :: src
     complex(mytype), dimension(:, :, :), intent(OUT) :: dst

     call transpose_x_to_y(src, dst, decomp_main)

  end subroutine transpose_x_to_y_complex_short

  subroutine transpose_x_to_y_complex_long(src, dst, decomp)

     implicit none

     complex(mytype), dimension(:, :, :), intent(IN) :: src
     complex(mytype), dimension(:, :, :), intent(OUT) :: dst
     TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
     integer :: istat, nsize
#endif
     if (dims(1) == 1) then
#if defined(_GPU)
        nsize = product(decomp%xsz)
        !$acc host_data use_device(src,dst)
        istat = cudaMemcpy(dst, src, nsize, cudaMemcpyDeviceToDevice)
        !$acc end host_data
        if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy")
#else
        dst = src
#endif
     else
        call transpose_x_to_y_complex(src, dst, decomp)
     end if

  end subroutine transpose_x_to_y_complex_long

  subroutine transpose_x_to_y_complex(src, dst, decomp)

     implicit none

     complex(mytype), dimension(:, :, :), intent(IN) :: src
     complex(mytype), dimension(:, :, :), intent(OUT) :: dst
     TYPE(DECOMP_INFO), intent(IN) :: decomp

     integer :: s1, s2, s3, d1, d2, d3
     integer :: ierror

#ifdef PROFILER
     if (decomp_profiler_transpose) call decomp_profiler_start("transp_x_y_c")
#endif

     s1 = SIZE(src, 1)
     s2 = SIZE(src, 2)
     s3 = SIZE(src, 3)
     d1 = SIZE(dst, 1)
     d2 = SIZE(dst, 2)
     d3 = SIZE(dst, 3)

     ! rearrange source array as send buffer
#if defined(_GPU)
     call mem_split_xy_complex(src, s1, s2, s3, work1_c_d, dims(1), &
                               decomp%x1dist, decomp)
#else
     call mem_split_xy_complex(src, s1, s2, s3, work1_c, dims(1), &
                               decomp%x1dist, decomp)
#endif

     ! define receive buffer
     ! transpose using MPI_ALLTOALL(V)
#ifdef EVEN
     call MPI_ALLTOALL(work1_c, decomp%x1count, &
                       complex_type, work2_c, decomp%y1count, &
                       complex_type, DECOMP_2D_COMM_COL, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALL")
#else

#if defined(_GPU)
#if defined(_NCCL)
     call decomp_2d_nccl_send_recv_col(work2_c_d, &
                                       work1_c_d, &
                                       decomp%x1disp, &
                                       decomp%x1cnts, &
                                       decomp%y1disp, &
                                       decomp%y1cnts, &
                                       dims(1), &
                                       decomp_buf_size)
#else
     call MPI_ALLTOALLV(work1_c_d, decomp%x1cnts, decomp%x1disp, &
                        complex_type, work2_c_d, decomp%y1cnts, decomp%y1disp, &
                        complex_type, DECOMP_2D_COMM_COL, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif
#else
     call MPI_ALLTOALLV(work1_c, decomp%x1cnts, decomp%x1disp, &
                        complex_type, work2_c, decomp%y1cnts, decomp%y1disp, &
                        complex_type, DECOMP_2D_COMM_COL, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif

#endif

     ! rearrange receive buffer
#if defined(_GPU)
     call mem_merge_xy_complex(work2_c_d, d1, d2, d3, dst, dims(1), &
                               decomp%y1dist, decomp)
#else
     call mem_merge_xy_complex(work2_c, d1, d2, d3, dst, dims(1), &
                               decomp%y1dist, decomp)
#endif

#ifdef PROFILER
     if (decomp_profiler_transpose) call decomp_profiler_end("transp_x_y_c")
#endif

     return
  end subroutine transpose_x_to_y_complex

  subroutine mem_split_xy_real(in, n1, n2, n3, out, iproc, dist, decomp)

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
        pos = m * decomp%x1count + 1
#else
        pos = decomp%x1disp(m) + 1
#endif

#if defined(_GPU)
        !$acc host_data use_device(in)
        istat = cudaMemcpy2D(out(pos), i2 - i1 + 1, in(i1, 1, 1), n1, i2 - i1 + 1, n2 * n3, cudaMemcpyDeviceToDevice)
        !$acc end host_data
        if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy2D")
#else
        do k = 1, n3
           do j = 1, n2
              do i = i1, i2
                 out(pos) = in(i, j, k)
                 pos = pos + 1
              end do
           end do
        end do
#endif
     end do

     return
  end subroutine mem_split_xy_real

  subroutine mem_split_xy_complex(in, n1, n2, n3, out, iproc, dist, decomp)

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
        pos = m * decomp%x1count + 1
#else
        pos = decomp%x1disp(m) + 1
#endif

#if defined(_GPU)
        !$acc host_data use_device(in)
        istat = cudaMemcpy2D(out(pos), i2 - i1 + 1, in(i1, 1, 1), n1, i2 - i1 + 1, n2 * n3, cudaMemcpyDeviceToDevice)
        !$acc end host_data
        if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy2D")
#else
        do k = 1, n3
           do j = 1, n2
              do i = i1, i2
                 out(pos) = in(i, j, k)
                 pos = pos + 1
              end do
           end do
        end do
#endif
     end do

     return
  end subroutine mem_split_xy_complex

  subroutine mem_merge_xy_real(in, n1, n2, n3, out, iproc, dist, decomp)

     implicit none

     integer, intent(IN) :: n1, n2, n3
     real(mytype), dimension(*), intent(IN) :: in
     real(mytype), dimension(n1, n2, n3), intent(OUT) :: out
     integer, intent(IN) :: iproc
     integer, dimension(0:iproc - 1), intent(IN) :: dist
     TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
     attributes(device) :: in
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
        pos = m * decomp%y1count + 1
#else
        pos = decomp%y1disp(m) + 1
#endif

#if defined(_GPU)
        !$acc host_data use_device(out)
        istat = cudaMemcpy2D(out(1, i1, 1), n1 * n2, in(pos), n1 * (i2 - i1 + 1), n1 * (i2 - i1 + 1), n3, cudaMemcpyDeviceToDevice)
        !$acc end host_data
        if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy2D")
#else
        do k = 1, n3
           do j = i1, i2
              do i = 1, n1
                 out(i, j, k) = in(pos)
                 pos = pos + 1
              end do
           end do
        end do
#endif
     end do

     return
  end subroutine mem_merge_xy_real

  subroutine mem_merge_xy_complex(in, n1, n2, n3, out, iproc, dist, decomp)

     implicit none

     integer, intent(IN) :: n1, n2, n3
     complex(mytype), dimension(*), intent(IN) :: in
     complex(mytype), dimension(n1, n2, n3), intent(OUT) :: out
     integer, intent(IN) :: iproc
     integer, dimension(0:iproc - 1), intent(IN) :: dist
     TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
     attributes(device) :: in
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
        pos = m * decomp%y1count + 1
#else
        pos = decomp%y1disp(m) + 1
#endif

#if defined(_GPU)
        !$acc host_data use_device(out)
        istat = cudaMemcpy2D(out(1, i1, 1), n1 * n2, in(pos), n1 * (i2 - i1 + 1), n1 * (i2 - i1 + 1), n3, cudaMemcpyDeviceToDevice)
        !$acc end host_data
        if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cudaMemcpy2D")
#else
        do k = 1, n3
           do j = i1, i2
              do i = 1, n1
                 out(i, j, k) = in(pos)
                 pos = pos + 1
              end do
           end do
        end do
#endif
     end do

     return
  end subroutine mem_merge_xy_complex
