!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

PROGRAM post

  USE decomp_2d
  USE decomp_2d_io
  USE variables
  USE param
  USE var
  USE MPI
  USE post_processing

  implicit none

  integer :: code,i,j,k,is,ii

  integer :: icrfile, file1, filen, ifile, ssfile, nt, iwrn, num, ttsize
  integer :: read_phi, read_u, read_ibm     !! 3D fields
  integer :: comp_visu, comp_post
  real(8) :: tstart,t1,trank,tranksum,ttotal,tremaining,telapsed,trstart, trend
  character(30) :: filename
  character(1) :: a

  TYPE(DECOMP_INFO) :: phG,ph1,ph2,ph3,ph4

  call ft_parameter(.true.)

  CALL MPI_INIT(code)
  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call init_coarser_mesh_statS(nstat,nstat,nstat,.true.) !start from 1 == true
  call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.) !start from 1 == true
  call init_coarser_mesh_statP(nprobe,nprobe,nprobe,.true.) !start from 1 == true
  call parameter()
  call init_variables
  call schemes()

  ux1=zero; 
  uy1=zero; 
  uz1=zero;
  phi1=zero; 
  diss1=zero; 
  pre1=zero; 

  read_phi=0; read_u=0; read_ibm=0
  open(10,file='post.prm',status='unknown',form='formatted')
  read (10,'(A1)') a
  read (10,'(A1)') a
  read (10,'(A1)') a
  read (10,*) file1
  read (10,*) filen
  read (10,*) icrfile
  read (10,*) ssfile
  read (10,'(A1)') a
  read (10,'(A1)') a
  read (10,'(A1)') a
  read (10,*) comp_post
  read (10,*) comp_visu
  close(10)

  nt = (filen-file1)/icrfile+1

  if(comp_post .eq. 1) then
     read_phi=1; read_ibm=1; read_u=1
  endif

  if(comp_visu .eq. 1) then
     call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.)    !start from 1 == true
     read_phi=1; read_u=1; read_ibm=1
  endif
  if ( iscalar.eq.0) read_phi=0
  if ( ivirt  .eq.0) read_ibm=0
  if ( read_ibm .eq. 1 ) then
     call decomp_2d_read_one(1,ep1,'./data/ibm0000')
  endif

  call init_post(ep1)

  ttsize=(read_phi*numscalar+read_u*3)*nx*ny*nz
  tstart=0.;t1=0.;trank=0.;tranksum=0.;ttotal=0.
  call cpu_time(tstart)

  do ii=1, nt
     call cpu_time(t1)
     ifile = (ii-1)*icrfile+file1
     write(filename,"(I4.4)") ifile
     t=dt*real(imodulo*ifile,mytype)
     itime=imodulo*ifile
     if (nrank==0) then
        write(*,*)'--------------------------------------------'
        write(*,*)'Snapshot',ifile, t
     endif

     !READ DATA
     call cpu_time(trstart)
     if ( read_phi .eq. 1 ) then
        do is=1, numscalar
           write(filename,"('./data/phi',I1.1,I4.4)") is, ifile
           call decomp_2d_read_one(1,phi1(:,:,:,is),filename)
        enddo
        call test_scalar_min_max(phi1)
     endif
     if ( read_u .eq. 1 ) then
        write(filename,"('./data/ux',I4.4)") ifile
        call decomp_2d_read_one(1,ux1,filename)
        write(filename,"('./data/uy',I4.4)") ifile
        call decomp_2d_read_one(1,uy1,filename)
        write(filename,"('./data/uz',I4.4)") ifile
        call decomp_2d_read_one(1,uz1,filename)
        call test_speed_min_max(ux1,uy1,uz1)
     endif
     call cpu_time(trend)

     if (ivisu.ne.0) then
        if (comp_visu .eq. 1) then
           call VISU_INSTA(ux1,uy1,uz1,phi1,ep1,.True.)
        endif
     endif
     if (ipost.ne.0) then
        if (comp_post .eq. 1) then
           call postprocessing(ux1,uy1,uz1,phi1,ep1)
        endif
     endif
     call cpu_time(trank)

     telapsed = (trank-tstart)/3600.
     tremaining  =  telapsed*(nt-ii)/(ii)

     if (nrank==0) then
        write(*,*)'Time per this snapshot (s):',real(trank-t1)
        write(*,*)'Reading speed (MB/s)',real((ttsize*1e-6)/(trend-trstart),4)
        write(*,"(' Remaining time:',I8,' h ',I2,' min')"), int(tremaining), int((tremaining-int(tremaining))*60.)
        write(*,"(' Elapsed time:',I8,' h ',I2,' min')"), int(telapsed), int((telapsed-int(telapsed))*60.)
     endif
  enddo

  call cpu_time(trank)
  ttotal=trank-tstart

  if (nrank==0) then
     write(*,*)'==========================================================='
     write(*,*)''
     write(*,*)'Post-processing finished successfully!'
     write(*,*)''
     write(*,*)'2DECOMP with p_row*p_col=',p_row,p_col
     write(*,*)''
     write(*,*)'nx*ny*nz=',nx*ny*nz
     write(*,*)'nx,ny,nz=',nx,ny,nz
     write(*,*)'dx,dy,dz=',dx,dy,dz
     write(*,*)''
     write(*,*)'Averaged time per snapshot (s):',real(ttotal/nt,4)
     write(*,*)'Total wallclock (s):',real(ttotal,4)
     write(*,*)'Total wallclock (m):',real(ttotal/60.,4)
     write(*,*)'Total wallclock (h):',real(ttotal/3600.,4)
     write(*,*)'Total wallclock (d):',real(ttotal*1.1574e-5,4)
     write(*,*)''
  endif

  call decomp_2d_finalize
  CALL MPI_FINALIZE(code)
end PROGRAM post
