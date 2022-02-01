!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module visu
  
  implicit none

  ! True to activate the XDMF output
  logical, save :: use_xdmf = .true.
  ! True to use the new enumeration
  logical, save :: filenamedigits = .false.
  ! output2D is defined in the input.i3d file
  !        0 for 3D output (default)
  !        1 for 2D output with X average
  !        2 for 2D output with Y average
  !        3 for 2D output with Z average
  integer, save :: output2D
  integer :: ioxdmf
  character(len=9) :: ifilenameformat = '(I3.3)'
  real, save :: tstart, tend

  character(len=*), parameter :: io_name = "solution-io"

  private
  public :: output2D, visu_init, visu_ready, visu_finalise, write_snapshot, end_snapshot, &
       write_field, io_name

contains

  !
  ! Initialize the visu module
  !
  subroutine visu_init()

    use MPI
    use param, only : ilmn, iscalar, ilast, ifirst, ioutput, istret
    use variables, only : numscalar, prec, nvisu
    use param, only : dx, dy, dz
    use decomp_2d, only : nrank, mytype, xszV, yszV, zszV, xsize, ysize, zsize
    use decomp_2d_io, only : decomp_2d_init_io, decomp_2d_open_io, decomp_2d_append_mode
    use decomp_2d_io, only : decomp_2d_register_variable

    
    implicit none

    ! Local variables
    integer :: noutput, nsnapout
    real(mytype) :: memout

    integer :: is

    ! HDD usage of visu module
    if (nrank==0) then
      noutput = (ilast - ifirst +1)/ioutput

      nsnapout = 4
      if (ilmn)         nsnapout = nsnapout + 1
      if (iscalar.ne.0) nsnapout = nsnapout + numscalar

      memout = prec * nsnapout * noutput
      if (output2D.eq.0) then
        memout = memout * xszV(1) * yszV(2) * zszV(3)
      else if (output2D.eq.1) then
        memout = memout *           yszV(2) * zszV(3)
      else if (output2D.eq.2) then
        memout = memout * xszV(1)           * zszV(3)
      else if (output2D.eq.3) then
        memout = memout * xszV(1) * yszV(2)
      endif
      write(*,*)'==========================================================='
      write(*,*)'Visu module requires ',real(memout*1e-9,4),'GB'
      write(*,*)'==========================================================='
    end if

    ! Safety check
    if (output2D < 0 .or. output2D > 3 &
        .or. (output2d == 2.and.istret /= 0)) then
      if (nrank.eq.0) write(*,*) "Visu module: incorrect value for output2D."
      call MPI_ABORT(MPI_COMM_WORLD, 0, noutput)
      stop
    endif

    call decomp_2d_init_io(io_name)

    !! Register variables
    call decomp_2d_register_variable(io_name, "ux", 1, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "uy", 1, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "uz", 1, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "pp", 1, 0, output2D, mytype)
    if (ilmn) then
       call decomp_2d_register_variable(io_name, "rho", 1, 0, output2D, mytype)
    endif
    if (iscalar.ne.0) then
       do is = 1, numscalar
          call decomp_2d_register_variable(io_name, "phi"//char(48+is), 1, 0, output2D, mytype)
       enddo
    endif
    
  end subroutine visu_init

  !
  ! Indicate visu ready for IO. This is only required for ADIOS2 backend, using MPIIO this
  ! subroutine doesn't do anything.
  ! XXX: Call after all visu initialisation (main visu + case visu)
  !
  subroutine visu_ready ()

    use decomp_2d_io, only : decomp_2d_open_io, decomp_2d_append_mode, decomp_2d_write_mode

    implicit none

    integer :: mode
    
#ifdef ADIOS2
    
    mode = decomp_2d_write_mode

    ! XXX: Currently opening BP4 files in append mode seems to corrupt data.
    ! if (.not.outloc_init) then
    !    if (irestart == 1) then
    !       !! Restarting - is the output already available to write to?
    !       inquire(file=gen_iodir_name("data", io_name), exist=dir_exists)
    !       if (dir_exists) then
    !          outloc_init = .true.
    !       end if
    !    end if
       
    !    if (.not.outloc_init) then !! Yes, yes, but check the restart check above.
    !       mode = decomp_2d_write_mode
    !    else
    !       mode = decomp_2d_append_mode
    !    end if
    !    outloc_init = .true.
    ! else
    !    mode = decomp_2d_append_mode
    ! end if

    call decomp_2d_open_io(io_name, "data", mode)
#endif
    
  end subroutine visu_ready
  
  !
  ! Finalise the visu module. When using the ADIOS2 backend this closes the IO which is held open
  ! for the duration of the simulation, otherwise does nothing.
  ! 
  subroutine visu_finalise()

    use decomp_2d_io, only : decomp_2d_close_io
    
    implicit none

#ifdef ADIOS2
    call decomp_2d_close_io(io_name, "data")
#endif
    
  end subroutine visu_finalise

  !
  ! Write a snapshot
  !
  subroutine write_snapshot(rho1, ux1, uy1, uz1, pp3, phi1, ep1, itime, num)

    use decomp_2d, only : transpose_z_to_y, transpose_y_to_x
    use decomp_2d, only : mytype, xsize, ysize, zsize
    use decomp_2d, only : nrank
    use decomp_2d_io, only : decomp_2d_start_io

    use param, only : nrhotime, ilmn, iscalar, ioutput, irestart

    use variables, only : sx, cifip6, cisip6, ciwip6, cifx6, cisx6, ciwx6
    use variables, only : sy, cifip6y, cisip6y, ciwip6y, cify6, cisy6, ciwy6
    use variables, only : sz, cifip6z, cisip6z, ciwip6z, cifz6, cisz6, ciwz6
    use variables, only : numscalar

    use var, only : pp1, ta1, di1, nxmsize
    use var, only : pp2, ppi2, dip2, ph2, nymsize
    use var, only : ppi3, dip3, ph3, nzmsize
    use var, only : npress

    use tools, only : rescale_pressure

    implicit none

    !! inputs
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ux1, uy1, uz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ep1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), nrhotime), intent(in) :: rho1
    real(mytype), dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize,npress), intent(in) :: pp3
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar), intent(in) :: phi1
    integer, intent(in) :: itime
    character(len=32), intent(out) :: num

    ! Local variables
    integer :: is
    integer :: ierr
    character(len=30) :: scname
    integer :: mode
    logical, save :: outloc_init = .false.
    logical :: dir_exists

    ! Update log file
    if (nrank.eq.0) then
      call cpu_time(tstart)
      print *,'Writing snapshots =>',itime/ioutput
    end if

#ifdef ADIOS2
    call decomp_2d_start_io(io_name, "data")
#endif
    
    ! Snapshot number
#ifndef ADIOS2
    if (filenamedigits) then
       ! New enumeration system, it works integrated with xcompact3d_toolbox
       write(num, ifilenameformat) itime
    else
       ! Classic enumeration system
       write(num, ifilenameformat) itime/ioutput
    endif
#else
    ! ADIOS2 is zero-indexed
    write(num, '(I0)') itime/ioutput - 1
#endif
    
    ! Write XDMF header
    if (use_xdmf) call write_xdmf_header(".", "snapshot", trim(num))

    ! Write velocity
    call write_field(ux1, ".", "ux", trim(num))
    call write_field(uy1, ".", "uy", trim(num))
    call write_field(uz1, ".", "uz", trim(num))

    ! Interpolate pressure
    !WORK Z-PENCILS
    call interzpv(ppi3,pp3(:,:,:,1),dip3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
    (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
    !WORK Y-PENCILS
    call transpose_z_to_y(ppi3,pp2,ph3) !nxm nym nz
    call interypv(ppi2,pp2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
            (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
    !WORK X-PENCILS
    call transpose_y_to_x(ppi2,pp1,ph2) !nxm ny nz
    call interxpv(ta1,pp1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
            nxmsize,xsize(1),xsize(2),xsize(3),1)


    ! Rescale pressure
    call rescale_pressure(ta1)

    ! Write pressure
    call write_field(ta1, ".", "pp", trim(num), .true., flush=.true.)

    ! LMN - write density
    if (ilmn) call write_field(rho1(:,:,:,1), ".", "rho", trim(num))

    ! Write scalars
    if (iscalar.ne.0) then
      do is = 1, numscalar
        write(scname,"('phi',I2.2)") is
        call write_field(phi1(:,:,:,is), ".", trim(scname), trim(num), .true.)
      enddo
    endif

  end subroutine write_snapshot

  subroutine end_snapshot(itime, num)

    use decomp_2d, only : nrank
    use decomp_2d_io, only : decomp_2d_end_io
    use param, only : istret, xlx, yly, zlz
    use variables, only : nx, ny, nz, beta
    use var, only : dt,t

    implicit none

    integer, intent(in) :: itime
    character(len=32), intent(in) :: num

    character(len=32) :: fmt2, fmt3, fmt4
    integer :: is
    integer :: ierr
    
    ! Write XDMF footer
    if (use_xdmf) call write_xdmf_footer()

    ! Add metadata if XDMF is not used
    if (.not. use_xdmf) then
      if (nrank==0) then

        write(fmt2,'("(A,I16)")')
        write(fmt3,'("(A,F16.4)")')
        write(fmt4,'("(A,F16.12)")')

        open(newunit=is,file="./data/snap"//trim(num)//".ini",action='write',status='replace')
        write(is,'(A)')'[domain]'
        write(is,fmt2) 'nx=      ',nx
        write(is,fmt2) 'ny=      ',ny
        write(is,fmt2) 'nz=      ',nz
        write(is,fmt2) 'istret=  ',istret
        write(is,fmt4) 'beta=    ',beta
        write(is,fmt3) 'Lx=      ',xlx
        write(is,fmt3) 'Ly=      ',yly
        write(is,fmt3) 'Lz=      ',zlz
        write(is,'(A)')'[time]'
        write(is,fmt2) 'itime=   ',itime
        write(is,fmt3) 'dt=      ',dt
        write(is,fmt3) 't =      ',t
        close(is)

      endif
    endif

#ifdef ADIOS2
    call decomp_2d_end_io(io_name, "data")
#endif
    
    ! Update log file
    if (nrank.eq.0) then
      call cpu_time(tend)
      write(*,'(" Time for writing snapshots (s): ",F12.8)') tend-tstart
    endif

  end subroutine end_snapshot

  !
  ! Write the header of the XDMF file
  ! Adapted from https://github.com/fschuch/Xcompact3d/blob/master/src/visu.f90
  !
  subroutine write_xdmf_header(pathname, filename, num)

    use variables, only : nvisu, yp
    use param, only : dx,dy,dz,istret
    use decomp_2d, only : mytype, nrank, xszV, yszV, zszV, ystV

    implicit none

    ! Arguments
    character(len=*), intent(in) :: pathname, filename, num

    ! Local variables
    integer :: i,k
    real(mytype) :: xp(xszV(1)), zp(zszV(3))

    if (nrank.eq.0) then
      OPEN(newunit=ioxdmf,file="./data/"//gen_snapshotname(pathname, filename, num, "xdmf"))

      write(ioxdmf,'(A22)')'<?xml version="1.0" ?>'
      write(ioxdmf,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
      write(ioxdmf,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
      write(ioxdmf,*)'<Domain>'
      if (istret.ne.0) then
        do i=1,xszV(1)
          xp(i) = real(i-1,mytype)*dx*nvisu
        enddo
        do k=1,zszV(3)
          zp(k) = real(k-1,mytype)*dz*nvisu
        enddo
        write(ioxdmf,*)'    <Topology name="topo" TopologyType="3DRectMesh"'
        if (output2D.eq.0) then
          write(ioxdmf,*)'        Dimensions="',zszV(3),yszV(2),xszV(1),'">'
        else if (output2D.eq.1) then
          write(ioxdmf,*)'        Dimensions="',zszV(3),yszV(2),1,'">'
        else if (output2D.eq.2) then
          write(ioxdmf,*)'        Dimensions="',zszV(3),1,xszV(1),'">'
        else if (output2D.eq.3) then
          write(ioxdmf,*)'        Dimensions="',1,yszV(2),xszV(1),'">'
        endif
        write(ioxdmf,*)'    </Topology>'
        write(ioxdmf,*)'    <Geometry name="geo" Type="VXVYVZ">'
        if (output2D.ne.1) then
          write(ioxdmf,*)'        <DataItem Dimensions="',xszV(1),'" NumberType="Float" Precision="4" Format="XML">'
          write(ioxdmf,*)'        ',xp(:)
        else
          write(ioxdmf,*)'        <DataItem Dimensions="1" NumberType="Float" Precision="4" Format="XML">'
          write(ioxdmf,*)'        ',xp(1)
        endif
        write(ioxdmf,*)'        </DataItem>'
        if (output2D.ne.2) then
          write(ioxdmf,*)'        <DataItem Dimensions="',yszV(2),'" NumberType="Float" Precision="4" Format="XML">'
          write(ioxdmf,*)'        ',yp(ystV(1)::nvisu)
        else
          write(ioxdmf,*)'        <DataItem Dimensions="1" NumberType="Float" Precision="4" Format="XML">'
          write(ioxdmf,*)'        ',yp(1)
        endif
        write(ioxdmf,*)'        </DataItem>'
        if (output2D.ne.3) then
          write(ioxdmf,*)'        <DataItem Dimensions="',zszV(3),'" NumberType="Float" Precision="4" Format="XML">'
          write(ioxdmf,*)'        ',zp(:)
        else
          write(ioxdmf,*)'        <DataItem Dimensions="1" NumberType="Float" Precision="4" Format="XML">'
          write(ioxdmf,*)'        ',zp(1)
        endif
        write(ioxdmf,*)'        </DataItem>'
        write(ioxdmf,*)'    </Geometry>'
      else
        write(ioxdmf,*)'    <Topology name="topo" TopologyType="3DCoRectMesh"'
        if (output2D.eq.0) then
          write(ioxdmf,*)'        Dimensions="',zszV(3),yszV(2),xszV(1),'">'
        else if (output2D.eq.1) then
          write(ioxdmf,*)'        Dimensions="',zszV(3),yszV(2),1,'">'
        else if (output2D.eq.2) then
          write(ioxdmf,*)'        Dimensions="',zszV(3),1,xszV(1),'">'
        else if (output2D.eq.3) then
          write(ioxdmf,*)'        Dimensions="',1,yszV(2),xszV(1),'">'
        endif
        write(ioxdmf,*)'    </Topology>'
        write(ioxdmf,*)'    <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
        write(ioxdmf,*)'        <!-- Origin -->'
        write(ioxdmf,*)'        <DataItem Format="XML" Dimensions="3">'
        write(ioxdmf,*)'        0.0 0.0 0.0'
        write(ioxdmf,*)'        </DataItem>'
        write(ioxdmf,*)'        <!-- DxDyDz -->'
        write(ioxdmf,*)'        <DataItem Format="XML" Dimensions="3">'
        if (output2D.eq.0) then
          write(ioxdmf,*)'        ',nvisu*dz,nvisu*dy,nvisu*dx
        else if (output2D.eq.1) then
          write(ioxdmf,*)'        ',dz,dy,1.
        else if (output2D.eq.2) then
          write(ioxdmf,*)'        ',dz,1.,dx
        else if (output2D.eq.3) then
          write(ioxdmf,*)'        ',1.,dy,dx
        endif
        write(ioxdmf,*)'        </DataItem>'
        write(ioxdmf,*)'    </Geometry>'
      endif
      write(ioxdmf,*)'    <Grid Name="'//num//'" GridType="Uniform">'
      write(ioxdmf,*)'        <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
      write(ioxdmf,*)'        <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
    endif
  end subroutine write_xdmf_header

  !
  ! Write the footer of the XDMF file
  ! Adapted from https://github.com/fschuch/Xcompact3d/blob/master/src/visu.f90
  !
  subroutine write_xdmf_footer()

    use decomp_2d, only : nrank

    implicit none

    if (nrank.eq.0) then
      write(ioxdmf,'(/)')
      write(ioxdmf,*)'    </Grid>'
      write(ioxdmf,*)'</Domain>'
      write(ioxdmf,'(A7)')'</Xdmf>'
      close(ioxdmf)
    endif

  end subroutine write_xdmf_footer

  !
  ! Write the given field for visualization
  ! Adapted from https://github.com/fschuch/Xcompact3d/blob/master/src/visu.f90
  !
  subroutine write_field(f1, pathname, filename, num, skip_ibm, flush)

    use mpi
    
    use var, only : ep1
    use var, only : zero, one
    use var, only : uvisu
    use param, only : iibm
    use decomp_2d, only : mytype, xsize, xszV, yszV, zszV
    use decomp_2d, only : nrank, fine_to_coarseV
    use decomp_2d_io, only : decomp_2d_write_one, decomp_2d_write_plane

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: f1
    character(len=*), intent(in) :: pathname, filename, num 
    logical, optional, intent(in) :: skip_ibm, flush

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: local_array
    logical :: mpiio, force_flush
    
    integer :: ierr

#ifndef ADIOS2
    mpiio = .true.
#else
    mpiio = .false.
#endif

    if (present(flush)) then
       force_flush = flush
    else
       force_flush = .false.
    end if
 
    if (use_xdmf) then
       if (nrank.eq.0) then
          write(ioxdmf,*)'        <Attribute Name="'//filename//'" Center="Node">'
#ifndef ADIOS2
          write(ioxdmf,*)'           <DataItem Format="Binary"'
#else
          write(ioxdmf,*)'           <DataItem Format="HDF"'
#endif
#ifdef DOUBLE_PREC
#ifdef SAVE_SINGLE
          if (output2D.eq.0) then
             write(ioxdmf,*)'            DataType="Float" Precision="4" Endian="little" Seek="0"'
          else
             write(ioxdmf,*)'            DataType="Float" Precision="8" Endian="little" Seek="0"'
          endif
#else
          write(ioxdmf,*)'            DataType="Float" Precision="8" Endian="little" Seek="0"'
#endif
#else
          write(ioxdmf,*)'            DataType="Float" Precision="4" Endian="little" Seek="0"'
#endif
          if (output2D.eq.0) then
             write(ioxdmf,*)'            Dimensions="',zszV(3),yszV(2),xszV(1),'">'
          else if (output2D.eq.1) then
             write(ioxdmf,*)'            Dimensions="',zszV(3),yszV(2),1,'">'
          else if (output2D.eq.2) then
             write(ioxdmf,*)'            Dimensions="',zszV(3),1,xszV(1),'">'
          else if (output2D.eq.3) then
             write(ioxdmf,*)'            Dimensions="',1,yszV(2),xszV(1),'">'
          endif
          write(ioxdmf,*)'              '//gen_h5path(gen_filename(pathname, filename, num, 'bin'), num)
          write(ioxdmf,*)'           </DataItem>'
          write(ioxdmf,*)'        </Attribute>'
       endif
    endif

    if ((iibm == 2) .and. .not.present(skip_ibm)) then
       local_array(:,:,:) = (one - ep1(:,:,:)) * f1(:,:,:)
    else
       local_array(:,:,:) = f1(:,:,:)
    endif
    if (output2D.eq.0) then
       if (mpiio .or. (iibm == 2) .or. force_flush) then
          !! XXX: This (re)uses a temporary array for data - need to force synchronous writes.
          uvisu = zero
          
          call fine_to_coarseV(1,local_array,uvisu)
          call decomp_2d_write_one(1,uvisu,"data",gen_filename(pathname, filename, num, 'bin'),2,io_name,&
               opt_deferred_writes=.false.)
       else
          call decomp_2d_write_one(1,f1,"data",gen_filename(pathname, filename, num, 'bin'),0,io_name)
       end if
    else
       call decomp_2d_write_plane(1,local_array,output2D,-1,"data",gen_filename(pathname, filename, num, 'bin'),io_name)
    endif

  end subroutine write_field

  function gen_snapshotname(pathname, varname, num, ext)
    character(len=*), intent(in) :: pathname, varname, num, ext
#ifndef ADIOS2
    character(len=(len(pathname) + 1 + len(varname) + 1 + len(num) + 1 + len(ext))) :: gen_snapshotname
    write(gen_snapshotname, "(A)") gen_filename(pathname, varname, num, ext)
#else
    character(len=(len(varname) + 1 + len(num) + 1 + len(ext))) :: gen_snapshotname
    write(gen_snapshotname, "(A)") varname//'-'//num//'.'//ext
#endif
  end function gen_snapshotname
  
  function gen_filename(pathname, varname, num, ext)

    character(len=*), intent(in) :: pathname, varname, num, ext
#ifndef ADIOS2
    character(len=(len(pathname) + 1 + len(varname) + 1 + len(num) + 1 + len(ext))) :: gen_filename
    write(gen_filename, "(A)") pathname//'/'//varname//'-'//num//'.'//ext
#else
    character(len=len(varname)) :: gen_filename
    write(gen_filename, "(A)") varname
#endif
    
  end function gen_filename

  function gen_h5path(filename, num)

    character(len=*), intent(in) :: filename, num
#ifndef ADIOS2
    character(len=*), parameter :: path_to_h5file = "./"
    character(len=(len(path_to_h5file) + len(filename))) :: gen_h5path
    write(gen_h5path, "(A)") path_to_h5file//filename
#else
    character(len=*), parameter :: path_to_h5file = "../data.hdf5:/Step"
    character(len=(len(path_to_h5file) + len(num) + 1+ len(filename))) :: gen_h5path
    write(gen_h5path, "(A)") path_to_h5file//num//"/"//filename
#endif
    
  end function gen_h5path
  
end module visu
