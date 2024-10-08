!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause
module particle

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank,nproc
  use mptool
  use constants, only : pi

  implicit none

  interface mclean
    module procedure mclean_mytype
    module procedure mclean_particle
  end interface mclean 
  !
  interface msize
    module procedure size_particle
    module procedure size_integer
  end interface msize
  !
  interface pa2a
    module procedure pa2a_particle
  end interface pa2a
  !
  interface mextend
     module procedure extend_particle
  end interface mextend
  !
  ! particles
  type partype

    private

    real(mytype) :: rho,mas,dimeter,re,vdiff,qe
    real(mytype) :: x(3),v(3),vf(3),f(3),b(3)
    real(mytype),allocatable,dimension(:,:) :: dx,dv
    integer :: id,rankinn,rank2go
    logical :: swap,new
    !+------------------+------------------------------------------+
    !|              rho | density                                  |
    !|              mas | mass                                     |
    !|          dimeter | dimeter                                  |
    !|               re | particle Reynolds number.                |
    !|            vdiff | velocity difference between fluid and    |
    !|                  | particle                                 |
    !|               qe | Electric charge                          |
    !|                x | spatial coordinates of particle          |
    !|                v | velocity of particle                     |
    !|               vf | velocity of fluids                       |
    !|                f | force on particle                        |
    !|                b | magnetic field on particle               |
    !|               dx | gradient of x,y,z to time, used for      |
    !|                  | temporal integration.                    |
    !|               dv | gradient of u,v,w to time, used for      |
    !|                  | temporal integration.                    |
    !|               id | the id of particle                       |
    !|          rankinn | the mpi rank which the particle is in    |
    !|          rank2go | the mpi rank which the particle will be  |
    !|             swap | a indicator to exchange particle         |
    !|              new | a indicator to a new particle            |
    !+------------------+------------------------------------------+
    
    contains

    private

    procedure :: init  => init_one_particle
    procedure :: reset => reset_one_particle
    procedure :: rep   => particle_reynolds_cal
    procedure :: force => particle_force_cal
    !
  end type partype

  integer :: numparticle,n_particle
  real(mytype) :: particletime,sub_time_step
  type(partype),allocatable,target :: partpack(:)
  real(mytype),allocatable,dimension(:) :: lxmin,lxmax,lymin,lymax,lzmin,lzmax
  !+------------------+-----------------------------------------------+
  !|      numparticle | number of particles in the domain             |
  !|       n_particle | total number of particles.                    |
  !|     particletime | current particle time .                       |
  !|    sub_time_step | particle advancing time step.                 |
  !|         partpack | pack of particles.                            |
  !|            lxmin |                                               |
  !|            lxmax |                                               |
  !|            lymin |                                               |
  !|            lymax |                                               |
  !|            lzmin |                                               |
  !|            lzmax | range of local domain.                        |
  !|                  |                                               |
  !+------------------+-----------------------------------------------+

  private

  public :: local_domain_size,particle_init,intt_particel,            &
            visu_particle,particle_checkpoint,n_particle

  contains

  !+-------------------------------------------------------------------+
  !| This subroutine is used to init a particle.                       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-Jun-2022  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine init_one_particle(pa)
    !
    use param, only : ntime
    !
    class(partype),target :: pa
    
    pa%swap=.false.
    pa%new =.true.

    pa%rankinn=nrank
    pa%rank2go=nrank

    pa%x=0.0; pa%v=0.0

    allocate(pa%dx(1:3,2:ntime),pa%dv(1:3,ntime))

    pa%dx=0.0
    pa%dv=0.0

    pa%dimeter=1.d-1
    pa%rho=10._mytype

    pa%qe=1.d-5

  end subroutine init_one_particle
  !+-------------------------------------------------------------------+
  !| The end of the subroutine initmesg.                               |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| This subroutine is used to reset a particle.                      |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-Jun-2022  | Created by J. Fang @ Imperial College              |
  !+-------------------------------------------------------------------+
  subroutine reset_one_particle(pa)
    !
    class(partype),target :: pa
    !
    pa%swap=.false.
    pa%new =.true.
    !
    pa%rankinn=nrank
    pa%rank2go=nrank
    !
    pa%x=0.0; pa%v=0.0
    !
    pa%dx=0.0
    pa%dv=0.0
    !
  end subroutine reset_one_particle
  !+-------------------------------------------------------------------+
  !| The end of the subroutine reset_one_particle.                     |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| This subroutine is used to calculate particle Reynolds number     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-Nov-2022  | Created by J. Fang @ Imperial College              |
  !+-------------------------------------------------------------------+
  subroutine particle_reynolds_cal(pa)
    !
    use param, only : re
    !
    class(partype),target :: pa
    !
    pa%vdiff = sqrt( (pa%vf(1)-pa%v(1))**2 + &
                     (pa%vf(2)-pa%v(2))**2 + &
                     (pa%vf(3)-pa%v(3))**2 )
    !
    pa%re = pa%dimeter*pa%vdiff*re
    !
  end subroutine particle_reynolds_cal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine particle_reynolds_cal.                  |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| This subroutine is used to calculate froce acting on a particle   |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-Nov-2022  | Created by J. Fang @ Imperial College              |
  !+-------------------------------------------------------------------+
  subroutine particle_force_cal(pa)
    !
    use param, only : re
    !
    class(partype),target :: pa
    !
    real(mytype) :: varc
    ! 
    varc=18.d0/(pa%rho*pa%dimeter**2*re)*(1.d0+0.15d0*pa%re**0.687d0)
    !
    pa%f(:) = varc*(pa%vf(:)-pa%v(:))
    !
  end subroutine particle_force_cal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine particle_reynolds_cal.                  |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| This subroutine is to report time cost for particles.             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 18-06-2022  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine partcle_report
    !
    integer :: ttal_particle
    !
    ! ttal_particle=psum(numparticle)
    !
    if(nrank==0) then
      ! write(*,*) 'Total number of particles:',ttal_particle
      ! write(*,*) 'Total time for particles :',real(part_time,4)
      ! write(*,*) '      time particles vel :',real(part_vel_time,4)
      ! write(*,*) '      time domain search :',real(part_dmck_time,4)
      ! write(*,*) '      time partical_swap :',real(part_comm_time,4)
      ! write(*,*) '           alltoall comm :',real(a2a_time,4)
      ! write(*,*) '           counting time :',real(count_time,4)
      ! write(*,*) '           table shareing:',real(table_share_time,4)
      ! write(*,*) '           data packing  :',real(data_pack_time,4)
      ! write(*,*) '           MPI Alltoall  :',real(mpi_comm_time,4)
      ! write(*,*) '           data unpacking:',real(data_unpack_time,4)
      ! write(*,*) '                  hdf5 io:',real(h5io_time,4)
    endif
    !
  end subroutine partcle_report
  !+-------------------------------------------------------------------+
  !| The end of the subroutine partcle_report.                         |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| This subroutine is to initilise particle module.                  |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-11-2021  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine particle_init
    !
    use param,     only : xlx,yly,zlz,irestart,t,dt
    use var,       only : itime
    !
    ! local data
    ! integer :: i,j,k,p
    ! real(mytype) :: dx,dy,dz
    !
    ! generate random particles within the local domain
    call particle_gen_random(partpack,numparticle)

    call partical_domain_check('bc_tgv',partpack)

    particletime=t
    sub_time_step=dt

  end subroutine particle_init
  !+-------------------------------------------------------------------+
  ! The end of the subroutine particle_init                            |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| This subroutine is to integrate particle coordinates in time.     |
  !+-------------------------------------------------------------------+
  !| only Euler scheme is used for now.                                |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16-11-2021  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine intt_particel(ux1,uy1,uz1,time1)
    !
    use variables 
    use param
    use decomp_2d, only : xsize
    !
    ! arguments
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1
    real(mytype),intent(in) :: time1
    !
    ! local data 
    integer :: p,psize,jpart,npart,total_num_part
    integer,save :: old_num_part=0
    type(partype),pointer :: pa
    real(mytype),save :: time0,tsro
    !
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: uxp,uyp,uzp
    real(mytype),allocatable,dimension(:,:,:),save :: ux0,uy0,uz0
    real(mytype),allocatable :: xcor(:,:),dxco(:,:,:),vcor(:,:),dvco(:,:,:)
    !
    logical,save :: firstcal=.true.
    !
    if(firstcal) then
      !
      allocate( ux0(xsize(1),xsize(2),xsize(3)), &
                uy0(xsize(1),xsize(2),xsize(3)), &
                uz0(xsize(1),xsize(2),xsize(3)) )
      !
      ux0=ux1
      uy0=uy1
      uz0=uz1
      !
      time0=time1
      !
      tsro=sub_time_step/dt
      !
      firstcal=.false.
      !
      return
      !
    endif
    !
    do while(particletime<time1)
      !
      particletime=particletime+sub_time_step
      !
      ! print*,'time0=',time0,'time1=',time1,'particletime=',particletime
      !
      uxp=linintp(time0,time1,ux0,ux1,particletime)
      uyp=linintp(time0,time1,uy0,uy1,particletime)
      uzp=linintp(time0,time1,uz0,uz1,particletime)
      !
      call particle_force(uxp,uyp,uzp)
      !
      allocate(xcor(3,numparticle),dxco(3,numparticle,ntime))
      allocate(vcor(3,numparticle),dvco(3,numparticle,ntime))
      !
      psize=msize(partpack)
      !
      npart=0
      do jpart=1,psize
        !
        pa=>partpack(jpart)
        !
        if(pa%new) cycle
        !
        npart=npart+1
        !
        xcor(:,npart)=pa%x(:)

        dxco(:,npart,1)=pa%v(:)

        dxco(:,npart,2:ntime)=pa%dx(:,2:ntime)

        vcor(:,npart)=pa%v(:)
        !
        dvco(:,npart,1)=pa%f(:)
        dvco(:,npart,2:ntime)=pa%dv(:,2:ntime)
        !
        if(npart==numparticle) exit
        !
      enddo
      !
      if (itimescheme.eq.1) then
         !>>> Euler
         vcor=gdt(itr)*tsro*dvco(:,:,1)+vcor
         !
         xcor=gdt(itr)*tsro*dxco(:,:,1)+xcor
         !
      elseif(itimescheme.eq.2) then
         !>>> Adam-Bashforth second order (AB2)
         !
         if(itime.eq.1 .and. irestart.eq.0) then
           ! Do first time step with Euler
           vcor=gdt(itr)*tsro*dvco(:,:,1)+vcor
           !
           xcor=gdt(itr)*tsro*dxco(:,:,1)+xcor
         else
           vcor=adt(itr)*dvco(:,:,1)+bdt(itr)*dvco(:,:,2)+vcor
           !
           xcor=adt(itr)*dxco(:,:,1)+bdt(itr)*dxco(:,:,2)+xcor
         endif
         dvco(:,:,2)=dvco(:,:,1)
         dxco(:,:,2)=dxco(:,:,1)
         !
      elseif(itimescheme.eq.3) then
         !>>> Adams-Bashforth third order (AB3)
         !
         ! Do first time step with Euler
         if(itime.eq.1.and.irestart.eq.0) then
            vcor=dt*tsro*dvco(:,:,1)+vcor
            !
            xcor=dt*tsro*dxco(:,:,1)+xcor
         elseif(itime.eq.2.and.irestart.eq.0) then
            ! Do second time step with AB2
            vcor=onepfive*dt*tsro*dvco(:,:,1)-half*dt*tsro*dvco(:,:,2)+vcor
            dvco(:,:,3)=dvco(:,:,2)
            !
            xcor=onepfive*dt*tsro*dxco(:,:,1)-half*dt*tsro*dxco(:,:,2)+xcor
            dxco(:,:,3)=dxco(:,:,2)
         else
            ! Finally using AB3
            vcor=adt(itr)*tsro*dvco(:,:,1)+bdt(itr)*tsro*dvco(:,:,2)+cdt(itr)*tsro*dvco(:,:,3)+vcor
            dvco(:,:,3)=dvco(:,:,2)
            !
            xcor=adt(itr)*tsro*dxco(:,:,1)+bdt(itr)*tsro*dxco(:,:,2)+cdt(itr)*tsro*dxco(:,:,3)+xcor
            dxco(:,:,3)=dxco(:,:,2)
         endif
         dvco(:,:,2)=dvco(:,:,1)
         dxco(:,:,2)=dxco(:,:,1)
         !
      elseif(itimescheme.eq.5) then
         !>>> Runge-Kutta (low storage) RK3
         if(itr.eq.1) then
            vcor=gdt(itr)*tsro*dvco(:,:,1)+vcor
            xcor=gdt(itr)*tsro*dxco(:,:,1)+xcor
         else
            vcor=adt(itr)*tsro*dvco(:,:,1)+bdt(itr)*tsro*dvco(:,:,2)+vcor
            xcor=adt(itr)*tsro*dxco(:,:,1)+bdt(itr)*tsro*dxco(:,:,2)+xcor
         endif
         dvco(:,:,2)=dvco(:,:,1)
         dxco(:,:,2)=dxco(:,:,1)
         !
      endif
      ! !
      ! put back from local array to particle array
      npart=0
      do jpart=1,psize
        !
        pa=>partpack(jpart)
        !
        if(pa%new) cycle
        !
        npart=npart+1
        !
        if(isnan(vcor(1,npart)) .or. isnan(vcor(2,npart)) .or. isnan(vcor(3,npart))) then
          print*, '!! particle velocity NaN @ intt_particel'
          print*, pa%v(:),pa%f(:)
          stop
        endif
        !
        pa%v(:)=vcor(:,npart)
        pa%x(:)=xcor(:,npart)
        ! pa%x(3)=0.d0 ! 2d x-y computation
        !
        pa%dv(:,2:ntime)=dvco(:,npart,2:ntime)
        pa%dx(:,2:ntime)=dxco(:,npart,2:ntime)
        !
        if(npart==numparticle) exit
        !
      enddo
      !
      deallocate(xcor,dxco,vcor,dvco)
      !
      ! call partical_domain_check('bc_tgv')
      call partical_domain_check('bc_tgv',partpack)
      !
      call partical_swap
      !
      total_num_part=psum(numparticle)
      !
      if(nrank==0) then
        if(total_num_part.ne.old_num_part) then
          print*,' ** number of particles changes from ',old_num_part,'->',total_num_part
        endif
        old_num_part=total_num_part
      endif
      !
    enddo
    !
    ux0=ux1
    uy0=uy1
    uz0=uz1
    !
    time0=time1

    !
  end subroutine intt_particel
  !+-------------------------------------------------------------------+
  ! The end of the subroutine intt_particel                            |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| This subroutine is to visulise particle data based on binary data |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 03-10-2024  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine visu_particle(itime)

    use param, only : ioutput
    use visu,  only : gen_snapshotname,int_to_str

    integer, intent(in) :: itime

    integer :: num,j,iovp,psize,total_num_part,p
    character(len=64) :: file2write

    real(mytype),allocatable :: xyz(:,:)

    ! Snapshot number
    num = itime/ioutput 

    if(numparticle>0) then
      !
      allocate(xyz(3,numparticle))
      !
      psize=msize(partpack)

      j=0
      do p=1,psize
        !
        if(partpack(p)%new) cycle
        !
        j=j+1

        xyz(1:3,j)  =partpack(p)%x(1:3)

      enddo

    endif

    total_num_part=psum(numparticle)
    !
    file2write='particle_coordinates_'//int_to_str(num)//'.bin'

    call pwrite('./data/'//trim(file2write),xyz)

    if(nrank==0) then

      OPEN(newunit=iovp,file=gen_snapshotname('./data/', 'snapshot_particle', num, "xdmf"))

      write(iovp,'(A)')'<?xml version="1.0" ?>'
      write(iovp,'(A)')'<Xdmf Version="3.0" xmlns:xi="http://www.w3.org/2001/XInclude">'
      write(iovp,'(A)')'  <Domain>'
      write(iovp,'(A)')'    <Grid GridType="Collection" CollectionType="Temporal">'

      write(iovp,'(A)')'      <Grid Name="ParticleData_t0" GridType="Uniform">'
      write(iovp,'(A,I0,A)')'        <Topology TopologyType="Polyvertex" NodesPerElement="',total_num_part,'"/>'
      write(iovp,'(A)')'        <Geometry GeometryType="XYZ">'
      write(iovp,'(A,I0,A)')'          <DataItem Format="binary" Dimensions="',total_num_part,' 3 " NumberType="Float" Precision="8">'
      write(iovp,'(A,A,A)')'            ',trim(file2write),' </DataItem>'
      write(iovp,'(A)')'        </Geometry>'
      write(iovp,'(A)')'      </Grid>'

      write(iovp,'(A)')'    </Grid>'
      write(iovp,'(A)')'  </Domain>'
      write(iovp,'(A)')'</Xdmf>'

      close(iovp)

      print*,'<< ./data/visu_particle.xdmf'

    endif

  end subroutine visu_particle
  !+-------------------------------------------------------------------+
  ! The end of the subroutine visu_particle                            |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| this function is to write a checkpoint file for restart           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 08-Oct-2024  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine particle_checkpoint(mode)

    character(len=*),intent(in) :: mode
    
    character(len=20) :: resfile

    real(mytype),allocatable,dimension(:,:) :: xp

    integer :: j,remainder

    resfile = "checkpoint-particles"

    if(mode=='write') then

      allocate(xp(1:3,numparticle))
  
      do j=1,numparticle
        xp(1:3,j)=partpack(j)%x(1:3)
      enddo
  
      call pwrite(resfile,xp)
  
      deallocate(xp)

    elseif(mode=='read') then


      numparticle=n_particle/nproc

      remainder=mod(n_particle,nproc)

      if(nrank<remainder) then
        numparticle=numparticle+1
      endif

      allocate(xp(1:3,numparticle))

      call pread(resfile,xp)

      allocate(partpack(1:numparticle))

      do j=1,numparticle
        !
        call partpack(j)%init()
        !
        partpack(j)%x(1:3)=xp(1:3,j)
        !
        partpack(j)%v(1) = 0.d0
        partpack(j)%v(2) = 0.d0
        partpack(j)%v(3) = 0.d0
        !
        partpack(j)%new  = .false.
        !
      enddo

      deallocate(xp)      

    endif

  end subroutine particle_checkpoint
  !+-------------------------------------------------------------------+
  ! The end of the subroutine particle_checkpoint                      |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| This subroutine is to get the fluids velocity on particles.       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 15-11-2021  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine field_var(ux1,uy1,uz1,B)
    !
    use MPI
    use param,     only : dx,dy,dz,istret,nclx,ncly,nclz,xlx,yly,zlz
    use variables, only : yp,ny,nz
    use decomp_2d, only : xsize,xstart,xend,update_halo
    use actuator_line_model_utils
    !
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),1:3),intent(in),optional :: B
    !
    ! local data
    integer :: jpart,npart,i,j,k,psize
    type(partype),pointer :: pa
    real(mytype) :: x1,y1,z1,x2,y2,z2
    real(mytype) :: test(xsize(1),xsize(2),xsize(3))
    real(mytype),allocatable,dimension(:,:,:) :: ux1_halo,uy1_halo,    &
                                                 uz1_halo,ux1_hal2,    &
                                                 uy1_hal2,uz1_hal2
    real(mytype),allocatable :: b_halo(:,:,:,:)
    !
    real(mytype),allocatable,save :: xx(:),yy(:),zz(:)
    logical,save :: firstcal=.true.
    !
    if(firstcal) then
      !
      allocate(xx(0:xsize(1)+1),yy(0:xsize(2)+1),zz(0:xsize(3)+1))
      !
      do i=0,xsize(1)+1
        xx(i)=real(i-1,mytype)*dx
      enddo
      do j=0,xsize(2)+1
        !
        if(j+xstart(2)-1>ny) then
          yy(j)=2.d0*yp(ny)-yp(ny-1)
        elseif(j+xstart(2)-1<1) then
          yy(j)=2.d0*yp(1)-yp(2)
        else
          yy(j)=yp(j+xstart(2)-1)
        endif
        !
      enddo
      !
      if(xsize(3)==1) then
        zz(:)=0.d0
      else
        do k=0,xsize(3)+1
          zz(k)=real((k+xstart(3)-2),mytype)*dz
        enddo
      endif
      !
      firstcal=.false.
      !
    endif
    !
    allocate( ux1_halo(1:xsize(1)+1,0:xsize(2)+1,0:xsize(3)+1),        &
              uy1_halo(1:xsize(1)+1,0:xsize(2)+1,0:xsize(3)+1),        &
              uz1_halo(1:xsize(1)+1,0:xsize(2)+1,0:xsize(3)+1) )
    ! 
    call pswap_yz(ux1,ux1_halo)
    call pswap_yz(uy1,uy1_halo)
    call pswap_yz(uz1,uz1_halo)
    !
    if(present(B)) then
      !
      allocate( b_halo(1:xsize(1)+1,0:xsize(2)+1,0:xsize(3)+1,1:3) )
      !
      call pswap_yz(b(:,:,:,1),b_halo(:,:,:,1))
      call pswap_yz(b(:,:,:,2),b_halo(:,:,:,2))
      call pswap_yz(b(:,:,:,3),b_halo(:,:,:,3))
      !
    endif
    !
    psize=msize(partpack)
    !
    npart=0
    !
    do jpart=1,psize
      !
      pa=>partpack(jpart)
      !
      if(pa%new) cycle
      !
      loopk: do k=1,xsize(3)
        !
        z1=zz(k)
        z2=zz(k+1)
        !
        do j=1,xsize(2)
          y1=yy(j)
          y2=yy(j+1)
          do i=1,xsize(1)
            x1=xx(i)
            x2=xx(i+1)
            !
            if( pa%x(1)>=x1 .and. pa%x(1)<=x2 .and. &
                pa%x(2)>=y1 .and. pa%x(2)<=y2 .and. &
                pa%x(3)>=z1 .and. pa%x(3)<=z2 ) then
              !
              ! locate the particle, do the interpolation
              ! print*,x1,x2,y1,y2,z1,z2
              pa%vf(1)=trilinear_interpolation( x1,y1,z1,            &
                                                x2,y2,z2,            &
                                            pa%x(1),pa%x(2),pa%x(3), &
                                            ux1_halo(i,j,k),     &
                                            ux1_halo(i+1,j,k),   &
                                            ux1_halo(i,j,k+1),   &
                                            ux1_halo(i+1,j,k+1), &
                                            ux1_halo(i,j+1,k),   &
                                            ux1_halo(i+1,j+1,k), &
                                            ux1_halo(i,j+1,k+1), &
                                            ux1_halo(i+1,j+1,k+1))
              !
              if(isnan(pa%vf(1))) then
                print*,' !! error in calculating vf x'
                print*,'        x1:',x1,y1,z1
                print*,'        x2:',x2,y2,z2
                print*,'      pa%x:',pa%x
                print*,'      rank:',nrank,'i,j,k',i,j,k
                print*,'  ux1_halo:',ux1_halo(i:i+1,j:j+1,k:k+1)
                stop
              endif
              !
              pa%vf(2)=trilinear_interpolation( x1,y1,z1,           &
                                               x2,y2,z2,            &
                                            pa%x(1),pa%x(2),pa%x(3),&
                                            uy1_halo(i,j,k),     &
                                            uy1_halo(i+1,j,k),   &
                                            uy1_halo(i,j,k+1),   &
                                            uy1_halo(i+1,j,k+1), &
                                            uy1_halo(i,j+1,k),   &
                                            uy1_halo(i+1,j+1,k), &
                                            uy1_halo(i,j+1,k+1), &
                                            uy1_halo(i+1,j+1,k+1)) 
              !
              if(isnan(pa%vf(2))) then
                print*,' !! error in calculating vf y'
                print*,'        x1:',x1,y1,z1
                print*,'        x2:',x2,y2,z2
                print*,'      pa%x:',pa%x
                print*,'      rank:',nrank,'i,j,k',i,j,k
                print*,'  uy1_halo:',uy1_halo(i:i+1,j:j+1,k:k+1)
                stop
              endif
              !
              pa%vf(3)=trilinear_interpolation( x1,y1,z1,            &
                                                x2,y2,z2,            &
                                            pa%x(1),pa%x(2),pa%x(3), &
                                            uz1_halo(i,j,k),     &
                                            uz1_halo(i+1,j,k),   &
                                            uz1_halo(i,j,k+1),   &
                                            uz1_halo(i+1,j,k+1), &
                                            uz1_halo(i,j+1,k),   &
                                            uz1_halo(i+1,j+1,k), &
                                            uz1_halo(i,j+1,k+1), &
                                            uz1_halo(i+1,j+1,k+1)) 
              !
              if(isnan(pa%vf(3))) then
                print*,' !! error in calculating vf z'
                print*,'        x1:',x1,y1,z1
                print*,'        x2:',x2,y2,z2
                print*,'      pa%x:',pa%x
                print*,'      rank:',nrank,'i,j,k',i,j,k
                print*,'  uz1_halo:',uz1_halo(i:i+1,j:j+1,k:k+1)
                stop
              endif
              !
              !
              if(present(B)) then
                !
                pa%b(1)=trilinear_interpolation( x1,y1,z1,            &
                                                  x2,y2,z2,            &
                                              pa%x(1),pa%x(2),pa%x(3), &
                                              b_halo(i,j,k,1),     &
                                              b_halo(i+1,j,k,1),   &
                                              b_halo(i,j,k+1,1),   &
                                              b_halo(i+1,j,k+1,1), &
                                              b_halo(i,j+1,k,1),   &
                                              b_halo(i+1,j+1,k,1), &
                                              b_halo(i,j+1,k+1,1), &
                                              b_halo(i+1,j+1,k+1,1))
                pa%b(2)=trilinear_interpolation( x1,y1,z1,            &
                                                 x2,y2,z2,            &
                                              pa%x(1),pa%x(2),pa%x(3),&
                                              b_halo(i,j,k,2),        &
                                              b_halo(i+1,j,k,2),   &
                                              b_halo(i,j,k+1,2),   &
                                              b_halo(i+1,j,k+1,2), &
                                              b_halo(i,j+1,k,2),   &
                                              b_halo(i+1,j+1,k,2), &
                                              b_halo(i,j+1,k+1,2), &
                                              b_halo(i+1,j+1,k+1,2) )
                pa%b(3)=trilinear_interpolation( x1,y1,z1,            &
                                                  x2,y2,z2,            &
                                              pa%x(1),pa%x(2),pa%x(3), &
                                              b_halo(i,j,k,3),     &
                                              b_halo(i+1,j,k,3),   &
                                              b_halo(i,j,k+1,3),   &
                                              b_halo(i+1,j,k+1,3), &
                                              b_halo(i,j+1,k,3),   &
                                              b_halo(i+1,j+1,k,3), &
                                              b_halo(i,j+1,k+1,3), &
                                              b_halo(i+1,j+1,k+1,3)) 
              if(isnan(pa%b(1)) .or. isnan(pa%b(2)) .or. isnan(pa%b(3))) then
                print*,' !! error in calculating b'
                print*,'        x1:',x1,y1,z1
                print*,'        x2:',x2,y2,z2
                print*,'      pa%x:',pa%x
                print*,'      rank:',nrank,'i,j,k',i,j,k
                print*,'  b_halo_x:',b_halo(i:i+1,j:j+1,k:k+1,1)
                print*,'  b_halo_y:',b_halo(i:i+1,j:j+1,k:k+1,2)
                print*,'  b_halo_z:',b_halo(i:i+1,j:j+1,k:k+1,3)
                stop
              endif
              !
              !
              endif
              !
              exit loopk
              !
            endif
            !
          enddo
        enddo
      enddo loopk
      !
      npart=npart+1
      !
      if(npart==numparticle) exit
      !
    enddo
    !
  end subroutine field_var
  !+-------------------------------------------------------------------+
  ! The end of the subroutine field_var                                |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to swap particle infomation between ranks      |
  !+-------------------------------------------------------------------+
  !| tecplot format for now                                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16-06-2022  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine partical_swap
    !
    use param,     only : nclx,ncly,nclz,xlx,yly,zlz
    !
    ! local data 
    integer :: p,psize,jrank,jpart,npart,n,newsize
    type(partype),pointer :: pa
    integer :: nsend(0:nproc-1),nrecv(0:nproc-1),nsend_total
    !+------------+--------------------------------+
    !| nsendtable | the table to record how much   |
    !|            | particle to send to a rank     |
    !+------------+--------------------------------+
    integer :: pr(0:nproc-1,1:numparticle)
    type(partype),allocatable :: pa2send(:),pa2recv(:)
    real(mytype) :: timebeg,tvar1,tvar11,tvar2,tvar3,tvar4
    !
    psize=msize(partpack)
    !
    nsend=0
    !
    n=0
    pr=0
    npart=0
    !
    do jpart=1,psize
      !
      pa=>partpack(jpart)
      !
      if(pa%new) cycle
      !
      ! to find out how many particle to send to which ranks
      if(pa%swap) then
        !
        n=n+1
        !
        nsend(pa%rank2go)=nsend(pa%rank2go)+1
        !
        pr(pa%rank2go,nsend(pa%rank2go))=jpart
        !
      endif
      !
      npart=npart+1
      !
      if(npart==numparticle) exit
      !
    enddo
    !
    nsend_total=n
    !
    ! synchronize recv table according to send table
    nrecv=ptabupd(nsend)
    !
    !
    ! to establish the buffer of storing particels about to send
    if(nsend_total>0) then
      !
      allocate(pa2send(1:nsend_total))
      !
      n=0
      do jrank=0,nproc-1
        !
        do jpart=1,nsend(jrank)
          !
          n=n+1
          !
          p=pr(jrank,jpart)
          !
          pa2send(n)=partpack(p)
          !
          call partpack(p)%reset()
          !
          numparticle=numparticle-1
          !
        enddo
        !
      enddo
      !
    endif 
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! swap particle among ranks
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call pa2a(pa2send,pa2recv,nsend,nrecv)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end of swap particle among ranks
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    call particle_add(partpack,pa2recv,n)
    !
    numparticle=numparticle+n
    !
  end subroutine partical_swap
  !+-------------------------------------------------------------------+
  ! The end of the subroutine partical_swap                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to add particles to the current particle arrary|
  !+-------------------------------------------------------------------+
  !| tecplot format for now                                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-06-2022  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine particle_add(particle_cur,particle_new,num_part_incr)
    !
    ! arguments
    type(partype),intent(inout),allocatable,target :: particle_cur(:)
    type(partype),intent(in),allocatable :: particle_new(:)
    integer,intent(out) :: num_part_incr
    !
    ! local data
    integer :: psize,newsize,n,jpart
    type(partype),pointer :: pa
    !
    psize=msize(particle_cur)
    !
    ! now add the received particle in to the array, dynamically
    if(numparticle+msize(particle_new)>psize) then
      !
      ! expand the particle array
      newsize=max(numparticle+msize(particle_new),numparticle+100)
      !
      call mextend(particle_cur,newsize)
      !
      psize=newsize
    endif
    !
    n=0
    do jpart=1,psize
      !
      ! print*,nrank,'|',jpart
      pa=>particle_cur(jpart)
      !
      ! the particle is free for re-assigning
      if(pa%new) then
        !
        if(n>=msize(particle_new)) exit
        !
        n=n+1
        !
        pa=particle_new(n)
        pa%new=.false.
        !
      endif
      !
    enddo
    !
    ! print*,nrank,'|',n,newsize
    num_part_incr=n
    !
  end subroutine particle_add
  !+-------------------------------------------------------------------+
  ! The end of the subroutine particle_add                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate particle array with random         |
  !| locations in a local domain.                                      |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 03-10-2024  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine particle_gen_random(particle_new,particle_size)
    !
    ! arguments
    type(partype),intent(out),allocatable :: particle_new(:)
    integer,intent(out) :: particle_size
    !
    ! local data
    integer :: p,j,local_size,total_particles
    real(mytype) :: dx,dy,dz,x,y,z,ran1,ran2,ran3
    !
    total_particles=512

    local_size=total_particles/nproc
    !
    allocate(particle_new(1:local_size))
    !
    p=0
    do j=1,local_size
      !
      call random_number(ran1)
      call random_number(ran2)
      call random_number(ran3)
      !
      x=(lxmax(nrank)-lxmin(nrank))*ran1+lxmin(nrank)
      y=(lymax(nrank)-lymin(nrank))*ran2+lymin(nrank)
      z=(lzmax(nrank)-lzmin(nrank))*ran3+lzmin(nrank)
      !
      ! print*,x,y,z
      ! print*,(y>=lymin(nrank) .and. y<=lymax(nrank))
      !
      if( x>=lxmin(nrank) .and. x<=lxmax(nrank) .and. &
          y>=lymin(nrank) .and. y<=lymax(nrank) .and. &
          z>=lzmin(nrank) .and. z<=lzmax(nrank) ) then
        !
        p=p+1
        !
        call particle_new(p)%init()
        !
        particle_new(p)%x(1)=x
        particle_new(p)%x(2)=y
        particle_new(p)%x(3)=z
        !
        particle_new(p)%v(1)   =0.d0
        particle_new(p)%v(2)   =0.d0
        particle_new(p)%v(3)   =0.d0
        !
        particle_new(p)%new=.false.
        !
      endif
      !
    enddo
    !
    call mclean(particle_new,p)
    !
    particle_size=p

    n_particle=psum(particle_size)

    if(nrank==0) print*,'** random particles are generated'
    !
  end subroutine particle_gen_random
  !+-------------------------------------------------------------------+
  !| The end of the subroutine particle_gen_random.                    |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| This subroutine is to calcualte the size and range of local domain|
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 15-11-2021  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine local_domain_size
    !
    use param,     only : dx,dy,dz,istret
    use variables, only : yp,ny
    use decomp_2d, only : xstart,xend
    use actuator_line_model_utils
    !
    integer :: nyr,nzr,jrank
    !
    allocate( lxmin(0:nproc-1),lxmax(0:nproc-1),  &
              lymin(0:nproc-1),lymax(0:nproc-1),  &
              lzmin(0:nproc-1),lzmax(0:nproc-1)   )
    !
    lxmin=0.d0; lxmax=0.d0
    lymin=0.d0; lymax=0.d0
    lzmin=0.d0; lzmax=0.d0
    !
    lxmin(nrank)=real(0,mytype)*dx
    lxmax(nrank)=real(xend(1),mytype)*dx
    !
    if (istret==0) then
      lymin(nrank)=real(xstart(2)-1,mytype)*dy
      lymax(nrank)=real(xend(2),mytype)*dy
    else
      lymin(nrank)=yp(xstart(2))
      nyr=min(ny,xend(2)+1)
      lymax(nrank)=yp(nyr)
    endif
    !
    lzmin(nrank)=real((xstart(3)-1),mytype)*dz
    nzr=xend(3)
    lzmax(nrank)=real(nzr,mytype)*dz
    !
    lxmin=psum(lxmin); lxmax=psum(lxmax)
    lymin=psum(lymin); lymax=psum(lymax)
    lzmin=psum(lzmin); lzmax=psum(lzmax)
    !
     ! do jrank=0,nproc-1
     !   print*,nrank,jrank,lxmin(jrank),lxmax(jrank),lymin(jrank),lymax(jrank),lzmin(jrank),lzmax(jrank)
     ! enddo
    !
    ! call mpistop
    !
  end subroutine local_domain_size
  !+-------------------------------------------------------------------+
  ! The end of the subroutine local_domain_size                        |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| This subroutine is to check if the particle is out of domain      |
  !+-------------------------------------------------------------------+
  !| tecplot format for now                                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16-06-2022  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine partical_domain_check(mode,particle)
    !
    use param,     only : nclx,ncly,nclz,xlx,yly,zlz

    !
    character(len=*),intent(in) :: mode
    type(partype),intent(in),allocatable,target :: particle(:)
    !
    ! local data 
    integer :: jpart,npart,psize,jrank,npcanc,npcanc_totl
    type(partype),pointer :: pa

    psize=msize(particle)

    if(mode=='out_disappear') then
      !
      npart=0
      npcanc=0
      do jpart=1,psize
        !
        pa=>particle(jpart)
        !
        if(pa%new) cycle
        !
        npart=npart+1
        if(npart>numparticle) exit
        ! if the particle is out of domain, mark it and subscribe the 
        ! total number of particles
        !
        if(nclx .and. (pa%x(1)>xlx .or. pa%x(1)<0)) then
          call pa%reset()
          npcanc=npcanc+1
          cycle
        endif
        !
        if(ncly .and. (pa%x(2)>yly .or. pa%x(2)<0)) then
          call pa%reset()
          npcanc=npcanc+1
          cycle
        endif
        !
        if(nclz .and. (pa%x(3)>zlz .or. pa%x(3)<0)) then
          call pa%reset()
          npcanc=npcanc+1
          cycle
        endif
        !
        if( pa%x(1)>=lxmin(nrank) .and. pa%x(1)<lxmax(nrank) .and. &
            pa%x(2)>=lymin(nrank) .and. pa%x(2)<lymax(nrank) .and. &
            pa%x(3)>=lzmin(nrank) .and. pa%x(3)<lzmax(nrank) ) then
          continue
        else
          !
          pa%swap=.true.
          !
          do jrank=0,nproc-1
            !
            ! to find which rank the particle are moving to and 
            ! mark
            if(jrank==nrank) cycle
            !
            if( pa%x(1)>=lxmin(jrank) .and. pa%x(1)<lxmax(jrank) .and. &
                pa%x(2)>=lymin(jrank) .and. pa%x(2)<lymax(jrank) .and. &
                pa%x(3)>=lzmin(jrank) .and. pa%x(3)<lzmax(jrank) ) then
              !
              pa%rank2go=jrank
              !
              exit
              !
            endif
            !
          enddo
          !
        endif
        !
      enddo
      !
      numparticle=numparticle-npcanc
      !
      npcanc_totl=psum(npcanc)
      if(nrank==0 .and. npcanc_totl>0) print*,' ** ',npcanc_totl,        &
                                     ' particles are moving out of domain'
      !
    elseif(mode=='bc_channel') then
      !
      npart=0
      npcanc=0
      do jpart=1,psize
        !
        pa=>particle(jpart)
        !
        if(pa%new) cycle
        !
        npart=npart+1
        !
        if(npart>numparticle) exit
        ! if the particle is out of domain, mark it and subscribe the 
        ! total number of particles
        !
        if(nclx) then
          !
          if(pa%x(1)>xlx) then
            pa%x(1)=pa%x(1)-xlx
          endif
          !
          if(pa%x(1)<0) then
            pa%x(1)=pa%x(1)+xlx
          endif
          !
        endif
        !
        if(ncly) then
          !
          if(pa%x(2)>yly) then
            pa%x(2)=pa%x(2)-yly
          endif 
          !
          if(pa%x(2)<0) then
            pa%x(2)=pa%x(2)+yly
          endif
          !
        else
          !
          ! reflect particles back in the domain.
          if(pa%x(2)>yly) then
            pa%x(2)=2.d0*yly-pa%x(2)
          endif
          if(pa%x(2)<0) then
            pa%x(2)=-pa%x(2)
          endif
          !
        endif
        !
        if(nclz) then
          !
          if(pa%x(3)>zlz) then
            pa%x(3)=pa%x(3)-zlz
          endif 
          !
          if(pa%x(3)<0) then
            pa%x(3)=pa%x(3)+zlz
          endif
          !
        endif
        !
        if(pa%x(1)>xlx .or. pa%x(1)<0 .or. &
           pa%x(2)>yly .or. pa%x(2)<0 .or. &
           pa%x(3)>zlz .or. pa%x(3)<0 ) then
            print*,' !! waring, the particle still moves out of domain 1'
            print*,nrank,jpart,'x:',pa%x(:),'v:',pa%v(:),'vf:',pa%vf(:),'f:',pa%f(:)
            stop
        endif
        !
        if( pa%x(1)>=lxmin(nrank) .and. pa%x(1)<lxmax(nrank) .and. &
            pa%x(2)>=lymin(nrank) .and. pa%x(2)<lymax(nrank) .and. &
            pa%x(3)>=lzmin(nrank) .and. pa%x(3)<lzmax(nrank) ) then
          continue
        else
          !
          pa%swap=.true.
          !
          do jrank=0,nproc-1
            !
            ! to find which rank the particle are moving to and 
            ! mark
            if(jrank==nrank) cycle
            !
            if( pa%x(1)>=lxmin(jrank) .and. pa%x(1)<lxmax(jrank) .and. &
                pa%x(2)>=lymin(jrank) .and. pa%x(2)<lymax(jrank) .and. &
                pa%x(3)>=lzmin(jrank) .and. pa%x(3)<lzmax(jrank) ) then
              !
              pa%rank2go=jrank
              !
              exit
              !
            endif
            !
          enddo
          !
        endif
        !
      enddo
      !
    elseif(mode=='bc_tgv') then
      !
      npart=0
      npcanc=0
      do jpart=1,psize
        !
        pa=>particle(jpart)
        !
        if(pa%new) cycle
        !
        npart=npart+1
        !
        if(npart>numparticle) exit
        ! if the particle is out of domain, mark it and subscribe the 
        ! total number of particles
        !
        if(nclx) then
          !
          if(pa%x(1)>xlx) then
            pa%x(1)=pa%x(1)-xlx
          endif
          !
          if(pa%x(1)<0) then
            pa%x(1)=pa%x(1)+xlx
          endif
          !
        else
          stop ' 1 @ partical_domain_check'
        endif
        !
        if(ncly) then
          !
          if(pa%x(2)>yly) then
            pa%x(2)=pa%x(2)-yly
          endif 
          !
          if(pa%x(2)<0) then
            pa%x(2)=pa%x(2)+yly
          endif
          !
        else
          print*,' ** particles x: ',pa%x(:)
          stop ' 2 @ partical_domain_check'
        endif
        !
        if(nclz) then
          !
          if(pa%x(3)>zlz) then
            pa%x(3)=pa%x(3)-zlz
          endif 
          !
          if(pa%x(3)<0) then
            pa%x(3)=pa%x(3)+zlz
          endif
          !
        else
          stop ' 3 @ partical_domain_check'
        endif
        !
        if(pa%x(1)>xlx .or. pa%x(1)<0 .or. &
           pa%x(2)>yly .or. pa%x(2)<0 .or. &
           pa%x(3)>zlz .or. pa%x(3)<0 ) then
            print*,' !! waring, the particle still moving out of domain 2'
            print*,nrank,jpart,'x:',pa%x(:),'v:',pa%v(:),'vf:',pa%vf(:),'f:',pa%f(:)
            stop
        endif
        !
        if( pa%x(1)>=lxmin(nrank) .and. pa%x(1)<lxmax(nrank) .and. &
            pa%x(2)>=lymin(nrank) .and. pa%x(2)<lymax(nrank) .and. &
            pa%x(3)>=lzmin(nrank) .and. pa%x(3)<lzmax(nrank) ) then
          continue
        else
          !
          pa%swap=.true.
          !
          do jrank=0,nproc-1
            !
            ! to find which rank the particle are moving to and 
            ! mark
            if(jrank==nrank) cycle
            !
            if( pa%x(1)>=lxmin(jrank) .and. pa%x(1)<lxmax(jrank) .and. &
                pa%x(2)>=lymin(jrank) .and. pa%x(2)<lymax(jrank) .and. &
                pa%x(3)>=lzmin(jrank) .and. pa%x(3)<lzmax(jrank) ) then
              !
              pa%rank2go=jrank
              !
              exit
              !
            endif
            !
          enddo
          !
        endif
        !
      enddo
      !
    else
      !
      stop ' !! mode not defined @ partical_domain_check!!'
      !
    endif
    !
  end subroutine partical_domain_check
  !+-------------------------------------------------------------------+
  ! The end of the subroutine partical_domain_check                    |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| This subroutine is to write particles using mpi-io.               |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 03-10-2024  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine dump_particle(particles)

    ! arguments
    type(partype),allocatable,intent(in) :: particles(:)

    ! local data
    character(len=5) :: num
    integer :: psize,total_num_part,p,j,k
    integer :: rank2coll
    character(len=64) :: file2write
    logical :: fexists
    real(mytype),allocatable :: xyz(:,:)


  end subroutine dump_particle
  !+-------------------------------------------------------------------+
  ! The end of the subroutine dump_particle                            |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| This subroutine is to calculate the force on particles.           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 15-11-2021  | Created by J. Fang                                  |
  !| 01-04-2023  | Add lorentz force by J. Fang                        |
  !+-------------------------------------------------------------------+
  subroutine particle_force(ux1,uy1,uz1)
    !
    use decomp_2d, only : xsize
    use var,       only : itime
    use param,     only : re,t,itr
    use mhd,       only : Bm,stuart
    !
    ! arguments
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1
    !
    ! local data
    integer :: psize,jpart,npart
    type(partype),pointer :: pa
    real(mytype) :: varc,vard,varq(3),maxforce,maxvdiff,vdiff,re_p,bmt(3)
    !
    logical,save :: firstcal=.true.
    !
    psize=msize(partpack)
    !
    ! print*,nrank,'-',numparticle,'-',psize
    !
    npart=0
    !
    maxforce=0.d0
    maxvdiff=0.d0
    !
    call field_var(ux1,uy1,uz1)
    ! print*,nrank,'|',numparticle
    !
    do jpart=1,psize
      !
      pa=>partpack(jpart)
      !
      if(pa%new) cycle
      !
      pa%v(:) = pa%vf(:)
      pa%f(:) = 0.d0
      !
    enddo
    !
    !
  end subroutine particle_force
  !+-------------------------------------------------------------------+
  ! The end of the subroutine particle_force                           |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| this function is to retune the size of a array                    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16-Jun-2022  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  pure function size_particle(var) result(nsize)
    !
    type(partype),allocatable,intent(in) :: var(:)
    integer :: nsize
    !
    if(allocated(var)) then
      nsize=size(var)
    else
      nsize=0
    endif
    !
    return
    !
  end function size_particle
  !
  pure function size_integer(var) result(nsize)
    !
    integer,allocatable,intent(in) :: var(:)
    integer :: nsize
    !
    if(allocated(var)) then
      nsize=size(var)
    else
      nsize=0
    endif
    !
    return
    !
  end function size_integer
  !+-------------------------------------------------------------------+
  !| The end of the subroutine size_particle.                          |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| this subroutine clean superfluous elements in a array             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 07-Nov-2018  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine mclean_mytype(var,n)
    !
    ! arguments
    real(mytype),allocatable,intent(inout) :: var(:)
    integer,intent(in) :: n
    !
    ! local data
    real(mytype),allocatable :: buffer(:)
    integer :: m
    logical :: lefc
    !
    if(.not.allocated(var)) return
    !
    if(n<=0) then
      deallocate(var)
      return
    endif
    !
    ! clean
    allocate(buffer(n))
    !
    buffer(1:n)=var(1:n)
    !
    deallocate(var)
    !
    call move_alloc(buffer,var)
    !
  end subroutine mclean_mytype
  !!
  subroutine mclean_particle(var,n)
    !
    ! arguments
    type(partype),allocatable,intent(inout) :: var(:)
    integer,intent(in) :: n
    !
    ! local data
    type(partype),allocatable :: buffer(:)
    integer :: m
    logical :: lefc
    !
    if(.not.allocated(var)) return
    !
    if(n<=0) then
      deallocate(var)
      return
    endif
    !
    ! clean
    allocate(buffer(n))
    !
    buffer(1:n)=var(1:n)
    !
    deallocate(var)
    !
    call move_alloc(buffer,var)
    !
  end subroutine mclean_particle
  !+-------------------------------------------------------------------+
  !| The end of the subroutine mclean.                                 |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| this subroutine is to expand an array.                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-Jun-2022  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine extend_particle(var,n)
    !
    ! arguments
    type(partype),allocatable,intent(inout) :: var(:)
    integer,intent(in) :: n
    !
    ! local data
    type(partype),allocatable :: buffer(:)
    integer :: m,jpart
    !
    if(.not. allocated(var)) then
      allocate(var(n))
      m=0
    else
      !
      m=size(var)
      !
      call move_alloc(var, buffer)
      !
      allocate(var(n))
      var(1:m)=buffer(1:m)
      !
    endif
    !
    ! initilise newly allocated particles
    do jpart=m+1,n
      call var(jpart)%init()
    enddo
    !
    return
    !
  end subroutine extend_particle
  !+-------------------------------------------------------------------+
  !| The end of the subroutine pa2a_particle.                          |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| this subroutine is to swap particles via alltoall.                |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-Jun-2022  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine pa2a_particle(datasend,datarecv,sendtabl,recvtabl)
    !
    use mpi
    !
    ! arguments
    type(partype),allocatable,intent(in) ::  datasend(:)
    type(partype),allocatable,intent(out) :: datarecv(:)
    integer,intent(in) :: sendtabl(0:),recvtabl(0:)
    !
    ! local data
    integer :: ierr,recvsize,jrank,jpart,jc,nindsize
    integer,allocatable :: senddispls(:),recvdispls(:)
    real(mytype),allocatable :: r8send(:,:),r8resv(:,:)
    !
    integer,save :: newtype
    !
    logical,save :: firstcal=.true.
    !
    if(firstcal) then
      call mpi_type_contiguous(6,mpi_real8,newtype,ierr)
      call mpi_type_commit(newtype,ierr)
      firstcal=.false.
    endif
    !
    allocate(senddispls(0:nproc-1),recvdispls(0:nproc-1))
    !
    senddispls=0
    recvdispls=0
    do jrank=1,nproc-1
      senddispls(jrank)=senddispls(jrank-1)+sendtabl(jrank-1)
      recvdispls(jrank)=recvdispls(jrank-1)+recvtabl(jrank-1)
    enddo
    recvsize=recvdispls(nproc-1)+recvtabl(nproc-1)
    !
    nindsize=msize(datasend)
    !
    allocate(r8send(6,nindsize))
    allocate(r8resv(6,recvsize))
    !
    r8resv=0.d0
    !
    do jpart=1,nindsize
      r8send(1,jpart)=datasend(jpart)%x(1)
      r8send(2,jpart)=datasend(jpart)%x(2)
      r8send(3,jpart)=datasend(jpart)%x(3)
      r8send(4,jpart)=datasend(jpart)%v(1)
      r8send(5,jpart)=datasend(jpart)%v(2)
      r8send(6,jpart)=datasend(jpart)%v(3)
    enddo
    !
    call mpi_alltoallv(r8send, sendtabl, senddispls, newtype, &
                       r8resv, recvtabl, recvdispls, newtype, &
                       mpi_comm_world, ierr)
    !
    allocate(datarecv(recvsize))
    !
    jc=0
    do jrank=0,nproc-1
      do jpart=1,recvtabl(jrank)
        !
        jc=jc+1
        !
        call datarecv(jc)%init()
        !
        datarecv(jc)%x(1)=r8resv(1,jc)
        datarecv(jc)%x(2)=r8resv(2,jc)
        datarecv(jc)%x(3)=r8resv(3,jc)
        datarecv(jc)%v(1)=r8resv(4,jc)
        datarecv(jc)%v(2)=r8resv(5,jc)
        datarecv(jc)%v(3)=r8resv(6,jc)
        !
        ! print*,nrank,'|',datarecv(jc)%x(1),datarecv(jc)%x(2),datarecv(jc)%x(3)
        !
      enddo
    enddo
    !
  end subroutine pa2a_particle
  !+-------------------------------------------------------------------+
  !| The end of the subroutine pa2a_particle.                          |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| this subroutine is to swap data in the y-z direction.             |
  !+-------------------------------------------------------------------+
  subroutine pswap_yz(varin,varout)
    !
    use variables, only : p_row, p_col
    use decomp_2d, only : xsize
    use param,     only : nclx,ncly,nclz
    use mpi
    !
    ! arguments
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: varin
    real(mytype),intent(out) :: varout(1:xsize(1)+1,0:xsize(2)+1,0:xsize(3)+1)
    !
    ! local data
    integer :: i,j,k,n,jrk,krk,ncou,mpitag,ierr
    integer :: status(mpi_status_size) 
    integer,allocatable,save :: mrank(:,:)
    integer,save :: upper,lower,front,bback
    logical,save :: init=.true.
    real(mytype),dimension(:,:),allocatable :: sbuf1,sbuf2,rbuf1,rbuf2
    !
    if(init) then
      !
      allocate(mrank(p_row,p_col))
      !
      n=-1
      do j=1,p_row
      do k=1,p_col
        !
        n=n+1
        mrank(j,k)=n
        !
        if(nrank==n) then
          jrk=j
          krk=k
        endif
        !
      enddo
      enddo
      !
      if(jrk==1) then
        upper=mrank(jrk+1,krk)
        !
        if(ncly) then
          lower=mrank(p_row,krk)
        else
          lower=mpi_proc_null
        endif
        !
      elseif(jrk==p_row) then
        !
        if(ncly) then
          upper=mrank(1,krk)
        else
          upper=mpi_proc_null
        endif
        !
        lower=mrank(jrk-1,krk)
      else
        upper=mrank(jrk+1,krk)
        lower=mrank(jrk-1,krk)
      endif
      !
      if(p_col==1) then
        front=mpi_proc_null
        bback=mpi_proc_null
      else
        !
        if(krk==1) then
          front=mrank(jrk,krk+1)
          !
          if(nclz) then
            bback=mrank(jrk,p_col)
          else
            bback=mpi_proc_null
          endif
          !
        elseif(krk==p_col) then
          if(nclz) then
            front=mrank(jrk,1)
          else
            front=mpi_proc_null
          endif
          !
          bback=mrank(jrk,krk-1)
        else
          front=mrank(jrk,krk+1)
          bback=mrank(jrk,krk-1)
        endif
      endif
      !
      init=.false.
      !
    endif
    !
    varout(1:xsize(1),1:xsize(2),1:xsize(3))=varin(1:xsize(1),1:xsize(2),1:xsize(3))
    !
    mpitag=1000
    !
    ! send & recv in the z direction
    !
    if(xsize(3)==1) then
      varout(1:xsize(1),1:xsize(2),0)=varout(1:xsize(1),1:xsize(2),1)
      varout(1:xsize(1),1:xsize(2),2)=varout(1:xsize(1),1:xsize(2),1)
    else
      !
      allocate(sbuf1(1:xsize(1),1:xsize(2)),sbuf2(1:xsize(1),1:xsize(2)),&
               rbuf1(1:xsize(1),1:xsize(2)),rbuf2(1:xsize(1),1:xsize(2)))
      !
      if(front .ne. mpi_proc_null) then
        sbuf1(1:xsize(1),1:xsize(2))=varout(1:xsize(1),1:xsize(2),xsize(3))
      endif
      if(bback .ne. mpi_proc_null) then
        sbuf2(1:xsize(1),1:xsize(2))=varout(1:xsize(1),1:xsize(2),1)
      endif
      !
      ncou=xsize(1)*xsize(2)
      !
      ! Message passing
      call mpi_sendrecv(sbuf1,ncou,mpi_real8,front, mpitag,             &
                        rbuf1,ncou,mpi_real8,bback, mpitag,             &
                                               mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      call mpi_sendrecv(sbuf2,ncou,mpi_real8,bback, mpitag,            &
                        rbuf2,ncou,mpi_real8,front, mpitag,            &
                                               mpi_comm_world,status,ierr)
      !
      if(bback .ne. mpi_proc_null) then
        varout(1:xsize(1),1:xsize(2),0)=rbuf1(1:xsize(1),1:xsize(2))
      endif
      if(front .ne. mpi_proc_null) then
        varout(1:xsize(1),1:xsize(2),xsize(3)+1)=rbuf2(1:xsize(1),1:xsize(2))
      endif
      !
      deallocate(sbuf1,sbuf2,rbuf1,rbuf2)
      !
    endif
    !
    ! end of Message passing in the z direction
    !
    ! 
    ! send & recv in the y direction
    allocate( sbuf1(1:xsize(1),0:xsize(3)+1),                          &
              sbuf2(1:xsize(1),0:xsize(3)+1),                          &
              rbuf1(1:xsize(1),0:xsize(3)+1),                          &
              rbuf2(1:xsize(1),0:xsize(3)+1) )

    if(upper .ne. mpi_proc_null) then
      sbuf1(1:xsize(1),0:xsize(3)+1)=varout(1:xsize(1),xsize(2),0:xsize(3)+1)
    endif
    if(lower .ne. mpi_proc_null) then
      sbuf2(1:xsize(1),0:xsize(3)+1)=varout(1:xsize(1),1,0:xsize(3)+1)
    endif
    !
    ncou=xsize(1)*(xsize(3)+2)
    !
    ! Message passing
    call mpi_sendrecv(sbuf1,ncou,mpi_real8,upper, mpitag,             &
                      rbuf1,ncou,mpi_real8,lower, mpitag,             &
                                             mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    call mpi_sendrecv(sbuf2,ncou,mpi_real8,lower, mpitag,            &
                      rbuf2,ncou,mpi_real8,upper, mpitag,            &
                                             mpi_comm_world,status,ierr)
    !
    if(upper .ne. mpi_proc_null) then
      varout(1:xsize(1),xsize(2)+1,0:xsize(3)+1)=rbuf2(1:xsize(1),0:xsize(3)+1)
    endif
    if(lower .ne. mpi_proc_null) then
      varout(1:xsize(1),0,0:xsize(3)+1)=rbuf1(1:xsize(1),0:xsize(3)+1)
    endif
    !
    deallocate(sbuf1,sbuf2,rbuf1,rbuf2)
    ! send & recv in the y direction
    !
    if(nclx) then
      varout(xsize(1)+1,:,:)=varout(1,:,:)
    endif
    !
    return
    !
  end subroutine pswap_yz
  !+-------------------------------------------------------------------+
  !| The end of the subroutine pswap_yz.                               |
  !+-------------------------------------------------------------------+

end module particle
!+---------------------------------------------------------------------+
!| The end of the module particle                                      |
!+---------------------------------------------------------------------+