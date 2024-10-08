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
  ! particles
  type partype
    !
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
    !
    contains
    !
    procedure :: init  => init_one_particle
    procedure :: reset => reset_one_particle
    procedure :: rep   => particle_reynolds_cal
    procedure :: force => particle_force_cal
    !
  end type partype

  integer :: numparticle,n_particle
  type(partype),allocatable,target :: partpack(:)
  real(mytype),allocatable,dimension(:) :: lxmin,lxmax,lymin,lymax,lzmin,lzmax
  !+------------------+-----------------------------------------------+
  !|      numparticle | number of particles in the domain             |
  !|       n_particle | total number of particles.                    |
  !|         partpack | pack of particles.                            |
  !|            lxmin |                                               |
  !|            lxmax |                                               |
  !|            lymin |                                               |
  !|            lymax |                                               |
  !|            lzmin |                                               |
  !|            lzmax | range of local domain.                        |
  !|                  |                                               |
  !+------------------+-----------------------------------------------+

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
    !
    pa%swap=.false.
    pa%new =.true.
    !
    pa%rankinn=nrank
    pa%rank2go=nrank
    !
    pa%x=0.0; pa%v=0.0
    !
    allocate(pa%dx(1:3,2:ntime),pa%dv(1:3,ntime))
    !
    pa%dx=0.0
    pa%dv=0.0
    !
    pa%dimeter=1.d-1
    pa%rho=10._mytype
    !
    pa%qe=1.d-5
    !
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
  !| This subroutine is used to calculate the Reynolds number of       |
  !| particles.                                                        |
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
  !| This subroutine is used to calculate the foce acting on a particle|                                        |
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

  end subroutine particle_init
  !+-------------------------------------------------------------------+
  ! The end of the subroutine particle_init                            |
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
    ! part_dmck_time=part_dmck_time+ptime()-timebeg
    ! print*,nrank,'|',numparticle
    ! do jpart=1,psize
    !   !
    !   pa=>particle(jpart)
    !   !
    !   if(pa%swap) then
    !     !
    !     write(*,'(3(A,1X,I0))')' ** particle',jpart,' moves from rank', &
    !                                  pa%rankinn,' to rank',pa%rank2go
    !     !
    !   endif
    !   !
    ! enddo
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

    ! if(numparticle>0) then
    !   !
    !   allocate(xyz(3,numparticle))
    !   !
    !   psize=msize(particles)

    !   j=0
    !   do p=1,psize
    !     !
    !     if(particles(p)%new) cycle
    !     !
    !     j=j+1

    !     xyz(1:3,j)  =particles(p)%x(1:3)

    !     k=k+3

    !   enddo

    !   rank2coll=nrank
    !   !
    ! else
    !   rank2coll=-1
    ! endif

    ! total_num_part=psum(numparticle)
    ! !
    ! write(num,'(I5.5)')part_file_numb

    ! file2write='./data/particle_coordinates_'//num//'.bin'

    ! call pwrite(trim(file2write),xyz)

  end subroutine dump_particle
  !+-------------------------------------------------------------------+
  ! The end of the subroutine partical_domain_check                    |
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
        partpack(j)%v(1)   =0.d0
        partpack(j)%v(2)   =0.d0
        partpack(j)%v(3)   =0.d0
        !
        partpack(j)%new=.false.
        !
      enddo

      deallocate(xp)      

    endif

  end subroutine particle_checkpoint
  !+-------------------------------------------------------------------+
  ! The end of the subroutine particle_checkpoint                      |
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

end module particle