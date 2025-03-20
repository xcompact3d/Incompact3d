!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

! the module defines the particle type
module particle_type
  
  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank

  implicit none

  ! type of particle
  type partype

    real(mytype) :: rho,mas,dimeter,re,vdiff,qe
    real(mytype) :: x(3),v(3),vf(3),f(3),b(3)
    real(mytype),allocatable,dimension(:,:) :: dx,dv
    integer :: id,rankinn,rank2go
    logical :: swap,inactive
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
    !|         inactive | a indicator to an inactive particle      |
    !+------------------+------------------------------------------+
    
    contains

    procedure :: init  => init_one_particle
    procedure :: reset => reset_one_particle
    procedure :: rep   => particle_reynolds_cal
    procedure :: force => particle_force_cal

  end type partype

  contains

  !+-------------------------------------------------------------------+
  !| This subroutine is used to init a particle.                       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-Jun-2022  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine init_one_particle(pa)

    use param, only : ntime

    class(partype) :: pa
    
    pa%swap=.false.
    pa%inactive =.true.

    pa%rankinn=nrank
    pa%rank2go=nrank

    pa%x=0.0; pa%v=0.0

    allocate(pa%dx(1:3,2:ntime),pa%dv(1:3,ntime))

    pa%dx=0.0
    pa%dv=0.0

    pa%dimeter=0.1_mytype
    pa%rho=10._mytype

    pa%qe=1.e-5_mytype

  end subroutine init_one_particle
  !+-------------------------------------------------------------------+
  !| The end of the subroutine initmesg.                               |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| This subroutine is to reset a particle to be new/inactivated .    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-Jun-2022  | Created by J. Fang @ Imperial College              |
  !+-------------------------------------------------------------------+
  subroutine reset_one_particle(pa)
    
    class(partype),target :: pa
    
    pa%swap=.false.
    pa%inactive =.true.
    
    pa%rankinn=nrank
    pa%rank2go=nrank
    
    pa%x=0.0; pa%v=0.0
    
    pa%dx=0.0
    pa%dv=0.0
    
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
    
    use param, only : re
    
    class(partype),target :: pa
    
    pa%vdiff = sqrt( (pa%vf(1)-pa%v(1))**2 + &
                     (pa%vf(2)-pa%v(2))**2 + &
                     (pa%vf(3)-pa%v(3))**2 )
    
    pa%re = pa%dimeter*pa%vdiff*re
    
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
    
    use param, only : re
    
    class(partype),target :: pa
    
    real(mytype) :: varc
     
    varc=18._mytype/(pa%rho*pa%dimeter**2*re)*(1._mytype+0.15_mytype*pa%re**0.687_mytype)
    
    pa%f(:) = varc*(pa%vf(:)-pa%v(:))
    
  end subroutine particle_force_cal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine particle_reynolds_cal.                  |
  !+-------------------------------------------------------------------+

end module particle_type

! this module contains some utility functions used for particle tracking 
module particle_utilities
  
  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nproc
  use particle_type
  
  implicit none

  interface mclean
    module procedure mclean_mytype
    module procedure mclean_particles
  end interface mclean 

  interface msize
    module procedure size_particle
    module procedure size_integer
  end interface msize

  interface pa2a
    module procedure pa2a_particle
  end interface pa2a

  interface mextend
     module procedure extend_particles
  end interface mextend

  contains

  !+-------------------------------------------------------------------+
  !| this function is to retune the size of a array                    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16-Jun-2022  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  ! The implementation of the functions is based on the PAMR code
  ! https://gitlab.com/fangjian/pamr/-/blob/master/src/memange.F90
  ! The original code is distributed under the APACHE-2.0 license reproduced below
  !
  ! Copyright 2018-2024 Jian Fang, STFC Daresbury Laboratory
  ! All Rights Reserved.
  !
  ! Licensed under the Apache License, Version 2.0 (the "License");
  ! you may not use this file except in compliance with the License.
  ! You may obtain a copy of the License at
  !
  ! https://gitlab.com/fangjian/pamr/-/blob/master/LICENSE
  !
  ! Unless required by applicable law or agreed to in writing, software
  ! distributed under the License is distributed on an
  ! "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
  ! either express or implied. See the License for the specific
  ! language governing permissions and limitations under the License.
  !+-------------------------------------------------------------------+
  pure function size_particle(var) result(nsize)
    
    type(partype),allocatable,intent(in) :: var(:)
    integer :: nsize
    
    if(allocated(var)) then
      nsize=size(var)
    else
      nsize=0
    endif
    
    return
    
  end function size_particle
  !
  pure function size_integer(var) result(nsize)
    
    integer,allocatable,intent(in) :: var(:)
    integer :: nsize
    
    if(allocated(var)) then
      nsize=size(var)
    else
      nsize=0
    endif
    
    return
    
  end function size_integer
  !+-------------------------------------------------------------------+
  !| The end of the subroutine size_particle.                          |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| this subroutine clean superfluous elements in a array             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| Created by J. Fang STFC Daresbury Laboratory                      |
  !+-------------------------------------------------------------------+
  ! The implementation of the functions is based on the PAMR code
  ! https://gitlab.com/fangjian/pamr/-/blob/master/src/memange.F90
  ! The original code is distributed under the APACHE-2.0 license reproduced below
  !
  ! Copyright 2018-2024 Jian Fang, STFC Daresbury Laboratory
  ! All Rights Reserved.
  !
  ! Licensed under the Apache License, Version 2.0 (the "License");
  ! you may not use this file except in compliance with the License.
  ! You may obtain a copy of the License at
  !
  ! https://gitlab.com/fangjian/pamr/-/blob/master/LICENSE
  !
  ! Unless required by applicable law or agreed to in writing, software
  ! distributed under the License is distributed on an
  ! "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
  ! either express or implied. See the License for the specific
  ! language governing permissions and limitations under the License.
  !+-------------------------------------------------------------------+
  subroutine mclean_mytype(var,n)
    
    ! arguments
    real(mytype),allocatable,intent(inout) :: var(:)
    integer,intent(in) :: n
    
    ! local data
    real(mytype),allocatable :: buffer(:)
    integer :: m
    logical :: lefc
    
    if(.not.allocated(var)) return
    
    if(n<=0) then
      deallocate(var)
      return
    endif
    
    ! clean
    allocate(buffer(n))
    
    buffer(1:n)=var(1:n)
    
    deallocate(var)
    
    call move_alloc(buffer,var)
    
  end subroutine mclean_mytype

  subroutine mclean_particles(var,n)
    
    ! arguments
    type(partype),allocatable,intent(inout) :: var(:)
    integer,intent(in) :: n
    
    ! local data
    type(partype),allocatable :: buffer(:)
    integer :: m
    logical :: lefc
    
    if(.not.allocated(var)) return
    
    if(n<=0) then
      deallocate(var)
      return
    endif
    
    ! clean
    allocate(buffer(n))
    
    buffer(1:n)=var(1:n)
    
    deallocate(var)
    
    call move_alloc(buffer,var)
    
  end subroutine mclean_particles
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
  ! The implementation of the functions is based on the PAMR code
  ! https://gitlab.com/fangjian/pamr/-/blob/master/src/memange.F90
  ! The original code is distributed under the APACHE-2.0 license reproduced below
  !
  ! Copyright 2018-2024 Jian Fang, STFC Daresbury Laboratory
  ! All Rights Reserved.
  !
  ! Licensed under the Apache License, Version 2.0 (the "License");
  ! you may not use this file except in compliance with the License.
  ! You may obtain a copy of the License at
  !
  ! https://gitlab.com/fangjian/pamr/-/blob/master/LICENSE
  !
  ! Unless required by applicable law or agreed to in writing, software
  ! distributed under the License is distributed on an
  ! "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
  ! either express or implied. See the License for the specific
  ! language governing permissions and limitations under the License.
  !+-------------------------------------------------------------------+
  subroutine extend_particles(var,n)
    
    ! arguments
    type(partype),allocatable,intent(inout) :: var(:)
    integer,intent(in) :: n
    
    ! local data
    type(partype),allocatable :: buffer(:)
    integer :: m,jpart
    
    if(.not. allocated(var)) then
      allocate(var(n))
      m=0
    else
    
      m=size(var)

      if(m==n) return
    
      call move_alloc(var, buffer)
    
      allocate(var(n))
      var(1:m)=buffer(1:m)
    
    endif
    
    ! initilise newly allocated particles
    do jpart=m+1,n
      call var(jpart)%init()
    enddo
    
    return
    
  end subroutine extend_particles
  !+-------------------------------------------------------------------+
  !| The end of the subroutine extend_particles.                       |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| this subroutine is to copy particle array.                        |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 18-Oct-2024  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  function particle_copy(vin) result(vout)

    ! arguments
    type(partype),allocatable,intent(in) :: vin(:)
    type(partype),allocatable :: vout(:)
    
    ! local data
    integer :: psize,big_number

    psize=msize(vin)

    if(psize>0) then
        allocate(vout(psize))
        vout=vin
    else
        allocate(vout(1))
        deallocate(vout)
    endif

    return

  end function particle_copy
  !+-------------------------------------------------------------------+
  !| The end of the subroutine particle_copy.                          |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| this subroutine is to swap particles via alltoall.                |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-Jun-2022  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine pa2a_particle(datasend,datarecv,sendtabl,recvtabl)
    
    use mpi
    use decomp_2d_constants, only : real_type
    
    ! arguments
    type(partype),allocatable,intent(in) ::  datasend(:)
    type(partype),allocatable,intent(out) :: datarecv(:)
    integer,intent(in) :: sendtabl(0:),recvtabl(0:)
    
    ! local data
    integer :: ierr,recvsize,jrank,jpart,jc,nindsize
    integer,allocatable :: senddispls(:),recvdispls(:)
    real(mytype),allocatable :: r8send(:,:),r8resv(:,:)
    
    integer,save :: newtype
    
    logical,save :: firstcal=.true.
    
    if(firstcal) then
      call mpi_type_contiguous(6,real_type,newtype,ierr)
      call mpi_type_commit(newtype,ierr)
      firstcal=.false.
    endif
    
    allocate(senddispls(0:nproc-1),recvdispls(0:nproc-1))
    
    senddispls=0
    recvdispls=0
    do jrank=1,nproc-1
      senddispls(jrank)=senddispls(jrank-1)+sendtabl(jrank-1)
      recvdispls(jrank)=recvdispls(jrank-1)+recvtabl(jrank-1)
    enddo
    recvsize=recvdispls(nproc-1)+recvtabl(nproc-1)
    
    nindsize=msize(datasend)
    
    allocate(r8send(6,nindsize))
    allocate(r8resv(6,recvsize))
    
    r8resv=0._mytype
    
    do jpart=1,nindsize
      r8send(1,jpart)=datasend(jpart)%x(1)
      r8send(2,jpart)=datasend(jpart)%x(2)
      r8send(3,jpart)=datasend(jpart)%x(3)
      r8send(4,jpart)=datasend(jpart)%v(1)
      r8send(5,jpart)=datasend(jpart)%v(2)
      r8send(6,jpart)=datasend(jpart)%v(3)
    enddo
    
    call mpi_alltoallv(r8send, sendtabl, senddispls, newtype, &
                       r8resv, recvtabl, recvdispls, newtype, &
                       mpi_comm_world, ierr)
    
    allocate(datarecv(recvsize))
    
    jc=0
    do jrank=0,nproc-1
      do jpart=1,recvtabl(jrank)
    
        jc=jc+1
    
        call datarecv(jc)%init()
    
        datarecv(jc)%x(1)=r8resv(1,jc)
        datarecv(jc)%x(2)=r8resv(2,jc)
        datarecv(jc)%x(3)=r8resv(3,jc)
        datarecv(jc)%v(1)=r8resv(4,jc)
        datarecv(jc)%v(2)=r8resv(5,jc)
        datarecv(jc)%v(3)=r8resv(6,jc)
    
      enddo
    enddo
    
  end subroutine pa2a_particle
  !+-------------------------------------------------------------------+
  !| The end of the subroutine pa2a_particle.                          |
  !+-------------------------------------------------------------------+
  
  !+-------------------------------------------------------------------+
  !| this subroutine is to swap data in the y-z direction.             |
  !+-------------------------------------------------------------------+
  subroutine pswap_yz(varin,varout)
    
    use decomp_2d, only : xsize,update_halo
    use param,     only : nclx,ncly,nclz
    
    ! arguments
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: varin
    real(mytype),intent(out) :: varout(1:xsize(1)+1,0:xsize(2)+1,0:xsize(3)+1)
    
    ! local data
    real(mytype),allocatable,dimension(:,:,:) :: var_halo

    call update_halo(varin,var_halo,1,opt_global=.true.,opt_pencil=1)

    
    varout(1:xsize(1),0:xsize(2)+1,0:xsize(3)+1)=var_halo(:,:,:)

    if(nclx) then
      varout(xsize(1)+1,:,:)=var_halo(1,:,:)
    endif

    return
    
  end subroutine pswap_yz
  !+-------------------------------------------------------------------+
  !| The end of the subroutine pswap_yz.                               |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| This subroutine is to swap particle infomation between ranks      |
  !+-------------------------------------------------------------------+
  !| tecplot format for now                                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16-06-2022  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine particle_swap(particles,num_active_particles)

    use mptool, only : ptabupd
    
    ! arguments
    type(partype),intent(inout),allocatable,target :: particles(:)
    integer,intent(inout),optional :: num_active_particles

    ! local data 
    integer :: p,psize,jrank,jpart,npart,n,newsize,npexit
    type(partype),pointer :: pa
    integer :: nsend(0:nproc-1),nrecv(0:nproc-1),nsend_total
    !+------------+--------------------------------+
    !| nsendtable | the table to record how much   |
    !|            | particle to send to a rank     |
    !+------------+--------------------------------+
    integer :: pr(0:nproc-1,1:msize(particles))
    type(partype),allocatable :: pa2send(:),pa2recv(:)
    real(mytype) :: timebeg,tvar1,tvar11,tvar2,tvar3,tvar4
    
    psize=msize(particles)
    
    nsend=0
    
    n=0
    pr=0
    npart=0

    if(present(num_active_particles)) then
        npexit=num_active_particles
    else
        npexit=psize
    endif
    
    do jpart=1,psize
    
      pa=>particles(jpart)
    
      if(pa%inactive) cycle
    
      ! to find out how many particle to send to which ranks
      if(pa%swap) then
    
        n=n+1
    
        nsend(pa%rank2go)=nsend(pa%rank2go)+1
    
        pr(pa%rank2go,nsend(pa%rank2go))=jpart
    
      endif
    
      npart=npart+1
    
      if(npart==npexit) exit

    enddo
    
    nsend_total=n
    
    ! synchronize recv table according to send table
    nrecv=ptabupd(nsend)
    
    ! to establish the buffer of storing particels about to send
    if(nsend_total>0) then
    
      allocate(pa2send(1:nsend_total))
    
      n=0
      do jrank=0,nproc-1
    
        do jpart=1,nsend(jrank)
    
          n=n+1
    
          p=pr(jrank,jpart)
    
          pa2send(n)=particles(p)
    
          call particles(p)%reset()

          if(present(num_active_particles)) num_active_particles=num_active_particles-1

        enddo
    
      enddo
    
    endif 
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! swap particle among ranks
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call pa2a(pa2send,pa2recv,nsend,nrecv)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end of swap particle among ranks
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    if(present(num_active_particles)) then
        call particle_add(particles,pa2recv,num_active_particles,n)
        num_active_particles=num_active_particles+n
    else
        call particle_add(particles,pa2recv)
    endif
    
  end subroutine particle_swap
  !+-------------------------------------------------------------------+
  ! The end of the subroutine particle_swap                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to remove all new particles from an array.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 22-10-2024  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine particle_array_clear(particles)
    
    use param,     only : nclx,ncly,nclz,xlx,yly,zlz
    
    ! arguments
    type(partype),intent(inout),allocatable,target :: particles(:)

    ! local data
    type(partype),allocatable :: buffer(:)
    integer :: n,j

    allocate(buffer(1:msize(particles)))

    n=0
    do j=1,msize(particles)

        if(particles(j)%inactive) cycle

        n=n+1

        buffer(n)=particles(j)

    enddo

    call mclean(buffer,n)

    deallocate(particles)

    call move_alloc(buffer,particles)

    return

  end subroutine particle_array_clear
  !+-------------------------------------------------------------------+
  ! The end of the subroutine particle_array_clear                     |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to add particles to the current particle arrary|
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-06-2022  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine particle_add(particle_cur,particle_new,num_active_particles,num_part_incr)
    
    ! arguments
    type(partype),intent(inout),allocatable,target :: particle_cur(:)
    type(partype),intent(in),allocatable :: particle_new(:)
    integer,intent(in),optional :: num_active_particles
    integer,intent(out),optional :: num_part_incr
    
    ! local data
    integer :: psize,newsize,n,jpart
    type(partype),pointer :: pa
    
    if(msize(particle_new)<=0) then
        if(present(num_part_incr)) num_part_incr=0
        return
    endif

    psize=msize(particle_cur)
    
    if(present(num_active_particles)) then

      ! now add the received particle in to the array, dynamically
      if(num_active_particles+msize(particle_new)>psize) then
    
        ! expand the particle array
        newsize=max(num_active_particles+msize(particle_new),num_active_particles+100)
      else
        newsize=psize
      endif

    else

        newsize=psize+msize(particle_new)

    endif

    call mextend(particle_cur,newsize)
    
    psize=newsize
    
    n=0
    do jpart=1,psize
    
      pa=>particle_cur(jpart)
    
      ! the particle is free for re-assigning
      if(pa%inactive) then
    
        if(n>=msize(particle_new)) exit
    
        n=n+1
    
        pa=particle_new(n)
        pa%inactive=.false.
    
      endif
    
    enddo
    
    if(present(num_part_incr)) num_part_incr=n
    
  end subroutine particle_add
  !+-------------------------------------------------------------------+
  ! The end of the subroutine particle_add                             |
  !+-------------------------------------------------------------------+
  
end module particle_utilities

module particle

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank,nproc,decomp_2d_abort
  use mptool
  use particle_type
  use particle_utilities
  use constants, only : pi

  implicit none

  integer :: n_local_particles,n_particles
  character(len=32) :: initype_particle
  character(len=16) :: bc_particle(6)
  real(mytype) :: particletime,sub_time_step,particle_inject_period,next_particle_inject_time
  type(partype),allocatable,target :: partpack(:),particles2inject(:)
  real(mytype),allocatable,dimension(:) :: lxmin,lxmax,lymin,lymax,lzmin,lzmax
  !+------------------+-----------------------------------------------+
  !|n_local_particles | number of particles in the domain             |
  !|      n_particles | total number of particles.                    |
  !| initype_particle | initilisation option for particle.            |
  !|      bc_particle | boundary condtion for particles.              |
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

  public :: local_domain_size,particle_init,intt_particles,           &
            particle_report,visu_particle,particle_checkpoint,        &
            n_particles,initype_particle,bc_particle,                 &
            particle_inject_period

  contains

  !+-------------------------------------------------------------------+
  !| This subroutine is to report time cost for particles.             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 18-06-2022  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine particle_report(reptyp)
    
    character(len=*),intent(in) :: reptyp

    integer :: ibc
    character(len=4) :: bcname(6)
    logical :: file_exists
    integer :: file_size

    if(nrank .ne. 0) return

    bcname(1) = 'xmin'
    bcname(2) = 'xmax'
    bcname(3) = 'ymin'
    bcname(4) = 'ymax'
    bcname(5) = 'zmin'
    bcname(6) = 'zmax'

    if(reptyp=='input') then

        do ibc=1,6
            if(trim(bc_particle(ibc)) == 'periodic') then
                print*,'** particles b.c. at ',bcname(ibc),': periodic'
            elseif(trim(bc_particle(ibc)) == 'reflective') then
                print*,'** particles b.c. at ',bcname(ibc),': reflective'
            elseif(trim(bc_particle(ibc)) == 'outflow') then
                print*,'** particles b.c. at ',bcname(ibc),': outflow'
            else
                print*,'!! bc_particle:',bc_particle
                call decomp_2d_abort(1, "bc_particle not defined")
            endif
        enddo

        if(trim(initype_particle)=='random') then
            write(*,'(A,I0,A)')' ** particle initilisation: generate ', &
                        n_particles,' random particles '
        elseif(trim(initype_particle)=='uniform') then
            write(*,'(A,I0,A)')' ** particle initilisation: generate ', &
                        n_particles,' uniformly distributed particles '
        else

            inquire(file=trim(initype_particle), exist=file_exists, size=file_size)
            if(file_exists) then
                write(*,*)'** read particles from file: ',trim(initype_particle)
            else
                call decomp_2d_abort(1, "initype_particle error @ particle_report:"//trim(initype_particle))
            endif
            
        endif

        if(particle_inject_period>0.0_mytype) then
            print*,'** particles are reguarly injected every t=',particle_inject_period
        endif

    endif

    
  end subroutine particle_report
  !+-------------------------------------------------------------------+
  !| The end of the subroutine particle_report.                        |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| This subroutine is to initilise particle module.                  |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-11-2021  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine particle_init(pxmin,pxmax,pymin,pymax,pzmin,pzmax)
    
    use param,     only : xlx,yly,zlz
    use var,       only : itime

    real(mytype),intent(in),optional :: pxmin,pxmax,pymin,pymax,pzmin,pzmax
    ! range of initial particles
    
    ! local data
    real(mytype) :: particle_range(6)
    logical :: file_exists
    integer :: psize

    particle_range(1) = 0.0_mytype
    particle_range(2) = xlx
    particle_range(3) = 0.0_mytype
    particle_range(4) = yly
    particle_range(5) = 0.0_mytype
    particle_range(6) = zlz

    if(present(pxmin)) particle_range(1) = pxmin
    if(present(pxmax)) particle_range(2) = pxmax
    if(present(pymin)) particle_range(3) = pymin
    if(present(pymax)) particle_range(4) = pymax
    if(present(pzmin)) particle_range(5) = pzmin
    if(present(pzmax)) particle_range(6) = pzmax

    ! generate random particles within the local domain
    if(trim(initype_particle)=='random') then 
        call particle_gen_random(partpack,psize,particle_range)
    elseif(trim(initype_particle)=='uniform') then 
        call particle_gen_uniform(partpack,psize,particle_range)
    else
        inquire(file=trim(initype_particle), exist=file_exists)
        if(file_exists) then
            call particle_read_bin(partpack,trim(initype_particle),psize)
        else
            call decomp_2d_abort(1,"initype_particle error @ particle_init")
        endif
    endif

    n_local_particles=psize

    n_particles=psum(n_local_particles)

    if(nrank==0) print*,'** ',n_particles,'particles are generated'

    if(particle_inject_period>0.0_mytype) then
        if(n_local_particles>0) particles2inject=partpack
        call particle_checkpoint(mode='write',filename='particle_inject.bin')

        next_particle_inject_time = 0.0_mytype + particle_inject_period
    endif

  end subroutine particle_init
  !+-------------------------------------------------------------------+
  ! The end of the subroutine particle_init                            |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| This subroutine is to add more particles to the domain.           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-11-2021  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine inject_particles(particles)
    
    ! local data
    type(partype),intent(in),allocatable :: particles(:)
    integer :: num_new_particle,n
    
    call particle_add(partpack,particles,n_local_particles,n)
    
    n_local_particles=n_local_particles+n

  end subroutine inject_particles
  !+-------------------------------------------------------------------+
  ! The end of the subroutine inject_particles                         |
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
  subroutine intt_particles(ux1,uy1,uz1,time1)
    
    use variables 
    use param
    use decomp_2d, only : xsize
    use, intrinsic :: ieee_arithmetic
    
    ! arguments
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1
    real(mytype),intent(in) :: time1
    
    ! local data 
    integer :: p,psize,jpart,npart,iptr
    integer,save :: old_num_part=0
    type(partype),pointer :: pa
    real(mytype),save :: time0,tsro
    !
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: uxp,uyp,uzp
    real(mytype),allocatable,dimension(:,:,:),save :: ux0,uy0,uz0
    real(mytype),allocatable :: xcor(:,:),dxco(:,:,:),vcor(:,:),dvco(:,:,:)
    
    logical,save :: firstcal=.true.,particle_inject=.false.
    
    if(firstcal) then

      allocate( ux0(xsize(1),xsize(2),xsize(3)), &
                uy0(xsize(1),xsize(2),xsize(3)), &
                uz0(xsize(1),xsize(2),xsize(3)) )

      ux0=ux1
      uy0=uy1
      uz0=uz1
      
      particletime=t

      sub_time_step=dt

      time0=time1

      tsro=sub_time_step/dt

      if(particle_inject_period>0.0_mytype) then
        particle_inject=.true.
      endif

      old_num_part=n_particles

      firstcal=.false.

    endif
    
    do while(particletime<time1)
    
      particletime=particletime+sub_time_step
      
      if(particle_inject .and. particletime>=next_particle_inject_time) then

        call inject_particles(particles2inject)

        next_particle_inject_time=next_particle_inject_time+particle_inject_period

      endif
    
      uxp=linintp(time0,time1,ux0,ux1,particletime)
      uyp=linintp(time0,time1,uy0,uy1,particletime)
      uzp=linintp(time0,time1,uz0,uz1,particletime)

      do iptr=1,iadvance_time

        call particle_force(uxp,uyp,uzp)

        psize=msize(partpack)
        allocate(xcor(3,psize),dxco(3,psize,ntime))
        allocate(vcor(3,psize),dvco(3,psize,ntime))

        npart=0

        do jpart=1,psize
  
          pa=>partpack(jpart)
  
          if(pa%inactive) cycle
  
          npart=npart+1
  
          if(npart>n_local_particles) exit
  
          xcor(:,npart)=pa%x(:)
  
          dxco(:,npart,1)=pa%v(:)
  
          dxco(:,npart,2:ntime)=pa%dx(:,2:ntime)
  
          vcor(:,npart)=pa%v(:)
  
          dvco(:,npart,1)=pa%f(:)
          dvco(:,npart,2:ntime)=pa%dv(:,2:ntime)
  
        enddo
    
        if (itimescheme.eq.1) then
           !>>> Euler
           vcor=gdt(iptr)*tsro*dvco(:,:,1)+vcor
           !
           xcor=gdt(iptr)*tsro*dxco(:,:,1)+xcor
           !
        elseif(itimescheme.eq.2) then
           !>>> Adam-Bashforth second order (AB2)
           !
           if(itime.eq.1 .and. irestart.eq.0) then
             ! Do first time step with Euler
             vcor=gdt(iptr)*tsro*dvco(:,:,1)+vcor
             !
             xcor=gdt(iptr)*tsro*dxco(:,:,1)+xcor
           else
             vcor=adt(iptr)*dvco(:,:,1)+bdt(iptr)*dvco(:,:,2)+vcor
             !
             xcor=adt(iptr)*dxco(:,:,1)+bdt(iptr)*dxco(:,:,2)+xcor
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
    
              xcor=dt*tsro*dxco(:,:,1)+xcor
           elseif(itime.eq.2.and.irestart.eq.0) then
              ! Do second time step with AB2
              vcor=onepfive*dt*tsro*dvco(:,:,1)-half*dt*tsro*dvco(:,:,2)+vcor
              dvco(:,:,3)=dvco(:,:,2)
    
              xcor=onepfive*dt*tsro*dxco(:,:,1)-half*dt*tsro*dxco(:,:,2)+xcor
              dxco(:,:,3)=dxco(:,:,2)
           else
              ! Finally using AB3
              vcor=adt(iptr)*tsro*dvco(:,:,1)+bdt(iptr)*tsro*dvco(:,:,2)+cdt(iptr)*tsro*dvco(:,:,3)+vcor
              dvco(:,:,3)=dvco(:,:,2)
    
              xcor=adt(iptr)*tsro*dxco(:,:,1)+bdt(iptr)*tsro*dxco(:,:,2)+cdt(iptr)*tsro*dxco(:,:,3)+xcor
              dxco(:,:,3)=dxco(:,:,2)
           endif
           dvco(:,:,2)=dvco(:,:,1)
           dxco(:,:,2)=dxco(:,:,1)
           !
        elseif(itimescheme.eq.5) then
           !>>> Runge-Kutta (low storage) RK3
           if(iptr.eq.1) then
              vcor=gdt(iptr)*tsro*dvco(:,:,1)+vcor
              xcor=gdt(iptr)*tsro*dxco(:,:,1)+xcor
           else
              vcor=adt(iptr)*tsro*dvco(:,:,1)+bdt(iptr)*tsro*dvco(:,:,2)+vcor
              xcor=adt(iptr)*tsro*dxco(:,:,1)+bdt(iptr)*tsro*dxco(:,:,2)+xcor
           endif
           dvco(:,:,2)=dvco(:,:,1)
           dxco(:,:,2)=dxco(:,:,1)
  
        endif

        ! put back from local array to particle array
        npart=0
        do jpart=1,psize
    
          pa=>partpack(jpart)
    
          if(pa%inactive) cycle
    
          npart=npart+1
  
          if(npart>n_local_particles) exit
    
          if(ieee_is_nan(vcor(1,npart)) .or. ieee_is_nan(vcor(2,npart)) .or. ieee_is_nan(vcor(3,npart))) then
            print*, npart,pa%x(:),pa%v(:)
            call decomp_2d_abort(1,"particle velocity NaN @ intt_particles")
          endif
    
          pa%v(:)=vcor(:,npart)
          pa%x(:)=xcor(:,npart)
          
          pa%dv(:,2:ntime)=dvco(:,npart,2:ntime)
          pa%dx(:,2:ntime)=dxco(:,npart,2:ntime)

        enddo

        deallocate(xcor,dxco,vcor,dvco)
        
        call particle_domain_check(partpack,n_local_particles)
        
        call particle_swap(partpack,n_local_particles)
      
      enddo

      n_particles=psum(n_local_particles)
  
      if(nrank==0) then
  
        if(n_particles.ne.old_num_part) then
          write(*,'(3(A,I0))')' ** number of particles changes from ',old_num_part, &
                                      ' -> ',n_particles,' : ',n_particles-old_num_part
          old_num_part=n_particles
        endif
        
      endif

    enddo
    
    ux0=ux1
    uy0=uy1
    uz0=uz1
    
    time0=time1
    
  end subroutine intt_particles
  !+-------------------------------------------------------------------+
  ! The end of the subroutine intt_particles                           |
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
    use utilities, only : gen_snapshotname,int_to_str

    integer, intent(in) :: itime

    integer :: num,j,iovp,psize,p
    character(len=64) :: file2write,fname

    real(mytype),allocatable :: xyz(:,:)

    ! Snapshot number
    num = itime/ioutput 

    psize=msize(partpack)
    
    allocate(xyz(3,n_local_particles))
    
    j=0
    do p=1,psize
    
      if(partpack(p)%inactive) cycle
    
      j=j+1

      xyz(1:3,j)  =partpack(p)%x(1:3)

      if(j==n_local_particles) exit

    enddo

    file2write='particle_coordinates_'//int_to_str(num)//'.bin'

    call pwrite('./data/'//trim(file2write),xyz)

    if(nrank==0) then

      fname = gen_snapshotname('./data', 'snapshot_particle', num, "xdmf")

      OPEN(newunit=iovp,file=trim(fname))

      write(iovp,'(A)')'<?xml version="1.0" ?>'
      write(iovp,'(A)')'<Xdmf Version="3.0" xmlns:xi="http://www.w3.org/2001/XInclude">'
      write(iovp,'(A)')'  <Domain>'
      write(iovp,'(A)')'    <Grid GridType="Collection" CollectionType="Temporal">'

      write(iovp,'(A)')'      <Grid Name="ParticleData_t0" GridType="Uniform">'
      write(iovp,'(A,I0,A)')'        <Topology TopologyType="Polyvertex" NodesPerElement="',n_particles,'"/>'
      write(iovp,'(A)')'        <Geometry GeometryType="XYZ">'
      write(iovp,'(2(A,I0),A)')'          <DataItem Format="binary" Dimensions="',n_particles,' 3 " NumberType="Float" Precision="',mytype,'">'
      write(iovp,'(A,A,A)')'            ',trim(file2write),' </DataItem>'
      write(iovp,'(A)')'        </Geometry>'
      write(iovp,'(A)')'      </Grid>'

      write(iovp,'(A)')'    </Grid>'
      write(iovp,'(A)')'  </Domain>'
      write(iovp,'(A)')'</Xdmf>'

      close(iovp)

      print*,'<< ',trim(fname)

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
  subroutine particle_checkpoint(mode,filename)

    character(len=*),intent(in) :: mode
    character(len=*),intent(in),optional :: filename

    real(mytype),allocatable,dimension(:,:) :: xp

    integer :: j,p,file_size,data_size,ninject,psize
    logical :: file_exists
    real(mytype) :: xpmin(3),xpmax(3)
    character(len=80) :: particle_res_file

    NAMELIST /ParTrack/ n_particles,particle_res_file,next_particle_inject_time

    if(mode=='write') then

      psize=msize(partpack)

      n_particles=psum(n_local_particles)

      allocate(xp(1:3,n_local_particles))
  
      j=0
      do p=1,psize
    
        if(partpack(p)%inactive) cycle
    
        j=j+1

        xp(1:3,j)  =partpack(p)%x(1:3)

        if(j==n_local_particles) exit

      enddo

      inquire(file=filename, exist=file_exists,size=file_size)

      if(file_exists) then
        call rename(filename,filename//'-bak')
      endif
  
      call pwrite(filename,xp)
  
      deallocate(xp)

      if(nrank == 0) then
        open (111,file='restart.info',action='write',access='append')
        write(111,'(A)')'&ParTrack'
        write(111,'(A)')'!========================='
        write(111,'(A,I13)') 'n_particles=  ',n_particles
        write(111,'(A)')     'particle_res_file=  "'//filename//'"'
        write(111,'(A,E20.13E2)') 'next_particle_inject_time=  ',next_particle_inject_time
        write(111,'(A)')'/End'
        write(111,'(A)')'!========================='
        close(111)
        print*,'<< restart.info'
      endif

    elseif(mode=='read') then

      open (111,file='restart.info',action='read')
      read(111, nml=ParTrack)
      close(111)
      if(nrank==0) print*,'>> restart.info'

      inquire(file=trim(particle_res_file), exist=file_exists,size=file_size)

      if(file_exists) then

        data_size=file_size/mytype/3

        if(data_size .ne. n_particles) then
            if(nrank==0) then
                print*,' !! warning !!, the datasize of file:',trim(particle_res_file),data_size, &
                          ' not consistent with n_particles',n_particles
            endif
        endif

        psize=numdist(n_particles)


        allocate(xp(1:3,psize))

        call pread(trim(particle_res_file),xp)

        allocate(partpack(1:psize))

        xpmin=1.e10_mytype
        xpmax=0.0_mytype

        do j=1,psize
    
          call partpack(j)%init()
    
          partpack(j)%x(1:3)=xp(1:3,j)
    
          partpack(j)%v(1) = 0._mytype
          partpack(j)%v(2) = 0._mytype
          partpack(j)%v(3) = 0._mytype
    
          partpack(j)%inactive  = .false.
    
          xpmin(1)=min(xpmin(1),partpack(j)%x(1))
          xpmin(2)=min(xpmin(2),partpack(j)%x(2))
          xpmin(3)=min(xpmin(3),partpack(j)%x(3))
          xpmax(1)=max(xpmax(1),partpack(j)%x(1))
          xpmax(2)=max(xpmax(2),partpack(j)%x(2))
          xpmax(3)=max(xpmax(3),partpack(j)%x(3))
        enddo

        deallocate(xp)

        n_local_particles=psize

        call particle_domain_check(partpack,n_local_particles)

        call particle_swap(partpack,n_local_particles)

        xpmin(1)=pmin(xpmin(1))
        xpmin(2)=pmin(xpmin(2))
        xpmin(3)=pmin(xpmin(3))
        xpmax(1)=pmax(xpmax(1))
        xpmax(2)=pmax(xpmax(2))
        xpmax(3)=pmax(xpmax(3))

        if(nrank==0) then
            print*,' ** min particle coordinates: ',xpmin
            print*,' ** max particle coordinates: ',xpmax
        endif

        if(particle_inject_period>0.0_mytype) then

            ! for case with injected particles, read the particle_inject.bin
            call particle_read_bin(particles=particles2inject,filename='particle_inject.bin',nump_read=ninject)

        endif

       else
         ! for the checkpoint file doesnot exist, re-intilise particles
         call decomp_2d_abort(1,"check point missing for particles, looking for: "//trim(particle_res_file))

       endif

    endif

  end subroutine particle_checkpoint
  !+-------------------------------------------------------------------+
  ! The end of the subroutine particle_checkpoint                      |
  !+-------------------------------------------------------------------+

  subroutine particle_read_bin(particles,filename,nump_read)

    ! arguments
    type(partype),intent(out),allocatable :: particles(:)
    character(len=*),intent(in) :: filename
    integer,intent(out) :: nump_read

    ! local data
    real(mytype),allocatable,dimension(:,:) :: xp
    integer :: j,file_size,data_size

    inquire(file=filename, size=file_size)

    data_size=file_size/mytype/3

    nump_read=numdist(data_size)

    allocate(xp(1:3,nump_read))

    call pread(filename,xp)

    allocate(particles(1:nump_read))

    do j=1,nump_read
    
      call particles(j)%init()
    
      particles(j)%x(1:3)=xp(1:3,j)
    
      particles(j)%v(1) = 0._mytype
      particles(j)%v(2) = 0._mytype
      particles(j)%v(3) = 0._mytype
    
      particles(j)%inactive  = .false.
    
    enddo

    deallocate(xp)

    call particle_domain_check(particles)

    call particle_swap(particles)

    call particle_array_clear(particles)

  end subroutine particle_read_bin

  !+-------------------------------------------------------------------+
  !| This subroutine is to get the fluids velocity on particles.       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 15-11-2021  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine field_var(ux1,uy1,uz1,bx,by,bz)
    
    use MPI
    use param,     only : dx,dy,dz,istret,nclx,ncly,nclz,xlx,yly,zlz
    use variables, only : yp,ny,nz
    use decomp_2d, only : xsize,xstart,xend,update_halo
    use actuator_line_model_utils, only: trilinear_interpolation
    
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in),optional :: bx,by,bz
    
    ! local data
    integer :: jpart,npart,i,j,k,psize
    type(partype),pointer :: pa
    real(mytype) :: x1,y1,z1,x2,y2,z2
    real(mytype) :: test(xsize(1),xsize(2),xsize(3))
    real(mytype),allocatable,dimension(:,:,:) :: ux1_halo,uy1_halo,    &
                                                 uz1_halo
    real(mytype),allocatable,dimension(:,:,:) :: bx_halo,by_halo,bz_halo
    
    real(mytype),allocatable,save :: xx(:),yy(:),zz(:)
    real(mytype) :: vec000(3),vec100(3),vec001(3),vec101(3), &
                    vec010(3),vec110(3),vec011(3),vec111(3)
    logical,save :: firstcal=.true.
    
    if(firstcal) then
    
      allocate(xx(0:xsize(1)+1),yy(0:xsize(2)+1),zz(0:xsize(3)+1))
    
      do i=0,xsize(1)+1
        xx(i)=real(i-1,mytype)*dx
      enddo
      do j=0,xsize(2)+1
    
        if(j+xstart(2)-1>ny) then
          yy(j)=2._mytype*yp(ny)-yp(ny-1)
        elseif(j+xstart(2)-1<1) then
          yy(j)=2._mytype*yp(1)-yp(2)
        else
          yy(j)=yp(j+xstart(2)-1)
        endif
    
      enddo
    
      if(xsize(3)==1) then
        zz(:)=0._mytype
      else
        do k=0,xsize(3)+1
          zz(k)=real((k+xstart(3)-2),mytype)*dz
        enddo
      endif
    
      firstcal=.false.
    
    endif
    
    allocate( ux1_halo(1:xsize(1)+1,0:xsize(2)+1,0:xsize(3)+1),        &
              uy1_halo(1:xsize(1)+1,0:xsize(2)+1,0:xsize(3)+1),        &
              uz1_halo(1:xsize(1)+1,0:xsize(2)+1,0:xsize(3)+1) )
    
    call pswap_yz(ux1,ux1_halo)
    call pswap_yz(uy1,uy1_halo)
    call pswap_yz(uz1,uz1_halo)

    if(present(bx) .and. present(by) .and. present(bz)) then
    
    allocate( bx_halo(1:xsize(1)+1,0:xsize(2)+1,0:xsize(3)+1),        &
              by_halo(1:xsize(1)+1,0:xsize(2)+1,0:xsize(3)+1),        &
              bz_halo(1:xsize(1)+1,0:xsize(2)+1,0:xsize(3)+1) )
    
      call pswap_yz(bx,bx_halo)
      call pswap_yz(by,by_halo)
      call pswap_yz(bz,bz_halo)
    
    endif
    
    psize=msize(partpack)
    
    npart=0
    
    do jpart=1,psize
    
      pa=>partpack(jpart)
    
      if(pa%inactive) cycle
    
      loopk: do k=1,xsize(3)
    
        z1=zz(k)
        z2=zz(k+1)
    
        do j=1,xsize(2)
          y1=yy(j)
          y2=yy(j+1)
          do i=1,xsize(1)
            x1=xx(i)
            x2=xx(i+1)
    
            if( pa%x(1)>=x1 .and. pa%x(1)<=x2 .and. &
                pa%x(2)>=y1 .and. pa%x(2)<=y2 .and. &
                pa%x(3)>=z1 .and. pa%x(3)<=z2 ) then

              vec000(1)=ux1_halo(i,  j,  k);   vec000(2)=uy1_halo(i,  j,  k);   vec000(3)=uz1_halo(i,  j,  k)
              vec100(1)=ux1_halo(i+1,j,  k);   vec100(2)=uy1_halo(i+1,j,  k);   vec100(3)=uz1_halo(i+1,j,  k)
              vec010(1)=ux1_halo(i,  j+1,k);   vec010(2)=uy1_halo(i,  j+1,k);   vec010(3)=uz1_halo(i,  j+1,k)
              vec110(1)=ux1_halo(i+1,j+1,k);   vec110(2)=uy1_halo(i+1,j+1,k);   vec110(3)=uz1_halo(i+1,j+1,k)
              vec001(1)=ux1_halo(i,  j,  k+1); vec001(2)=uy1_halo(i,  j,  k+1); vec001(3)=uz1_halo(i,  j,  k+1)
              vec101(1)=ux1_halo(i+1,j,  k+1); vec101(2)=uy1_halo(i+1,j,  k+1); vec101(3)=uz1_halo(i+1,j,  k+1)
              vec011(1)=ux1_halo(i,  j+1,k+1); vec011(2)=uy1_halo(i,  j+1,k+1); vec011(3)=uz1_halo(i,  j+1,k+1)
              vec111(1)=ux1_halo(i+1,j+1,k+1); vec111(2)=uy1_halo(i+1,j+1,k+1); vec111(3)=uz1_halo(i+1,j+1,k+1)
    
              ! locate the particle, do the interpolation
              pa%vf=trilinear_interpolation( x1,y1,z1,x2,y2,z2,           &
                                             pa%x(1),pa%x(2),pa%x(3),     &
                                             vec000,vec100,vec001,vec101, &
                                             vec010,vec110,vec011,vec111)
  
    
              if(present(bx) .and. present(by) .and. present(bz)) then
                !
                vec000(1)=bx_halo(i,  j,  k);   vec000(2)=by_halo(i,  j,  k);   vec000(3)=bz_halo(i,  j,  k)
                vec100(1)=bx_halo(i+1,j,  k);   vec100(2)=by_halo(i+1,j,  k);   vec100(3)=bz_halo(i+1,j,  k)
                vec010(1)=bx_halo(i,  j+1,k);   vec010(2)=by_halo(i,  j+1,k);   vec010(3)=bz_halo(i,  j+1,k)
                vec110(1)=bx_halo(i+1,j+1,k);   vec110(2)=by_halo(i+1,j+1,k);   vec110(3)=bz_halo(i+1,j+1,k)
                vec001(1)=bx_halo(i,  j,  k+1); vec001(2)=by_halo(i,  j,  k+1); vec001(3)=bz_halo(i,  j,  k+1)
                vec101(1)=bx_halo(i+1,j,  k+1); vec101(2)=by_halo(i+1,j,  k+1); vec101(3)=bz_halo(i+1,j,  k+1)
                vec011(1)=bx_halo(i,  j+1,k+1); vec011(2)=by_halo(i,  j+1,k+1); vec011(3)=bz_halo(i,  j+1,k+1)
                vec111(1)=bx_halo(i+1,j+1,k+1); vec111(2)=by_halo(i+1,j+1,k+1); vec111(3)=bz_halo(i+1,j+1,k+1)

                pa%b=trilinear_interpolation( x1,y1,z1,x2,y2,z2,            &
                                              pa%x(1),pa%x(2),pa%x(3),      &
                                              vec000,vec100,vec001,vec101,  &
                                              vec010,vec110,vec011,vec111)
    
              endif
    
              exit loopk

            endif
    
          enddo
        enddo
      enddo loopk
    
      npart=npart+1
    
      if(npart==n_local_particles) exit

    enddo
    
  end subroutine field_var
  !+-------------------------------------------------------------------+
  ! The end of the subroutine field_var                                |
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
  subroutine particle_gen_random(particle_new,particle_size,particle_range)
    
    ! arguments
    type(partype),intent(out),allocatable :: particle_new(:)
    integer,intent(out) :: particle_size
    real(mytype),intent(in) :: particle_range(6)
    
    ! local data
    integer :: p,j,local_size,code,n,i
    integer,allocatable :: seed(:)
    real(mytype) :: dx,dy,dz,x,y,z,ran1,ran2,ran3,plx,ply,plz
    

    plx=particle_range(2)-particle_range(1)
    ply=particle_range(4)-particle_range(3)
    plz=particle_range(6)-particle_range(5)

    local_size=numdist(n_particles)
    
    allocate(particle_new(1:local_size))
    
    call system_clock(count=code)
    call random_seed(size = n)
    allocate(seed(n))
    seed=code+63946*(nrank+1)
    call random_seed(put = seed)

    p=0
    do j=1,local_size
    
      call random_number(ran1)
      call random_number(ran2)
      call random_number(ran3)
    
      x=plx*ran1+particle_range(1)
      y=ply*ran2+particle_range(3)
      z=plz*ran3+particle_range(5)
    
      if( x>=lxmin(nrank) .and. x<lxmax(nrank) .and. &
          y>=lymin(nrank) .and. y<lymax(nrank) .and. &
          z>=lzmin(nrank) .and. z<lzmax(nrank) ) then
    
        p=p+1
    
        call particle_new(p)%init()
    
        particle_new(p)%x(1)=x
        particle_new(p)%x(2)=y
        particle_new(p)%x(3)=z
    
        particle_new(p)%v(1)   =0._mytype
        particle_new(p)%v(2)   =0._mytype
        particle_new(p)%v(3)   =0._mytype
    
        particle_new(p)%inactive=.false.
    
      endif
    
    enddo
    
    call mclean(particle_new,p)
    
    particle_size=p

    deallocate(seed)
    
  end subroutine particle_gen_random

  subroutine particle_gen_uniform(particle_new,particle_size,particle_range)
    
    use param,     only : nclx,ncly,nclz,xlx,yly,zlz

    ! arguments
    type(partype),intent(out),allocatable :: particle_new(:)
    integer,intent(out) :: particle_size
    real(mytype),intent(in) :: particle_range(6)
    
    ! local data
    integer :: p,i,j,k,n,local_size,code,npx,npy,npz
    integer :: max_n_particles
    integer,allocatable :: seed(:)
    real(mytype) :: x,y,z,detalx,detaly,detalz
    real(mytype) :: plx,ply,plz
    

    plx=particle_range(2)-particle_range(1)
    ply=particle_range(4)-particle_range(3)
    plz=particle_range(6)-particle_range(5)

    call particle_dis_xyz(n_particles,plx,ply,plz,npx,npy,npz)
    
    if(nrank==0) then
        print*,'** actual number of particle in 3 direction',npx,npy,npz,'=',npx*npy*npz
    endif

    detalx = (particle_range(2)-particle_range(1))/real(npx,mytype)
    detaly = (particle_range(4)-particle_range(3))/real(npy,mytype)
    detalz = (particle_range(6)-particle_range(5))/real(npz,mytype)
    
    max_n_particles=npx*npy*npz

    allocate(particle_new(max_n_particles))
    
    p=0
    do k=1,npz
    do j=1,npy
    do i=1,npx
    
      x=particle_range(1)+0.5_mytype*detalx+detalx*real(i-1,mytype)
      y=particle_range(3)+0.5_mytype*detaly+detaly*real(j-1,mytype)
      z=particle_range(5)+0.5_mytype*detalz+detalz*real(k-1,mytype)
    
      if( x>=lxmin(nrank) .and. x<lxmax(nrank) .and. &
          y>=lymin(nrank) .and. y<lymax(nrank) .and. &
          z>=lzmin(nrank) .and. z<lzmax(nrank) ) then
    
        p=p+1
    
        call particle_new(p)%init()
    
        particle_new(p)%x(1) = x
        particle_new(p)%x(2) = y
        particle_new(p)%x(3) = z
    
        particle_new(p)%v(1) = 0._mytype
        particle_new(p)%v(2) = 0._mytype
        particle_new(p)%v(3) = 0._mytype
    
        particle_new(p)%inactive=.false.
    
      endif
    
    enddo
    enddo
    enddo
    
    call mclean(particle_new,p)
    
    particle_size=p

  end subroutine particle_gen_uniform
  !+-------------------------------------------------------------------+
  !| The end of the subroutine particle_gen.                           |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| This subroutine is to distribute particles evenly in 3 directions |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-10-2024  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine particle_dis_xyz(nptotal,plx,ply,plz,npx,npy,npz)

    ! arguments
    integer,intent(in) ::  nptotal
    real(mytype),intent(in) :: plx,ply,plz
    integer,intent(out) ::  npx,npy,npz

    ! local data
    real(mytype) :: n1


    if(abs(plx)<epsilon(1.0_mytype) .and. abs(ply)<epsilon(1.0_mytype) .and. abs(plz)<epsilon(1.0_mytype)) then

        npx=1; npy=1; npz=1

    elseif(abs(plx)<epsilon(1.0_mytype) .and. abs(ply)<epsilon(1.0_mytype)) then

        npx=1
        npy=1
        npz=nptotal

    elseif(abs(plx)<epsilon(1.0_mytype) .and. abs(plz)<epsilon(1.0_mytype)) then

        npx=1
        npy=nptotal
        npz=1

    elseif(abs(ply)<epsilon(1.0_mytype) .and. abs(plz)<epsilon(1.0_mytype)) then

        npx=nptotal
        npy=1
        npz=1

    elseif(abs(plx)<epsilon(1.0_mytype)) then

        n1=sqrt(real(nptotal,mytype)/(ply*plz))
        
        npx=1 
        npy=int(n1*ply)
        npz=int(n1*plz)

    elseif(abs(ply)<epsilon(1.0_mytype)) then

        n1=sqrt(real(nptotal,mytype)/(plx*plz))
        
        npx=int(n1*plx)
        npy=1
        npz=int(n1*plz)

    elseif(abs(plz)<epsilon(1.0_mytype)) then

        n1=sqrt(real(nptotal,mytype)/(plx*ply))
        
        npx=int(n1*plx)
        npy=int(n1*ply)
        npz=1

    else

        n1=(real(nptotal,mytype)/(plx*ply*plz))**(1.0/3.0)
        npx=int(n1*plx)
        npy=int(n1*ply)
        npz=int(n1*plz)

    endif

    npx = max(npx,1)
    npy = max(npy,1)
    npz = max(npz,1)

  end subroutine particle_dis_xyz
  !+-------------------------------------------------------------------+
  !| The end of the subroutine particle_dis_xyz.                       |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| This subroutine is to calcualte the size and range of local domain|
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 15-11-2021  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine local_domain_size
    
    use param,     only : dx,dy,dz,istret
    use variables, only : yp,ny
    use decomp_2d, only : xstart,xend
    use actuator_line_model_utils
    
    integer :: nyr,nzr,jrank
    
    allocate( lxmin(0:nproc-1),lxmax(0:nproc-1),  &
              lymin(0:nproc-1),lymax(0:nproc-1),  &
              lzmin(0:nproc-1),lzmax(0:nproc-1)   )
    
    lxmin=0._mytype; lxmax=0._mytype
    lymin=0._mytype; lymax=0._mytype
    lzmin=0._mytype; lzmax=0._mytype
    
    lxmin(nrank)=real(0,mytype)*dx
    lxmax(nrank)=real(xend(1),mytype)*dx
    
    if (istret==0) then
      lymin(nrank)=real(xstart(2)-1,mytype)*dy
      lymax(nrank)=real(xend(2),mytype)*dy
    else
      lymin(nrank)=yp(xstart(2))
      nyr=min(ny,xend(2)+1)
      lymax(nrank)=yp(nyr)
    endif
    
    lzmin(nrank)=real((xstart(3)-1),mytype)*dz
    nzr=xend(3)
    lzmax(nrank)=real(nzr,mytype)*dz
    
    lxmin=psum(lxmin); lxmax=psum(lxmax)
    lymin=psum(lymin); lymax=psum(lymax)
    lzmin=psum(lzmin); lzmax=psum(lzmax)
    
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
  subroutine particle_domain_check(particle,num_active_particles)
    
    use param,     only : xlx,yly,zlz

    
    type(partype),intent(in),allocatable,target :: particle(:)
    integer,intent(inout),optional :: num_active_particles
    
    ! local data 
    integer :: jpart,npart,psize,jrank,counter,npcanc_totl,nbc,npexit
    type(partype),pointer :: pa

    psize=msize(particle)

    if(present(num_active_particles)) then
        npexit=num_active_particles
    else
        npexit=psize
    endif

    npart=0
    counter=0
    do jpart=1,psize

      pa=>particle(jpart)

      if(pa%inactive) cycle

      npart=npart+1

      if(npart>npexit) exit

      if(pa%x(1)<0._mytype) then

        ! xmin face
        call particle_bc(face=1,bctype=bc_particle(1),particle=pa,particle_deduce=counter)

      elseif(pa%x(1)>xlx) then

        ! xmax face
        call particle_bc(face=2,bctype=bc_particle(2),particle=pa,particle_deduce=counter)
      
      endif

      if(pa%x(2)<0._mytype) then

        ! ymin face
        call particle_bc(face=3,bctype=bc_particle(3),particle=pa,particle_deduce=counter)

      elseif(pa%x(2)>yly) then

        ! ymax face
        call particle_bc(face=4,bctype=bc_particle(4),particle=pa,particle_deduce=counter)
      
      endif

      if(pa%x(3)<0._mytype) then

        ! zmin face
        call particle_bc(face=5,bctype=bc_particle(5),particle=pa,particle_deduce=counter)

      elseif(pa%x(3)>zlz) then

        ! zmax face
        call particle_bc(face=6,bctype=bc_particle(6),particle=pa,particle_deduce=counter)
      
      endif

      if(pa%inactive) cycle 
      ! if the particle is inactive due the implementation of boundary condtion.

      ! checking local domain boundary 
      if( pa%x(1)>=lxmin(nrank) .and. pa%x(1)<lxmax(nrank) .and. &
          pa%x(2)>=lymin(nrank) .and. pa%x(2)<lymax(nrank) .and. &
          pa%x(3)>=lzmin(nrank) .and. pa%x(3)<lzmax(nrank) ) then
          continue
      else
    
        pa%swap=.true.
    
        do jrank=0,nproc-1
    
          ! to find which rank the particle are moving to and 
          ! mark
          if(jrank==nrank) cycle
    
          if( pa%x(1)>=lxmin(jrank) .and. pa%x(1)<lxmax(jrank) .and. &
              pa%x(2)>=lymin(jrank) .and. pa%x(2)<lymax(jrank) .and. &
              pa%x(3)>=lzmin(jrank) .and. pa%x(3)<lzmax(jrank) ) then
    
            pa%rank2go=jrank
    
            exit
    
          endif
    
        enddo
    
      endif

    enddo

    if(present(num_active_particles)) then
      num_active_particles=num_active_particles-counter
    endif
    
  end subroutine particle_domain_check
  !+-------------------------------------------------------------------+
  ! The end of the subroutine particle_domain_check                    |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| This subroutine is to apply b.c. to particles.                    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 03-10-2024  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine particle_bc(face,bctype,particle,particle_deduce)

    use param,     only : xlx,yly,zlz

    ! arguments
    integer,intent(in) :: face
    character(len=*),intent(in) :: bctype
    type(partype),intent(inout) :: particle
    integer,intent(inout) :: particle_deduce

    ! local data
    real(mytype), dimension(6), save :: bcord, lenpe
    logical,save :: firstcal=.true.
    integer :: idir
    
    if(firstcal) then

        bcord(1)=0.0_mytype
        bcord(2)=xlx
        bcord(3)=0.0_mytype
        bcord(4)=yly
        bcord(5)=0.0_mytype
        bcord(6)=zlz

        ! lenpe is to get the particle back to the domain for periodic boundaries.
        ! defined as the distance (+ or -) between the two paring periodic boundaries.
        lenpe(1)=xlx
        lenpe(2)=-xlx
        lenpe(3)=yly
        lenpe(4)=-yly
        lenpe(5)=zlz
        lenpe(6)=-zlz

        firstcal=.false.
    endif

    idir=(face+1)/2

    if(idir>3) then
        print*,idir,face
        call decomp_2d_abort(1,"idir error @ particle_bc")
    endif

    if(bctype=='periodic') then
        particle%x(idir)=particle%x(idir)+lenpe(face)
    elseif(bctype=='reflective') then
        particle%x(idir)=-particle%x(idir)+2.0_mytype*bcord(face)
    elseif(bctype=='outflow') then
        call particle%reset()
        particle_deduce=particle_deduce+1
    else
        call decomp_2d_abort(1, "!! bc not defined @ particle_bc")
    endif

  end subroutine particle_bc
  !+-------------------------------------------------------------------+
  ! The end of the subroutine particle_bc                              |
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
    
    use decomp_2d, only : xsize
    use var,       only : itime
    use param,     only : re,t,itr
    use mhd,       only : Bm,stuart
    
    ! arguments
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1
    
    ! local data
    integer :: psize,jpart,npart
    type(partype),pointer :: pa
    real(mytype) :: varc,vard,varq(3),maxforce,maxvdiff,vdiff,re_p,bmt(3)
    
    logical,save :: firstcal=.true.
    
    psize=msize(partpack)
    
    npart=0
    
    maxforce=0._mytype
    maxvdiff=0._mytype
    
    call field_var(ux1,uy1,uz1)
    

    do jpart=1,psize
    
      pa=>partpack(jpart)
    
      if(pa%inactive) cycle
    
      pa%v(:) = pa%vf(:)
      pa%f(:) = 0._mytype
    
    enddo
    
    
  end subroutine particle_force
  !+-------------------------------------------------------------------+
  ! The end of the subroutine particle_force                           |
  !+-------------------------------------------------------------------+

end module particle
!+---------------------------------------------------------------------+
!| The end of the module particle                                      |
!+---------------------------------------------------------------------+
