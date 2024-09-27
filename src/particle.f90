!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause
module particle

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank,nproc
  use mptool
  use constants, only : pi

  implicit none

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
    ! !
    integer :: num_particle,n
    real(mytype),allocatable,dimension(:) :: x_p

    num_particle=100*(nrank+1)

    allocate(x_p(num_particle))

    do n=1,num_particle
      call random_number(x_p(n))
      x_p(n)=x_p(n)+nrank*100._mytype

      ! print*,nrank,'|',n,'-',x_p(n)
    enddo
    
    ! call pwrite('particle_init.bin',x_p)

    call pread('particle_init.bin',x_p,num_particle)

    print*,nrank,'|',x_p(1)

    ! if(irestart==0) then
    !   !
    !   call particle_gen(particle,numparticle)
    !   !
    !   particle_file_numb=0
    !   !
    !   call h5write_particle()
    !   !
    ! else
    !   ! call h5read_particle(particle,numparticle)
    !   call particle_gen(particle,numparticle)
    !   !
    !   particle_file_numb=0
    !   !
    !   call h5write_particle()
    !   !
    ! endif
    ! !
    ! ! call partical_domain_check('bc_tgv')
    ! ! call partical_domain_check('bc_channel')
    ! call partical_domain_check('out_disappear')
    ! !
    ! !
    ! ! call partical_swap
    ! !
    ! ! call write_particle()
    ! !
    ! part_time=0.d0
    ! part_comm_time=0.d0
    ! part_vel_time=0.d0
    ! part_dmck_time=0.d0
    ! a2a_time=0.d0
    ! count_time=0.d0
    ! table_share_time=0.d0
    ! data_pack_time=0.d0
    ! data_unpack_time=0.d0
    ! mpi_comm_time=0.d0
    ! h5io_time=0.d0
    ! !
    ! particletime=t
    ! sub_time_step=dt
    ! sub_time_step=0.1d0*dt
    !
    ! lorentzforce=.true.
    !
    stop
    !
  end subroutine particle_init
  !+-------------------------------------------------------------------+
  ! The end of the subroutine particle_init                            |
  !+-------------------------------------------------------------------+
end module particle