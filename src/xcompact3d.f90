!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

program xcompact3d

  use var
  use case

  use transeq, only : calculate_transeq_rhs
  use time_integrators, only : int_time
  use navier, only : velocity_to_momentum, momentum_to_velocity, pre_correc, &
       calc_divu_constraint, solve_poisson, cor_vel
  use tools, only : restart, simu_stats, apply_spatial_filter, read_inflow
  use turbine, only : compute_turbines
  use ibm_param
  use ibm, only : body
  use genepsi, only : genepsi3d
  use ellipsoid_utils, only: lin_step, ang_step, QuaternionNorm
  use forces, only : force, init_forces, iforces,update_forces, xld,xrd,yld,yud,zld,zrd,torque_calc,nvol
  implicit none
  real(mytype)  :: dummy,drag(10),lift(10),lat(10),grav_effy(10),grav_effx(10),grav_effz(10),xtorq(10),ytorq(10),ztorq(10),maxrad
  integer :: iounit,ierr,i
  real, dimension(100) :: x
  character(len=30) :: filename!, filename2



  call init_xcompact3d()

  iounit = 135
  !Print forces out on ellip
  if ((nrank==0).and.(force_csv.eq.1)) then 
   open(unit=20, file='force_out.dat', status='unknown',form='formatted')
   ! if (ierr /= 0) then
   !    print *, 'Error opening file.'
   !    stop
   ! end if
   write(*,*) 'Outputting forces' 
  end if

  if (nrank==0) then
   do i = 1,nbody
      write(filename,"('body.dat',I1.1)") i
      open(unit=11+i, file=filename, status='unknown', form='formatted')
   enddo
  endif
!   do i = 1,100
!    x(i) = i
!   enddo
!   open(unit=3, file='testcsv.dat', status='new',action='write',iostat=ierr)

!   do i = 1,100
!    write(3,*) x(i)
!   enddo

  

  do itime=ifirst,ilast
     !t=itime*dt
     t=t0 + (itime0 + itime + 1 - ifirst)*dt
     call simu_stats(2)

     if (iturbine.ne.0) call compute_turbines(ux1, uy1, uz1)

     if (iin.eq.3.and.mod(itime,ntimesteps)==1) then
        call read_inflow(ux_inflow,uy_inflow,uz_inflow,itime/ntimesteps)
     endif

     if ((itype.eq.itype_abl.or.iturbine.ne.0).and.(ifilter.ne.0).and.(ilesmod.ne.0)) then
        call filter(C_filter)
        call apply_spatial_filter(ux1,uy1,uz1,phi1)
     endif


     do itr=1,iadvance_time

        call set_fluid_properties(rho1,mu1)
        call boundary_conditions(rho1,ux1,uy1,uz1,phi1,ep1)

        if (imove.eq.1) then ! update epsi for moving objects
          if ((iibm.eq.2).or.(iibm.eq.3)) then
             call genepsi3d(ep1)
             do i = 1,nobjmax 
               maxrad = max(shape(i,1),shape(i,2),shape(i,3))
               if (iforces.eq.1) then
                  xld(i) = position(i,1) - maxrad * ra(i) * cvl_scalar
                  xrd(i) = position(i,1) + maxrad * ra(i) * cvl_scalar
                  yld(i) = position(i,2) - maxrad * ra(i) * cvl_scalar
                  yud(i) = position(i,2) + maxrad * ra(i) * cvl_scalar
                  zld(i) = position(i,3) - maxrad * ra(i) * cvl_scalar
                  zrd(i) = position(i,3) + maxrad * ra(i) * cvl_scalar
                  ! write(*,*) "CV bounds = ", xld(i), xrd(i), yld(i), yud(i), zld(i), zrd(i)
                  
               endif
            enddo
            if (itime.eq.ifirst) then 
               call init_forces()
            else 
               call update_forces()
            endif
          else if (iibm.eq.1) then
             call body(ux1,uy1,uz1,ep1)
          endif
        endif
        call calculate_transeq_rhs(drho1,dux1,duy1,duz1,dphi1,rho1,ux1,uy1,uz1,ep1,phi1,divu3)
#ifdef DEBG
        call check_transients()
#endif
        
        if (ilmn) then
           !! XXX N.B. from this point, X-pencil velocity arrays contain momentum (LMN only).
           call velocity_to_momentum(rho1,ux1,uy1,uz1)
        endif

        call int_time(rho1,ux1,uy1,uz1,phi1,drho1,dux1,duy1,duz1,dphi1)
        call pre_correc(ux1,uy1,uz1,ep1)

        call calc_divu_constraint(divu3,rho1,phi1)
        call solve_poisson(pp3,px1,py1,pz1,rho1,ux1,uy1,uz1,ep1,drho1,divu3)
        call cor_vel(ux1,uy1,uz1,px1,py1,pz1)

        if (ilmn) then
           call momentum_to_velocity(rho1,ux1,uy1,uz1)
           !! XXX N.B. from this point, X-pencil velocity arrays contain velocity (LMN only).
           !! Note - all other solvers work on velocity always
        endif
        
        call test_flow(rho1,ux1,uy1,uz1,phi1,ep1,drho1,divu3)

        !Add force calculation here
      !   if (nrank.eq.0) then 
      !   write(*,*) 'Going to call force from xcompact3d, itr = ', itr
      !   endif 
        call force(ux1,uy1,uz1,ep1,drag,lift,lat,1)
        grav_effx = grav_x*(rho_s-1.0)
        grav_effy = grav_y*(rho_s-1.0)
        grav_effz = grav_z*(rho_s-1.0)
        do i = 1,nbody
         linearForce(i,:) = [drag(i)-grav_effx(i), lift(i)-grav_effy(i), lat(i)-grav_effz(i)]
        enddo
        if (nozdrift==1) then
            linearForce(:,3)=zero
        endif

        if (bodies_fixed==1) then
            linearForce(:,:)=zero
        endif

        if ((nrank==0).and.(force_csv.eq.1)) then
         ! open(unit=20, file='force_out.dat', action='write')
         write(20, *) linearForce(1,1), linearForce(1,2), linearForce(1,3)
         write(*,*) 'Writing forces', linearForce(1,1), linearForce(1,2), linearForce(1,3)
         flush(20)
        endif 
        
        
        if (torques_flag.eq.1) then 
         call torque_calc(ux1,uy1,uz1,ep1,xtorq,ytorq,ztorq,1)
        endif 
        if (orientations_free.eq.1) then 
         do i = 1,nvol 
            torque(i,:) = [xtorq(i), ytorq(i), ztorq(i)]
         enddo
         if (ztorq_only.eq.1) then
            torque(:,1) = zero
            torque(:,2) = zero
         endif
        else 
         torque(:,:) = zero
        endif
      !   if (nrank==0) then

      !   if (bodies_fixed==0) then 
        do i = 1,nvol

         call lin_step(position(i,:),linearVelocity(i,:),linearForce(i,:),ellip_m(i),dt,position_1,linearVelocity_1)
         call ang_step(orientation(i,:),angularVelocity(i,:),torque(i,:),inertia(i,:,:),dt,orientation_1,angularVelocity_1)
        
         position(i,:) = position_1
         linearVelocity(i,:) = linearVelocity_1

         orientation(i,:) = orientation_1
         angularVelocity(i,:) = angularVelocity_1
        enddo


      if (nrank==0) then
         do i = 1,nbody
            write(11+i ,*) t, position(i,1), position(i,2), position(i,3), orientation(i,1), orientation(i,2), orientation(i,3), orientation(i,4), linearVelocity(i,1), linearVelocity(i,2), linearVelocity(i,3), angularVelocity(i,2), angularVelocity(i,3), angularVelocity(i,4), linearForce(i,1), linearForce(i,2), linearForce(i,3), torque(i,1), torque(i,2), torque(i,3)
            flush(11+i)
         enddo
      endif

         if ((nrank==0).and.(mod(itime,ilist)==0)) then 
            do i = 1,nbody
               write(*,*) "Body", i
               write(*,*) "Position =         ", position(i,:)
               write(*,*) "Orientation =      ", orientation(i,:)
               write(*,*) "Linear velocity =  ", linearVelocity(i,:)
               write(*,*) "Angular velocity = ", angularVelocity(i,:)
               write(*,*) "Linear Force = ", linearForce(i,:)
               write(*,*) "Torque = ", torque(i,:)
            enddo
         ! call QuaternionNorm(angularVelocity,dummy)

         ! write(*,*) 'Norm of angvel = ', dummy
         endif   

      !   endif 

      !   if (nrank==0) then 
      !    write(*,*) 'Centroid position is ', position
      !    write(*,*) 'Orientation is ', orientation
      !   end if

     enddo !! End sub timesteps

     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,rho1,drho1,mu1,1)

     call simu_stats(3)

     call postprocessing(rho1,ux1,uy1,uz1,pp3,phi1,ep1)

  enddo !! End time loop

  close(iounit)

  call finalise_xcompact3d()

end program xcompact3d
!########################################################################
!########################################################################
subroutine init_xcompact3d()

  use MPI
  use decomp_2d
  use decomp_2d_io, only : decomp_2d_io_init
  USE decomp_2d_poisson, ONLY : decomp_2d_poisson_init
  use case
  use sandbox, only : init_sandbox
  use forces

  use var

  use navier, only : calc_divu_constraint
  use tools, only : test_speed_min_max, test_scalar_min_max, &
       restart, &
       simu_stats, compute_cfldiff, &
       init_inflow_outflow

  use param, only : ilesmod, jles,itype
  use param, only : irestart

  use variables, only : nx, ny, nz, nxm, nym, nzm
  use variables, only : p_row, p_col
  use variables, only : nstat, nvisu, nprobe, ilist

  use les, only: init_explicit_les
  use turbine, only: init_turbines

  use visu, only : visu_init, visu_ready

  use genepsi, only : genepsi3d, epsi_init, param_assign
  use ibm, only : body

  use probes, only : init_probes
!   use case, only : param_assign
  implicit none

  integer :: ierr

  integer :: nargin, FNLength, status, DecInd
  logical :: back
  character(len=80) :: InputFN, FNBase
    
  !! Initialise MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  ! Handle input file like a boss -- GD
  nargin=command_argument_count()
  if (nargin <1) then
     InputFN='input.i3d'
     if (nrank==0) write(*,*) 'Xcompact3d is run with the default file -->', trim(InputFN)
  elseif (nargin >= 1) then
     call get_command_argument(1,InputFN,FNLength,status)
     back=.true.
     FNBase=inputFN((index(InputFN,'/',back)+1):len(InputFN))
     DecInd=index(FNBase,'.',back)
     if (DecInd >1) then
        FNBase=FNBase(1:(DecInd-1))
     end if
     if (nrank==0) write(*,*) 'Xcompact3d is run with the provided file -->', trim(InputFN)
  endif

#ifdef ADIOS2
  if (nrank .eq. 0) then
     print *, " WARNING === WARNING === WARNING === WARNING === WARNING"
     print *, " WARNING: Running Xcompact3d with ADIOS2"
     print *, "          this is currently experimental"
     print *, "          for safety of results it is recommended"
     print *, "          to run the default build as this feature"
     print *, "          is developed. Thank you for trying it."
     print *, " WARNING === WARNING === WARNING === WARNING === WARNING"
  endif
#endif
  
  call parameter(InputFN)

  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call decomp_2d_io_init()
  call init_coarser_mesh_statS(nstat,nstat,nstat,.true.)    !start from 1 == true
  call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.)    !start from 1 == true
  call init_coarser_mesh_statP(nprobe,nprobe,nprobe,.true.) !start from 1 == true
  !div: nx ny nz --> nxm ny nz --> nxm nym nz --> nxm nym nzm
  call decomp_info_init(nxm, nym, nzm, ph1)
  call decomp_info_init(nxm, ny, nz, ph4)
  !gradp: nxm nym nzm -> nxm nym nz --> nxm ny nz --> nx ny nz
  call decomp_info_init(nxm, ny, nz, ph2)
  call decomp_info_init(nxm, nym, nz, ph3)

  call init_variables()
  call param_assign()

  call schemes()

  call decomp_2d_poisson_init()
  call decomp_info_init(nxm,nym,nzm,phG)

  if (ilesmod.ne.0) then
     if (jles.gt.0)  call init_explicit_les()
  endif

  if ((iibm.eq.2).or.(iibm.eq.3)) then
   !   call boundary_conditions()
     call genepsi3d(ep1)
  else if (iibm.eq.1) then
     call epsi_init(ep1)
     call body(ux1,uy1,uz1,ep1)
  endif

  if (iforces.eq.1) then
   !   call init_forces()
     if (irestart==1) then
        call restart_forces(0)
     endif
  endif

  !####################################################################
  ! initialise visu
  if (ivisu.ne.0) then
     call visu_init()
     call visu_case_init() !! XXX: If you get error about uninitialised IO, look here.
                           !! Ensures additional case-specific variables declared for IO
     call visu_ready()
  end if
  ! compute diffusion number of simulation
  call compute_cfldiff()
  !####################################################################
  if (irestart==0) then
     call init(rho1,ux1,uy1,uz1,ep1,phi1,drho1,dux1,duy1,duz1,dphi1,pp3,px1,py1,pz1)
     itime = 0
     call preprocessing(rho1,ux1,uy1,uz1,pp3,phi1,ep1)
  else
     itr=1
     if (itype == itype_sandbox) then
        call init_sandbox(ux1,uy1,uz1,ep1,phi1,1)
     end if
     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,rho1,drho1,mu1,0)
  endif

  if ((ioutflow.eq.1).or.(iin.eq.3)) then
     call init_inflow_outflow()
  end if

  if ((iibm.eq.2).or.(iibm.eq.3)) then
     call genepsi3d(ep1)
  else if ((iibm.eq.1).or.(iibm.eq.3)) then
     call body(ux1,uy1,uz1,ep1)
  endif

  if (mod(itime, ilist) == 0 .or. itime == ifirst) then
     call test_speed_min_max(ux1,uy1,uz1)
     if (iscalar==1) call test_scalar_min_max(phi1)
  endif

  call simu_stats(1)

  call calc_divu_constraint(divu3, rho1, phi1)

  call init_probes()

  if (iturbine.ne.0) call init_turbines(ux1, uy1, uz1)

  if (itype==2) then
     if(nrank.eq.0)then
        open(42,file='time_evol.dat',form='formatted')
     endif
  endif
  if (itype==5) then
     if(nrank.eq.0)then
        open(38,file='forces.dat',form='formatted')
     endif
  endif

endsubroutine init_xcompact3d
!########################################################################
!########################################################################
subroutine finalise_xcompact3d()

  use MPI
  use decomp_2d
  use decomp_2d_io, only : decomp_2d_io_finalise

  use tools, only : simu_stats
  use param, only : itype, jles, ilesmod
  use probes, only : finalize_probes
  use visu, only : visu_finalise
  use les, only: finalise_explicit_les

  implicit none

  integer :: ierr
  
  if (itype==2) then
     if(nrank.eq.0)then
        close(42)
     endif
  endif
  if (itype==5) then
     if(nrank.eq.0)then
        close(38)
     endif
  endif
  
  call simu_stats(4)
  call finalize_probes()
  call visu_finalise()
  if (ilesmod.ne.0) then
     if (jles.gt.0) call finalise_explicit_les()
  endif
  call decomp_2d_io_finalise()
  call decomp_2d_finalize
  CALL MPI_FINALIZE(ierr)

endsubroutine finalise_xcompact3d

subroutine check_transients()

  use decomp_2d, only : mytype
  use mpi
  use var
  
  implicit none

  real(mytype) :: dep, dep1
  integer :: code
   
  dep=maxval(abs(dux1))
  call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
  if (nrank == 0) write(*,*)'## MAX dux1 ', dep1
 
  dep=maxval(abs(duy1))
  call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
  if (nrank == 0) write(*,*)'## MAX duy1 ', dep1
 
  dep=maxval(abs(duz1))
  call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
  if (nrank == 0) write(*,*)'## MAX duz1 ', dep1
  
end subroutine check_transients
