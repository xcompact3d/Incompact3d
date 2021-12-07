module x3d_unit_testing_tools

   use MPI
   use iso_c_binding

   implicit none

   logical, save :: check_error

   private
   public :: log, check_error, one_time_step

contains

   ! Save given error values in a text file
   subroutine log(filename, l1, l2, linf, time)

      use decomp_2d, only : mytype, nrank, nproc, real_type, decomp_2d_abort

      implicit none

      character(len=*), intent(in) :: filename
      real(mytype), intent(in) :: l1, l2, linf, time

      logical :: file_found
      integer :: iounit, code
      real(mytype), dimension(nproc) :: alltime

      alltime(:) = 0.d0
      alltime(nrank+1) = time
      call MPI_ALLREDUCE(MPI_IN_PLACE,alltime,nproc,real_type,MPI_SUM,MPI_COMM_WORLD,code)
      if (code /= 0) call decomp_2d_abort(code, "MPI_ALLREDUCE")

      if (nrank == 0) then
         write(*,*) trim(filename)//" : ", l1, l2, linf, alltime
         inquire(file=filename//".dat", exist=file_found)
         if (file_found) then
            open(newunit=iounit, file=trim(filename)//".dat", form='formatted', status='old', access='append')
         else
            open(newunit=iounit, file=trim(filename)//".dat", form='formatted', status='new')
         endif
         write(iounit, *) l1, l2, linf, alltime
         close(iounit)
      endif

   end subroutine log

   ! Performs one time step, duplicated from src/xcompact3d.f90
   ! TODO: avoid duplicate
   subroutine one_time_step()

      use var
      use case

      use param, only : catching_signal
      use transeq, only : calculate_transeq_rhs
      use time_integrators, only : int_time
      use navier, only : velocity_to_momentum, momentum_to_velocity, pre_correc, &
           calc_divu_constraint, solve_poisson, cor_vel
      use tools, only : restart, simu_stats, apply_spatial_filter, read_inflow
      use turbine, only : compute_turbines
      use ibm, only : body, imove
      use genepsi, only : genepsi3d

      implicit none

      !t=itime*dt
      t=t0 + (itime0 + itime + 1 - ifirst)*dt
      call simu_stats(2)

      if (iturbine /= 0) call compute_turbines(ux1, uy1, uz1)

      if (iin == 3.and.mod(itime,ntimesteps)==1) then
         call read_inflow(ux_inflow,uy_inflow,uz_inflow,itime/ntimesteps)
      endif

      if ((itype == itype_abl.or.iturbine /= 0).and.(ifilter /= 0).and.(ilesmod /= 0)) then
         call filter(C_filter)
         call apply_spatial_filter(ux1,uy1,uz1,phi1)
      endif

      do itr=1,iadvance_time

         call set_fluid_properties(rho1,mu1)
         call boundary_conditions(rho1,ux1,uy1,uz1,phi1,ep1)

         if (imove == 1) then ! update epsi for moving objects
            if ((iibm == 2).or.(iibm == 3)) then
               call genepsi3d(ep1)
            else if (iibm == 1) then
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

      enddo !! End sub timesteps

      call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,rho1,drho1,mu1,1)

      call simu_stats(3)

      call postprocessing(rho1,ux1,uy1,uz1,pp3,phi1,ep1)

   end subroutine one_time_step

end module x3d_unit_testing_tools
