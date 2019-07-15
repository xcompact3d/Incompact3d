PROGRAM incompact3d

  USE decomp_2d
  USE decomp_2d_poisson, ONLY : decomp_2d_poisson_init
  use decomp_2d_io
  USE variables
  USE ibm
  USE param
  USE var
  USE MPI
  USE derivX
  USE derivZ
  USE case
  USE simulation_stats
  USE forces

  USE transeq, ONLY : calculate_transeq_rhs

  implicit none

  integer :: ierr

  integer :: nargin, FNLength, status, DecInd
  logical :: back
  character(len=80) :: InputFN, FNBase

  !!-------------------------------------------------------------------------------
  !! Initialisation
  !!-------------------------------------------------------------------------------

  !! Initialise MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  ! Handle input file like a boss -- GD
  nargin=command_argument_count()
  if (nargin <1) then
     InputFN='input.i3d'
     if (nrank==0) print*, 'Incompact3d is run with the default file -->', InputFN
  elseif (nargin.ge.1) then
     if (nrank==0) print*, 'Program is run with the provided file -->', InputFN

     call get_command_argument(1,InputFN,FNLength,status)
     back=.true.
     FNBase=inputFN((index(InputFN,'/',back)+1):len(InputFN))
     DecInd=index(FNBase,'.',back)
     if (DecInd >1) then
        FNBase=FNBase(1:(DecInd-1))
     end if
  endif

  call parameter(InputFN)

  call decomp_2d_init(nx,ny,nz,p_row,p_col)
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

  call schemes()

  !if (nrank==0) call stabiltemp()

  call decomp_2d_poisson_init()
  call decomp_info_init(nxm,nym,nzm,phG)

  if (ilesmod.ne.0) then
     call filter(0.45_mytype)
     if (jles.le.3)  call init_explicit_les()
  endif

  if (irestart==0) then
     call init(rho1,ux1,uy1,uz1,ep1,phi1,drho1,dux1,duy1,duz1,dphi1,pp3,px1,py1,pz1)
     CALL visu(rho1, ux1, uy1, uz1, pp3(:,:,:,1),phi1, 0)
  else
     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,0)
  endif

  if (iibm.eq.2) then
     call genepsi3d(ep1)
  else if (iibm.eq.1) then
     call body(ux1,uy1,uz1,ep1,0)
  endif

  if (iforces) then
     call init_forces()
     if (irestart==1) call restart_forces(0)
  endif

  call test_speed_min_max(ux1,uy1,uz1)
  if (iscalar==1) call test_scalar_min_max(phi1)

  call simu_stats(1)

  call calc_divu_constraint(divu3, rho1, phi1)
  !!-------------------------------------------------------------------------------
  !! End initialisation
  !!-------------------------------------------------------------------------------
  if(nrank.eq.0)then
     open(42,file='time_evol.dat',form='formatted')
  endif

  do itime=ifirst,ilast
     t=itime*dt
     call simu_stats(2)

     call postprocessing(ux1,uy1,uz1,pp3,phi1,ep1)

     do itr=1,iadvance_time

        !!-------------------------------------------------------------------------
        !! Initialise timestep
        !!-------------------------------------------------------------------------
        call boundary_conditions(rho1,ux1,uy1,uz1,phi1,ep1)
        !!-------------------------------------------------------------------------
        !! End initialise timestep
        !!-------------------------------------------------------------------------

        CALL calculate_transeq_rhs(drho1,dux1,duy1,duz1,dphi1,rho1,ux1,uy1,uz1,ep1,phi1,divu3)

        !!-------------------------------------------------------------------------
        !! Time integrate transport equations
        !!-------------------------------------------------------------------------

        !! XXX N.B. from this point, X-pencil velocity arrays contain momentum.
        call velocity_to_momentum(rho1, ux1, uy1, uz1)

        call intt(rho1,ux1,uy1,uz1,phi1,drho1,dux1,duy1,duz1,dphi1)
        call pre_correc(ux1,uy1,uz1,ep1)
        !!-------------------------------------------------------------------------
        !! End time integrate transport equations
        !!-------------------------------------------------------------------------

        if (iibm==1) then !solid body old school
           call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,1)
           call body(ux1,uy1,uz1,ep1,1)
           call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,2)
        endif

        !!-------------------------------------------------------------------------
        !! Poisson solver and velocity correction
        !!-------------------------------------------------------------------------
        call calc_divu_constraint(divu3, rho1, phi1)
        call solve_poisson(pp3, px1, py1, pz1, rho1, ux1, uy1, uz1, ep1, drho1, divu3)
        call corpg(ux1,uy1,uz1,px1,py1,pz1)

        call momentum_to_velocity(rho1, ux1, uy1, uz1)
        !! XXX N.B. from this point, X-pencil velocity arrays contain velocity.

        !!-------------------------------------------------------------------------
        !! End Poisson solver and velocity correction
        !!-------------------------------------------------------------------------

        if (mod(itime,10)==0) then
           call divergence(dv3,rho1,ux1,uy1,uz1,ep1,drho1,divu3,2)
           call test_speed_min_max(ux1,uy1,uz1)
           if (iscalar==1) call test_scalar_min_max(phi1)
        endif

     enddo !! End sub timesteps

     if(iforces) then
        call force(ux1,uy1,ep1,ta1,tb1,tc1,td1,di1,&
             ux2,uy2,ta2,tb2,tc2,td2,di2)
        if (mod(itime,icheckpoint).eq.0) then
           call restart_forces(1)
        endif
     endif

     !!----------------------------------------------------------------------------
     !! Post-processing / IO
     !!----------------------------------------------------------------------------

     if (mod(itime,icheckpoint).eq.0) then
        call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,1)
     endif

     call simu_stats(3)

     CALL visu(rho1, ux1, uy1, uz1, pp3(:,:,:,1), phi1, itime)
     !!----------------------------------------------------------------------------
     !! End post-processing / IO
     !!----------------------------------------------------------------------------

  enddo !! End time loop

  !!-------------------------------------------------------------------------------
  !! End simulation
  !!-------------------------------------------------------------------------------
  if(nrank.eq.0)then
     close(42)
  endif
  call simu_stats(4)
  call decomp_2d_finalize
  CALL MPI_FINALIZE(ierr)

END PROGRAM incompact3d
