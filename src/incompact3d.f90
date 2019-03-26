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

#ifdef FORCES
  USE forces
#endif
  implicit none

  integer :: code,nlock,i,j,k,ii,bcx,bcy,bcz,fh,ierror
  real(mytype) :: x,y,z,timeleft

  integer :: ErrFlag, nargin, FNLength, status, DecInd, output_counter
  logical :: back
  character(len=80) :: InputFN, FNBase
  character(len=20) :: filename

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Initialisation
  !!-------------------------------------------------------------------------------

  !! Initialise MPI
  call MPI_INIT(code)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)

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
     call init(rho1,ux1,uy1,uz1,ep1,phi1,drho1,dux1,duy1,duz1,dphi1)
     CALL visu(rho1, ux1, uy1, uz1, pp3, phi1, 0)
  else
     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3,phi1,dphi1,px1,py1,pz1,0)
  endif

  if (iibm.eq.2) then
     call genepsi3d(ep1)
  else if (iibm.eq.1) then
     call body(ux1,uy1,uz1,ep1,0)
  endif

#ifdef FORCES
  call init_forces()
  if (irestart==1) call restart_forces(0)
#endif

  call test_speed_min_max(ux1,uy1,uz1)
  if (iscalar==1) call test_scalar_min_max(phi1)

  call simu_stats(1)

  call calc_divu_constraint(divu3, rho1)
  !!-------------------------------------------------------------------------------
  !! End initialisation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(nrank.eq.0)then
     open(42,file='time_evol.dat',form='formatted')
  endif

  do itime=ifirst,ilast
     t=itime*dt
     call simu_stats(2)
     
     call postprocessing(ux1,uy1,uz1,phi1,ep1)
     
     do itr=1,iadvance_time

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Initialise timestep
        !!-------------------------------------------------------------------------
        call boundary_conditions(ux1,uy1,uz1,phi1,ep1)
        !!-------------------------------------------------------------------------
        !! End initialise timestep
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL calculate_transeq_rhs(drho1,dux1,duy1,duz1,dphi1,rho1,ux1,uy1,uz1,ep1,phi1,divu3)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Time integrate transport equations
        !!-------------------------------------------------------------------------
        if (ilmn) then
           !! XXX N.B. from this point, X-pencil velocity arrays contain momentum.
           call primary_to_conserved(rho1, ux1)
           call primary_to_conserved(rho1, uy1)
           call primary_to_conserved(rho1, uz1)
        endif
        
        call intt(rho1,ux1,uy1,uz1,phi1,drho1,dux1,duy1,duz1,dphi1)
        call pre_correc(ux1,uy1,uz1,ep1)
        !!-------------------------------------------------------------------------
        !! End time integrate transport equations
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (iibm==1) then !solid body old school
           call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,1)
           call body(ux1,uy1,uz1,ep1,1)
           call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,2)
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Poisson solver and velocity correction
        !!-------------------------------------------------------------------------
        call calc_divu_constraint(divu3, rho1)
        call solve_poisson(pp3, rho1, ux1, uy1, uz1, ep1, drho1, divu3)
        call gradp(px1,py1,pz1,pp3)
        call corpg(ux1,uy1,uz1,px1,py1,pz1)

        if (ilmn) then
           !! XXX N.B. from this point, X-pencil velocity arrays contain velocity.
           call conserved_to_primary(rho1, ux1)
           call conserved_to_primary(rho1, uy1)
           call conserved_to_primary(rho1, uz1)
        endif
        !!-------------------------------------------------------------------------
        !! End Poisson solver and velocity correction
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (mod(itime,10)==0) then
           call divergence(dv3,rho1,ux1,uy1,uz1,ep1,drho1,divu3,2)
           call test_speed_min_max(ux1,uy1,uz1)
           if (iscalar==1) call test_scalar_min_max(phi1)
        endif

     enddo !! End sub timesteps
#ifdef FORCES
     call force(ux1,uy1,ep1,ta1,tb1,tc1,td1,di1,&
          ux2,uy2,ta2,tb2,tc2,td2,di2)
     if (mod(itime,icheckpoint).eq.0) then
        call restart_forces(1)
     endif
#endif

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !! Post-processing / IO
     !!----------------------------------------------------------------------------

     if (mod(itime,icheckpoint).eq.0) then
        call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3,phi1,dphi1,px1,py1,pz1,1)
     endif

     call simu_stats(3)

     CALL visu(rho1, ux1, uy1, uz1, pp3, phi1, itime)
     !!----------------------------------------------------------------------------
     !! End post-processing / IO
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  enddo !! End time loop

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! End simulation
  !!-------------------------------------------------------------------------------
  if(nrank.eq.0)then
     close(42)
  endif
  call simu_stats(4)
  call decomp_2d_finalize
  CALL MPI_FINALIZE(code)

END PROGRAM incompact3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  SUBROUTINE: calculate_transeq_rhs
!!      AUTHOR: Paul Bartholomew
!! DESCRIPTION: Calculates the right hand sides of all transport
!!              equations - momentum, scalar transport, etc.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calculate_transeq_rhs(drho1,dux1,duy1,duz1,dphi1,rho1,ux1,uy1,uz1,ep1,phi1,divu3)

  USE decomp_2d, ONLY : mytype, xsize, zsize
  USE variables, ONLY : numscalar
  USE param, ONLY : ntime, ilmn, nrhotime
  USE transeq, ONLY : momentum_rhs_eq, continuity_rhs_eq, scalar
  
  IMPLICIT NONE

  !! Inputs
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ux1, uy1, uz1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime), INTENT(IN) :: rho1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), numscalar), INTENT(IN) :: phi1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ep1
  REAL(mytype), DIMENSION(zsize(1), zsize(2), zsize(3)), INTENT(IN) :: divu3

  !! Outputs
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime) :: drho1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime, numscalar) :: dphi1
  
  !! Momentum equations
  CALL momentum_rhs_eq(dux1,duy1,duz1,rho1,ux1,uy1,uz1,ep1,phi1,divu3)

  !! Scalar equations
  !! XXX Not yet LMN!!!
  CALL scalar(dphi1, ux1, uy1, uz1, phi1)

  !! Other (LMN, ...)
  IF (ilmn) THEN
     CALL continuity_rhs_eq(drho1, rho1, ux1, uy1, uz1, divu3)
  ENDIF
  
END SUBROUTINE calculate_transeq_rhs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  SUBROUTINE: solve_poisson
!!      AUTHOR: Paul Bartholomew
!! DESCRIPTION: Takes the intermediate momentum field as input,
!!              computes div and solves pressure-Poisson equation.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE solve_poisson(pp3, rho1, ux1, uy1, uz1, ep1, drho1, divu3)

  USE decomp_2d, ONLY : mytype, xsize, zsize, ph1
  USE decomp_2d_poisson, ONLY : poisson
  USE var, ONLY : nzmsize
  USE param, ONLY : ntime, nrhotime

  IMPLICIT NONE

  INTEGER :: nlock

  !! Inputs
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ux1, uy1, uz1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ep1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime), INTENT(IN) :: rho1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime), INTENT(IN) :: drho1
  REAL(mytype), DIMENSION(zsize(1), zsize(2), zsize(3)), INTENT(IN) :: divu3

  !! Outputs
  REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize) :: pp3

  nlock = 1 !! Corresponds to computing div(u*)
  CALL divergence(pp3,rho1,ux1,uy1,uz1,ep1,drho1,divu3,nlock)

  CALL poisson(pp3)

END SUBROUTINE solve_poisson

SUBROUTINE visu(rho1, ux1, uy1, uz1, pp3, phi1, itime)

  USE decomp_2d, ONLY : mytype, xsize, ysize, zsize
  USE decomp_2d, ONLY : fine_to_coarseV, transpose_z_to_y, transpose_y_to_x
  USE decomp_2d_io, ONLY : decomp_2d_write_one
  USE param, ONLY : ivisu, ioutput, nrhotime, ilmn, iscalar
  USE variables, ONLY : sx, cifip6, cisip6, ciwip6, cifx6, cisx6, ciwx6
  USE variables, ONLY : sy, cifip6y, cisip6y, ciwip6y, cify6, cisy6, ciwy6
  USE variables, ONLY : sz, cifip6z, cisip6z, ciwip6z, cifz6, cisz6, ciwz6
  USE variables, ONLY : numscalar
  USE var, ONLY : uvisu
  USE var, ONLY : pp1, ta1, di1, nxmsize
  USE var, ONLY : pp2, ppi2, dip2, ph2, nymsize
  USE var, ONLY : ppi3, dip3, ph3, nzmsize

  IMPLICIT NONE

  CHARACTER(len=30) :: filename
  
  !! Inputs
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ux1, uy1, uz1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime), INTENT(IN) :: rho1
  REAL(mytype), DIMENSION(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize), INTENT(IN) :: pp3
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), numscalar), INTENT(IN) :: phi1
  INTEGER, INTENT(IN) :: itime

  INTEGER :: is

  IF ((ivisu.NE.0).AND.(MOD(itime, ioutput).EQ.0)) THEN
     !! Write velocity
     uvisu=0.
     call fine_to_coarseV(1,ux1,uvisu)
990  format('ux',I3.3)
     write(filename, 990) itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)
     
     uvisu=0.
     call fine_to_coarseV(1,uy1,uvisu)
991  format('uy',I3.3)
     write(filename, 991) itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)
     
     uvisu=0.
     call fine_to_coarseV(1,uz1,uvisu)
992  format('uz',I3.3)
     write(filename, 992) itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)

     !WORK Z-PENCILS
     call interzpv(ppi3,pp3,dip3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
          (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
     !WORK Y-PENCILS
     call transpose_z_to_y(ppi3,pp2,ph3) !nxm nym nz
     call interypv(ppi2,pp2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
          (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
     !WORK X-PENCILS
     call transpose_y_to_x(ppi2,pp1,ph2) !nxm ny nz
     call interxpv(ta1,pp1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
          nxmsize,xsize(1),xsize(2),xsize(3),1)

     uvisu=0._mytype
     call fine_to_coarseV(1,ta1,uvisu)
993  format('pp',I3.3)
     write(filename, 993) itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)

     !! LMN - write out density
     IF (ilmn) THEN
        uvisu=0.
        call fine_to_coarseV(1,rho1(:,:,:,1),uvisu)
994     format('rho',I3.3)
        write(filename, 994) itime/ioutput
        call decomp_2d_write_one(1,uvisu,filename,2)
     ENDIF

     !! Scalars
     IF (iscalar.NE.0) THEN
995     format('phi',I1.1,I3.3)
        DO is = 1, numscalar
           uvisu=0.
           call fine_to_coarseV(1,phi1(:,:,:,is),uvisu)
           write(filename, 995) is, itime/ioutput
           call decomp_2d_write_one(1,uvisu,filename,2)
        ENDDO
     ENDIF
  ENDIF
END SUBROUTINE visu

SUBROUTINE intt(rho1, ux1, uy1, uz1, phi1, drho1, dux1, duy1, duz1, dphi1)

  USE decomp_2d, ONLY : mytype, xsize
  USE param, ONLY : ntime, nrhotime, ilmn, iscalar
  USE variables, ONLY : numscalar

  IMPLICIT NONE

  !! INPUT/OUTPUT
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime) :: drho1, dux1, duy1, duz1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime, numscalar) :: dphi1

  !! OUTPUT
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), numscalar) :: phi1

  !! LOCAL
  INTEGER :: is
  
  CALL int_time_momentum(ux1, uy1, uz1, dux1, duy1, duz1)

  IF (iscalar.NE.0) THEN
     DO is = 1, numscalar
        CALL int_time(phi1(:,:,:,is), dphi1(:,:,:,:,is))
     ENDDO
  ENDIF

  IF (ilmn) THEN
     CALL int_time_continuity(rho1, drho1)
  ENDIF
  
ENDSUBROUTINE intt
