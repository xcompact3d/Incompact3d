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
     call init(ux1,uy1,uz1,ep1,phi1,dux1,duy1,duz1,phis1,phiss1)
  else
     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3,phi1,px1,py1,pz1,0)
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
        call boundary_conditions(ux1,uy1,uz1,phi1)
        !!-------------------------------------------------------------------------
        !! End initialise timestep
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL calculate_transeq_rhs(dux1,duy1,duz1,ux1,uy1,uz1,ep1,phi1)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Time integrate transport equations
        !!-------------------------------------------------------------------------
        call int_time_momentum(ux1,uy1,uz1,dux1,duy1,duz1)
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
        call solve_poisson(pp3, ux1, uy1, uz1, ep1)
        call gradp(px1,py1,pz1,pp3)
        call corpg(ux1,uy1,uz1,px1,py1,pz1)
        !!-------------------------------------------------------------------------
        !! End Poisson solver and velocity correction
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (mod(itime,10)==0) then
           call divergence(ux1,uy1,uz1,ep1,dv3,2)
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
        call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3,phi1,px1,py1,pz1,1)
     endif

     call simu_stats(3)

     CALL visu(ux1, uy1, uz1, itime)
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
SUBROUTINE calculate_transeq_rhs(dux1,duy1,duz1,ux1,uy1,uz1,ep1,phi1)

  USE decomp_2d, ONLY : mytype, xsize
  USE variables, ONLY : numscalar
  USE param, ONLY : ntime
  USE transeq, ONLY : momentum_rhs_eq
  
  IMPLICIT NONE

  !! Inputs
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ux1, uy1, uz1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), numscalar), INTENT(IN) :: phi1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ep1

  !! Outputs
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1
  
  !! Momentum equations
  CALL momentum_rhs_eq(dux1,duy1,duz1,ux1,uy1,uz1,ep1,phi1)

  !! Scalar equations

  !! Other (LMN, ...)
  
END SUBROUTINE calculate_transeq_rhs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  SUBROUTINE: solve_poisson
!!      AUTHOR: Paul Bartholomew
!! DESCRIPTION: Takes the intermediate momentum field as input,
!!              computes div and solves pressure-Poisson equation.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE solve_poisson(pp3, ux1, uy1, uz1, ep1)

  USE decomp_2d, ONLY : mytype, xsize, ph1
  USE decomp_2d_poisson, ONLY : poisson
  USE var, ONLY : nzmsize

  IMPLICIT NONE

  INTEGER :: nlock

  !! Inputs
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ux1, uy1, uz1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ep1

  !! Outputs
  REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize) :: pp3

  nlock = 1 !! Corresponds to computing div(u*)
  CALL divergence(ux1,uy1,uz1,ep1,pp3,nlock)

  CALL poisson(pp3)

END SUBROUTINE solve_poisson

SUBROUTINE visu(ux1, uy1, uz1, itime)

  USE decomp_2d, ONLY : mytype, xsize, ysize, zsize
  USE decomp_2d, ONLY : fine_to_coarseV, transpose_x_to_y, transpose_y_to_z, &
       transpose_z_to_y, transpose_y_to_x
  USE decomp_2d_io, ONLY : decomp_2d_write_one
  USE param, ONLY : ivisu, ioutput
  USE variables, ONLY : sx, ffx, fsx, fwx, ffxp, fsxp, fwxp
  USE variables, ONLY : sy, ffy, fsy, fwy, ffyp, fsyp, fwyp, ppy
  USE variables, ONLY : sz, ffz, fsz, fwz, ffzp, fszp, fwzp
  USE variables, ONLY : derx, dery, derz
  USE var, ONLY : uvisu
  USE var, ONLY : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
  USE var, ONLY : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3

  IMPLICIT NONE

  CHARACTER(len=30) :: filename
  
  !! Inputs
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ux1, uy1, uz1
  INTEGER, INTENT(IN) :: itime

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

     !! Calculate and write vorticity

     !x-derivatives
     call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
     call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
     call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
     
     !y-derivatives
     call transpose_x_to_y(ux1,td2)
     call transpose_x_to_y(uy1,te2)
     call transpose_x_to_y(uz1,tf2)
     call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
     call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
     call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
     
     !!z-derivatives
     call transpose_y_to_z(td2,td3)
     call transpose_y_to_z(te2,te3)
     call transpose_y_to_z(tf2,tf3)
     call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
     call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
     call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
     
     !!all back to x-pencils
     call transpose_z_to_y(ta3,td2)
     call transpose_z_to_y(tb3,te2)
     call transpose_z_to_y(tc3,tf2)
     call transpose_y_to_x(td2,tg1)
     call transpose_y_to_x(te2,th1)
     call transpose_y_to_x(tf2,ti1)
     call transpose_y_to_x(ta2,td1)
     call transpose_y_to_x(tb2,te1)
     call transpose_y_to_x(tc2,tf1)
     
     !du/dx=ta1 du/dy=td1 and du/dz=tg1
     !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
     !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1
     
     di1(:,:,:)=sqrt((tf1(:,:,:)-th1(:,:,:))**2+(tg1(:,:,:)-tc1(:,:,:))**2+&
          (tb1(:,:,:)-td1(:,:,:))**2)
     uvisu=0.
     call fine_to_coarseV(1,di1,uvisu)
993  format('vort',I3.3)
     write(filename, 993) itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)

     !! TODO Write pressure
  ENDIF
END SUBROUTINE visu
