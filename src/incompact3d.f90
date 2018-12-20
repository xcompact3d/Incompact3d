PROGRAM incompact3d

  USE decomp_2d
  USE decomp_2d_poisson
  use decomp_2d_io
  USE variables
  USE param
  USE var
  USE MPI
  USE derivX
  USE derivZ
#ifdef POST
  USE post_processing
#endif
#ifdef FORCES
  USE forces
#endif
  implicit none

  integer :: code,nlock,i,j,k,ii,bcx,bcy,bcz,fh,ierror
  real(mytype) :: x,y,z,timeleft
  real(8) :: tstart,t1,trank,tranksum,ttotal,tremaining,telapsed

  integer :: ErrFlag, nargin, FNLength, status, DecInd, output_counter
  logical :: back
  character(len=80) :: InputFN, FNBase
  character(len=20) :: filename

  ! Handle input file like a boss -- GD
  nargin=command_argument_count()
  if (nargin <1) then
     InputFN='input.i3d'
     print*, 'Incompact3d is run with the default file -->', InputFN
  elseif (nargin.ge.1) then
     print*, 'Program is run with the provided file -->', InputFN
     
     call get_command_argument(1,InputFN,FNLength,status)
     back=.true.
     FNBase=inputFN((index(InputFN,'/',back)+1):len(InputFN))
     DecInd=index(FNBase,'.',back)
     if (DecInd >1) then
        FNBase=FNBase(1:(DecInd-1))
     end if
 endif
 
 call parameter(InputFN)
  
  call MPI_INIT(code)
  
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

#ifdef IMPLICIT
  iimplicit=1
  if (nrank==0) print *,'--SEMI IMPLICIT CODE (IN BETA)-------------------'
#endif

  call schemes()

#ifdef ELES
  call filter()
#endif

  !if (nrank==0) call stabiltemp()

  call decomp_2d_poisson_init()
  call decomp_info_init(nxm,nym,nzm,phG)

#ifdef ELES
  call init_explicit_les()
#endif

  if (irestart==0) then
     call init(ux1,uy1,uz1,ep1,phi1,dux1,duy1,duz1,phis1,phiss1)
  else
     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3,phi1,px1,py1,pz1,0)
  endif


#ifdef IBM
  if(ilag.eq.1)  call genepsi3d(ep1)
  if(ivirt.eq.1) call body(ux1,uy1,uz1,ep1,0)
#endif

#ifdef POST
  call init_post(ep1)
#endif

#ifdef FORCES
  call init_forces()
  if (irestart==1) call restart_forces(0)
#endif

  if (irestart==0) then
     itime=0
#ifdef VISU
     call VISU_INSTA (ux1,uy1,uz1,phi1,ep1,diss1,.false.)
#endif
#ifdef POST
     call postprocessing(ux1,uy1,uz1,phi1,ep1)
#endif
  endif

  call test_speed_min_max(ux1,uy1,uz1)
  if (iscalar==1) call test_scalar_min_max(phi1)


  tstart=zero;t1=zero;trank=zero;tranksum=zero;ttotal=zero
  call cpu_time(tstart)

 

  do itime=ifirst,ilast
     t=itime*dt
     call cpu_time(t1)
     if (nrank==0) then
        print *,'-----------------------------------------------------------'
        write(*,"(' Time step =',i7,'/',i7,', Time unit =',F9.4)") itime,ilast,t
     endif

     do itr=1,iadvance_time

        call boundary_conditions(ux1,uy1,uz1,phi1,ep1)
        call momentum_rhs_eq(dux1,duy1,duz1,ux1,uy1,uz1,ep1,phi1)
        call int_time_momentum(ux1,uy1,uz1,dux1,duy1,duz1)
        call pre_correc(ux1,uy1,uz1,ep1)
#ifdef IBM
        if (ivirt==1) then !solid body old school
           !we are in X-pencil
           call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,1)
           call body(ux1,uy1,uz1,ep1,1)
           call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,2)
        endif
#endif
        call divergence(ux1,uy1,uz1,ep1,pp3,1)
        call poisson(pp3)
        call gradp(px1,py1,pz1,pp3)
        call corpg(ux1,uy1,uz1,px1,py1,pz1)

        if (mod(itime,itest)==0) then
           call divergence(ux1,uy1,uz1,ep1,dv3,2)
           call test_speed_min_max(ux1,uy1,uz1)
           if (iscalar==1) call test_scalar_min_max(phi1)
        endif

     enddo

#ifdef STATS
     call CONVERGENCE_STATISTIC(ux1,ep1,u1sum_tik,u1sum_tak,tsum)
     call CONVERGENCE_STATISTIC2(ux1,ep1,tik1,tik2,tak1,tak2)
#endif
#ifdef FORCES
     call force(ux1,uy1,ep1,ta1,tb1,tc1,td1,di1,&
          ux2,uy2,ta2,tb2,tc2,td2,di2)
     if (mod(itime,icheckpoint).eq.0) then
        call restart_forces(1)
     endif
#endif
#ifdef POST
     call postprocessing(ux1,uy1,uz1,phi1,ep1)
     if (nprobes.ne.0) call write_probes(ux1,uy1,uz1,phi1)
#endif

     call cpu_time(trank)
     if (nrank==0) print *,'Time per this time step (s):',real(trank-t1)
     if (nrank==0) print *,'Snapshot current/final ',int(itime/ioutput),int(ilast/ioutput)

#ifdef VISU
     if (mod(itime,ioutput).eq.0) then
        call VISU_INSTA (ux1,uy1,uz1,phi1,ep1,.false.)

        if (save_pre.eq.1.OR.save_prem.eq.1) call VISU_PRE (pp3,ta1,tb1,di1,&
             ta2,tb2,di2,ta3,di3,nxmsize,nymsize,nzmsize,uvisu,pre1)
     endif
#endif

#ifdef STATS
     if (t.ge.callstat) then
        call cpu_time(t1)
        call OVERALL_STATISTIC(ux1,uy1,uz1,phi1,pre1,diss1,ta1,&
             u1sum,v1sum,w1sum,u2sum,v2sum,w2sum,&
             u3sum,v3sum,w3sum,u4sum,v4sum,w4sum,&
             uvsum,vwsum,vwsum,disssum,tsum,&
             psum,ppsum,upsum,vpsum,wpsum)

        if (save_dudx.eq.1) call EXTRA_STAT (ux1,uy1,uz1,uvisu,tsum,dudxsum,utmapsum)
        call cpu_time(trank)
        if (nrank==0) print *,'Time in statistics (s):',real(trank-t1,4)

     endif
#endif

     if (mod(itime,icheckpoint).eq.0) then
        call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3,phi1,px1,py1,pz1,1)
     endif

     call cpu_time(trank)

     telapsed = (trank-tstart)/thirtysixthousand
     tremaining  = telapsed*(ilast-itime)/(itime-ifirst)

     if (nrank==0) then
        write(*,"(' Remaining time:',I8,' h ',I2,' min')") int(tremaining), int((tremaining-int(tremaining))*sixty)
        write(*,"(' Elapsed time:  ',I8,' h ',I2,' min')") int(telapsed), int((telapsed-int(telapsed))*sixty)
     endif
  enddo

  call cpu_time(trank); ttotal=trank-tstart

  if (nrank==0) then
     print *,'==========================================================='
     print *,''
     print *,'Good job! Incompact3d finished successfully!'
     print *,''
     print *,'2DECOMP with p_row*p_col=',p_row,p_col
     print *,''
     print *,'nx*ny*nz=',nx*ny*nz
     print *,'nx,ny,nz=',nx,ny,nz
     print *,'dx,dy,dz=',dx,dy,dz
     print *,''
     print *,'Averaged time per step (s):',real(ttotal/(ilast-ifirst),4)
     print *,'Total wallclock (s):',real(ttotal,4)
     print *,'Total wallclock (m):',real(ttotal/sixty,4)
     print *,'Total wallclock (h):',real(ttotal/thirtysixthousand,4)
     print *,'Total wallclock (d):',real(ttotal*1.1574e-5,4)
     print *,''
  endif

  call decomp_2d_finalize
  CALL MPI_FINALIZE(code)

end PROGRAM incompact3d
