module actuator_line_model

    use decomp_2d, only: mytype, nrank
    use variables, only : ilist
    use param, only: itime
    use actuator_line_model_utils
    use airfoils
    use actuator_line_element
    use actuator_line_turbine
    use actuator_line_controller
    use actuator_line_write_output
    use dynstall

    implicit none

    type(ActuatorLineType), allocatable, save :: Actuatorline(:)
    type(TurbineType), allocatable, save :: Turbine(:)
    integer, save :: Ntur, Nal
    real(mytype), save :: DeltaT, Visc, ctime
    logical,save :: actuator_line_model_writeFlag=.true.

    public :: actuator_line_model_init, actuator_line_model_write_output, &
              actuator_line_model_update, actuator_line_model_compute_forces

contains

    !*******************************************************************************
    !
    subroutine actuator_line_model_init(Nturbines,Nactuatorlines,turbines_file,actuatorlines_file,dt)
    !
    !*******************************************************************************

      use param, only: irestart, iturboutput

      implicit none
      integer :: Nturbines, Nactuatorlines
      character(len=80), dimension(100), intent(in) :: turbines_file, actuatorlines_file
      real(mytype), intent(in) :: dt
      integer :: itur,ial
      character(1000) :: ReadLine
      character(100) :: filename, dir

      if (nrank==0) then
         write(*,*) '==========================================================='
         write(*,*) 'Initializing the Actuator Line Model'
         write(*,*) 'Developed by G. Deskos 2015-2017'
         write(*,*) '==========================================================='
      endif

      ! Specify turbines
      Ntur=Nturbines
      if (nrank==0) then
         write(*,*) 'Number of turbines : ', Ntur
      endif
      call get_turbine_options(turbines_file)
      if (Ntur>0) then
         if (irestart==1) then
            ! Read the checkpoint information and rotate actuator lines accordingly
            write(filename,"('restart',I7.7'.alm')") itime
            open(17,file=filename)
            read(17,*)
            do itur=1,Ntur
               read(17,'(A)') ReadLine
               read(ReadLine,*) turbine(itur)%AzimAngle, turbine(itur)%angularVel, turbine(itur)%cbp_old
            enddo
            close(17)
            ! Read statistics
            dir=adjustl(trim(dirname(itime/iturboutput)))
            do itur=1,Ntur
               open(17,File=trim(dir)//'_'//trim(turbine(itur)%name)//'.stat')
               read(17,*)
               read(17,*) turbine(itur)%T_ave, turbine(itur)%P_ave, turbine(itur)%Torque_ave 
               close(17)
            enddo
         endif
         do itur=1,Ntur
            call set_turbine_geometry(Turbine(itur))
         enddo

      endif

      ! Speficy actuator lines
      Nal=NActuatorlines
      if (nrank==0) then
         write(*,*) 'Number of actuator lines : ', Nal
      endif
      call get_actuatorline_options(actuatorlines_file)
      if (Nal>0) then
         do ial=1,Nal
            call set_actuatorline_geometry(actuatorline(ial))
         enddo
      endif

      ! Initialize time
      ctime=0.0
      DeltaT=dt

    end subroutine actuator_line_model_init

    !*******************************************************************************
    !
    subroutine actuator_line_model_write_output(dump_no)
    !
    !*******************************************************************************

      implicit none
      integer, intent(in) :: dump_no
      integer :: itur,ial,iblade
      character(len=100) :: dir

      dir=adjustl(trim(dirname(dump_no)))

      if (Ntur>0) then
         do itur=1,Ntur
            call actuator_line_turbine_write_output(turbine(itur),dir)
            call actuator_line_turbine_write_statistics(turbine(itur),dir)
            do iblade=1,turbine(itur)%NBlades
               call actuator_line_element_write_output(turbine(itur)%Blade(iblade),dir)
               if (turbine(itur)%Blade(iblade)%do_Sheng_stall) then
                  call dynamic_stall_write_output(turbine(itur)%Blade(iblade),dir)
               endif
            enddo
            if (turbine(itur)%Has_Tower) then
               call actuator_line_element_write_output(turbine(itur)%tower,dir)
            endif
         enddo
      endif

      if (Nal>0) then
         do ial=1,Nal
            call actuator_line_element_write_output(actuatorline(ial),dir)
            if (actuatorline(ial)%do_Sheng_stall) then
               call dynamic_stall_write_output(actuatorline(ial),dir)
            endif
         enddo
      endif

    end subroutine actuator_line_model_write_output

    !*******************************************************************************
    !
    subroutine get_turbine_options(turbines_path)
    !
    !*******************************************************************************

      use decomp_2d, only: mytype
      use param, only: u1,u2
      use param, only: zero, zpthree, zptwoone, zpfive, one, twentyone
      use constants

      implicit none
      character(len=80),dimension(100),intent(in) :: turbines_path
      integer :: i,j,k
      integer :: nfoils
      character(len=100) :: name, blade_geom, tower_geom, dynstall_param_file, AeroElastInputFile, AeroElastSolverFile, list_controller_file
      character(len=100), dimension(MaxNAirfoils) :: afname
      real(mytype), dimension(3) :: origin
      integer :: numblades,numfoil,towerFlag, TypeFlag, OperFlag, RotFlag, AddedMassFlag, DynStallFlag, EndEffectsFlag
      integer :: TipCorr, RootCorr, RandomWalkForcingFlag, AeroElastFlag, AeroElastModel, ConstantCirculationFlag
      real(mytype) :: toweroffset,tower_drag,tower_lift,tower_strouhal, uref, tsr, ShenC1, ShenC2, GammaCirc
      real(mytype) :: BladeInertia, GeneratorInertia, GBRatio, GBEfficiency, RatedGenSpeed
      real(mytype) :: RatedLimitGenTorque, CutInGenSpeed, Region2StartGenSpeed, Region2EndGenSpeed,Kgen
      real(mytype) :: RatedPower, MaximumTorque
      real(mytype) :: yaw_angle, shaft_tilt_angle, blade_cone_angle
      NAMELIST/TurbineSpecs/name,origin,numblades,blade_geom,numfoil,afname,towerFlag,towerOffset, &
         tower_geom,tower_drag,tower_lift,tower_strouhal,TypeFlag, OperFlag, tsr, uref,RotFlag, AddedMassFlag, &
         RandomWalkForcingFlag, DynStallFlag,dynstall_param_file,EndEffectsFlag,TipCorr, RootCorr,ShenC1, ShenC2, &
         ConstantCirculationFlag, GammaCirc, &
         yaw_angle, shaft_tilt_angle, blade_cone_angle,AeroElastFlag, AeroElastModel, AeroElastInputFile, AeroElastSolverFile, &
         BladeInertia, GeneratorInertia, GBRatio, GBEfficiency, RatedGenSpeed, RatedLimitGenTorque, CutInGenSpeed, &
         Region2StartGenSpeed,Region2EndGenSpeed,Kgen,RatedPower,MaximumTorque,list_controller_file

      if (nrank==0) then
         write(*,*) 'Loading the turbine options ...'
      endif

      ! Allocate arrays
      allocate(Turbine(Ntur))

      !==========================================
      ! Get turbine options and initialize them
      !==========================================
      do i=1,Ntur
         ! Set some default parameters
         numfoil=0
         towerFlag=0
         tower_lift=zpthree
         tower_drag=one
         tower_strouhal=zptwoone
         AddedMassFlag=0
         RandomWalkForcingFlag=0
         DynStallFlag=0
         EndEffectsFlag=0
         ConstantCirculationFlag=0
         TipCorr=0
         RootCorr=0
         ShenC1=0.125_mytype
         ShenC2=twentyone
         yaw_angle=zero
         shaft_tilt_angle=zero
         blade_cone_angle=zero

         open(100,File=turbines_path(i))
         read(100,nml=TurbineSpecs)
         Turbine(i)%name=name
         Turbine(i)%ID=i
         Turbine(i)%origin=origin
         Turbine(i)%NBlades=numblades
         Turbine(i)%blade_geom_file=blade_geom

         ! Allocate Blades
         allocate(Turbine(i)%Blade(Turbine(i)%NBlades))

         ! Count how many Airfoil Sections are available
         nfoils = numfoil

         if (nfoils==0) then
            write(*,*) 'You need to provide at least one static_foils_data entry for the computation of the blade forces'
            stop
         else
            if (nrank==0) then
               write(*,*) 'Loaded a total of : ', nfoils, 'airfoil data sets'
               do k=1,nfoils
                  write(*,*) 'Airfoil :', k, 'is ', afname(k)
               enddo
            endif
         endif

         ! Assign variables on the blade level
         do j=1,Turbine(i)%NBlades
            ! Assign the blade inertia
            Turbine(i)%Blade(j)%Inertia=BladeInertia

            ! Allocate the memory of the Airfoils
            Turbine(i)%Blade(j)%NAirfoilData=nfoils
            allocate(Turbine(i)%Blade(j)%AirfoilData(nfoils))

            ! Read and Store Airfoils
            do k=1,Turbine(i)%Blade(j)%NAirfoilData
               Turbine(i)%Blade(j)%AirfoilData(k)%afname=afname(k)
               call airfoil_init_data(Turbine(i)%Blade(j)%AirfoilData(k))
            enddo
         enddo

         ! Tower
         if (towerFlag==1) then
            Turbine(i)%Has_Tower=.true.
            Turbine(i)%TowerOffset=toweroffset
            Turbine(i)%Tower%geom_file=tower_geom
            Turbine(i)%TowerDrag=tower_drag
            Turbine(i)%TowerLift=tower_lift
            Turbine(i)%TowerStrouhal=tower_strouhal
         endif

         ! Check the type of turbine (choose between Horizontal and Vertical Axis turbines)
         if (TypeFlag==1) then
            Turbine(i)%Type='Horizontal_Axis'
            Turbine(i)%RotN=[one,zero,zero]
            Turbine(i)%shaft_tilt_angle=shaft_tilt_angle
            Turbine(i)%yaw_angle=yaw_angle
            Turbine(i)%blade_cone_angle=blade_cone_angle
         else if (TypeFlag==2) then
            write(*,*) 'Not ready yet'
            stop
         else
            write(*,*) 'You should not be here'
            stop
         endif

         ! Get Operation Options
         if (OperFlag==1) then
            Turbine(i)%Is_constant_rotation_operated=.true.
            Turbine(i)%Uref=uref
            Turbine(i)%TSR=tsr
         else if (OperFlag==2) then
            Turbine(i)%Is_NRELController=.true.
            ! Assign Uref and TSR (this applies only to the first time step)
            Turbine(i)%Uref=zpfive*(u1+u2)
            Turbine(i)%TSR=tsr
            ! Assign the blade inertia
            do j=1,Turbine(i)%NBlades
               Turbine(i)%Blade(j)%Inertia=BladeInertia
            enddo
            ! Initialize Contoller
            call init_controller(Turbine(i)%Controller,GeneratorInertia,GBRatio,GBEfficiency,&
               RatedGenSpeed,RatedLimitGenTorque,CutInGenSpeed,&
               Region2StartGenSpeed,Region2EndGenSpeed,Kgen,RatedPower,&
               MaximumTorque)
            Turbine(i)%Controller%IStatus=0
         else if (OperFlag==3) then
            Turbine(i)%Is_ListController= .true.
            ! Assign Uref and TSR (this applies only to the first time step)
            Turbine(i)%Uref=uref
            Turbine(i)%TSR=tsr
            ! Read the list_controller file
            call read_list_controller_file(list_controller_file,turbine(i))
         else if (OperFlag==4) then
            Turbine(i)%Is_upstreamvel_controlled=.true.
            ! Assign Uref and TSR (this applies only to the first time step)
            Turbine(i)%Uref=uref
            Turbine(i)%TSR=tsr
         else
            write(*,*) 'You should not be here. Select other turbine operation option'
            stop
         endif

         if (RotFlag==1) then
            Turbine(i)%IsClockwise=.true.
         else if (RotFlag==2) then
            Turbine(i)%IsCounterClockwise=.true.
            Turbine(i)%RotN=-Turbine(i)%RotN
         else
            write(*,*) 'You should not be here. The options are clockwise and counterclockwise'
            stop
         endif

         ! Get Unsteady Effects Modelling Options
         if (RandomWalkForcingFlag==1) then
            do j=1,Turbine(i)%NBlades
               Turbine(i)%Blade(j)%do_random_walk_forcing=.true.
            enddo
         endif

         if (AddedMassFlag==1) then
            do j=1,Turbine(i)%NBlades
               Turbine(i)%Blade(j)%do_added_mass=.true.
            enddo
         endif

         if (DynStallFlag>0) then
            if (DynStallFlag==1) then ! Sheng et al. model
               do j=1,Turbine(i)%NBlades
                  Turbine(i)%Blade(j)%do_Sheng_stall=.true.
                  Turbine(i)%Blade(j)%DynStallFile=dynstall_param_file
               enddo
            endif
            if (DynStallFlag==2) then ! Legacy LB model
               do j=1,Turbine(i)%NBlades
                  Turbine(i)%Blade(j)%do_LB_stall=.true.
                  Turbine(i)%Blade(j)%DynStallFile=dynstall_param_file
               enddo
            endif
         endif

         if (ConstantCirculationFlag==1) then
            do j=1,Turbine(i)%NBlades
               Turbine(i)%Blade(j)%Is_constant_circulation=.true.
               Turbine(i)%Blade(j)%GammaCirc=GammaCirc
            enddo
         endif

         if (EndEffectsFlag>0) then
            Turbine(i)%Has_BladeEndEffectModelling=.true.
            if (EndEffectsFlag==1) then
               Turbine(i)%EndEffectModel_is_Glauert=.true.
               if (TipCorr==1) Turbine(i)%do_tip_correction=.true.
               if (RootCorr==1) Turbine(i)%do_root_correction=.true.
            else if (EndEffectsFlag==2) then
               Turbine(i)%EndEffectModel_is_Shen=.true.
               Turbine(i)%ShenCoeff_c1=ShenC1
               Turbine(i)%ShenCoeff_c2=ShenC2
               if (TipCorr==1) Turbine(i)%do_tip_correction=.true.
               if (RootCorr==1) Turbine(i)%do_root_correction=.true.
            endif
         endif

         if (AeroElastFlag==1) then
            Turbine(i)%do_aeroelasticity=.true.
         endif

         ! close turbine file
         close(100)
      enddo

    end subroutine get_turbine_options

    !*******************************************************************************
    !
    subroutine get_actuatorline_options(actuatorline_path)
    !
    !*******************************************************************************

      use constants

      implicit none
      integer :: i,k
      character(len=80), dimension(100), intent(in) :: actuatorline_path
      character(len=100) :: name, actuatorline_geom, dynstall_param_file
      character(len=100*MaxNAirfoils) :: afname
      real(mytype), dimension(3) :: origin
      real(mytype) :: PitchStartTime, PitchEndTime, PitchAngleInit, PitchAmp, AngularPitchFreq
      integer :: numfoil, AddedMassFlag, DynStallFlag, EndEffectsFlag, RandomWalkForcingFlag, PitchControlFlag
      NAMELIST/ActuatorLineSpecs/name, origin, actuatorline_geom, numfoil, afname, AddedMassFlag, &
         DynStallFlag, PitchControlFlag, PitchStartTime, &
         PitchEndTime, PitchAngleInit, PitchAmp, AngularPitchFreq, dynstall_param_file

      if (nrank==0) then
         write(*,*) 'Loading the actuator line options ...'
      endif

       ! Allocate arrays
      allocate(actuatorline(Nal))

      !==========================================
      ! Get actuator line options and initialize them
      !==========================================
      do i=1,Nal
         open(200,File=actuatorline_path(i))
         read(200,nml=ActuatorLineSpecs)
         close(200)
         Actuatorline(i)%name=name
         Actuatorline(i)%COR=origin
         Actuatorline(i)%geom_file=actuatorline_geom

         ! Count how many Airfoil Sections are available
         Actuatorline(i)%NAirfoilData=numfoil

         ! Allocate the memory of the Airfoils
         allocate(Actuatorline(i)%AirfoilData(Actuatorline(i)%NAirfoilData))

         ! Read and Store Airfoils
         do k=1,Actuatorline(i)%NAirfoilData
            Actuatorline(i)%AirfoilData(k)%afname=afname
            call airfoil_init_data(Actuatorline(i)%AirfoilData(k))
         enddo

         ! Get Dynamic Loads Modelling Options
         if (AddedMassFlag==1) then
            Actuatorline%do_added_mass=.true.
         endif

         if (DynStallFlag>0) then
            if (DynStallFlag==1) then
               Actuatorline%do_Sheng_stall=.true.
               Actuatorline%DynStallFile=dynstall_param_file
            endif
            if (DynstallFlag==2) then
               Actuatorline%do_LB_stall=.true.
               Actuatorline%DynStallFile=dynstall_param_file
            endif
         endif

         ! Get Pitching Opions
         if (PitchControlFlag==1) then ! Sinusoidal
            Actuatorline%pitch_control=.true.
            Actuatorline(i)%pitch_start_time=PitchStartTime
            Actuatorline(i)%pitch_end_time=PitchEndTime
            Actuatorline(i)%pitch_angle_init=PitchAngleInit
            Actuatorline(i)%pitchAmp=PitchAmp
            Actuatorline(i)%angular_pitch_freq=AngularPitchFreq
         else if (PitchControlFlag==2) then ! Ramping up
            Actuatorline%pitch_control=.true.
            write(*,*) 'Ramping up not implemented yet'
            stop
         endif

      enddo

    end subroutine get_actuatorline_options

    !*******************************************************************************
    !
    subroutine actuator_line_model_update(current_time,dt)
    !
    !*******************************************************************************

      use decomp_2d, only: mytype
      use param, only: zpfive, two, sixty, onehundredeighty
      use dbg_schemes, only: sin_prec

      implicit none
      real(mytype), intent(inout) :: current_time, dt
      integer :: i,j,k, Nstation
      real(mytype) :: theta, pitch_angle, deltapitch, pitch_angle_old
      real(mytype) :: WSRotorAve, Omega

      !==========================================
      ! This routine updates the location of the actuator lines
      !==========================================

      ctime=current_time
      DeltaT=dt

      if (Ntur>0) then
         do i=1,Ntur
            if (Turbine(i)%Is_constant_rotation_operated) then
               ! Computes the rigid-body velocity
               theta=Turbine(i)%angularVel*DeltaT
               Turbine(i)%AzimAngle=Turbine(i)%AzimAngle+theta
               call rotate_turbine(Turbine(i),Turbine(i)%RotN,theta)

               ! Computes the displacement on the turbine
               !if (turbine(i)%do_aeroelasticity) then
               !    call actuator_line_beam_solve(turbine(i)%beam,DeltaT)
               !endif

               ! Computes the new velocity
               call Compute_Turbine_RotVel(Turbine(i))

            else if (Turbine(i)%Is_NRELController) then
               ! First do control
               call compute_rotor_upstream_velocity(Turbine(i))
               call operate_controller(Turbine(i)%Controller,ctime,Turbine(i)%NBlades,Turbine(i)%angularVel)
               Turbine(i)%deltaOmega=(Turbine(i)%Torque-Turbine(i)%Controller%GearBoxRatio*Turbine(i)%Controller%GenTrq)/(Turbine(i)%IRotor+Turbine(i)%Controller%GearBoxRatio**2.*Turbine(i)%Controller%IGenerator)*DeltaT

               ! Then calculate the angular velocity and compute DeltaTheta and AzimAngle
               Turbine(i)%angularVel=Turbine(i)%angularVel+Turbine(i)%deltaOmega
               theta=Turbine(i)%angularVel*DeltaT
               Turbine(i)%AzimAngle=Turbine(i)%AzimAngle+theta

               ! Then do pitch control
               do j=1,Turbine(i)%NBlades
                  if (Turbine(i)%IsClockwise) then
                     Turbine(i)%cbp=-Turbine(i)%Controller%PitCom(j)
                  else
                     stop
                  endif
                  deltapitch=Turbine(i)%cbp-Turbine(i)%cbp_old
                  if (nrank==0.and.mod(itime,ilist)==0) write(*,*) 'Doing Pitch control', -deltapitch*onehundredeighty/pi
                  call pitch_actuator_line(Turbine(i)%Blade(j),deltapitch)
               enddo
               Turbine(i)%cbp_old=Turbine(i)%cbp

               call rotate_turbine(Turbine(i),Turbine(i)%RotN,theta)
               call Compute_Turbine_RotVel(Turbine(i))

               ! After you do both variable speed and pitch control update the status of the controller
               Turbine(i)%Controller%IStatus=Turbine(i)%Controller%IStatus+1

            else if (Turbine(i)%Is_ListController) then
               ! Compute the rotor averaged wind speed
               call compute_rotor_upstream_velocity(Turbine(i))
               WSRotorAve=Turbine(i)%Uref

               ! Compute Omega and pitch by interpolating from the list
               call from_list_controller(Omega,pitch_angle,turbine(i),WSRotorAve)
               Turbine(i)%angularVel=Omega/(sixty/(two*pi)) ! Convert to rad/s
               theta=Turbine(i)%angularVel*DeltaT
               Turbine(i)%AzimAngle=Turbine(i)%AzimAngle+theta
               call rotate_turbine(Turbine(i),Turbine(i)%RotN,theta)
               call Compute_Turbine_RotVel(Turbine(i))

               ! Do pitch control
               if (Turbine(i)%IsClockwise) then
                  Turbine(i)%cbp = pitch_angle*pi/onehundredeighty ! Convert to rad
               else
                  Turbine(i)%cbp = pitch_angle*pi/onehundredeighty ! Convert to rad
               endif
               deltapitch = Turbine(i)%cbp - Turbine(i)%cbp_old
               do j=1,Turbine(i)%NBlades
                  call pitch_actuator_line(Turbine(i)%Blade(j),deltapitch)
               enddo
               Turbine(i)%cbp_old=Turbine(i)%cbp

            else if (Turbine(i)%Is_upstreamvel_controlled) then
               ! Compute the rotor averaged wind speed
               call compute_rotor_upstream_velocity(Turbine(i))
               Turbine(i)%angularVel=Turbine(i)%Uref*Turbine(i)%TSR/Turbine(i)%Rmax
               theta=Turbine(i)%angularVel*DeltaT
               Turbine(i)%AzimAngle=Turbine(i)%AzimAngle+theta
               call rotate_turbine(Turbine(i),Turbine(i)%RotN,theta)
               call Compute_Turbine_RotVel(Turbine(i))
            endif
         enddo
         ! DO FARM_LEVEL CONTROL
      endif

      if (Nal>0) then
         do i=1,Nal
            if (ActuatorLine(i)%pitch_control.and.ctime > ActuatorLine(i)%pitch_start_time.and.ctime < ActuatorLine(i)%pitch_end_time) then
               ! Do harmonic pitch control for all elements (stations) of the actuator line
               Nstation=ActuatorLine(i)%NElem+1
               do j=1,Nstation
                  ActuatorLine(i)%pitch(j)=actuatorline(i)%pitchAmp*sin_prec(actuatorline(i)%angular_pitch_freq*(ctime-ActuatorLine(i)%pitch_start_time))
               enddo
               if (nrank==0.and.mod(itime,ilist)==0) then        
                  write(*,*) 'Harmonic pitch :'
                  write(*,*) 'Current pitch angle : ', sum(ActuatorLine(i)%pitch)/(ActuatorLine%Nelem+1)
               endif
               call pitch_actuator_line(actuatorline(i),ActuatorLine(i)%pitch(j))
            endif
         enddo
      endif

      return
      
    end subroutine actuator_line_model_update

    !*******************************************************************************
    !
    subroutine actuator_line_model_compute_forces
    !
    !*******************************************************************************

      use param, only: rho_air, itime, initstat

      implicit none
      integer :: i,j

      ! Get into each Turbine and Compute the Forces blade by blade and element by element
      if (Ntur>0) then
         do i=1,Ntur
            ! First compute the end effects on the turbine
            if (Turbine(i)%Has_BladeEndEffectModelling) then
               call Compute_Turbine_EndEffects(Turbine(i))
            endif

            ! Then compute the coefficients
            do j=1,Turbine(i)%Nblades
               call Compute_ActuatorLine_Forces(Turbine(i)%Blade(j),rho_air,visc,deltaT,ctime)
            enddo
            call Compute_performance(Turbine(i), rho_air)

            ! Tower
            if (Turbine(i)%has_tower) then
               call Compute_Tower_Forces(Turbine(i)%Tower,rho_air,visc,ctime,Turbine(i)%TowerLift,Turbine(i)%TowerDrag,Turbine(i)%TowerStrouhal)
            endif
         enddo

         ! Get statistics
         if (itime>=initstat) then
            call actuator_line_statistics()
         endif
      endif

      if (Nal>0) then
         do i=1,Nal
            call Compute_ActuatorLine_Forces(ActuatorLine(i),rho_air,visc,deltaT,ctime)
         enddo
      endif

      return

    end subroutine actuator_line_model_compute_forces

    !*******************************************************************************
    !
    subroutine actuator_line_statistics()
    !
    !*******************************************************************************

      implicit none
      integer :: itur

      do itur=1,Ntur
         turbine(itur)%P_ave=turbine(itur)%P_ave+Turbine(itur)%Power
         turbine(itur)%T_ave=turbine(itur)%T_ave+Turbine(itur)%Thrust
         turbine(itur)%Torque_ave=turbine(itur)%Torque_ave+Turbine(itur)%Torque
      enddo

      return

    end subroutine actuator_line_statistics

    !*******************************************************************************
    !
    subroutine actuator_line_model_write_restart()
    !
    !*******************************************************************************

      implicit none
      integer :: itur
      character(len=30) :: filename

      if (Ntur>0) then
         write(filename,"('restart',I7.7'.alm')") itime
         open(2021,file=filename)
         write(2021,*) 'Azimuthal angle, Angular velocity, Collective blade pitch'
         do itur=1,Ntur
            write(2021,*) turbine(itur)%AzimAngle, turbine(itur)%angularVel, turbine(itur)%cbp
         enddo
         close(2021)
      endif

      return 

    end subroutine actuator_line_model_write_restart
    
end module actuator_line_model
