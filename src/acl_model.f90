module actuator_line_model

   ! Use the actuator_line Modules
   use decomp_2d, only: mytype, nrank
   use actuator_line_model_utils
   use airfoils
   use actuator_line_element
   use actuator_line_turbine
   use actuator_line_controller
   use actuator_line_write_output
   use dynstall

   implicit none

   type(ActuatorLineType), allocatable, save :: Actuatorline(:)
   type(TurbineType), allocatable, save :: Turbine(:) ! Turbine
   integer,save :: Ntur, Nal ! Number of the turbines
   real(mytype),save :: DeltaT, Visc, ctime
   logical,save :: actuator_line_model_writeFlag=.true.


   public :: actuator_line_model_init, actuator_line_model_write_output, &
      actuator_line_model_update, actuator_line_model_compute_forces

contains

   subroutine actuator_line_model_init(Nturbines,Nactuatorlines,turbines_file,actuatorlines_file,dt)

      use param, only: ialmrestart, filealmrestart

      implicit none
      integer :: Nturbines, Nactuatorlines
      character(len=80),dimension(100),intent(in) :: turbines_file, actuatorlines_file
      real(mytype), intent(in) :: dt
      integer :: itur,ial
      character(1000) :: ReadLine

      if (nrank==0) then
         write(6,*) '====================================================='
         write(6,*) 'Initializing the Actuator Line Model'
         write(6,*) 'Developed by G. Deskos 2015-2017'
         write(6,*) 'gd1414@ic.ac.uk'
         write(6,*) '====================================================='
      end if

      !### Specify Turbines
      Ntur=Nturbines
      if (nrank==0) then
         write(6,*) 'Number of turbines : ', Ntur
      endif
      call get_turbine_options(turbines_file)
      if (Ntur>0) then
         if (ialmrestart==1) then
            ! Read the checkpoint information and rotate actuator lines accordingly
            open(17,file=filealmrestart)
            ! Read the azimuthal angle
            do itur=1,Ntur
               read(17,'(A)') ReadLine ! Turb ....
               read(ReadLine,*) turbine(itur)%AzimAngle, turbine(itur)%angularVel,turbine(itur)%cbp
            end do
            close(17)
         endif
         do itur=1,Ntur
            call set_turbine_geometry(Turbine(itur))
         end do
      endif

      !### Speficy Actuator Lines
      Nal=NActuatorlines
      if (nrank==0) then
         write(6,*) 'Number of actuator lines : ', Nal
      endif
      call get_actuatorline_options(actuatorlines_file)
      if(Nal>0) then
         do ial=1,Nal
            call set_actuatorline_geometry(actuatorline(ial))
         end do
      endif

      !##### Initialize Time
      ctime=0.0
      DeltaT=dt

   end subroutine actuator_line_model_init
   !############################################################################
   !############################################################################
   subroutine actuator_line_model_write_output(dump_no)

      implicit none
      integer,intent(in) :: dump_no
      integer :: itur,ial,iblade
      character(len=100) :: dir

      !call system('mkdir -p ALM/'//adjustl(trim(dirname(dump_no))))

      !dir='ALM/'//adjustl(trim(dirname(dump_no)))

      dir=adjustl(trim(dirname(dump_no)))

      if (Ntur>0) then
         do itur=1,Ntur
            call actuator_line_turbine_write_output(turbine(itur),dir)
            do iblade=1,turbine(itur)%NBlades
               call actuator_line_element_write_output(turbine(itur)%Blade(iblade),dir)

               if (turbine(itur)%Blade(iblade)%do_Sheng_stall) then
                  call dynamic_stall_write_output(turbine(itur)%Blade(iblade),dir)
               end if

            end do
            if(turbine(itur)%Has_Tower) then
               call actuator_line_element_write_output(turbine(itur)%tower,dir)
            endif
         end do
      endif

      if (Nal>0) then
         do ial=1,Nal
            call actuator_line_element_write_output(actuatorline(ial),dir)
            if (actuatorline(ial)%do_Sheng_stall) then
               call dynamic_stall_write_output(actuatorline(ial),dir)
            end if
         end do
      endif

   end subroutine actuator_line_model_write_output

   subroutine get_turbine_options(turbines_path)

      use decomp_2d, only : mytype
      use param, only: u1,u2
      use param, only: zero, zpthree, zptwoone, zpfive, one, twentyone
      use constants
      implicit none

      character(len=80),dimension(100),intent(in) :: turbines_path
      integer :: i,j,k
      ! integer, parameter :: MaxReadLine = 1000
      integer :: nfoils
      !-------------------------------------
      ! Dummy variables
      !-------------------------------------
      character(len=100) :: name, blade_geom, tower_geom, dynstall_param_file, AeroElastInputFile, AeroElastSolverFile, list_controller_file
      character(len=100),dimension(MaxNAirfoils) :: afname
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
         write(6,*) 'Loading the turbine options ...'
      endif

      ! Allocate Turbines Arrays
      Allocate(Turbine(Ntur))
      ! ==========================================
      ! Get Turbines' options and INITIALIZE THEM
      ! ==========================================
      do i=1, Ntur
         ! Set some default parameters
         !++++++++++++++++++++++++++++++++
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
         !+++++++++++++++++++++++++++++++++
         open(100,File=turbines_path(i))
         read(100,nml=TurbineSpecs)
         Turbine(i)%name=name
         Turbine(i)%ID=i
         Turbine(i)%origin=origin
         Turbine(i)%NBlades=numblades
         Turbine(i)%blade_geom_file=blade_geom

         ! Allocate Blades
         Allocate(Turbine(i)%Blade(Turbine(i)%NBlades))

         ! Count how many Airfoil Sections are available
         nfoils = numfoil

         if(nfoils==0) then
            write(6,*) "You need to provide at least one static_foils_data entry for the computation of the blade forces"
            stop
         else
            if (nrank==0) then
               write(6,*) "Loaded a total of : ", nfoils, "airfoil data sets"
               do k=1,nfoils
                  write(6,*) "Airfoil :", k, "is ", afname(k)
               enddo
            endif
         end if

         ! Assign variables on the blade level
         do j=1,Turbine(i)%NBlades

            ! Assign the blade inertia
            Turbine(i)%Blade(j)%Inertia=BladeInertia

            ! Allocate the memory of the Airfoils
            Turbine(i)%Blade(j)%NAirfoilData=nfoils
            Allocate(Turbine(i)%Blade(j)%AirfoilData(nfoils))

            do k=1, Turbine(i)%Blade(j)%NAirfoilData
               Turbine(i)%Blade(j)%AirfoilData(k)%afname=afname(k)
               ! Read and Store Airfoils
               call airfoil_init_data(Turbine(i)%Blade(j)%AirfoilData(k))
            end do
         end do

         ! ## Tower ##
         if (towerFlag==1) then
            Turbine(i)%Has_Tower=.true.
            Turbine(i)%TowerOffset=toweroffset
            Turbine(i)%Tower%geom_file=tower_geom
            Turbine(i)%TowerDrag=tower_drag
            Turbine(i)%TowerLift=tower_lift
            Turbine(i)%TowerStrouhal=tower_strouhal
         endif

         !#############2  Get turbine_specs #################
         ! Check the typ of Turbine (choose between Horizontal and Vertical Axis turbines)
         if(TypeFlag==1) then
            Turbine(i)%Type='Horizontal_Axis'
            Turbine(i)%RotN=[one,zero,zero]
            Turbine(i)%shaft_tilt_angle=shaft_tilt_angle
            Turbine(i)%yaw_angle=yaw_angle
            Turbine(i)%blade_cone_angle=blade_cone_angle
         elseif(TypeFlag==2) then
            !        call get_option(trim(turbine_path(i))//"/type/Vertical_Axis/axis_of_rotation",Turbine(i)%RotN)
            !        call get_option(trim(turbine_path(i))//"/type/Vertical_Axis/distance_from_axis",Turbine(i)%dist_from_axis)
            write(6,*) 'Not ready yet'
            stop
         else
            write(6,*) "You should not be here"
            stop
         end if

         !##############3 Get Operation Options ######################
         if (OperFlag==1) then
            Turbine(i)%Is_constant_rotation_operated= .true.
            Turbine(i)%Uref=uref
            Turbine(i)%TSR=tsr
         else if(OperFlag==2) then
            Turbine(i)%Is_NRELController = .true.
            Turbine(i)%Uref=zpfive*(u1+u2)
            Turbine(i)%TSR=tsr
            ! Assign the Uref and TSR to compute the optimum tip-speed ratio (This applies only to the first time step)
            ! Controller Variables
            !--------------------------
            do j=1,Turbine(i)%NBlades
               ! Assign the blade inertia
               Turbine(i)%Blade(j)%Inertia=BladeInertia
            enddo
            ! Initialize Contoller
            call init_controller(Turbine(i)%Controller,GeneratorInertia,GBRatio,GBEfficiency,&
               RatedGenSpeed,RatedLimitGenTorque,CutInGenSpeed,&
               Region2StartGenSpeed,Region2EndGenSpeed,Kgen,RatedPower,&
               MaximumTorque)

            Turbine(i)%Controller%IStatus=0
         else if(OperFlag==3) then
            ! Assign the Uref and TSR to compute the optimum tip-speed ratio (This applies only to the first time step)
            Turbine(i)%Uref=uref
            Turbine(i)%TSR=tsr
            Turbine(i)%Is_ListController= .true.
            !> Read the list_controller file
            call read_list_controller_file(list_controller_file,turbine(i))
         else if(OperFlag==4) then
            Turbine(i)%TSR=tsr
            Turbine(i)%Is_upstreamvel_controlled=.true.
            Turbine(i)%Uref=uref
         else
            write(*,*) "Only constant_rotation (1) and control_based (2) is used"
            stop
         endif

         if(RotFlag==1) then
            Turbine(i)%IsClockwise=.true.
         elseif(RotFlag==2) then
            Turbine(i)%IsCounterClockwise=.true.
            Turbine(i)%RotN=-Turbine(i)%RotN
         else
            write(6,*) "You should not be here. The options are clockwise and counterclockwise"
            stop
         endif

         !##################4 Get Unsteady Effects Modelling Options ##################
         if(RandomWalkForcingFlag==1) then
            do j=1,Turbine(i)%NBlades
               Turbine(i)%Blade(j)%do_random_walk_forcing=.true.
            end do
         endif

         if(AddedMassFlag==1) then
            do j=1,Turbine(i)%NBlades
               Turbine(i)%Blade(j)%do_added_mass=.true.
            end do
         endif

         if(DynStallFlag>0) then
            if(DynStallFlag==1) then ! Do Sheng et al. modelling
               do j=1,Turbine(i)%NBlades
                  Turbine(i)%Blade(j)%do_Sheng_stall=.true.
                  Turbine(i)%Blade(j)%DynStallFile=dynstall_param_file
               end do
            endif
            if(DynStallFlag==2) then ! Do the legacy LB model
               do j=1,Turbine(i)%NBlades
                  Turbine(i)%Blade(j)%do_LB_stall=.true.
                  Turbine(i)%Blade(j)%DynStallFile=dynstall_param_file
               end do
            endif
         endif

         if (ConstantCirculationFlag==1) then
            do j=1,Turbine(i)%NBlades
               Turbine(i)%Blade(j)%Is_constant_circulation=.true.
               Turbine(i)%Blade(j)%GammaCirc=GammaCirc
            end do
         endif

         if (EndEffectsFlag>0) then
            Turbine(i)%Has_BladeEndEffectModelling=.true.
            if(EndEffectsFlag==1) then
               Turbine(i)%EndEffectModel_is_Glauret=.true.
               if(TipCorr==1) Turbine(i)%do_tip_correction=.true.
               if (RootCorr==1) Turbine(i)%do_root_correction=.true.
            else if(EndEffectsFlag==2) then
               Turbine(i)%EndEffectModel_is_Shen=.true.
               Turbine(i)%ShenCoeff_c1=ShenC1
               Turbine(i)%ShenCoeff_c2=ShenC2
               if(TipCorr==1) Turbine(i)%do_tip_correction=.true.
               if (RootCorr==1) Turbine(i)%do_root_correction=.true.
            endif
         endif

         if (AeroElastFlag==1) then
            Turbine(i)%do_aeroelasticity=.true.
         endif


      end do

   end subroutine get_turbine_options

   subroutine get_actuatorline_options(actuatorline_path)

      use constants
      implicit none

      integer :: i,k
      ! integer, parameter :: MaxReadLine = 1000
      character(len=80), dimension(100),intent(in) :: actuatorline_path
      !-------------------------------------
      ! Dummy variables
      !-------------------------------------
      character(len=100) :: name, actuatorline_geom, dynstall_param_file
      character(len=100*MaxNAirfoils) :: afname
      real(mytype), dimension(3) :: origin
      real(mytype) :: PitchStartTime, PitchEndTime, PitchAngleInit, PitchAmp, AngularPitchFreq
      integer :: numfoil, AddedMassFlag, DynStallFlag, EndEffectsFlag, RandomWalkForcingFlag, PitchControlFlag
      NAMELIST/ActuatorLineSpecs/name,origin, actuatorline_geom,numfoil,afname, AddedMassFlag, &
         DynStallFlag,EndEffectsFlag, PitchControlFlag,PitchStartTime, &
         PitchEndTime, PitchAngleInit, PitchAmp, AngularPitchFreq, dynstall_param_file

      if (nrank==0) then
         write(6,*) 'Loading the actuator line options ...'
      endif

       ! Allocate Turbines Arrays
      allocate(actuatorline(Nal))

      ! ==========================================
      ! Get Actuator lines' options and INITIALIZE THEM
      ! ==========================================
      do i=1, Nal
         open(200,File=actuatorline_path(i))
         read(200,nml=ActuatorLineSpecs)
         close(200)
         Actuatorline(i)%name=name
         Actuatorline(i)%COR=origin
         Actuatorline(i)%geom_file=actuatorline_geom

         ! Count how many Airfoil Sections are available
         Actuatorline(i)%NAirfoilData=numfoil
         ! Allocate the memory of the Airfoils
         Allocate(Actuatorline(i)%AirfoilData(Actuatorline(i)%NAirfoilData))

         do k=1, Actuatorline(i)%NAirfoilData

            Actuatorline(i)%AirfoilData(k)%afname=afname

            ! Read and Store Airfoils
            call airfoil_init_data(Actuatorline(i)%AirfoilData(k))
         end do

         !##################4 Get Dynamic Loads Modelling Options ##################
         if(AddedMassFlag==1) then
            Actuatorline%do_added_mass=.true.
         endif

         if(DynStallFlag>0) then
            if(DynStallFlag==1) then
               Actuatorline%do_Sheng_stall=.true.
               Actuatorline%DynStallFile=dynstall_param_file
            endif
            if(DynstallFlag==2) then
               Actuatorline%do_LB_stall=.true.
               Actuatorline%DynStallFile=dynstall_param_file
            endif
         endif

         !    !##################4 Get Pitching Opions ##################
         if(PitchControlFlag==1) then ! Sinusoidal
            Actuatorline%pitch_control=.true.
            Actuatorline(i)%pitch_start_time=PitchStartTime
            Actuatorline(i)%pitch_end_time=PitchEndTime

            Actuatorline(i)%pitch_angle_init=PitchAngleInit
            Actuatorline(i)%pitchAmp=PitchAmp
            Actuatorline(i)%angular_pitch_freq=AngularPitchFreq
         elseif(PitchControlFlag==2) then ! Ramping up
            Actuatorline%pitch_control=.true.
            write(*,*) 'Ramping up not implemented yet'
            stop
         endif

      end do

   end subroutine get_actuatorline_options

   subroutine actuator_line_model_update(current_time,dt)

      use decomp_2d, only: mytype
      use param, only : zpfive, two, sixty, onehundredeighty
      use dbg_schemes, only: sin_prec

      implicit none
      real(mytype),intent(inout) :: current_time, dt
      integer :: i,j,k, Nstation
      real(mytype) :: theta, pitch_angle, deltapitch, pitch_angle_old
      real(mytype) :: WSRotorAve,Omega
      ! This routine updates the location of the actuator lines

      ctime=current_time
      DeltaT=dt

      if (Ntur>0) then
         do i=1,Ntur
            if(Turbine(i)%Is_constant_rotation_operated) then
               theta=Turbine(i)%angularVel*DeltaT
               Turbine(i)%AzimAngle=Turbine(i)%AzimAngle+theta
               ! Computes the rigid-body velocity
               call rotate_turbine(Turbine(i),Turbine(i)%RotN,theta)
               ! Computes the displacement on the turbine
               !if(turbine(i)%do_aeroelasticity) then
               !    call actuator_line_beam_solve(turbine(i)%beam,DeltaT)
               !endif
               ! Computes the new velocity
               call Compute_Turbine_RotVel(Turbine(i))
            else if(Turbine(i)%Is_NRELController) then
               !if(nrank==0) write(*,*) 'Entering the control-based operation for turbine', Turbine(i)%name
               ! First do control
               call compute_rotor_upstream_velocity(Turbine(i))
               call operate_controller(Turbine(i)%Controller,ctime,Turbine(i)%NBlades,Turbine(i)%angularVel)
               Turbine(i)%deltaOmega=(Turbine(i)%Torque-Turbine(i)%Controller%GearBoxRatio*Turbine(i)%Controller%GenTrq)/(Turbine(i)%IRotor+Turbine(i)%Controller%GearBoxRatio**2.*Turbine(i)%Controller%IGenerator)*DeltaT
               Turbine(i)%angularVel=Turbine(i)%angularVel+Turbine(i)%deltaOmega
               ! Then Calculate the angular velocity and compute the DeltaTheta  and AzimAngle
               theta=Turbine(i)%angularVel*DeltaT
               Turbine(i)%AzimAngle=Turbine(i)%AzimAngle+theta

               ! Then do picth control (if not zero)
               do j=1,Turbine(i)%NBlades
                  if (Turbine(i)%IsClockwise) then
                     Turbine(i)%cbp=-Turbine(i)%Controller%PitCom(j)
                  else
                     stop
                  endif
                  deltapitch=Turbine(i)%cbp-Turbine(i)%cbp_old
                  if(nrank==0) write(*,*) 'Doing Pitch control', -deltapitch*onehundredeighty/pi
                  call pitch_actuator_line(Turbine(i)%Blade(j),deltapitch)
               enddo
               Turbine(i)%cbp_old=Turbine(i)%cbp

               call rotate_turbine(Turbine(i),Turbine(i)%RotN,theta)
               call Compute_Turbine_RotVel(Turbine(i))

               ! After you do both variable speed and pitch control update the status of the controller
               Turbine(i)%Controller%IStatus=Turbine(i)%Controller%IStatus+1

            else if(Turbine(i)%Is_ListController) then

               if(nrank==0) write(*,*) 'Entering the List-controlled operation for the turbine', Turbine(i)%name
               !> Compute the rotor averaged wind speed
               call compute_rotor_upstream_velocity(Turbine(i))
               Turbine(i)%Uref=WSRotorAve
               !> Compute Omega and pitch by interpolating from the list
               call from_list_controller(Omega,pitch_angle,turbine(i),WSRotorAve)
               Turbine(i)%angularVel=Omega/(sixty/(two*pi))  !Translate in Rad/sec
               theta=Turbine(i)%angularVel*DeltaT
               Turbine(i)%AzimAngle=Turbine(i)%AzimAngle+theta
               call rotate_turbine(Turbine(i),Turbine(i)%RotN,theta)
               call Compute_Turbine_RotVel(Turbine(i))

               !> Do pitch control
               ! Then do picth control (if not zero)
               if (Turbine(i)%IsClockwise) then
                  deltapitch=pitch_angle
               else
                  deltapitch=pitch_angle
               endif

               do j=1,Turbine(i)%NBlades
                  call pitch_actuator_line(Turbine(i)%Blade(j),deltapitch)
               enddo
            else if (Turbine(i)%Is_upstreamvel_controlled) then
               call compute_rotor_upstream_velocity(Turbine(i))
               !Turbine(i)%angularVel=Turbine(i)%Uref*Turbine(i)%TSR/Turbine(i)%Rmax
               Turbine(i)%angularVel=sqrt(abs(Turbine(i)%Torque)/(zpfive*0.432_mytype*Turbine(i)%A*Turbine(i)%Rmax**3./Turbine(i)%TSR**3.))
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
            if(ActuatorLine(i)%pitch_control.and.ctime > ActuatorLine(i)%pitch_start_time.and.ctime < ActuatorLine(i)%pitch_end_time) then
               !> Do harmonic pitch control for all elements of the actuator line
               Nstation=ActuatorLine(i)%NElem+1
               do j=1,Nstation
                  ActuatorLine(i)%pitch(j)=actuatorline(i)%pitchAmp*sin_prec(actuatorline(i)%angular_pitch_freq*(ctime-ActuatorLine(i)%pitch_start_time))
               end do
               if(nrank==0) then
                  write(*,*) '-----------------------'
                  write(*,*) ' Harmonic pitch :'
                  write(*,*) '-----------------------'
                  write(*,*) 'Current pitch angle : ', sum(ActuatorLine(i)%pitch)/(ActuatorLine%Nelem+1)
               endif

               call pitch_actuator_line(actuatorline(i),ActuatorLine(i)%pitch(j))
            endif
         enddo
      endif

      return
   end subroutine actuator_line_model_update

   subroutine actuator_line_model_compute_forces

      use param, only: rho_air
      implicit none

      integer :: i,j
      ! Zero the Source Term at each time step

      !write(*,*) 'Entering the actuator_line_model_compute_forces'

      if (Ntur>0) then

         ! Get into each Turbine and Compute the Forces blade by blade and element by element
         do i=1,Ntur

            ! First compute the end effects on the turbine and
            if (Turbine(i)%Has_BladeEndEffectModelling) then
               call Compute_Turbine_EndEffects(Turbine(i))
            endif

            ! Then compute the coefficients
            do j=1,Turbine(i)%Nblades
               call Compute_ActuatorLine_Forces(Turbine(i)%Blade(j),rho_air,visc,deltaT,ctime)
            end do
            call Compute_performance(Turbine(i), rho_air)

            ! Tower
            if(Turbine(i)%has_tower) then
               call Compute_Tower_Forces(Turbine(i)%Tower,rho_air,visc,ctime,Turbine(i)%TowerLift,Turbine(i)%TowerDrag,Turbine(i)%TowerStrouhal)
            endif

         end do
      end if

      if (Nal>0) then
         do i=1,Nal
            call Compute_ActuatorLine_Forces(ActuatorLine(i),rho_air,visc,deltaT,ctime)
         end do
      end if

      !write(*,*) 'Exiting actuator_line_model_compute_forces'
      return

   end subroutine actuator_line_model_compute_forces

   subroutine actuator_line_statistics()

      implicit none
      integer :: itur

      do itur=1,Ntur
         turbine(itur)%CP_ave    =turbine(itur)%CP_ave+Turbine(itur)%CP
         turbine(itur)%CT_ave    =turbine(itur)%CT_ave+Turbine(itur)%CT
         turbine(itur)%Torque_ave=turbine(itur)%CP_ave+Turbine(itur)%Torque
         if(nrank==0) then
            print *, 'Writing statistics for alm'
            print *, 'Turbine name , CT_ave, CP_ave, Torque_ave'
            print *, turbine(itur)%name, turbine(itur)%CT_ave, turbine(itur)%CP_ave, turbine(itur)%Torque_ave
         endif

      enddo

      return

   end subroutine actuator_line_statistics

end module actuator_line_model
