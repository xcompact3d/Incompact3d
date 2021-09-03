module actuator_line_turbine
    
    use decomp_2d, only: mytype, nrank
    use actuator_line_model_utils
    use Airfoils
    use actuator_line_element
    !use actuator_line_beam_model
    use actuator_line_controller

    implicit none

type TurbineType 
    character(len=100) :: name
    character(len=100) :: blade_geom_file
    character(len=100) :: type
    integer :: ID
    integer :: NBlades
    real(mytype), dimension(3) :: RotN, origin ! Rotational vectors in the normal and perpendicular directions
    real(mytype) :: shaft_tilt_angle, blade_cone_angle, yaw_angle 
    real(mytype) :: Rmax ! Reference radius, velocity, viscosity
    real(mytype) :: IRotor ! Inertia of the Rotor
    real(mytype) :: A ! Rotor area
    real(mytype) :: Torque, angularVel,deltaOmega,TSR,Uref ! Torque and rotation for the shaft  
    real(mytype) :: Ux_upstream=0.0
    real(mytype) :: Uy_upstream=0.0
    real(mytype) :: Uz_upstream=0.0
    real(mytype) :: AzimAngle=0.0
    real(mytype) :: dist_from_axis=0.0
    real(mytype) :: cbp=0.0, cbp_old=0.0 ! Collective blade pitch
    integer :: No_rev=0.0
    logical :: Is_constant_rotation_operated = .false. ! For a constant rotational velocity (in Revolutions Per Minute)
    logical :: Is_NRELController = .false. ! Active control-based rotational velocity using the NREL controller
    logical :: Is_ListController = .false. ! Active control-based rotational velocity using the rotor averaged velocity and
                                           ! interpolation lists
    logical :: Is_upstreamvel_controlled=.false. 

    logical :: do_aeroelasticity=.true.    ! Flag
    !type(BeamType) :: beam         ! Elastic structural beam model for the blades
    

    integer :: ContEntries
    real(mytype), allocatable :: ContWindSpeed(:), ContOmega(:), ContPitch(:)

    type(ControllerType) :: controller     ! Contoller
    logical :: IsClockwise = .false.
    logical :: IsCounterClockwise = .false. 
    logical :: Has_Tower=.false.
    real(mytype) :: Towerheight, TowerOffset, TowerDrag, TowerLift, TowerStrouhal
    logical :: Has_BladeEndEffectModelling=.false.
    logical :: do_tip_correction=.false.
    logical :: do_root_correction=.false.
    logical :: EndEffectModel_is_Glauret=.false.
    logical :: EndEffectModel_is_Shen=.false.
    real(mytype) :: ShenCoeff_c1, ShenCoeff_c2
    type(ActuatorLineType), allocatable :: Blade(:)
    type(AirfoilType), allocatable :: AirfoilData(:)
    type(ActuatorLineType) :: tower
    real(mytype) :: CP  ! Power coefficient 
    real(mytype) :: CTR ! Torque coefficient 
    real(mytype) :: CFx ! Fx coefficient 
    real(mytype) :: CFy ! Fy coefficient 
    real(mytype) :: CFz ! Fz coefficient 
    real(mytype) :: CT  ! Thrust coefficient    
    real(mytype) :: Thrust, Power ! Absolute values for Thrust and Power

    ! Rotor Statistics
    real(mytype) :: CT_ave=0.0_mytype, CP_ave=0.0_mytype, Torque_ave=0.0_mytype

end type TurbineType
    

contains
    
    subroutine set_turbine_geometry(turbine)
    use param, only: ialmrestart
    implicit none
    type(TurbineType),intent(inout) :: turbine
    real(mytype), allocatable :: rR(:),ctoR(:),pitch(:),thick(:)
    real(mytype) :: SVec(3), theta
    integer :: Nstations, iblade, Istation,ielem

    call read_actuatorline_geometry(turbine%blade_geom_file,turbine%Rmax,SVec,rR,ctoR,pitch,thick,Nstations)
    ! Make sure that the spanwise is [0 0 1]
    Svec = [0.0d0,0.0d0,1.0d0]
    ! Make sure that origin is [0,0,0] : we set everything to origin 0 and then translate the
    ! turbine to the actual origin(this is for simplicity)
    theta=2*pi/turbine%Nblades
    do iblade=1,turbine%Nblades
    call allocate_actuatorline(Turbine%blade(iblade),Nstations)
    turbine%blade(iblade)%name=trim(turbine%name)//'_blade_'//int2str(iblade)
    
    turbine%blade(iblade)%COR(1:3)=turbine%origin(1:3)
    turbine%blade(iblade)%L=turbine%Rmax
    turbine%blade(iblade)%NElem=Nstations-1 
    
    do istation=1,Nstations
    turbine%blade(iblade)%QCx(istation)=rR(istation)*turbine%Rmax*Svec(1)!+turbine%blade(iblade)%COR(1)+turbine%dist_from_axis
    turbine%blade(iblade)%QCy(istation)=rR(istation)*turbine%Rmax*Svec(2)!+turbine%blade(iblade)%COR(2)
    turbine%blade(iblade)%QCz(istation)=rR(istation)*turbine%Rmax*Svec(3)!+turbine%blade(iblade)%COR(3)
    if(turbine%IsCounterClockwise) then
        turbine%blade(iblade)%tx(istation)=sin(pitch(istation)/180.0*pi)    
        turbine%blade(iblade)%ty(istation)=-cos(pitch(istation)/180.0*pi)    
        turbine%blade(iblade)%tz(istation)= 0.0
        turbine%blade(iblade)%C(istation)=ctoR(istation)*turbine%Rmax
        turbine%blade(iblade)%thick(istation)=thick(istation)
        turbine%blade(iblade)%pitch(istation)=pitch(istation)/180.0*pi
    elseif(turbine%IsClockwise) then
        turbine%blade(iblade)%tx(istation)=sin(pitch(istation)/180.0*pi)    
        turbine%blade(iblade)%ty(istation)=cos(pitch(istation)/180.0*pi)    
        turbine%blade(iblade)%tz(istation)= 0.0
        turbine%blade(iblade)%C(istation)=ctoR(istation)*turbine%Rmax
        turbine%blade(iblade)%thick(istation)=thick(istation)
        turbine%blade(iblade)%pitch(istation)=pitch(istation)/180.0*pi
        turbine%blade(iblade)%FlipN = .true.
    endif
    !### Do the blade cone angle ###
    ! Rotate coordinates (around y)
    call QuatRot(turbine%blade(iblade)%QCx(istation),turbine%blade(iblade)%QCy(istation),turbine%blade(iblade)%QCz(istation),turbine%blade_cone_angle*pi/180.0d0,&
                0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,0.d0,turbine%blade(iblade)%QCx(istation),turbine%blade(iblade)%QCy(istation),turbine%blade(iblade)%QCz(istation))
    ! Rotate tangential vectors (around y)
    call QuatRot(turbine%blade(iblade)%tx(istation),turbine%blade(iblade)%ty(istation),turbine%blade(iblade)%tz(istation),turbine%blade_cone_angle*pi/180.0d0,&
                0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,0.d0,turbine%blade(iblade)%tx(istation),turbine%blade(iblade)%ty(istation),turbine%blade(iblade)%tz(istation)) 
    
    ! Translate to the COR of each turbine
    turbine%blade(iblade)%QCx(istation)=turbine%blade(iblade)%QCx(istation)+turbine%blade(iblade)%COR(1)
    turbine%blade(iblade)%QCy(istation)=turbine%blade(iblade)%QCy(istation)+turbine%blade(iblade)%COR(2)
    turbine%blade(iblade)%QCz(istation)=turbine%blade(iblade)%QCz(istation)+turbine%blade(iblade)%COR(3)

    end do
    
    call rotate_actuatorline(turbine%blade(iblade),turbine%blade(iblade)%COR,turbine%RotN,(iblade-1)*theta+turbine%AzimAngle)   
 
    call make_actuatorline_geometry(turbine%blade(iblade))
    
    ! Populate element Airfoils 
    call populate_blade_airfoils(turbine%blade(iblade)%NElem,turbine%Blade(iblade)%NAirfoilData,turbine%Blade(iblade)%EAirfoil,turbine%Blade(iblade)%AirfoilData,turbine%Blade(iblade)%ETtoC)
    
    turbine%Blade(iblade)%EAOA_LAST(:)=-666
    turbine%Blade(iblade)%Eepsilon(:)=0.0
    turbine%Blade(iblade)%EEndeffects_factor(:)=1.0
    
    ! Initialise Dynamic stall model
    do ielem=1,turbine%blade(iblade)%Nelem
    if(turbine%blade(iblade)%do_Sheng_stall) then
    call dystl_init(turbine%blade(iblade)%EDynstall(ielem),turbine%blade(iblade)%DynStallFile)
    endif
    if(turbine%blade(iblade)%do_lb_stall) then
    call dystl_init_lb(turbine%blade(iblade)%ELBstall(ielem),turbine%blade(iblade)%DynStallFile)
    endif
    end do 
    end do
   
    ! Rotate the turbine according to the tilt and yaw angle
    ! Yaw
    call rotate_turbine(turbine,(/0.0d0,1.0d0,0.0d0/),turbine%yaw_angle*pi/180.0d0)
    ! Tilt
    call rotate_turbine(turbine,(/0.0d0,0.0d0,1.0d0/),-turbine%shaft_tilt_angle*pi/180.0d0)
   
    ! Set the rotational axis
    call QuatRot(turbine%RotN(1),turbine%RotN(2),turbine%RotN(3),turbine%yaw_angle*pi/180.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,0.d0,&
            turbine%RotN(1),turbine%RotN(2),turbine%RotN(3))
    call QuatRot(turbine%RotN(1),turbine%RotN(2),turbine%RotN(3),-turbine%shaft_tilt_angle*pi/180.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.d0,&
            turbine%RotN(1),turbine%RotN(2),turbine%RotN(3))
    !if(turbine%do_aeroelasticity) then
    !call actuator_line_beam_model_init(turbine%beam,turbine%blade,turbine%NBlades)
    !endif
    
    if (nrank==0) then        
    write(6,*) 'Turbine Name : ', adjustl(turbine%name)
    write(6,*) '-------------------------------------------------------------------'
    write(6,*) 'Number of Blades : ', turbine%Nblades
    write(6,*) 'Origin           : ', turbine%origin
    write(6,*) 'Axis of Rotation : ', turbine%RotN
    write(6,*) '-------------------------------------------------------------------'
    end if
     
    !=========================================================
    ! Create a Tower
    !=========================================================
    if(turbine%has_Tower) then
    call read_actuatorline_geometry(turbine%tower%geom_file,turbine%Towerheight,SVec,rR,ctoR,pitch,thick,Nstations)
    ! Make sure that the spanwise is [0 0 1]
    Svec = (/0.0,1.0,0.0/)
    
    call allocate_actuatorline(Turbine%Tower,Nstations)
    turbine%tower%name=trim(turbine%name)//'_tower'
   
    turbine%Tower%COR=turbine%origin
    turbine%Tower%NElem=Nstations-1  
    turbine%Tower%L=turbine%Towerheight

    do istation=1,Nstations
    turbine%Tower%QCx(istation)= turbine%Tower%COR(1) + turbine%TowerOffset  
    turbine%Tower%QCy(istation)= rR(istation)*turbine%Towerheight*Svec(2) 
    turbine%Tower%QCz(istation)= turbine%Tower%COR(3) 
    turbine%Tower%tx(istation)= 1.0    
    turbine%Tower%ty(istation)= 0.0    
    turbine%Tower%tz(istation)= 0.0
    turbine%Tower%C(istation)=ctoR(istation)*turbine%Towerheight
    turbine%Tower%thick(istation)=thick(istation)
    turbine%Tower%pitch(istation)=pitch(istation)/180.0*pi
    enddo
    
    call make_actuatorline_geometry(turbine%tower)
    
    turbine%tower%EAOA_LAST(:)=-666
    
    !Set the tower body velocity to zero
    turbine%tower%EVbx(:)=0.0
    turbine%tower%EVby(:)=0.0
    turbine%tower%EVbz(:)=0.0
    turbine%tower%EObx(:)=0.0
    turbine%tower%EOby(:)=0.0
    turbine%tower%EObz(:)=0.0
    turbine%tower%Eepsilon(:)=0.0
    endif
     
    !========================================================
    !Compute a number of global parameters for the turbine
    !========================================================
    turbine%angularVel=turbine%Uref*turbine%TSR/turbine%Rmax
    turbine%A=pi*turbine%Rmax**2
    
    turbine%IRotor=0. 
    do iblade=1,turbine%NBlades
    turbine%IRotor=turbine%IRotor+turbine%Blade(iblade)%Inertia
    enddo   
 
    call Compute_Turbine_RotVel(turbine)
    
    return

    end subroutine set_turbine_geometry
   
    subroutine compute_performance(turbine, rho_air)

    implicit none
    type(TurbineType), intent(inout) :: turbine
    real(mytype), intent(in) :: rho_air
    real(mytype) :: Torque_i,FX_i,FY_i,FZ_i,fx_tot,fy_tot,fz_tot,torq_tot 
    real(mytype) :: xe,ye,ze,o1,o2,o3,fx,fy,fz,trx,try,trz,te,ms,sxe,sye,sze
    real(mytype) :: rotx,roty,rotz, rot_vel_mod
    integer :: iblade, ielem
    !write(*,*) 'In calculate_performance'

    RotX=turbine%RotN(1)
    RotY=turbine%RotN(2)
    RotZ=turbine%RotN(3)
    rot_vel_mod = sqrt(RotX*RotX + RotY*RotY + RotZ*RotZ)

    ! Compute Torque for each Blade
        Fx_tot=0.
        Fy_tot=0.
        Fz_tot=0.
        Torq_tot=0

        do iblade=1,turbine%Nblades
        Fx_i=0.
        Fy_i=0.
        Fz_i=0.
        Torque_i=0.

        do ielem=1,turbine%blade(iblade)%NElem
            
            xe=turbine%blade(iblade)%PEX(ielem)
            ye=turbine%blade(iblade)%PEY(ielem)
            ze=turbine%blade(iblade)%PEZ(ielem)
            o1=turbine%blade(iblade)%COR(1)
            o2=turbine%blade(iblade)%COR(2)
            o3=turbine%blade(iblade)%COR(3)
            fx=turbine%blade(iblade)%EFX(ielem)
            fy=turbine%blade(iblade)%EFY(ielem)
            fz=turbine%blade(iblade)%EFZ(ielem)
            ms=turbine%blade(iblade)%EMS(ielem)
            sxe=turbine%blade(iblade)%sex(ielem)
            sye=turbine%blade(iblade)%sey(ielem)
            sze=turbine%blade(iblade)%sez(ielem)

            call cross(xe-o1,ye-o2,ze-o3,fx,fy,fz,trx,try,trz)
            te=(trx*RotX+try*RotY+trz*RotZ)+ms*(sxe*RotX+sye*RotY+sze*RotZ)
            
            Fx_i=Fx_i+fx
            Fy_i=Fy_i+fy
            Fz_i=Fz_i+fz
            Torque_i=Torque_i+te
        end do
            
            Fx_tot=Fx_tot+Fx_i
            Fy_tot=Fy_tot+Fy_i
            Fz_tot=Fz_tot+Fz_i
            Torq_tot=torq_tot+Torque_i
    end do

  
    ! Coefficients and Absolute values
    turbine%CFx=FX_tot/(0.5*rho_air*turbine%A*turbine%Uref**2)
    turbine%CFy=FY_tot/(0.5*rho_air*turbine%A*turbine%Uref**2)
    turbine%CFz=Fz_tot/(0.5*rho_air*turbine%A*turbine%Uref**2)
    turbine%Thrust=sqrt(FX_tot**2.0+FY_tot**2.0+FZ_tot**2.0)
    turbine%CT=sqrt(turbine%CFx**2.0+turbine%CFy**2.0+turbine%CFz**2.0)
    !turbine%Thrust=(Fx_tot*RotX +Fy_tot*RotY +Fz_tot*RotZ)/rot_vel_mod
    !turbine%CT=turbine%Thrust/(0.5*rho_air*turbine%A*turbine%Uref**2)
    turbine%torque=Torq_tot
    turbine%CTR=Torq_tot/(0.5*rho_air*turbine%A*turbine%Rmax*turbine%Uref**2.0)
    turbine%Power=abs(Torq_tot)*turbine%angularVel
    turbine%CP= abs(turbine%CTR)*turbine%TSR
    
    ! PRINT ON SCREEN
    if(nrank==0) then
        write(6,*) "Turbine : ",   turbine%name
        write(6,*) "======================================="
        write(6,*) "Thrust coeff : ",   turbine%CT
        write(6,*) "Power  coeff  : ",  turbine%CP
        write(6,*) "Thrust  : ",  turbine%Thrust/1000., "kN"
        write(6,*) "Torque  : ",  turbine%Torque/1000., "kN m"
        write(6,*) "Power   : ",  turbine%Power/1000000., "MW"
    endif

    end subroutine compute_performance
    
    subroutine Compute_Turbine_EndEffects(turbine)
    
    implicit none
    type(TurbineType),intent(inout) :: turbine
    integer :: iblade,ielem
    real(mytype) ::g1,alpha,pitch,F,Froot,Ftip,rtip, rroot, phi, sphi,axis_mag
    
    
    do iblade=1,turbine%Nblades
    ! Compute angle phi at the tip 
    do ielem=1,turbine%blade(iblade)%Nelem
        
        phi=turbine%blade(iblade)%EAOA(ielem)+turbine%blade(iblade)%Epitch(ielem)
        !if (nrank==0) print *, "AoA = ", turbine%blade(iblade)%EAOA(ielem), "Pitch =", turbine%blade(iblade)%EPitch(ielem), "phi =", phi
        rroot=turbine%blade(iblade)%ERdist(ielem)/turbine%Rmax
        rtip=(turbine%Rmax-turbine%blade(iblade)%ERdist(ielem))/turbine%Rmax
        
        ! Set them into 1.0 before the calculation
        Ftip=1.0
        Froot=1.0
        if (turbine%do_tip_correction) then 
            if (turbine%EndEffectModel_is_Glauret) then
                g1=1.0
            else if (turbine%EndEffectModel_is_Shen) then
                g1=dexp(-turbine%ShenCoeff_c1*(turbine%NBlades*turbine%TSR-turbine%ShenCoeff_c2))+0.1
            else
                write(*,*) "Only Glauret and ShenEtAl2005 are available at the moment"
                stop
            endif
            Ftip=2.0/pi*dacos(dexp(-g1*turbine%Nblades/2.0*(1.0/rroot-1.0)/dsin(phi)))
            if (dabs(dexp(-g1*turbine%Nblades/2.0*(1.0/rroot-1.0)/dsin(phi)))>1.) then
                if (nrank==0) print *, "Something went wrong with the tip correction model -- phi =", phi
                Ftip=1.
            endif
        endif
        if (turbine%do_root_correction) then
            if (turbine%EndEffectModel_is_Glauret) then
                g1=1.0
            else if (turbine%EndEffectModel_is_Shen) then
                g1=dexp(-turbine%ShenCoeff_c1*(turbine%NBlades*turbine%TSR-turbine%ShenCoeff_c2))+0.1
            else
                write(*,*) "Only Glauret and ShenEtAl2005 are available at the moment"
                stop
            endif
            Froot=2.0/pi*dacos(dexp(-g1*turbine%Nblades/2.0*(1.0/rtip-1.0)/dsin(phi)))
            if (dabs(dexp(-g1*turbine%Nblades/2.0*(1.0/rtip-1.0)/dsin(phi)))>1.) then
                if (nrank==0) print *, "Something went wrong with the root correction model -- phi =", phi
                Froot=1.
            endif
        endif

            turbine%blade(iblade)%EEndeffects_factor(ielem)=Ftip*Froot 
        end do
    end do

    end subroutine Compute_Turbine_EndEffects

    subroutine Compute_Turbine_RotVel(turbine)
    implicit none
    type(TurbineType),intent(inout) :: turbine
    integer :: iblade,ielem
    real(mytype) :: wRotX,wRotY,wRotZ,Rx,Ry,Rz,ublade,vblade,wblade
    
    !========================================================
    ! Compute Element local rotational velocity
    !========================================================
    wRotX=turbine%angularVel*turbine%RotN(1)
    wRotY=turbine%angularVel*turbine%RotN(2)
    wRotZ=turbine%angularVel*turbine%RotN(3)
    
    do iblade=1,turbine%NBlades
    do ielem=1,turbine%Blade(iblade)%Nelem
     
    Rx=-turbine%blade(iblade)%COR(1)+turbine%Blade(iblade)%PEx(ielem);
    Ry=-turbine%blade(iblade)%COR(2)+turbine%Blade(iblade)%PEy(ielem);
    Rz=-turbine%blade(iblade)%COR(3)+turbine%Blade(iblade)%PEz(ielem);

   ! ! Find the cross product Ublade = Omega x R
    call cross(wRotX,wRotY,wRotZ,Rx,Ry,Rz,ublade,vblade,wblade)
    
    turbine%Blade(iblade)%EVbx(ielem)=ublade
    turbine%Blade(iblade)%EVby(ielem)=vblade
    turbine%Blade(iblade)%EVbz(ielem)=wblade
    turbine%Blade(iblade)%EObx(ielem)=wRotX
    turbine%Blade(iblade)%EOby(ielem)=wRotY
    turbine%Blade(iblade)%EObz(ielem)=wRotZ
    
    end do
    end do
    
    end subroutine Compute_Turbine_RotVel

    subroutine rotate_turbine(turbine,Axis,theta)
    
    implicit none
    type(TurbineType),intent(inout) :: turbine
    real(mytype),intent(in) :: Axis(3)
    real(mytype),intent(in) :: theta
    real(mytype) :: nrx,nry,nrz,px,py,pz 
    integer :: j,ielem,i
    real(mytype) :: vrx,vry,vrz,VMag
    real(mytype) :: xtmp,ytmp,ztmp, txtmp, tytmp, tztmp
    ! Rotates data in blade arrays. Rotate element end geometry and recalculate element geometry.

        
    nrx=Axis(1)
    nry=Axis(2)
    nrz=Axis(3) 
    px =turbine%origin(1)
    py =turbine%origin(2)
    pz =turbine%origin(3)


    do j=1,turbine%NBlades
        do ielem=1,turbine%Blade(j)%Nelem+1
        ! Blade end locations (quarter chord). xBE(MaxSegEnds)
        xtmp=turbine%Blade(j)%QCx(ielem)
        ytmp=turbine%Blade(J)%QCy(ielem)
        ztmp=turbine%Blade(j)%QCz(ielem)
        
        Call QuatRot(xtmp,ytmp,ztmp,theta,nrx,nry,nrz,px,py,pz,vrx,vry,vrz)
        turbine%Blade(j)%QCx(ielem)=vrx                                       
        turbine%Blade(j)%QCy(ielem)=vry                                       
        turbine%Blade(j)%QCz(ielem)=vrz                                  
        
        txtmp=turbine%Blade(j)%tx(ielem)
        tytmp=turbine%Blade(j)%ty(ielem)
        tztmp=turbine%Blade(j)%tz(ielem)
        
        ! Tangent vectors (rotated assuming the origina as 0)
        Call QuatRot(txtmp,tytmp,tztmp,theta,nrx,nry,nrz,0.0d0,0.0d0,0.0d0,vrx,vry,vrz)
        VMag=sqrt(vrx**2+vry**2+vrz**2)
        turbine%Blade(j)%tx(ielem)=vrx/VMag                                      
        turbine%Blade(j)%ty(ielem)=vry/VMag                                  
        turbine%Blade(j)%tz(ielem)=vrz/VMag                                       
  
        end do
        
        call make_actuatorline_geometry(turbine%Blade(j))
    end do 
    
    end subroutine rotate_turbine  
   
    subroutine read_list_controller_file(FN,turbine) 
    
    implicit none
    character(len=100),intent(in)  :: FN ! FileName of the geometry file
    type(TurbineType),intent(inout) :: turbine
    integer :: i
    character(1000) :: ReadLine
    
    open(22,file=FN)

    ! Read the Number of List Controller Entries 
    read(22,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) turbine%ContEntries

    allocate(turbine%ContWindSpeed(turbine%ContEntries),turbine%ContOmega(turbine%ContEntries), turbine%ContPitch(turbine%ContEntries))
    ! Read the stations specs
    do i=1,turbine%ContEntries
    
    read(22,'(A)') ReadLine ! Blade ....

    read(ReadLine,*) turbine%ContWindSpeed(i), turbine%ContOmega(i), turbine%ContPitch(i)

    end do
    
    close(22)

    end subroutine read_list_controller_file 

    subroutine from_list_controller(Omega,pitch,turbine,WindSpeed)

        implicit none
        type(TurbineType),intent(inout) :: turbine
        real(mytype),intent(in) :: WindSpeed
        real(mytype),intent(out) :: Omega,pitch
        real(mytype) :: mindiff, inter
        integer :: i,j,inear
        integer :: ilower,iupper

        
        if(WindSpeed<=turbine%ContWindSpeed(1)) then
            Omega=turbine%ContOmega(1) ! First entry
            Pitch=turbine%ContPitch(1) ! First entry
        else if(WindSpeed>=turbine%ContWindSpeed(turbine%ContEntries)) then
            Omega=turbine%ContOmega(turbine%ContEntries) ! Last entry
            Pitch=turbine%ContPitch(turbine%ContEntries) ! Last entry
        
        else !> Find the nearest list value
        mindiff=1e6
        do i=1,turbine%ContEntries
            if(abs(WindSpeed-turbine%ContWindSpeed(i))<mindiff) then
               mindiff=abs(WindSpeed-turbine%ContWindSpeed(i)) 
               inear=i
            endif 
        enddo
        
        !> Find the upper and lower indices
        if(turbine%ContWindSpeed(inear)<WindSpeed) then
            ilower=inear
            iupper=inear+1
        else
            ilower=inear-1
            iupper=inear
        endif

        !> Interpolate
        inter=(WindSpeed-turbine%ContWindSpeed(ilower))/(turbine%ContWindSpeed(iupper)-turbine%ContWindSpeed(ilower))
        Omega=turbine%ContOmega(ilower)+inter*(turbine%ContOmega(iupper)-turbine%ContOmega(ilower))
        Pitch=turbine%ContPitch(ilower)+inter*(turbine%ContPitch(iupper)-turbine%ContPitch(ilower))

        endif
        
        return

    end subroutine from_list_controller

    subroutine compute_rotor_upstream_velocity(turbine)

use actuator_line_model_utils ! used only for the trilinear interpolation
USE param 
USE decomp_2d
use var
use MPI

        implicit none
        type(TurbineType),intent(inout) ::turbine
        integer :: iblade, ielem,Nelem
        real(mytype) :: Rupstream(3)
        real(mytype) :: Ux,Uy,Uz,Phixy,Phixz
        real(mytype) :: Ux_part, Uy_part, Uz_part, Phixy_part, Phixz_part
	    real(mytype) :: ymin,ymax,zmin,zmax
	    real(mytype) :: xmesh,ymesh,zmesh
	    real(mytype) :: dist, min_dist 
        integer :: min_i,min_j,min_k
        integer :: i,j,k,code

        Ux=0.
        Uy=0.
        Uz=0.
       
        ! Velocity is calculated at a probe point (closest) at one D upstream the turbine

        Rupstream(:)=turbine%origin(:)-turbine%rotN(:)*2.*turbine%Rmax       
        if (istret.eq.0) then 
        ymin=(xstart(2)-1)*dy-dy/2.0 ! Add -dy/2.0 overlap
        ymax=(xend(2)-1)*dy+dy/2.0   ! Add +dy/2.0 overlap
        else
        ymin=yp(xstart(2))
        ymax=yp(xend(2))
        endif
        
        zmin=(xstart(3)-1)*dz-dz/2.0 ! Add a -dz/2.0 overlap
        zmax=(xend(3)-1)*dz+dz/2.0   ! Add a +dz/2.0 overlap
         
        
        min_dist=1e6
        if((Rupstream(2)>=ymin).and.(Rupstream(2)<=ymax).and.(Rupstream(3)>=zmin).and.(Rupstream(3)<=zmax)) then
            !write(*,*) 'Warning: I own this node'
            do k=xstart(3),xend(3)
            zmesh=(k-1)*dz 
            do j=xstart(2),xend(2)
            if (istret.eq.0) ymesh=(j-1)*dy
            if (istret.ne.0) ymesh=yp(j)
            do i=xstart(1),xend(1)
            xmesh=(i-1)*dx
            dist = sqrt((Rupstream(1)-xmesh)**2.+(Rupstream(2)-ymesh)**2.+(Rupstream(3)-zmesh)**2.) 
            
            if (dist<min_dist) then
                min_dist=dist
                min_i=i
                min_j=j
                min_k=k
            endif

            enddo
            enddo
            enddo
            
            Ux_part=ux1(min_i,min_j,min_k)
            Uy_part=uy1(min_i,min_j,min_k)
            Uz_part=uz1(min_i,min_j,min_k)
            Phixy_part=0.0
            Phixz_part=0.0    
        else
            Ux_part=0.0
            Uy_part=0.0
            Uz_part=0.0
            Phixy_part=0.0
            Phixz_part=0.0
            !write(*,*) 'Warning: I do not own this node' 
        endif
           
        call MPI_ALLREDUCE(Ux_part,Ux,1,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,code)
        if (code.ne.0) call decomp_2d_abort(code, "MPI_ALLREDUCE")
        call MPI_ALLREDUCE(Uy_part,Uy,1,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,code)
        if (code.ne.0) call decomp_2d_abort(code, "MPI_ALLREDUCE")
        call MPI_ALLREDUCE(Uz_part,Uz,1,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,code)
        if (code.ne.0) call decomp_2d_abort(code, "MPI_ALLREDUCE")
        
        Turbine%Ux_upstream=Ux
        Turbine%Uy_upstream=Uy
        Turbine%Uz_upstream=Uz

        Turbine%Uref=sqrt(Ux**2.0+Uy**2.0+Uz**2.0)
        return
    
    end subroutine Compute_Rotor_upstream_Velocity

end module actuator_line_turbine 
