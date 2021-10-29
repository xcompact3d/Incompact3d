module actuator_line_element

    use decomp_2d, only: mytype
    use actuator_line_model_utils 
    use airfoils
    use dynstall_legacy
    use dynstall

type ActuatorLineType
    integer :: NElem                            ! Number of Elements of the Blade
    character(len=100) :: name                  ! Actuator line name
    character(len=100) :: geom_file             ! Actuator line file name (is not used for the turbines)
    character(len=100) :: dynstallfile          ! Dynstallfile to load options

    logical :: Is_constant_circulation=.false.  ! Is constant circulation flag for verification purposes

    ! Station parameters
    logical :: FlipN =.false.                   ! Flip Normal
    real(mytype), allocatable :: QCx(:)         ! Blade quarter-chord line x coordinates at element ends (stations)
    real(mytype), allocatable :: QCy(:)         ! Blade quarter-chord line y coordinates at element ends (stations)
    real(mytype), allocatable :: QCz(:)         ! Blade quarter-chord line z coordinates at element ends (stations)
    real(mytype), allocatable :: tx(:)          ! Blade unit tangent vector (rearward chord line direction) x-comp at element ends
    real(mytype), allocatable :: ty(:)          ! Blade unit tangent vector (rearward chord line direction) y-comp at element ends 
    real(mytype), allocatable :: tz(:)          ! Blade unit tangent vector (rearward chord line direction) z-comp at element ends  
    real(mytype), allocatable :: C(:)           ! Blade chord length at element ends
    real(mytype), allocatable :: thick(:)       ! Blade thickness at element ends
    real(mytype), allocatable :: pitch(:)       ! Blade station pitch at element ends

    ! Element parameters 
    real(mytype), allocatable :: PEx(:)         ! Element centre x coordinates
    real(mytype), allocatable :: PEy(:)         ! Element centre y coordinates
    real(mytype), allocatable :: PEz(:)         ! Element centre z coordinates
    real(mytype), allocatable :: tEx(:)         ! Element unit tangent vector (rearward chord line direction) x-component
    real(mytype), allocatable :: tEy(:)         ! Element unit tangent vector (rearward chord line direction) y-component
    real(mytype), allocatable :: tEz(:)         ! Element unit tangent vector (rearward chord line direction) z-component
    real(mytype), allocatable :: nEx(:)         ! Element unit normal vector x-component
    real(mytype), allocatable :: nEy(:)         ! Element unit normal vector y-component
    real(mytype), allocatable :: nEz(:)         ! Element unit normal vector z-component
    real(mytype), allocatable :: sEx(:)         ! Element unit spanwise vector x-component 
    real(mytype), allocatable :: sEy(:)         ! Element unit spanwise vector y-component
    real(mytype), allocatable :: sEz(:)         ! Element unit spanwise vector z-component
    real(mytype), allocatable :: EC(:)          ! Element chord length
    real(mytype), allocatable :: EDS(:)         ! Element spanwise distance (length)
    real(mytype), allocatable :: EArea(:)       ! Element Area
    real(mytype), allocatable :: Eepsilon(:)    ! Element Force Projection Parameter
    real(mytype), allocatable :: ERdist(:)      ! Element Distance from the origin 
    real(mytype), allocatable :: ETtoC(:)       ! Element thickness to Chord ratio
    
    ! Angle of Attack, Pitch, local Reynolds number and relative velocity
    real(mytype), allocatable :: Epitch(:)      ! Element pitch angle
    real(mytype), allocatable :: EAOA(:)        ! Element Current angle of Attack (used in added mass terms)
    real(mytype), allocatable :: EAOAdot(:)     ! Element AOA rate of change
    real(mytype), allocatable :: EUn(:)         ! Element Current normal velocity (used in added mass terms)
    real(mytype), allocatable :: EUndot(:)      ! Element Current normal velocity (used in added mass terms)
    real(mytype), allocatable :: EAOA_LAST(:)   ! Element Last angle of Attack (used in added mass terms)
    real(mytype), allocatable :: EUn_LAST(:)    ! Element Last normal velocity (used in added mass terms)
    real(mytype), allocatable :: ERe(:)         ! Element local Reynolds number
    real(mytype), allocatable :: EUr(:)         ! Element local velocity magnitude
    
    ! Velocity of the Fluid at the actuator Line Locations
    real(mytype), allocatable :: EVx(:)         ! Element Local fluid Velocity in the global x-direction
    real(mytype), allocatable :: EVy(:)         ! Element Local fluid Velocity in the global y-direction
    real(mytype), allocatable :: EVz(:)         ! Element Local fluid Velocity in the global z-direction
     
    ! Body Velocity of the Actuator Line 
    real(mytype), allocatable :: EVbx(:)        ! Element Local body Velocity in the global x-direction
    real(mytype), allocatable :: EVby(:)        ! Element Local body Velocity in the global y-direction
    real(mytype), allocatable :: EVbz(:)        ! Element Local body Velocity in the global z-direction
    real(mytype), allocatable :: EObx(:)        ! Element Local body angular Velocity in the global x-direction
    real(mytype), allocatable :: EOby(:)        ! Element Local body angular Velocity in the global y-direction
    real(mytype), allocatable :: EObz(:)        ! Element Local body angular Velocity in the global z-direction
    
    ! Element Forces CD, CL CM25 
    real(mytype), allocatable :: ECD(:)         ! Element Drag Coefficient
    real(mytype), allocatable :: ECL(:)         ! Element Lift Coefficient 
    real(mytype), allocatable :: ECLcirc(:)     ! Element Circulation Lift Coefficient
    real(mytype), allocatable :: ECM(:)         ! Element Moment Coefficient
    real(mytype), allocatable :: ECN(:)         ! Element Normal Force Coefficient
    real(mytype), allocatable :: ECT(:)         ! Element Tangential Force Coefficient 
    
    ! Element Circulation Gamma 
    real(mytype), allocatable :: EGamma(:)      ! Element Circulation (Gamma)

    ! Element Forces in the nts direction
    real(mytype), allocatable :: EFn(:)         ! Element Force in the normal direction
    real(mytype), allocatable :: EFt(:)         ! Element Force in the tangential direction (rearward chord line direction) 
    real(mytype), allocatable :: EMS(:)         ! Element Moment 

    ! Element Forces and Torque in the xyz direction
    real(mytype), allocatable :: EFx(:)         ! Element Force in the global x-direction
    real(mytype), allocatable :: EFy(:)         ! Element Force in the global y-direction
    real(mytype), allocatable :: EFz(:)         ! Element Force in the global z-direction

    ! Influence Matrix for rbf interpolations
    real(mytype), allocatable :: A_rbf(:,:)     

    real(mytype), allocatable :: EEndeffects_factor(:) ! End effects factor for the blade (initialize as one)

    ! Element Airfoil Data
    type(AirfoilType), allocatable :: EAirfoil(:)      ! Element Airfoil 
    integer :: NAirfoilData
    type(AirfoilType), allocatable :: AirfoilData(:)   ! Element Airfoil 
    type(DS_Type), allocatable :: EDynstall(:)         ! Element Dynamic Stall Model: Sheng et al. (2008)
    type(LB_Type), allocatable :: ELBStall(:)          ! Element Dynamic Stall Model: classic Leishman-Beddoes (taken from CACTUS)
    
    ! Forces and Torques on the ActuatorLine 
    real(mytype) :: Fx      ! Element Force in the global x-direction
    real(mytype) :: Fy      ! Element Force in the global y-direction
    real(mytype) :: Fz      ! Element Force in the global z-direction

    real(mytype) :: Area    ! Effective Airfoil Area
    real(mytype) :: L 
    real(mytype) :: Inertia ! Moment of inertia

    ! Degrees of Freedom    
    real(mytype) :: COR(3)       ! Center of Rotation
    real(mytype) :: SpanWise(3)  ! Point of Rotation
    real(mytype) :: GammaCirc    ! GammaCirc (for the constant circulation simulations)
   
    ! Unsteady Loading
    logical :: do_added_mass=.false.
    logical :: do_Sheng_stall=.false.
    logical :: do_LB_stall=.false.
    logical :: do_DynStall_AlphaEquiv=.false.
    logical :: do_random_walk_forcing=.false.

    ! Blade/Actuator_line Pitch
    logical :: pitch_control=.false.
    real(mytype) :: pitch_start_time
    real(mytype) :: pitch_end_time

    ! Harmonic Pitch control parameters
    real(mytype) :: angular_pitch_freq
    real(mytype) :: pitch_angle_init
    real(mytype) :: pitchAmp

end type ActuatorLineType

contains

    !*******************************************************************************
    !
    subroutine set_actuatorline_geometry(actuatorline)
    !
    !*******************************************************************************

      use param, only: zero, one
      use dbg_schemes, only: cos_prec, sin_prec

      implicit none
      type(ActuatorLineType), intent(inout) :: actuatorline
      real(mytype), allocatable :: rR(:),ctoR(:),pitch(:),thick(:)
      real(mytype) :: SVec(3), length
      integer :: Nstations, Istation, ielem

      if (nrank==0) then
         write(*,*) 'Actuator line name : ', actuatorline%name
      endif

      call read_actuatorline_geometry(actuatorline%geom_file,length,SVec,rR,ctoR,pitch,thick,Nstations)
    
      call allocate_actuatorline(actuatorline,Nstations)
   
      actuatorline%SpanWise=SVec
      actuatorline%L=length

      ! The directions of vectors etc are just hacked...
      actuatorline%Nelem=Nstations-1
      do istation=1,Nstations
         actuatorline%QCx(istation)=rR(istation)*length*Svec(1)+actuatorline%COR(1)
         actuatorline%QCy(istation)=rR(istation)*length*Svec(2)+actuatorline%COR(2)
         actuatorline%QCz(istation)=rR(istation)*length*Svec(3)+actuatorline%COR(3)
    
         if (actuatorline%pitch_control) then
            actuatorline%pitch(istation)=actuatorline%pitch_angle_init*conrad
         else
            actuatorline%pitch(istation)=pitch(istation)*conrad
         endif
 
         actuatorline%tx(istation)= cos_prec(actuatorline%pitch(istation))
         actuatorline%ty(istation)=-sin_prec(actuatorline%pitch(istation))
         actuatorline%tz(istation)= zero
         actuatorline%C(istation)=ctoR(istation)*length
         actuatorline%thick(istation)=thick(istation)
      enddo
      actuatorline%FlipN=.false.

      call make_actuatorline_geometry(actuatorline)
    
      ! Compute the area
      do ielem=1,actuatorline%Nelem
         actuatorline%Area=actuatorline%Area+actuatorline%EArea(ielem)
      enddo
    
      ! Initial values for the body linear and angular velocities
      actuatorline%EVbx(:)=zero
      actuatorline%EVby(:)=zero
      actuatorline%EVbz(:)=zero
      actuatorline%EObx(:)=zero
      actuatorline%EOby(:)=zero
      actuatorline%EObz(:)=zero
      actuatorline%Eepsilon(:)=zero
      actuatorline%EEndeffects_factor(:)=one
    
      ! Populate element airfoils
      call populate_blade_airfoils(actuatorline%NElem,actuatorline%NAirfoilData,actuatorline%EAirfoil,actuatorline%AirfoilData,actuatorline%ETtoC)
    
      actuatorline%EAOA_LAST(:)=-666.0_mytype
    
      ! If Sheng dynamic stall is enabled
      if (actuatorline%do_Sheng_stall) then
         do ielem=1,actuatorline%Nelem
            call dystl_init(actuatorline%EDynstall(ielem),actuatorline%dynstallfile)
         enddo
      endif
    
      ! If LB dynamic stall is enabled
      if (actuatorline%do_lb_stall) then
         do ielem=1,actuatorline%Nelem
            call dystl_init_lb(actuatorline%ELBStall(ielem),actuatorline%dynstallfile)
         enddo
      endif

    end subroutine set_actuatorline_geometry
    
    !*******************************************************************************
    !
    subroutine Compute_ActuatorLine_Forces(act_line,rho_air,visc,dt,time)
    !
    !*******************************************************************************

      use dbg_schemes, only: sqrt_prec, sin_prec, cos_prec
      use param, only: zpone, zpfive, zero, one, two, four, eight

      implicit none
      type(ActuatorLineType), intent(inout) :: act_line
      real(mytype), intent(in) :: rho_air,visc,dt,time
      real(mytype) :: wRotX,wRotY,wRotZ,ub,vb,wb,u,v,w
      real(mytype) :: xe,ye,ze,nxe,nye,nze,txe,tye,tze,sxe,sye,sze,ElemArea,ElemChord
      real(mytype) :: urdn,urdc,ur,alpha,ds
      real(mytype) :: CL,CD,CN,CT,CM25,MS,FN,FT,FX,Fy,Fz
      real(mytype) :: dal,dUn
      real(mytype) :: CLstat, CDstat
      real(mytype) :: CLdyn,CDdyn, CM25dyn, CNAM,CTAM,CMAM
      real(mytype) :: rand(3000), freq, Strouhal ! To add random walk on the lift/drag coefficients
      integer :: ielem
  
      call random_number(rand)
     
      do ielem=1,act_line%NElem
         ! Assign global values to local values (to make life easier)
         nxe=act_line%nEx(ielem)
         nye=act_line%nEy(ielem)
         nze=act_line%nEz(ielem)
         txe=act_line%tEx(ielem)
         tye=act_line%tEy(ielem)
         tze=act_line%tEz(ielem)
         sxe=act_line%sEx(ielem)
         sye=act_line%sEy(ielem)
         sze=act_line%sEz(ielem)
         ElemArea=act_line%EArea(ielem)
         ElemChord=act_line%EC(ielem)
         u=act_line%EVx(ielem)
         v=act_line%EVy(ielem)
         w=act_line%EVz(ielem)

         ! Solid body rotation of the elements
         ub=act_line%EVbx(ielem)
         vb=act_line%EVby(ielem)
         wb=act_line%EVbz(ielem)
         wRotx=act_line%EObx(ielem)
         wRoty=act_line%EOby(ielem)
         wRotz=act_line%EObz(ielem)

         ! Calculate element normal and tangential velocity components
         urdn=nxe*(u-ub)+nye*(v-vb)+nze*(w-wb) ! Normal
         urdc=txe*(u-ub)+tye*(v-vb)+tze*(w-wb) ! Tangential
    
         ur=sqrt_prec(urdn**2.0+urdc**2.0)
         act_line%EUr(ielem)=ur
    
         ! This is the dynamic angle of attack
         alpha=atan2(urdn,urdc)
         act_line%ERe(ielem)=ur*ElemChord/visc
    
         act_line%EAOA(ielem)=alpha
         act_line%EUn(ielem)=urdn
    
         if (act_line%EAOA_Last(ielem)<-665.0_mytype) then
            dal=zero
            dUn=zero
         else
            dal=act_line%EAOA(ielem)-act_line%EAOA_Last(ielem)
            dUn=act_line%EUn(ielem)-act_line%EUn_LAST(ielem)
         endif
         act_line%EAOAdot(ielem)=dal/dt
         act_line%EUndot(ielem)=dUn/dt
         
         ! Compute the airfoil coefficients
         call compute_StaticLoads(act_line%EAirfoil(ielem),alpha,act_line%ERe(ielem),CN,CT,CM25,CL,CD)
 
         ! Correct for dynamic stall
         if (act_line%do_Sheng_stall) then
            call DynstallCorrect(act_line%EDynstall(ielem),act_line%EAirfoil(ielem),time,dt,act_line%EUr(ielem),&
                   ElemChord,alpha,act_line%ERe(ielem),CLdyn,CDdyn,CM25dyn)
            CL=CLdyn
            CD=CDdyn
            CM25=CM25dyn
         endif
   
         if (act_line%do_LB_stall) then
            CLstat=CL
            CDstat=CD
            call LB_DynStall(act_line%EAirfoil(ielem),act_line%ELBStall(ielem),CLstat,CDstat,alpha,alpha,act_line%ERe(ielem),CLdyn,CDdyn)
            CL=CLdyn
            CD=CDdyn
            ds=two*act_line%EUr(ielem)*dt/ElemChord
            call LB_UpdateStates(act_line%ELBStall(ielem),act_line%EAirfoil(ielem),act_line%ERe(ielem),ds)
         endif
    
         ! Correct for added mass
         if (act_line%do_added_mass) then
            CNAM=-pi*Elemchord*act_line%EUndot(ielem)/(eight*ur**2)
            CTAM= pi*Elemchord*act_line%EAOAdot(ielem)*urdn/(eight*ur*2)
            CMAM=-CNAM/four-urdn*urdc/(eight*ur**2)
            CT=CT+CTAM
            CN=CN+CNAM
            CL=CN*cos_prec(alpha)-CT*sin_prec(alpha)
            CD=CN*sin_prec(alpha)+CT*cos_prec(alpha)
            CM25=CM25+CMAM
         endif
   
         ! Apply random walk on the lift and drag forces
         if (act_line%do_random_walk_forcing) then
            Strouhal=0.17_mytype
            freq=Strouhal*ur/max(ElemChord,0.0001_mytype)
            CL=CL*(one+zpone*sin_prec(two*pi*freq*time)+0.05_mytype*(-one+two*rand(ielem)))
            CD=CD*(one+0.05_mytype*(-one+two*rand(ielem)))
         endif

         if (act_line%Is_constant_circulation) then
            CL=two*act_line%GammaCirc/(ur*ElemChord)
            CD=zero
            MS=zero
         endif

         ! Tangential and normal coeffs
         CN=CL*cos_prec(alpha)+CD*sin_prec(alpha)
         CT=-CL*sin_prec(alpha)+CD*cos_prec(alpha)

         ! Apply coeffs to calculate tangential and normal Forces
         FN=zpfive*rho_air*CN*ElemArea*ur**2.0
         FT=zpfive*rho_air*CT*ElemArea*ur**2.0
         MS=zpfive*rho_air*CM25*ElemChord*ElemArea*ur**2.0
    
         ! Correct for tip/root loss effects
         FN=FN*act_line%EEndeffects_factor(ielem)
         FT=FT*act_line%EEndeffects_factor(ielem)
   
         ! Compute forces in the X, Y, Z axis and torque
         FX=FN*nxe+FT*txe
         FY=FN*nye+FT*tye
         FZ=FN*nze+FT*tze
         
         !==========================================
         ! Assign the derived types
         !==========================================
    
         ! Local Load Coefficients
         act_line%ECD(ielem)=CD
         act_line%ECL(ielem)=CL
         act_line%ECM(ielem)=CM25
         act_line%ECN(ielem)=CN
         act_line%ECT(ielem)=CT
    
         ! Local Coordinate-system Forces
         act_line%EFn(ielem)=FN
         act_line%EFt(ielem)=FT
         act_line%EMS(ielem)=MS
    
         ! Global Forces and Torques
         act_line%EFx(ielem)=FX
         act_line%EFy(ielem)=FY
         act_line%EFz(ielem)=FZ
    
         ! Set the AOA_LAST before exiting the routine
         act_line%EAOA_LAST(ielem)=alpha
         act_line%EUn_last(ielem)=urdn
         act_line%EGamma(ielem)=CL*ElemChord*ur/two
      enddo

    end subroutine compute_Actuatorline_Forces
    
    !*******************************************************************************
    !
    subroutine compute_Tower_Forces(tower,rho_air,visc,time,CL,CD,Str)
    !
    !*******************************************************************************

      use param, only: zero, zpfive, one, two
      use dbg_schemes, only: sin_prec, cos_prec

      implicit none
      type(ActuatorLineType), intent(inout) :: tower
      real(mytype), intent(in) :: rho_air,visc, time, CL, CD, Str
      real(mytype) :: rand(3000)
      real(mytype) :: xe,ye,ze,nxe,nye,nze,txe,tye,tze,sxe,sye,sze,ElemArea
      real(mytype) :: u,v,w,ub,vb,wb,urdn,urdc,ur,Diameter,freq,alpha
      real(mytype) :: CN,CT,FN,FT,FX,Fy,Fz
      integer :: ielem
    
      call random_number(rand)

      do ielem=1,tower%NElem
         xe=tower%PEX(ielem)
         ye=tower%PEY(ielem)
         ze=tower%PEZ(ielem)
         nxe=tower%nEx(ielem)
         nye=tower%nEy(ielem)
         nze=tower%nEz(ielem)
         txe=tower%tEx(ielem)
         tye=tower%tEy(ielem)
         tze=tower%tEz(ielem)
         sxe=tower%sEx(ielem)
         sye=tower%sEy(ielem)
         sze=tower%sEz(ielem)
         Diameter=tower%EC(ielem)
         ElemArea=tower%EArea(ielem)
         u=tower%EVx(ielem)
         v=tower%EVy(ielem)
         w=tower%EVz(ielem)
         ub=zero
         vb=zero
         wb=zero
           
         ! Calculate element normal and tangential velocity components
         urdn=nxe*(u-ub)+nye*(v-vb)+nze*(w-wb) ! Normal
         urdc=txe*(u-ub)+tye*(v-vb)+tze*(w-wb) ! Tangential

         ur=sqrt(urdn**2.0+urdc**2.0)
         tower%EUr(ielem)=ur
         alpha=atan2(urdn,urdc)
         tower%EAOA(ielem)=alpha
         tower%ERE(ielem)=ur*Diameter/visc

         ! Apply random walk on the lift and drag forces
         freq=Str*ur/max(Diameter,0.0001_mytype)
         tower%ECL(ielem)=CL*sin_prec(two*freq*pi*time)
         tower%ECL(ielem)=tower%ECL(ielem)*(one+0.05_mytype*(-one+two*rand(ielem)))
         tower%ECD(ielem)=CD

         ! Tangential and normal coeffs
         CN=tower%ECL(ielem)*cos_prec(alpha)+tower%ECD(ielem)*sin_prec(alpha)
         CT=-tower%ECL(ielem)*sin_prec(alpha)+tower%ECD(ielem)*cos_prec(alpha)

         ! Apply Coeffs to calculate tangential and normal Forces
         FN=zpfive*rho_air*CN*ElemArea*ur**2.0
         FT=zpfive*rho_air*CT*ElemArea*ur**2.0

         ! Compute forces in the X, Y, Z axis and torque
         FX=FN*nxe+FT*txe
         FY=FN*nye+FT*tye
         FZ=FN*nze+FT*tze

         ! Assign the derived types
         tower%ECN(ielem)=CN
         tower%ECT(ielem)=CT
         tower%EFN(ielem)=FN
         tower%EFT(ielem)=FT
         tower%EMS(ielem)=zero
         tower%EFX(ielem)=FX
         tower%EFY(ielem)=FY
         tower%EFZ(ielem)=FZ
      enddo

    end subroutine compute_Tower_Forces

    !*******************************************************************************
    !
    subroutine pitch_actuator_line(act_line,pitch_angle)
    !
    !*******************************************************************************

      use dbg_schemes, only: sqrt_prec

      implicit none
      type(ActuatorLineType), intent(inout) :: act_line
      real(mytype), intent(inout) :: pitch_angle ! Pitch in degrees
      real(mytype) :: S, sx, sy, sz, tx,ty,tz, new_tx, new_ty, new_tz, tmag
      integer :: istation, Nstation
    
      ! Change the pitch angle by changing n,t and s unit vectors
      Nstation=act_line%Nelem+1
   
      ! Compute the blade-wise unit vector
      S=sqrt_prec((act_line%QCX(NStation)-act_line%COR(1))**2. + &
                  (act_line%QCY(NStation)-act_line%COR(2))**2. + &
                  (act_line%QCZ(NStation)-act_line%COR(3))**2.)
       
      sx=(act_line%QCX(NStation)-act_line%COR(1))/S
      sy=(act_line%QCY(NStation)-act_line%COR(2))/S
      sz=(act_line%QCZ(NStation)-act_line%COR(3))/S
    
      do istation=1,Nstation
         tx=act_line%tx(istation)
         ty=act_line%ty(istation)
         tz=act_line%tz(istation)

         call QuatRot(tx,ty,tz,pitch_angle,sx,sy,sz,0.0_mytype,0.0_mytype,0.0_mytype,new_tx,new_ty,new_tz)
         tmag=sqrt(new_tx**2.+new_ty**2.+new_tz**2.)
         act_line%tx(istation)=new_tx/tmag
         act_line%ty(istation)=new_ty/tmag
         act_line%tz(istation)=new_tz/tmag

         act_line%pitch(istation)=act_line%pitch(istation)-pitch_angle
      enddo

      call make_actuatorline_geometry(act_line)

      return

    end subroutine pitch_actuator_line

    !*******************************************************************************
    !
    subroutine populate_blade_airfoils(NElem,NData,EAirfoil,AirfoilData,ETtoC)
    !
    !*******************************************************************************

      use param, only: zero
      use dbg_schemes, only: abs_prec

      implicit none
      integer, intent(in) :: NElem,NData
      type(AirfoilType), dimension(NElem) :: EAirfoil
      type(AirfoilType), dimension(NData) :: AirfoilData
      real(mytype), dimension(NElem), intent(in) :: ETtoC
      real(mytype), dimension(NData) :: Thicks, diffthicks
      integer :: ielem,idata,imin,imax,iint

      !==========================================
      ! This routine initialises the airfoil struct for the blades by interpolating from the data
      ! We need to interpolate from two or more
      !==========================================
    
      do ielem=1,NElem
         EAirfoil(ielem)%afname = int2str(ielem)
         call allocate_airfoil(EAirfoil(ielem),MaxAOAVals,MaxReVals)
    
         Thicks(:)=zero
         diffthicks(:)=zero
         imin=0
         imax=0
         iint=0
    
         do idata=1,NData
            Thicks(idata)=AirfoilData(idata)%tc
            diffthicks(idata)=abs_prec(AirfoilData(idata)%tc-ETtoC(ielem))
         enddo
         imin=minloc(Thicks,1)
         imax=maxloc(Thicks,1)
         iint=minloc(diffthicks,1)
        
         if (ETtoC(ielem)>=Thicks(imax)) then
            call copy_airfoil_values(EAirfoil(ielem),AirfoilData(imax))
         else if (ETtoC(ielem)<=Thicks(imin)) then
            call copy_airfoil_values(EAirfoil(ielem),AirfoilData(imin))
         else
            call copy_airfoil_values(EAirfoil(ielem),AirfoilData(iint))
         endif
      enddo

    end subroutine populate_blade_airfoils
    
    !*******************************************************************************
    !
    subroutine rotate_actuatorline(actuatorline,origin,rotN,theta)
    !
    !*******************************************************************************

      use dbg_schemes, only: sqrt_prec
      use param, only : zero

      implicit none
      type(ActuatorLineType), intent(inout) :: actuatorline
      real(mytype), intent(in) :: origin(3), rotN(3)
      real(mytype) :: theta,nrx,nry,nrz,px,py,pz
      real(mytype) :: vrx,vry,vrz,VMag
      real(mytype) :: xtmp,ytmp,ztmp, txtmp, tytmp, tztmp
      integer :: ielem
      
      ! Specify the rotation axis and the normal vector of rotation 
      nrx=rotN(1)
      nry=rotN(2)
      nrz=rotN(3)
    
      px=origin(1)
      py=origin(2)
      pz=origin(3)

      do ielem=1,actuatorline%Nelem+1
         ! Blade end locations (quarter chord, stations)
         xtmp=actuatorline%QCx(ielem)
         ytmp=actuatorline%QCy(ielem)
         ztmp=actuatorline%QCz(ielem)
    
         call QuatRot(xtmp,ytmp,ztmp,theta,nrx,nry,nrz,px,py,pz,vrx,vry,vrz)
         actuatorline%QCx(ielem)=vrx
         actuatorline%QCy(ielem)=vry
         actuatorline%QCz(ielem)=vrz
    
         ! Tangent vectors
         txtmp=actuatorline%tx(ielem)
         tytmp=actuatorline%ty(ielem)
         tztmp=actuatorline%tz(ielem)
    
         call QuatRot(txtmp,tytmp,tztmp,theta,nrx,nry,nrz,zero,zero,zero,vrx,vry,vrz)
         VMag=sqrt_prec(vrx**2+vry**2+vrz**2)
         actuatorline%tx(ielem)=vrx/VMag
         actuatorline%ty(ielem)=vry/VMag
         actuatorline%tz(ielem)=vrz/VMag
      enddo

    end subroutine rotate_actuatorline
       
    !*******************************************************************************
    !
    subroutine read_actuatorline_geometry(FN,Rmax,SpanwiseVec,rR,ctoR,pitch,thick,Nstations)
    !
    !*******************************************************************************
    
      implicit none
      character(len=100), intent(in) :: FN ! Filename of the geometry file
      real(mytype), dimension(3) :: SpanwiseVec
      real(mytype), allocatable, intent(out) :: rR(:),ctoR(:),pitch(:),thick(:)
      real(mytype), intent(out) :: Rmax
      integer, intent(out) :: Nstations
      integer :: i
      character(1000) :: ReadLine
    
      open(15,file=FN)

      ! Read the Number of Blades
      read(15,'(A)') ReadLine
      read(ReadLine(index(ReadLine,':')+1:),*) Rmax
    
      ! Read Spanwise actuator line axis
      read(15,'(A)') ReadLine
      read(ReadLine(index(ReadLine,':')+1:),*) SpanwiseVec(1), SpanwiseVec(2), SpanwiseVec(3)
    
      ! Read the stations specs
      read(15,'(A)') ReadLine
      read(ReadLine(index(ReadLine,':')+1:),*) Nstations

      allocate(rR(Nstations),ctoR(Nstations),pitch(Nstations),thick(Nstations))
      do i=1,NStations
         read(15,'(A)') ReadLine
         read(ReadLine,*) rR(i), ctoR(i), pitch(i), thick(i)
      enddo
    
      close(15)

    end subroutine read_actuatorline_geometry

    !*******************************************************************************
    !
    subroutine make_actuatorline_geometry(blade)
    !
    !*******************************************************************************

      use param, only: two
      use dbg_schemes, only: sqrt_prec

      implicit none
      type(ActuatorLineType), intent(inout) :: blade ! For simplity. In fact this is an actuator line.
      integer :: nbe, nej, j
      real(mytype) :: sEM, tEM, nEM
      real(mytype) :: sE(3),tE(3),normE(3),P1(3),P2(3),P3(3),P4(3),V1(3),V2(3),V3(3),V4(3),A1(3),A2(3)

      ! Calculates element geometry from element end (stations) geometry
      nbe=blade%NElem

      do j=1,nbe
         nej=1+j

         ! Element center locations
         blade%PEx(nej-1)=(blade%QCx(nej)+blade%QCx(nej-1))/two
         blade%PEy(nej-1)=(blade%QCy(nej)+blade%QCy(nej-1))/two
         blade%PEz(nej-1)=(blade%QCz(nej)+blade%QCz(nej-1))/two
         blade%ERdist(nej-1)=sqrt((blade%PEX(nej-1)-blade%COR(1))**2 +(blade%PEY(nej-1)-blade%COR(2))**2+(blade%PEZ(nej-1)-blade%COR(3))**2) ! Element center length

         ! Set spanwise and tangential vectors
         sE=(/blade%QCx(nej)-blade%QCx(nej-1),blade%QCy(nej)-blade%QCy(nej-1),blade%QCz(nej)-blade%QCz(nej-1)/) ! nominal element spanwise direction set opposite to QC line
         sEM=sqrt(dot_product(sE,sE))
         blade%EDS(nej-1) = sEM
         sE=sE/sEM

         tE=(/blade%tx(nej)+blade%tx(nej-1),blade%ty(nej)+blade%ty(nej-1),blade%tz(nej)+blade%tz(nej-1)/)/2.0

         ! Force tE normal to sE
         tE=tE-dot_product(tE,sE)*sE
         tEM=sqrt(dot_product(tE,tE))
         tE=tE/tEM

         blade%sEx(nej-1)=sE(1)
         blade%sEy(nej-1)=sE(2)
         blade%sEz(nej-1)=sE(3)
         blade%tEx(nej-1)=tE(1)
         blade%tEy(nej-1)=tE(2)
         blade%tEz(nej-1)=tE(3)

         ! Calc normal vector
         call cross(sE(1),sE(2),sE(3),tE(1),tE(2),tE(3),normE(1),normE(2),normE(3))
         nEM=sqrt(dot_product(normE,normE))
         normE=normE/nEM
         blade%nEx(nej-1)=normE(1)
         blade%nEy(nej-1)=normE(2)
         blade%nEz(nej-1)=normE(3)

         if (blade%FlipN) then
            blade%nEx(nej-1)=-blade%nEx(nej-1)
            blade%nEy(nej-1)=-blade%nEy(nej-1)
            blade%nEz(nej-1)=-blade%nEz(nej-1)
            blade%sEx(nej-1)=-blade%sEx(nej-1)
            blade%sEy(nej-1)=-blade%sEy(nej-1)
            blade%sEz(nej-1)=-blade%sEz(nej-1)
         endif

         ! Calc element area and chord
         P1=(/blade%QCx(nej-1)-0.25*blade%C(nej-1)*blade%tx(nej-1),blade%QCy(nej-1)-0.25*blade%C(nej-1)*blade%ty(nej-1),blade%QCz(nej-1)-0.25*blade%C(nej-1)*blade%tz(nej-1)/)
         P2=(/blade%QCx(nej-1)+0.75*blade%C(nej-1)*blade%tx(nej-1),blade%QCy(nej-1)+0.75*blade%C(nej-1)*blade%ty(nej-1),blade%QCz(nej-1)+0.75*blade%C(nej-1)*blade%tz(nej-1)/)
         P3=(/blade%QCx(nej)+0.75*blade%C(nej)*blade%tx(nej),blade%QCy(nej)+0.75*blade%C(nej)*blade%ty(nej),blade%QCz(nej)+0.75*blade%C(nej)*blade%tz(nej)/)
         P4=(/blade%QCx(nej)-0.25*blade%C(nej)*blade%tx(nej),blade%QCy(nej)-0.25*blade%C(nej)*blade%ty(nej),blade%QCz(nej)-0.25*blade%C(nej)*blade%tz(nej)/)

         V1=P2-P1
         V2=P3-P2
         V3=P4-P3
         V4=P1-P4

         ! Calc quad area from two triangular facets
         call cross(V1(1),V1(2),V1(3),V2(1),V2(2),V2(3),A1(1),A1(2),A1(3))
         A1=A1/2.0
         call cross(V3(1),V3(2),V3(3),V4(1),V4(2),V4(3),A2(1),A2(2),A2(3))
         A2=A2/2.0
         blade%EArea(nej-1)=sqrt(dot_product(A1,A1))+sqrt(dot_product(A2,A2))

         ! Calc average element chord from area and span
         blade%EC(nej-1)=blade%EArea(nej-1)/sEM
         blade%ETtoC(nej-1)=0.5*(blade%thick(nej)+blade%thick(nej-1))
         blade%Epitch(nej-1)=0.5*(blade%pitch(nej)+blade%pitch(nej-1))   
      enddo

      return

    end subroutine make_actuatorline_geometry
    
    !*******************************************************************************
    !
    subroutine allocate_actuatorline(actuatorline,NStations)
    !
    !*******************************************************************************

      implicit none
    
      type(ActuatorLineType) :: actuatorline
      integer, intent(in) :: Nstations
      integer :: NElem
    
      Nelem=Nstations-1
      actuatorline%Nelem = Nelem

      allocate(actuatorline%QCx(NElem+1))
      allocate(actuatorline%QCy(NElem+1))
      allocate(actuatorline%QCz(NElem+1))
      allocate(actuatorline%tx(NElem+1))
      allocate(actuatorline%ty(NElem+1))
      allocate(actuatorline%tz(NElem+1))
      allocate(actuatorline%C(NElem+1))
      allocate(actuatorline%thick(NElem+1))
      allocate(actuatorline%pitch(NElem+1))
      allocate(actuatorline%PEx(NElem))
      allocate(actuatorline%PEy(NElem))
      allocate(actuatorline%PEz(NElem))
      allocate(actuatorline%tEx(NElem))
      allocate(actuatorline%tEy(NElem))
      allocate(actuatorline%tEz(NElem))
      allocate(actuatorline%nEx(NElem))
      allocate(actuatorline%nEy(NElem))
      allocate(actuatorline%nEz(NElem))
      allocate(actuatorline%sEx(NElem))
      allocate(actuatorline%sEy(NElem))
      allocate(actuatorline%sEz(NElem))
      allocate(actuatorline%EC(NElem))
      allocate(actuatorline%EDS(NElem))
      allocate(actuatorline%EArea(NElem))
      allocate(actuatorline%ETtoC(NElem))
      allocate(actuatorline%Eepsilon(NElem))
      allocate(actuatorline%EAirfoil(Nelem))
      allocate(actuatorline%EDynstall(Nelem))
      allocate(actuatorline%ELBstall(Nelem))
      allocate(actuatorline%ERdist(Nelem))
      allocate(actuatorline%EVx(NElem))
      allocate(actuatorline%EVy(NElem))
      allocate(actuatorline%EVz(NElem))
      allocate(actuatorline%EVbx(NElem))
      allocate(actuatorline%EVby(NElem))
      allocate(actuatorline%EVbz(NElem))
      allocate(actuatorline%EObx(NElem))
      allocate(actuatorline%EOby(NElem))
      allocate(actuatorline%EObz(NElem))
      allocate(actuatorline%Epitch(Nelem))
      allocate(actuatorline%EAOA(Nelem))
      allocate(actuatorline%EAOAdot(Nelem))
      allocate(actuatorline%EUndot(Nelem))
      allocate(actuatorline%ERE(Nelem))
      allocate(actuatorline%EUr(Nelem))
      allocate(actuatorline%EUn(Nelem))
      allocate(actuatorline%EAOA_LAST(Nelem))
      allocate(actuatorline%EUn_LAST(Nelem))
      allocate(actuatorline%ECD(Nelem))
      allocate(actuatorline%ECL(Nelem))
      allocate(actuatorline%ECLcirc(Nelem))
      allocate(actuatorline%ECM(Nelem))
      allocate(actuatorline%EGamma(Nelem))
      allocate(actuatorline%ECN(Nelem))
      allocate(actuatorline%ECT(Nelem))
      allocate(actuatorline%EFn(NElem))
      allocate(actuatorline%EFt(NElem))
      allocate(actuatorline%EMS(NElem))
      allocate(actuatorline%EFx(NElem))
      allocate(actuatorline%EFy(NElem))
      allocate(actuatorline%EFz(NElem))
      allocate(actuatorline%A_rbf(NElem,NElem))
      allocate(actuatorline%EEndeffects_factor(NElem))

      actuatorline%EAOA=0.0_mytype

    end subroutine allocate_actuatorline

end module actuator_line_element
