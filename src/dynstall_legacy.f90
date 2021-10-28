module dynstall_legacy

    use airfoils

    implicit none
       
type LB_Type
    logical :: StallFlag = .false.
    real(mytype) :: dp
    real(mytype) :: dF
    real(mytype) :: dCNv
    real(mytype) :: sLEv
    integer :: LESepState
    real(mytype) :: CLRef 
    real(mytype) :: CLRefLE
    real(mytype) :: CLA
    real(mytype) :: CLCritP
    real(mytype) :: CLCritN
    integer :: CLRateFlag
    real(mytype) :: Fstat
    real(mytype) :: F
    real(mytype) :: cv
    real(mytype) :: dcv
    real(mytype) :: CLRef_Last
    real(mytype) :: CLRefLE_Last
    real(mytype) :: Fstat_Last
    real(mytype) :: cv_Last
    ! Additional LB diagnostic output
    integer, dimension(9) :: LB_LogicOutputs
    integer :: Logic_W(9)=[1,1,1,1,1,3,2,1,1]    
end type LB_Type

    Type(LB_Type) :: lb

contains

    !*******************************************************************************
    !
    subroutine dystl_init_LB(lb,dynstallfile)
    !
    !*******************************************************************************
        
      use param, only: zero, one

      implicit none
      type(LB_Type) :: lb
      character :: dynstallfile*80
      real(mytype) :: CLcritp, CLcritn, CLalpha
      NAMELIST/LBParam/CLcritp,CLcritn,CLalpha

      lb%StallFlag = .true.
      lb%dp=zero
      lb%dF=zero
      lb%dCNv=zero
      lb%LESepState=zero
      lb%sLEv=zero
      lb%CLRef_Last=zero
      lb%CLRefLE_Last=zero
      lb%Fstat_Last=one
      lb%cv_Last=zero
      lb%LB_LogicOutputs(:)=0
        
      open(30,file=dynstallfile) 
      read(30,nml=LBParam)
      close(30)

      lb%CLCritP=CLcritp
      lb%CLCritN=CLcritn
      lb%CLA=CLalpha

      return

    end subroutine dystl_init_LB

    !*******************************************************************************
    !
    subroutine LB_EvalIdealCL(AOA,AOA0,CLa,RefFlag,CLID)
    !
    !*******************************************************************************

      use decomp_2d, only: mytype
      use param, only: one, two, thirty
      use dbg_schemes, only: sin_prec
    
      implicit none
      real(mytype) :: AOA, AOA0, CLa
      integer :: RefFlag
      real(mytype) :: CLID
      real(mytype) :: IDS, aID, d1, CLaI, aIDc

      ! AOA inputs in radians
      ! AOA0 is zero lift AOA
      ! RefFlag defines whether to output reference CL or ideal CL
      ! CLa is reference lift slope (per radian) to be used for reference CL (ideal CLa is 2*pi)

      aID=AOA-AOA0
      call Force180(aID)
      ! reflect function across axis for abs(aID)>90
      IDS=1
      if (aID>pi/two) then
         aID=pi-aID
         IDS=-one
      else if (aID<-pi/two) then
         aID=-pi-aID
         IDS=-one
      endif

      ! If RefFlag is 1, output reference CL, otherwise round off ideal CL at high AOA
      if (RefFlag==1) then
         CLID=IDS*CLa*aID
      else
         ! round off the ideal CL after cutoff AOA
         aIDc=thirty*conrad
         d1=1.8_mytype
         CLaI=two*pi
         if (abs(aID)<aIDc) then
            CLID=IDS*CLaI*aID
         else
            CLID=IDS*(CLaI*(aIDc-one/d1*sin_prec(d1*aIDc))+CLaI/d1*sin_prec(d1*aID))
         endif
      endif
    
    end subroutine LB_EvalIdealCL

    !*******************************************************************************
    !
    subroutine Force180(a)
    !
    !*******************************************************************************
        
      use dbg_schemes, only: exp_prec
      use param, only: two
    
      implicit none
      real(mytype) :: a

      ! alpha in radians
      if (a>pi) then
         a=a-two*pi
      else if (a<-pi) then
         a=a+two*pi
      endif

    end subroutine Force180
    
    !*******************************************************************************
    !
    subroutine LB_UpdateStates(lb,airfoil,Re,ds)
    !
    !*******************************************************************************
        
      use decomp_2d, only: mytype
      use param, only: half, one, two, three, four, eleven
      use dbg_schemes, only: exp_prec

      implicit none
      type(LB_Type), intent(inout) :: lb
      type(AirfoilType), intent(in) :: airfoil
      real(mytype), intent(in) :: Re, ds
      integer :: i, nei, j, IsBE
      real(mytype) :: Tf,TfRef,Tp,TvRef, Tv, TvL

      ! Set model parameters. All of these are potentially a function of Mach
      ! and are set to low mach values...
      Tp=four                   ! time constant on LE pressure response to change in CL
      TfRef=three               ! time constant on TE separation point travel
      TvRef=6.3_mytype          ! time constant on LE vortex lift indicial function
      TvL=eleven                ! Characteristic LE vortex travel time

      ! Eval LE separation state
      ! Note: CLCrit is critical (ideal) CL value for LE separation. This is
      ! approximately equal to the CL that would exist at the angle of attack at
      ! max CL if the CL curve had remained linear
    
      if (lb%LESepState==0 .and. (lb%CLRefLE>lb%CLCritP .or. lb%CLRefLE<lb%CLCritN)) then
         ! In LE separation state
         lb%LESepState=1
         lb%sLEv=0 ! reset leading edge vortex time counter

         ! Set logic state flags (for model diagnosis output)
         lb%LB_LogicOutputs(5)=1
      else if (lb%LESepState==1 .and. (lb%CLRefLE<lb%CLCritP .and. lb%CLRefLE>lb%CLCritN)) then
         ! Out of LE separation state
         lb%LESepState=0
         lb%sLEv=0 ! reset leading edge vortex time counter

         ! Set logic state flags (for model diagnosis output)
         lb%LB_LogicOutputs(5)=2
      endif

      ! Set time constants based on LE separation state and TE separation point
      ! location. Different time constants for abs(CL) increasing vs. decreasing
      if (lb%LESepState==1) then
         if (lb%sLEv<TvL) then
            if (lb%CLRateFlag>0) then
               !Tf=3.0*TfRef ! original
               Tf=four*TfRef
               Tv=TvRef

               ! Set logic state flags (for model diagnosis output)
               lb%LB_LogicOutputs(8)=1
            else
               Tf=half*TfRef
               Tv=half*TvRef

               ! Set logic state flags (for model diagnosis output)
               lb%LB_LogicOutputs(8)=2
            endif

            ! Set logic state flags (for model diagnosis output)
            lb%LB_LogicOutputs(7)=1
         else if (lb%sLEv<two*TvL) then
            if (lb%CLRateFlag>0) then
               ! orig 
               !Tf=1.0/3.0*TfRef
               !Tv=1.0/4.0*TvRef
               Tf=two*TfRef
               Tv=TvRef

               ! Set logic state flags (for model diagnosis output)
               lb%LB_LogicOutputs(8)=3
            else
               Tf=half*TfRef
               Tv=half*TvRef

               ! Set logic state flags (for model diagnosis output)
               lb%LB_LogicOutputs(8)=4                                                        
            endif

            ! Set logic state flags (for model diagnosis output)
            lb%LB_LogicOutputs(7)=2                                                
         else
            ! orig
            !Tf=4.0*TfRef
            !Tv=0.9*TvRef
            Tf=TfRef
            Tv=TvRef

            ! Set logic state flags (for model diagnosis output)
            lb%LB_LogicOutputs(7)=3
         endif

        ! Set logic state flags (for model diagnosis output)
        lb%LB_LogicOutputs(6)=1
      else
         if (lb%F>0.7) then
            Tf=TfRef

            ! Set logic state flags (for model diagnosis output)
            lb%LB_LogicOutputs(7)=4
         else
            Tf=two*TfRef

            ! Set logic state flags (for model diagnosis output)
            lb%LB_LogicOutputs(7)=5
         endif
         Tv=TvRef

         ! Set logic state flags (for model diagnosis output)
         lb%LB_LogicOutputs(6)=2
      endif

      ! Update LE vortex time counter if in LE separation state
      if (lb%LESepState==1) then
         lb%sLEv=lb%sLEv+ds

         ! Set logic state flags (for model diagnosis output)
         lb%LB_LogicOutputs(9)=1
      endif

      ! Update states, first order lag equations, exponential recursion form (midpoint rule version)
      lb%dp=lb%dp*exp_prec(-ds/Tp)+(lb%CLRef-lb%CLRef_Last)*exp_prec(-ds/(two*Tp))
      lb%dF=lb%dF*exp_prec(-ds/Tf)+(lb%Fstat-lb%Fstat_Last)*exp_prec(-ds/(two*Tf))
      lb%dCNv=lb%dCNv*exp(-ds/Tv)+lb%dcv*exp(-ds/(two*Tv))

      ! Update lagged values
      lb%CLRef_Last=lb%CLRef
      lb%CLRefLE_Last=lb%CLRefLE
      lb%Fstat_Last=lb%Fstat 
      lb%cv_Last=lb%cv

    end subroutine LB_UpdateStates

    !*******************************************************************************
    !
    subroutine LB_LogicChecksum(LBCheck)
    !
    !*******************************************************************************
    
      implicit none
      integer :: LBCheck
      integer :: Loop

      ! Calculates a checksum for the logic states in the LB model
      LBCheck=0
      do Loop=1,9
         LBCheck=LBCheck+lb%LB_LogicOutputs(Loop)*lb%Logic_W(Loop)
      enddo

    end subroutine LB_LogicChecksum

    !*******************************************************************************
    !
    subroutine LB_DynStall(airfoil,lb,CLstat,CDstat,alphaL,alpha5,Re,CL,CD)
    !
    !*******************************************************************************

      use decomp_2d, only: mytype
      use param, only: zero, zpone, zptwofive, one, two, four, fifty
      use dbg_schemes, only: abs_prec, sqrt_prec, sin_prec, cos_prec, exp_prec

      implicit none
      type(AirfoilType) :: airfoil       ! Airfoil structure
      type(LB_Type) :: lb                ! Leishmann-Beddoes model structure
      real(mytype) :: CLstat, CDstat, alphaL, alpha5, Re, CL, CD
      real(mytype) :: AOA0, CLID, Trans, dCLRefLE, dAOARefLE, AOARefLE, CLstatF, C, C1, CLIDF 
      real(mytype) :: CLRatio, CLsep, CLF, dCDF, KD, CLa, NOF, dCLv, dCDv, acut, CLCritP, CLCritN

      !==========================================
      ! Routine that computes the Leishmann-Beddoes dynamic stall model
      ! with incompressible reduction and returns corrected values for 
      ! CL and CD having taken into account the dynamic stall effects
      !==========================================

      ! Airfoil data
      AOA0=airfoil%alzer
    
      ! Model constants
      KD=zpone          ! Trailing Edge separation drag factor

      ! Evaluate the ideal CL curve at current AOA
      call LB_EvalIdealCL(alphaL,AOA0,lb%CLa,1,lb%CLRef) 
      call LB_EvalIdealCL(alphaL,AOA0,lb%CLa,0,CLID)
    
      ! calc lagged ideal CL for comparison with critical LE separation CL
      Trans=(cos(alphaL-AOA0))**2 ! fair effect to zero at 90 deg. AOA...
      dCLRefLE=Trans*lb%dp ! dp is lagged CLRef change
      dAOARefLE=dCLRefLE/CLa

      ! define reference LE CL and AOA
      lb%CLRefLE=lb%CLRef-dCLRefLE
      if (lb%CLRefLE*(lb%CLRefLE-lb%CLRefLE_Last) > 0) then
         lb%CLRateFlag=1
      else
         lb%CLRateFlag=0
      endif
      AOARefLE=alphaL-dAOARefLE
      call Force180(AOARefLE)

      ! calc effective static TE separation point using effective LE AOA
      call EvalStaticCoeff(Re,AOARefLE*condeg,CLstatF,C,C1,airfoil)
      call LB_EvalIdealCL(AOARefLE,AOA0,CLa,0,CLIDF)
      if (abs_prec(CLIDF)<0.001_mytype) then
         CLRatio=999_mytype
      else
         CLRatio=CLstatF/CLIDF;
      endif

      if (CLRatio > zptwofive) then
         lb%Fstat=min((sqrt(four*CLRatio)-one)**2,one)

         ! Test logic
         lb%LB_LogicOutputs(1)=1
      else
         lb%Fstat=0

         ! Test logic
         lb%LB_LogicOutputs(1)=2
      endif
      ! calc lagged Fstat to represent dynamic TE separation point
      lb%F=lb%Fstat-lb%dF
      ! force limits on lagged F (needed due to discretization error...)
      lb%F=min(max(lb%F,zero),one)

      ! Calc dynamic CL due to TE separation as fairing between fully attached and fully separated predictions from the Kirchoff approximation at current AOA
      if (abs(CLID)<0.001_mytype) then
         CLRatio=999_mytype
      else
         CLRatio=CLstat/CLID
      endif

      if (CLRatio > one) then
         CLID=CLstat

         ! Test logic
         lb%LB_LogicOutputs(2)=1
      endif

      if (CLRatio > zptwofive) then
         CLsep=CLID/four

         ! Test logic
         lb%LB_LogicOutputs(3)=1
      else
         CLsep=CLstat

         ! Test logic
         lb%LB_LogicOutputs(3)=2
      endif
      CLF=CLsep+CLID*zptwofive*(lb%F+two*sqrt_prec(lb%F))
      dCDF=KD*(CLstat-CLF)*sign(one,CLstat)

      ! LE vortex lift component, dCNv is a lagged change in the added normal force due
      ! to LE vortex shedding. Assumed to affect lift coeff as an added circulation...
      dCLv=lb%dCNv*cos_prec(alpha5)
      dCDv=lb%dCNv*sin_prec(alpha5)
      ! vortex feed is given by the rate at which lift (circulation) is being shed due to dynamic separation. Lift component due to separation is defined by the
      ! difference between the ideal lift and the lift including dynamic separation effects.
      lb%cv=CLID-CLF
      lb%dcv=lb%cv-lb%cv_Last
      ! If the sign of dcv is opposite the reference LE CL, set to zero to disallow negative vorticity from shedding from the leading edge. Also, limit the model 
      ! at AOA>acut or if the magnitude of the reference CL is decreasing...
      acut=fifty*conrad
      if (sign(one,lb%dcv*lb%CLRefLE)<0 .OR. abs(alphaL-AOA0)>acut .OR. lb%CLRateFlag<0) then
         lb%dcv=zero

         ! Test logic
         lb%LB_LogicOutputs(4)=1
      endif

      ! Total lift and drag
      CL=CLF+dCLv
      CD=CDstat+dCDF+dCDv
    
      return

    end subroutine LB_DynStall

end module dynstall_legacy
