module dynstall
    
    use decomp_2d, only: nrank, mytype
    use airfoils

    implicit none
    real(mytype), parameter :: SMALL= TINY(1.0) ! The number slighty greater than unity in single precision.
    
type DS_Type 
    ! These parameters follow the formulation of turbinesFOAM
    logical :: StallFlag = .false.
    logical :: StallFlag_prev = .false. 
    logical :: do_calcAlphaEquiv=.false.
    real(mytype) :: deltaS,time,time_prev
    real(mytype) :: X,X_prev,Y,Y_prev,Z,Z_prev
    real(mytype) :: Mach 
    real(mytype) :: T1,T2,T3
    real(mytype) :: Tp, Tf, Tv, Tvl ! Time periods
    real(mytype) :: A1,A2,A3,b1,b2,b3,K1,K2,r,r0,alphaDS0DiffDeg
    real(mytype) :: alphaEquiv, etaL, etaL_prev
    real(mytype) :: D,D_prev,DF,DF_prev
    real(mytype) :: CV,CV_prev,CNV,CNV_prev
    real(mytype) :: fprime,fprime_prev,fDoublePrime
    real(mytype) :: eta,S1,S2,CD0
    real(mytype) :: tau, tau_prev
    integer :: nNewTimes
    real(mytype) :: f ! not needed
    real(mytype) :: fcrit
    real(mytype) :: CMC
    real(mytype) :: CNC, CNAlpha
    real(mytype) :: lambdaL,lambdaL_prev
    real(mytype) :: lambdaM,lambdaM_prev
    real(mytype) :: TI,E0, AlphaZeroLift
    real(mytype) :: speed_of_sound
    real(mytype) :: H,H_prev,J,J_prev
    real(mytype) :: CNI,CMI,CNP,CNP_prev,CN1
    real(mytype) :: CN_prev,CNprime
    real(mytype) :: DP,DP_prev
    real(mytype) :: AlphaPrime,AlphaPrime_prev
    real(mytype) :: DAlpha, DAlpha_prev, TAlpha
    real(mytype) :: Alpha, Alpha_prev, deltaAlpha,deltaAlpha_prev, AlphaCrit
    real(mytype) :: alpha1, alphaDS0, alphaSS, dalphaDS
    real(mytype) :: Re
    real(mytype) :: Vx,CNF,CT,CN,cmFitExponent,K0,CM
    real(mytype), allocatable :: ReList(:), CNAlphaList(:), CD0List(:), CN1List(:), alpha1List(:), S1List(:), S2List(:)
    real(mytype), allocatable :: alphaSSList(:), K1List(:), K2List(:)
end type DS_Type

contains

    !*******************************************************************************
    !
    subroutine dystl_init(ds,dynstallfile)
    !
    !*******************************************************************************
        
      implicit none
      type(DS_Type), intent(inout) :: ds
      character(len=100), intent(in) :: dynstallfile
      character(1000) :: ReadLine
      integer :: i, Nstall
      real(mytype) :: Tp, Tf, TAlpha, alphaDS0DiffDeg, r0, Tv, Tvl, B1, B2, eta, E0, CNAlpha,CD0,alpha1,S1,S2,alphaSS,K1,K2 
      NAMELIST/DynstallParam/Tp, Tf, TAlpha, alphaDS0DiffDeg, r0, Tv, Tvl, B1, B2, eta, E0, CNAlpha, CD0,alpha1,S1,S2,alphaSS,K1,K2

      ! Default values for the parameters
      Tp=1.7
      Tf=3.0
      TAlpha=6.25
      alphaDS0DiffDeg=3.8
      r0=0.01
      Tv=11.0
      Tvl=8.7
      B1=0.5
      B2=0.2
      eta=0.98
      E0=0.16
      CNAlpha=6.498  
      CD0=0.01   
      alpha1=10.92  
      S1=0.0483 
      S2=0.0894 
      alphaSS=12.49 
      K1=2.77e-8 
      K2=-1.43e-6

      ! Read from file
      open(30,file=dynstallfile) 
      read(30,nml=DynstallParam)
      close(30)
      ds%Tp=Tp
      ds%Tf=Tf
      ds%TAlpha=TAlpha
      ds%alphaDS0DiffDeg=alphaDS0DiffDeg
      ds%r0=r0
      ds%Tv=Tv
      ds%Tvl=Tvl
      ds%B1=B1
      ds%B2=B2
      ds%eta=eta
      ds%E0=E0
      ds%CNAlpha=CNAlpha  
      ds%CD0=CD0   
      ds%Alpha1=Alpha1  
      ds%S1=S1 
      ds%S2=S2 
      ds%AlphaSS=alphaSS 
      ds%K1=K1
      ds%K2=K2
            
      ! Convert degrees to radians
      ds%alpha1=ds%alpha1*conrad
      ds%alphaSS=ds%alphaSS*conrad

      ! Define CN1 using fcrit=0.7 (Hardcoded)
      ds%CN1=ds%CNAlpha*ds%alpha1*(1.+sqrt(0.7)/2.0)**2.

      ! Here we need to compute S1, S2, K1, K2
        
      ! Set all other initial conditions
      ds%X=0.0
      ds%X_prev=0.0
      ds%Y=0.0
      ds%Y_prev=0.0
      ds%Z=0.0
      ds%Z_prev=0.0
      ds%Alpha=0.0
      ds%Alpha_prev=0.0
      ds%deltaAlpha=0.0
      ds%deltaAlpha_prev=0.0
      ds%time_prev=0.0
      ds%D=0.0
      ds%D_prev=0.0
      ds%DP=0.0
      ds%DP_prev=0.0
      ds%CNP=0.0
      ds%CNP_prev=0.0
      ds%fprime=1.0
      ds%fprime_prev=1.0
      ds%DF=0.0
      ds%DF_prev=0.0
      ds%CV=0.0
      ds%CV_prev=0.0
      ds%CNV=0.0
      ds%CNV_prev=0.0
      ds%eta=0.95
      ds%tau=0.0
      ds%tau_prev=0.0
      ds%nNewTimes=0
      ds%K0=1e-6
      ds%K1=0.0
      ds%K2=0.0
      ds%cmFitExponent=2
      ds%CM=0.0
      ds%Re=0.0
      ds%speed_of_sound=343.0
      ds%etaL=0.0
      ds%etaL_prev=0.0
      ds%A1=0.165
      ds%A2=0.335
      ds%A3=0.5
      ds%T1=20.0
      ds%T2=4.5
      ds%H=0.0
      ds%H_prev=0.0
      ds%lambdaL=0.0
      ds%lambdaL_prev=0.0
      ds%J=0.0
      ds%J_prev=0.0
      ds%lambdaM=0.0
      ds%lambdaM_prev=0.0
      ds%CMC=0.0
      ds%CMI=0.0
      ds%fcrit=0.7
      ds%AlphaPrime=0.0
      ds%AlphaCrit=17.0
      ds%r=0.0
      ds%DAlpha=0.0
      ds%DAlpha_prev=0.0
 
      return

    end subroutine dystl_init
   
    !*******************************************************************************
    !
    subroutine calcAlphaEquiv(ds)
    !
    !*******************************************************************************
        
      implicit none
      type(DS_Type), intent(inout) :: ds
      real(mytype) :: beta

      ! Calculates the equivalent angle of attack after having applied the defeciency functions
      ds%T3=1.25*ds%Mach
      beta=1-ds%Mach**2
        
      ds%X=ds%X_prev*exp(-beta*ds%deltaS/ds%T1)+ds%A1*(ds%etaL-ds%etaL_prev)*exp(-beta*ds%deltaS/(2.0*ds%T1))
      ds%Y=ds%Y_prev*exp(-beta*ds%deltaS/ds%T2)+ds%A2*(ds%etaL-ds%etaL_prev)*exp(-beta*ds%deltaS/(2.0*ds%T2))
      ds%Z=ds%Z_prev*exp(-beta*ds%deltaS/ds%T3)+ds%A3*(ds%etaL-ds%etaL_prev)*exp(-beta*ds%deltaS/(2.0*ds%T3))

      ds%alphaEquiv = ds%alpha-ds%X-ds%Y-ds%Z

      if (abs(ds%alphaEquiv)>2*pi) then
         ds%alphaEquiv=mod(ds%alphaEquiv,2*pi)
      endif

    end subroutine calcAlphaEquiv
    
    !*******************************************************************************
    !
    subroutine DynstallCorrect(dynstall,airfoil,time,dt,Urel,chord,alpha,Re,CLdyn,CDdyn,CM25dyn)
    !
    !*******************************************************************************
        
      implicit none
      type(DS_Type), intent(inout) :: dynstall
      type(AirfoilType), intent(in) :: airfoil
      real(mytype), intent(in) :: time,dt,Urel,chord,Re
      real(mytype), intent(inout) :: alpha
      real(mytype), intent(out) :: CLdyn,CDdyn,CM25dyn
      real(mytype) :: mach

      ! update previous values if time has changed
      dynstall%nNewTimes=dynstall%nNewTimes+1

      ! Set previous angle equal to the current one if first time step with dynstall 
      if (dynstall%nNewTimes <=1) then
         dynstall%alpha_Prev=alpha
      endif
    
      dynstall%Alpha=alpha
      dynstall%mach=urel/dynstall%speed_of_sound 
      dynstall%Re=Re
      dynstall%deltaAlpha=dynstall%Alpha-dynstall%alpha_Prev
      dynstall%deltaS=2.0*Urel*dt/chord
 
      dynstall%alphaEquiv = dynstall%Alpha
      dynstall%AlphaZeroLift=airfoil%alzer*conrad ! Make it in rad
      ! Ok
      !if (nrank==0) write(*,*) 180*dynstall%Alpha/pi, dynstall%mach, dynstall%Re, dynstall%deltaAlpha, dynstall%deltaS
      ! Ok

      ! Evaluate static coefficient data if if has changed from the Reynolds number correction
        
      !call EvalStaticData(dynstall)
        
      call calcUnsteady(dynstall,chord,Urel,dt)
        
      call calcSeparated(dynstall)
    
      ! Modify Coefficients
      CLdyn=dynstall%CN*cos(dynstall%alpha)+dynstall%CT*sin(dynstall%alpha)
      CDdyn=dynstall%CN*sin(dynstall%alpha)-dynstall%CT*cos(dynstall%alpha)+dynstall%CD0
      CM25dyn=dynstall%CM
          
      call update_DynStall(dynstall,time)
      
      return

    end subroutine DynstallCorrect
    
    !*******************************************************************************
    !
    subroutine EvalStaticData(ds)
    !
    !*******************************************************************************

      implicit none
      type(DS_Type), intent(inout) :: ds
      integer :: ilist,iup,ilow,imin
      real(mytype) ::minRediff, Rediff

      minRediff=1e12
      do ilist=1,size(ds%Relist)
         Rediff=abs(ds%ReList(ilist)-ds%Re)
         if (Rediff<minRediff) then
            imin=ilist
            minRediff=Rediff
         endif
      enddo

      if (ds%Re<=ds%ReList(1)) then
         iup=1
         ilow=1
      else if (ds%Re>=ds%ReList(size(ds%Relist))) then
         iup=size(ds%Relist)
         ilow=size(ds%Relist)
      else if (ds%Re>ds%ReList(imin).and.ds%Re<ds%ReList(size(ds%ReList))) then
         ilow=imin
         iup=imin+1
      else if (ds%Re<ds%ReList(imin).and.ds%Re>ds%ReList(1)) then
         ilow=imin-1
         iup=imin
      endif 
        
      ! Define the Static Coefficients
      if (iup==ilow) then
         ds%CNAlpha=ds%CNAlphaList(ilow)
         ds%CD0=ds%CD0List(ilow)
         ds%CN1=ds%CN1List(ilow)
         ds%alpha1=ds%Alpha1List(ilow)
         ds%S1=ds%S1List(ilow)
         ds%S2=ds%S2List(ilow)
         ds%alphaSS=ds%alphaSSList(ilow)
         ds%K1=ds%K1List(ilow)
         ds%K2=ds%K2List(ilow)
      else
         ds%CNAlpha=ds%CNAlphalist(ilow)+(ds%Re-ds%ReList(ilow))*(ds%CNAlphalist(iup)-ds%CNAlphaList(ilow))/(ds%ReList(iup)-ds%ReList(ilow))
         ds%CD0=ds%CD0List(ilow)+(ds%Re-ds%ReList(ilow))*(ds%CD0List(iup)-ds%CD0List(ilow))/(ds%ReList(iup)-ds%ReList(ilow))
         ds%alpha1=ds%alpha1list(ilow)+(ds%Re-ds%ReList(ilow))*(ds%alpha1list(iup)-ds%alpha1List(ilow))/(ds%ReList(iup)-ds%ReList(ilow))
         ds%S1=ds%S1List(ilow)+(ds%Re-ds%ReList(ilow))*(ds%S1list(iup)-ds%S1List(ilow))/(ds%ReList(iup)-ds%ReList(ilow))
         ds%S2=ds%S2List(ilow)+(ds%Re-ds%ReList(ilow))*(ds%S2List(iup)-ds%S2List(ilow))/(ds%ReList(iup)-ds%ReList(ilow))
         ds%alphaSS=ds%alphaSSList(ilow)+(ds%Re-ds%ReList(ilow))*(ds%alphaSSlist(iup)-ds%alphaSSList(ilow))/(ds%ReList(iup)-ds%ReList(ilow))
         ds%K1=ds%K1List(ilow)+(ds%Re-ds%ReList(ilow))*(ds%K1list(iup)-ds%K1List(ilow))/(ds%ReList(iup)-ds%ReList(ilow))
         ds%K2=ds%K2List(ilow)+(ds%Re-ds%ReList(ilow))*(ds%K2list(iup)-ds%K2List(ilow))/(ds%ReList(iup)-ds%ReList(ilow))
      endif

      ds%alpha1=ds%alpha1*conrad
      ds%alphaSS=ds%alphaSS*conrad
      ! Define CN1 using fcrit=0.7 (Hardcoded)
      ds%CN1=ds%CNAlpha*ds%alpha1*(1.+sqrt(0.7)/2.0)**2.

      ! Here we need to compute S1, S2, K1, K2
       
      return

    end subroutine EvalStaticData

    !*******************************************************************************
    !
    subroutine calcUnsteady(ds,chord,Ur,dt)
    !
    !*******************************************************************************

      implicit none
      type(DS_Type), intent(inout) :: ds
      real(mytype), intent(in) :: chord,Ur,dt
      real(mytype) :: kAlpha, dAlphaDS
    
      ! Calculate the circulatory normal force coefficient
      ds%CNC=ds%CNAlpha*(ds%alphaEquiv-ds%AlphaZeroLift) ! Here CNAlpha is the Normal Force/ alpha (rad) slope

      ! Calculate the impulsive normal force coefficient
      ds%lambdaL=(pi/4.0)*(ds%alpha+chord/(4.0*Ur)*ds%deltaAlpha/dt)

      ds%TI=chord/ds%speed_of_sound*(1.0+3.0*ds%mach)/4.0

      ds%H=ds%H_prev*exp(-dt/ds%TI)+(ds%lambdaL-ds%lambdaL_prev)*exp(-dt/(2.0*ds%TI))
        
      ds%CNI=4.0/ds%mach*ds%H
    
      ! Calculate the impulsive moment coefficient
      ds%lambdaM=3*pi/16.0*(ds%alpha+chord/(4.0*Ur)*ds%deltaAlpha/dt)+pi/16.*chord/Ur*ds%deltaAlpha/dt
        
      ds%J=ds%J_prev*exp(-dt/ds%TI)+(ds%lambdaM-ds%lambdaM_prev)*exp(-dt/(2.0*ds%TI))
        
      ds%CMI=-4.0/ds%mach*ds%J

      ! Calculate total normal force and pitching moment coefficient
      ds%CNP = ds%CNC + ds%CNI

      ! Apply first-order lag to normal force coefficient
      ds%DP=ds%DP_prev*exp(-ds%deltaS/ds%Tp)+(ds%CNP-ds%CNP_prev)*exp(-ds%deltaS/(2.0*ds%Tp))
        
      ds%CNprime=ds%CNP-ds%DP

      ! Calculate lagged angle of attack
      ds%DAlpha=ds%DAlpha_prev*exp(-ds%deltaS/ds%TAlpha)+(ds%alpha-ds%alpha_prev)*exp(-ds%deltaS/(2.0*ds%TAlpha))
        
      ds%AlphaPrime=ds%alpha-ds%DAlpha
        
      ! Calculate reduced pitch rate
      ds%r=ds%deltaAlpha/dt*chord/(2.*Ur)

      ! Calculate alphaDS0
      dAlphaDS=ds%alphaDS0DiffDeg*conrad
      ds%alphaDS0=ds%alphaSS + dAlphaDS
        
      if (abs(ds%r)>=ds%r0) then
         ds%alphaCrit=ds%alphaDS0
      else
         ds%alphaCrit=ds%alphaSS+(ds%alphaDS0-ds%alphaSS)*abs(ds%r)/ds%r0
      endif

      if (abs(ds%AlphaPrime)>ds%AlphaCrit) then
         ds%StallFlag=.true.
         if (nrank==0) write(*,*) 'Section is Stalled'
      endif

      return

    end subroutine calcUnsteady

    !*******************************************************************************
    !
    subroutine calcSeparated(ds)
    !
    !*******************************************************************************

      implicit none
      type(DS_type), intent(inout) :: ds
      real(mytype) :: f, Tf,Tv, Tst, KN, m, cmf, cpv,cmv

      ! Calculate trailing-edge separation point
      if (abs(ds%alphaPrime)<ds%alpha1) then
         ds%fprime=1.0-0.4*exp((abs(ds%alphaPrime)-ds%alpha1)/ds%S1) 
      else
         ds%fprime=0.02+0.58*exp((ds%alpha1-abs(ds%alphaPrime))/ds%S2)
      endif

      ! Calculate vortex tracking time
      if (ds%StallFlag_prev.eqv..false.) then
         ds%tau=0.0
      else
         if (ds%tau==ds%tau_prev) then
            ds%tau=ds%tau_prev+ds%deltaS
         endif
      endif
        
      if (ds%tau>ds%Tf) then
         Tf=ds%Tf/2.
      else
         Tf=ds%Tf
      endif
        
      ! Calculate dynamic separation point
      ds%DF=ds%DF_prev*exp(-ds%deltaS/Tf)+(ds%fprime-ds%fprime_prev)*exp(-ds%deltaS/(2.*Tf))
      ds%fDoublePrime=ds%fprime-ds%DF

      ! Calculate vortex modulation parameter
      if (ds%tau>=0.0.and.ds%tau<=ds%Tv) then
         ds%Vx=sin(pi*ds%tau/(2.0*ds%Tv))**1.5
      else if (ds%tau>ds%Tv) then
         ds%Vx=cos(pi*(ds%tau-ds%Tv)/ds%Tvl)**2.0
      endif
        
      if (ds%tau>ds%Tvl) then
         ds%Vx=0.0
      endif

      ! Calculate normal force coefficient including dynamic separation point
      ds%CNF=ds%CNAlpha*(ds%alphaEquiv-ds%AlphaZeroLift)*((1.0+sqrt(ds%fDoublePrime))/2.0)**2 !+ds%CNI

      ! Calculate tangential force coefficient
      ds%CT=ds%eta*ds%CNAlpha*(ds%alphaEquiv-ds%AlphaZeroLift)**2.*(sqrt(ds%fDoublePrime)-ds%E0)
        
      ! Calculate static trailing-edge separation point
      if (abs(ds%alpha)<abs(ds%alpha1)) then
         f=1.0-0.4*exp((abs(ds%alpha)-ds%alpha1)/ds%S1)
      else
         f=0.02+0.58*exp((ds%alpha1-abs(ds%alpha))/ds%S2)
      endif

      ! Calculate vortex lift contribution
      ds%CNV=ds%B1*(ds%fDoublePrime-f)*ds%VX

      ! Total normal force coefficient is the combination of that from 
      ! circulatory effects, impulsive effects, dynamic separation, and vortex lift
      ds%CN=ds%CNF+ds%CNV
      
      ! Calculate moment coefficient
      m=ds%cmFitExponent
      cmf=(ds%K0+ds%K1*(1-ds%fDoublePrime)+ds%K2*sin(pi*ds%fDoublePrime**2))*ds%CNC
      ! + moment coefficient at Zero lift angle of attack
      cmv = ds%B2*(1.0-cos(pi*ds%tau/ds%Tv))*ds%CNV
      ds%CM=cmf+cmv+ds%CMI
        
      return

    end subroutine calcSeparated
    
    !*******************************************************************************
    !
    subroutine update_dynstall(ds,time)
    !
    !*******************************************************************************

      implicit none
      type(DS_type), intent(inout) ::ds
      real(mytype), intent(in) :: time

      ds%time_prev=           time
      ds%Alpha_prev=          ds%Alpha
      ds%X_prev=              ds%X
      ds%Y_prev=              ds%Y
      ds%Z_prev=              ds%Z
      ds%deltaAlpha_prev =    ds%deltaAlpha
      ds%D_prev=              ds%D
      ds%DP_prev=             ds%DP
      ds%CNP_prev=            ds%CNP
      ds%DF_prev=             ds%DF
      ds%fprime_prev=         ds%fprime
      ds%CV_prev=             ds%CV
      ds%CNV_prev=            ds%CNV
      ds%StallFlag_prev=      ds%StallFlag
      ds%tau_prev=            ds%tau
      ds%etaL_prev=           ds%etaL
      ds%H_prev=              ds%H
      ds%lambdaL_prev=        ds%lambdaL
      ds%J_prev=              ds%J
      ds%lambdaM_prev=        ds%lambdaM
      ds%AlphaPrime_prev=     ds%AlphaPrime
      ds%DAlpha_prev=         ds%DAlpha

      return

    end subroutine update_dynstall

end module
