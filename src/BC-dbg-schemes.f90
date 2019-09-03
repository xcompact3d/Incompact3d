!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!!        FILE: BC-dbg-schemes.f90
!!!      AUTHOR: Felipe Schuch
!!!    MODIFIED: Paul Bartholomew
!!! DESCRIPTION: This module sets up debugging for schemes.
!!!   CHANGELOG: [2019-02-19] Making module private by default
!!               [2019-02-19] Turning file into a module
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module dbg_schemes

  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  integer :: nprobes

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_dbg, boundary_conditions_dbg, postprocess_dbg

contains
  !********************************************************************
  subroutine init_dbg (ux1,uy1,uz1,ep1,phi1)

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    integer :: code, ierror

    call debug_schemes()
    call MPI_ABORT(MPI_COMM_WORLD,code,ierror)

    return
  end subroutine init_dbg
  !********************************************************************
  subroutine boundary_conditions_dbg (ux,uy,uz,phi)

    USE param
    USE variables
    USE decomp_2d

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    return
  end subroutine boundary_conditions_dbg

  function sin_prec(x) result(y)
    USE decomp_2d, only : mytype
    real(mytype), intent(in) :: x
    real(mytype) :: y
#ifdef DOUBLE_PREC
    y = dsin(x)
#else
    y = sin(x)
#endif
  end function sin_prec

  function cos_prec(x) result(y)
    USE decomp_2d, only : mytype
    real(mytype), intent(in) :: x
    real(mytype) :: y
#ifdef DOUBLE_PREC
    y = dcos(x)
#else
    y = cos(x)
#endif
  end function cos_prec

  function abs_prec(x) result(y)
    USE decomp_2d, only : mytype
    real(mytype), intent(in) :: x
    real(mytype) :: y
#ifdef DOUBLE_PREC
    y = dabs(x)
#else
    y = abs(x)
#endif
  end function abs_prec

  subroutine debug_schemes()

    USE param
    USE variables
    USE decomp_2d
    USE var, only : nxmsize, nymsize, nzmsize

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: fx1, ffx1, ffxp1, fiffx1, fiffxp1, fxp1, dfdx1, dfdxp1, dfdxx1, dfdxxp1, di1, rand1
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: fy2, ffy2, ffyp2, fiffy2, fiffyp2, fyp2, dfdy2, dfdyp2, dfdyy2, dfdyyp2, di2, rand2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: fz3, ffz3, ffzp3, fiffz3, fiffzp3, fzp3, dfdz3, dfdzp3, dfdzz3, dfdzzp3, di3, rand3
    real(mytype), save, allocatable, dimension(:,:,:) :: test1,test11,test2,test22,test3,test33
    real(mytype) :: x,x1,y,y1,z,z1
    integer :: i,j,k
    character(len=30) :: filename

    allocate(test1(nxmsize,xsize(2),xsize(3)))
    allocate(test11(nxmsize,xsize(2),xsize(3)))
    allocate(test2(ysize(1),nymsize,ysize(3)))
    allocate(test22(ysize(1),nymsize,ysize(3)))
    allocate(test3(zsize(1),zsize(2),nzmsize))
    allocate(test33(zsize(1),zsize(2),nzmsize))

    nclx1=1
    nclxn=1
    ncly1=1
    nclyn=1
    nclz1=1
    nclzn=1
    call schemes()

    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             x = real(i-1,mytype)*dx*four*pi
             fx1(i,j,k)   = sin_prec(x)  !odd
             fxp1(i,j,k)  = cos_prec(x) !even
          enddo
       enddo
    enddo
    call derx (dfdx1 ,fx1 ,di1,sx,ffx ,fsx ,fwx ,xsize(1),xsize(2),xsize(3),0)
    call derx (dfdxp1,fxp1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
    call derxx (dfdxx1 ,fx1 ,di1,sx,sfx ,ssx ,swx ,xsize(1),xsize(2),xsize(3),0)
    call derxx (dfdxxp1,fxp1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
    if (nrank.eq.0) then
       write(filename,"('schemes_x',I1.1,I1.1,I1.1,I4.4)") jLES,nclx1,nclxn,nx
       open(67,file=trim(filename),status='unknown',form='formatted')
       do i=1,xsize(1)
          x = real(i-1,mytype)*dx
          write(67,'(9E14.6)') x,&
               four*pi*cos_prec(four*pi*x),dfdx1(i,1,1),&
               -four*pi*sin_prec(four*pi*x),dfdxp1(i,1,1),&
               -sixteen*pi*pi*sin_prec(four*pi*x),dfdxx1(i,1,1),&
               -sixteen*pi*pi*cos_prec(four*pi*x),dfdxxp1(i,1,1)
       enddo
       close(67)
    endif
    call derxvp(test1,fx1,di1,sx,cfx6,csx6,cwx6,xsize(1),nxmsize,xsize(2),xsize(3),0)
    call interxvp(test11,fxp1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)
    if (nrank.eq.0) then
       write(filename,"('schemes_xVP',I1.1,I1.1,I1.1,I4.4)") jLES,nclx1,nclxn,nx
       open(68,file=trim(filename),status='unknown',form='formatted')
       do i=1,nxmsize
          x1 = real(i-half,mytype)*dx
          write(68,'(5E14.6)') x1,&
               four*pi*cos_prec(four*pi*x1),test1(i,1,1),&
               cos_prec(four*pi*x1),test11(i,1,1)
       enddo
       close(68)
    endif
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,nxmsize
             x1 = real(i-0.5,mytype)*dx*four*pi
             test1(i,j,k) = cos_prec(x1)
          enddo
       enddo
    enddo
    call derxpv(fx1,test1,di1,sx,cfip6,csip6,cwip6,cfx6,csx6,cwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)
    call interxpv(fxp1,test1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)
    if (nrank.eq.0) then
       write(filename,"('schemes_xPV',I1.1,I1.1,I1.1,I4.4)") jLES,nclx1,nclxn,nx
       open(69,file=trim(filename),status='unknown',form='formatted')
       do i=1,nxmsize
          x = real(i-1.,mytype)*dx
          write(69,'(5E14.6)') x,&
               -four*pi*sin_prec(four*pi*x),fx1(i,1,1),&
               cos_prec(four*pi*x),fxp1(i,1,1)
       enddo
       close(69)
    endif

    ! FILTER
    call random_number(rand1)
    call filter(0.45_mytype)
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             x = real(i-1,mytype)*dx*four*pi
             ffx1(i,j,k)   = sin_prec(x)+sin_prec(ten*x)+0.25*rand1(i,j,k) !odd
             ffxp1(i,j,k)  = cos_prec(x)+cos_prec(ten*x)+0.25*rand1(i,j,k) !even
          enddo
       enddo
    enddo

    call filx (fiffx1  ,ffx1  ,di1,fisx,fiffx ,fifsx ,fifwx ,xsize(1),xsize(2),xsize(3),0)
    call filx (fiffxp1 ,ffxp1 ,di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1)
    if (nrank.eq.0) then
       write(filename,"('filter_x',I1.1,I1.1,I1.1,I4.4)") jLES,nclx1,nclxn,nx
       open(70,file=trim(filename),status='unknown',form='formatted')
       do i=1,xsize(2)
          x = real(i-1,mytype)*dx
          write(70,'(9E14.6)') x,&
               ffx1(i,1,1),fiffx1(i,1,1),&
               ffxp1(i,1,1),fiffxp1(i,1,1)
       enddo
       close(70)
    endif

    do k=1,ysize(3)
       do j=1,ysize(2)
          y = real(j-1,mytype)*dy*4*pi
          do i=1,ysize(1)
             fy2(i,j,k) = sin_prec(y)
             fyp2(i,j,k) = cos_prec(y)
          enddo
       enddo
    enddo
    call dery (dfdy2 ,fy2 ,di2,sy,ffy ,fsy ,fwy ,ppy,ysize(1),ysize(2),ysize(3),0)
    call dery (dfdyp2,fyp2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    call deryy (dfdyy2 ,fy2 ,di2,sy,sfy ,ssy ,swy ,ysize(1),ysize(2),ysize(3),0)
    call deryy (dfdyyp2,fyp2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
    if (nrank.eq.0) then
       write(filename,"('schemes_y',I1.1,I1.1,I1.1,I4.4)") jLES,ncly1,nclyn,ny
       open(67,file=trim(filename),status='unknown',form='formatted')
       do j=1,ysize(2)
          y = real(j-1,mytype)*dy
          write(67,'(9E14.6)') y,&
               four*pi*cos_prec(four*pi*y),dfdy2(1,j,1),&
               -four*pi*sin_prec(four*pi*y),dfdyp2(1,j,1),&
               -sixteen*pi*pi*sin_prec(four*pi*y),dfdyy2(1,j,1),&
               -sixteen*pi*pi*cos_prec(four*pi*y),dfdyyp2(1,j,1)
       enddo
       close(67)
    endif
    call deryvp(test2,fy2,di2,sy,cfy6,csy6,cwy6,ppyi,ysize(1),ysize(2),nymsize,ysize(3),0)
    call interyvp(test22,fyp2,di2,sy,cifyp6,cisyp6,ciwyp6,ysize(1),ysize(2),nymsize,ysize(3),1)
    if (nrank.eq.0) then
       write(filename,"('schemes_yVP',I1.1,I1.1,I1.1,I4.4)") jLES,ncly1,nclyn,ny
       open(68,file=trim(filename),status='unknown',form='formatted')
       do j=1,nymsize
          y1 = real(j-half,mytype)*dy
          write(68,'(5E14.6)') y1,&
               four*pi*cos_prec(four*pi*y1),test2(1,j,1),&
               cos_prec(four*pi*y1),test22(1,j,1)
       enddo
       close(68)
    endif
    do k=1,ysize(3)
       do j=1,nymsize
          y1 = real(j-half,mytype)*dy*4*pi
          do i=1,ysize(1)
             test2(i,j,k) = cos_prec(y1)
          enddo
       enddo
    enddo
    call derypv(fy2,test2,di2,sy,cfip6y,csip6y,cwip6y,cfy6,csy6,cwy6,ppy,&
         ysize(1),nymsize,ysize(2),ysize(3),1)
    call interypv(fyp2,test2,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
         ysize(1),nymsize,ysize(2),ysize(3),1)
    if (nrank.eq.0) then
       write(filename,"('schemes_yPV',I1.1,I1.1,I1.1,I4.4)") jLES,ncly1,nclyn,ny
       open(69,file=trim(filename),status='unknown',form='formatted')
       do j=1,ysize(2)
          y = real(j-1,mytype)*dy
          write(69,'(5E14.6)') y,&
               -four*pi*sin_prec(four*pi*y),fy2(1,j,1),&
               cos_prec(four*pi*y),fyp2(1,j,1)
       enddo
       close(69)
    endif

    ! FILTER
    call filter(0.45_mytype)
    do k=1,ysize(3)
       do j=1,ysize(2)
          y = real(j-1,mytype)*dy*four*pi
          do i=1,ysize(1)
             ffy2(i,j,k)  = sin_prec(y)
             ffyp2(i,j,k) = cos_prec(y)
          enddo
       enddo
    enddo

    call fily (fiffy2  ,ffy2  ,di2,fisy,fiffy ,fifsy ,fifwy ,ppy,ysize(1),ysize(2),ysize(3),0)
    call fily (fiffyp2 ,ffyp2 ,di2,fisy,fiffyp,fifsyp,fifwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    if (nrank.eq.0) then
       write(filename,"('filter_y',I1.1,I1.1,I1.1,I4.4)") jLES,ncly1,nclyn,ny
       open(70,file=trim(filename),status='unknown',form='formatted')
       do j=1,ysize(2)
          y = real(j-1,mytype)*dy
          write(70,'(9E14.6)') y,&
               ffy2(1,j,1),fiffy2(1,j,1),&
               ffyp2(1,j,1),fiffyp2(1,j,1)
       enddo
       close(70)
    endif

    do k=1,zsize(3)
       z = real(k-1,mytype)*dz*4*pi
       do j=1,zsize(2)
          do i=1,zsize(1)
             fz3(i,j,k) = sin_prec(z)
             fzp3(i,j,k) = cos_prec(z)
          enddo
       enddo
    enddo
    call derz (dfdz3 ,fz3 ,di3,sz,ffz ,fsz ,fwz ,zsize(1),zsize(2),zsize(3),0)
    call derz (dfdzp3,fzp3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
    call derzz (dfdzz3 ,fz3 ,di3,sz,sfz ,ssz ,swz ,zsize(1),zsize(2),zsize(3),0)
    call derzz (dfdzzp3,fzp3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
    if (nrank.eq.0) then
       write(filename,"('schemes_z',I1.1,I1.1,I1.1,I4.4)") jLES,nclz1,nclzn,nz
       open(67,file=trim(filename),status='unknown',form='formatted')
       do k=1,zsize(3)
          z = real(k-1,mytype)*dz
          write(67,'(9E14.6)') z,&
               four*pi*cos_prec(four*pi*z),dfdz3(1,1,k),&
               -four*pi*sin_prec(four*pi*z),dfdzp3(1,1,k),&
               -sixteen*pi*pi*sin_prec(four*pi*z),dfdzz3(1,1,k),&
               -sixteen*pi*pi*cos_prec(four*pi*z),dfdzzp3(1,1,k)
       enddo
       close(67)
    endif
    call derzvp(test3,fz3,di3,sz,cfz6,csz6,cwz6,zsize(1),zsize(2),zsize(3),nzmsize,0)
    call interzvp(test33,fzp3,di3,sz,cifzp6,ciszp6,ciwzp6,zsize(1),zsize(2),zsize(3),nzmsize,1)
    if (nrank.eq.0) then
       write(filename,"('schemes_zVP',I1.1,I1.1,I1.1,I4.4)") jLES,nclz1,nclzn,nz
       open(68,file=trim(filename),status='unknown',form='formatted')
       do k=1,nzmsize
          z1 = real(k-half,mytype)*dz
          write(68,'(5E14.6)') z1,&
               four*pi*cos_prec(four*pi*z1),test3(1,1,k),&
               cos_prec(four*pi*z1),test33(1,1,k)
       enddo
       close(68)
    endif

    do k=1,nzmsize
       z1 = real(k-half,mytype)*dz*4*pi
       do j=1,zsize(2)
          do i=1,zsize(1)
             test3(i,j,k) = cos_prec(z1)
          enddo
       enddo
    enddo
    call derzpv(fz3,test3,di3,sz,cfip6z,csip6z,cwip6z,cfz6,csz6,cwz6,zsize(1),zsize(2),nzmsize,zsize(3),1)
    call interzpv(fzp3,test3,di3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,zsize(1),zsize(2),nzmsize,zsize(3),1)
    if (nrank.eq.0) then
       write(filename,"('schemes_zPV',I1.1,I1.1,I1.1,I4.4)") jLES,nclz1,nclzn,nz
       open(69,file=trim(filename),status='unknown',form='formatted')
       do k=1,zsize(3)
          z = real(k-1,mytype)*dz
          write(69,'(5E14.6)') z,&
               -four*pi*sin_prec(four*pi*z),fz3(1,1,k),&
               cos_prec(four*pi*z),fzp3(1,1,k)
       enddo
       close(69)
    endif

    ! FILTER
    call filter(0.45_mytype)
    do k=1,zsize(3)
       z = real(k-1,mytype)*dz*four*pi
       do j=1,zsize(2)
          do i=1,zsize(1)
             ffz3(i,j,k)  = sin_prec(z)
             ffzp3(i,j,k) = cos_prec(z)
          enddo
       enddo
    enddo
    call filz (fiffz3  ,ffz3  ,di3,fisz,fiffz ,fifsz ,fifwz ,zsize(1),zsize(2),zsize(3),0)
    call filz (fiffzp3 ,ffzp3 ,di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1)
    if (nrank.eq.0) then
       write(filename,"('filter_z',I1.1,I1.1,I1.1,I4.4)") jLES,nclz1,nclzn,nz
       open(70,file=trim(filename),status='unknown',form='formatted')
       do k=1,zsize(3)
          z = real(k-1,mytype)*dz
          write(70,'(9E14.6)') z,&
               ffz3 (1,1,k),fiffz3 (1,1,k),&
               ffzp3(1,1,k),fiffzp3(1,1,k)
       enddo
       close(70)
    endif

    !###############################################################
    !###############################################################

    nclx1=2
    nclxn=2
    ncly1=2
    nclyn=2
    nclz1=2
    nclzn=2
    call schemes()

    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             x = real(i-1,mytype)*dx*four*pi
             fx1(i,j,k) = sin_prec(x)  !odd
          enddo
       enddo
    enddo
    call derx (dfdx1 ,fx1 ,di1,sx,ffx ,fsx ,fwx ,xsize(1),xsize(2),xsize(3),0)
    call derxx (dfdxx1 ,fx1 ,di1,sx,sfx ,ssx ,swx ,xsize(1),xsize(2),xsize(3),0)
    if (nrank.eq.0) then
       write(filename,"('schemes_x',I1.1,I1.1,I1.1,I4.4)") jLES,nclx1,nclxn,nx
       open(67,file=trim(filename),status='unknown',form='formatted')
       do i=1,xsize(1)
          x = real(i-1,mytype)*dx
          write(67,'(5E14.6)') x,&
               four*pi*cos_prec(four*pi*x),dfdx1(i,1,1),&
               -sixteen*pi*pi*sin_prec(four*pi*x),dfdxx1(i,1,1)
       enddo
       close(67)
    endif

    ! FILTER
    call random_number(rand1)
    call filter(-0.45_mytype)
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             x = real(i-1,mytype)*dx*four*pi
             ffx1(i,j,k)   = sin_prec(x)+sin_prec(two*x)+rand1(i,j,k) !odd
             ffxp1(i,j,k)  = cos_prec(x)+cos_prec(two*x)+rand1(i,j,k) !even
          enddo
       enddo
    enddo

    call filx (fiffx1  ,ffx1  ,di1,fisx,fiffx ,fifsx ,fifwx ,xsize(1),xsize(2),xsize(3),0)
    call filx (fiffxp1 ,ffxp1 ,di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1)
    if (nrank.eq.0) then
       write(filename,"('filter_x',I1.1,I1.1,I1.1,I4.4)") jLES,nclx1,nclxn,nx
       open(70,file=trim(filename),status='unknown',form='formatted')
       do i=1,xsize(2)
          x = real(i-1,mytype)*dx
          write(70,'(9E14.6)') x,&
               ffx1(i,1,1),fiffx1(i,1,1),&
               ffxp1(i,1,1),fiffxp1(i,1,1)
       enddo
       close(70)
    endif

    do k=1,ysize(3)
       do j=1,ysize(2)
          y = real(j-1,mytype)*dy*4*pi
          do i=1,ysize(1)
             fy2(i,j,k) = sin_prec(y)
          enddo
       enddo
    enddo
    call dery (dfdy2 ,fy2 ,di2,sy,ffy ,fsy ,fwy ,ppy,ysize(1),ysize(2),ysize(3),0)
    call deryy (dfdyy2 ,fy2 ,di2,sy,sfy ,ssy ,swy ,ysize(1),ysize(2),ysize(3),0)
    if (nrank.eq.0) then
       write(filename,"('schemes_y',I1.1,I1.1,I1.1,I4.4)") jLES,ncly1,nclyn,ny
       open(67,file=trim(filename),status='unknown',form='formatted')
       do j=1,ysize(2)
          y = real(j-1,mytype)*dy
          write(67,'(5E14.6)') y,&
               four*pi*cos_prec(four*pi*y),dfdy2(1,j,1),&
               -sixteen*pi*pi*sin_prec(four*pi*y),dfdyy2(1,j,1)
       enddo
       close(67)
    endif
    ! FILTER
    call filter(0.45_mytype)
    do k=1,ysize(3)
       do j=1,ysize(2)
          y = real(j-1,mytype)*dy*four*pi
          do i=1,ysize(1)
             ffy2(i,j,k)  = sin_prec(y)
             ffyp2(i,j,k) = cos_prec(y)
          enddo
       enddo
    enddo

    call fily (fiffy2  ,ffy2  ,di2,fisy,fiffy ,fifsy ,fifwy ,ppy,ysize(1),ysize(2),ysize(3),0)
    call fily (fiffyp2 ,ffyp2 ,di2,fisy,fiffyp,fifsyp,fifwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    if (nrank.eq.0) then
       write(filename,"('filter_y',I1.1,I1.1,I1.1,I4.4)") jLES,ncly1,nclyn,ny
       open(70,file=trim(filename),status='unknown',form='formatted')
       do j=1,ysize(2)
          y = real(j-1,mytype)*dy
          write(70,'(9E14.6)') y,&
               ffy2(1,j,1),fiffy2(1,j,1),&
               ffyp2(1,j,1),fiffyp2(1,j,1)
       enddo
       close(70)
    endif
    do k=1,zsize(3)
       z = real(k-1,mytype)*dz*4*pi
       do j=1,zsize(2)
          do i=1,zsize(1)
             fz3(i,j,k) = sin_prec(z)
          enddo
       enddo
    enddo
    call derz (dfdz3 ,fz3 ,di3,sz,ffz ,fsz ,fwz ,zsize(1),zsize(2),zsize(3),0)
    call derzz (dfdzz3 ,fz3 ,di3,sz,sfz ,ssz ,swz ,zsize(1),zsize(2),zsize(3),0)
    if (nrank.eq.0) then
       write(filename,"('schemes_z',I1.1,I1.1,I1.1,I4.4)") jLES,nclz1,nclzn,nz
       open(67,file=trim(filename),status='unknown',form='formatted')
       do k=1,zsize(3)
          z = real(k-1,mytype)*dz
          write(67,'(5E14.6)') z,&
               four*pi*cos_prec(four*pi*z),dfdz3(1,1,k),&
               -sixteen*pi*pi*sin_prec(four*pi*z),dfdzz3(1,1,k)
       enddo
       close(67)
    endif
    ! FILTER
    call filter(0.45_mytype)
    do k=1,zsize(3)
       z = real(k-1,mytype)*dz*four*pi
       do j=1,zsize(2)
          do i=1,zsize(1)
             ffz3(i,j,k)  = sin_prec(z)
             ffzp3(i,j,k) = cos_prec(z)
          enddo
       enddo
    enddo
    call filz (fiffz3  ,ffz3  ,di3,fisz,fiffz ,fifsz ,fifwz ,zsize(1),zsize(2),zsize(3),0)
    call filz (fiffzp3 ,ffzp3 ,di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1)
    if (nrank.eq.0) then
       write(filename,"('filter_z',I1.1,I1.1,I1.1,I4.4)") jLES,nclz1,nclzn,nz
       open(70,file=trim(filename),status='unknown',form='formatted')
       do k=1,zsize(3)
          z = real(k-1,mytype)*dz
          write(70,'(9E14.6)') z,&
               ffz3 (1,1,k),fiffz3 (1,1,k),&
               ffzp3(1,1,k),fiffzp3(1,1,k)
       enddo
       close(70)
    endif

    !###############################################################
    !###############################################################

    nclx1=2
    nclxn=1
    ncly1=2
    nclyn=1
    nclz1=2
    nclzn=1
    call schemes()

    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             x = real(i-1,mytype)*dx*four*pi
             fx1(i,j,k) = sin_prec(x)  !odd
             fxp1(i,j,k) = cos_prec(x) !even
          enddo
       enddo
    enddo
    call derx (dfdx1 ,fx1 ,di1,sx,ffx ,fsx ,fwx ,xsize(1),xsize(2),xsize(3),0)
    call derx (dfdxp1,fxp1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
    call derxx (dfdxx1 ,fx1 ,di1,sx,sfx ,ssx ,swx ,xsize(1),xsize(2),xsize(3),0)
    call derxx (dfdxxp1,fxp1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
    if (nrank.eq.0) then
       write(filename,"('schemes_x',I1.1,I1.1,I1.1,I4.4)") jLES,nclx1,nclxn,nx
       open(67,file=trim(filename),status='unknown',form='formatted')
       do i=1,xsize(1)
          x = real(i-1,mytype)*dx
          write(67,'(9E14.6)') x,&
               four*pi*cos_prec(four*pi*x),dfdx1(i,1,1),&
               -four*pi*sin_prec(four*pi*x),dfdxp1(i,1,1),&
               -sixteen*pi*pi*sin_prec(four*pi*x),dfdxx1(i,1,1),&
               -sixteen*pi*pi*cos_prec(four*pi*x),dfdxxp1(i,1,1)
       enddo
       close(67)
    endif
    ! FILTER
    call random_number(rand1)
    call filter(0.45_mytype)
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             x = real(i-1,mytype)*dx*four*pi
             ffx1(i,j,k)   = sin_prec(x)+sin_prec(ten*x)+0.25*rand1(i,j,k) !odd
             ffxp1(i,j,k)  = cos_prec(x)+cos_prec(ten*x)+0.25*rand1(i,j,k) !even
          enddo
       enddo
    enddo

    call filx (fiffx1  ,ffx1  ,di1,fisx,fiffx ,fifsx ,fifwx ,xsize(1),xsize(2),xsize(3),0)
    call filx (fiffxp1 ,ffxp1 ,di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1)
    if (nrank.eq.0) then
       write(filename,"('filter_x',I1.1,I1.1,I1.1,I4.4)") jLES,nclx1,nclxn,nx
       open(70,file=trim(filename),status='unknown',form='formatted')
       do i=1,xsize(2)
          x = real(i-1,mytype)*dx
          write(70,'(9E14.6)') x,&
               ffx1(i,1,1),fiffx1(i,1,1),&
               ffxp1(i,1,1),fiffxp1(i,1,1)
       enddo
       close(70)
    endif

    do k=1,ysize(3)
       do j=1,ysize(2)
          y = real(j-1,mytype)*dy*4*pi
          do i=1,ysize(1)
             fy2(i,j,k) = sin_prec(y)
             fyp2(i,j,k) = cos_prec(y)
          enddo
       enddo
    enddo
    call dery (dfdy2 ,fy2 ,di2,sy,ffy ,fsy ,fwy ,ppy,ysize(1),ysize(2),ysize(3),0)
    call dery (dfdyp2,fyp2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    call deryy (dfdyy2 ,fy2 ,di2,sy,sfy ,ssy ,swy ,ysize(1),ysize(2),ysize(3),0)
    call deryy (dfdyyp2,fyp2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
    if (nrank.eq.0) then
       write(filename,"('schemes_y',I1.1,I1.1,I1.1,I4.4)") jLES,ncly1,nclyn,ny
       open(67,file=trim(filename),status='unknown',form='formatted')
       do j=1,ysize(2)
          y = real(j-1,mytype)*dy
          write(67,'(9E14.6)') y,&
               four*pi*cos_prec(four*pi*y),dfdy2(1,j,1),&
               -four*pi*sin_prec(four*pi*y),dfdyp2(1,j,1),&
               -sixteen*pi*pi*sin_prec(four*pi*y),dfdyy2(1,j,1),&
               -sixteen*pi*pi*cos_prec(four*pi*y),dfdyyp2(1,j,1)
       enddo
       close(67)
    endif
    ! FILTER
    call filter(0.45_mytype)
    do k=1,ysize(3)
       do j=1,ysize(2)
          y = real(j-1,mytype)*dy*four*pi
          do i=1,ysize(1)
             ffy2(i,j,k)  = sin_prec(y)
             ffyp2(i,j,k) = cos_prec(y)
          enddo
       enddo
    enddo

    call fily (fiffy2  ,ffy2  ,di2,fisy,fiffy ,fifsy ,fifwy ,ppy,ysize(1),ysize(2),ysize(3),0)
    call fily (fiffyp2 ,ffyp2 ,di2,fisy,fiffyp,fifsyp,fifwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    if (nrank.eq.0) then
       write(filename,"('filter_y',I1.1,I1.1,I1.1,I4.4)") jLES,ncly1,nclyn,ny
       open(70,file=trim(filename),status='unknown',form='formatted')
       do j=1,ysize(2)
          y = real(j-1,mytype)*dy
          write(70,'(9E14.6)') y,&
               ffy2(1,j,1),fiffy2(1,j,1),&
               ffyp2(1,j,1),fiffyp2(1,j,1)
       enddo
       close(70)
    endif

    do k=1,zsize(3)
       z = real(k-1,mytype)*dz*4*pi
       do j=1,zsize(2)
          do i=1,zsize(1)
             fz3(i,j,k) = sin_prec(z)
             fzp3(i,j,k) = cos_prec(z)
          enddo
       enddo
    enddo
    call derz (dfdz3 ,fz3 ,di3,sz,ffz ,fsz ,fwz ,zsize(1),zsize(2),zsize(3),0)
    call derz (dfdzp3,fzp3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
    call derzz (dfdzz3 ,fz3 ,di3,sz,sfz ,ssz ,swz ,zsize(1),zsize(2),zsize(3),0)
    call derzz (dfdzzp3,fzp3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
    if (nrank.eq.0) then
       write(filename,"('schemes_z',I1.1,I1.1,I1.1,I4.4)") jLES,nclz1,nclzn,nz
       open(67,file=trim(filename),status='unknown',form='formatted')
       do k=1,zsize(3)
          z = real(k-1,mytype)*dz
          write(67,'(9E14.6)') z,&
               four*pi*cos_prec(four*pi*z),dfdz3(1,1,k),&
               -four*pi*sin_prec(four*pi*z),dfdzp3(1,1,k),&
               -sixteen*pi*pi*sin_prec(four*pi*z),dfdzz3(1,1,k),&
               -sixteen*pi*pi*cos_prec(four*pi*z),dfdzzp3(1,1,k)
       enddo
       close(67)
    endif
    ! FILTER
    call filter(0.45_mytype)
    do k=1,zsize(3)
       z = real(k-1,mytype)*dz*four*pi
       do j=1,zsize(2)
          do i=1,zsize(1)
             ffz3(i,j,k)  = sin_prec(z)
             ffzp3(i,j,k) = cos_prec(z)
          enddo
       enddo
    enddo
    call filz (fiffz3  ,ffz3  ,di3,fisz,fiffz ,fifsz ,fifwz ,zsize(1),zsize(2),zsize(3),0)
    call filz (fiffzp3 ,ffzp3 ,di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1)
    if (nrank.eq.0) then
       write(filename,"('filter_z',I1.1,I1.1,I1.1,I4.4)") jLES,nclz1,nclzn,nz
       open(70,file=trim(filename),status='unknown',form='formatted')
       do k=1,zsize(3)
          z = real(k-1,mytype)*dz
          write(70,'(9E14.6)') z,&
               ffz3 (1,1,k),fiffz3 (1,1,k),&
               ffzp3(1,1,k),fiffzp3(1,1,k)
       enddo
       close(70)
    endif

    !###############################################################
    !###############################################################

    nclx1=1
    nclxn=2
    ncly1=1
    nclyn=2
    nclz1=1
    nclzn=2
    call schemes()

    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             x = real(i-1,mytype)*dx*four*pi
             fx1(i,j,k) = sin_prec(x)  !odd
             fxp1(i,j,k) = cos_prec(x) !even
          enddo
       enddo
    enddo
    call derx (dfdx1 ,fx1 ,di1,sx,ffx ,fsx ,fwx ,xsize(1),xsize(2),xsize(3),0)
    call derx (dfdxp1,fxp1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
    call derxx (dfdxx1 ,fx1 ,di1,sx,sfx ,ssx ,swx ,xsize(1),xsize(2),xsize(3),0)
    call derxx (dfdxxp1,fxp1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
    if (nrank.eq.0) then
       write(filename,"('schemes_x',I1.1,I1.1,I1.1,I4.4)") jLES,nclx1,nclxn,nx
       open(67,file=trim(filename),status='unknown',form='formatted')
       do i=1,xsize(1)
          x = real(i-1,mytype)*dx
          write(67,'(9E14.6)') x,&
               four*pi*cos_prec(four*pi*x),dfdx1(i,1,1),&
               -four*pi*sin_prec(four*pi*x),dfdxp1(i,1,1),&
               -sixteen*pi*pi*sin_prec(four*pi*x),dfdxx1(i,1,1),&
               -sixteen*pi*pi*cos_prec(four*pi*x),dfdxxp1(i,1,1)
       enddo
       close(67)
    endif
    ! FILTER
    call random_number(rand1)
    call filter(0.45_mytype)
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             x = real(i-1,mytype)*dx*four*pi
             ffx1(i,j,k)   = sin_prec(x)!+sin_prec(ten*x)+0.25*rand1(i,j,k) !odd
             ffxp1(i,j,k)  = cos_prec(x)!+cos_prec(ten*x)+0.25*rand1(i,j,k) !even
          enddo
       enddo
    enddo

    call filx (fiffx1  ,ffx1  ,di1,fisx,fiffx ,fifsx ,fifwx ,xsize(1),xsize(2),xsize(3),0)
    call filx (fiffxp1 ,ffxp1 ,di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1)
    if (nrank.eq.0) then
       write(filename,"('filter_x',I1.1,I1.1,I1.1,I4.4)") jLES,nclx1,nclxn,nx
       open(70,file=trim(filename),status='unknown',form='formatted')
       do i=1,xsize(2)
          x = real(i-1,mytype)*dx
          write(70,'(9E14.6)') x,&
               ffx1(i,1,1),fiffx1(i,1,1),&
               ffxp1(i,1,1),fiffxp1(i,1,1)
       enddo
       close(70)
    endif

    do k=1,ysize(3)
       do j=1,ysize(2)
          y = real(j-1,mytype)*dy*4*pi
          do i=1,ysize(1)
             fy2(i,j,k) = sin_prec(y)
             fyp2(i,j,k) = cos_prec(y)
          enddo
       enddo
    enddo
    call dery (dfdy2 ,fy2 ,di2,sy,ffy ,fsy ,fwy ,ppy,ysize(1),ysize(2),ysize(3),0)
    call dery (dfdyp2,fyp2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    call deryy (dfdyy2 ,fy2 ,di2,sy,sfy ,ssy ,swy ,ysize(1),ysize(2),ysize(3),0)
    call deryy (dfdyyp2,fyp2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
    if (nrank.eq.0) then
       write(filename,"('schemes_y',I1.1,I1.1,I1.1,I4.4)") jLES,ncly1,nclyn,ny
       open(67,file=trim(filename),status='unknown',form='formatted')
       do j=1,ysize(2)
          y = real(j-1,mytype)*dy
          write(67,'(9E14.6)') y,&
               four*pi*cos_prec(four*pi*y),dfdy2(1,j,1),&
               -four*pi*sin_prec(four*pi*y),dfdyp2(1,j,1),&
               -sixteen*pi*pi*sin_prec(four*pi*y),dfdyy2(1,j,1),&
               -sixteen*pi*pi*cos_prec(four*pi*y),dfdyyp2(1,j,1)
       enddo
       close(67)
    endif
    ! FILTER
    call filter(0.45_mytype)
    do k=1,ysize(3)
       do j=1,ysize(2)
          y = real(j-1,mytype)*dy*four*pi
          do i=1,ysize(1)
             ffy2(i,j,k)  = sin_prec(y)
             ffyp2(i,j,k) = cos_prec(y)
          enddo
       enddo
    enddo

    call fily (fiffy2  ,ffy2  ,di2,fisy,fiffy ,fifsy ,fifwy ,ppy,ysize(1),ysize(2),ysize(3),0)
    call fily (fiffyp2 ,ffyp2 ,di2,fisy,fiffyp,fifsyp,fifwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    if (nrank.eq.0) then
       write(filename,"('filter_y',I1.1,I1.1,I1.1,I4.4)") jLES,ncly1,nclyn,ny
       open(70,file=trim(filename),status='unknown',form='formatted')
       do j=1,ysize(2)
          y = real(j-1,mytype)*dy
          write(70,'(9E14.6)') y,&
               ffy2(1,j,1),fiffy2(1,j,1),&
               ffyp2(1,j,1),fiffyp2(1,j,1)
       enddo
       close(70)
    endif

    do k=1,zsize(3)
       z = real(k-1,mytype)*dz*4*pi
       do j=1,zsize(2)
          do i=1,zsize(1)
             fz3(i,j,k) = sin_prec(z)
             fzp3(i,j,k) = cos_prec(z)
          enddo
       enddo
    enddo
    call derz (dfdz3 ,fz3 ,di3,sz,ffz ,fsz ,fwz ,zsize(1),zsize(2),zsize(3),0)
    call derz (dfdzp3,fzp3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
    call derzz (dfdzz3 ,fz3 ,di3,sz,sfz ,ssz ,swz ,zsize(1),zsize(2),zsize(3),0)
    call derzz (dfdzzp3,fzp3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
    if (nrank.eq.0) then
       write(filename,"('schemes_z',I1.1,I1.1,I1.1,I4.4)") jLES,nclz1,nclzn,nz
       open(67,file=trim(filename),status='unknown',form='formatted')
       do k=1,zsize(3)
          z = real(k-1,mytype)*dz
          write(67,'(9E14.6)') z,&
               four*pi*cos_prec(four*pi*z),dfdz3(1,1,k),&
               -four*pi*sin_prec(four*pi*z),dfdzp3(1,1,k),&
               -sixteen*pi*pi*sin_prec(four*pi*z),dfdzz3(1,1,k),&
               -sixteen*pi*pi*cos_prec(four*pi*z),dfdzzp3(1,1,k)
       enddo
       close(67)
    endif
    ! FILTER
    call filter(0.45_mytype)
    do k=1,zsize(3)
       z = real(k-1,mytype)*dz*four*pi
       do j=1,zsize(2)
          do i=1,zsize(1)
             ffz3(i,j,k)  = sin_prec(z)
             ffzp3(i,j,k) = cos_prec(z)
          enddo
       enddo
    enddo
    call filz (fiffz3  ,ffz3  ,di3,fisz,fiffz ,fifsz ,fifwz ,zsize(1),zsize(2),zsize(3),0)
    call filz (fiffzp3 ,ffzp3 ,di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1)
    if (nrank.eq.0) then
       write(filename,"('filter_z',I1.1,I1.1,I1.1,I4.4)") jLES,nclz1,nclzn,nz
       open(70,file=trim(filename),status='unknown',form='formatted')
       do k=1,zsize(3)
          z = real(k-1,mytype)*dz
          write(70,'(9E14.6)') z,&
               ffz3 (1,1,k),fiffz3 (1,1,k),&
               ffzp3(1,1,k),fiffzp3(1,1,k)
       enddo
       close(70)
    endif

    !###############################################################
    !###############################################################
    deallocate(test1)
    deallocate(test11)
    deallocate(test2)
    deallocate(test22)
    deallocate(test3)
    deallocate(test33)
    nclx1=0
    nclxn=0
    ncly1=0
    nclyn=0
    nclz1=0
    nclzn=0
    nclx=.true.
    ncly=.true.
    nclz=.true.
    nxm=nx
    nym=ny
    nzm=nz
    nxmsize=xsize(1)
    nymsize=ysize(2)
    nzmsize=zsize(3)
    dx=xlx/real(nxm,mytype)
    dy=yly/real(nym,mytype)
    dz=zlz/real(nzm,mytype)
    dx2=dx*dx
    dy2=dy*dy
    dz2=dz*dz
    call schemes()
    allocate(test1(nxmsize,xsize(2),xsize(3)))
    allocate(test11(nxmsize,xsize(2),xsize(3)))
    allocate(test2(ysize(1),nymsize,ysize(3)))
    allocate(test22(ysize(1),nymsize,ysize(3)))
    allocate(test3(zsize(1),zsize(2),nzmsize))
    allocate(test33(zsize(1),zsize(2),nzmsize))

    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             x = real(i-1,mytype)*dx*four*pi
             fx1(i,j,k) = sin_prec(x)  !odd
             fxp1(i,j,k) = cos_prec(x)
          enddo
       enddo
    enddo
    call derx (dfdx1 ,fx1 ,di1,sx,ffx ,fsx ,fwx ,xsize(1),xsize(2),xsize(3),0)
    call derxx (dfdxx1 ,fx1 ,di1,sx,sfx ,ssx ,swx ,xsize(1),xsize(2),xsize(3),0)
    if (nrank.eq.0) then
       write(filename,"('schemes_x',I1.1,I1.1,I1.1,I4.4)") jLES,nclx1,nclxn,nx
       open(67,file=trim(filename),status='unknown',form='formatted')
       do i=1,xsize(1)
          x = real(i-1,mytype)*dx
          write(67,'(5E14.6)') x,&
               four*pi*cos_prec(four*pi*x),dfdx1(i,1,1),&
               -sixteen*pi*pi*sin_prec(four*pi*x),dfdxx1(i,1,1)
       enddo
       close(67)
    endif
    call derxvp(test1,fx1,di1,sx,cfx6,csx6,cwx6,xsize(1),nxmsize,xsize(2),xsize(3),0)
    call interxvp(test11,fxp1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)
    if (nrank.eq.0) then
       write(filename,"('schemes_xVP',I1.1,I1.1,I1.1,I4.4)") jLES,nclx1,nclxn,nx
       open(68,file=trim(filename),status='unknown',form='formatted')
       do i=1,nxmsize
          x1 = real(i-half,mytype)*dx
          write(68,'(5E14.6)') x1,&
               four*pi*cos_prec(four*pi*x1),test1(i,1,1),&
               cos_prec(four*pi*x1),test11(i,1,1)
       enddo
       close(68)
    endif
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,nxmsize
             x1 = real(i-0.5,mytype)*dx*four*pi
             test1(i,j,k) = cos_prec(x1)
          enddo
       enddo
    enddo
    call derxpv(fx1,test1,di1,sx,cfip6,csip6,cwip6,cfx6,csx6,cwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)
    call interxpv(fxp1,test11,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)
    if (nrank.eq.0) then
       write(filename,"('schemes_xPV',I1.1,I1.1,I1.1,I4.4)") jLES,nclx1,nclxn,nx
       open(69,file=trim(filename),status='unknown',form='formatted')
       do i=1,nxmsize
          x = real(i-1.,mytype)*dx
          write(69,'(5E14.6)') x,&
               -four*pi*sin_prec(four*pi*x),fx1(i,1,1),&
               cos_prec(four*pi*x),fxp1(i,1,1)
       enddo
       close(69)
    endif

    ! FILTER
    call random_number(rand1)
    call filter(-0.45_mytype)
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             x = real(i-1,mytype)*dx*four*pi
             ffx1(i,j,k)   = sin_prec(x)+sin_prec(ten*x)+rand1(i,j,k) !odd
             ffxp1(i,j,k)  = cos_prec(x)+cos_prec(ten*x)+rand1(i,j,k) !even
          enddo
       enddo
    enddo

    call filx (fiffx1  ,ffx1  ,di1,fisx,fiffx ,fifsx ,fifwx ,xsize(1),xsize(2),xsize(3),0)
    call filx (fiffxp1 ,ffxp1 ,di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1)
    if (nrank.eq.0) then
       write(filename,"('filter_x',I1.1,I1.1,I1.1,I4.4)") jLES,nclx1,nclxn,nx
       open(70,file=trim(filename),status='unknown',form='formatted')
       do i=1,xsize(2)
          x = real(i-1,mytype)*dx
          write(70,'(9E14.6)') x,&
               ffx1(i,1,1),fiffx1(i,1,1),&
               ffxp1(i,1,1),fiffxp1(i,1,1)
       enddo
       close(70)
    endif

    do k=1,ysize(3)
       do j=1,ysize(2)
          y = real(j-1,mytype)*dy*4*pi
          do i=1,ysize(1)
             fy2(i,j,k) = sin_prec(y)
             fyp2(i,j,k) = cos_prec(y)
          enddo
       enddo
    enddo
    call dery (dfdy2 ,fy2 ,di2,sy,ffy ,fsy ,fwy ,ppy,ysize(1),ysize(2),ysize(3),0)
    call deryy (dfdyy2 ,fy2 ,di2,sy,sfy ,ssy ,swy ,ysize(1),ysize(2),ysize(3),0)
    if (nrank.eq.0) then
       write(filename,"('schemes_y',I1.1,I1.1,I1.1,I4.4)") jLES,ncly1,nclyn,ny
       open(67,file=trim(filename),status='unknown',form='formatted')
       do j=1,ysize(2)
          y = real(j-1,mytype)*dy
          write(67,'(5E14.6)') y,&
               four*pi*cos_prec(four*pi*y),dfdy2(1,j,1),&
               -sixteen*pi*pi*sin_prec(four*pi*y),dfdyy2(1,j,1)
       enddo
       close(67)
    endif

    call deryvp(test2,fy2,di2,sy,cfy6,csy6,cwy6,ppyi,ysize(1),ysize(2),nymsize,ysize(3),0)
    call interyvp(test22,fyp2,di2,sy,cifyp6,cisyp6,ciwyp6,ysize(1),ysize(2),nymsize,ysize(3),1)
    if (nrank.eq.0) then
       write(filename,"('schemes_yVP',I1.1,I1.1,I1.1,I4.4)") jLES,ncly1,nclyn,ny
       open(68,file=trim(filename),status='unknown',form='formatted')
       do j=1,nymsize
          y1 = real(j-half,mytype)*dy
          write(68,'(5E14.6)') y1,&
               four*pi*cos_prec(four*pi*y1),test2(1,j,1),&
               cos_prec(four*pi*y1),test22(1,j,1)
       enddo
       close(68)
    endif
    do k=1,ysize(3)
       do j=1,nymsize
          y1 = real(j-0.5,mytype)*dy*4*pi
          do i=1,ysize(1)
             test2(i,j,k) = cos_prec(y1)
          enddo
       enddo
    enddo
    call derypv(fy2,test2,di2,sy,cfip6y,csip6y,cwip6y,cfy6,csy6,cwy6,ppy,&
         ysize(1),nymsize,ysize(2),ysize(3),1)
    call interypv(fyp2,test2,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
         ysize(1),nymsize,ysize(2),ysize(3),1)
    if (nrank.eq.0) then
       write(filename,"('schemes_yPV',I1.1,I1.1,I1.1,I4.4)") jLES,ncly1,nclyn,ny
       open(69,file=trim(filename),status='unknown',form='formatted')
       do j=1,ysize(2)
          y = real(j-1,mytype)*dy
          write(69,'(5E14.6)') y,&
               -four*pi*sin_prec(four*pi*y),fy2(1,j,1),&
               cos_prec(four*pi*y),fyp2(1,j,1)
       enddo
       close(69)
    endif
    ! FILTER
    call filter(0.45_mytype)
    do k=1,ysize(3)
       do j=1,ysize(2)
          y = real(j-1,mytype)*dy*four*pi
          do i=1,ysize(1)
             ffy2(i,j,k)  = sin_prec(y)
             ffyp2(i,j,k) = cos_prec(y)
          enddo
       enddo
    enddo

    call fily (fiffy2  ,ffy2  ,di2,fisy,fiffy ,fifsy ,fifwy ,ppy,ysize(1),ysize(2),ysize(3),0)
    call fily (fiffyp2 ,ffyp2 ,di2,fisy,fiffyp,fifsyp,fifwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    if (nrank.eq.0) then
       write(filename,"('filter_y',I1.1,I1.1,I1.1,I4.4)") jLES,ncly1,nclyn,ny
       open(70,file=trim(filename),status='unknown',form='formatted')
       do j=1,ysize(2)
          y = real(j-1,mytype)*dy
          write(70,'(9E14.6)') y,&
               ffy2(1,j,1),fiffy2(1,j,1),&
               ffyp2(1,j,1),fiffyp2(1,j,1)
       enddo
       close(70)
    endif
    do k=1,zsize(3)
       z = real(k-1,mytype)*dz*4*pi
       do j=1,zsize(2)
          do i=1,zsize(1)
             fz3(i,j,k) = sin_prec(z)
             fzp3(i,j,k) = cos_prec(z)
          enddo
       enddo
    enddo
    call derz (dfdz3 ,fz3 ,di3,sz,ffz ,fsz ,fwz ,zsize(1),zsize(2),zsize(3),0)
    call derzz (dfdzz3 ,fz3 ,di3,sz,sfz ,ssz ,swz ,zsize(1),zsize(2),zsize(3),0)
    if (nrank.eq.0) then
       write(filename,"('schemes_z',I1.1,I1.1,I1.1,I4.4)") jLES,nclz1,nclzn,nz
       open(67,file=trim(filename),status='unknown',form='formatted')
       do k=1,zsize(3)
          z = real(k-1,mytype)*dz
          write(67,'(5E14.6)') z,&
               four*pi*cos_prec(four*pi*z),dfdz3(1,1,k),&
               -sixteen*pi*pi*sin_prec(four*pi*z),dfdzz3(1,1,k)
       enddo
       close(67)
    endif
    call derzvp(test3,fz3,di3,sz,cfz6,csz6,cwz6,zsize(1),zsize(2),zsize(3),nzmsize,0)
    call interzvp(test33,fzp3,di3,sz,cifzp6,ciszp6,ciwzp6,zsize(1),zsize(2),zsize(3),nzmsize,1)
    if (nrank.eq.0) then
       write(filename,"('schemes_zVP',I1.1,I1.1,I1.1,I4.4)") jLES,nclz1,nclzn,nz
       open(68,file=trim(filename),status='unknown',form='formatted')
       do k=1,nzmsize
          z1 = real(k-half,mytype)*dz
          write(68,'(5E14.6)') z1,&
               four*pi*cos_prec(four*pi*z1),test3(1,1,k),&
               cos_prec(four*pi*z1),test33(1,1,k)
       enddo
       close(68)
    endif

    do k=1,nzmsize
       z1 = real(k-half,mytype)*dz*4*pi
       do j=1,zsize(2)
          do i=1,zsize(1)
             test3(i,j,k) = cos_prec(z1)
          enddo
       enddo
    enddo
    call derzpv(fz3,test3,di3,sz,cfip6z,csip6z,cwip6z,cfz6,csz6,cwz6,zsize(1),zsize(2),nzmsize,zsize(3),1)
    call interzpv(fzp3,test3,di3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,zsize(1),zsize(2),nzmsize,zsize(3),1)
    if (nrank.eq.0) then
       write(filename,"('schemes_zPV',I1.1,I1.1,I1.1,I4.4)") jLES,nclz1,nclzn,nz
       open(69,file=trim(filename),status='unknown',form='formatted')
       do k=1,zsize(3)
          z = real(k-1,mytype)*dz
          write(69,'(5E14.6)') z,&
               -four*pi*sin_prec(four*pi*z),fz3(1,1,k),&
               cos_prec(four*pi*z),fzp3(1,1,k)
       enddo
       close(69)
    endif
    ! FILTER
    call filter(0.45_mytype)
    do k=1,zsize(3)
       z = real(k-1,mytype)*dz*four*pi
       do j=1,zsize(2)
          do i=1,zsize(1)
             ffz3(i,j,k)  = sin_prec(z)
             ffzp3(i,j,k) = cos_prec(z)
          enddo
       enddo
    enddo
    call filz (fiffz3  ,ffz3  ,di3,fisz,fiffz ,fifsz ,fifwz ,zsize(1),zsize(2),zsize(3),0)
    call filz (fiffzp3 ,ffzp3 ,di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1)
    if (nrank.eq.0) then
       write(filename,"('filter_z',I1.1,I1.1,I1.1,I4.4)") jLES,nclz1,nclzn,nz
       open(70,file=trim(filename),status='unknown',form='formatted')
       do k=1,zsize(3)
          z = real(k-1,mytype)*dz
          write(70,'(9E14.6)') z,&
               ffz3 (1,1,k),fiffz3 (1,1,k),&
               ffzp3(1,1,k),fiffzp3(1,1,k)
       enddo
       close(70)
    endif

    stop 'stop debug_schemes'
  end subroutine debug_schemes

  subroutine init_post(ep1)

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ep1

  end subroutine init_post

  subroutine postprocess_dbg(ux1,uy1,uz1,phi1,ep1)

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ep1

  end subroutine postprocess_dbg

  subroutine write_probes(ux1,uy1,uz1,phi1)

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ux1, uy1, uz1
    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),numscalar) :: phi1

  end subroutine write_probes

end module dbg_schemes
