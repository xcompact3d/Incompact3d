module decomp_2d_poisson

  use decomp_2d
  use decomp_2d_fft

  use param
  use variables

  implicit none

  private        ! Make everything private unless declared public

!  real(mytype), private, parameter :: PI = 3.14159265358979323846_mytype

#ifdef DOUBLE_PREC
  real(mytype), parameter :: epsilon = 1.e-16_mytype
#else
  real(mytype), parameter :: epsilon = 1.e-8_mytype
#endif

  ! boundary conditions
  integer, save :: bcx, bcy, bcz

  ! decomposition object for physical space
  TYPE(DECOMP_INFO), save :: ph
  
  ! decomposition object for spectral space
  TYPE(DECOMP_INFO), save :: sp

  ! store sine/cosine factors
  real(mytype), save, allocatable, dimension(:) :: az,bz
  real(mytype), save, allocatable, dimension(:) :: ay,by
  real(mytype), save, allocatable, dimension(:) :: ax,bx

  ! wave numbers
  complex(mytype), save, allocatable, dimension(:,:,:) :: kxyz
  !wave numbers for stretching in a pentadiagonal matrice
  complex(mytype), save, allocatable, dimension(:,:,:,:) :: a,a2,a3
  ! work arrays, 
  ! naming convention: cw (complex); rw (real); 
  !                    1 = X-pencil; 2 = Y-pencil; 3 = Z-pencil
  real(mytype), allocatable, dimension(:,:,:) :: rw1,rw1b,rw2,rw2b,rw3
  complex(mytype), allocatable, dimension(:,:,:) :: cw1,cw1b,cw2,cw22,cw2b,cw2c

  ! underlying FFT library only needs to be initialised once
  logical, save :: fft_initialised = .false.

  ABSTRACT INTERFACE
     SUBROUTINE poisson_xxx(rhs)
       use decomp_2d, only : mytype
       real(mytype), dimension(:,:,:), intent(INOUT) :: rhs
     END SUBROUTINE poisson_xxx
  END INTERFACE
  PROCEDURE (poisson_xxx), POINTER :: poisson

  public :: decomp_2d_poisson_init,decomp_2d_poisson_finalize,poisson
contains



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialise Poisson solver for given boundary conditions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_poisson_init()

    implicit none

    integer :: nx, ny, nz, i

    if (nclx) then
       bcx=0
    else
       bcx=1
    endif
    if (ncly) then
       bcy=0
    else
       bcy=1
    endif
    if (nclz) then
       bcz=0
    else
       bcz=1
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Top level wrapper
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (bcx==0 .and. bcy==0 .and. bcz==0) then
       poisson => poisson_000
    else if (bcx==1 .and. bcy==0 .and. bcz==0) then
       poisson => poisson_100
    else if (bcx==0 .and. bcy==1 .and. bcz==0) then
       poisson => poisson_010
    else if (bcx==1 .and. bcy==1) then   ! 110 & 111
       poisson => poisson_11x
    else
       stop 'boundary condition not supported'
    end if

    nx = nx_global
    ny = ny_global
    nz = nz_global

    ! pressure-grid having 1 fewer point for non-periodic directions
    if (bcx==1) nx=nx-1
    if (bcy==1) ny=ny-1
    if (bcz==1) nz=nz-1

#ifdef DEBG 
  if (nrank .eq. 0) print *,'# decomp_2d_poisson_init start'
#endif

    allocate(ax(nx),bx(nx))
    allocate(ay(ny),by(ny))
    allocate(az(nz),bz(nz))
    call abxyz(ax,ay,az,bx,by,bz,nx,ny,nz,bcx,bcy,bcz)

#ifdef DEBG 
  if (nrank .eq. 0) print *,'# decomp_2d_poisson_init decomp_info_init'
#endif

    call decomp_info_init(nx, ny, nz, ph)
    call decomp_info_init(nx, ny, nz/2+1, sp)

#ifdef DEBG 
  if (nrank .eq. 0) print *,'# decomp_2d_poisson_init decomp_info_init ok'
#endif

    ! allocate work space
    if (bcx==0 .and. bcy==0 .and. bcz==0) then
       allocate(cw1(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
       allocate(kxyz(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
       allocate(a(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a2(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a3(sp%yst(1):sp%yen(1),ny,sp%yst(3):sp%yen(3),5))
    else if (bcx==1 .and. bcy==0 .and. bcz==0) then
       allocate(cw1(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
                 sp%xst(3):sp%xen(3)))
       allocate(cw1b(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
                 sp%xst(3):sp%xen(3)))
       allocate(rw1(ph%xst(1):ph%xen(1),ph%xst(2):ph%xen(2), &
            ph%xst(3):ph%xen(3)))
       allocate(rw1b(ph%xst(1):ph%xen(1),ph%xst(2):ph%xen(2), &
            ph%xst(3):ph%xen(3)))
       allocate(rw2(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
       allocate(kxyz(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
       allocate(a(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a2(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a3(sp%yst(1):sp%yen(1),ny,sp%yst(3):sp%yen(3),5))
    else if (bcx==0 .and. bcy==1 .and. bcz==0) then
       allocate(rw2(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
       allocate(rw2b(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
       allocate(cw1(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
                 sp%xst(3):sp%xen(3)))
       allocate(cw2(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
                 sp%yst(3):sp%yen(3)))
       allocate(cw22(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
       allocate(cw2b(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
                 sp%yst(3):sp%yen(3)))
       allocate(cw2c(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
                 sp%yst(3):sp%yen(3)))
       allocate(kxyz(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
       allocate(a(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a2(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a3(sp%yst(1):sp%yen(1),ny,sp%yst(3):sp%yen(3),5))
    else if (bcx==1 .and. bcy==1) then
       allocate(cw1(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
                 sp%xst(3):sp%xen(3)))
       allocate(cw1b(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
                 sp%xst(3):sp%xen(3)))
       allocate(cw2(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
                 sp%yst(3):sp%yen(3)))
       allocate(cw22(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
       allocate(cw2b(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
                 sp%yst(3):sp%yen(3)))
       allocate(cw2c(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
                 sp%yst(3):sp%yen(3)))
       allocate(rw1(ph%xst(1):ph%xen(1),ph%xst(2):ph%xen(2), &
            ph%xst(3):ph%xen(3)))
       allocate(rw1b(ph%xst(1):ph%xen(1),ph%xst(2):ph%xen(2), &
            ph%xst(3):ph%xen(3)))
       allocate(rw2(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
       allocate(rw2b(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
       if (bcz==1) then  
          allocate(rw3(ph%zsz(1),ph%zsz(2),ph%zsz(3)))
       end if
       allocate(kxyz(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))    
       allocate(a(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a2(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a3(sp%yst(1):sp%yen(1),nym,sp%yst(3):sp%yen(3),5))      
    end if 

#ifdef DEBG 
  if (nrank .eq. 0) print *,'# decomp_2d_poisson_init before waves'
#endif

    call waves()

#ifdef DEBG 
  if (nrank .eq. 0) print *,'# decomp_2d_poisson_init end'
#endif

    return
  end subroutine decomp_2d_poisson_init


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Release memory used by Poisson solver
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_poisson_finalize

    implicit none

    deallocate(ax,bx,ay,by,az,bz)

    call decomp_info_finalize(ph)
    call decomp_info_finalize(sp)

    call decomp_2d_fft_finalize
    fft_initialised = .false.

    deallocate(kxyz)

    if (bcx==0 .and. bcy==0 .and. bcz==0) then
       deallocate(cw1)
       deallocate(a,a2,a3)
    else if (bcx==1 .and. bcy==0 .and. bcz==0) then
       deallocate(cw1,cw1b,rw1,rw1b,rw2)
       deallocate(a,a2,a3)
    else if (bcx==0 .and. bcy==1 .and. bcz==0) then
       deallocate(cw1,cw2,cw2b,rw2,rw2b)
       deallocate(a,a2,a3)
    else if (bcx==1 .and. bcy==1) then
       deallocate(cw1,cw1b,cw2,cw2b,rw1,rw1b,rw2,rw2b)
       deallocate(a,a2,a3)
       if (bcz==1) then
          deallocate(rw3)
       end if
    end if

    return
  end subroutine decomp_2d_poisson_finalize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Solving 3D Poisson equation with periodic B.C in all 3 dimensions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine poisson_000(rhs)

    use derivX
    use derivY
    use derivZ

    ! right-hand-side of Poisson as input
    ! solution of Poisson as output
    real(mytype), dimension(:,:,:), intent(INOUT) :: rhs

    integer, dimension(3) :: fft_start, fft_end, fft_size

    complex(mytype) :: xyzk

    complex(mytype) :: ytt,xtt,ztt,yt1,xt1,yt2,xt2
    complex(mytype) :: xtt1,ytt1,ztt1,zt1,zt2


    real(mytype) :: tmp1, tmp2,x ,y, z
    
    integer :: nx,ny,nz, i,j,k

    nx = nx_global
    ny = ny_global
    nz = nz_global

    if (.not. fft_initialised) then
       call decomp_2d_fft_init(PHYSICAL_IN_Z)
       fft_initialised = .true.
    end if

    ! compute r2c transform 
    call decomp_2d_fft_3d(rhs,cw1)

    ! normalisation
    cw1 = cw1 / real(nx, kind=mytype) /real(ny, kind=mytype) &
         / real(nz, kind=mytype)

    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)

             ! post-processing in spectral space

             ! POST PROCESSING IN Z
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*bz(k)+tmp2*az(k), &
                  tmp2*bz(k)-tmp1*az(k), kind=mytype)

             ! POST PROCESSING IN Y
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*by(j)+tmp2*ay(j), &
                  tmp2*by(j)-tmp1*ay(j), kind=mytype)
             if (j.gt.(ny/2+1)) cw1(i,j,k)=-cw1(i,j,k)

             ! POST PROCESSING IN X
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*bx(i)+tmp2*ax(i), &
                  tmp2*bx(i)-tmp1*ax(i), kind=mytype)
             if (i.gt.(nx/2+1)) cw1(i,j,k)=-cw1(i,j,k)

             ! Solve Poisson
             tmp1=real(kxyz(i,j,k), kind=mytype)
             tmp2=aimag(kxyz(i,j,k))
             ! CANNOT DO A DIVISION BY ZERO
             if ((tmp1.lt.epsilon).or.(tmp2.lt.epsilon)) then
                cw1(i,j,k)=0._mytype
!                print *,'DIV 0',i,j,k,epsilon
             else
                cw1(i,j,k)=cmplx( real(cw1(i,j,k), kind=mytype) / (-tmp1), &
                     aimag(cw1(i,j,k))/(-tmp2), kind=mytype)
             end if

           !Print result in spectal space after Poisson
      !     if (abs(out(i,j,k)) > 1.0e-4) then
      !        write(*,*) 'AFTER',i,j,k,out(i,j,k),xyzk
      !     end if

             ! post-processing backward
             
             ! POST PROCESSING IN Z
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*bz(k)-tmp2*az(k), &
                  -tmp2*bz(k)-tmp1*az(k), kind=mytype)

             ! POST PROCESSING IN Y
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*by(j)+tmp2*ay(j), &
                  tmp2*by(j)-tmp1*ay(j), kind=mytype)
             if (j.gt.(ny/2+1)) cw1(i,j,k)=-cw1(i,j,k)

             ! POST PROCESSING IN X
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*bx(i)+tmp2*ax(i), &
                  -tmp2*bx(i)+tmp1*ax(i), kind=mytype)
             if (i.gt.(nx/2+1)) cw1(i,j,k)=-cw1(i,j,k)

          end do
       end do
    end do
             
    ! compute c2r transform
    call decomp_2d_fft_3d(cw1,rhs)
    
 !   call decomp_2d_fft_finalize

    return
  end subroutine poisson_000


  subroutine poisson_100(rhs)

    implicit none

    real(mytype), dimension(:,:,:), intent(INOUT) :: rhs

    complex(mytype) :: xyzk
    real(mytype) :: tmp1, tmp2, tmp3, tmp4
    real(mytype) :: xx1,xx2,xx3,xx4,xx5,xx6,xx7,xx8
    
    integer :: nx,ny,nz, i,j,k, itmp

100 format(1x,a8,3I4,2F12.6)

    nx = nx_global - 1
    ny = ny_global
    nz = nz_global

    ! rhs is in Z-pencil but requires global operations in X
    call transpose_z_to_y(rhs,rw2,ph)
    call transpose_y_to_x(rw2,rw1,ph)
    do k=ph%xst(3),ph%xen(3)
       do j=ph%xst(2),ph%xen(2)
          do i=1,nx/2
             rw1b(i,j,k)=rw1(2*(i-1)+1,j,k)
          enddo
          do i=nx/2+1,nx
             rw1b(i,j,k)=rw1(2*nx-2*i+2,j,k)
          enddo
       enddo
    end do

    call transpose_x_to_y(rw1b,rw2,ph)
    call transpose_y_to_z(rw2,rhs,ph)

    if (.not. fft_initialised) then
       call decomp_2d_fft_init(PHYSICAL_IN_Z,nx,ny,nz)
       fft_initialised = .true.
    end if

    ! compute r2c transform 
    call decomp_2d_fft_3d(rhs,cw1)

    ! normalisation
    cw1 = cw1 / real(nx, kind=mytype) /real(ny, kind=mytype) &
         / real(nz, kind=mytype)
#ifdef DEBUG
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             if (abs(cw1(i,j,k)) > 1.0e-4) then
                write(*,100) 'START',i,j,k,cw1(i,j,k)
             end if
          end do
       end do
    end do
#endif

    ! post-processing in spectral space

    ! POST PROCESSING IN Z
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*bz(k)+tmp2*az(k), &
                  tmp2*bz(k)-tmp1*az(k), kind=mytype)
#ifdef DEBUG
             if (abs(cw1(i,j,k)) > 1.0e-4) &
                  write(*,100) 'after z',i,j,k,cw1(i,j,k)
#endif
          end do
       end do
    end do

    ! POST PROCESSING IN Y
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*by(j)+tmp2*ay(j), &
                  tmp2*by(j)-tmp1*ay(j), kind=mytype)
             if (j.gt.(ny/2+1)) cw1(i,j,k)=-cw1(i,j,k)
#ifdef DEBUG
             if (abs(cw1(i,j,k)) > 1.0e-4) &
                  write(*,100) 'after y',i,j,k,cw1(i,j,k)
#endif
          end do
       end do
    end do

    ! POST PROCESSING IN X
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          cw1b(1,j,k)=cw1(1,j,k)
          do i = 2,nx
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             tmp3 = real(cw1(nx-i+2,j,k), kind=mytype)
             tmp4 = aimag(cw1(nx-i+2,j,k))
             xx1=tmp1*bx(i)/2._mytype
             xx2=tmp1*ax(i)/2._mytype
             xx3=tmp2*bx(i)/2._mytype
             xx4=tmp2*ax(i)/2._mytype
             xx5=tmp3*bx(i)/2._mytype
             xx6=tmp3*ax(i)/2._mytype
             xx7=tmp4*bx(i)/2._mytype
             xx8=tmp4*ax(i)/2._mytype
             cw1b(i,j,k) = cmplx(xx1+xx4+xx5-xx8,-xx2+xx3+xx6+xx7, &
                  kind=mytype)  
          end do
       end do
    end do
#ifdef DEBUG
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             if (abs(cw1b(i,j,k)) > 1.0e-4) then
                write(*,100) 'after x',i,j,k,cw1b(i,j,k)
             end if
          end do
       end do
    end do
#endif

    ! Solve Poisson
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             !tmp1=real(zk2(k)+yk2(j)+xk2(i), kind=mytype)
             !tmp2=aimag(zk2(k)+yk2(j)+xk2(i))
             tmp1=real(kxyz(i,j,k), kind=mytype)
             tmp2=aimag(kxyz(i,j,k))
             !xyzk=cmplx(tmp1,tmp2, kind=mytype)
             ! CANNOT DO A DIVISION BY ZERO
             if ((abs(tmp1).lt.epsilon).and.(abs(tmp2).lt.epsilon)) then    
                cw1b(i,j,k)=cmplx(0._mytype,0._mytype, kind=mytype)
             end if
             if ((abs(tmp1).lt.epsilon).and.(abs(tmp2).ge.epsilon)) then
                cw1b(i,j,k)=cmplx(0._mytype, &
                     aimag(cw1b(i,j,k))/(-tmp2), kind=mytype)
             end if
             if ((abs(tmp1).ge.epsilon).and.(abs(tmp2).lt.epsilon)) then    
                cw1b(i,j,k)=cmplx( real(cw1b(i,j,k), kind=mytype) &
                     /(-tmp1), 0._mytype, kind=mytype)
             end if
             if ((abs(tmp1).ge.epsilon).and.(abs(tmp2).ge.epsilon)) then
                cw1b(i,j,k)=cmplx( real(cw1b(i,j,k), kind=mytype) &
                     /(-tmp1), &
                     aimag(cw1b(i,j,k))/(-tmp2), kind=mytype)
             end if
#ifdef DEBUG
             if (abs(cw1b(i,j,k)) > 1.0e-4) &
                  write(*,100) 'AFTER',i,j,k,cw1b(i,j,k)
#endif
          end do
       end do
    end do
    
    ! post-processing backward

    ! POST PROCESSING IN X
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          cw1(1,j,k)=cw1b(1,j,k)
          do i = 2,nx 
             tmp1 = real(cw1b(i,j,k), kind=mytype)
             tmp2 = aimag(cw1b(i,j,k))
             tmp3 = real(cw1b(nx-i+2,j,k), kind=mytype)
             tmp4 = aimag(cw1b(nx-i+2,j,k))
             xx1=tmp1*bx(i)
             xx2=tmp1*ax(i)
             xx3=tmp2*bx(i)
             xx4=tmp2*ax(i)
             xx5=tmp3*bx(i)
             xx6=tmp3*ax(i)
             xx7=tmp4*bx(i)
             xx8=tmp4*ax(i)
             cw1(i,j,k) = cmplx(xx1-xx4+xx6+xx7,-(-xx2-xx3+xx5-xx8), &
                  kind=mytype)        
          end do
       end do
    end do
#ifdef DEBUG
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             if (abs(cw1(i,j,k)) > 1.0e-4) then
                write(*,100) 'AFTER X',i,j,k,cw1(i,j,k)
             end if
          end do
       end do
    end do
#endif

    ! POST PROCESSING IN Y
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*by(j)-tmp2*ay(j), &
                  tmp2*by(j)+tmp1*ay(j), kind=mytype)
             if (j.gt.(ny/2+1)) cw1(i,j,k)=-cw1(i,j,k)
#ifdef DEBUG
             if (abs(cw1(i,j,k)) > 1.0e-4) &
                  write(*,100) 'AFTER Y',i,j,k,cw1(i,j,k)
#endif
          end do
       end do
    end do

    ! POST PROCESSING IN Z
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*bz(k)-tmp2*az(k), &
                  tmp2*bz(k)+tmp1*az(k), kind=mytype)
#ifdef DEBUG
             if (abs(cw1(i,j,k)) > 1.0e-4) &
                  write(*,100) 'END',i,j,k,cw1(i,j,k)
#endif
          end do
       end do
    end do

    ! compute c2r transform
    call decomp_2d_fft_3d(cw1,rhs)

    ! rhs is in Z-pencil but requires global operations in X
    call transpose_z_to_y(rhs,rw2,ph)
    call transpose_y_to_x(rw2,rw1,ph)
    do k=ph%xst(3),ph%xen(3)
       do j=ph%xst(2),ph%xen(2)
          do i=1,nx/2
             rw1b(2*i-1,j,k)=rw1(i,j,k)
          enddo
          do i=1,nx/2
             rw1b(2*i,j,k)=rw1(nx-i+1,j,k)
          enddo
       enddo
    end do
    call transpose_x_to_y(rw1b,rw2,ph)
    call transpose_y_to_z(rw2,rhs,ph)
    
  !  call decomp_2d_fft_finalize

    return
  end subroutine poisson_100

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Solving 3D Poisson equation: Neumann in Y; periodic in X & Z
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine poisson_010(rhs)

    implicit none

    real(mytype), dimension(:,:,:), intent(INOUT) :: rhs

    complex(mytype) :: xyzk
    real(mytype) :: tmp1, tmp2, tmp3, tmp4
    real(mytype) :: xx1,xx2,xx3,xx4,xx5,xx6,xx7,xx8

    integer :: nx,ny,nz, i,j,k

100 format(1x,a8,3I4,2F12.6)

    nx = nx_global
    ny = ny_global - 1
    nz = nz_global

    ! rhs is in Z-pencil but requires global operations in Y
    call transpose_z_to_y(rhs,rw2,ph)
    do k=ph%yst(3),ph%yen(3)
       do i=ph%yst(1),ph%yen(1)
          do j=1,ny/2
             rw2b(i,j,k)=rw2(i,2*(j-1)+1,k)
          enddo
          do j=ny/2+1,ny
             rw2b(i,j,k)=rw2(i,2*ny-2*j+2,k)
          enddo
       enddo
    end do
    call transpose_y_to_z(rw2b,rhs,ph)

    if (.not. fft_initialised) then
       call decomp_2d_fft_init(PHYSICAL_IN_Z,nx,ny,nz)
       fft_initialised = .true.
    end if
    ! compute r2c transform 
    call decomp_2d_fft_3d(rhs,cw1)

    ! normalisation
    cw1 = cw1 / real(nx, kind=mytype) /real(ny, kind=mytype) &
         / real(nz, kind=mytype)
#ifdef DEBUG
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             if (abs(cw1(i,j,k)) > 1.0e-4) then
                write(*,100) 'START',i,j,k,cw1(i,j,k)
             end if
          end do
       end do
    end do
#endif

    ! post-processing in spectral space

    ! POST PROCESSING IN Z
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*bz(k)+tmp2*az(k), &
                  tmp2*bz(k)-tmp1*az(k), kind=mytype)
#ifdef DEBUG
             if (abs(cw1(i,j,k)) > 1.0e-4) &
                  write(*,100) 'after z',i,j,k,cw1(i,j,k)
#endif
          end do
       end do
    end do

    ! POST PROCESSING IN X
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*bx(i)+tmp2*ax(i), &
                  tmp2*bx(i)-tmp1*ax(i), kind=mytype)
             if (i.gt.(nx/2+1)) cw1(i,j,k)=-cw1(i,j,k)
#ifdef DEBUG
             if (abs(cw1(i,j,k)) > 1.0e-4) &
                  write(*,100) 'after x',i,j,k,cw1(i,j,k)
#endif
          end do
       end do
    end do
    
    ! POST PROCESSING IN Y
    ! NEED TO BE IN Y PENCILS!!!!!!!!!!!!!!!
    call transpose_x_to_y(cw1,cw2,sp)

    do k = sp%yst(3), sp%yen(3)
       do i = sp%yst(1), sp%yen(1)
          cw2b(i,1,k)=cw2(i,1,k)
          do j = 2,ny      
             tmp1 = real(cw2(i,j,k), kind=mytype)
             tmp2 = aimag(cw2(i,j,k))
             tmp3 = real(cw2(i,ny-j+2,k), kind=mytype)
             tmp4 = aimag(cw2(i,ny-j+2,k))
             xx1=tmp1*by(j)/2._mytype
             xx2=tmp1*ay(j)/2._mytype
             xx3=tmp2*by(j)/2._mytype
             xx4=tmp2*ay(j)/2._mytype
             xx5=tmp3*by(j)/2._mytype
             xx6=tmp3*ay(j)/2._mytype
             xx7=tmp4*by(j)/2._mytype
             xx8=tmp4*ay(j)/2._mytype
             cw2b(i,j,k) = cmplx(xx1+xx4+xx5-xx8,-xx2+xx3+xx6+xx7, &
                  kind=mytype)
          end do
       end do
    end do
#ifdef DEBUG
    do k = sp%yst(3), sp%yen(3)
       do j = sp%yst(2), sp%yen(2)
          do i = sp%yst(1), sp%yen(1)
             if (abs(cw2b(i,j,k)) > 1.0e-4) then
                write(*,100) 'after y',i,j,k,cw2b(i,j,k)
                print *,kxyz(i,j,k)
             end if
          end do
       end do
    end do
#endif

    if (istret==0) then 

    ! Solve Poisson
    ! doing wave number division in Y-pencil
    do k = sp%yst(3), sp%yen(3)
       do j = sp%yst(2), sp%yen(2)
          do i = sp%yst(1), sp%yen(1)
             !tmp1=real(zk2(k)+yk2(j)+xk2(i), kind=mytype)
             !tmp2=aimag(zk2(k)+yk2(j)+xk2(i))
             tmp1=real(kxyz(i,j,k), kind=mytype)
             tmp2=aimag(kxyz(i,j,k))
             !xyzk=cmplx(tmp1,tmp2, kind=mytype)
             !CANNOT DO A DIVISION BY ZERO
             if ((abs(tmp1).lt.epsilon).and.(abs(tmp2).lt.epsilon)) then    
                cw2b(i,j,k)=cmplx(0._mytype,0._mytype, kind=mytype)
             end if
             if ((abs(tmp1).lt.epsilon).and.(abs(tmp2).ge.epsilon)) then
                cw2b(i,j,k)=cmplx(0._mytype, &
                     aimag(cw2b(i,j,k))/(-tmp2), kind=mytype)
             end if
             if ((abs(tmp1).ge.epsilon).and.(abs(tmp2).lt.epsilon)) then    
                cw2b(i,j,k)=cmplx( real(cw2b(i,j,k), kind=mytype) &
                     /(-tmp1), 0._mytype, kind=mytype)
             end if
             if ((abs(tmp1).ge.epsilon).and.(abs(tmp2).ge.epsilon)) then
                cw2b(i,j,k)=cmplx( real(cw2b(i,j,k), kind=mytype) &
                     /(-tmp1), &
                     aimag(cw2b(i,j,k))/(-tmp2), kind=mytype)
             end if
          end do
       end do
    end do

    else
       
       call matrice_refinement()
!       do k = sp%yst(3), sp%yen(3)
!          do j = 1,ny/2
!             do i = sp%yst(1), sp%yen(1)
!                print *,i,j,k,a(i,j,k,3)
!!                if (nrank.le.1) print *,i,j,k,a(i,j,k,3)
!!                if (nrank.gt.1) print *,i+4,j,k,a(i,j,k,3)
!             enddo
!          enddo
!       enddo
     

       if (istret.ne.3) then
          cw2(:,:,:)=0.;cw2c(:,:,:)=0.
          do k = sp%yst(3), sp%yen(3)
          do j = 1,ny/2
          do i = sp%yst(1), sp%yen(1)
             cw2(i,j,k)=cw2b(i,2*j-1,k) 
             cw2c(i,j,k)=cw2b(i,2*j,k)
          enddo
          enddo
          enddo

 !   do k = sp%yst(3), sp%yen(3)
 !      do j = 1,ny/2
 !         do i = sp%yst(1), sp%yen(1)
 !            if (abs(cw2(i,j,k)) > 1.0e-4) then
 !               write(*,*) 'before IN',i,j,k,cw2(i,j,k)!*2.
 !!            end if
  !        end do
  !     end do
  !  end do
          
         call inversion5_v1(a,cw2,sp)
         call inversion5_v1(a2,cw2c,sp)

!         cw2(1,1,1)=cw2(1,1,1)*0.5
         
 
!   do k = sp%yst(3), sp%yen(3)
!       do j = 1,ny/2
!          do i = sp%yst(1), sp%yen(1)
!             if (abs(cw2c(i,j,k)) > 1.0e-4) then
!                write(*,*) 'after IN',i,j,k,cw2c(i,j,k)!*2.
!             end if
!          end do
!       end do
!    end do

          cw2b(:,:,:)=0.
          do k=sp%yst(3), sp%yen(3)
          do j=1,ny-1,2
          do i=sp%yst(1), sp%yen(1)
             cw2b(i,j,k)=cw2(i,(j+1)/2,k)
          enddo
          enddo
          do j=2,ny,2
          do i=sp%yst(1), sp%yen(1)
             cw2b(i,j,k)=cw2c(i,j/2,k)
          enddo
          enddo
          enddo
          !do k=sp%yst(3), sp%yen(3)
          !do i=sp%yst(1), sp%yen(1)
          !   if ((xkx(i)==0).and.(zkz(k)==0)) then
          !   !   cw2b(i,1,1)=0.
          !   !   cw2b(i,ny,1)=0.
          !   endif
          !enddo
          !enddo
       else
          do k = sp%yst(3), sp%yen(3)
          do j = 1,ny
          do i = sp%yst(1), sp%yen(1)
             cw2(i,j,k)=cw2b(i,j,k) 
          enddo
          enddo
          enddo
          call inversion5_v2(a3,cw2,sp)
          do k = sp%yst(3), sp%yen(3)
          do j = 1,ny
          do i = sp%yst(1), sp%yen(1)
             cw2b(i,j,k)=cw2(i,j,k) 
          enddo
          enddo
          enddo
       endif

    endif

!    print *,nrank, sp%yst(3),sp%yen(3),sp%yst(1),sp%yen(1)

!we are in Y pencil
    do k = sp%yst(3), sp%yen(3)  
    do i = sp%yst(1), sp%yen(1)
       if ((i==nx/2+1).and.(k==nz/2+1)) then
          cw2b(i,:,k)=0.
       endif
    enddo
    enddo
#ifdef DEBUG
    do k = sp%yst(3), sp%yen(3)
       do j = sp%yst(2), sp%yen(2)
          do i = sp%yst(1), sp%yen(1)
             if (abs(cw2b(i,j,k)) > 1.0e-4) then
                write(*,100) 'AFTER',i,j,k,cw2b(i,j,k)
                print *,kxyz(i,j,k)
             end if
          end do
       end do
    end do
#endif

    ! post-processing backward

    ! POST PROCESSING IN Y
    do k = sp%yst(3), sp%yen(3)
       do i = sp%yst(1), sp%yen(1)
          cw2(i,1,k)=cw2b(i,1,k)
          do j = 2,ny 
             tmp1 = real(cw2b(i,j,k), kind=mytype)
             tmp2 = aimag(cw2b(i,j,k))
             tmp3 = real(cw2b(i,ny-j+2,k), kind=mytype)
             tmp4 = aimag(cw2b(i,ny-j+2,k))
             xx1=tmp1*by(j)
             xx2=tmp1*ay(j)
             xx3=tmp2*by(j)
             xx4=tmp2*ay(j)
             xx5=tmp3*by(j)
             xx6=tmp3*ay(j)
             xx7=tmp4*by(j)
             xx8=tmp4*ay(j)
             cw2(i,j,k) = cmplx(xx1-xx4+xx6+xx7,-(-xx2-xx3+xx5-xx8), &
                  kind=mytype) 
          end do
       end do
    end do
           
    ! Back to X-pencil
    call transpose_y_to_x(cw2,cw1,sp)
#ifdef DEBUG
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             if (abs(cw1(i,j,k)) > 1.0e-4) then
                write(*,100) 'AFTER Y',i,j,k,cw1(i,j,k)
             end if
          end do
       end do
    end do
#endif
    
    ! POST PROCESSING IN X
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*bx(i)-tmp2*ax(i), &
                  tmp2*bx(i)+tmp1*ax(i), kind=mytype)
             if (i.gt.(nx/2+1)) cw1(i,j,k)=-cw1(i,j,k)
#ifdef DEBUG
             if (abs(cw1(i,j,k)) > 1.0e-4) &
                  write(*,100) 'AFTER X',i,j,k,cw1(i,j,k)
#endif
          end do
       end do
    end do

    ! POST PROCESSING IN Z
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*bz(k)-tmp2*az(k), &
                  tmp2*bz(k)+tmp1*az(k), kind=mytype)
#ifdef DEBUG
             if (abs(cw1(i,j,k)) > 1.0e-4) &
                  write(*,100) 'END',i,j,k,cw1(i,j,k)
#endif
          end do
       end do
    end do

    ! compute c2r transform, back to physical space
    call decomp_2d_fft_3d(cw1,rhs)

    ! rhs is in Z-pencil but requires global operations in Y
    call transpose_z_to_y(rhs,rw2,ph)
    do k=ph%yst(3),ph%yen(3)
       do i=ph%yst(1),ph%yen(1)
          do j=1,ny/2
             rw2b(i,2*j-1,k)=rw2(i,j,k)
          enddo
          do j=1,ny/2
             rw2b(i,2*j,k)=rw2(i,ny-j+1,k)
          enddo
       enddo
    end do
    call transpose_y_to_z(rw2b,rhs,ph)

  !  call decomp_2d_fft_finalize

    return
  end subroutine poisson_010


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Solving 3D Poisson equation: Neumann in X, Y; Neumann/periodic in Z
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine poisson_11x(rhs)

    implicit none

    real(mytype), dimension(:,:,:), intent(INOUT) :: rhs

    complex(mytype) :: xyzk
    real(mytype) :: tmp1, tmp2, tmp3, tmp4
    real(mytype) :: xx1,xx2,xx3,xx4,xx5,xx6,xx7,xx8

    integer :: nx,ny,nz, i,j,k

100 format(1x,a8,3I4,2F12.6)

    nx = nx_global - 1
    ny = ny_global - 1

    if (bcz==1) then	
       nz = nz_global - 1
    else if (bcz==0) then
       nz = nz_global
    end if
   
    if (bcz==1) then  
       do j=1,ph%zsz(2)
          do i=1,ph%zsz(1)
             do k=1,nz/2
                rw3(i,j,k)=rhs(i,j,2*(k-1)+1)
             end do
             do k=nz/2+1,nz
                rw3(i,j,k)=rhs(i,j,2*nz-2*k+2)
             end do
          end do
       end do
       call transpose_z_to_y(rw3,rw2,ph)
    else if (bcz==0) then     
       call transpose_z_to_y(rhs,rw2,ph)
    end if
    

    do k=ph%yst(3),ph%yen(3)
       do i=ph%yst(1),ph%yen(1)
          do j=1,ny/2
             rw2b(i,j,k)=rw2(i,2*(j-1)+1,k)
          end do
          do j=ny/2+1,ny
             rw2b(i,j,k)=rw2(i,2*ny-2*j+2,k)
          end do
       end do
    end do

    ! the global operations in X
    call transpose_y_to_x(rw2b,rw1,ph)

    do k=ph%xst(3),ph%xen(3)
       do j=ph%xst(2),ph%xen(2)
          do i=1,nx/2
             rw1b(i,j,k)=rw1(2*(i-1)+1,j,k)
          end do
          do i=nx/2+1,nx
             rw1b(i,j,k)=rw1(2*nx-2*i+2,j,k)
          end do
       end do
    end do

    ! back to Z-pencil
    call transpose_x_to_y(rw1b,rw2,ph)
    call transpose_y_to_z(rw2,rhs,ph)

    if (.not. fft_initialised) then
       call decomp_2d_fft_init(PHYSICAL_IN_Z,nx,ny,nz)
       fft_initialised = .true.
    end if

    ! compute r2c transform 

    call decomp_2d_fft_3d(rhs,cw1)

    ! normalisation
    cw1 = cw1 / real(nx, kind=mytype) /real(ny, kind=mytype) &
         / real(nz, kind=mytype)
#ifdef DEBUG
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             if (abs(cw1(i,j,k)) > 1.0e-4) then
                write(*,100) 'START',i,j,k,cw1(i,j,k)
             end if
          end do
       end do
    end do
#endif

    ! post-processing in spectral space

    ! POST PROCESSING IN Z
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*bz(k)+tmp2*az(k), &
                  tmp2*bz(k)-tmp1*az(k), kind=mytype)
#ifdef DEBUG
             if (abs(cw1(i,j,k)) > 1.0e-4) &
                  write(*,100) 'after z',i,j,k,cw1(i,j,k)
#endif
          end do
       end do
    end do

    ! POST PROCESSING IN Y
    ! WE HAVE TO BE IN Y PENCILS
    call transpose_x_to_y(cw1,cw2,sp)
    do k = sp%yst(3), sp%yen(3)
       do i = sp%yst(1), sp%yen(1)
          cw2b(i,1,k)=cw2(i,1,k)
          do j = 2,ny 
             tmp1 = real(cw2(i,j,k), kind=mytype)
             tmp2 = aimag(cw2(i,j,k))
             tmp3 = real(cw2(i,ny-j+2,k), kind=mytype)
             tmp4 = aimag(cw2(i,ny-j+2,k))
             xx1=tmp1*by(j)/2._mytype
             xx2=tmp1*ay(j)/2._mytype
             xx3=tmp2*by(j)/2._mytype
             xx4=tmp2*ay(j)/2._mytype
             xx5=tmp3*by(j)/2._mytype
             xx6=tmp3*ay(j)/2._mytype
             xx7=tmp4*by(j)/2._mytype
             xx8=tmp4*ay(j)/2._mytype 
             cw2b(i,j,k) = cmplx(xx1+xx4+xx5-xx8,-xx2+xx3+xx6+xx7, &
                  kind=mytype)  
          end do
       end do
    end do
    
    ! back to X-pencil
    call transpose_y_to_x(cw2b,cw1,sp)
#ifdef DEBUG
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             if (abs(cw1(i,j,k)) > 1.0e-4) then
                write(*,100) 'after y',i,j,k,cw1(i,j,k)
             end if
          end do
       end do
    end do
#endif
    
    ! POST PROCESSING IN X
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          cw1b(1,j,k)=cw1(1,j,k)
          do i = 2,nx
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             tmp3 = real(cw1(nx-i+2,j,k), kind=mytype)
             tmp4 = aimag(cw1(nx-i+2,j,k))
             xx1=tmp1*bx(i)/2._mytype
             xx2=tmp1*ax(i)/2._mytype
             xx3=tmp2*bx(i)/2._mytype
             xx4=tmp2*ax(i)/2._mytype
             xx5=tmp3*bx(i)/2._mytype
             xx6=tmp3*ax(i)/2._mytype
             xx7=tmp4*bx(i)/2._mytype
             xx8=tmp4*ax(i)/2._mytype
             cw1b(i,j,k) = cmplx(xx1+xx4+xx5-xx8,-xx2+xx3+xx6+xx7, &
                  kind=mytype)  
          end do
       end do
    end do
    
#ifdef DEBUG
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             if (abs(cw1b(i,j,k)) > 1.0e-4) then
                write(*,*) 'BEFORE',i,j,k,cw1b(i,j,k)
             end if
          end do
       end do
    end do
#endif

    if (istret==0) then

    ! Solve Poisson
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             !tmp1=real(zk2(k)+yk2(j)+xk2(i), kind=mytype)
             !tmp2=aimag(zk2(k)+yk2(j)+xk2(i))
             tmp1=real(kxyz(i,j,k), kind=mytype)
             tmp2=aimag(kxyz(i,j,k))
             !xyzk=cmplx(tmp1,tmp2, kind=mytype)
             !CANNOT DO A DIVISION BY ZERO
             if ((abs(tmp1).lt.epsilon).and.(abs(tmp2).lt.epsilon)) then    
                cw1b(i,j,k)=cmplx(0._mytype,0._mytype, kind=mytype)
             end if
             if ((abs(tmp1).lt.epsilon).and.(abs(tmp2).ge.epsilon)) then
                cw1b(i,j,k)=cmplx(0._mytype, &
                     aimag(cw1b(i,j,k))/(-tmp2), kind=mytype)
             end if
             if ((abs(tmp1).ge.epsilon).and.(abs(tmp2).lt.epsilon)) then    
                cw1b(i,j,k)=cmplx( real(cw1b(i,j,k), kind=mytype) &
                     /(-tmp1), 0._mytype, kind=mytype)
             end if
             if ((abs(tmp1).ge.epsilon).and.(abs(tmp2).ge.epsilon)) then
                cw1b(i,j,k)=cmplx( real(cw1b(i,j,k), kind=mytype) &
                     /(-tmp1), &
                     aimag(cw1b(i,j,k))/(-tmp2), kind=mytype)
             end if
          end do
       end do
    end do

    else
       call matrice_refinement()
! the stretching is only working in Y pencils
       call transpose_x_to_y(cw1b,cw2b,sp)
       !we are now in Y pencil
       
       if (istret.ne.3) then
          cw2(:,:,:)=0.;cw2c(:,:,:)=0.
          do k = sp%yst(3), sp%yen(3)
          do j = 1,ny/2
          do i = sp%yst(1), sp%yen(1)
             cw2(i,j,k)=cw2b(i,2*j-1,k) 
             cw2c(i,j,k)=cw2b(i,2*j,k)
          enddo
          enddo
          enddo
          call inversion5_v1(a,cw2,sp)
          call inversion5_v1(a2,cw2c,sp)

          cw2b(:,:,:)=0.
          do k=sp%yst(3), sp%yen(3)
          do j=1,ny-1,2
          do i=sp%yst(1), sp%yen(1)
             cw2b(i,j,k)=cw2(i,(j+1)/2,k)
          enddo
          enddo
          do j=2,ny,2
          do i=sp%yst(1), sp%yen(1)
             cw2b(i,j,k)=cw2c(i,j/2,k)
          enddo
          enddo
          enddo
       else
          cw2(:,:,:)=0.
          do k = sp%yst(3), sp%yen(3)
          do j = sp%yst(2), sp%yen(2)
          do i = sp%yst(1), sp%yen(1)
             cw2(i,j,k)=cw2b(i,j,k) 
          enddo
          enddo
          enddo

          call inversion5_v2(a3,cw2,sp)

          do k = sp%yst(3), sp%yen(3)
          do j = sp%yst(2), sp%yen(2)
          do i = sp%yst(1), sp%yen(1)
             cw2b(i,j,k)=cw2(i,j,k) 
          enddo
          enddo
          enddo
       endif
!we have to go back in X pencils
       call transpose_y_to_x(cw2b,cw1b,sp)
    endif

#ifdef DEBUG
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             if (abs(cw1b(i,j,k)) > 1.0e-6) then
                write(*,*) 'AFTER',i,j,k,cw1b(i,j,k)
             end if
          end do
       end do
    end do
#endif
!stop
    ! post-processing backward
    
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          cw1(1,j,k)=cw1b(1,j,k)
          do i = 2,nx
             tmp1 = real(cw1b(i,j,k), kind=mytype)
             tmp2 = aimag(cw1b(i,j,k))
             tmp3 = real(cw1b(nx-i+2,j,k), kind=mytype)
             tmp4 = aimag(cw1b(nx-i+2,j,k))
             xx1=tmp1*bx(i)
             xx2=tmp1*ax(i)
             xx3=tmp2*bx(i)
             xx4=tmp2*ax(i)
             xx5=tmp3*bx(i)
             xx6=tmp3*ax(i)
             xx7=tmp4*bx(i)
             xx8=tmp4*ax(i)
             cw1(i,j,k) = cmplx(xx1-xx4+xx6+xx7,-(-xx2-xx3+xx5-xx8), &
                  kind=mytype)        
          end do
       end do
    end do
#ifdef DEBUG
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             if (abs(cw1(i,j,k)) > 1.0e-4) then
                write(*,100) 'AFTER X',i,j,k,cw1(i,j,k)
             end if
          end do
       end do
    end do
#endif

    ! POST PROCESSING IN Y
    ! NEED to be in Y-pencil
    call transpose_x_to_y(cw1,cw2,sp)
    do k = sp%yst(3), sp%yen(3)
       do i = sp%yst(1), sp%yen(1)
          cw2b(i,1,k)=cw2(i,1,k)
          do j = 2,ny
             tmp1 = real(cw2(i,j,k), kind=mytype)
             tmp2 = aimag(cw2(i,j,k))
             tmp3 = real(cw2(i,ny-j+2,k), kind=mytype)
             tmp4 = aimag(cw2(i,ny-j+2,k))
             xx1=tmp1*by(j)
             xx2=tmp1*ay(j)
             xx3=tmp2*by(j)
             xx4=tmp2*ay(j)
             xx5=tmp3*by(j)
             xx6=tmp3*ay(j)
             xx7=tmp4*by(j)
             xx8=tmp4*ay(j)
             cw2b(i,j,k) = cmplx(xx1-xx4+xx6+xx7,-(-xx2-xx3+xx5-xx8), &
                  kind=mytype)        
          end do
       end do
    end do
#ifdef DEBUG
    do k = sp%yst(3), sp%yen(3)
       do j = sp%yst(2), sp%yen(2)
          do i = sp%yst(1), sp%yen(1)
             if (abs(cw2b(i,j,k)) > 1.0e-4) then
                write(*,100) 'AFTER Y',i,j,k,cw2b(i,j,k)
             end if
          end do
       end do
    end do
#endif
    ! back to X-pencil
    call transpose_y_to_x(cw2b,cw1,sp)
    
    ! POST PROCESSING IN Z
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*bz(k)-tmp2*az(k), &
                  tmp2*bz(k)+tmp1*az(k), kind=mytype)
#ifdef DEBUG
             if (abs(cw1(i,j,k)) > 1.0e-4) &
                  write(*,100) 'END',i,j,k,cw1(i,j,k)
#endif
          end do
       end do
    end do

    ! compute c2r transform, back to physical space
    call decomp_2d_fft_3d(cw1,rhs)

    if (bcz==1) then 
       do j=1,ph%zsz(2)
          do i=1,ph%zsz(1)
             do k=1,nz/2
                rw3(i,j,2*k-1)=rhs(i,j,k)
             end do
             do k=1,nz/2
                rw3(i,j,2*k)=rhs(i,j,nz-k+1)
             end do
          end do
       end do
       call transpose_z_to_y(rw3,rw2,ph)
    else if (bcz==0) then 
       call transpose_z_to_y(rhs,rw2,ph)   
    end if
    
    do k=ph%yst(3),ph%yen(3)
       do i=ph%yst(1),ph%yen(1)
          do j=1,ny/2
             rw2b(i,2*j-1,k)=rw2(i,j,k)
          end do
          do j=1,ny/2
             rw2b(i,2*j,k)=rw2(i,ny-j+1,k)
          end do
       enddo
    end do
    call transpose_y_to_x(rw2b,rw1,ph)
    do k=ph%xst(3),ph%xen(3)
       do j=ph%xst(2),ph%xen(2)
          do i=1,nx/2
             rw1b(2*i-1,j,k)=rw1(i,j,k)
          enddo
          do i=1,nx/2
             rw1b(2*i,j,k)=rw1(nx-i+1,j,k)
          enddo
       enddo
    end do
    call transpose_x_to_y(rw1b,rw2,ph)
    call transpose_y_to_z(rw2,rhs,ph)

  !  call decomp_2d_fft_finalize

   

    return
  end subroutine poisson_11x


  
  subroutine abxyz(ax,ay,az,bx,by,bz,nx,ny,nz,bcx,bcy,bcz)
	
    use param

    implicit none

    integer, intent(IN) :: nx,ny,nz
    integer, intent(IN) :: bcx,bcy,bcz
    real(mytype), dimension(:), intent(OUT) :: ax,bx
    real(mytype), dimension(:), intent(OUT) :: ay,by
    real(mytype), dimension(:), intent(OUT) :: az,bz

    integer :: i,j,k

    if (bcx==0) then
       do i=1,nx
          ax(i) = sin(real(i-1, kind=mytype)*PI/real(nx, kind=mytype))
          bx(i) = cos(real(i-1, kind=mytype)*PI/real(nx, kind=mytype))
       end do
    else if (bcx==1) then
       do i=1,nx
          ax(i) = sin(real(i-1, kind=mytype)*PI/2.0_mytype/ &
               real(nx, kind=mytype))
          bx(i) = cos(real(i-1, kind=mytype)*PI/2.0_mytype/ &
               real(nx, kind=mytype))
       end do
    end if

    if (bcy==0) then
       do j=1,ny
          ay(j) = sin(real(j-1, kind=mytype)*PI/real(ny, kind=mytype))
          by(j) = cos(real(j-1, kind=mytype)*PI/real(ny, kind=mytype))
       end do
    else if (bcy==1) then
       do j=1,ny
          ay(j) = sin(real(j-1, kind=mytype)*PI/2.0_mytype/ &
               real(ny, kind=mytype))
          by(j) = cos(real(j-1, kind=mytype)*PI/2.0_mytype/ &
               real(ny, kind=mytype))
       end do
    end if

    if (bcz==0) then
       do k=1,nz
          az(k) = sin(real(k-1, kind=mytype)*PI/real(nz, kind=mytype))
          bz(k) = cos(real(k-1, kind=mytype)*PI/real(nz, kind=mytype))
       end do
    else if (bcz==1) then
       do k=1,nz
          az(k) = sin(real(k-1, kind=mytype)*PI/2.0_mytype/ &
               real(nz, kind=mytype))
          bz(k) = cos(real(k-1, kind=mytype)*PI/2.0_mytype/ &
               real(nz, kind=mytype))
       end do
    end if

    return
  end subroutine abxyz

! ***********************************************************
!
subroutine waves ()
!
!***********************************************************

USE derivX 
USE derivY 
USE derivZ 
USE param
USE decomp_2d
USE variables
use decomp_2d_fft

implicit none

integer :: i,j,k
real(mytype) :: w,wp,w1,w1p 
complex(mytype) :: xyzk
complex(mytype) :: ytt,xtt,ztt,yt1,xt1,yt2,xt2
complex(mytype) :: xtt1,ytt1,ztt1,zt1,zt2,tmp1,tmp2,tmp3
complex(mytype) :: tmp4,tmp5,tmp6

xkx(:)=0. ; xk2(:)=0. ; yky(:)=0. ; yk2(:)=0.
zkz(:)=0. ; zk2(:)=0.

!WAVE NUMBER IN X
if (bcx==0) then
   do i=1,nx/2+1
      w=2.*pi*(i-1)/nx
      wp=acix6*2.*dx*sin(w/2.)+(bcix6*2.*dx)*sin(3./2.*w)
      wp=wp/(1.+2.*alcaix6*cos(w))
      xkx(i)=cmplx(nx*wp/xlx,nx*wp/xlx, kind=mytype)
      exs(i)=cmplx(nx*w/xlx,nx*w/xlx, kind=mytype)
      xk2(i)=cmplx((nx*wp/xlx)**2,(nx*wp/xlx)**2, kind=mytype)
   enddo
   do i=nx/2+2,nx
      xkx(i)=xkx(nx-i+2)
      exs(i)=exs(nx-i+2)
      xk2(i)=xk2(nx-i+2)
   enddo
else
   do i=1,nx
      w=2.*pi*0.5*(i-1)/nxm
      wp=acix6*2.*dx*sin(w/2.)+(bcix6*2.*dx)*sin(3./2.*w)
      wp=wp/(1.+2.*alcaix6*cos(w))
      xkx(i)=cmplx(nxm*wp/xlx,nxm*wp/xlx, kind=mytype)
      exs(i)=cmplx(nxm*w/xlx,nxm*w/xlx, kind=mytype)
      xk2(i)=cmplx((nxm*wp/xlx)**2,(nxm*wp/xlx)**2, kind=mytype)
   enddo
   xkx(1)=0.
   exs(1)=0.
   xk2(1)=0.
endif

!WAVE NUMBER IN Y
if (bcy==0) then
   do j=1,ny/2+1
      w=2.*pi*(j-1)/ny
      wp=aciy6*2.*dy*sin(w/2.)+(bciy6*2.*dy)*sin(3./2.*w)
      wp=wp/(1.+2.*alcaiy6*cos(w))
      if (istret==0) yky(j)=cmplx(ny*wp/yly,ny*wp/yly, kind=mytype)
      if (istret.ne.0) yky(j)=cmplx(ny*wp,ny*wp, kind=mytype)
      eys(j)=cmplx(ny*w/yly,ny*w/yly, kind=mytype)
      yk2(j)=cmplx((ny*wp/yly)**2,(ny*wp/yly)**2, kind=mytype)
   enddo
   do j=ny/2+2,ny
      yky(j)=yky(ny-j+2)
      eys(j)=eys(ny-j+2)
      yk2(j)=yk2(ny-j+2)
   enddo
else
   do j=1,ny
      w=2.*pi*0.5*(j-1)/nym
      wp=aciy6*2.*dy*sin(w/2.)+(bciy6*2.*dy)*sin(3./2.*w)
      wp=wp/(1.+2.*alcaiy6*cos(w))
      if (istret==0) yky(j)=cmplx(nym*wp/yly,nym*wp/yly, kind=mytype)
      if (istret.ne.0) yky(j)=cmplx(nym*wp,nym*wp, kind=mytype)
      eys(j)=cmplx(nym*w/yly,nym*w/yly, kind=mytype)
      yk2(j)=cmplx((nym*wp/yly)**2,(nym*wp/yly)**2, kind=mytype)
   enddo
   yky(1)=0.
   eys(1)=0.
   yk2(1)=0.
endif

!WAVE NUMBER IN Z
if (bcz==0) then
   do k=1,nz/2+1
      w=2.*pi*(k-1)/nz
      wp=aciz6*2.*dz*sin(w/2.)+(bciz6*2.*dz)*sin(3./2.*w)
      wp=wp/(1.+2.*alcaiz6*cos(w))
      zkz(k)=cmplx(nz*wp/zlz,nz*wp/zlz, kind=mytype)
      ezs(k)=cmplx(nz*w/zlz,nz*w/zlz, kind=mytype)
      zk2(k)=cmplx((nz*wp/zlz)**2,(nz*wp/zlz)**2, kind=mytype)
   enddo
else
   do k=1,nz/2+1
      w=2.*pi*0.5*(k-1)/nzm
      w1=2.*pi*0.5*(nzm-k+1)/nzm
      wp=aciz6*2.*dz*sin(w/2.)+(bciz6*2.*dz)*sin(3./2.*w)
      wp=wp/(1.+2.*alcaiz6*cos(w))
      w1p=aciz6*2.*dz*sin(w1/2.)+(bciz6*2.*dz)*sin(3./2.*w1)
      w1p=w1p/(1.+2.*alcaiz6*cos(w1))     
      zkz(k)=cmplx(nzm*wp/zlz,-nzm*w1p/zlz, kind=mytype)
      ezs(k)=cmplx(nzm*w/zlz,nzm*w1/zlz, kind=mytype)
      zk2(k)=cmplx((nzm*wp/zlz)**2,(nzm*w1p/zlz)**2, kind=mytype)
   enddo
endif
!
!if (nrank==0) then
!   do i=1,nx
!      print *,i,ezs(i)
!   enddo
!endif
!stop

if ((bcx==0).and.(bcz==0).and.bcy.ne.0) then
do k = sp%yst(3), sp%yen(3)
do j = sp%yst(2), sp%yen(2)
do i = sp%yst(1), sp%yen(1)
   xtt=cmplx((bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx/2.)+&
        cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx/2.)+&
        dicix6*2.*cos(real(exs(i), kind=mytype)*7.*dx/2.)),&
        (bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx/2.)+&
        cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx/2.)+&
        dicix6*2.*cos(real(exs(i), kind=mytype)*7.*dx/2.)), kind=mytype)
   ytt=cmplx((biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dy/2.)+&
        ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dy/2.)+&
        diciy6*2.*cos(real(eys(j), kind=mytype)*7.*dy/2.)),&
        (biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dy/2.)+&
        ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dy/2.)+&
        diciy6*2.*cos(real(eys(j), kind=mytype)*7.*dy/2.)), kind=mytype)
   ztt=cmplx((biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dz/2.)+&
        ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dz/2.)+&
        diciz6*2.*cos(real(ezs(k), kind=mytype)*7.*dz/2.)),&
        (biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dz/2.)+&
        ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dz/2.)+&
        diciz6*2.*cos(real(ezs(k), kind=mytype)*7.*dz/2.)), kind=mytype)
   xtt1=cmplx((aicix6*2.*cos(real(exs(i), kind=mytype)*dx/2.)),&
        (aicix6*2.*cos(real(exs(i), kind=mytype)*dx/2.)), kind=mytype)
   ytt1=cmplx((aiciy6*2.*cos(real(eys(j), kind=mytype)*dy/2.)),&
        (aiciy6*2.*cos(real(eys(j), kind=mytype)*dy/2.)), kind=mytype)
   ztt1=cmplx((aiciz6*2.*cos(real(ezs(k), kind=mytype)*dz/2.)),&
        (aiciz6*2.*cos(real(ezs(k), kind=mytype)*dz/2.)), kind=mytype)
   xt1=cmplx((1.+2.*ailcaix6*cos(real(exs(i), kind=mytype)*dx)),&
        (1.+2.*ailcaix6*cos(real(exs(i), kind=mytype)*dx)), kind=mytype)
   yt1=cmplx((1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dy)),&
        (1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dy)), kind=mytype)
   zt1=cmplx((1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dz)),&
        (1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dz)), kind=mytype)
   xt2=xk2(i)*((((ytt1+ytt)/yt1)*((ztt1+ztt)/zt1))**2)
   yt2=yk2(j)*((((xtt1+xtt)/xt1)*((ztt1+ztt)/zt1))**2)
   zt2=zk2(k)*((((xtt1+xtt)/xt1)*((ytt1+ytt)/yt1))**2)
   xyzk=xt2+yt2+zt2
   kxyz(i,j,k)=xyzk
!   print *,i,j,k, kxyz(i,j,k)
enddo
enddo
enddo
else
if (bcz==0) then
      do k = sp%xst(3),sp%xen(3)
      do j = sp%xst(2),sp%xen(2)
      do i = sp%xst(1),sp%xen(1)
         xtt=cmplx((bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx/2.)+&
              cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx/2.)+&
              dicix6*2.*cos(real(exs(i), kind=mytype)*7.*dx/2.)),&
              (bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx/2.)+&
              cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx/2.)+&
              dicix6*2.*cos(real(exs(i), kind=mytype)*7.*dx/2.)), kind=mytype)
         ytt=cmplx((biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dy/2.)+&
              ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dy/2.)+&
              diciy6*2.*cos(real(eys(j), kind=mytype)*7.*dy/2.)),&
              (biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dy/2.)+&
              ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dy/2.)+&
              diciy6*2.*cos(real(eys(j), kind=mytype)*7.*dy/2.)), kind=mytype)
         ztt=cmplx((biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dz/2.)+&
              ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dz/2.)+&
              diciz6*2.*cos(real(ezs(k), kind=mytype)*7.*dz/2.)),&
              (biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dz/2.)+&
              ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dz/2.)+&
              diciz6*2.*cos(real(ezs(k), kind=mytype)*7.*dz/2.)), kind=mytype)
         xtt1=cmplx((aicix6*2.*cos(real(exs(i), kind=mytype)*dx/2.)),&
              (aicix6*2.*cos(real(exs(i), kind=mytype)*dx/2.)), kind=mytype)
         ytt1=cmplx((aiciy6*2.*cos(real(eys(j), kind=mytype)*dy/2.)),&
              (aiciy6*2.*cos(real(eys(j), kind=mytype)*dy/2.)), kind=mytype)
         ztt1=cmplx((aiciz6*2.*cos(real(ezs(k), kind=mytype)*dz/2.)),&
              (aiciz6*2.*cos(real(ezs(k), kind=mytype)*dz/2.)), kind=mytype)
         xt1=cmplx((1.+2.*ailcaix6*cos(real(exs(i), kind=mytype)*dx)),&
              (1.+2.*ailcaix6*cos(real(exs(i))*dx)), kind=mytype)
         yt1=cmplx((1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dy)),&
              (1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dy)), kind=mytype)
         zt1=cmplx((1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dz)),&
              (1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dz)), kind=mytype)
         xt2=xk2(i)*((((ytt1+ytt)/yt1)*((ztt1+ztt)/zt1))**2)
         yt2=yk2(j)*((((xtt1+xtt)/xt1)*((ztt1+ztt)/zt1))**2)
         zt2=zk2(k)*((((xtt1+xtt)/xt1)*((ytt1+ytt)/yt1))**2)
         xyzk=xt2+yt2+zt2
         kxyz(i,j,k)=xyzk
!   print *,i,j,k, kxyz(i,j,k)
      enddo
      enddo
      enddo
   else
      do k = sp%xst(3),sp%xen(3)
      do j = sp%xst(2),sp%xen(2)
      do i = sp%xst(1),sp%xen(1)  
         xtt=cmplx((bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx/2.)+&
              cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx/2.)+&
              dicix6*2.*cos(real(exs(i), kind=mytype)*7.*dx/2.)),&
              (bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx/2.)+&
              cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx/2.)+&
              dicix6*2.*cos(real(exs(i), kind=mytype)*7.*dx/2.)), kind=mytype)
         ytt=cmplx((biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dy/2.)+&
              ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dy/2.)+&
              diciy6*2.*cos(real(eys(j), kind=mytype)*7.*dy/2.)),&
              (biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dy/2.)+&
              ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dy/2.)+&
              diciy6*2.*cos(real(eys(j), kind=mytype)*7.*dy/2.)), kind=mytype)
         !
         ztt=cmplx((biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dz/2.)+&
              ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dz/2.)+&
              diciz6*2.*cos(real(ezs(k), kind=mytype)*7.*dz/2.)),&
              (biciz6*2.*cos(aimag(ezs(k))*3.*dz/2.)+&
              ciciz6*2.*cos(aimag(ezs(k))*5.*dz/2.)+&
              diciz6*2.*cos(aimag(ezs(k))*7.*dz/2.)), kind=mytype)
         !
         xtt1=cmplx((aicix6*2.*cos(real(exs(i), kind=mytype)*dx/2.)),&
              (aicix6*2.*cos(real(exs(i), kind=mytype)*dx/2.)), kind=mytype)
         ytt1=cmplx((aiciy6*2.*cos(real(eys(j), kind=mytype)*dy/2.)),&
              (aiciy6*2.*cos(real(eys(j), kind=mytype)*dy/2.)), kind=mytype)
         !
         ztt1=cmplx((aiciz6*2.*cos(real(ezs(k), kind=mytype)*dz/2.)),&
              (aiciz6*2.*cos(aimag(ezs(k))*dz/2.)), kind=mytype)
         !
         xt1=cmplx((1.+2.*ailcaix6*cos(real(exs(i), kind=mytype)*dx)),&
              (1.+2.*ailcaix6*cos(real(exs(i), kind=mytype)*dx)), kind=mytype)
         yt1=cmplx((1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dy)),&
              (1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dy)), kind=mytype)
         zt1=cmplx((1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dz)),&
              (1.+2.*ailcaiz6*cos(aimag(ezs(k))*dz)), kind=mytype)
         
         tmp1=cmplx(real(ztt1+ztt, kind=mytype)/real(zt1, kind=mytype),&
              aimag(ztt1+ztt)/aimag(zt1), kind=mytype)
         tmp2=cmplx(real(ytt1+ytt, kind=mytype)/real(yt1, kind=mytype),&
              real(ytt1+ytt, kind=mytype)/real(yt1, kind=mytype), kind=mytype)
         tmp3=cmplx(real(xtt1+xtt, kind=mytype)/real(xt1, kind=mytype),&
              real(xtt1+xtt, kind=mytype)/real(xt1, kind=mytype), kind=mytype)
         
         tmp4=cmplx((real(tmp1, kind=mytype)*real(tmp2, kind=mytype))**2,(aimag(tmp1)*aimag(tmp2))**2, kind=mytype)
         tmp5=cmplx((real(tmp1, kind=mytype)*real(tmp3, kind=mytype))**2,(aimag(tmp1)*aimag(tmp3))**2, kind=mytype)
         tmp6=cmplx((real(tmp3, kind=mytype)*real(tmp2, kind=mytype))**2,(aimag(tmp3)*aimag(tmp2))**2, kind=mytype)
         
         tmp1=cmplx(real(tmp4, kind=mytype)*real(xk2(i), kind=mytype),aimag(tmp4)*aimag(xk2(i)), kind=mytype)
         tmp2=cmplx(real(tmp5, kind=mytype)*real(yk2(j), kind=mytype),aimag(tmp5)*aimag(yk2(j)), kind=mytype)
         tmp3=cmplx(real(tmp6, kind=mytype)*real(zk2(k), kind=mytype),aimag(tmp6)*aimag(zk2(k)), kind=mytype)

         xyzk=tmp1+tmp2+tmp3
         kxyz(i,j,k)=xyzk
!         print *,i,j,k,zt1,yt1
      enddo
      enddo
      enddo
   endif
endif


!          do k=1,1!nz
!          do j=1,ny
!          do i=1,1!!nx
!             print *,j,a(i,j,k,3),kxyz(i,j,k)
!          enddo
!          enddo
!          enddo

end subroutine waves

!**************************************************************************
!
subroutine matrice_refinement()
!
!**************************************************************************

USE decomp_2d
USE variables
USE param
USE var
USE MPI
USE derivX 
USE derivY 
USE derivZ 

implicit none

integer :: i,j,k

complex(mytype),dimension(sp%yst(1):sp%yen(1)) :: transx
complex(mytype),dimension(sp%yst(2):sp%yen(2)) :: transy
complex(mytype),dimension(sp%yst(3):sp%yen(3)) :: transz
real(mytype) :: xa0,xa1 
complex(mytype) :: ytt,xtt,ztt,yt1,xt1,yt2,xt2
complex(mytype) :: xtt1,ytt1,ztt1,zt1,zt2,tmp1,tmp2,tmp3

do i = sp%yst(1),sp%yen(1)
   xtt=cmplx((bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx/2.)+&
        cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx/2.)+&
        dicix6*2.*cos(real(exs(i), kind=mytype)*7.*dx/2.)),&
        (bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx/2.)+&
        cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx/2.)+&
        dicix6*2.*cos(real(exs(i), kind=mytype)*7.*dx/2.)), kind=mytype)
   xtt1=cmplx((aicix6*2.*cos(real(exs(i), kind=mytype)*dx/2.)),&
        (aicix6*2.*cos(real(exs(i), kind=mytype)*dx/2.)), kind=mytype)
   xt1=cmplx((1.+2.*ailcaix6*cos(real(exs(i), kind=mytype)*dx)),&
        (1.+2.*ailcaix6*cos(real(exs(i), kind=mytype)*dx)), kind=mytype)
   transx(i)=cmplx(real(xtt1+xtt, kind=mytype)/real(xt1, kind=mytype),&
        real(xtt1+xtt, kind=mytype)/real(xt1, kind=mytype), kind=mytype)!(xtt+xtt)/xt1
enddo
do j = sp%yst(2),sp%yen(2)
   ytt=cmplx((biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dy/2.)+&
        ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dy/2.)+&
        diciy6*2.*cos(real(eys(j), kind=mytype)*7.*dy/2.)),&
        (biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dy/2.)+&
        ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dy/2.)+&
        diciy6*2.*cos(real(eys(j), kind=mytype)*7.*dy/2.)), kind=mytype)
   ytt1=cmplx((aiciy6*2.*cos(real(eys(j), kind=mytype)*dy/2.)),&
        (aiciy6*2.*cos(real(eys(j), kind=mytype)*dy/2.)), kind=mytype)
   yt1=cmplx((1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dy)),&
        (1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dy)), kind=mytype)
   transy(j)=cmplx(real(ytt1+ytt, kind=mytype)/real(yt1, kind=mytype),&
        real(ytt1+ytt, kind=mytype)/real(yt1, kind=mytype), kind=mytype)!(ytt+ytt)/yt1
enddo
if (bcz==0) then
   do k = sp%yst(3),sp%yen(3)
      ztt=cmplx((biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dz/2.)+&
           ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dz/2.)+&
           diciz6*2.*cos(real(ezs(k), kind=mytype)*7.*dz/2.)),&
           (biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dz/2.)+&
           ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dz/2.)+&
           diciz6*2.*cos(real(ezs(k), kind=mytype)*7.*dz/2.)), kind=mytype)
      ztt1=cmplx((aiciz6*2.*cos(real(ezs(k), kind=mytype)*dz/2.)),&
           (aiciz6*2.*cos(real(ezs(k), kind=mytype)*dz/2.)), kind=mytype)
      zt1=cmplx((1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dz)),&
           (1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dz)), kind=mytype)
      transz(k)=cmplx(real(ztt1+ztt, kind=mytype)/real(zt1, kind=mytype),&
           aimag(ztt1+ztt)/aimag(zt1), kind=mytype)!(ztt+ztt)/zt1
   enddo
else
   do k = sp%yst(3),sp%yen(3)
      ztt=cmplx((biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dz/2.)+&
           ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dz/2.)),&
           (biciz6*2.*cos(aimag(ezs(k))*3.*dz/2.)+&
           ciciz6*2.*cos(aimag(ezs(k))*5.*dz/2.)), kind=mytype)
      ztt1=cmplx((aiciz6*2.*cos(real(ezs(k), kind=mytype)*dz/2.)),&
           (aiciz6*2.*cos(aimag(ezs(k))*dz/2.)), kind=mytype)
      zt1=cmplx((1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dz)),&
           (1.+2.*ailcaiz6*cos(aimag(ezs(k))*dz)), kind=mytype)
      transz(k)=cmplx(real(ztt1+ztt, kind=mytype)/real(zt1, kind=mytype),&
           aimag(ztt1+ztt)/aimag(zt1), kind=mytype)!(ztt+ztt)/zt1
   enddo
endif

if ((istret==1).or.(istret==2)) then

  
   xa0=alpha/pi+1./2./beta/pi
   if (istret==1) xa1=1./4./beta/pi
   if (istret==2) xa1=-1./4./beta/pi
!
!construction of the pentadiagonal matrice
!
   do k=sp%yst(3),sp%yen(3)
   do j=1,ny/2
   do i=sp%yst(1),sp%yen(1)
      cw22(i,j,k)=cmplx(real(yky(2*j-1), kind=mytype)*real(transx(i), kind=mytype)*real(transz(k), kind=mytype),&
           aimag(yky(2*j-1))*aimag(transx(i))*aimag(transz(k)), kind=mytype)
      cw2(i,j,k)=cmplx(real(yky(2*j), kind=mytype)*real(transx(i), kind=mytype)*real(transz(k), kind=mytype),&
           aimag(yky(2*j))*aimag(transx(i))*aimag(transz(k)), kind=mytype)
   enddo
   enddo
   enddo


   

!main diagonal 
   do k=sp%yst(3),sp%yen(3)
   do j=2,ny/2-1
   do i=sp%yst(1),sp%yen(1)
      a(i,j,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(2*j-1), kind=mytype)*real(transy(2*j-1), kind=mytype)&
           *real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
           -(real(zk2(k), kind=mytype)*real(transy(2*j-1), kind=mytype)*real(transy(2*j-1), kind=mytype)*&
           real(transx(i), kind=mytype)*real(transx(i), kind=mytype))&
           -real(cw22(i,j,k), kind=mytype)*real(cw22(i,j,k), kind=mytype)*xa0*xa0-&
           xa1*xa1*(real(cw22(i,j,k), kind=mytype)*real(cw22(i,j-1,k), kind=mytype)+real(cw22(i,j,k), kind=mytype)*&
           real(cw22(i,j+1,k), kind=mytype)),&
           -(aimag(xk2(i))*aimag(transy(2*j-1))*aimag(transy(2*j-1))*aimag(transz(k))*aimag(transz(k)))&
           -(aimag(zk2(k))*aimag(transy(2*j-1))*aimag(transy(2*j-1))*aimag(transx(i))*aimag(transx(i)))&
           -aimag(cw22(i,j,k))*aimag(cw22(i,j,k))*xa0*xa0-&
           xa1*xa1*(aimag(cw22(i,j,k))*aimag(cw22(i,j-1,k))+aimag(cw22(i,j,k))*aimag(cw22(i,j+1,k))), kind=mytype)
      a2(i,j,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(2*j), kind=mytype)*real(transy(2*j), kind=mytype)*&
           real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
           -(real(zk2(k), kind=mytype)*real(transy(2*j), kind=mytype)*real(transy(2*j), kind=mytype)*&
           real(transx(i), kind=mytype)*real(transx(i), kind=mytype))&
           -real(cw2(i,j,k), kind=mytype)*real(cw2(i,j,k), kind=mytype)*xa0*xa0-&
           xa1*xa1*(real(cw2(i,j,k), kind=mytype)*real(cw2(i,j-1,k), kind=mytype)+real(cw2(i,j,k), kind=mytype)*&
           real(cw2(i,j+1,k), kind=mytype)),&
           -(aimag(xk2(i))*aimag(transy(2*j))*aimag(transy(2*j))*aimag(transz(k))*aimag(transz(k)))&
           -(aimag(zk2(k))*aimag(transy(2*j))*aimag(transy(2*j))*aimag(transx(i))*aimag(transx(i)))&
           -aimag(cw2(i,j,k))*aimag(cw2(i,j,k))*xa0*xa0-&
           xa1*xa1*(aimag(cw2(i,j,k))*aimag(cw2(i,j-1,k))+aimag(cw2(i,j,k))*aimag(cw2(i,j+1,k))), kind=mytype)
   enddo
   enddo
   do i=sp%yst(1),sp%yen(1)
      a(i,1,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(1), kind=mytype)*real(transy(1), kind=mytype)*&
           real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
           -(real(zk2(k), kind=mytype)*real(transy(1), kind=mytype)*real(transy(1), kind=mytype)*&
           real(transx(i), kind=mytype)*real(transx(i), kind=mytype))&
           -real(cw22(i,1,k), kind=mytype)*real(cw22(i,1,k), kind=mytype)*xa0*xa0-xa1*xa1*(real(cw22(i,1,k), kind=mytype)*&
           real(cw22(i,2,k), kind=mytype)),&
           -(aimag(xk2(i))*aimag(transy(1))*aimag(transy(1))*aimag(transz(k))*aimag(transz(k)))&
           -(aimag(zk2(k))*aimag(transy(1))*aimag(transy(1))*aimag(transx(i))*aimag(transx(i)))&
           -aimag(cw22(i,1,k))*aimag(cw22(i,1,k))*xa0*xa0-xa1*xa1*(aimag(cw22(i,1,k))*aimag(cw22(i,2,k))), kind=mytype)
      a(i,ny/2,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(ny-2), kind=mytype)*real(transy(ny-2), kind=mytype)&
           *real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
           -(real(zk2(k), kind=mytype)*real(transy(ny-2), kind=mytype)*real(transy(ny-2), kind=mytype)*&
           real(transx(i), kind=mytype)*real(transx(i), kind=mytype))&
           -real(cw22(i,ny/2,k), kind=mytype)*real(cw22(i,ny/2,k), kind=mytype)*xa0*xa0-&
           xa1*xa1*(real(cw22(i,ny/2,k), kind=mytype)*real(cw22(i,ny/2-1,k), kind=mytype)),&
           -(aimag(xk2(i))*aimag(transy(ny-2))*aimag(transy(ny-2))*aimag(transz(k))*aimag(transz(k)))&
           -(aimag(zk2(k))*aimag(transy(ny-2))*aimag(transy(ny-2))*aimag(transx(i))*aimag(transx(i)))&
           -aimag(cw22(i,ny/2,k))*aimag(cw22(i,ny/2,k))*xa0*xa0-&
           xa1*xa1*(aimag(cw22(i,ny/2,k))*aimag(cw22(i,ny/2-1,k))), kind=mytype)
      a2(i,1,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(2), kind=mytype)*real(transy(2), kind=mytype)*&
           real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
           -(real(zk2(k), kind=mytype)*real(transy(2), kind=mytype)*real(transy(2), kind=mytype)*&
           real(transx(i), kind=mytype)*real(transx(i), kind=mytype))&
           -real(cw2(i,1,k), kind=mytype)*real(cw2(i,1,k), kind=mytype)*(xa0-xa1)*(xa0+xa1)-xa1*xa1*(real(cw2(i,1,k), kind=mytype)*&
           real(cw2(i,2,k), kind=mytype)),&
           -(aimag(xk2(i))*aimag(transy(2))*aimag(transy(2))*aimag(transz(k))*aimag(transz(k)))&
           -(aimag(zk2(k))*aimag(transy(2))*aimag(transy(2))*aimag(transx(i))*aimag(transx(i)))&
           -aimag(cw2(i,1,k))*aimag(cw2(i,1,k))*(xa0-xa1)*(xa0+xa1)-xa1*xa1*(aimag(cw2(i,1,k))*aimag(cw2(i,2,k))), kind=mytype)
      a2(i,ny/2,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(ny-1), kind=mytype)*real(transy(ny-1), kind=mytype)*&
           real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
           -(real(zk2(k), kind=mytype)*real(transy(ny-1), kind=mytype)*real(transy(ny-1), kind=mytype)*&
           real(transx(i), kind=mytype)*real(transx(i), kind=mytype))&
           -real(cw2(i,ny/2,k), kind=mytype)*real(cw2(i,ny/2,k), kind=mytype)*(xa0+xa1)*(xa0+xa1)-&
           xa1*xa1*(real(cw2(i,ny/2,k), kind=mytype)*real(cw2(i,ny/2-1,k), kind=mytype)),&
           -(aimag(xk2(i))*aimag(transy(ny-1))*aimag(transy(ny-1))*aimag(transz(k))*aimag(transz(k)))&
           -(aimag(zk2(k))*aimag(transy(ny-1))*aimag(transy(ny-1))*aimag(transx(i))*aimag(transx(i)))&
           -aimag(cw2(i,ny/2,k))*aimag(cw2(i,ny/2,k))*(xa0+xa1)*(xa0+xa1)-&
           xa1*xa1*(aimag(cw2(i,ny/2,k))*aimag(cw2(i,ny/2-1,k))), kind=mytype)
   enddo
   enddo
   
 

 

!sup diag +1
   do k=sp%yst(3),sp%yen(3)
   do j=2,ny/2-1
   do i=sp%yst(1),sp%yen(1)   
      a(i,j,k,4)=cmplx(xa0*xa1*(real(cw22(i,j,k), kind=mytype)*real(cw22(i,j+1,k), kind=mytype)+real(cw22(i,j+1,k), kind=mytype)*&
           real(cw22(i,j+1,k), kind=mytype)),&
           xa0*xa1*(aimag(cw22(i,j,k))*aimag(cw22(i,j+1,k))+aimag(cw22(i,j+1,k))*aimag(cw22(i,j+1,k))), kind=mytype)
      a2(i,j,k,4)=cmplx(xa0*xa1*(real(cw2(i,j,k), kind=mytype)*real(cw2(i,j+1,k), kind=mytype)+real(cw2(i,j+1,k), kind=mytype)*&
           real(cw2(i,j+1,k), kind=mytype)),&
           xa0*xa1*(aimag(cw2(i,j,k))*aimag(cw2(i,j+1,k))+aimag(cw2(i,j+1,k))*aimag(cw2(i,j+1,k))), kind=mytype)
   enddo
   enddo
   do i=sp%yst(1),sp%yen(1)   
      a(i,1,k,4)=2.*cmplx((xa0*xa1*(real(cw22(i,1,k), kind=mytype)*real(cw22(i,2,k), kind=mytype)+real(cw22(i,2,k), kind=mytype)*&
           real(cw22(i,2,k), kind=mytype))),&
           (xa0*xa1*(aimag(cw22(i,1,k))*aimag(cw22(i,2,k))+aimag(cw22(i,2,k))*aimag(cw22(i,2,k)))), kind=mytype)
      a2(i,1,k,4)=cmplx((xa0-xa1)*xa1*(real(cw2(i,1,k), kind=mytype)*real(cw2(i,2,k), kind=mytype))&
           +xa0*xa1*(real(cw2(i,2,k), kind=mytype)*real(cw2(i,2,k), kind=mytype)),&
           (xa0-xa1)*xa1*(aimag(cw2(i,1,k))*aimag(cw2(i,2,k)))&
           +xa0*xa1*(aimag(cw2(i,2,k))*aimag(cw2(i,2,k))), kind=mytype)
      a2(i,ny/2-1,k,4)=cmplx(xa0*xa1*(real(cw2(i,ny/2-1,k), kind=mytype)*real(cw2(i,ny/2,k), kind=mytype))&
           +(xa0+xa1)*xa1*(real(cw2(i,ny/2,k), kind=mytype)*real(cw2(i,ny/2,k), kind=mytype)),&
           xa0*xa1*(aimag(cw2(i,ny/2-1,k))*aimag(cw2(i,ny/2,k)))&
           +(xa0+xa1)*xa1*(aimag(cw2(i,ny/2,k))*aimag(cw2(i,ny/2,k))), kind=mytype)
      a2(i,ny/2,k,4)=0.
   enddo
   enddo



!sup diag +2
   do k=sp%yst(3),sp%yen(3)
   do i=sp%yst(1),sp%yen(1)   
   do j=1,ny/2-2
      a(i,j,k,5)=cmplx(-real(cw22(i,j+1,k), kind=mytype)*real(cw22(i,j+2,k), kind=mytype)*xa1*xa1,&
           -aimag(cw22(i,j+1,k))*aimag(cw22(i,j+2,k))*xa1*xa1, kind=mytype)
      a2(i,j,k,5)=cmplx(-real(cw2(i,j+1,k), kind=mytype)*real(cw2(i,j+2,k), kind=mytype)*xa1*xa1,&
           -aimag(cw2(i,j+1,k))*aimag(cw2(i,j+2,k))*xa1*xa1, kind=mytype)
   enddo
   a(i,1,k,5)=cmplx(real(a(i,1,k,5), kind=mytype)*2.,aimag(a(i,1,k,5))*2., kind=mytype)
   a(i,ny/2-1,k,5)=0.
   a(i,ny/2,k,5)=0.
   a2(i,ny/2-1,k,5)=0.
   a2(i,ny/2,k,5)=0. 
   enddo
   enddo



!inf diag -1
   do k=sp%yst(3),sp%yen(3)
   do i=sp%yst(1),sp%yen(1)   
   do j=2,ny/2
      a(i,j,k,2)=cmplx(xa0*xa1*(real(cw22(i,j,k), kind=mytype)*real(cw22(i,j-1,k), kind=mytype)+real(cw22(i,j-1,k), kind=mytype)*&
           real(cw22(i,j-1,k), kind=mytype)),&
           xa0*xa1*(aimag(cw22(i,j,k))*aimag(cw22(i,j-1,k))+aimag(cw22(i,j-1,k))*aimag(cw22(i,j-1,k))), kind=mytype)
      a2(i,j,k,2)=cmplx(xa0*xa1*(real(cw2(i,j,k), kind=mytype)*real(cw2(i,j-1,k), kind=mytype)+real(cw2(i,j-1,k), kind=mytype)*&
           real(cw2(i,j-1,k), kind=mytype)),&
           xa0*xa1*(aimag(cw2(i,j,k))*aimag(cw2(i,j-1,k))+aimag(cw2(i,j-1,k))*aimag(cw2(i,j-1,k))), kind=mytype)
   enddo
   a(i,1,k,2)=0.
   a2(i,1,k,2)=0.
   a2(i,2,k,2)=cmplx(xa0*xa1*(real(cw2(i,2,k), kind=mytype)*real(cw2(i,1,k), kind=mytype))&
        +(xa0+xa1)*xa1*(real(cw2(i,1,k), kind=mytype)*real(cw2(i,1,k), kind=mytype)),&
        xa0*xa1*(aimag(cw2(i,2,k))*aimag(cw2(i,1,k)))&
        +(xa0+xa1)*xa1*(aimag(cw2(i,1,k))*aimag(cw2(i,1,k))), kind=mytype)
   a2(i,ny/2,k,2)=cmplx((xa0+xa1)*xa1*(real(cw2(i,ny/2,k), kind=mytype)*real(cw2(i,ny/2-1,k), kind=mytype))&
        +xa0*xa1*(real(cw2(i,ny/2-1,k), kind=mytype)*real(cw2(i,ny/2-1,k), kind=mytype)),&
        (xa0+xa1)*xa1*(aimag(cw2(i,ny/2,k))*aimag(cw2(i,ny/2-1,k)))&
        +xa0*xa1*(aimag(cw2(i,ny/2-1,k))*aimag(cw2(i,ny/2-1,k))), kind=mytype)
   enddo
   enddo
!inf diag -2
   do k=sp%yst(3),sp%yen(3)
   do i=sp%yst(1),sp%yen(1)  
   do j=3,ny/2
      a(i,j,k,1)=cmplx(-real(cw22(i,j-1,k), kind=mytype)*real(cw22(i,j-2,k), kind=mytype)*xa1*xa1,&
           -aimag(cw22(i,j-1,k))*aimag(cw22(i,j-2,k))*xa1*xa1, kind=mytype)
      a2(i,j,k,1)=cmplx(-real(cw2(i,j-1,k), kind=mytype)*real(cw2(i,j-2,k), kind=mytype)*xa1*xa1,&
           -aimag(cw2(i,j-1,k))*aimag(cw2(i,j-2,k))*xa1*xa1, kind=mytype)
   enddo
   a(i,1,k,1)=0.
   a(i,2,k,1)=0.
   a2(i,1,k,1)=0.
   a2(i,2,k,1)=0.
   enddo
   enddo
!not to have a singular matrice
   do k=sp%yst(3),sp%yen(3)
   do i=sp%yst(1),sp%yen(1)
      if ((real(xk2(i), kind=mytype)==0.).and.(real(zk2(k), kind=mytype)==0)) then
         a(i,1,k,3)=cmplx(1.,1., kind=mytype)
         a(i,1,k,4)=0.
         a(i,1,k,5)=0.
      endif
   enddo
   enddo

else
   xa0=alpha/pi+1./2./beta/pi
   xa1=-1./4./beta/pi 
!
!construction of the pentadiagonal matrice
!   
   do k=sp%yst(3),sp%yen(3)
   do j=1,nym
   do i=sp%yst(1),sp%yen(1)
      cw22(i,j,k)=cmplx(real(yky(j), kind=mytype)*real(transx(i), kind=mytype)*real(transz(k), kind=mytype),&
           aimag(yky(j))*aimag(transx(i))*aimag(transz(k)), kind=mytype)
   enddo
   enddo
   enddo

!main diagonal 
   do k=sp%yst(3),sp%yen(3)
   do j=2,nym-1
   do i=sp%yst(1),sp%yen(1)
      a3(i,j,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(j), kind=mytype)*real(transy(j), kind=mytype)*&
           real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
           -(real(zk2(k), kind=mytype)*real(transy(j), kind=mytype)*real(transy(j), kind=mytype)*&
           real(transx(i), kind=mytype)*real(transx(i), kind=mytype))&
           -real(cw22(i,j,k), kind=mytype)*real(cw22(i,j,k), kind=mytype)*xa0*xa0-&
           xa1*xa1*(real(cw22(i,j,k), kind=mytype)*real(cw22(i,j-1,k), kind=mytype)+real(cw22(i,j,k), kind=mytype)*&
           real(cw22(i,j+1,k), kind=mytype)),&
           -(aimag(xk2(i))*aimag(transy(j))*aimag(transy(j))*aimag(transz(k))*aimag(transz(k)))&
           -(aimag(zk2(k))*aimag(transy(j))*aimag(transy(j))*aimag(transx(i))*aimag(transx(i)))&
           -aimag(cw22(i,j,k))*aimag(cw22(i,j,k))*xa0*xa0-&
           xa1*xa1*(aimag(cw22(i,j,k))*aimag(cw22(i,j-1,k))+aimag(cw22(i,j,k))*aimag(cw22(i,j+1,k))), kind=mytype)
   enddo
   enddo
   enddo

   do k=sp%yst(3),sp%yen(3)
   do i=sp%yst(1),sp%yen(1)
      a3(i,1,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(1), kind=mytype)*real(transy(1), kind=mytype)*&
           real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
           -(real(zk2(k), kind=mytype)*real(transy(1), kind=mytype)*real(transy(1), kind=mytype)*real(transx(i), kind=mytype)*&
           real(transx(i), kind=mytype))&
           -real(cw22(i,1,k), kind=mytype)*real(cw22(i,1,k), kind=mytype)*xa0*xa0-xa1*xa1*(real(cw22(i,1,k), kind=mytype)*&
           real(cw22(i,2,k), kind=mytype)),&
           -(aimag(xk2(i))*aimag(transy(1))*aimag(transy(1))*aimag(transz(k))*aimag(transz(k)))&
           -(aimag(zk2(k))*aimag(transy(1))*aimag(transy(1))*aimag(transx(i))*aimag(transx(i)))&
           -aimag(cw22(i,1,k))*aimag(cw22(i,1,k))*xa0*xa0-xa1*xa1*(aimag(cw22(i,1,k))*aimag(cw22(i,2,k))), kind=mytype)
      a3(i,nym,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(nym), kind=mytype)*real(transy(nym), kind=mytype)*&
           real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
           -(real(zk2(k), kind=mytype)*real(transy(nym), kind=mytype)*real(transy(nym), kind=mytype)*real(transx(i), kind=mytype)*&
           real(transx(i), kind=mytype))&
           -real(cw22(i,nym,k), kind=mytype)*real(cw22(i,nym,k), kind=mytype)*xa0*xa0-&
           xa1*xa1*(real(cw22(i,nym,k), kind=mytype)*real(cw22(i,nym-1,k), kind=mytype)),&
           -(aimag(xk2(i))*aimag(transy(nym))*aimag(transy(nym))*aimag(transz(k))*aimag(transz(k)))&
           -(aimag(zk2(k))*aimag(transy(nym))*aimag(transy(nym))*aimag(transx(i))*aimag(transx(i)))&
           -aimag(cw22(i,nym,k))*aimag(cw22(i,nym,k))*xa0*xa0-&
           xa1*xa1*(aimag(cw22(i,nym,k))*aimag(cw22(i,nym-1,k))), kind=mytype)
   enddo
   enddo

   


!sup diag +1
   do k=sp%yst(3),sp%yen(3)
   do i=sp%yst(1),sp%yen(1)
   do j=2,nym-1
      a3(i,j,k,4)=cmplx(xa0*xa1*(real(cw22(i,j,k), kind=mytype)*real(cw22(i,j+1,k), kind=mytype)+real(cw22(i,j+1,k), kind=mytype)*&
           real(cw22(i,j+1,k), kind=mytype)),&
           xa0*xa1*(aimag(cw22(i,j,k))*aimag(cw22(i,j+1,k))+aimag(cw22(i,j+1,k))*aimag(cw22(i,j+1,k))), kind=mytype)
   enddo
   a3(i,1,k,4)=cmplx((xa0*xa1*(real(cw22(i,1,k), kind=mytype)*real(cw22(i,2,k), kind=mytype)+real(cw22(i,2,k), kind=mytype)*&
        real(cw22(i,2,k), kind=mytype))),&
        (xa0*xa1*(aimag(cw22(i,1,k))*aimag(cw22(i,2,k))+aimag(cw22(i,2,k))*aimag(cw22(i,2,k)))), kind=mytype)
   enddo
   enddo
!sup diag +2
   do k=sp%yst(3),sp%yen(3)
   do i=sp%yst(1),sp%yen(1)
   do j=1,nym-2
      a3(i,j,k,5)=cmplx(-real(cw22(i,j+1,k), kind=mytype)*real(cw22(i,j+2,k), kind=mytype)*xa1*xa1,&
           -aimag(cw22(i,j+1,k))*aimag(cw22(i,j+2,k))*xa1*xa1, kind=mytype)
   enddo
   !a3(i,1,k,5)=a3(i,1,k,5)*2.
   !a3(i,1,k,5)=0.
   a3(i,nym-1,k,5)=0.
   a3(i,nym,k,5)=0.
   enddo
   enddo


!inf diag -1
   do k=sp%yst(3),sp%yen(3)
   do i=sp%yst(1),sp%yen(1)
   do j=2,nym
      a3(i,j,k,2)=cmplx(xa0*xa1*(real(cw22(i,j,k), kind=mytype)*real(cw22(i,j-1,k), kind=mytype)+real(cw22(i,j-1,k), kind=mytype)*&
           real(cw22(i,j-1,k), kind=mytype)),&
           xa0*xa1*(aimag(cw22(i,j,k))*aimag(cw22(i,j-1,k))+aimag(cw22(i,j-1,k))*aimag(cw22(i,j-1,k))), kind=mytype)
   enddo
   a3(i,1,k,2)=0.
   enddo
   enddo
!inf diag -2
   do k=sp%yst(3),sp%yen(3)
   do i=sp%yst(1),sp%yen(1)
   do j=3,nym
      a3(i,j,k,1)=cmplx(-real(cw22(i,j-1,k), kind=mytype)*real(cw22(i,j-2,k), kind=mytype)*xa1*xa1,&
           -aimag(cw22(i,j-1,k))*aimag(cw22(i,j-2,k))*xa1*xa1, kind=mytype)
   enddo
   a3(i,1,k,1)=0.
   a3(i,2,k,1)=0.
   enddo
   enddo 

!not to have a singular matrice
!   do k=sp%yst(3),sp%yen(3)
!   do i=sp%yst(1),sp%yen(1)
!      if ((xkx(i)==0.).and.(zkz(k)==0)) then
if (nrank==0) then
   a3(1,1,1,3)=cmplx(1.,1., kind=mytype)
   a3(1,1,1,4)=0.
   a3(1,1,1,5)=0.
endif
!      endif
!   enddo
!   enddo
endif


   


return
end subroutine matrice_refinement

end module decomp_2d_poisson

