module flow_type
  use decomp_2d, only : mytype
  integer :: initstats1,initstats2

end module flow_type

subroutine ft_parameter(arg)

  USE param
  USE variables
  USE flow_type
  USE complex_geometry
  USE decomp_2d, only : nrank
  implicit none

  logical,intent(in) :: arg
  character :: a

  nclx1 = 0 !Boundary condition in x=0  (0: Periodic, 1:Free-slip, 2: Dirichlet)
  nclxn = 0 !Boundary condition in x=Lx (0: Periodic, 1:Free-slip, 2: Dirichlet)
  ncly1 = 2 !Boundary condition in y=0  (0: Periodic, 1:Free-slip, 2: Dirichlet)
  nclyn = 2 !Boundary condition in y=Ly (0: Periodic, 1:Free-slip, 2: Dirichlet)
  nclz1 = 0 !Boundary condition in z=0  (0: Periodic, 1:Free-slip, 2: Dirichlet)
  nclzn = 0 !Boundary condition in z=Lz (0: Periodic, 1:Free-slip, 2: Dirichlet)

  open(10,file='BC-Periodic-hill.prm',status='unknown',form='formatted')
  read (10,*) a ! 
  read (10,*) a ! INCOMPACT 3D computational parameters
  read (10,*) a !
  read (10,*) nx
  read (10,*) ny
  read (10,*) nz
  read (10,*) nphi
  read (10,*) p_row
  read (10,*) p_col
  if (arg) then
     close(10)
     return
  endif
  read (10,*) a ! 
  read (10,*) a ! INCOMPACT 3D Flow parameters
  read (10,*) a !
  read (10,*) xlx
  read (10,*) yly
  read (10,*) zlz
  read (10,*) re
  read (10,*) noise
  read (10,*) dt
  read (10,*) a !
  read (10,*) a ! INCOMPACT3D Flow configuration
  read (10,*) a !
  read (10,*) iin
  read (10,*) ifirst
  read (10,*) ilast
  read (10,*) nscheme
  read (10,*) istret
  read (10,*) beta
  read (10,*) a !
  read (10,*) a ! INCOMPACT 3D File parameters
  read (10,*) a !
  read (10,*) ilit
  read (10,*) icheckpoint
  read (10,*) imodulo
  read (10,*) initstats1
  read (10,*) initstats2
  read (10,*) a !
  read (10,*) a ! NUMERICAL DISSIPATION
  read (10,*) a !
  read (10,*) jLES
  read (10,*) fpi2
  read (10,*) a !
  read (10,*) a ! INCOMPACT 3D Body (old school)
  read (10,*) a !
  read (10,*) ivirt
  read (10,*) a !
  read (10,*) a ! INCOMPACT 3D Forcing with Lagrangian Polynomials
  read (10,*) a !
  read (10,*) ilag
  read (10,*) npif 
  read (10,*) izap  
  read (10,*) nraf
  read (10,*) nobjmax
  close(10)

  if (nrank==0) then
     print *,'=======================Periodic hill======================='
     write(*,"(' initstats1         : ',I15)") initstats1
     write(*,"(' initstats2         : ',I15)") initstats2
     print *,'==========================================================='
  endif
  return
end subroutine ft_parameter
!********************************************************************
subroutine geomcomplex(epsi,nxi,nxf,ny,nyi,nyf,nzi,nzf,dx,yp,dz,remp)
  use param, only : xlx,zero,one,two,three,nine,fourteen,twenty,twentyeight
  use decomp_2d, only : mytype
  implicit none
  !
  real(mytype),dimension(nxi:nxf,nyi:nyf,nzi:nzf) :: epsi
  real(mytype),dimension(ny) :: yp
  integer                    :: nxi,nxf,ny,nyi,nyf,nzi,nzf
  real(mytype)               :: dx,dz
  real(mytype)               :: remp
  integer                    :: i,j,k
  real(mytype)               :: xm,ym
  real(mytype)               :: zeromach
  !
  real(mytype), dimension(nxi:nxf) :: dune
  real(mytype) :: y_bump
  !
  zeromach=one
  do while ((one + zeromach / two) .gt. one)
     zeromach = zeromach/two
  end do
  zeromach = 1.0e1*zeromach
  !
  y_bump=zero
  dune=zero
  do i=nxi,nxf
     xm=real(i-1,mytype)*dx
     if (xm.gt.xlx/two) then
        xm = (xlx-xm)*twentyeight
     else
        xm = xm*twentyeight
     endif
     if ((xm.ge.zero).and.(xm.le.nine)) then
        y_bump=min(28.,28+0.006775070969851*xm**two-2.124527775800E-03*xm**three)
     endif
     if ((xm.ge.9.).and.(xm.le.fourteen)) then
        y_bump=  2.507355893131E+01      +9.754803562315E-01*xm&
             -1.016116352781E-01*xm**two +1.889794677828E-03*xm**three
     endif
     if ((xm.ge.14.).and.(xm.le.twenty)) then
        y_bump=  2.579601052357E+01      +8.206693007457E-01*xm &
             -9.055370274339E-02*xm**two +1.626510569859E-03*xm**three
     endif
     if ((xm.ge.20.).and.(xm.le.30._mytype)) then
        y_bump=  4.046435022819E+01      -1.379581654948E+00*xm &
             +1.945884504128E-02*xm**two -2.070318932190E-04*xm**three
     endif
     if ((xm.ge.30.).and.(xm.le.40._mytype)) then
        y_bump=  1.792461334664E+01      +8.743920332081E-01*xm &
             -5.567361123058E-02*xm**two +6.277731764683E-04*xm**three
     endif
     if ((xm.ge.40.).and.(xm.le.54._mytype)) then
        y_bump=max(0.,5.639011190988E+01  -2.010520359035E+00*xm &
             +1.644919857549E-02*xm**two  +2.674976141766E-05*xm**three)
     endif
     dune(i)=y_bump/twentyeight
  enddo

  do k=nzi,nzf
     do j=nyi,nyf
        ym=yp(j)
        do i=nxi,nxf
           if (ym-dune(i).le.zeromach) then
              epsi(i,j,k)=remp
           endif
        enddo
     enddo
  enddo
  return
end subroutine geomcomplex
!********************************************************************
subroutine init (ux1,uy1,uz1,ep1,phi1,dux1,duy1,duz1,phis1,phiss1)

  USE decomp_2d
  USE decomp_2d_io
  USE variables
  USE param
  USE MPI

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1,phis1,phiss1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1

  real(mytype) :: y,r,um,r3,x,z,h,ct
  real(mytype) :: cx0,cy0,cz0,hg,lg
  integer :: k,j,i,ijk,fh,ierror,ii,is,code
  integer (kind=MPI_OFFSET_KIND) :: disp

  if (iscalar==1) then

     phi1 = one !change as much as you want

     !do not delete this
     phis1=phi1
     phiss1=phis1

  endif
  ux1=zero;uy1=zero;uz1=zero
  if (iin.ne.0) then
     call system_clock(count=code)
     if (iin.eq.2) code=0
     call random_seed(size = ii)
     call random_seed(put = code+63946*nrank*(/ (i - 1, i = 1, ii) /))

     call random_number(ux1)
     call random_number(uy1)
     call random_number(uz1)
  endif

  !modulation of the random noise
  do k=1,xsize(3)
     do j=1,xsize(2)
        do i=1,xsize(1)
           ux1(i,j,k)=noise*(two*ux1(i,j,k)-one)
           uy1(i,j,k)=noise*(two*uy1(i,j,k)-one)
           uz1(i,j,k)=noise*(two*uz1(i,j,k)-one)
        enddo
     enddo
  enddo

  !initial velocity profile
  do k=1,xsize(3)
     do j=1,xsize(2)
        if (istret.eq.0) y=real(j+xstart(2)-1-1,mytype)*dy
        if (istret.ne.0) y=yp(j+xstart(2)-1)
        if (y.lt.yly-two) then
           do i=1,xsize(1)
              ux1(i,j,k) = zero
              uy1(i,j,k) = zero
              uz1(i,j,k) = zero
           enddo
        else
           do i=1,xsize(1)
              ux1(i,j,k)=(ux1(i,j,k)+one)*(one-(y+one-yly)**two)
              uy1(i,j,k)=(uy1(i,j,k))*(one-(y+one-yly)**two)
              uz1(i,j,k)=(uz1(i,j,k))*(one-(y+one-yly)**two)
           enddo
        endif
     enddo
  enddo

  !INIT FOR G AND U=MEAN FLOW + NOISE
  do k=1,xsize(3)
     do j=1,xsize(2)
        do i=1,xsize(1)
           ux1(i,j,k)=ux1(i,j,k)+bxx1(j,k)
           uy1(i,j,k)=uy1(i,j,k)+bxy1(j,k)
           uz1(i,j,k)=uz1(i,j,k)+bxz1(j,k)
           dux1(i,j,k,1)=ux1(i,j,k)
           duy1(i,j,k,1)=uy1(i,j,k)
           duz1(i,j,k,1)=uz1(i,j,k)
           dux1(i,j,k,2)=dux1(i,j,k,1)
           duy1(i,j,k,2)=duy1(i,j,k,1)
           duz1(i,j,k,2)=duz1(i,j,k,1)
        enddo
     enddo
  enddo

#ifdef DEBG
  if (nrank .eq. 0) print *,'# init end ok'
#endif

  return
end subroutine init
!********************************************************************
subroutine boundary_conditions (ux,uy,uz,phi,ep1)

  USE param
  USE variables
  USE decomp_2d

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,ep1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi

  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: gx

  ux = ux*(one-ep1)
  call transpose_x_to_y(ux,gx)
  call channel(gx,(two/three)*(two/yly))
  call transpose_y_to_x(gx,ux)

  return
end subroutine boundary_conditions
#ifdef POST
!********************************************************************
module post_processing

  USE decomp_2d
  USE variables
  USE param
  USE flow_type

  implicit none
  !
  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character
  !
  !probes
  integer :: nprobes, ntimes1, ntimes2
  integer, save, allocatable, dimension(:) :: rankprobes, nxprobes, nyprobes, nzprobes
  !
  real(mytype),save,allocatable,dimension(:,:) :: usum,vsum,wsum,uusum,uvsum,uwsum,vvsum,vwsum,wwsum

contains

  !############################################################################
  subroutine init_post(ep1)

    USE MPI

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ep1
    real(mytype) :: x, xprobes, yprobes, zprobes
    integer :: i,j,k,code
    character :: a

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init_post start'
#endif

    allocate(usum(zsize(1),zsize(2)),vsum(zsize(1),zsize(2)),wsum(zsize(1),zsize(2)))
    allocate(uusum(zsize(1),zsize(2)),uvsum(zsize(1),zsize(2)),uwsum(zsize(1),zsize(2)))
    allocate(vvsum(zsize(1),zsize(2)),vwsum(zsize(1),zsize(2)),wwsum(zsize(1),zsize(2)))
    usum=zero;vsum=zero;wsum=zero
    uusum=zero;uvsum=zero;uwsum=zero
    vvsum=zero;vwsum=zero;wwsum=zero
    ntimes1 = 0
    ntimes2 = 0
    nprobes  = 0

    !probes
    !WORK X-PENCILS
    open(10,file='probes.prm',status='unknown',form='formatted')
    read (10,*) nprobes
    read (10,*) a
    if (nprobes .gt. 0) then
       allocate(nxprobes(nprobes), nyprobes(nprobes), nzprobes(nprobes), rankprobes(nprobes))
       rankprobes(:)=0
       do i=1, nprobes
          read (10,*) xprobes, yprobes, zprobes
          !x
          if (nclx) then
             nxprobes(i)=int(xprobes/dx)
          else
             nxprobes(i)=int(xprobes/dx+1)
          end if
          !y
          if (ncly) then
             nyprobes(i)=int(yprobes/dy)
          else
             nyprobes(i)=int(yprobes/dy+1)
          end if
          !z
          if (nclz) then
             nzprobes(i)=int(zprobes/dz)
          else
             nzprobes(i)=int(zprobes/dz+1)
          end if
          if       (xstart(1) .le. nxprobes(i) .and. nxprobes(i) .le. xend(1)) then
             if    (xstart(2) .le. nyprobes(i) .and. nyprobes(i) .le. xend(2)) then
                if (xstart(3) .le. nzprobes(i) .and. nzprobes(i) .le. xend(3)) then
                   rankprobes(i)=1
                endif
             endif
          endif
       enddo
    endif
    close(10)

#ifdef DEBG 
    if (nrank .eq. 0) print *,'# init_post ok'
#endif

  end subroutine init_post
  !############################################################################
  subroutine postprocessing(ux1,uy1,uz1,phi1,ep1) !By Felipe Schuch

    USE decomp_2d_io

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1
    !
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2, uy2, uz2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ux3, uy3, uz3, temp3
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) ::  uprime, vprime, wprime
    !
    real(mytype),dimension(zsize(1),zsize(2)) :: um, vm, wm
    !
    integer :: i,j,k
    character(len=30) :: filename

    if (itime.ge.initstats1) then

       call transpose_x_to_y(ux1,ux2)
       call transpose_x_to_y(uy1,uy2)
       call transpose_x_to_y(uz1,uz2)

       call transpose_y_to_z(ux2,ux3)
       call transpose_y_to_z(uy2,uy3)
       call transpose_y_to_z(uz2,uz3)

       call mean_plane_z(ux3,zsize(1),zsize(2),zsize(3),um)
       call mean_plane_z(uy3,zsize(1),zsize(2),zsize(3),vm)
       call mean_plane_z(uz3,zsize(1),zsize(2),zsize(3),wm)

       usum = usum + um
       vsum = vsum + vm
       wsum = wsum + wm

       ntimes1 = ntimes1 + 1

       if (itime.ge.initstats2) then

          uprime = zero
          vprime = zero
          wprime = zero

          do k=1,zsize(3)
             do j=1,zsize(2)
                do i=1,zsize(1)
                   uprime(i,j,k) = ux3(i,j,k) - usum(i,j)/real(ntimes1,mytype)
                   vprime(i,j,k) = uy3(i,j,k) - vsum(i,j)/real(ntimes1,mytype)
                   wprime(i,j,k) = uz3(i,j,k) - wsum(i,j)/real(ntimes1,mytype)
                enddo
             enddo
          enddo

          do k=1,zsize(3)
             do j=1,zsize(2)
                do i=1,zsize(1)
                   uusum(i,j) = uusum(i,j) + uprime(i,j,k)*uprime(i,j,k)/real(nz,mytype)
                   uvsum(i,j) = uvsum(i,j) + uprime(i,j,k)*vprime(i,j,k)/real(nz,mytype)
                   uwsum(i,j) = uwsum(i,j) + uprime(i,j,k)*wprime(i,j,k)/real(nz,mytype)
                   vvsum(i,j) = vvsum(i,j) + vprime(i,j,k)*vprime(i,j,k)/real(nz,mytype)
                   vwsum(i,j) = vwsum(i,j) + vprime(i,j,k)*wprime(i,j,k)/real(nz,mytype)
                   wwsum(i,j) = wwsum(i,j) + wprime(i,j,k)*wprime(i,j,k)/real(nz,mytype)
                enddo
             enddo
          enddo

          ntimes2 = ntimes2 + 1
       endif

       if (mod(itime,imodulo).eq.0) then !write results
          temp3(:,:,1) = usum/real(ntimes1,mytype)
          write(filename,"('./data/usum',I4.4)") itime/imodulo
          call decomp_2d_write_plane(3,temp3,3,1,filename)
          temp3(:,:,1) = vsum/real(ntimes1,mytype)
          write(filename,"('./data/vsum',I4.4)") itime/imodulo
          call decomp_2d_write_plane(3,temp3,3,1,filename)
          temp3(:,:,1) = wsum/real(ntimes1,mytype)
          write(filename,"('./data/wsum',I4.4)") itime/imodulo
          call decomp_2d_write_plane(3,temp3,3,1,filename)
          temp3(:,:,1) = uusum/real(ntimes2,mytype)
          write(filename,"('./data/uusum',I4.4)") itime/imodulo
          call decomp_2d_write_plane(3,temp3,3,1,filename)
          temp3(:,:,1) = uvsum/real(ntimes2,mytype)
          write(filename,"('./data/uvsum',I4.4)") itime/imodulo
          call decomp_2d_write_plane(3,temp3,3,1,filename)
          temp3(:,:,1) = uwsum/real(ntimes2,mytype)
          write(filename,"('./data/uwsum',I4.4)") itime/imodulo
          call decomp_2d_write_plane(3,temp3,3,1,filename)
          temp3(:,:,1) = vvsum/real(ntimes2,mytype)
          write(filename,"('./data/vvsum',I4.4)") itime/imodulo
          call decomp_2d_write_plane(3,temp3,3,1,filename)
          temp3(:,:,1) = vwsum/real(ntimes2,mytype)
          write(filename,"('./data/vwsum',I4.4)") itime/imodulo
          call decomp_2d_write_plane(3,temp3,3,1,filename)
          temp3(:,:,1) = wwsum/real(ntimes2,mytype)
          write(filename,"('./data/wwsum',I4.4)") itime/imodulo
          call decomp_2d_write_plane(3,temp3,3,1,filename)
       endif
    endif
    return
  end subroutine postprocessing
  !############################################################################
  subroutine write_probes(ux1,uy1,uz1,phi1) !By Felipe Schuch

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ux1, uy1, uz1
    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),nphi) :: phi1

    integer :: i
    character(len=30) :: filename
    FS = 1+3+nphi !Number of columns
    write(fileformat, '( "(",I4,"(E14.6),A)" )' ) FS
    FS = FS*14+1  !Line width

    do i=1, nprobes
       if (rankprobes(i) .eq. 1) then
          write(filename,"('./out/probe',I4.4)") i
          open(67,file=trim(filename),status='unknown',form='formatted'&
               ,access='direct',recl=FS)
          write(67,fileformat,rec=itime) t,&                         !1
               ux1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !2
               uy1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !3
               uz1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !4
               phi1(nxprobes(i),nyprobes(i),nzprobes(i),:),&         !nphi
               NL                                                    !+1
          close(67)
       endif
    enddo

  end subroutine write_probes
  !############################################################################
end module post_processing
#endif

