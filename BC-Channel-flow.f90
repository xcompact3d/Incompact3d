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

  open(10,file='BC-Channel-flow.prm',status='old',form='formatted')
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
  xlx = xlx*pi
  read (10,*) yly
  read (10,*) zlz
  zlz = zlz*pi
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
  read (10,*) isave
  read (10,*) imodulo
  read (10,*) wrotation
  read (10,*) irotation
  read (10,*) initstats1
  read (10,*) initstats2
  read (10,*) a !
  read (10,*) a ! NUMERICAL DISSIPATION
  read (10,*) a !
  read (10,*) jLES
  read (10,*) fpi2

  if (nrank==0) then
     print *,'==================Turbulent channel flow==================='
     write(*,"(' irotation          : ',I15)") irotation
     write(*,"(' wrotation          : ',F15.8)") wrotation
     write(*,"(' initstats1         : ',I15)") initstats1
     write(*,"(' initstats2         : ',I15)") initstats2
     print *,'==========================================================='
  endif
  return
end subroutine ft_parameter
!********************************************************************
subroutine init (ux1,uy1,uz1,ep1,phi1,gx1,gy1,gz1,phis1,hx1,hy1,hz1,phiss1)

  USE decomp_2d
  USE decomp_2d_io
  USE variables
  USE param
  USE MPI

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gx1,gy1,gz1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: hx1,hy1,hz1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1,phis1,phiss1

  real(mytype) :: y,r,um,r3,x,z,h,ct
  real(mytype) :: cx0,cy0,cz0,hg,lg
  integer :: k,j,i,ijk,fh,ierror,ii,is,code
  integer (kind=MPI_OFFSET_KIND) :: disp

  integer, dimension (:), allocatable :: seed

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

  !modulation of the random noise + initial velocity profile
  do k=1,xsize(3)
     do j=1,xsize(2)
        if (istret.eq.0) y=real(j+xstart(2)-1-1,mytype)*dy-yly/two
        if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/two
        um=exp(-zptwo*y*y)
        do i=1,xsize(1)
           ux1(i,j,k)=noise*um*(two*ux1(i,j,k)-one)+one-y*y
           uy1(i,j,k)=noise*um*(two*uy1(i,j,k)-one)
           uz1(i,j,k)=noise*um*(two*uz1(i,j,k)-one)
        enddo
     enddo
  enddo

  !INIT FOR G AND U=MEAN FLOW + NOISE
  do k=1,xsize(3)
     do j=1,xsize(2)
        do i=1,xsize(1)
           ux1(i,j,k)=ux1(i,j,k)+bxx1(j,k)
           uy1(i,j,k)=uy1(i,j,k)+bxy1(j,k)
           uz1(i,j,k)=uz1(i,j,k)+bxz1(j,k)
           gx1(i,j,k)=ux1(i,j,k)
           gy1(i,j,k)=uy1(i,j,k)
           gz1(i,j,k)=uz1(i,j,k)
           hx1(i,j,k)=gx1(i,j,k)
           hy1(i,j,k)=gy1(i,j,k)
           hz1(i,j,k)=gz1(i,j,k)
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
!!$  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ut

  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: gx

  call transpose_x_to_y(ux,gx)
  call channel(gx,two/three)
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
  integer, save :: nprobes, ntimes1, ntimes2
  integer, save, allocatable, dimension(:) :: rankprobes, nxprobes, nyprobes, nzprobes
  !
  real(mytype),save,allocatable,dimension(:) :: usum,vsum,wsum,uusum,uvsum,uwsum,vvsum,vwsum,wwsum

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

    allocate(usum(ysize(2)),vsum(ysize(2)),wsum(ysize(2)))
    allocate(uusum(ysize(2)),uvsum(ysize(2)),uwsum(ysize(2)))
    allocate(vvsum(ysize(2)),vwsum(ysize(2)),wwsum(ysize(2)))
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

    USE MPI

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1
    !
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2, uy2, uz2
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) ::  uprime, vprime, wprime
    !
    real(mytype),dimension(ysize(2)) :: um, vm, wm
    real(mytype),dimension(ysize(2)) :: um1,vm1,wm1
    real(mytype),dimension(ysize(2)) :: uum, uvm, uwm, vvm, vwm, wwm
    real(mytype),dimension(ysize(2)) :: uum1,uvm1,uwm1,vvm1,vwm1,wwm1
    !
    integer :: i,j,k,code
    character(len=30) :: filename

    if (itime.ge.initstats1) then

       call transpose_x_to_y(ux1,ux2)
       call transpose_x_to_y(uy1,uy2)
       call transpose_x_to_y(uz1,uz2)

       um = zero; vm = zero; wm = zero
       um1= zero; vm1= zero; wm1= zero
       do j=1,ysize(2)
          um(j) = sum(ux2(:,j,:))
          vm(j) = sum(uy2(:,j,:))
          wm(j) = sum(uz2(:,j,:))
       enddo

       call MPI_ALLREDUCE(um,um1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       call MPI_ALLREDUCE(vm,vm1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       call MPI_ALLREDUCE(wm,wm1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)

       usum = usum + um1/real(nx*nz,mytype)
       vsum = vsum + vm1/real(nx*nz,mytype)
       wsum = wsum + wm1/real(nx*nz,mytype)

       ntimes1 = ntimes1 + 1

       if (itime.ge.initstats2) then

          do k=1,ysize(3)
             do j=1,ysize(2)
                do i=1,ysize(1)
                   uprime(i,j,k) = ux2(i,j,k) - usum(j)/real(ntimes1,mytype)
                   vprime(i,j,k) = uy2(i,j,k) - vsum(j)/real(ntimes1,mytype)
                   wprime(i,j,k) = uz2(i,j,k) - wsum(j)/real(ntimes1,mytype)
                enddo
             enddo
          enddo

          uum=zero ;uvm=zero ;uwm=zero ;vvm=zero ;vwm=zero ;wwm=zero
          uum1=zero;uvm1=zero;uwm1=zero;vvm1=zero;vwm1=zero;wwm1=zero
          do k=1,ysize(3)
             do j=1,ysize(2)
                do i=1,ysize(1)
                   uum(j) = uum(j) + uprime(i,j,k)*uprime(i,j,k)
                   uvm(j) = uvm(j) + uprime(i,j,k)*vprime(i,j,k)
                   uwm(j) = uwm(j) + uprime(i,j,k)*wprime(i,j,k)
                   vvm(j) = vvm(j) + vprime(i,j,k)*vprime(i,j,k)
                   vwm(j) = vwm(j) + vprime(i,j,k)*wprime(i,j,k)
                   wwm(j) = wwm(j) + wprime(i,j,k)*wprime(i,j,k)
                enddo
             enddo
          enddo

          call MPI_ALLREDUCE(uum,uum1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
          call MPI_ALLREDUCE(uvm,uvm1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
          call MPI_ALLREDUCE(uwm,uwm1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
          call MPI_ALLREDUCE(vvm,vvm1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
          call MPI_ALLREDUCE(vwm,vwm1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
          call MPI_ALLREDUCE(wwm,wwm1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)

          uusum = uusum + uum1/real(nx*nz,mytype)
          uvsum = uvsum + uvm1/real(nx*nz,mytype)
          uwsum = uwsum + uwm1/real(nx*nz,mytype)
          vvsum = vvsum + vvm1/real(nx*nz,mytype)
          vwsum = vwsum + vwm1/real(nx*nz,mytype)
          wwsum = wwsum + wwm1/real(nx*nz,mytype)

          ntimes2 = ntimes2 + 1
       endif

       if (mod(itime,imodulo).eq.0) then !write results
          if (nrank.eq.0) then
             write(filename,"('./out/stats',I4.4)") itime/imodulo
             open(67,file=trim(filename),status='unknown',form='formatted')
             do j=1,ysize(2)
                write(67,'(10E16.8)') yp(j),&
                     usum(j)/real(ntimes1,mytype),&
                     vsum(j)/real(ntimes1,mytype),&
                     wsum(j)/real(ntimes1,mytype),&
                     uusum(j)/real(ntimes2,mytype),&
                     uvsum(j)/real(ntimes2,mytype),&
                     uwsum(j)/real(ntimes2,mytype),&
                     vvsum(j)/real(ntimes2,mytype),&
                     vwsum(j)/real(ntimes2,mytype),&
                     wwsum(j)/real(ntimes2,mytype)
             enddo
             close(67)
          endif
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
