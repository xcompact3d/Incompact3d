module flow_type
  use decomp_2d, only : mytype
  
  integer :: ncil
  real(mytype) :: ra
  real(mytype),allocatable,dimension(:) :: cex,cey

end module flow_type

#ifdef FORCES
#include "forces.f90"
#endif

subroutine ft_parameter(arg)

  USE param
  USE variables
  USE flow_type
  USE complex_geometry
#ifdef FORCES
  USE forces
#endif
  USE decomp_2d, only : nrank
  implicit none

  logical,intent(in) :: arg
  integer :: i
  character :: a

  nclx1 = 2 !Boundary condition in x=0  (0: Periodic, 1:Free-slip, 2: Dirichlet)
  nclxn = 2 !Boundary condition in x=Lx (0: Periodic, 1:Free-slip, 2: Dirichlet)
  ncly1 = 1 !Boundary condition in y=0  (0: Periodic, 1:Free-slip, 2: Dirichlet)
  nclyn = 1 !Boundary condition in y=Ly (0: Periodic, 1:Free-slip, 2: Dirichlet)
  nclz1 = 0 !Boundary condition in z=0  (0: Periodic, 1:Free-slip, 2: Dirichlet)
  nclzn = 0 !Boundary condition in z=Lz (0: Periodic, 1:Free-slip, 2: Dirichlet)

  open(10,file='BC-Cylinder.prm',status='unknown',form='formatted')
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
  read (10,*) angle
  read (10,*) u1
  read (10,*) u2
  read (10,*) noise
  read (10,*) noise1
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
  read (10,*) a !
  read (10,*) a ! NUMERICAL DISSIPATION
  read (10,*) a !
  read (10,*) jLES
  read (10,*) fpi2
  read (10,*) a !
  read (10,*) a ! INCOMPACT 3D Body (old school)
  read (10,*) a !
  read (10,*) ivirt
  read (10,*) ncil
  allocate(cex(ncil),cey(ncil))
  read (10,*) cex
  read (10,*) cey
  read (10,*) ra
  read (10,*) a !
  read (10,*) a ! INCOMPACT 3D Forcing with Lagrangian Polynomials
  read (10,*) a !
  read (10,*) ilag
  read (10,*) npif 
  read (10,*) izap  
  read (10,*) nraf
  read (10,*) nobjmax
#ifdef FORCES
  read (10,*) a !
  read (10,*) a !INCOMPACT 3D Forces - Drag and Lift coefficients
  read (10,*) a !
  read (10,*) nvol
  allocate(xld(nvol),xrd(nvol),yld(nvol),yud(nvol))
  read (10,*) xld
  read (10,*) xrd
  read (10,*) yld
  read (10,*) yud
#endif
  close(10)

  if (nrank==0) then
     print *,'=======================Cylinder============================'
     do i=1,ncil
        write(*,"(' Cylinder           : #',I1)") i
        write(*,"(' cex, cey, ra       : (',F6.2,',',F6.2,',',F6.2,')')") cex(i), cey(i), ra
     enddo
     print *,'==========================================================='
  endif
  return
end subroutine ft_parameter
!********************************************************************
subroutine geomcomplex(epsi,nxi,nxf,ny,nyi,nyf,nzi,nzf,dx,yp,dz,remp)
  use flow_type, only : cex,cey,ra,ncil
  use decomp_2d, only : mytype
  use param, only : zero, one, two
  implicit none
  !
  real(mytype),dimension(nxi:nxf,nyi:nyf,nzi:nzf) :: epsi
  real(mytype),dimension(ny) :: yp
  integer                    :: nxi,nxf,ny,nyi,nyf,nzi,nzf
  real(mytype)               :: dx,dz
  real(mytype)               :: remp
  integer                    :: i,ic,j,k
  real(mytype)               :: xm,ym,r
  real(mytype)               :: zeromach

  zeromach=one
  do while ((one + zeromach / two) .gt. one)
     zeromach = zeromach/two
  end do
  zeromach = 1.0e1*zeromach
  
  do ic=1, ncil
     do k=nzi,nzf
        do j=nyi,nyf
           ym=yp(j)
           do i=nxi,nxf
              xm=real(i-1,mytype)*dx
              r=sqrt((xm-cex(ic))**two+(ym-cey(ic))**two)
              if (r-ra .gt. zeromach) cycle
              epsi(i,j,k)=remp
           enddo
        enddo
     enddo
  enddo
  !
  return
end subroutine geomcomplex
!********************************************************************
subroutine boundary_conditions (ux,uy,uz,phi,ep1)

  USE param
  USE variables
  USE decomp_2d

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,ep1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi

  call inflow (ux,uy,uz,phi)
  call outflow (ux,uy,uz,phi)

  return
end subroutine boundary_conditions
!********************************************************************
subroutine inflow (ux,uy,uz,phi)

  USE param
  USE variables
  USE decomp_2d

  implicit none

  integer  :: j,k,is
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi

  call random_number(bxo)
  call random_number(byo)
  call random_number(bzo)
  do k=1,xsize(3)
     do j=1,xsize(2)
        bxx1(j,k)=one+bxo(j,k)*noise1
        bxy1(j,k)=zero+byo(j,k)*noise1
        bxz1(j,k)=zero+bzo(j,k)*noise1
     enddo
  enddo

  if (iscalar.eq.1) then
     do is=1, nphi
        do k=1,xsize(3)
           do j=1,xsize(2)
              phi(1,j,k,is)=cp(is)
           enddo
        enddo
     enddo
  endif

  return
end subroutine inflow
!********************************************************************
subroutine outflow (ux,uy,uz,phi)

  USE param
  USE variables
  USE decomp_2d
  USE MPI

  implicit none

  integer :: i,j,k,is,code
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi
  real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx,cz,uxmin,uxmax,uxmin1,uxmax1

  udx=one/dx; udy=one/dy; udz=one/dz; uddx=half/dx; uddy=half/dy; uddz=half/dz

  uxmax=-1609.
  uxmin=1609.
  do k=1,xsize(3)
     do j=1,xsize(2)
        if (ux(nx-1,j,k).gt.uxmax) uxmax=ux(nx-1,j,k)
        if (ux(nx-1,j,k).lt.uxmin) uxmin=ux(nx-1,j,k)
     enddo
  enddo

  call MPI_ALLREDUCE(uxmax,uxmax1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
  call MPI_ALLREDUCE(uxmin,uxmin1,1,real_type,MPI_MIN,MPI_COMM_WORLD,code)

  if (u1.eq.0) cx=(half*(uxmax1+uxmin1))*gdt(itr)*udx
  if (u1.eq.1) cx=uxmax1*gdt(itr)*udx
  if (u1.eq.2) cx=u2*gdt(itr)*udx    !works better

  do k=1,xsize(3)
     do j=1,xsize(2)
        bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
        bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
        bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
     enddo
  enddo

  if (iscalar==1) then
     if (u2.eq.0.) cx=(half*(uxmax1+uxmin1))*gdt(itr)*udx
     if (u2.eq.1.) cx=uxmax1*gdt(itr)*udx
     if (u2.eq.2.) cx=u2*gdt(itr)*udx    !works better
     do k=1,xsize(3)
        do j=1,xsize(2)
           phi(nx,j,k,:)=phi(nx,j,k,:)-cx*(phi(nx,j,k,:)-phi(nx-1,j,k,:))
        enddo
     enddo
  endif

  if (nrank==0) write(*,*) "Outflow velocity ux nx=n min max=",real(uxmin1,4),real(uxmax1,4)

  return
end subroutine outflow
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

  ux1=zero; uy1=zero; uz1=zero

  if (iin.ne.0) then
     call system_clock(count=code)
     if (iin.eq.2) code=0
     call random_seed(size = ii)
     call random_seed(put = code+63946*nrank*(/ (i - 1, i = 1, ii) /))

     call random_number(ux1)
     call random_number(uy1)
     call random_number(uz1)


     do k=1,xsize(3)
        do j=1,xsize(2)
           do i=1,xsize(1)
              ux1(i,j,k)=noise*(ux1(i,j,k)-0.5)
              uy1(i,j,k)=noise*(uy1(i,j,k)-0.5)
              uz1(i,j,k)=noise*(uz1(i,j,k)-0.5)
           enddo
        enddo
     enddo

     !modulation of the random noise
     do k=1,xsize(3)
        do j=1,xsize(2)
           if (istret.eq.0) y=(j+xstart(2)-1-1)*dy-yly/2.
           if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/2.
           um=exp(-zptwo*y*y)
           do i=1,xsize(1)
              ux1(i,j,k)=um*ux1(i,j,k)
              uy1(i,j,k)=um*uy1(i,j,k)
              uz1(i,j,k)=um*uz1(i,j,k)
           enddo
        enddo
     enddo
  endif

  !INIT FOR G AND U=MEAN FLOW + NOISE
  do k=1,xsize(3)
     do j=1,xsize(2)
        do i=1,xsize(1)
           ux1(i,j,k)=ux1(i,j,k)+one
           uy1(i,j,k)=uy1(i,j,k)
           uz1(i,j,k)=uz1(i,j,k)
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
