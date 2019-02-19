
!********************************************************************
subroutine geomcomplex(epsi,nxi,nxf,ny,nyi,nyf,nzi,nzf,dx,yp,dz,remp)
  use decomp_2d, only : mytype
  use param, only : zero, one, two
  use ibm
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

  do k=nzi,nzf
     do j=nyi,nyf
        ym=yp(j)
        do i=nxi,nxf
           xm=real(i-1,mytype)*dx
           r=sqrt((xm-cex)**two+(ym-cey)**two)
           if (r-ra .gt. zeromach) cycle
           epsi(i,j,k)=remp
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
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

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
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

  call random_number(bxo)
  call random_number(byo)
  call random_number(bzo)
  do k=1,xsize(3)
     do j=1,xsize(2)
        bxx1(j,k)=one+bxo(j,k)*inflow_noise
        bxy1(j,k)=zero+byo(j,k)*inflow_noise
        bxz1(j,k)=zero+bzo(j,k)*inflow_noise
     enddo
  enddo

  if (iscalar.eq.1) then
     do is=1, numscalar
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
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
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
subroutine init (ux1,uy1,uz1,ep1,phi1,phis1,phiss1)

  USE decomp_2d
  USE decomp_2d_io
  USE variables
  USE param
  USE MPI

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1,phis1,phiss1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1

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
              ux1(i,j,k)=init_noise*(ux1(i,j,k)-0.5)
              uy1(i,j,k)=init_noise*(uy1(i,j,k)-0.5)
              uz1(i,j,k)=init_noise*(uz1(i,j,k)-0.5)
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
module post_processing

  USE decomp_2d
  USE variables
  USE param

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
    USE decomp_2d
    USE decomp_2d_io
    USE var, only : umean,vmean,wmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
    USE var, only : uvisu
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    character(len=30) :: filename

    if (itime.ge.ioutput) then
    
       !x-derivatives
       call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
       call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
       call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
       !y-derivatives
       call transpose_x_to_y(ux1,td2)
       call transpose_x_to_y(uy1,te2)
       call transpose_x_to_y(uz1,tf2)
       call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
       call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
       call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
       !!z-derivatives
       call transpose_y_to_z(td2,td3)
       call transpose_y_to_z(te2,te3)
       call transpose_y_to_z(tf2,tf3)
       call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
       call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
       call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
       !!all back to x-pencils
       call transpose_z_to_y(ta3,td2)
       call transpose_z_to_y(tb3,te2)
       call transpose_z_to_y(tc3,tf2)
       call transpose_y_to_x(td2,tg1)
       call transpose_y_to_x(te2,th1)
       call transpose_y_to_x(tf2,ti1)
       call transpose_y_to_x(ta2,td1)
       call transpose_y_to_x(tb2,te1)
       call transpose_y_to_x(tc2,tf1)
       !du/dx=ta1 du/dy=td1 and du/dz=tg1
       !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
       !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1
       !############################################################################
       !VORTICITY
       di1(:,:,:)=0.
       di1(:,:,:)=sqrt((tf1(:,:,:)-th1(:,:,:))**2+(tg1(:,:,:)-tc1(:,:,:))**2+&
            (tb1(:,:,:)-td1(:,:,:))**2)
       uvisu=0.
       call fine_to_coarseV(1,di1,uvisu)
990    format('vort',I3.3)
       write(filename, 990) itime/ioutput
       call decomp_2d_write_one(1,uvisu,filename,2)
       !############################################################################

       !############################################################################
       !VELOCITY
       uvisu=0.
       call fine_to_coarseV(1,ux1,uvisu)
993    format('ux',I3.3)
       write(filename, 993) itime/ioutput
       call decomp_2d_write_one(1,uvisu,filename,2)
       uvisu=0.
       call fine_to_coarseV(1,uy1,uvisu)
994    format('uy',I3.3)
       write(filename, 994) itime/ioutput
       call decomp_2d_write_one(1,uvisu,filename,2)
       uvisu=0.
       call fine_to_coarseV(1,uz1,uvisu)
995    format('uz',I3.3)
       write(filename, 995) itime/ioutput
       call decomp_2d_write_one(1,uvisu,filename,2)
       !############################################################################

       !############################################################################
       !PRESSURE       
       !NOT READY
       
 endif   
    
    if (itime.ge.initstat) then
       !umean=ux1
       call fine_to_coarseS(1,ux1,tmean)
       umean(:,:,:)=umean(:,:,:)+tmean(:,:,:)
       
       !vmean=uy1
       call fine_to_coarseS(1,uy1,tmean)
       vmean(:,:,:)=vmean(:,:,:)+tmean(:,:,:)

       !wmean=uz1
       call fine_to_coarseS(1,uz1,tmean)
       wmean(:,:,:)=wmean(:,:,:)+tmean(:,:,:)

       !uumean=ux1*ux1
       ta1(:,:,:)=ux1(:,:,:)*ux1(:,:,:)
       call fine_to_coarseS(1,ta1,tmean)
       uumean(:,:,:)=uumean(:,:,:)+tmean(:,:,:)

       !vvmean=uy1*uy1
       ta1(:,:,:)=uy1(:,:,:)*uy1(:,:,:)
       call fine_to_coarseS(1,ta1,tmean)
       vvmean(:,:,:)=vvmean(:,:,:)+tmean(:,:,:)

       !wwmean=uz1*uz1
       ta1(:,:,:)=uz1(:,:,:)*uz1(:,:,:)
       call fine_to_coarseS(1,ta1,tmean)
       wwmean(:,:,:)=wwmean(:,:,:)+tmean(:,:,:)

       !uvmean=ux1*uy1
       ta1(:,:,:)=ux1(:,:,:)*uy1(:,:,:)
       call fine_to_coarseS(1,ta1,tmean)
       uvmean(:,:,:)=uvmean(:,:,:)+tmean(:,:,:)

       !uwmean=ux1*uz1
       ta1(:,:,:)=ux1(:,:,:)*uz1(:,:,:)
       call fine_to_coarseS(1,ta1,tmean)
       uwmean(:,:,:)=uwmean(:,:,:)+tmean(:,:,:)

       !vwmean=uy1*uz1
       ta1(:,:,:)=uy1(:,:,:)*uz1(:,:,:)
       call fine_to_coarseS(1,ta1,tmean)
       vwmean(:,:,:)=vwmean(:,:,:)+tmean(:,:,:)

       if (mod(itime,icheckpoint)==0) then
          if (nrank==0) print *,'===========================================================<<<<<'
          if (nrank==0) print *,'Writing stat file',itime
          write(filename,"('umean.dat',I7.7)") itime
          call decomp_2d_write_one(1,umean,filename,1)
          write(filename,"('vmean.dat',I7.7)") itime
          call decomp_2d_write_one(1,vmean,filename,1)
          write(filename,"('wmean.dat',I7.7)") itime
          call decomp_2d_write_one(1,wmean,filename,1)
          write(filename,"('uumean.dat',I7.7)") itime
          call decomp_2d_write_one(1,uumean,filename,1)
          write(filename,"('vvmean.dat',I7.7)") itime
          call decomp_2d_write_one(1,vvmean,filename,1)
          write(filename,"('wwmean.dat',I7.7)") itime
          call decomp_2d_write_one(1,wwmean,filename,1)
          write(filename,"('uvmean.dat',I7.7)") itime
          call decomp_2d_write_one(1,uvmean,filename,1)
          write(filename,"('uwmean.dat',I7.7)") itime
          call decomp_2d_write_one(1,uwmean,filename,1)
          write(filename,"('vwmean.dat',I7.7)") itime
          call decomp_2d_write_one(1,vwmean,filename,1)  
          if (nrank==0) print *,'write stat done!'
          if (nrank==0) print *,'===========================================================<<<<<'
          if (nrank==0) then
             write(filename,"('umean.dat',I7.7)") itime-icheckpoint
             call system ("rm " //filename)
             write(filename,"('vmean.dat',I7.7)") itime-icheckpoint
             call system ("rm " //filename)
             write(filename,"('wmean.dat',I7.7)") itime-icheckpoint
             call system ("rm " //filename)
             write(filename,"('uumean.dat',I7.7)") itime-icheckpoint
             call system ("rm " //filename)
             write(filename,"('vvmean.dat',I7.7)") itime-icheckpoint
             call system ("rm " //filename)
             write(filename,"('wwmean.dat',I7.7)") itime-icheckpoint
             call system ("rm " //filename)
             write(filename,"('uvmean.dat',I7.7)") itime-icheckpoint
             call system ("rm " //filename)
             write(filename,"('uwmean.dat',I7.7)") itime-icheckpoint
             call system ("rm " //filename)
             write(filename,"('vwmean.dat',I7.7)") itime-icheckpoint
             call system ("rm " //filename)
          endif
       endif

    endif

    return
  end subroutine postprocessing
  !############################################################################
  subroutine write_probes(ux1,uy1,uz1,phi1) !By Felipe Schuch

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ux1, uy1, uz1
    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),numscalar) :: phi1

    integer :: i
    character(len=30) :: filename
    FS = 1+3+numscalar !Number of columns
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
               phi1(nxprobes(i),nyprobes(i),nzprobes(i),:),&         !numscalar
               NL                                                    !+1
          close(67)
       endif
    enddo

  end subroutine write_probes
  !############################################################################
end module post_processing

