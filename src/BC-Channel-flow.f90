!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!!        FILE: BC-Channel-flow.f90
!!!      AUTHOR: ??
!!!    MODIFIED: Paul Bartholomew
!!! DESCRIPTION: This module describes the channel flow.
!!!   CHANGELOG: [2019-02-19] Making module private by default
!!               [2019-02-19] Turning file into a module
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module channel

  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  !probes
  integer, save :: nprobes, ntimes1, ntimes2
  integer, save, allocatable, dimension(:) :: rankprobes, nxprobes, nyprobes, nzprobes

  real(mytype),save,allocatable,dimension(:) :: usum,vsum,wsum,uusum,uvsum,uwsum,vvsum,vwsum,wwsum

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_channel, boundary_conditions_channel, postprocessing_channel, momentum_forcing_channel

contains

  subroutine init_channel (ux1,uy1,uz1,ep1,phi1,dux1,duy1,duz1,dphi1)

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime,numscalar) :: dphi1

    real(mytype) :: y,r,um,r3,x,z,h,ct
    real(mytype) :: cx0,cy0,cz0,hg,lg
    integer :: k,j,i,ijk,fh,ierror,ii,is,it,code
    integer (kind=MPI_OFFSET_KIND) :: disp

    integer, dimension (:), allocatable :: seed

    if (iscalar==1) then

       phi1(:,:,:,:) = zero !change as much as you want

       do is = 1, numscalar
          dphi1(:,:,:,1,is) = phi1(:,:,:,is)
          do it = 2, ntime
             dphi1(:,:,:,it,is) = dphi1(:,:,:,it - 1,is)
          enddo
       enddo

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
             ux1(i,j,k)=init_noise*um*(two*ux1(i,j,k)-one)+one-y*y
             uy1(i,j,k)=init_noise*um*(two*uy1(i,j,k)-one)
             uz1(i,j,k)=init_noise*um*(two*uz1(i,j,k)-one)
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
             dux1(i,j,k,1)=ux1(i,j,k)
             duy1(i,j,k,1)=uy1(i,j,k)
             duz1(i,j,k,1)=uz1(i,j,k)
             do is = 2, ntime
                dux1(i,j,k,is)=dux1(i,j,k,is - 1)
                duy1(i,j,k,is)=duy1(i,j,k,is - 1)
                duz1(i,j,k,is)=duz1(i,j,k,is - 1)
             enddo
          enddo
       enddo
    enddo

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init end ok'
#endif

    return
  end subroutine init_channel
  !********************************************************************
  subroutine boundary_conditions_channel (ux,uy,uz,phi)

    USE param
    USE variables
    USE decomp_2d

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
!!$  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ut

    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: gx
    real(mytype) :: x, y, z
    integer :: i, j, k, is

    call transpose_x_to_y(ux,gx)
    call channel_flrt(gx,two/three)
    call transpose_y_to_x(gx,ux)

    if (iscalar.ne.0) then
       if (nclxS1.eq.2) then
          i = 1
          phi(i,:,:,:) = zero
       endif
       if (nclxSn.eq.2) then
          i = xsize(1)
          phi(i,:,:,:) = phi(i - 1,:,:,:)
       endif

       if ((nclyS1.eq.2).and.(xstart(2).eq.1)) then
          !! Generate a hot patch on bottom boundary
          do k = 1, xsize(3)
             z = real(k + xstart(3) - 2, mytype) * dz - half * zlz
             if (abs(z).lt.zlz/four) then
                j = 1
                do i = 1, xsize(1)
                   x = real(i + xstart(1) - 2, mytype) * dx
                   if ((x.gt.0.1*xlx).and.(x.lt.0.3*xlx)) then
                      do is = 1, numscalar
                         phi(i, j, k, is) = one
                      enddo
                   else
                      do is = 1, numscalar
                         phi(i, j, k, is) = zero
                      enddo
                   endif
                enddo
             endif
          enddo
       endif
    endif

    return
  end subroutine boundary_conditions_channel

  !********************************************************************
  !
  subroutine channel_flrt (ux,constant)
    !
    !********************************************************************

    USE decomp_2d
    USE decomp_2d_poisson
    USE variables
    USE param
    USE var
    USE MPI

    implicit none

    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux
    real(mytype) :: constant

    integer :: j,i,k,code
    real(mytype) :: can,ut3,ut,ut4

    ut3=zero
    do k=1,ysize(3)
       do i=1,ysize(1)
          ut=zero
          do j=1,ny-1
             if (istret.eq.0) then
                ut=ut+dy*(ux(i,j+1,k)-half*(ux(i,j+1,k)-ux(i,j,k)))
             else
                ut=ut+(yp(j+1)-yp(j))*(ux(i,j+1,k)-half*(ux(i,j+1,k)-ux(i,j,k)))
             endif
          enddo
          ut=ut/yly
          ut3=ut3+ut
       enddo
    enddo
    ut3=ut3/(real(nx*nz,mytype))

    call MPI_ALLREDUCE(ut3,ut4,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

    can=-(constant-ut4)

    if (nrank==0) print *,nrank,'UT',ut4,can

    do k=1,ysize(3)
       do i=1,ysize(1)
          do j=2,ny-1
             ux(i,j,k)=ux(i,j,k)-can
          enddo
       enddo
    enddo

    return
  end subroutine channel_flrt
  !********************************************************************

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
  subroutine postprocessing_channel(ux1,uy1,uz1,pp3,phi1,ep1) !By Felipe Schuch

    USE MPI
    USE decomp_2d
    USE decomp_2d_io
    USE var, only : umean,vmean,wmean,pmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
    USE var, only : ta1, pp1, di1
    USE var, only : ppi3, dip3
    USE var, only : pp2, ppi2, dip2
    
    USE var, ONLY : nxmsize, nymsize, nzmsize
    USE param, ONLY : npress

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(in) :: pp3
    character(len=30) :: filename

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

       !pmean

       call interzpv(ppi3,pp3,dip3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
            (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
       call transpose_z_to_y(ppi3,pp2,ph3) !nxm nym nz
       call interypv(ppi2,pp2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
            (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
       call transpose_y_to_x(ppi2,pp1,ph2) !nxm ny nz
       call interxpv(ta1,pp1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
            nxmsize,xsize(1),xsize(2),xsize(3),1)
       
       call fine_to_coarseS(1,pp1,tmean)
       pmean(:,:,:)=pmean(:,:,:)+tmean(:,:,:)

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
          write(filename,"('pmean.dat',I7.7)") itime
          call decomp_2d_write_one(1,pmean,filename,1)
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
             write(filename,"('pmean.dat',I7.7)") itime-icheckpoint
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
  end subroutine postprocessing_channel
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  SUBROUTINE: momentum_forcing
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Applies rotation for t < spinup_time.
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE momentum_forcing_channel(dux1, duy1, ux1, uy1)

    IMPLICIT NONE

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ux1, uy1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1

    if (itime.lt.spinup_time) then
       if (nrank==0) print *,'Rotating turbulent channel at speed ',wrotation
       dux1(:,:,:,1) = dux1(:,:,:,1) - wrotation*uy1(:,:,:)
       duy1(:,:,:,1) = duy1(:,:,:,1) + wrotation*ux1(:,:,:)
    endif

  ENDSUBROUTINE momentum_forcing_channel

end module channel
