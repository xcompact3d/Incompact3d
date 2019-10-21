!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!!        FILE: BC-TBL.f90
!!!      AUTHOR: ??
!!!    MODIFIED: Arash Hamzehloo
!!! DESCRIPTION: This module describes the turbulent boundary layer.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module tbl

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
  PUBLIC :: init_tbl, boundary_conditions_tbl, postprocess_tbl

contains

  subroutine init_tbl (ux1,uy1,uz1,ep1,phi1)

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: y,r,um,r3,x,z,h,ct
    real(mytype) :: cx0,cy0,cz0,hg,lg
    integer :: k,j,i,fh,ierror,ii,is,it,code
    integer (kind=MPI_OFFSET_KIND) :: disp

    integer, dimension (:), allocatable :: seed

    if (iscalar==1) then

       phi1(:,:,:,:) = zero !change as much as you want

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
          um=exp(-16*y*y)
          do i=1,xsize(1)
             ux1(i,j,k)=init_noise*um*(two*ux1(i,j,k)-one)
             uy1(i,j,k)=init_noise*um*(two*uy1(i,j,k)-one)
             uz1(i,j,k)=init_noise*um*(two*uz1(i,j,k)-one)
          enddo
       enddo
    enddo

    call ecoule(ux1,uy1,uz1,phi1)

    !INIT FOR G AND U=MEAN FLOW + NOISE
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             ux1(i,j,k)=ux1(i,j,k)+bxx1(j,k)
             uy1(i,j,k)=uy1(i,j,k)+bxy1(j,k)
             uz1(i,j,k)=uz1(i,j,k)+bxz1(j,k)
          enddo
       enddo
    enddo

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init end ok'
#endif

    return
  end subroutine init_tbl
  !********************************************************************
  subroutine boundary_conditions_tbl (ux,uy,uz,phi)

    USE param
    USE variables
    USE decomp_2d

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    real(mytype) :: x, y, z, um
    real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx

    integer :: i, j, k, is

    !INFLOW

    call ecoule(ux,uy,uz,phi)

    call random_number(bxo)
    call random_number(byo)
    call random_number(bzo)

    if (iin.eq.1) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             if (istret.eq.0) y=(j+xstart(2)-1-1)*dy-1.
             if (istret.ne.0) y=yp(j+xstart(2)-1)-1.
             um=exp(-16*y*y)
             bxx1(j,k)=bxx1(j,k)+bxo(j,k)*inflow_noise*um*0.
             bxy1(j,k)=bxy1(j,k)+byo(j,k)*inflow_noise*um*0.
             bxz1(j,k)=bxz1(j,k)+bzo(j,k)*inflow_noise*um*0.
          enddo
       enddo
       !if (iscalar==1) then
       !do k=1,xsize(3)
       !do j=1,xsize(2)
       !bxo(j,k)=phi(xsize(1)/4,j,k)
       !enddo
       !enddo
       !endif
    endif

    !OUTFLOW

    udx=one/dx
    udy=one/dy
    udz=one/dz
    uddx=half/dx
    uddy=half/dy
    uddz=half/dz

    do k=1,xsize(3)
       do j=1,xsize(2)

          cx=ux(nx,j,k)*gdt(itr)*udx

          if (cx.LT.0.0) cx=0.0
          bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
          bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
          bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
       enddo
    enddo

    if (nclyn == 2) THEN
       do k = 1, xsize(3)
          do i = 1, xsize(1)
             byxn(i, k) = ux(i, xsize(2) - 1, k)
             byyn(i, k) = zero
             byzn(i, k) = uz(i, xsize(2) - 1, k)
          enddo
       enddo
    endif

    !SCALAR
    if (iscalar.ne.0) then
       if (nclxS1.eq.2) then
          i = 1
          phi(i,:,:,:) = zero
       endif
       if (nclxSn.eq.2) then
          i = xsize(1)
          phi(i,:,:,:) = phi(i - 1,:,:,:)
       endif
    endif

    return
  end subroutine boundary_conditions_tbl

  !********************************************************************
  subroutine ecoule(ux1,uy1,uz1,phi1)


    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: eta_bl, f_bl, g_bl, x_bl,h_bl
    real(mytype) :: delta_int, delta_eta, eps_eta


    real(mytype) :: x, y, z
    integer :: i, j, k, is

    do k=1,xsize(3)
       do j=1,xsize(2)
          if (istret.eq.0) y=(j+xstart(2)-1-1)*dy
          if (istret.ne.0) y=yp(j+xstart(2)-1)

          eta_bl=y*4.91/9.0

          !OLD POLYNOMIAL FITTING

          delta_eta=0.0
          eps_eta=0.0
          delta_int=0.2

          if (eta_bl .ge. (7.5/9.0)) then
             delta_eta=eta_bl-7.5/9.0
             eta_bl=7.5/9.0
             eps_eta=0.00015
          end if

          f_bl=1678.64209592595*eta_bl**14-11089.6925017429*eta_bl**13 &
               +31996.4350140670*eta_bl**12-52671.5249779799*eta_bl**11 &
               +54176.1691167667*eta_bl**10-35842.8204706097*eta_bl**9  &
               +15201.3088871240*eta_bl**8 -4080.17137935648*eta_bl**7  &
               +702.129634528103*eta_bl**6 -56.2063925805318*eta_bl**5  &
               -17.0181128273914*eta_bl**4 +0.819582894357566*eta_bl**3  &
               -0.0601348202321954*eta_bl**2 +2.98973991270405*eta_bl**1

          f_bl=f_bl+(1-exp(-delta_eta/delta_int))*eps_eta


          if (eta_bl .ge. (7.15/9.0)) then
             delta_int=0.8
             delta_eta=eta_bl-7.15/9.0
             eta_bl=7.15/9.0
             eps_eta=0.0005
          end if

          g_bl=4924.05284779754*eta_bl**14-34686.2970972733*eta_bl**13 &
               +108130.253843618*eta_bl**12-195823.099139525*eta_bl**11 &
               +227305.908339065*eta_bl**10-176106.001047617*eta_bl**9  &
               +92234.5885895112*eta_bl**8 -32700.3687158807*eta_bl**7  &
               +7923.51008739107*eta_bl**6 -1331.09245288739*eta_bl**5  &
               +130.109496961069*eta_bl**4 -7.64507811014497*eta_bl**3  &
               +6.94303207046209*eta_bl**2 -0.00209716712558639*eta_bl**1 ! &

          g_bl=g_bl+(1-exp(-delta_eta/delta_int))*eps_eta


          x_bl=1.0/(4.91**2*xnu)

          bxx1(j,k)=f_bl/1.0002014996204402/1.0000000359138641 !To assure 1.0 in infinity
          bxy1(j,k)=g_bl*sqrt(xnu/x_bl)/1.000546554
          bxz1(j,k)=0.0

       enddo
    enddo

    !STORE VALUE F_BL_INF G_BL_INF (ONLY ONE MORE TIME)------------------

    y=yly
    eta_bl=y*4.91/9.0  !The 9 is due to interpolation

    delta_eta=0.0
    eps_eta=0.0
    delta_int=0.2

    if (eta_bl .ge. (7.5/9.0)) then
       delta_eta=eta_bl-7.5/9.0
       eta_bl=7.5/9.0
       eps_eta=0.00015
    end if

    f_bl_inf=1678.64209592595*eta_bl**14-11089.6925017429*eta_bl**13 &
         +31996.4350140670*eta_bl**12-52671.5249779799*eta_bl**11 &
         +54176.1691167667*eta_bl**10-35842.8204706097*eta_bl**9  &
         +15201.3088871240*eta_bl**8 -4080.17137935648*eta_bl**7  &
         +702.129634528103*eta_bl**6 -56.2063925805318*eta_bl**5  &
         -17.0181128273914*eta_bl**4 +0.819582894357566*eta_bl**3  &
         -0.0601348202321954*eta_bl**2 +2.98973991270405*eta_bl**1


    f_bl_inf=f_bl_inf+(1-exp(-delta_eta/delta_int))*eps_eta
    f_bl_inf=f_bl_inf/1.0002014996204402/1.0000000359138641 !To assure 1.0 in infinity

    if (eta_bl .ge. (7.15/9.0)) then
       delta_int=0.8
       delta_eta=eta_bl-7.15/9.0
       eta_bl=7.15/9.0
       eps_eta=0.0005
    end if

    g_bl_inf=4924.05284779754*eta_bl**14-34686.2970972733*eta_bl**13 &
         +108130.253843618*eta_bl**12-195823.099139525*eta_bl**11 &
         +227305.908339065*eta_bl**10-176106.001047617*eta_bl**9  &
         +92234.5885895112*eta_bl**8 -32700.3687158807*eta_bl**7  &
         +7923.51008739107*eta_bl**6 -1331.09245288739*eta_bl**5  &
         +130.109496961069*eta_bl**4 -7.64507811014497*eta_bl**3  &
         +6.94303207046209*eta_bl**2 -0.00209716712558639*eta_bl**1


    g_bl_inf=g_bl_inf+(1-exp(-delta_eta/delta_int))*eps_eta
    g_bl_inf=g_bl_inf/1.000546554

    return
  end subroutine ecoule

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
  subroutine postprocess_tbl(ux1,uy1,uz1,pp3,phi1,ep1) !By Felipe Schuch

    USE MPI
    USE decomp_2d
    USE decomp_2d_io
    USE var, only : umean,vmean,wmean,pmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
    USE var, only : phimean, phiphimean
    USE var, only : ta1, pp1, di1
    USE var, only : ppi3, dip3
    USE var, only : pp2, ppi2, dip2

    USE var, ONLY : nxmsize, nymsize, nzmsize
    USE param, ONLY : npress

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(in) :: pp3
    character(len=30) :: filename

    integer :: is

    return
  end subroutine postprocess_tbl
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

end module tbl
