!################################################################################
!This file is part of Xcompact3d.
!
!Xcompact3d
!Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!    Xcompact3d is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation.
!
!    Xcompact3d is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the code.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    We kindly request that you cite Xcompact3d/Incompact3d in your
!    publications and presentations. The following citations are suggested:
!
!    1-Laizet S. & Lamballais E., 2009, High-order compact schemes for
!    incompressible flows: a simple and efficient method with the quasi-spectral
!    accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
!
!    2-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence
!    problems with up to 0(10^5) computational cores, Int. J. of Numerical
!    Methods in Fluids, vol 67 (11), pp 1735-1757
!################################################################################

module tbl

  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_tbl, boundary_conditions_tbl, postprocess_tbl, tbl_flrt

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

       phi1(:,:,:,:) = 0.25 !change as much as you want
          if ((nclyS1.eq.2).and.(xstart(2).eq.1)) then
             !! Generate a hot patch on bottom boundary
             phi1(:,1,:,:) = one
          endif
          if ((nclySn.eq.2).and.(xend(2).eq.ny)) THEN
             phi1(:,xsize(2),:,:) = 0.25
          endif

    endif
    ux1=zero;uy1=zero;uz1=zero

    !a blasius profile is created in ecoule and then duplicated for the all domain
    call blasius()

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

    !INFLOW with an update of bxx1, byy1 and bzz1 at the inlet

    call blasius()
    !INLET FOR SCALAR, TO BE CONSISTENT WITH INITIAL CONDITION
    if (iscalar.eq.1) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             phi(1,:,:,:)=0.25
             if ((xstart(2).eq.1)) then
                phi(:,1,:,:) = one
             endif
             if ((xend(2).eq.ny)) THEN
                phi(:,xsize(2),:,:) = 0.25
             endif
          enddo
       enddo
    endif

    !OUTFLOW based on a 1D convection equation

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
          if (iscalar.eq.1) phi(nx,:,:,:) =  phi(nx,:,:,:) - cx*(phi(nx,:,:,:)-phi(nx-1,:,:,:))
          enddo
    enddo

    !! Bottom Boundary
    if (ncly1 == 2) THEN
      do k = 1, xsize(3)
        do i = 1, xsize(1)
          byx1(i, k) = zero
          byy1(i, k) = zero
          byz1(i, k) = zero
        enddo
      enddo
    endif
    !! Top Boundary
    if (nclyn == 2) then
       do k = 1, xsize(3)
          do i = 1, xsize(1)
             byxn(i, k) = ux(i, xsize(2) - 1, k)
             byyn(i, k) = uy(i, xsize(2) - 1, k)
             byzn(i, k) = uz(i, xsize(2) - 1, k)
          enddo
       enddo
    endif

    !SCALAR   
    if (itimescheme.ne.7) then
    if (iscalar.ne.0) then
          if ((nclyS1.eq.2).and.(xstart(2).eq.1)) then
             !! Generate a hot patch on bottom boundary
             phi(1,1,:,:) = one
          endif
          if ((nclySn.eq.2).and.(xend(2).eq.ny)) THEN
             phi(1,xsize(2),:,:) = phi(1,xsize(2)-1,:,:)
          endif
    endif
    endif

    !update of the flow rate (what is coming in the domain is getting out)
    call tbl_flrt(ux,uy,uz)

    return
  end subroutine boundary_conditions_tbl

  !********************************************************************

  !********************************************************************
!
subroutine tbl_flrt (ux1,uy1,uz1)
  !
  !********************************************************************

  USE decomp_2d
  USE decomp_2d_poisson
  USE variables
  USE param
  USE MPI

  implicit none
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2

  integer :: j,i,k,code
  real(mytype) :: can,ut1,ut2,ut3,ut4,utt1,utt2,utt3,utt4,udif

  ux1(1,:,:)=bxx1(:,:)
  ux1(nx,:,:)=bxxn(:,:)

  call transpose_x_to_y(ux1,ux2)
  call transpose_x_to_y(uy1,uy2)
  ! Flow rate at the inlet
  ut1=zero;utt1=zero
  if (ystart(1)==1) then !! CPUs at the inlet
    do k=1,ysize(3)
      do j=1,ysize(2)-1
        ut1=ut1+(yp(j+1)-yp(j))*(ux2(1,j+1,k)-half*(ux2(1,j+1,k)-ux2(1,j,k)))
      enddo
    enddo
    ! ut1=ut1/real(ysize(3),mytype)
  endif
  call MPI_ALLREDUCE(ut1,utt1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
  utt1=utt1/real(nz,mytype) !! Volume flow rate per unit spanwise dist
  ! Flow rate at the outlet
  ut2=zero;utt2=zero
  if (yend(1)==nx) then !! CPUs at the outlet
    do k=1,ysize(3)
      do j=1,ysize(2)-1
        ut2=ut2+(yp(j+1)-yp(j))*(ux2(ysize(1),j+1,k)-half*(ux2(ysize(1),j+1,k)-ux2(ysize(1),j,k)))
      enddo
    enddo
    ! ut2=ut2/real(ysize(3),mytype)
  endif

  call MPI_ALLREDUCE(ut2,utt2,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
  utt2=utt2/real(nz,mytype) !! Volume flow rate per unit spanwise dist

  ! Flow rate at the top and bottom
  ut3=zero
  ut4=zero
  do k=1,ysize(3)
    do i=1,ysize(1)
      ut3=ut3+uy2(i,1,k)
      ut4=ut4+uy2(i,ny,k)
    enddo
  enddo
  call MPI_ALLREDUCE(ut3,utt3,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
  call MPI_ALLREDUCE(ut4,utt4,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
  utt3=utt3/(real(nx*nz,mytype))*xlx  !!! Volume flow rate per unit spanwise dist
  utt4=utt4/(real(nx*nz,mytype))*xlx  !!! Volume flow rate per unit spanwise dist

  !! velocity correction
  udif=(utt1-utt2+utt3-utt4)/yly
  if (nrank==0 .and. mod(itime,1)==0) then
    write(*,"(' Mass balance: L-BC, R-BC,',2f12.6)") utt1,utt2
    write(*,"(' Mass balance: B-BC, T-BC, Crr-Vel',3f11.5)") utt3,utt4,udif
  endif
  ! do k=1,xsize(3)
  !   do j=1,xsize(2)
  !     ux1(nx,i,k)=ux1(nx,i,k)+udif
  !   enddo
  ! enddo
  do k=1,xsize(3)
    do j=1,xsize(2)
      bxxn(j,k)=bxxn(j,k)+udif
    enddo
  enddo


  return
end subroutine tbl_flrt

!********************************************************************
  subroutine blasius()


    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none

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
  end subroutine blasius

  !############################################################################
  subroutine postprocess_tbl(ux1,uy1,uz1,ep1) !By Felipe Schuch

    USE MPI
    USE decomp_2d
    USE decomp_2d_io
    USE var, only : umean,vmean,wmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
    USE var, only : uvisu
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    character(len=30) :: filename

       !! Write vorticity as an example of post processing
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

    di1(:,:,:)=sqrt((tf1(:,:,:)-th1(:,:,:))**2+(tg1(:,:,:)-tc1(:,:,:))**2+&
         (tb1(:,:,:)-td1(:,:,:))**2)
    if (iibm==2) then
       di1(:,:,:) = (one - ep1(:,:,:)) * di1(:,:,:)
    endif
    uvisu=0.
    call fine_to_coarseV(1,di1,uvisu)
994 format('vort',I3.3)
    write(filename, 994) itime/ioutput
    call decomp_2d_write_one(1,uvisu,filename,2)

    return
  end subroutine postprocess_tbl

end module tbl
