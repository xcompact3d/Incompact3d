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

module stats

  implicit none

  private
  public overall_statistic

contains

  subroutine overall_statistic(ux1,uy1,uz1,phi1,pp3,ep1)

    use param
    use variables
    use decomp_2d
    use decomp_2d_io

    use var, only : nxmsize, nymsize, nzmsize
    use var, only : ppi3, dip3
    use var, only : pp2, ppi2, dip2
    use var, only : pp1, ta1, di1

    use var, only : pmean
    use var, only : umean, uumean
    use var, only : vmean, vvmean
    use var, only : wmean, wwmean
    use var, only : uvmean, uwmean
    use var, only : vwmean
    use var, only : phimean, phiphimean

    implicit none

    !! Inputs
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar),intent(in) :: phi1
    real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3

    !! Locals
    integer :: is
    character(len=30) :: filename

    if (itime.lt.initstat) then
       return
    endif

    !! Mean pressure
    !WORK Z-PENCILS
    call interzpv(ppi3,pp3(:,:,:,1),dip3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
         (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
    !WORK Y-PENCILS
    call transpose_z_to_y(ppi3,pp2,ph3) !nxm nym nz
    call interypv(ppi2,pp2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
         (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
    !WORK X-PENCILS
    call transpose_y_to_x(ppi2,pp1,ph2) !nxm ny nz
    call interxpv(ta1,pp1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)
    call update_average_scalar(pmean, ta1, ep1)

    !! Mean velocity
    call update_average_vector(umean, vmean, wmean, &
                               ux1, uy1, uz1, ep1)

    !! Second-order velocity moments
    call update_variance_vector(uumean, vvmean, wwmean, uvmean, uwmean, vwmean, &
                                ux1, uy1, uz1, ep1)

    !! Scalar statistics
    if (iscalar==1) then
       do is=1, numscalar
          !pmean=phi1
          call update_average_scalar(phimean(:,:,:,is), phi1(:,:,:,is), ep1)

          !phiphimean=phi1*phi1
          call update_average_scalar(phiphimean(:,:,:,is), phi1(:,:,:,is)*phi1(:,:,:,is), ep1)
       enddo
    endif

    if (mod(itime,icheckpoint)==0) then

       if (nrank==0) then
          print *,'===========================================================<<<<<'
          print *,'Writing stat file',itime
       endif

       write(filename,"('pmean.dat',I7.7)") itime
       call decomp_2d_write_one(1,pmean,filename,1)

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

       write(filename,"('kmean.dat',I7.7)") itime
       call decomp_2d_write_one(1,half*(uumean+vvmean+wwmean),filename,1)

       if (iscalar==1) then
          do is=1, numscalar
             write(filename,"('phi',I2.2,'mean.dat',I7.7)") is, itime
             call decomp_2d_write_one(1,phimean(:,:,:,is),filename,1)
             write(filename,"('phiphi',I2.2,'mean.dat',I7.7)") is, itime
             call decomp_2d_write_one(1,phiphimean(:,:,:,is),filename,1)
          enddo
       endif

       if (nrank==0) then
          print *,'write stat done!'
          print *,'===========================================================<<<<<'
       endif
    endif

  end subroutine overall_statistic

  !
  ! Basic function, can be applied to arrays
  !
  elemental real(mytype) function one_minus_ep1(var, ep1)

    use decomp_2d, only : mytype
    use param, only : iibm, one

    implicit none

    ! inputs
    real(mytype), intent(in) :: var, ep1

    if (iibm==2) then
      one_minus_ep1 = (one - ep1) * var
    else
      one_minus_ep1 = var
    endif

  end function one_minus_ep1

  !
  ! Update um, the average of ux
  !
  subroutine update_average_scalar(um, ux, ep)

    use decomp_2d, only : mytype, xsize, xstS, xenS, fine_to_coarseS
    use var, only : di1, tmean

    implicit none

    ! inputs
    real(mytype), dimension(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)), intent(inout) :: um
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux, ep

    di1 = one_minus_ep1(ux, ep)
    call fine_to_coarseS(1, di1, tmean)
    um = um + tmean

  end subroutine update_average_scalar

  !
  ! Update (um, vm, wm), the average of (ux, uy, uz)
  !
  subroutine update_average_vector(um, vm, wm, ux, uy, uz, ep)

    use decomp_2d, only : mytype, xsize, xstS, xenS

    implicit none

    ! inputs
    real(mytype), dimension(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)), intent(inout) :: um, vm, wm
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux, uy, uz, ep

    call update_average_scalar(um, ux, ep)
    call update_average_scalar(vm, uy, ep)
    call update_average_scalar(wm, uz, ep)

  end subroutine update_average_vector

  !
  ! Update (uum, vvm, wwm, uvm, uwm, vwm), the variance of the vector (ux, uy, uz)
  !
  subroutine update_variance_vector(uum, vvm, wwm, uvm, uwm, vwm, ux, uy, uz, ep)

    use decomp_2d, only : mytype, xsize, xstS, xenS

    implicit none

    ! inputs
    real(mytype), dimension(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)), intent(inout) :: uum, vvm, wwm, uvm, uwm, vwm
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux, uy, uz, ep

    call update_average_scalar(uum, ux*ux, ep)
    call update_average_scalar(vvm, uy*uy, ep)
    call update_average_scalar(wwm, uz*uz, ep)
    call update_average_scalar(uvm, ux*uy, ep)
    call update_average_scalar(uwm, ux*uz, ep)
    call update_average_scalar(vwm, uy*uz, ep)

  end subroutine update_variance_vector

endmodule stats

! !############################################################################
! subroutine CONVERGENCE_STATISTIC(ux1,ep1,u1sum_tik,u1sum_tak,tsum)

!   USE param
!   USE variables
!   USE decomp_2d
!   USE decomp_2d_io

!   implicit none

!   real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,ep1
!   real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: u1sum_tik,u1sum_tak,tsum
!   character(len=30) :: filename

!   if (iibm.ne.0) then
!      call fine_to_coarseS(1,(1-ep1)*ux1,tsum)
!   else
!      call fine_to_coarseS(1,ux1,tsum)
!   endif

!   u1sum_tik=u1sum_tik+tsum

!   if (mod(itime,ntik)==0) then
!      if (nrank.eq.0) print *,'Saving uxmean tik =>',(itime/(ntik/2)-1)
!      write(filename,"('uxmean',I4.4)") (itime/(ntik/2)-1)
!      call decomp_2d_write_one(1,u1sum_tik/ntik,filename,1)
!      u1sum_tik = zero
!   endif

!   if (itime.ge.(ntik/2)) then

!      u1sum_tak=u1sum_tak+tsum

!      if (itime.gt.(ntik/2).AND.mod(itime+ntik/2,ntik)==0) then
!         if (nrank.eq.0) print *,'Saving uxmean tak =>',(itime/(ntik/2)-1)
!         write(filename,"('uxmean',I4.4)") (itime/(ntik/2) -1)
!         call decomp_2d_write_one(1,u1sum_tak/ntik,filename,1)
!         u1sum_tak = zero
!      endif

!   endif

! end subroutine CONVERGENCE_STATISTIC
! !############################################################################
! subroutine CONVERGENCE_STATISTIC2(ux1,ep1,tik1,tik2,tak1,tak2)

!   USE param
!   USE variables
!   USE decomp_2d
!   USE decomp_2d_io

!   implicit none

!   real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,ep1
!   real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: tik1,tik2,tak1,tak2,tsum
!   real(mytype) :: rms1
!   character(len=30) :: filename

!   rms1 = zero

!   if (iibm.ne.0) then
!      call fine_to_coarseS(1,(1-ep1)*ux1,tsum)
!   else
!      call fine_to_coarseS(1,ux1,tsum)
!   endif

!   tik1 = tik1 + tsum
!   tak1 = tak1 + tsum

!   if (mod(itime,ntik)==0) then
!      tik2 = tik1/ntik
!      tik1 = zero
!   end if

!   if (itime.eq.(ntik/2.)) tak1 = zero

!   if ((itime.gt.(ntik/2.)).AND.(mod(itime+ntik/2,ntik)==0)) then
!      tak2 = tak1/ntik
!      tak1 = zero
!   end if

!   if (itime.ge.(int((3./2.)*ntik))) then

!      if ((mod(itime,ntik)==0).OR.(mod(itime+ntik/2,ntik)==0)) then

!         call RMS(tik2,tak2,rms1)

!         if (nrank .eq. 0) then

!            print *,'RMS=',rms1

!            write(filename,"('rms',I8.8)") itime
!            open(67,file=trim(filename),status='unknown',form='formatted')
!            write(67,"(2E14.6,I14)") t,rms1,itime
!            close(67)

!            rms1=zero

!         end if
!      endif
!   end if

! end subroutine CONVERGENCE_STATISTIC2
! !############################################################################
! subroutine RMS(meanA,meanB,rms1)

!   USE param
!   USE variables
!   USE decomp_2d
!   USE MPI

!   implicit none

!   real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: meanA, meanB
!   real(mytype) :: rms0,rms1
!   character(len=30) :: filename
!   integer :: i,j,k,code

!   rms1 = zero
!   do k = 1, xszS(3)
!      do j = 1, xszS(2)
!         do i = 1, xszS(1)
!            rms1=rms1+(meanB(i,j,k)-meanA(i,j,k))**2
!         enddo
!      enddo
!   enddo

!   rms0 = zero
!   rms0 = sqrt(rms1/(nx*ny*nz))
!   rms1 = zero

!   call MPI_REDUCE(rms0,rms1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
!   return

! end subroutine RMS
! !############################################################################
! subroutine EXTRA_STAT (ux1,uy1,uz1,uvisu,tsum,dudxsum,utmapsum)

!   USE param
!   USE variables
!   USE decomp_2d
!   USE decomp_2d_io

!   implicit none

!   real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
!   real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu,tsum,dudxsum,utmapsum

!   real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,td1,tf1,di1 !ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,phim1,temp1
!   real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tc2,td2,tf2 !ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,phim2,di2,temp2
!   !real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: !ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phim3,temp3



!   !x-derivatives
!   call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
!   !call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
!   !call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
!   !y-derivatives
!   call transpose_x_to_y(ux1,td2)
!   !call transpose_x_to_y(uy1,te2)
!   call transpose_x_to_y(uz1,tf2)
!   call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
!   !call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
!   call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
!   !!z-derivatives
!   !call transpose_y_to_z(td2,td3)
!   !call transpose_y_to_z(te2,te3)
!   !call transpose_y_to_z(tf2,tf3)
!   !call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
!   !call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
!   !call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
!   !!all back to x-pencils
!   !call transpose_z_to_y(ta3,td2)
!   !call transpose_z_to_y(tb3,te2)
!   !call transpose_z_to_y(tc3,tf2)
!   !call transpose_y_to_x(td2,tg1)
!   !call transpose_y_to_x(te2,th1)
!   !call transpose_y_to_x(tf2,ti1)
!   call transpose_y_to_x(ta2,td1)
!   !call transpose_y_to_x(tb2,te1)
!   call transpose_y_to_x(tc2,tf1)
!   !du/dx=ta1 du/dy=td1 and du/dz=tg1
!   !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
!   !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1


!   !FIRST DERIVATIVE
!   call fine_to_coarseS(1,ta1,tsum)
!   dudxsum = dudxsum + tsum
!   if (mod(itime,imodulo)==0) then
!      call decomp_2d_write_one(1,dudx,'dudx.sum',1)
!   endif

!   !FRICTION VELOCITY OVER ITERATION
!   di1 = zero
!   di1 = sqrt(sqrt( td1**2 + tf1**2 )*xnu) !du/dy**2 + dw/dy**2

!   call fine_to_coarseS(1,ta1,tsum)
!   utmapsum=utmapsum+tsum
!   if (mod(itime,imodulo)==0) then
!      call decomp_2d_write_plane(1,utmapsum,2,1,'utmap.sum')
!   endif



! end subroutine EXTRA_STAT
