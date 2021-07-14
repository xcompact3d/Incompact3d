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

  ! Experimental
  ! .false. requires nstat=1 and xenS(1)=xstS(1)
  logical, parameter :: flag_3D_IO = .true.

  private
  public overall_statistic

contains

  !
  ! Initialize to zero all statistics
  !
  subroutine init_statistic

    use param, only : zero, iscalar
    use var, only : tmean
    use var, only : pmean, ppmean
    use var, only : umean, uumean
    use var, only : vmean, vvmean
    use var, only : wmean, wwmean
    use var, only : uvmean, uwmean
    use var, only : vwmean
    use var, only : phimean, phiphimean

    implicit none

    tmean = zero
    pmean = zero
    umean = zero
    uumean = zero
    vmean = zero
    vvmean = zero
    wmean = zero
    wwmean = zero
    uvmean = zero
    uwmean = zero
    vwmean = zero
    ppmean = zero
    if (iscalar==1) then
      phimean = zero
      phiphimean = zero
    endif

  end subroutine init_statistic

  !
  ! Read all statistics if possible
  !
  subroutine restart_statistic

    use param, only : initstat, irestart, ifirst, zero
    use variables, only : nstat
    use var, only : tmean

    implicit none

    ! No reading for statistics when nstat > 1 or no restart
    if (nstat.gt.1 .or. irestart.eq.0) then
       call init_statistic()
       initstat = ifirst
       return
    endif

    ! Temporary array
    tmean = zero

    ! Read all statistics
    call read_or_write_all_stats(.true.)

  end subroutine restart_statistic

  !
  ! Statistics: perform all IO
  !
  subroutine read_or_write_all_stats(flag_read)

    use param, only : iscalar, itime
    use variables, only : numscalar
    use decomp_2d, only : nrank
    use var, only : pmean, ppmean
    use var, only : umean, uumean
    use var, only : vmean, vvmean
    use var, only : wmean, wwmean
    use var, only : uvmean, uwmean
    use var, only : vwmean
    use var, only : phimean, phiphimean

    implicit none

    ! Argument
    logical, intent(in) :: flag_read

    ! Local variables
    integer :: is
    character(len=30) :: filename

    if (nrank==0) then
      print *,'==========================================================='
      if (flag_read) then
        print *,'Reading stat file', itime
      else
        print *,'Writing stat file', itime
      endif
    endif

    write(filename,"('pmean.dat',I7.7)") itime
    call read_or_write_one_stat(flag_read, filename, pmean)
    write(filename,"('umean.dat',I7.7)") itime
    call read_or_write_one_stat(flag_read, filename, umean)
    write(filename,"('vmean.dat',I7.7)") itime
    call read_or_write_one_stat(flag_read, filename, vmean)
    write(filename,"('wmean.dat',I7.7)") itime
    call read_or_write_one_stat(flag_read, filename, wmean)

    write(filename,"('uumean.dat',I7.7)") itime
    call read_or_write_one_stat(flag_read, filename, uumean)
    write(filename,"('vvmean.dat',I7.7)") itime
    call read_or_write_one_stat(flag_read, filename, vvmean)
    write(filename,"('wwmean.dat',I7.7)") itime
    call read_or_write_one_stat(flag_read, filename, wwmean)
    write(filename,"('ppmean.dat',I7.7)") itime
    call read_or_write_one_stat(flag_read, filename, ppmean)

    write(filename,"('uvmean.dat',I7.7)") itime
    call read_or_write_one_stat(flag_read, filename, uvmean)
    write(filename,"('uwmean.dat',I7.7)") itime
    call read_or_write_one_stat(flag_read, filename, uwmean)
    write(filename,"('vwmean.dat',I7.7)") itime
    call read_or_write_one_stat(flag_read, filename, vwmean)

    if (iscalar==1) then
       do is=1, numscalar
          write(filename,"('phi',I2.2,'mean.dat',I7.7)") is, itime
          call read_or_write_one_stat(flag_read, filename, phimean(:,:,:,is))
          write(filename,"('phiphi',I2.2,'mean.dat',I7.7)") is, itime
          call read_or_write_one_stat(flag_read, filename, phiphimean(:,:,:,is))
       enddo
    endif

    if (nrank==0) then
      if (flag_read) then
        print *,'Read stat done!'
      else
        print *,'Write stat done!'
      endif
      print *,'==========================================================='
    endif

  end subroutine read_or_write_all_stats

  !
  ! Statistics: perform one IO
  !
  subroutine read_or_write_one_stat(flag_read, filename, array)

    use decomp_2d, only : mytype, xstS, xenS
    use decomp_2d_io, only : decomp_2d_read_one, decomp_2d_write_one

    implicit none

    ! Arguments
    logical, intent(in) :: flag_read
    character(len=*), intent(in) :: filename
    real(mytype), dimension(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)), intent(inout) :: array

    if (flag_read) then
      if (flag_3D_IO) then
        ! There was a check for nvisu = 1 before
        call decomp_2d_read_one(1, array, filename)
      else
        call decomp_2d_read_plane(1, array, filename)
      endif
    else
      if (flag_3D_IO) then
        call decomp_2d_write_one(1, array, filename, 1)
      else
        call decomp_2d_write_plane(1, array, 1, 1, filename)
      endif
    endif

  end subroutine read_or_write_one_stat

  !
  ! Statistics : Intialize, update and perform IO
  !
  subroutine overall_statistic(ux1,uy1,uz1,phi1,pp3,ep1)

    use param
    use variables
    use decomp_2d
    use decomp_2d_io
    use tools, only : rescale_pressure

    use var, only : nxmsize, nymsize, nzmsize
    use var, only : ppi3, dip3
    use var, only : pp2, ppi2, dip2
    use var, only : pp1, ta1, di1

    use var, only : pmean, ppmean
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
    elseif (itime.eq.initstat) then
       call init_statistic()
    elseif (itime.eq.ifirst) then
        call restart_statistic()
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
    ! Convert to physical pressure
    call rescale_pressure(ta1)
    call update_average_scalar(pmean, ta1, ep1)
    call update_average_scalar(ppmean, ta1**2, ep1)

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
          call update_average_scalar(phiphimean(:,:,:,is), phi1(:,:,:,is)**2, ep1)
       enddo
    endif

    ! Write all statistics
    if (mod(itime,icheckpoint)==0) then
       call read_or_write_all_stats(.false.)
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
    use param, only : itime, initstat
    use var, only : di1, tmean

    implicit none

    ! inputs
    real(mytype), dimension(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)), intent(inout) :: um
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux, ep

    di1 = one_minus_ep1(ux, ep)
    if (flag_3D_IO) then
      call fine_to_coarseS(1, di1, tmean)
    else
      tmean(xstS(1),:,:) = sum(di1(:,:,:), dim=1)/real(xsize(1),kind=mytype)
    endif
    um = um + (tmean - um) / real(itime-initstat+1, kind=mytype)

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
