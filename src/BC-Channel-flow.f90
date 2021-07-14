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
module channel

  use decomp_2d
  use variables
  use param

  implicit none

  integer, save :: iochannel
  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_channel, boundary_conditions_channel, postprocess_channel, &
            visu_channel, momentum_forcing_channel, &
            geomcomplex_channel, finalize_channel

contains
  !############################################################################
  subroutine init_channel (ux1,uy1,uz1,ep1,phi1)

    use decomp_2d_io
    use MPI

    implicit none

    ! Arguments
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(out) :: ux1,uy1,uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar), intent(out) :: phi1

    ! Local variables
    logical :: read_from_file
    real(mytype) :: y, um
    integer :: k, j, i, is, code

    !
    ! Simulation can start from a 3D snapshot
    !
    if (nrank.eq.0) inquire(file="channel_init_ux", exist=read_from_file)
    call MPI_BCAST(read_from_file,1,MPI_LOGICAL,0,MPI_COMM_WORLD,code)
    if (code.ne.0) call decomp_2d_abort(code, "MPI_BCAST")
    if (read_from_file) then

       if (nrank.eq.0) print *, "Channel: init from snapshot."
       call decomp_2d_read_one(1,ux1,"channel_init_ux")
       call decomp_2d_read_one(1,uy1,"channel_init_uy")
       call decomp_2d_read_one(1,uz1,"channel_init_uz")
       if (iscalar.eq.1) then
          call decomp_2d_read_one(1,phi1(:,:,:,1),"channel_init_t")
          if (numscalar.gt.1) then
             phi1(:,:,:,2:numscalar) = zero
          endif
       endif

    else

       if (iscalar==1) then
          if (nrank.eq.0) print *,'Imposing linear temperature profile'
          do k=1,xsize(3)
             do j=1,xsize(2)
                if (istret.eq.0) y=real(j+xstart(2)-2,mytype)*dy
                if (istret.ne.0) y=yp(j+xstart(2)-1)
                do i=1,xsize(1)
                   phi1(i,j,k,:) = one - y/yly
                enddo
             enddo
          enddo

          if ((nclyS1.eq.2).and.(xstart(2).eq.1)) then
             !! Generate a hot patch on bottom boundary
             phi1(:,1,:,:) = one
          endif
          if ((nclySn.eq.2).and.(xend(2).eq.ny)) then
             phi1(:,xsize(2),:,:) = zero
          endif
       endif

       ux1=zero;uy1=zero;uz1=zero
       if (iin.ne.0) then
          call system_clock(count=code)
          if (iin.eq.2) code=0
          call random_seed(size = j)
          call random_seed(put = code+63946*nrank*(/ (i - 1, i = 1, j) /))

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
             enddo
          enddo
       enddo

    endif

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init end ok'
#endif

  end subroutine init_channel
  !############################################################################
  !############################################################################
  subroutine boundary_conditions_channel (ux,uy,uz,phi)

    use var, only : di2

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    if (.not.cpg) then ! if not constant pressure gradient
        call channel_cfr(ux,two/three)
    end if

    if (iscalar.ne.0) then
       if (iimplicit.le.0) then
          if ((nclyS1.eq.2).and.(xstart(2).eq.1)) then
             !! Generate a hot patch on bottom boundary
             phi(:,1,:,:) = one
          endif
          if ((nclySn.eq.2).and.(xend(2).eq.ny)) THEN
             phi(:,xsize(2),:,:) = zero
          endif
       else
          !
          ! Implicit boundary conditions are usually given in input file
          ! It is possible to modify g_sc here
          ! It is not possible to modify alpha_sc and beta_sc here
          !
          ! Bottom temperature if alpha_sc(:,1)=1 and beta_sc(:,1)=0 (default)
          !if (nclyS1.eq.2) g_sc(:,1) = one
          ! Top temperature if alpha_sc(:,2)=1 and beta_sc(:,2)=0 (default)
          !if (nclySn.eq.2) g_sc(:,2) = zero
       endif
    endif

  end subroutine boundary_conditions_channel
  !############################################################################
  !!
  !! Compute average of given array on the current CPU
  !!
  !! MPI call should follow to compute the average over the domain
  !!
  !############################################################################
  function channel_local_average(array)

    implicit none

    ! Argument
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: array
    ! Output
    real(mytype) :: channel_local_average
    ! Local variables
    integer :: i, j, k, jloc

    channel_local_average = zero
    do k = 1, xsize(3)
       do jloc = 1, xsize(2)
          j = jloc + xstart(2) - 1
          do i = 1, xsize(1)
            channel_local_average = channel_local_average + array(i,jloc,k) / ppy(j)
          enddo
       enddo
    enddo

  end function channel_local_average
  !############################################################################
  !!
  !!  SUBROUTINE: channel_cfr
  !!      AUTHOR: Kay Schäfer
  !! DESCRIPTION: Inforces constant flow rate without need of data transposition
  !!
  !############################################################################
  subroutine channel_cfr (ux, constant)

    use MPI

    implicit none

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux
    real(mytype), intent(in) :: constant

    integer :: code, i, j, k, jloc
    real(mytype) :: can, ub, coeff

    coeff = dy / (yly * real(xsize(1) * zsize(3), kind=mytype))

    ub = channel_local_average(ux) * coeff

    call MPI_ALLREDUCE(MPI_IN_PLACE,ub,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    if (code.ne.0) call decomp_2d_abort(code, "MPI_ALLREDUCE")

    can = - (constant - ub)

    if (nrank==0) print *, nrank, 'UT', ub, can

    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          ux(i,j,k) = ux(i,j,k) - can
        enddo
      enddo
    enddo

  end subroutine channel_cfr
  !############################################################################
  !!
  !! Monitor average and variance of the velocity and scalar
  !!
  !############################################################################
  subroutine postprocess_channel(ux1,uy1,uz1,pp3,phi1,ep1)

    use var, ONLY : di1, nzmsize
    use MPI

    implicit none

    ! Arguments
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3

    ! Local variables
    integer :: code, is, iref
    real(mytype) :: array(2*(3+numscalar)), coeff
    character(len=100) :: fileformat

    ! Init.
    array = zero
    coeff = dy / (yly * real(xsize(1) * zsize(3), kind=mytype))

    ! Average
    iref = 0
    array(iref+1) = channel_local_average(ux1)
    array(iref+2) = channel_local_average(uy1)
    array(iref+3) = channel_local_average(uz1)
    do is = 1, numscalar
      array(iref+3+is) = channel_local_average(phi1(:,:,:,is))
    enddo

    ! Variance
    iref = 3 + numscalar
    di1 = ux1**2; array(iref+1) = channel_local_average(di1)
    di1 = uy1**2; array(iref+2) = channel_local_average(di1)
    di1 = uz1**2; array(iref+3) = channel_local_average(di1)
    do is = 1, numscalar
      di1 = phi1(:,:,:,is)**2
      array(iref+3+is) = channel_local_average(di1)
    enddo

    ! Rescaling and MPI communication
    array = array * coeff
    if (nrank.eq.0) then
      call MPI_REDUCE(MPI_IN_PLACE, array, 2*(3+numscalar), real_type, MPI_SUM, 0, MPI_COMM_WORLD, code)
    else
      call MPI_REDUCE(array, array, 2*(3+numscalar), real_type, MPI_SUM, 0, MPI_COMM_WORLD, code)
    endif
    if (code.ne.0) call decomp_2d_abort(code, "MPI_REDUCE")

    ! Compute variance and log
    if (nrank.eq.0) then
      ! Compute variance
      do is = 1, 3 + numscalar
        array(3+numscalar+is) = array(3+numscalar+is) - array(is)**2
      enddo
      ! Print header at the first time step
      if (itime.eq.ifirst) then
        open(newunit=iochannel, file='channel.dat', form='formatted')
        write(iochannel,*) "# <u>, <v>, <w>, <T>, <u'²>, <v'²>, <w'²>, <T'²>"
      endif
      ! Use format to avoid printing all the digits
      write(fileformat, '( "(",I4,"(E18.10),A)" )' ) 2*(3+numscalar)
      write(iochannel, fileformat) array
    endif

  end subroutine postprocess_channel
  !############################################################################
  !!
  !!  SUBROUTINE: visu_channel
  !!      AUTHOR: FS
  !! DESCRIPTION: Performs channel-specific visualization
  !!
  !############################################################################
  subroutine visu_channel(ux1, uy1, uz1, pp3, phi1, ep1, num)

    use var, only : ux2, uy2, uz2, ux3, uy3, uz3
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    use var, ONLY : nzmsize
    use visu, only : write_field
    
    use ibm_param, only : ubcx,ubcy,ubcz

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    character(len=32), intent(in) :: num

    ! Write vorticity as an example of post processing

    ! Perform communications if needed
    if (sync_vel_needed) then
      call transpose_x_to_y(ux1,ux2)
      call transpose_x_to_y(uy1,uy2)
      call transpose_x_to_y(uz1,uz2)
      call transpose_y_to_z(ux2,ux3)
      call transpose_y_to_z(uy2,uy3)
      call transpose_y_to_z(uz2,uz3)
      sync_vel_needed = .false.
    endif

    !x-derivatives
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz)
    !y-derivatives
    call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
    call dery (tb2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
    call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
    !!z-derivatives
    call derz (ta3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
    call derz (tb3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcy)
    call derz (tc3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcz)
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

    !Q=-0.5*(ta1**2+te1**2+di1**2)-td1*tb1-tg1*tc1-th1*tf1
    di1 = zero
    di1(:,:,:) = - 0.5*(ta1(:,:,:)**2 + te1(:,:,:)**2 + ti1(:,:,:)**2) &
                 - td1(:,:,:) * tb1(:,:,:) &
                 - tg1(:,:,:) * tc1(:,:,:) &
                 - th1(:,:,:) * tf1(:,:,:)
    call write_field(di1, ".", "critq", trim(num))

  end subroutine visu_channel
  !############################################################################
  !############################################################################
  !!
  !!  SUBROUTINE: momentum_forcing
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Applies rotation for t < spinup_time.
  !!
  !############################################################################
  subroutine momentum_forcing_channel(dux1, duy1, ux1, uy1)

    implicit none

    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1
    real(mytype), intent(inout), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1

    if (cpg) then
        !! fcpg: add constant pressure gradient in streamwise direction
        dux1(:,:,:,1) = dux1(:,:,:,1) + fcpg !* (re/re_cent)**2
    endif

    if (irestart.eq.0 .and. itime.lt.spinup_time) then
       if (nrank==0) print *,'Rotating turbulent channel at speed ',wrotation
       dux1(:,:,:,1) = dux1(:,:,:,1) - wrotation*uy1(:,:,:)
       duy1(:,:,:,1) = duy1(:,:,:,1) + wrotation*ux1(:,:,:)
    endif

  end subroutine momentum_forcing_channel
  !############################################################################
  !############################################################################
  subroutine geomcomplex_channel(epsi,nxi,nxf,ny,nyi,nyf,nzi,nzf,yp,remp)

    use ibm

    implicit none

    integer                    :: nxi,nxf,ny,nyi,nyf,nzi,nzf
    real(mytype),dimension(nxi:nxf,nyi:nyf,nzi:nzf) :: epsi
    real(mytype),dimension(ny) :: yp
    real(mytype)               :: remp
    integer                    :: j
    real(mytype)               :: ym
    real(mytype)               :: zeromach
    real(mytype)               :: h

    epsi(:,:,:) = zero
    h = (yly - two) / two

    zeromach=one
    do while ((one + zeromach / two) .gt. one)
       zeromach = zeromach/two
    end do
    zeromach = 1.0e1*zeromach

    do j=nyi,nyf
       ym=yp(j)
       if ((ym.le.h).or.(ym.ge.(h+two))) then
          epsi(:,j,:)=remp
       endif
    enddo

    return
  end subroutine geomcomplex_channel
  !############################################################################
  !!
  !! Finalize the channel case
  !!
  !############################################################################
  subroutine finalize_channel()

    implicit none

    if (nrank.eq.0) then
      close(iochannel)
    endif

  end subroutine finalize_channel

end module channel
