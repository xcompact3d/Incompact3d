!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module les

  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
  use visu, only :  output2D
  use utilities, only : gen_filename

  character(len=*), parameter :: io_turb = "turb-io", &
       turb_dir = "turb-data"
contains

  subroutine init_explicit_les
    !================================================================================
    !
    !  SUBROUTINE: init_explicit_les
    ! DESCRIPTION: Initialises the explicit LES parameters
    !      AUTHOR: G. Deskos <g.deskos14@imperial.ac.uk>
    !
    !================================================================================

    USE param
    USE variables
    use decomp_2d_io, only : decomp_2d_init_io, decomp_2d_register_variable, decomp_2d_open_io, decomp_2d_write_mode

    implicit none

    integer :: j

    !you can specify other metric e.g. max(dx,dy,dz) or min(dx,dy,dz)

    do j = 1, ny - 1
       del(j) = (dx * (yp(j + 1) - yp(j)) * dz)**(one / three)
    enddo
    del(ny) = del(ny - 1)

    if(nrank==0) then

       write(*, *) ' '
       write(*, *) '++++++++++++++++++++++++++++++++'
       write(*, *) 'LES Modelling'
       if(jles==1) then
          write(*, *) ' Classic Smagorinsky is used ... '
          write(*, *) ' Smagorinsky constant = ', smagcst
          if (itype==itype_abl) then
             write(*,*) ' with Mason and Thomson damping function ... '
             write(*,*) ' n = ', nSmag
          endif
       else if (jles==2) then
          write(*, *) ' WALE SGS model is used ... '
          write(*, *) ' Max value for the WALE constant  = ', walecst
       else if (jles==3) then
          write(*, *) ' Dynamic Smagorinsky is used ... '
          write(*, *) ' Max value for the dynamic constant field = ', maxdsmagcst
       endif
       write(*, *) '++++++++++++++++++++++++++++++++'
       !if (nrank==0) print *, "Del y min max= ", minval(del), maxval(del)
       write(*, *) '++++++++++++++++++++++++++++++++'
       write(*, *) ' '
    endif

    !! Initialise visualisation output
    call decomp_2d_init_io(io_turb)

    if (jles .eq. 1) then ! Smagorinsky
       call decomp_2d_register_variable(io_turb, "nut_smag", 1, 0, output2D, mytype)
    else if (jles .eq. 2) then ! WALE
       call decomp_2d_register_variable(io_turb, "nut_wale", 1, 0, output2D, mytype)
    else if (jles .eq. 3) then ! Lilly-style Dynamic Smagorinsky
       call decomp_2d_register_variable(io_turb, "nut_smag", 1, 0, output2D, mytype) ! Written by call to smag in dynsmag
       call decomp_2d_register_variable(io_turb, "dsmagcst_final", 1, 0, output2D, mytype)
       call decomp_2d_register_variable(io_turb, "nut_dynsmag", 1, 0, output2D, mytype)
    end if

#ifdef ADIOS2
    call decomp_2d_open_io(io_turb, turb_dir, decomp_2d_write_mode)
#endif

  end subroutine init_explicit_les
  subroutine finalise_explicit_les()

    use decomp_2d_io, only : decomp_2d_close_io
    implicit none
    
#ifdef ADIOS2
    call decomp_2d_close_io(io_turb, turb_dir)
#endif

  end subroutine finalise_explicit_les
  
  !************************************************************
  subroutine compute_SGS(sgsx1,sgsy1,sgsz1,ux1,uy1,uz1,phi1,ep1)
    !================================================================================
    !
    !  SUBROUTINE: compute_SGS
    ! DESCRIPTION: computes the SGS terms (divergence of the SGS stresses) used in the
    !              momentum equation
    !      AUTHOR: G. Deskos <g.deskos14@imperial.ac.uk>
    !
    !================================================================================

    USE param
    USE variables
    USE decomp_2d_io
    use var, only: nut1
    use var, only: ta1, tb1, tc1
    use abl, only: wall_sgs_slip, wall_sgs_noslip
    implicit none

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ux1, uy1, uz1, ep1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar), intent(in) :: phi1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(out) :: sgsx1, sgsy1, sgsz1
    !real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: w_fluxx1, w_fluxy1, w_fluxz1
    !real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: w_sgsx1, w_sgsy1, w_gsz1

    ! Calculate eddy-viscosity
    if(jles.eq.1) then ! Smagorinsky
       call smag(nut1,ux1,uy1,uz1)

    elseif(jles.eq.2) then ! Wall-adapting local eddy-viscosity (WALE) model
       call wale(nut1,ux1,uy1,uz1)

    elseif(jles.eq.3) then ! Lilly-style Dynamic Smagorinsky
       call dynsmag(nut1,ux1,uy1,uz1,ep1)

    endif

    ! we can store wallflux in ta1, tb1, tc1
    ! this will be done in all call to wall_sgs_ for the ABL
    if(iconserv.eq.0) then ! Non-conservative form for calculating the divergence of the SGS stresses
       call sgs_mom_nonconservative(sgsx1,sgsy1,sgsz1,ux1,uy1,uz1,nut1,ep1)

       ! SGS correction for ABL
       if(itype.eq.itype_abl) then
          ! No-slip wall
          if (ncly1==2) then 
             call wall_sgs_noslip(ux1,uy1,uz1,nut1,ta1,tb1,tc1)
             if(xstart(2)==1) then
               sgsx1(:,2,:) = -ta1(:,2,:)
               sgsy1(:,2,:) = -tb1(:,2,:)
               sgsz1(:,2,:) = -tc1(:,2,:)
             endif
          ! Slip wall
          elseif (ncly1==1) then
             call wall_sgs_slip(ux1,uy1,uz1,phi1,nut1,ta1,tb1,tc1)
             if(xstart(2)==1) then
                sgsx1(:,1,:) = ta1(:,1,:)
                sgsy1(:,1,:) = tb1(:,1,:)
                sgsz1(:,1,:) = tc1(:,1,:)
             endif
          endif
       endif

    elseif (iconserv.eq.1) then ! Conservative form for calculating the divergence of the SGS stresses (used with wall functions)
       call sgs_mom_conservative(sgsx1,sgsy1,sgsz1,ux1,uy1,uz1,nut1)

    endif

    return

  end subroutine compute_SGS


  subroutine smag(nut1,ux1,uy1,uz1)
    !================================================================================
    !
    !  SUBROUTINE: smag
    ! DESCRIPTION: Calculates the eddy-viscosity nut according to the standard
    !              Smagorinsky model
    !      AUTHOR: G. Deskos <g.deskos14@imperial.ac.uk>
    !
    !================================================================================

    use MPI
    USE param
    USE variables
    USE decomp_2d_io
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,di2
    USE var, only : ux3,uy3,uz3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
    USE var, only : sxx1,syy1,szz1,sxy1,sxz1,syz1,srt_smag
    USE var, only : gxx1,gyx1,gzx1,gxy2,gyy2,gzy2,gxz3,gyz3,gzz3
    USE var, only : gxy1,gyy1,gzy1,gxz2,gyz2,gzz2,gxz1,gyz1,gzz1
    USE var, only : sxx2,syy2,szz2,sxy2,sxz2,syz2,srt_smag2,nut2
    USE var, only : sxx3,syy3,szz3,sxy3,sxz3,syz3
    USE ibm_param

    implicit none

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: nut1
    real(mytype) :: smag_constant, y, length
    real(mytype) :: nutmin, nutmax
    real(mytype) :: srtmin, srtmax

    integer :: i, j, k, ierr
    character(len = 30) :: filename


    ! INFO about the auxillary arrays
    !--------------------------------------------------------
    ! gxx= dux/dx; gyx=duy/dx; gzx=duz/dx;
    ! gxy= dux/dy; gyy=duy/dy; gzy=duz/dy;
    ! gxz= dux/dz; gyz=duy/dz; gzz=duz/dz

    call derx (gxx1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)
    call derx (gyx1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
    call derx (gzx1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz)

    sxx1(:,:,:) = gxx1(:,:,:)

    !WORK Y-PENCILS
    call transpose_x_to_y(ux1,ux2)
    call transpose_x_to_y(uy1,uy2)
    call transpose_x_to_y(uz1,uz2)
    call transpose_x_to_y(gyx1,ta2)

    call dery (gxy2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
    call dery (gyy2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
    call dery (gzy2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)

    sxy2(:,:,:)=half*(gxy2(:,:,:)+ta2(:,:,:))
    syy2(:,:,:)=gyy2(:,:,:)

    !WORK Z-PENCILS
    call transpose_y_to_z(ux2,ux3)
    call transpose_y_to_z(uy2,uy3)
    call transpose_y_to_z(uz2,uz3)
    call transpose_y_to_z(gzy2,ta3)

    call derz(gxz3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
    call derz(gyz3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcz)
    call derz(gzz3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcz)

    szz3(:,:,:)=gzz3(:,:,:)
    syz3(:,:,:)=half*(gyz3(:,:,:)+ta3(:,:,:))

    !WORK Y-PENCILS
    call transpose_z_to_y(syz3,syz2)
    call transpose_z_to_y(szz3,szz2)

    call transpose_z_to_y(gxz3,gxz2)
    call transpose_z_to_y(gyz3,gyz2)
    call transpose_z_to_y(gzz3,gzz2)

    !WORK X-PENCILS
    call transpose_y_to_x(sxy2,sxy1)
    call transpose_y_to_x(syy2,syy1)
    call transpose_y_to_x(syz2,syz1)
    call transpose_y_to_x(szz2,szz1)

    call transpose_y_to_x(gxy2,gxy1)
    call transpose_y_to_x(gyy2,gyy1)
    call transpose_y_to_x(gzy2,gzy1)
    call transpose_y_to_x(gxz2,gxz1)
    call transpose_y_to_x(gyz2,gyz1)
    call transpose_y_to_x(gzz2,gzz1)

    sxz1(:,:,:)=half*(gzx1(:,:,:)+gxz1(:,:,:))

    srt_smag = zero
    srt_smag = sxx1 * sxx1 + syy1 * syy1 + szz1 * szz1 + two * sxy1 * sxy1 + two * sxz1 * sxz1 + two * syz1 * syz1

    nut1 = zero; nut2 = zero
    call transpose_x_to_y(srt_smag, srt_smag2)
    do k = 1, ysize(3)
       do j = 1, ysize(2)
          do i = 1, ysize(1)
             if(itype.eq.itype_abl) then
                !Mason and Thomson damping coefficient
                if (istret == 0) y=real(j+ystart(2)-1-1,mytype)*dy
                if (istret /= 0) y=yp(j+ystart(2)-1)
                smag_constant=(smagcst**(-nSmag)+(k_roughness*(y/del(j)+z_zero/del(j)))**(-nSmag))**(-one/nSmag)
                length=smag_constant*del(j)
             else
                length=smagcst*del(j)
             endif
             !Calculate eddy visc nu_t
             nut2(i, j, k) = ((length)**two) * sqrt(two * srt_smag2(i, j, k))
          enddo
       enddo
    enddo
    call transpose_y_to_x(nut2, nut1)


    if (mod(itime,ilist)==0) then 
      srtmin = minval(srt_smag)
      srtmax = maxval(srt_smag)
      nutmin = minval(nut1)
      nutmax = maxval(nut1)
      srtmin = zero
      srtmax = zero
      nutmin = zero
      nutmax = zero

      call MPI_ALLREDUCE(MPI_IN_PLACE,srtmin,1,real_type,MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,nutmin,1,real_type,MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,srtmax,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,nutmax,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)

      if (mod(itime, ilist) == 0) then
         if (nrank==0) then 
            write(*,*) "smag srt_smag min max= ", srtmin, srtmax  
            write(*,*) "smag nut1     min max= ", nutmin, nutmax
         endif
      endif
    endif

    if (mod(itime, ioutput).eq.0) then

#ifdef ADIOS2
       if (jles /= 3) then
          call decomp_2d_start_io(io_turb, turb_dir)
       end if
#endif
       call decomp_2d_write_one(1, nut1, turb_dir, &
            gen_filename(".", "nut_smag", itime / ioutput, "bin"), &
            2, io_turb)
#ifdef ADIOS2
       if (jles /= 3) then
          call decomp_2d_end_io(io_turb, turb_dir)
       end if
#endif

    endif

  end subroutine smag

  subroutine dynsmag(nut1,ux1,uy1,uz1,ep1)
    !================================================================================
    !
    !  SUBROUTINE: dynsmag
    ! DESCRIPTION: Calculates the eddy-viscosity nut according to the Lilly-Germano
    !              dynamic Smagorinsky model
    !      AUTHOR: G. Deskos <g.deskos14@imperial.ac.uk>
    !
    !================================================================================

    USE param
    USE variables
    USE decomp_2d_io
    USE MPI
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,di2
    USE var, only : ux3,uy3,uz3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
    USE var, only : sxx1,syy1,szz1,sxy1,sxz1,syz1,srt_smag
    USE var, only : sxx2,syy2,szz2,sxy2,sxz2,syz2,srt_smag2,nut2
    USE var, only : sxx3,syy3,szz3,sxy3,sxz3,syz3
    use tools, only : mean_plane_z
    USE ibm_param
    USE param, only : zero
    
    implicit none

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: nut1

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux1f, uy1f, uz1f
    real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: ux2f, uy2f, uz2f, ep2
    real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: ux3f, uy3f, uz3f

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: uxx1, uyy1, uzz1, uxy1, uxz1, uyz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: uxx1f, uyy1f, uzz1f, uxy1f, uxz1f, uyz1f
    real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: uxx2f, uyy2f, uzz2f, uxy2f, uxz2f, uyz2f
    real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: uxx3f, uyy3f, uzz3f, uxy3f, uxz3f, uyz3f


    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: sxx1f, syy1f, szz1f, sxy1f, sxz1f, syz1f
    real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: syy2f, szz2f, sxy2f, syz2f
    real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: szz3f, syz3f

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: gxx1f, gyx1f, gzx1f, gxz1f
    real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: gxy2f, gyy2f, gzy2f, gxz2f
    real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: gxz3f, gyz3f, gzz3f

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: axx1, ayy1, azz1, axy1, axz1, ayz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: axx1f, ayy1f, azz1f, axy1f, axz1f, ayz1f
    real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: axx2f, ayy2f, azz2f, axy2f, axz2f, ayz2f
    real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: axx3f, ayy3f, azz3f, axy3f, axz3f, ayz3f

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: bbxx1, bbyy1, bbzz1, bbxy1, bbxz1, bbyz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: lxx1, lyy1, lzz1, lxy1, lxz1, lyz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: mxx1, myy1, mzz1, mxy1, mxz1, myz1

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: smagC1, smagC1f, dsmagcst1
    real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: smagC2, smagC2f, dsmagcst2
    real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: smagC3, smagC3f, dsmagcst3

    integer :: i,j,k

    character(len = 30) :: filename

    if((iibm==1).or.(iibm==2).or.(iibm==3)) then
       ta1 = ux1 * (one - ep1)
       tb1 = uy1 * (one - ep1)
       tc1 = uz1 * (one - ep1)
    else
       ta1 = ux1
       tb1 = uy1
       tc1 = uz1
    endif

    if (itime.eq.ifirst) call smag(nut1,ux1,uy1,uz1)

    uxx1 = ta1 * ta1
    uyy1 = tb1 * tb1
    uzz1 = tc1 * tc1
    uxy1 = ta1 * tb1
    uxz1 = ta1 * tc1
    uyz1 = tb1 * tc1

    ! Initialise the filter
    call filter(zero)

    call filx(ux1f, ta1, di1,fisx,fiffx ,fifsx ,fifwx ,xsize(1),xsize(2),xsize(3),0,ubcx) !ux1
    call filx(uy1f, tb1, di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1,ubcy) !uy1
    call filx(uz1f, tc1, di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1,ubcz) !uz1

    call filx(uxx1f, uxx1, di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1,ubcx*ubcx) !ux1*ux1
    call filx(uyy1f, uyy1, di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1,ubcy*ubcy) !uy1*uy1
    call filx(uzz1f, uzz1, di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1,ubcz*ubcz) !uz1*uz1

    call filx(uxy1f, uxy1, di1,fisx,fiffx ,fifsx ,fifwx ,xsize(1),xsize(2),xsize(3),0,ubcx*ubcy) !ux1*uy1
    call filx(uxz1f, uxz1, di1,fisx,fiffx ,fifsx ,fifwx ,xsize(1),xsize(2),xsize(3),0,ubcx*ubcz) !ux1*uz1
    call filx(uyz1f, uyz1, di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1,ubcy*ubcz) !uy1*uz1

    if (mod(itime, ilist) == 0) then
       if (nrank==0) write(*,*) "filx ux= ", maxval(ta1), maxval(ux1f), maxval(ta1) - maxval(ux1f)
    endif

   if((iibm==1).or.(iibm==2).or.(iibm==3)) then
       ux1f = ux1f * (one - ep1)
       uy1f = uy1f * (one - ep1)
       uz1f = uz1f * (one - ep1)
       uxx1f = uxx1f * (one - ep1)
       uyy1f = uyy1f * (one - ep1)
       uzz1f = uzz1f * (one - ep1)
       uxy1f = uxy1f * (one - ep1)
       uxz1f = uxz1f * (one - ep1)
       uyz1f = uyz1f * (one - ep1)
       call transpose_x_to_y(ep1, ep2)
    endif

    call transpose_x_to_y(ux1f, ta2)
    call transpose_x_to_y(uy1f, tb2)
    call transpose_x_to_y(uz1f, tc2)
    call transpose_x_to_y(uxx1f, td2)
    call transpose_x_to_y(uyy1f, te2)
    call transpose_x_to_y(uzz1f, tf2)
    call transpose_x_to_y(uxy1f, tg2)
    call transpose_x_to_y(uxz1f, th2)
    call transpose_x_to_y(uyz1f, ti2)

    call fily(ux2f, ta2, di2,fisy,fiffyp,fifsyp,fifwyp,ysize(1),ysize(2),ysize(3),1,ubcx) !ux2
    call fily(uy2f, tb2, di2,fisy,fiffy ,fifsy ,fifwy ,ysize(1),ysize(2),ysize(3),0,ubcy) !uy2
    call fily(uz2f, tc2, di2,fisy,fiffyp,fifsyp,fifwyp,ysize(1),ysize(2),ysize(3),1,ubcz) !uz2

    call fily(uxx2f, td2, di2,fisy,fiffyp,fifsyp,fifwyp,ysize(1),ysize(2),ysize(3),1,ubcx*ubcx) !ux2*ux2
    call fily(uyy2f, te2, di2,fisy,fiffyp,fifsyp,fifwyp,ysize(1),ysize(2),ysize(3),1,ubcy*ubcy) !uy2*uy2
    call fily(uzz2f, tf2, di2,fisy,fiffyp,fifsyp,fifwyp,ysize(1),ysize(2),ysize(3),1,ubcz*ubcz) !uz2*uz2

    call fily(uxy2f, tg2, di2,fisy,fiffy ,fifsy ,fifwy ,ysize(1),ysize(2),ysize(3),0,ubcx*ubcy) !ux2*uy2
    call fily(uxz2f, th2, di2,fisy,fiffyp,fifsyp,fifwyp,ysize(1),ysize(2),ysize(3),1,ubcx*ubcz) !ux2*uz2
    call fily(uyz2f, ti2, di2,fisy,fiffy ,fifsy ,fifwy ,ysize(1),ysize(2),ysize(3),0,ubcy*ubcz) !uy2*uz2

    if (mod(itime, ilist) == 0) then
       if (nrank==0) write(*,*) "fily ux= ", maxval(ta2), maxval(ux2f), maxval(ta2) - maxval(ux2f)
    endif

    if((iibm==1).or.(iibm==2).or.(iibm==3)) then   
       ux2f = ux2f * (one - ep2)
       uy2f = uy2f * (one - ep2)
       uz2f = uz2f * (one - ep2)
       uxx2f = uxx2f * (one - ep2)
       uyy2f = uyy2f * (one - ep2)
       uzz2f = uzz2f * (one - ep2)
       uxy2f = uxy2f * (one - ep2)
       uxz2f = uxz2f * (one - ep2)
       uyz2f = uyz2f * (one - ep2)
    endif

    ta2 = zero; tb2 = zero; tc2 = zero
    td2 = zero; te2 = zero; tf2 = zero
    tg2 = zero; th2 = zero; ti2 = zero

    call transpose_y_to_z(ux2f, ta3)
    call transpose_y_to_z(uy2f, tb3)
    call transpose_y_to_z(uz2f, tc3)
    call transpose_y_to_z(uxx2f, td3)
    call transpose_y_to_z(uyy2f, te3)
    call transpose_y_to_z(uzz2f, tf3)
    call transpose_y_to_z(uxy2f, tg3)
    call transpose_y_to_z(uxz2f, th3)
    call transpose_y_to_z(uyz2f, ti3)

    call filz(ux3f, ta3, di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1,ubcx) !ux3
    call filz(uy3f, tb3, di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1,ubcy) !uy3
    call filz(uz3f, tc3, di3,fisz,fiffz ,fifsz ,fifwz ,zsize(1),zsize(2),zsize(3),0,ubcz) !uz3

    call filz(uxx3f, td3, di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1,ubcx*ubcx) !ux3*ux3
    call filz(uyy3f, te3, di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1,ubcy*ubcy) !uy3*uy3
    call filz(uzz3f, tf3, di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1,ubcz*ubcz) !uz3*uz3

    call filz(uxy3f, tg3, di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1,ubcx*ubcy) !ux3*uy3
    call filz(uxz3f, th3, di3,fisz,fiffz ,fifsz ,fifwz ,zsize(1),zsize(2),zsize(3),0,ubcx*ubcz) !ux3*uz3
    call filz(uyz3f, ti3, di3,fisz,fiffz ,fifsz ,fifwz ,zsize(1),zsize(2),zsize(3),0,ubcy*ubcz) !uy3*uz3

    if (mod(itime, ilist) == 0) then
       if (nrank==0) write(*,*) "filz ux= ", maxval(ta3), maxval(ux3f), maxval(ta3) - maxval(ux3f)
    endif

    ta3 = zero; tb3 = zero; tc3 = zero
    td3 = zero; te3 = zero; tf3 = zero
    tg3 = zero; th3 = zero; ti3 = zero

    ux2f = zero;uy2f = zero;uz2f = zero

    call transpose_z_to_y(ux3f, ux2f)
    call transpose_z_to_y(uy3f, uy2f)
    call transpose_z_to_y(uz3f, uz2f)
    call transpose_z_to_y(uxx3f, uxx2f)
    call transpose_z_to_y(uyy3f, uyy2f)
    call transpose_z_to_y(uzz3f, uzz2f)
    call transpose_z_to_y(uxy3f, uxy2f)
    call transpose_z_to_y(uxz3f, uxz2f)
    call transpose_z_to_y(uyz3f, uyz2f)

    ux1f = zero;uy1f = zero;uz1f = zero

    call transpose_y_to_x(ux2f, ux1f)
    call transpose_y_to_x(uy2f, uy1f)
    call transpose_y_to_x(uz2f, uz1f)
    call transpose_y_to_x(uxx2f, uxx1f)
    call transpose_y_to_x(uyy2f, uyy1f)
    call transpose_y_to_x(uzz2f, uzz1f)
    call transpose_y_to_x(uxy2f, uxy1f)
    call transpose_y_to_x(uxz2f, uxz1f)
    call transpose_y_to_x(uyz2f, uyz1f)

    if((iibm==1).or.(iibm==2).or.(iibm==3)) then
       ux1f = ux1f * (one - ep1)
       uy1f = uy1f * (one - ep1)
       uz1f = uz1f * (one - ep1)
       uxx1f = uxx1f * (one - ep1)
       uyy1f = uyy1f * (one - ep1)
       uzz1f = uzz1f * (one - ep1)
       uxy1f = uxy1f * (one - ep1)
       uxz1f = uxz1f * (one - ep1)
       uyz1f = uyz1f * (one - ep1)
    endif

    !Lij tensor OK
    lxx1 = uxx1f - ux1f * ux1f
    lyy1 = uyy1f - uy1f * uy1f
    lzz1 = uzz1f - uz1f * uz1f
    lxy1 = uxy1f - ux1f * uy1f
    lxz1 = uxz1f - ux1f * uz1f
    lyz1 = uyz1f - uy1f * uz1f

    call derx (gxx1f, ux1f, di1, sx, ffx, fsx, fwx, xsize(1), xsize(2), xsize(3), 0,ubcx)
    call derx (gyx1f, uy1f, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1,ubcy)
    call derx (gzx1f, uz1f, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1,ubcz)

    sxx1f = gxx1f

    !WORK Y-PENCILS
    call transpose_x_to_y(ux1f, ux2f)
    call transpose_x_to_y(uy1f, uy2f)
    call transpose_x_to_y(uz1f, uz2f)
    call transpose_x_to_y(gyx1f, ta2)

    call dery (gxy2f, ux2f, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1,ubcx)
    call dery (gyy2f, uy2f, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0,ubcy)
    call dery (gzy2f, uz2f, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1,ubcz)

    sxy2f = half * (gxy2f + ta2)
    syy2f = gyy2f

    !WORK Z-PENCILS
    call transpose_y_to_z(ux2f, ux3f)
    call transpose_y_to_z(uy2f, uy3f)
    call transpose_y_to_z(uz2f, uz3f)
    call transpose_y_to_z(gzy2f, ta3)

    call derz(gxz3f, ux3f, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1,ubcx)
    call derz(gyz3f, uy3f, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1,ubcy)
    call derz(gzz3f, uz3f, di3, sz, ffz, fsz, fwz, zsize(1), zsize(2), zsize(3), 0,ubcz)

    szz3f = gzz3f
    syz3f = half * (gyz3f + ta3)

    !WORK Y-PENCILS
    call transpose_z_to_y(syz3f, syz2f)
    call transpose_z_to_y(szz3f, szz2f)
    call transpose_z_to_y(gxz3f, gxz2f)

    !WORK X-PENCILS
    call transpose_y_to_x(sxy2f, sxy1f)
    call transpose_y_to_x(syy2f, syy1f)
    call transpose_y_to_x(syz2f, syz1f)
    call transpose_y_to_x(szz2f, szz1f)
    call transpose_y_to_x(gxz2f, gxz1f)

    sxz1f = half * (gzx1f + gxz1f)

    !Bij tensor with u test filtered OK
    bbxx1 = -two * sqrt(two * (sxx1f * sxx1f + syy1f * syy1f + szz1f * szz1f + two * sxy1f * sxy1f + two * sxz1f * sxz1f + two * syz1f * syz1f)) * sxx1f
    bbyy1 = -two * sqrt(two * (sxx1f * sxx1f + syy1f * syy1f + szz1f * szz1f + two * sxy1f * sxy1f + two * sxz1f * sxz1f + two * syz1f * syz1f)) * syy1f
    bbzz1 = -two * sqrt(two * (sxx1f * sxx1f + syy1f * syy1f + szz1f * szz1f + two * sxy1f * sxy1f + two * sxz1f * sxz1f + two * syz1f * syz1f)) * szz1f

    bbxy1 = -two * sqrt(two * (sxx1f * sxx1f + syy1f * syy1f + szz1f * szz1f + two * sxy1f * sxy1f + two * sxz1f * sxz1f + two * syz1f * syz1f)) * sxy1f
    bbxz1 = -two * sqrt(two * (sxx1f * sxx1f + syy1f * syy1f + szz1f * szz1f + two * sxy1f * sxy1f + two * sxz1f * sxz1f + two * syz1f * syz1f)) * sxz1f
    bbyz1 = -two * sqrt(two * (sxx1f * sxx1f + syy1f * syy1f + szz1f * szz1f + two * sxy1f * sxy1f + two * sxz1f * sxz1f + two * syz1f * syz1f)) * syz1f

    call smag(nut1,ux1,uy1,uz1) ! to update the Aij tensor

    !Aij tensor with u
    axx1 = -two * sqrt(two * (sxx1 * sxx1 + syy1 * syy1 + szz1 * szz1 + two * sxy1 * sxy1 + two * sxz1 * sxz1 + two * syz1 * syz1)) * sxx1
    ayy1 = -two * sqrt(two * (sxx1 * sxx1 + syy1 * syy1 + szz1 * szz1 + two * sxy1 * sxy1 + two * sxz1 * sxz1 + two * syz1 * syz1)) * syy1
    azz1 = -two * sqrt(two * (sxx1 * sxx1 + syy1 * syy1 + szz1 * szz1 + two * sxy1 * sxy1 + two * sxz1 * sxz1 + two * syz1 * syz1)) * szz1

    axy1 = -two * sqrt(two * (sxx1 * sxx1 + syy1 * syy1 + szz1 * szz1 + two * sxy1 * sxy1 + two * sxz1 * sxz1 + two * syz1 * syz1)) * sxy1
    axz1 = -two * sqrt(two * (sxx1 * sxx1 + syy1 * syy1 + szz1 * szz1 + two * sxy1 * sxy1 + two * sxz1 * sxz1 + two * syz1 * syz1)) * sxz1
    ayz1 = -two * sqrt(two * (sxx1 * sxx1 + syy1 * syy1 + szz1 * szz1 + two * sxy1 * sxy1 + two * sxz1 * sxz1 + two * syz1 * syz1)) * syz1

    if((iibm==1).or.(iibm==2).or.(iibm==3)) then
       bbxx1 = bbxx1 * (one - ep1)
       bbyy1 = bbyy1 * (one - ep1)
       bbzz1 = bbzz1 * (one - ep1)
       bbxy1 = bbxy1 * (one - ep1)
       bbxz1 = bbxz1 * (one - ep1)
       bbyz1 = bbyz1 * (one - ep1)
       axx1 = axx1 * (one - ep1)
       ayy1 = ayy1 * (one - ep1)
       azz1 = azz1 * (one - ep1)
       axy1 = axy1 * (one - ep1)
       axz1 = axz1 * (one - ep1)
       ayz1 = ayz1 * (one - ep1)
    endif

    call transpose_x_to_y(axx1, ta2)
    call transpose_x_to_y(ayy1, tb2)
    call transpose_x_to_y(azz1, tc2)
    call transpose_x_to_y(axy1, td2)
    call transpose_x_to_y(axz1, te2)
    call transpose_x_to_y(ayz1, tf2)

    do j = 1, ysize(2)
       ta2(:, j, :) = ta2(:, j, :) * (del(j))**two
       tb2(:, j, :) = tb2(:, j, :) * (del(j))**two
       tc2(:, j, :) = tc2(:, j, :) * (del(j))**two
       td2(:, j, :) = td2(:, j, :) * (del(j))**two
       te2(:, j, :) = te2(:, j, :) * (del(j))**two
       tf2(:, j, :) = tf2(:, j, :) * (del(j))**two
    enddo

    call transpose_y_to_x(ta2, axx1)
    call transpose_y_to_x(tb2, ayy1)
    call transpose_y_to_x(tc2, azz1)
    call transpose_y_to_x(td2, axy1)
    call transpose_y_to_x(te2, axz1)
    call transpose_y_to_x(tf2, ayz1)

    call transpose_x_to_y(bbxx1, ta2)
    call transpose_x_to_y(bbyy1, tb2)
    call transpose_x_to_y(bbzz1, tc2)
    call transpose_x_to_y(bbxy1, td2)
    call transpose_x_to_y(bbxz1, te2)
    call transpose_x_to_y(bbyz1, tf2)

    do j = 1, ysize(2)
       ta2(:, j, :) = ta2(:, j, :) * (two * del(j))**two
       tb2(:, j, :) = tb2(:, j, :) * (two * del(j))**two
       tc2(:, j, :) = tc2(:, j, :) * (two * del(j))**two
       td2(:, j, :) = td2(:, j, :) * (two * del(j))**two
       te2(:, j, :) = te2(:, j, :) * (two * del(j))**two
       tf2(:, j, :) = tf2(:, j, :) * (two * del(j))**two
    enddo

    call transpose_y_to_x(ta2, bbxx1)
    call transpose_y_to_x(tb2, bbyy1)
    call transpose_y_to_x(tc2, bbzz1)
    call transpose_y_to_x(td2, bbxy1)
    call transpose_y_to_x(te2, bbxz1)
    call transpose_y_to_x(tf2, bbyz1)

    !Need to filter Aij components

    call filx(axx1f, axx1, di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1,zero)
    call filx(ayy1f, ayy1, di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1,zero)
    call filx(azz1f, azz1, di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1,zero)

    call filx(axy1f, axy1, di1,fisx,fiffx ,fifsx ,fifwx ,xsize(1),xsize(2),xsize(3),0,zero)
    call filx(axz1f, axz1, di1,fisx,fiffx ,fifsx ,fifwx ,xsize(1),xsize(2),xsize(3),0,zero)
    call filx(ayz1f, ayz1, di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1,zero)

    if (mod(itime, ilist) == 0) then
       if (nrank==0) write(*,*) "filx axx1= ", maxval(axx1), maxval(axx1f), maxval(axx1) - maxval(axx1f)
    endif

    if((iibm==1).or.(iibm==2).or.(iibm==3)) then
       axx1f = axx1f * (one - ep1)
       ayy1f = ayy1f * (one - ep1)
       azz1f = azz1f * (one - ep1)
       axy1f = axy1f * (one - ep1)
       axz1f = axz1f * (one - ep1)
       ayz1f = ayz1f * (one - ep1)
    endif

    call transpose_x_to_y(axx1f, ta2)
    call transpose_x_to_y(ayy1f, tb2)
    call transpose_x_to_y(azz1f, tc2)
    call transpose_x_to_y(axy1f, td2)
    call transpose_x_to_y(axz1f, te2)
    call transpose_x_to_y(ayz1f, tf2)

    call fily(axx2f, ta2, di2,fisy,fiffyp,fifsyp,fifwyp,ysize(1),ysize(2),ysize(3),1,zero)
    call fily(ayy2f, tb2, di2,fisy,fiffyp,fifsyp,fifwyp,ysize(1),ysize(2),ysize(3),1,zero)
    call fily(azz2f, tc2, di2,fisy,fiffyp,fifsyp,fifwyp,ysize(1),ysize(2),ysize(3),1,zero)

    call fily(axy2f, td2, di2,fisy,fiffy ,fifsy ,fifwy ,ysize(1),ysize(2),ysize(3),0,zero)
    call fily(axz2f, te2, di2,fisy,fiffyp,fifsyp,fifwyp,ysize(1),ysize(2),ysize(3),1,zero)
    call fily(ayz2f, tf2, di2,fisy,fiffy ,fifsy ,fifwy ,ysize(1),ysize(2),ysize(3),0,zero)

    if (mod(itime, ilist) == 0) then
       if (nrank==0) write(*,*) "fily axx2= ", maxval(ta2), maxval(axx2f), maxval(ta2) - maxval(axx2f)
    endif

    if((iibm==1).or.(iibm==2).or.(iibm==3)) then
       axx2f = axx2f * (one - ep2)
       ayy2f = ayy2f * (one - ep2)
       azz2f = azz2f * (one - ep2)
       axy2f = axy2f * (one - ep2)
       axz2f = axz2f * (one - ep2)
       ayz2f = ayz2f * (one - ep2)
    endif

    ta2 = zero; tb2 = zero; tc2 = zero
    td2 = zero; te2 = zero; tf2 = zero

    call transpose_y_to_z(axx2f, ta3)
    call transpose_y_to_z(ayy2f, tb3)
    call transpose_y_to_z(azz2f, tc3)
    call transpose_y_to_z(axy2f, td3)
    call transpose_y_to_z(axz2f, te3)
    call transpose_y_to_z(ayz2f, tf3)


    call filz(axx3f, ta3, di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1,zero)
    call filz(ayy3f, tb3, di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1,zero)
    call filz(azz3f, tc3, di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1,zero)

    call filz(axy3f, td3, di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1,zero)
    call filz(axz3f, te3, di3,fisz,fiffz ,fifsz ,fifwz ,zsize(1),zsize(2),zsize(3),0,zero)
    call filz(ayz3f, tf3, di3,fisz,fiffz ,fifsz ,fifwz ,zsize(1),zsize(2),zsize(3),0,zero)

    if (mod(itime, ilist) == 0) then
       if (nrank==0) write(*,*) "filz axx3= ", maxval(ta3), maxval(axx3f), maxval(ta3) - maxval(axx3f)
    endif

    ta3 = zero; tb3 = zero; tc3 = zero
    td3 = zero; te3 = zero; tf3 = zero

    call transpose_z_to_y(axx3f, axx2f)
    call transpose_z_to_y(ayy3f, ayy2f)
    call transpose_z_to_y(azz3f, azz2f)
    call transpose_z_to_y(axy3f, axy2f)
    call transpose_z_to_y(axz3f, axz2f)
    call transpose_z_to_y(ayz3f, ayz2f)

    call transpose_y_to_x(axx2f, axx1f)
    call transpose_y_to_x(ayy2f, ayy1f)
    call transpose_y_to_x(azz2f, azz1f)
    call transpose_y_to_x(axy2f, axy1f)
    call transpose_y_to_x(axz2f, axz1f)
    call transpose_y_to_x(ayz2f, ayz1f)

    if((iibm==1).or.(iibm==2).or.(iibm==3)) then
       axx1f = axx1f * (one - ep1)
       ayy1f = ayy1f * (one - ep1)
       azz1f = azz1f * (one - ep1)
       axy1f = axy1f * (one - ep1)
       axz1f = axz1f * (one - ep1)
       ayz1f = ayz1f * (one - ep1)
    endif

    !Mij tensor OK !!!M_{ij} = B_{ij} - Aij
    mxx1 = bbxx1 - axx1f
    myy1 = bbyy1 - ayy1f
    mzz1 = bbzz1 - azz1f
    mxy1 = bbxy1 - axy1f
    mxz1 = bbxz1 - axz1f
    myz1 = bbyz1 - ayz1f

    !Lij deviator
    lxx1 = lxx1 - (lxx1 + lyy1 + lzz1) / three
    lyy1 = lyy1 - (lxx1 + lyy1 + lzz1) / three
    lzz1 = lzz1 - (lxx1 + lyy1 + lzz1) / three

    if((iibm==1).or.(iibm==2).or.(iibm==3)) then
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                if (ep1(i, j, k) .eq. one) then
                   ta1(i, j, k) = zero
                   tb1(i, j, k) = one
                endif
             enddo
          enddo
       enddo
    endif

    !MODEL OF LILLY (1992)
    smagC1 = (lxx1 * mxx1 + lyy1 * myy1 + lzz1 * mzz1 + two * (lxy1 * mxy1 + lxz1 * mxz1 + lyz1 * myz1)) / &
         (mxx1 * mxx1 + myy1 * myy1 + mzz1 * mzz1 + two * (mxy1 * mxy1 + mxz1 * mxz1 + myz1 * myz1)) !l/M

    ! Limiter for the dynamic Smagorinsky constant
    do k = 1, xsize(3)
       do j = 1, xsize(2)
          do i = 1, xsize(1)
             if (smagC1(i, j, k).gt. maxdsmagcst) smagC1(i, j, k) = maxdsmagcst!zero
             if (smagC1(i, j, k).lt. 0.0) smagC1(i, j, k) = zero
          enddo
       enddo
    enddo

    !FILTERING THE NON-CONSTANT CONSTANT
    call filx(smagC1f, smagC1, di1,fisx,fiffx ,fifsx ,fifwx ,xsize(1),xsize(2),xsize(3),0,zero)

    call transpose_x_to_y(smagC1f, ta2)
    call fily(smagC2f, ta2, di2,fisy,fiffy ,fifsy ,fifwy ,ysize(1),ysize(2),ysize(3),0,zero)

    call transpose_y_to_z(smagC2f, ta3)
    call filz(smagC3f, ta3, di3,fisz,fiffz ,fifsz ,fifwz ,zsize(1),zsize(2),zsize(3),0,zero)

    if (mod(itime, ilist) == 0) then
       if (nrank==0) write(*,*) "filx smagC1= ", maxval(smagC1), maxval(smagC1f), maxval(smagC1) - maxval(smagC1f)
       if (nrank==0) write(*,*) "fily smagC1= ", maxval(ta2), maxval(smagC2f), maxval(ta2) - maxval(smagC2f)
       if (nrank==0) write(*,*) "filz smagC1= ", maxval(ta3), maxval(smagC3f), maxval(ta3) - maxval(smagC3f)
    endif


    dsmagcst3 = zero
    call mean_plane_z(smagC3f, zsize(1), zsize(2), zsize(3), dsmagcst3(:, :, 1))

    do k = 2, zsize(3)
       dsmagcst3(:, :, k) = dsmagcst3(:, :, 1)
    enddo

    call transpose_z_to_y(dsmagcst3, dsmagcst2)
    !call transpose_z_to_y(smagC3f,dsmagcst2)

    call transpose_y_to_x(dsmagcst2, dsmagcst1)

    ! do ijk = 1, nvect1 !ERIC LIMITEUR SI BESOIN
    !   if (dsmagcst1(ijk, 1, 1).gt. maxdsmagcst) dsmagcst1(ijk, 1, 1) =zero
    !   if (dsmagcst1(ijk, 1, 1).lt. 0.0) dsmagcst1(ijk, 1, 1) = zero
    ! enddo

    nut1 = zero; nut2 = zero

    ! Use the standard smagorinsky to calculate the srt_smag
#ifdef ADIOS2
    if (mod(itime, ioutput) == 0) then
       ! Start IO for writing in SMAG
       call decomp_2d_start_io(io_turb, turb_dir)
    end if
#endif

    call transpose_x_to_y(srt_smag, srt_smag2)
    call transpose_x_to_y(dsmagcst1, dsmagcst2)
    do k = 1, ysize(3)
       do j = 1, ysize(2)
          do i = 1, ysize(1)
             nut2(i, j, k) = dsmagcst2(i, j, k) * ((del(j))**two) * sqrt(two * srt_smag2(i, j, k))
          enddo
       enddo
    enddo
    call transpose_y_to_x(nut2, nut1)

    if (mod(itime,ilist)==0) then
       !if (nrank==0) write(*,*) "dsmagc init   min max= ", minval(smagC1), maxval(smagC1)
       if (nrank==0) write(*,*) "dsmagc final  min max= ", minval(dsmagcst1), maxval(dsmagcst1)
       if (nrank==0) write(*,*) "dsmag nut1    min max= ", minval(nut1), maxval(nut1)
    endif

    if (mod(itime, ioutput) == 0) then

      ! write(filename, "('./data/dsmagcst_initial',I4.4)") itime / imodulo
      ! call decomp_2d_write_one(1, smagC1, filename, 2)

       call decomp_2d_write_one(1, dsmagcst1, turb_dir, gen_filename(".", "dsmagcst_final", itime / ioutput, "bin"), 2, io_turb)
       call decomp_2d_write_one(1, nut1, turb_dir, gen_filename(".", "nut_dynsmag", itime / ioutput, "bin"), 2, io_turb)
#ifdef ADIOS2
       call decomp_2d_end_io(io_turb, turb_dir)
#endif
    endif

  end subroutine dynsmag

  subroutine wale(nut1,ux1,uy1,uz1)
  !================================================================================
  !
  !  SUBROUTINE: wale
  ! DESCRIPTION: Calculates the eddy-viscosity nut according to the wall-adapting
  !              local eddy-viscosity (WALE) model:
  !
  !              Nicoud, F. and Ducros, F., 1999. Subgrid-scale stress modelling
  !              based on the square of the velocity gradient tensor.
  !              Flow, turbulence and Combustion, 62(3), pp.183-200.
  !
  !      AUTHOR: Arash Hamzehloo <a.hamzehloo@imperial.ac.uk>
  !
  !================================================================================

  USE param
  USE variables
  USE decomp_2d_io
  USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
  USE var, only : ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,di2
  USE var, only : ux3,uy3,uz3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
  USE var, only : sxx1,syy1,szz1,sxy1,sxz1,syz1
  USE var, only : gxx1,gyx1,gzx1,gxy2,gyy2,gzy2,gxz3,gyz3,gzz3
  USE var, only : gxy1,gyy1,gzy1,gxz2,gyz2,gzz2,gxz1,gyz1,gzz1
  USE var, only : sxx2,syy2,szz2,sxy2,sxz2,syz2,srt_smag2,nut2
  USE var, only : sxx3,syy3,szz3,sxy3,sxz3,syz3
  USE var, only : sdxx1,sdyy1,sdzz1,sdxy1,sdxz1,sdyz1
  USE var, only : sdxx2,sdyy2,sdzz2,sdxy2,sdxz2,sdyz2
  USE var, only : sdxx3,sdyy3,sdzz3,sdxy3,sdxz3,sdyz3
  USE var, only : srt_wale,srt_wale2,srt_wale3,srt_wale4
  USE ibm_param
  
  implicit none

  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: nut1

  integer :: i, j, k
  character(len = 30) :: filename


  ! INFO about the auxillary arrays
  !--------------------------------------------------------
  ! gxx= dux/dx; gyx=duy/dx; gzx=duz/dx;
  ! gxy= dux/dy; gyy=duy/dy; gzy=duz/dy;
  ! gxz= dux/dz; gyz=duy/dz; gzz=duz/dz

  call derx (gxx1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)
  call derx (gyx1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
  call derx (gzx1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz)

  sxx1(:,:,:) = gxx1(:,:,:)

  !WORK Y-PENCILS
  call transpose_x_to_y(ux1,ux2)
  call transpose_x_to_y(uy1,uy2)
  call transpose_x_to_y(uz1,uz2)
  call transpose_x_to_y(gyx1,ta2)

  call dery (gxy2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
  call dery (gyy2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
  call dery (gzy2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)

  syy2(:,:,:)=gyy2(:,:,:)
  sxy2(:,:,:)=half*(gxy2(:,:,:)+ta2(:,:,:))


  !WORK Z-PENCILS
  call transpose_y_to_z(ux2,ux3)
  call transpose_y_to_z(uy2,uy3)
  call transpose_y_to_z(uz2,uz3)
  call transpose_y_to_z(gzy2,ta3)

  call derz(gxz3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
  call derz(gyz3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcy)
  call derz(gzz3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcz)

  szz3(:,:,:)=gzz3(:,:,:)
  syz3(:,:,:)=half*(gyz3(:,:,:)+ta3(:,:,:))

  !WORK Y-PENCILS
  call transpose_z_to_y(syz3,syz2)
  call transpose_z_to_y(szz3,szz2)

  call transpose_z_to_y(sdyz3,sdyz2)
  call transpose_z_to_y(sdzz3,sdzz2)

  call transpose_z_to_y(gxz3,gxz2)
  call transpose_z_to_y(gyz3,gyz2)
  call transpose_z_to_y(gzz3,gzz2)

  !WORK X-PENCILS
  call transpose_y_to_x(sxy2,sxy1)
  call transpose_y_to_x(syy2,syy1)
  call transpose_y_to_x(syz2,syz1)
  call transpose_y_to_x(szz2,szz1)

  call transpose_y_to_x(sdxy2,sdxy1)
  call transpose_y_to_x(sdyy2,sdyy1)
  call transpose_y_to_x(sdyz2,sdyz1)
  call transpose_y_to_x(sdzz2,sdzz1)

  call transpose_y_to_x(gxy2,gxy1)
  call transpose_y_to_x(gyy2,gyy1)
  call transpose_y_to_x(gzy2,gzy1)
  call transpose_y_to_x(gxz2,gxz1)
  call transpose_y_to_x(gyz2,gyz1)
  call transpose_y_to_x(gzz2,gzz1)

  sxz1(:,:,:)=half*(gzx1(:,:,:)+gxz1(:,:,:))

  sdxx1(:,:,:)=gxx1(:,:,:)*gxx1(:,:,:)+gxy1(:,:,:)*gyx1(:,:,:)+&
       gxz1(:,:,:)*gzx1(:,:,:)-(one/three)*(gxx1(:,:,:)*gxx1(:,:,:)+&
       gyy1(:,:,:)*gyy1(:,:,:)+gzz1(:,:,:)*gzz1(:,:,:)+&
       two*gxy1(:,:,:)*gyx1(:,:,:)+two*gxz1(:,:,:)*gzx1(:,:,:)+&
       two*gzy1(:,:,:)*gyz1(:,:,:))

  sdxy1(:,:,:)=half*(gxx1(:,:,:)*gxy1(:,:,:)+gxy1(:,:,:)*gyy1(:,:,:)+&
       gxz1(:,:,:)*gzy1(:,:,:)+gyx1(:,:,:)*gxx1(:,:,:)+&
       gyy1(:,:,:)*gyx1(:,:,:)+gyz1(:,:,:)*gzx1(:,:,:))

  sdxz1(:,:,:)=half*(gxx1(:,:,:)*gxz1(:,:,:)+gxy1(:,:,:)*gyz1(:,:,:)+&
       gxz1(:,:,:)*gzz1(:,:,:)+gzx1(:,:,:)*gxx1(:,:,:)+&
       gzy1(:,:,:)*gyx1(:,:,:)+gzz1(:,:,:)*gzx1(:,:,:))

  sdyy1(:,:,:)=gyx1(:,:,:)*gxy1(:,:,:)+gyy1(:,:,:)*gyy1(:,:,:)+&
       gyz1(:,:,:)*gzy1(:,:,:)-(one/three)*(gxx1(:,:,:)*gxx1(:,:,:)+&
       gyy1(:,:,:)*gyy1(:,:,:)+gzz1(:,:,:)*gzz1(:,:,:)+&
       two*gxy1(:,:,:)*gyx1(:,:,:)+two*gxz1(:,:,:)*gzx1(:,:,:)+&
       two*gzy1(:,:,:)*gyz1(:,:,:))

  sdyz1(:,:,:)=half*(gyx1(:,:,:)*gxz1(:,:,:)+gyy1(:,:,:)*gyz1(:,:,:)+&
       gyz1(:,:,:)*gzz1(:,:,:)+gzx1(:,:,:)*gxy1(:,:,:)+&
       gzy1(:,:,:)*gyy1(:,:,:)+gzz1(:,:,:)*gzy1(:,:,:))

  sdzz1(:,:,:)=gzx1(:,:,:)*gxz1(:,:,:)+gzy1(:,:,:)*gyz1(:,:,:)+&
       gzz1(:,:,:)*gzz1(:,:,:)-(one/three)*(gxx1(:,:,:)*gxx1(:,:,:)+&
       gyy1(:,:,:)*gyy1(:,:,:)+gzz1(:,:,:)*gzz1(:,:,:)+&
       two*gxy1(:,:,:)*gyx1(:,:,:)+two*gxz1(:,:,:)*gzx1(:,:,:)+&
       two*gzy1(:,:,:)*gyz1(:,:,:))

  srt_wale = zero
  srt_wale2 = zero
  srt_wale3 = zero
  srt_wale4 = zero

  srt_wale = sxx1 * sxx1 + syy1 * syy1 + szz1 * szz1 + two * sxy1 * sxy1 &
           + two * sxz1 * sxz1 + two * syz1 * syz1

  srt_wale3 = sdxx1 * sdxx1 + sdyy1 * sdyy1 + sdzz1 * sdzz1 + two * sdxy1 * sdxy1 &
            + two * sdxz1 * sdxz1 + two * sdyz1 * sdyz1

  nut1 = zero; nut2 = zero
  call transpose_x_to_y(srt_wale, srt_wale2)
  call transpose_x_to_y(srt_wale3, srt_wale4)
  do k = 1, ysize(3)
     do j = 1, ysize(2)
        do i = 1, ysize(1)
           nut2(i, j, k) = ((walecst * del(j))**two) * ((srt_wale4(i, j, k)**(three/two))/((srt_wale2(i, j, k)**(five/two)) &
                         +(srt_wale4(i, j, k)**(five/four))))
        enddo
     enddo
  enddo
  call transpose_y_to_x(nut2, nut1)

  if (mod(itime, ilist) == 0) then
     if (nrank==0) write(*,*) "WALE SS min max= ", minval(srt_wale), maxval(srt_wale)
     if (nrank==0) write(*,*) "WALE SdSd min max= ", minval(srt_wale3), maxval(srt_wale3)
     if (nrank==0) write(*,*) "WALE nut1     min max= ", minval(nut1), maxval(nut1)
  endif

  if (mod(itime, ioutput).eq.0) then

#ifdef ADIOS2
     call decomp_2d_start_io(io_turb, turb_dir)
#endif
     call decomp_2d_write_one(1, nut1, turb_dir, gen_filename(".", "nut_wale", itime / ioutput, "bin"), 2, io_turb)
#ifdef ADIOS2
     call decomp_2d_end_io(io_turb, turb_dir)
#endif

  endif

end subroutine wale


  subroutine sgs_mom_nonconservative(sgsx1,sgsy1,sgsz1,ux1,uy1,uz1,nut1,ep1)
    !================================================================================
    !
    !  SUBROUTINE: sgs_mom_nonconservative
    ! DESCRIPTION: Calculates the divergence of the sub-grid-scale stresses
    !              using a non-conservative formulation
    !      AUTHOR: G. Deskos <g.deskos14@imperial.ac.uk>
    !
    !================================================================================

    USE param
    USE variables
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
    USE var, only : ux3,uy3,uz3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
    USE var, only : sgsx2,sgsy2,sgsz2,nut2
    USE var, only : sgsx3,sgsy3,sgsz3,nut3
    USE var, only : sxx1,sxy1,sxz1,syy1,syz1,szz1
    USE var, only : sxy2,syy2,syz2,sxz2,szz2,sxz3,syz3,szz3
    USE param, only : mytype, zero
    USE ibm_param
    
    implicit none

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)),intent(in) :: ux1, uy1, uz1, nut1, ep1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)),intent(out):: sgsx1, sgsy1, sgsz1

    integer :: i, j, k, ijk, nvect1

    ta1 = zero; ta2 = zero; ta3 = zero
    sgsx1=zero; sgsy1=zero; sgsz1=zero
    sgsx2=zero; sgsy2=zero; sgsz2=zero
    sgsx3=zero; sgsy3=zero; sgsz3=zero
    !WORK X-PENCILS
    call derx (ta1,nut1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,zero)

    call derxx (td1,ux1,di1,sx,sfx ,ssx ,swx ,xsize(1),xsize(2),xsize(3),0,zero)
    call derxx (te1,uy1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1,zero)
    call derxx (tf1,uz1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1,zero)

    sgsx1 = td1 * nut1 + two * sxx1 * ta1
    sgsy1 = te1 * nut1 + two * sxy1 * ta1
    sgsz1 = tf1 * nut1 + two * sxz1 * ta1


    !WORK Y-PENCILS
    call transpose_x_to_y(sgsx1, sgsx2)
    call transpose_x_to_y(sgsy1, sgsy2)
    call transpose_x_to_y(sgsz1, sgsz2)
    call transpose_x_to_y(sxz1, sxz2)
    call transpose_x_to_y(sxy1, sxy2)
    call transpose_x_to_y(syy1, syy2)
    call transpose_x_to_y(syz1, syz2)
    call transpose_x_to_y(szz1, szz2)
    call transpose_x_to_y(nut1, nut2)
    call transpose_x_to_y(ux1, ux2)
    call transpose_x_to_y(uy1, uy2)
    call transpose_x_to_y(uz1, uz2)

    call dery (ta2, nut2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1,zero)

    !-->for ux
    td2 = zero
    iimplicit = -iimplicit
    if (istret.ne.0) then
       call deryy (td2, ux2, di2, sy, sfyp, ssyp, swyp, ysize(1), ysize(2), ysize(3), 1,ubcx)
       call dery (te2, ux2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1,ubcx)
       do k = 1, ysize(3)
          do j = 1, ysize(2)
             do i = 1, ysize(1)
                td2(i, j, k) = td2(i, j, k) * pp2y(j) - pp4y(j) * te2(i, j, k)
             enddo
          enddo
       enddo
    else
       call deryy (td2, ux2, di2, sy, sfyp, ssyp, swyp, ysize(1), ysize(2), ysize(3), 1,ubcx)
    endif

    !-->for uy
    te2 = zero
    if (istret.ne.0) then
       call deryy (te2, uy2, di2, sy, sfy, ssy, swy, ysize(1), ysize(2), ysize(3), 0,ubcy)
       call dery (tf2, uy2, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0,ubcy)
       do k = 1, ysize(3)
          do j = 1, ysize(2)
             do i = 1, ysize(1)
                te2(i, j, k) = te2(i, j, k) * pp2y(j) - pp4y(j) * tf2(i, j, k)
             enddo
          enddo
       enddo
    else
       call deryy (te2, uy2, di2, sy, sfy, ssy, swy, ysize(1), ysize(2), ysize(3), 0,ubcy)
    endif

    !-->for uz
    tf2 = zero
    if (istret.ne.0) then
       call deryy (tf2, uz2, di2, sy, sfyp, ssyp, swyp, ysize(1), ysize(2), ysize(3), 1,ubcz)
       call dery (tj2, uz2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1,ubcz)
       do k = 1, ysize(3)
          do j = 1, ysize(2)
             do i = 1, ysize(1)
                tf2(i, j, k) = tf2(i, j, k) * pp2y(j) - pp4y(j) * tj2(i, j, k)
             enddo
          enddo
       enddo
    else
       call deryy (tf2, uz2, di2, sy, sfyp, ssyp, swyp, ysize(1), ysize(2), ysize(3), 1,ubcz)
    endif
    iimplicit = -iimplicit

    sgsx2 = sgsx2 + nut2 * td2 + two * sxy2 * ta2
    sgsy2 = sgsy2 + nut2 * te2 + two * syy2 * ta2
    sgsz2 = sgsz2 + nut2 * tf2 + two * syz2 * ta2

    !WORK Z-PENCILS
    call transpose_y_to_z(sgsx2, sgsx3)
    call transpose_y_to_z(sgsy2, sgsy3)
    call transpose_y_to_z(sgsz2, sgsz3)
    call transpose_y_to_z(sxz2, sxz3)
    call transpose_y_to_z(syz2, syz3)
    call transpose_y_to_z(szz2, szz3)
    call transpose_y_to_z(nut2, nut3)
    call transpose_y_to_z(ux2, ux3)
    call transpose_y_to_z(uy2, uy3)
    call transpose_y_to_z(uz2, uz3)

    call derz (ta3, nut3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1, zero)

    call derzz (td3, ux3, di3, sz, sfzp, sszp, swzp, zsize(1), zsize(2), zsize(3), 1, ubcx)
    call derzz (te3, uy3, di3, sz, sfzp, sszp, swzp, zsize(1), zsize(2), zsize(3), 1, ubcy)
    call derzz (tf3, uz3, di3, sz, sfz, ssz, swz, zsize(1), zsize(2), zsize(3), 0, ubcz)

    sgsx3 = sgsx3 + nut3 * td3 + two * sxz3 * ta3
    sgsy3 = sgsy3 + nut3 * te3 + two * syz3 * ta3
    sgsz3 = sgsz3 + nut3 * tf3 + two * szz3 * ta3

    call transpose_z_to_y(sgsx3, sgsx2)
    call transpose_z_to_y(sgsy3, sgsy2)
    call transpose_z_to_y(sgsz3, sgsz2)

    call transpose_y_to_x(sgsx2, sgsx1)
    call transpose_y_to_x(sgsy2, sgsy1)
    call transpose_y_to_x(sgsz2, sgsz1)

    if((iibm==1).or.(iibm==2).or.(iibm==3)) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             do i=1,xsize(1)
                if(ep1(i,j, k).eq.1) then
                   sgsx1(i,j,k) = zero
                   sgsy1(i,j,k) = zero
                   sgsz1(i,j,k) = zero
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine sgs_mom_nonconservative

  !************************************************************
  subroutine sgs_scalar_nonconservative(sgsphi1,nut1,phi1,is)

    USE param
    USE variables

    USE var, only: di1,tb1,di2,tb2,di3,tb3,tc1,tc2,tc3
    USE abl, only: wall_sgs_slip_scalar

    implicit none

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: sgsphi1, phi1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: dphidy1
    real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: phi2, sgsphi2
    real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: phi3, sgsphi3

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: nut1, dnut1
    real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: nut2, dnut2
    real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: nut3, dnut3

    real(mytype) :: Pr
    integer :: is, i, j, k

    sgsphi1 = zero; sgsphi2 = zero; sgsphi3 = zero

    call derxS (dnut1, nut1, di1, sx, ffxpS, fsxpS, fwxpS, xsize(1), xsize(2), xsize(3), 1, zero)
    call transpose_x_to_y(nut1, nut2)
    call deryS (dnut2, nut2, di2, sy, ffypS, fsypS, fwypS, ppy, ysize(1), ysize(2), ysize(3), 1, zero)
    call transpose_y_to_z(nut2, nut3)
    call derzS (dnut3, nut3, di3, sz, ffzpS, fszpS, fwzpS, zsize(1), zsize(2), zsize(3), 1, zero)

    ! kappat = nut/Pr
    Pr = Sc(is)

    call derxS (tb1, phi1, di1, sx, ffxpS, fsxpS, fwxpS, xsize(1), xsize(2), xsize(3), 1, zero)
    call derxxS(tc1, phi1, di1, sx, sfxpS, ssxpS, swxpS, xsize(1), xsize(2), xsize(3), 1, zero)
    sgsphi1 = tb1 * (dnut1/Pr) + tc1 * (nut1/Pr)

    call transpose_x_to_y(phi1, phi2)
    call transpose_x_to_y(sgsphi1, sgsphi2)

    call deryS (tb2, phi2, di2, sy, ffypS, fsypS, fwypS, ppy, ysize(1), ysize(2), ysize(3), 1, zero)
    iimplicit = - iimplicit
    call deryyS(tc2, phi2, di2, sy, sfypS, ssypS, swypS, ysize(1), ysize(2), ysize(3), 1, zero)
    iimplicit = - iimplicit
    if (istret.ne.0) then
       do k = 1, ysize(3)
       do j = 1, ysize(2)
       do i = 1, ysize(1)
          tc2(i,j,k) = tc2(i,j,k)*pp2y(j)-pp4y(j)*tb2(i,j,k)
       enddo
       enddo
       enddo
    endif
    sgsphi2 = sgsphi2 + tb2 * (dnut2/Pr) + tc2 * (nut2/Pr)

    call transpose_y_to_z(phi2, phi3)
    call transpose_y_to_z(sgsphi2, sgsphi3)

    call derzS (tb3, phi3, di3, sz, ffzpS, fszpS, fwzpS, zsize(1), zsize(2), zsize(3), 1, zero)
    call derzzS(tc3, phi3, di3, sz, sfzpS, sszpS, swzpS, zsize(1), zsize(2), zsize(3), 1, zero)
    sgsphi3 = sgsphi3 + tb3 * (dnut3/Pr) + tc3 * (nut3/Pr)

    call transpose_z_to_y(sgsphi3, sgsphi2)
    call transpose_y_to_x(sgsphi2, sgsphi1)

    ! SGS correction for ABL
    if (itype.eq.itype_abl.and.is==1.and.ibuoyancy.eq.1) then
       call transpose_y_to_x(tb2,dphidy1)
       call wall_sgs_slip_scalar(sgsphi1,nut1,dphidy1)
    endif

  end subroutine sgs_scalar_nonconservative

  !************************************************************
  subroutine sgs_mom_conservative(sgsx1,sgsy1,sgsz1,ux1,uy1,uz1,nut1)

    USE param
    USE variables
    use MPI
    USE var, only : ta1,tb1,tc1,di1
    use var, only : td1, te1, tf1, tg1, th1, ti1
    USE var, only : ta2,tb2,tc2,di2
    use var, only :      te2, tf2, tg2, th2, ti2
    USE var, only : ta3,tb3,tc3,di3
    use var, only :           tf3,      th3, ti3
    USE var, only : sgsx2,sgsy2,sgsz2
    USE var, only : sgsx3,sgsy3,sgsz3
    USE var, only : sxx1,sxy1,sxz1,syy1,syz1,szz1
    USE abl, only : wall_sgs_noslip
    use ibm_param, only: ubcx, ubcy, ubcz
    
    implicit none 

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)),intent(in) :: ux1, uy1, uz1, nut1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)),intent(out) :: sgsx1, sgsy1, sgsz1

    !                                                         ta1        tb1        tc1 
    !real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: wallsgsx1, wallsgsy1, wallsgsz1
    !
    !                                                         td1  te1  tf1  tg1  th1  ti1         
    !real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: txx, txy, txz, tyy, tyz, tzz 
    
    integer :: i, j, k, code, ierr
    
    ! Construct stress tensor
    td1 = 2.0*nut1*sxx1
    te1 = 2.0*nut1*sxy1
    tf1 = 2.0*nut1*sxz1
    tg1 = 2.0*nut1*syy1
    th1 = 2.0*nut1*syz1
    ti1 = 2.0*nut1*szz1   
    
    ! Add wall model for ABL
    if (itype.eq.itype_abl) then
      if (ncly1==2) then 
        call wall_sgs_noslip(ux1,uy1,uz1,nut1,ta1,tb1,tc1)
        if (xstart(2)==1) then
          td1(:,2,:) = 0.
          te1(:,2,:) = - ta1(:,2,:)! txy1(:,2,:)
          tf1(:,2,:) = 0.
          tg1(:,2,:) = 0.
          th1(:,2,:) = - tc1(:,2,:)! tyz1(:,2,:)
          ti1(:,2,:) = 0.
        endif
      elseif (ncly1==1) then 
        write(*,*) 'Simulation stopped: slip bottom wall not supported with iconserv=1'
        call MPI_ABORT(MPI_COMM_WORLD,code,ierr); stop 
      endif
    endif
    
    ! Compute derivatives
    ta1 = zero; ta2 = zero; ta3 = zero
    tb1 = zero; tb2 = zero; tb3 = zero
    tc1 = zero; tc2 = zero; tc3 = zero
    sgsx1=zero; sgsy1=zero; sgsz1=zero
    sgsx2=zero; sgsy2=zero; sgsz2=zero
    sgsx3=zero; sgsy3=zero; sgsz3=zero

    ! WORK X-PENCILS
    call derx (sgsx1,td1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,zero)
    call derx (sgsy1,te1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,zero)
    call derx (sgsz1,tf1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,zero)

    !call filter(0.48d0)
    !call filx(taf1,ta1,di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1,zero)
    !call filx(tbf1,tb1,di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1,zero)
    !call filx(tcf1,tc1,di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1,zero)

    ! Stored into var immediatly not here
    !sgsx1 = ta1
    !sgsy1 = tb1
    !sgsz1 = tc1

    ! WORK Y-PENCILS
    call transpose_x_to_y(te1, te2)
    call transpose_x_to_y(tg1, tg2)
    call transpose_x_to_y(th1, th2)
    call transpose_x_to_y(tf1, tf2)
    call transpose_x_to_y(ti1, ti2)
    call transpose_x_to_y(sgsx1, sgsx2)
    call transpose_x_to_y(sgsy1, sgsy2)
    call transpose_x_to_y(sgsz1, sgsz2)

    call dery (ta2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,zero)
    call dery (tb2,tg2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,zero)
    call dery (tc2,th2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,zero)

    !call fily(taf2,ta2,di2,fisy,fiffyp,fifsyp,fifwyp,ysize(1),ysize(2),ysize(3),1,zero)
    !call fily(tbf2,tb2,di2,fisy,fiffyp,fifsyp,fifwyp,ysize(1),ysize(2),ysize(3),1,zero)
    !call fily(tcf2,tc2,di2,fisy,fiffyp,fifsyp,fifwyp,ysize(1),ysize(2),ysize(3),1,zero)

    sgsx2 = sgsx2 + ta2
    sgsy2 = sgsy2 + tb2
    sgsz2 = sgsz2 + tc2

    ! WORK Z-PENCILS
    call transpose_y_to_z(sgsx2, sgsx3)
    call transpose_y_to_z(sgsy2, sgsy3)
    call transpose_y_to_z(sgsz2, sgsz3)
    call transpose_y_to_z(tf2, tf3)
    call transpose_y_to_z(th2, th3)
    call transpose_y_to_z(ti2, ti3)

    call derz (ta3, tf3, di3, sz, ffz, fsz, fwz, zsize(1), zsize(2), zsize(3), 0, zero)
    call derz (tb3, th3, di3, sz, ffz, fsz, fwz, zsize(1), zsize(2), zsize(3), 0, zero)
    call derz (tc3, ti3, di3, sz, ffz, fsz, fwz, zsize(1), zsize(2), zsize(3), 0, zero)

    !call filz(taf3,ta3,di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1,zero)
    !call filz(tbf3,tb3,di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1,zero)
    !call filz(tcf3,tc3,di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1,zero)

    sgsx3 = sgsx3 + ta3
    sgsy3 = sgsy3 + tb3
    sgsz3 = sgsz3 + tc3

    call transpose_z_to_y(sgsx3, sgsx2)
    call transpose_z_to_y(sgsy3, sgsy2)
    call transpose_z_to_y(sgsz3, sgsz2)

    call transpose_y_to_x(sgsx2, sgsx1)
    call transpose_y_to_x(sgsy2, sgsy1)
    call transpose_y_to_x(sgsz2, sgsz1)

  end subroutine sgs_mom_conservative

end module les
