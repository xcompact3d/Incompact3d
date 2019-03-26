MODULE transeq

  PRIVATE
  PUBLIC :: momentum_rhs_eq, continuity_rhs_eq, scalar
  
CONTAINS

  subroutine momentum_rhs_eq(dux1,duy1,duz1,rho1,ux1,uy1,uz1,ep1,phi1,divu3)

    USE param
    USE variables
    USE decomp_2d
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : rho2,ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
    USE var, only : rho3,ux3,uy3,uz3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
    USE var, only : sgsx1,sgsy1,sgsz1

    USE case, ONLY : momentum_forcing

    implicit none

    !! INPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    real(mytype),intent(in),dimension(zsize(1),zsize(2),zsize(3)) :: divu3

    !! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1

    integer :: ijk,nvect1,nvect2,nvect3,i,j,k,is

    !SKEW SYMMETRIC FORM
    !WORK X-PENCILS
    ta1(:,:,:) = rho1(:,:,:,1) * ux1(:,:,:) * ux1(:,:,:)
    tb1(:,:,:) = rho1(:,:,:,1) * ux1(:,:,:) * uy1(:,:,:)
    tc1(:,:,:) = rho1(:,:,:,1) * ux1(:,:,:) * uz1(:,:,:)

    call derx (td1,ta1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
    call derx (te1,tb1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
    call derx (tf1,tc1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)

    ta1(:,:,:) = td1(:,:,:) + rho1(:,:,:,1) * ux1(:,:,:) * ta1(:,:,:)
    tb1(:,:,:) = te1(:,:,:) + rho1(:,:,:,1) * ux1(:,:,:) * tb1(:,:,:)
    tc1(:,:,:) = tf1(:,:,:) + rho1(:,:,:,1) * ux1(:,:,:) * tc1(:,:,:)

    if (ilmn) then
       !! Quasi-skew symmetric terms
       call derx (td1,rho1(:,:,:,1),di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
       ta1(:,:,:) = ta1(:,:,:) + ux1(:,:,:) * ux1(:,:,:) * td1(:,:,:)
       tb1(:,:,:) = tb1(:,:,:) + uy1(:,:,:) * ux1(:,:,:) * td1(:,:,:)
       tc1(:,:,:) = tc1(:,:,:) + uz1(:,:,:) * ux1(:,:,:) * td1(:,:,:)
    endif

    call transpose_x_to_y(ux1,ux2)
    call transpose_x_to_y(uy1,uy2)
    call transpose_x_to_y(uz1,uz2)
    call transpose_x_to_y(ta1,ta2)
    call transpose_x_to_y(tb1,tb2)
    call transpose_x_to_y(tc1,tc2)

    if (ilmn) then
       call transpose_x_to_y(rho1(:,:,:,1),rho2)
    else
       rho2(:,:,:) = one
    endif
    
    !WORK Y-PENCILS
    td2(:,:,:) = rho2(:,:,:) * ux2(:,:,:) * uy2(:,:,:)
    te2(:,:,:) = rho2(:,:,:) * uy2(:,:,:) * uy2(:,:,:)
    tf2(:,:,:) = rho2(:,:,:) * uz2(:,:,:) * uy2(:,:,:)

    call dery (tg2,td2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
    call dery (th2,te2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    call dery (ti2,tf2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
    call dery (td2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    call dery (te2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
    call dery (tf2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)

    ta2(:,:,:) = ta2(:,:,:) + (tg2(:,:,:) + rho2(:,:,:) * uy2(:,:,:) * td2(:,:,:))
    tb2(:,:,:) = tb2(:,:,:) + (th2(:,:,:) + rho2(:,:,:) * uy2(:,:,:) * te2(:,:,:))
    tc2(:,:,:) = tc2(:,:,:) + (ti2(:,:,:) + rho2(:,:,:) * uy2(:,:,:) * tf2(:,:,:))

    if (ilmn) then
       !! Quasi-skew symmetric terms
       call dery (te2,rho2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
       ta2(:,:,:) = ta2(:,:,:) + ux2(:,:,:) * uy2(:,:,:) * te2(:,:,:)
       tb2(:,:,:) = tb2(:,:,:) + uy2(:,:,:) * uy2(:,:,:) * te2(:,:,:)
       tc2(:,:,:) = tc2(:,:,:) + uz2(:,:,:) * uy2(:,:,:) * te2(:,:,:)
    endif

    call transpose_y_to_z(ux2,ux3)
    call transpose_y_to_z(uy2,uy3)
    call transpose_y_to_z(uz2,uz3)
    call transpose_y_to_z(ta2,ta3)
    call transpose_y_to_z(tb2,tb3)
    call transpose_y_to_z(tc2,tc3)

    if (ilmn) then
       call transpose_y_to_z(rho2,rho3)
    else
       rho3(:,:,:) = one
    endif

    !WORK Z-PENCILS
    td3(:,:,:) = rho3(:,:,:) * ux3(:,:,:) * uz3(:,:,:)
    te3(:,:,:) = rho3(:,:,:) * uy3(:,:,:) * uz3(:,:,:)
    tf3(:,:,:) = rho3(:,:,:) * uz3(:,:,:) * uz3(:,:,:)

    call derz (tg3,td3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
    call derz (th3,te3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
    call derz (ti3,tf3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
    call derz (td3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
    call derz (te3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
    call derz (tf3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)

    ta3(:,:,:) = ta3(:,:,:) + (tg3(:,:,:) + rho3(:,:,:) * uz3(:,:,:) * td3(:,:,:))
    tb3(:,:,:) = tb3(:,:,:) + (th3(:,:,:) + rho3(:,:,:) * uz3(:,:,:) * te3(:,:,:))
    tc3(:,:,:) = tc3(:,:,:) + (ti3(:,:,:) + rho3(:,:,:) * uz3(:,:,:) * tf3(:,:,:))

    if (ilmn) then
       !! Quasi-skew symmetric terms
       call derz (tf3,rho3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
       ta3(:,:,:) = ta3(:,:,:) + ux3(:,:,:) * uz3(:,:,:) * tf3(:,:,:)
       tb3(:,:,:) = tb3(:,:,:) + uy3(:,:,:) * uz3(:,:,:) * tf3(:,:,:)
       tc3(:,:,:) = tc3(:,:,:) + uz3(:,:,:) * uz3(:,:,:) * tf3(:,:,:)

       !! Add the additional divu terms
       ta3(:,:,:) = ta3(:,:,:) + rho3(:,:,:) * ux3(:,:,:) * divu3(:,:,:)
       tb3(:,:,:) = tb3(:,:,:) + rho3(:,:,:) * uy3(:,:,:) * divu3(:,:,:)
       tc3(:,:,:) = tc3(:,:,:) + rho3(:,:,:) * uz3(:,:,:) * divu3(:,:,:)
    endif

    !! Skew symmetric - need to multiply by half
    ta3(:,:,:) = half * ta3(:,:,:)
    tb3(:,:,:) = half * tb3(:,:,:)
    tc3(:,:,:) = half * tc3(:,:,:)

    !ALL THE CONVECTIVE TERMS ARE IN TA3, TB3 and TC3
    td3 = ta3 
    te3 = tb3 
    tf3 = tc3 

    !DIFFUSIVE TERMS IN Z
    call derzz (ta3,ux3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
    call derzz (tb3,uy3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
    call derzz (tc3,uz3,di3,sz,sfz ,ssz ,swz ,zsize(1),zsize(2),zsize(3),0)

    !WORK Y-PENCILS
    call transpose_z_to_y(ta3,ta2)
    call transpose_z_to_y(tb3,tb2)
    call transpose_z_to_y(tc3,tc2)
    call transpose_z_to_y(td3,td2)
    call transpose_z_to_y(te3,te2)
    call transpose_z_to_y(tf3,tf2)

    tg2 = td2 
    th2 = te2 
    ti2 = tf2 

    !DIFFUSIVE TERMS IN Y
    !-->for ux
    call deryy (td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
    if (istret.ne.0) then
       call dery (te2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
       do k = 1,ysize(3)
          do j = 1,ysize(2)
             do i = 1,ysize(1)
                td2(i,j,k) = td2(i,j,k)*pp2y(j)-pp4y(j)*te2(i,j,k)
             enddo
          enddo
       enddo
    endif

    !-->for uy
    call deryy (te2,uy2,di2,sy,sfy,ssy,swy,ysize(1),ysize(2),ysize(3),0)
    if (istret.ne.0) then
       call dery (tf2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
       do k = 1,ysize(3)
          do j = 1,ysize(2)
             do i = 1,ysize(1)
                te2(i,j,k) = te2(i,j,k)*pp2y(j)-pp4y(j)*tf2(i,j,k)
             enddo
          enddo
       enddo
    endif

    !-->for uz
    call deryy (tf2,uz2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
    if (istret.ne.0) then
       call dery (tj2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
       do k = 1,ysize(3)
          do j = 1,ysize(2)
             do i = 1,ysize(1)
                tf2(i,j,k) = tf2(i,j,k)*pp2y(j)-pp4y(j)*tj2(i,j,k)
             enddo
          enddo
       enddo
    endif

    ta2 = ta2 + td2 
    tb2 = tb2 + te2 
    tc2 = tc2 + tf2 

    !WORK X-PENCILS
    call transpose_y_to_x(ta2,ta1)
    call transpose_y_to_x(tb2,tb1)
    call transpose_y_to_x(tc2,tc1) !diff
    call transpose_y_to_x(tg2,td1)
    call transpose_y_to_x(th2,te1)
    call transpose_y_to_x(ti2,tf1) !conv

    tg1 = td1 
    th1 = te1 
    ti1 = tf1 

    !DIFFUSIVE TERMS IN X
    call derxx (td1,ux1,di1,sx,sfx ,ssx ,swx ,xsize(1),xsize(2),xsize(3),0)
    call derxx (te1,uy1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
    call derxx (tf1,uz1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)

    ta1 = ta1 + td1 
    tb1 = tb1 + te1 
    tc1 = tc1 + tf1

    ! di1 =  zero
    ! do is = 1, numscalar
    !    di1 = di1  + phi1(:,:,:,is)*ri(is) !Mod. by Ricardo
    ! enddo

    !FINAL SUM: DIFF TERMS + CONV TERMS
    dux1(:,:,:,1) = xnu*ta1(:,:,:) - tg1(:,:,:) ! + di1(:,:,:)*anglex  !+x
    duy1(:,:,:,1) = xnu*tb1(:,:,:) - th1(:,:,:) ! - di1(:,:,:)*angley  !+y
    duz1(:,:,:,1) = xnu*tc1(:,:,:) - ti1(:,:,:) ! !+- di1       !+z

    if (ilmn) then
       call momentum_full_viscstress_tensor(dux1(:,:,:,1), duy1(:,:,:,1), duz1(:,:,:,1), divu3)
    endif

    dux1(:,:,:,1) = xnu*ta1(:,:,:) - tg1(:,:,:) + di1(:,:,:)*anglex  !+x
    duy1(:,:,:,1) = xnu*tb1(:,:,:) - th1(:,:,:) - di1(:,:,:)*angley  !+y
    duz1(:,:,:,1) = xnu*tc1(:,:,:) - ti1(:,:,:) !+- di1       !+z


    ! If LES modelling is enabled, add the SGS stresses
    if (ilesmod.ne.0.and.jles.le.3.) then
        ! Wall model for LES
        if (iwall.eq.1) then 
        call compute_SGS(sgsx1,sgsy1,sgsz1,ux1,uy1,uz1,ep1,1)
        else
        call compute_SGS(sgsx1,sgsy1,sgsz1,ux1,uy1,uz1,ep1,0)
        endif
        ! Calculate SGS stresses (conservative/non-conservative formulation)
        dux1(:,:,:,1) = dux1(:,:,:,1) + sgsx1(:,:,:)
        duy1(:,:,:,1) = duy1(:,:,:,1) + sgsy1(:,:,:)
        duz1(:,:,:,1) = duz1(:,:,:,1) + sgsz1(:,:,:)
    endif

    call momentum_forcing(dux1, duy1, duz1, rho1, ux1, uy1, uz1)

    if (itrip == 1) then
       call tripping(tb1,td1)
       if (nrank == 0) print *,'TRIPPING KTH STYLE!!'
    endif

  end subroutine momentum_rhs_eq
  !************************************************************

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  SUBROUTINE: momentum_full_viscstress_tensor
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: In an incompressible flow the viscous stress
  !!              tensor reduces to
  !!                d2u^j / {dx^i}^2
  !!              however if div(u) != 0 we have
  !!                d2u^j / {dx^i}^2 + 1/3 d/dx^j div(u)
  !!              and further if \mu != const.
  !!                \mu (d2u^j / {dx^i}^2 + 1/3 d/dx^j div(u))
  !!                  + d\mu/dx^i (du^j/dx^i + du^i/dx^j
  !!                  - 2/3 div(u) \delta^{ij})
  !!              This subroutine computes the additional
  !!              contributions not accounted for in the
  !!              incompressible solver.
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine momentum_full_viscstress_tensor(ta1, tb1, tc1, divu3)

    USE param
    USE variables
    USE decomp_2d
    USE var, only : td1,te1,tf1,tg1,di1
    USE var, only : ta2,tb2,tc2,td2,di2
    USE var, only : tc3,di3

    IMPLICIT NONE

    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ta1, tb1, tc1
    REAL(mytype), DIMENSION(zsize(1), zsize(2), zsize(3)), INTENT(IN) :: divu3

    INTEGER :: i, j, k
    REAL(mytype) :: one_third

    one_third = one / three
    
    call derz (tc3,divu3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
    call transpose_z_to_y(tc3, tc2)
    call transpose_z_to_y(divu3, ta2)

    call dery(tb2,ta2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
    call transpose_y_to_x(tb2, te1)
    call transpose_y_to_x(tc2, tf1)
    call transpose_y_to_x(ta2, tg1)

    call derx(td1,tg1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)

    ta1(:,:,:) = ta1(:,:,:) + one_third * xnu * td1(:,:,:)
    tb1(:,:,:) = tb1(:,:,:) + one_third * xnu * te1(:,:,:)
    tc1(:,:,:) = tc1(:,:,:) + one_third * xnu * tf1(:,:,:)

  end subroutine momentum_full_viscstress_tensor
  
  !************************************************************
  subroutine scalar(dphi1, ux1, uy1, uz1, phi1)

    USE param
    USE variables
    USE decomp_2d
    
    USE var, ONLY : ta1,tb1,tc1,td1,di1
    USE var, ONLY : uy2,phi2,ta2,tb2,tc2,td2,di2
    USE var, ONLY : uz3,phi3,ta3,tb3,di3

    implicit none

    !! INPUTS
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    !! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime,numscalar) :: dphi1

    !! LOCALS
    integer :: i, j, k, is

    !!=====================================================================
    !! XXX It is assumed that ux,uy,uz are already updated in all pencils!
    !!=====================================================================
    do is = 1, numscalar

       !X PENCILS
       call derxS (tb1,phi1(:,:,:,is),di1,sx,ffxpS,fsxpS,fwxpS,xsize(1),xsize(2),xsize(3),1)
       tb1(:,:,:) = ux1(:,:,:) * tb1(:,:,:)
       call derxxS (ta1,phi1(:,:,:,is),di1,sx,sfxpS,ssxpS,swxpS,xsize(1),xsize(2),xsize(3),1)
       call transpose_x_to_y(phi1(:,:,:,is),phi2(:,:,:,is))

       !Y PENCILS
       call deryS (tb2,phi2(:,:,:,is),di2,sy,ffypS,fsypS,fwypS,ppy,ysize(1),ysize(2),ysize(3),1)
       tb2(:,:,:) = uy2(:,:,:) * tb2(:,:,:)
       call deryyS (ta2,phi2(:,:,:,is),di2,sy,sfypS,ssypS,swypS,ysize(1),ysize(2),ysize(3),1)
       if (istret.ne.0) then
          call deryS (tc2,phi2(:,:,:,is),di2,sy,ffypS,fsypS,fwypS,ppy,ysize(1),ysize(2),ysize(3),1)
          do k = 1,ysize(3)
             do j = 1,ysize(2)
                do i = 1,ysize(1)
                   ta2(i,j,k) = ta2(i,j,k)*pp2y(j)-pp4y(j)*tc2(i,j,k)
                enddo
             enddo
          enddo
       endif
       call transpose_y_to_z(phi2(:,:,:,is),phi3(:,:,:,is))

       !Z PENCILS
       call derzS (tb3,phi3(:,:,:,is),di3,sz,ffzpS,fszpS,fwzpS,zsize(1),zsize(2),zsize(3),1)
       tb3 = tb3*uz3
       call derzzS (ta3,phi3(:,:,:,is),di3,sz,sfzpS,sszpS,swzpS,zsize(1),zsize(2),zsize(3),1)

       call transpose_z_to_y(ta3,tc2)
       call transpose_z_to_y(tb3,td2)

       !Y PENCILS ADD TERMS
       tc2 = tc2+ta2
       td2 = td2+tb2

       call transpose_y_to_x(tc2,tc1)
       call transpose_y_to_x(td2,td1)

       !X PENCILS ADD TERMS
       ta1 = ta1+tc1 !SECOND DERIVATIVE
       tb1 = tb1+td1 !FIRST DERIVATIVE

       dphi1(:,:,:,1,is) = (xnu/sc(is))*ta1(:,:,:) - tb1(:,:,:)

    end do !loop numscalar

  end subroutine scalar

  SUBROUTINE continuity_rhs_eq(drho1, rho1, ux1, uy1, uz1, divu3)

    USE decomp_2d, ONLY : mytype, xsize, ysize, zsize
    USE decomp_2d, ONLY : transpose_z_to_y, transpose_y_to_x
    USE param, ONLY : ntime, nrhotime
    USE variables

    USE var, ONLY : ta1, di1
    USE var, ONLY : rho2, uy2, ta2, tb2, di2
    USE var, ONLY : rho3, uz3, ta3, di3
    
    IMPLICIT NONE

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), INTENT(IN), DIMENSION(zsize(1), zsize(2), zsize(3)) :: divu3
    
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime) :: drho1

    !! XXX All variables up to date - no need to transpose

    CALL derz (ta3, rho3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1)
    ta3(:,:,:) = uz3(:,:,:) * ta3(:,:,:) + rho3(:,:,:) * divu3(:,:,:)
    
    CALL transpose_z_to_y(ta3, tb2)
    CALL dery (ta2, rho2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1)
    ta2(:,:,:) = uy2(:,:,:) * ta2(:,:,:) + tb2(:,:,:)

    CALL transpose_y_to_x(ta2, ta1)
    CALL derx (drho1(:,:,:,1), rho1(:,:,:,1), &
         di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1)
    drho1(:,:,:,1) = -(ux1(:,:,:) * drho1(:,:,:,1) + ta1(:,:,:))

  ENDSUBROUTINE continuity_rhs_eq
  
END MODULE transeq
