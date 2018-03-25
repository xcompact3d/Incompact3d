subroutine convdiff(ux1,uy1,uz1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
     ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
     ux3,uy3,uz3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phi1,ep1,nut1)

  USE param
  USE variables
  USE decomp_2d

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,sumphi,nut1,ep1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ux3,uy3,uz3
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3

  integer :: ijk,nvect1,nvect2,nvect3,i,j,k,is

#ifdef ELES
  !############################## EXPLICIT LES MODELLING #######

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gxx1,gyx1,gzx1,gxy1,gyy1,gzy1,gxz1,gyz1,gzz1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: asxx1,asyy1,aszz1,asxy1,asxz1,asyz1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: sxx1,syy1,szz1,sxy1,sxz1,syz1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: srt_smag,sgsx1,sgsy1,sgsz1

  sgsx1=zero; sgsy1=zero; sgsz1=zero; nut1=zero; srt_smag=zero


     if (jLES==2) then !SMAG

        call smag(ux1,uy1,uz1,gxx1,gyx1,gzx1,gxy1,gyy1,gzy1,gxz1,gyz1,gzz1,&
             sxx1,syy1,szz1,sxy1,sxz1,syz1,srt_smag,nut1,ta2,ta3,di1,di2,di3)

     elseif (jLES == 3) then !WALE

        call smag(ux1,uy1,uz1,gxx1,gyx1,gzx1,gxy1,gyy1,gzy1,gxz1,gyz1,gzz1,&
             sxx1,syy1,szz1,sxy1,sxz1,syz1,srt_smag,nut1,ta2,ta3,di1,di2,di3)

        call wale(gxx1,gyx1,gzx1,gxy1,gyy1,gzy1,gxz1,gyz1,gzz1,srt_smag,nut1)

     elseif (jLES==4) then !DSMAG

        call smag(ux1,uy1,uz1,gxx1,gyx1,gzx1,gxy1,gyy1,gzy1,gxz1,gyz1,gzz1,&
             sxx1,syy1,szz1,sxy1,sxz1,syz1,srt_smag,nut1,ta2,ta3,di1,di2,di3)

        call dynsmag(ux1,uy1,uz1,ep1,sxx1,syy1,szz1,sxy1,sxz1,syz1,&
             srt_smag,nut1,di1,ta1,tb1,tc1,td1,ta2,tb2,tc2,td2,te2,tf2,&
             tg2,th2,ti2,di2,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3)
     endif

     call lesdiff(ux1,uy1,uz1,sxx1,syy1,szz1,sxy1,sxz1,syz1,nut1,&
          sgsx1,sgsy1,sgsz1,ta1,td1,te1,tf1,di1,ta2,td2,te2,tf2,tj2,di2,&
          ta3,td3,te3,tf3,di3)


  ta1 = zero; tb1 = zero; tc1 = zero
  td1 = zero; te1 = zero; tf1 = zero

  ta2 = zero; tb2 = zero; tc2 = zero
  td2 = zero; te2 = zero; tf2 = zero
  tj2 = zero

  ta3 = zero; tb3 = zero; tc3 = zero
  td3 = zero; te3 = zero; tf3 = zero

  !########################## ENDING EXPLIICIT LES TERMS ##################################
#endif

  nvect1 = xsize(1)*xsize(2)*xsize(3)
  nvect2 = ysize(1)*ysize(2)*ysize(3)
  nvect3 = zsize(1)*zsize(2)*zsize(3)

  !  !u.rot(u)
  !  !WORK X-PENCILS
  !   call derx (ta1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  !   call derx (tb1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  !   call transpose_x_to_y(ux1,ux2)
  !   call transpose_x_to_y(uy1,uy2)
  !   call transpose_x_to_y(uz1,uz2)
  !   call transpose_x_to_y(ta1,ta2)
  !   call transpose_x_to_y(tb1,tb2)
  !  !WORK Y-PENCILS
  !   call dery (tc2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) 
  !   call dery (td2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) 
  !   call transpose_y_to_z(ux2,ux3)
  !   call transpose_y_to_z(uy2,uy3)
  !   call transpose_y_to_z(uz2,uz3)
  !   call transpose_y_to_z(ta2,ta3)
  !   call transpose_y_to_z(tb2,tb3)
  !   call transpose_y_to_z(tc2,tc3)
  !   call transpose_y_to_z(td2,td3)
  !  !WORK Z-PENCILS
  !   call derz (te3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
  !   call derz (tf3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
  !   do ijk = 1,nvect3
  !      ta3(ijk,1,1) = uz3(ijk,1,1)*(te3(ijk,1,1)-tb3(ijk,1,1))-&
  !           uy3(ijk,1,1)*(ta3(ijk,1,1)-tc3(ijk,1,1))
  !      tb3(ijk,1,1) = ux3(ijk,1,1)*(ta3(ijk,1,1)-tc3(ijk,1,1))-&
  !           uz3(ijk,1,1)*(td3(ijk,1,1)-tf3(ijk,1,1))
  !      tc3(ijk,1,1) = uy3(ijk,1,1)*(td3(ijk,1,1)-tf3(ijk,1,1))-&
  !           ux3(ijk,1,1)*(te3(ijk,1,1)-tb3(ijk,1,1))
  !   enddo
  !SKEW!
  !WORK X-PENCILS
  do ijk = 1,nvect1
     ta1(ijk,1,1) = ux1(ijk,1,1)*ux1(ijk,1,1)
     tb1(ijk,1,1) = ux1(ijk,1,1)*uy1(ijk,1,1)
     tc1(ijk,1,1) = ux1(ijk,1,1)*uz1(ijk,1,1)
  enddo

  call derx (td1,ta1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  call derx (te1,tb1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
  call derx (tf1,tc1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
  call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
  call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)

  do ijk = 1,nvect1
     ta1(ijk,1,1) = half*td1(ijk,1,1)+half*ux1(ijk,1,1)*ta1(ijk,1,1)
     tb1(ijk,1,1) = half*te1(ijk,1,1)+half*ux1(ijk,1,1)*tb1(ijk,1,1)
     tc1(ijk,1,1) = half*tf1(ijk,1,1)+half*ux1(ijk,1,1)*tc1(ijk,1,1)
  enddo

  call transpose_x_to_y(ux1,ux2)
  call transpose_x_to_y(uy1,uy2)
  call transpose_x_to_y(uz1,uz2)
  call transpose_x_to_y(ta1,ta2)
  call transpose_x_to_y(tb1,tb2)
  call transpose_x_to_y(tc1,tc2)

  !WORK Y-PENCILS
  do ijk = 1,nvect2
     td2(ijk,1,1) = ux2(ijk,1,1)*uy2(ijk,1,1)
     te2(ijk,1,1) = uy2(ijk,1,1)*uy2(ijk,1,1)
     tf2(ijk,1,1) = uz2(ijk,1,1)*uy2(ijk,1,1)
  enddo

  call dery (tg2,td2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
  call dery (th2,te2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
  call dery (ti2,tf2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
  call dery (td2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
  call dery (te2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
  call dery (tf2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)

  do ijk = 1,nvect2
     ta2(ijk,1,1) = ta2(ijk,1,1)+half*tg2(ijk,1,1)+half*uy2(ijk,1,1)*td2(ijk,1,1)
     tb2(ijk,1,1) = tb2(ijk,1,1)+half*th2(ijk,1,1)+half*uy2(ijk,1,1)*te2(ijk,1,1)
     tc2(ijk,1,1) = tc2(ijk,1,1)+half*ti2(ijk,1,1)+half*uy2(ijk,1,1)*tf2(ijk,1,1)
  enddo
  call transpose_y_to_z(ux2,ux3)
  call transpose_y_to_z(uy2,uy3)
  call transpose_y_to_z(uz2,uz3)
  call transpose_y_to_z(ta2,ta3)
  call transpose_y_to_z(tb2,tb3)
  call transpose_y_to_z(tc2,tc3)

  !WORK Z-PENCILS
  do ijk = 1,nvect3
     td3(ijk,1,1) = ux3(ijk,1,1)*uz3(ijk,1,1)
     te3(ijk,1,1) = uy3(ijk,1,1)*uz3(ijk,1,1)
     tf3(ijk,1,1) = uz3(ijk,1,1)*uz3(ijk,1,1)
  enddo

  call derz (tg3,td3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
  call derz (th3,te3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
  call derz (ti3,tf3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
  call derz (td3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
  call derz (te3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
  call derz (tf3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)

  do ijk = 1,nvect3
     ta3(ijk,1,1) = ta3(ijk,1,1)+half*tg3(ijk,1,1)+half*uz3(ijk,1,1)*td3(ijk,1,1)
     tb3(ijk,1,1) = tb3(ijk,1,1)+half*th3(ijk,1,1)+half*uz3(ijk,1,1)*te3(ijk,1,1)
     tc3(ijk,1,1) = tc3(ijk,1,1)+half*ti3(ijk,1,1)+half*uz3(ijk,1,1)*tf3(ijk,1,1)
  enddo

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
  if (istret.ne.0) then
     call deryy (td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
     call dery (te2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
     do k = 1,ysize(3)
        do j = 1,ysize(2)
           do i = 1,ysize(1)
              td2(i,j,k) = td2(i,j,k)*pp2y(j)-pp4y(j)*te2(i,j,k)
           enddo
        enddo
     enddo
  else
     call deryy (td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
  endif

  !-->for uy
  if (istret.ne.0) then
     call deryy (te2,uy2,di2,sy,sfy,ssy,swy,ysize(1),ysize(2),ysize(3),0)
     call dery (tf2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
     do k = 1,ysize(3)
        do j = 1,ysize(2)
           do i = 1,ysize(1)
              te2(i,j,k) = te2(i,j,k)*pp2y(j)-pp4y(j)*tf2(i,j,k)
           enddo
        enddo
     enddo
  else
     call deryy (te2,uy2,di2,sy,sfy,ssy,swy,ysize(1),ysize(2),ysize(3),0)
  endif

  !-->for uz
  if (istret.ne.0) then
     call deryy (tf2,uz2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
     call dery (tj2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
     do k = 1,ysize(3)
        do j = 1,ysize(2)
           do i = 1,ysize(1)
              tf2(i,j,k) = tf2(i,j,k)*pp2y(j)-pp4y(j)*tj2(i,j,k)
           enddo
        enddo
     enddo
  else
     call deryy (tf2,uz2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
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

  sumphi =  zero

  do is = 1, nphi
     sumphi = sumphi  + phi1(:,:,:,is)*ri(is) !Mod. by Ricardo
  enddo

  !FINAL SUM: DIFF TERMS + CONV TERMS
  ta1 = xnu*ta1 - tg1 + sumphi*anglex  !+x
  tb1 = xnu*tb1 - th1 - sumphi*angley  !+y
  tc1 = xnu*tc1 - ti1 !+- sumphi       !+z

  if (itime.lt.irotation) then
     if (nrank==0) print *,'Rotating turbulent channel!'
     ta1 = ta1 - wrotation*uy1
     tb1 = tb1 + wrotation*ux1
  endif

#ifdef ELES
  ta1 = ta1 + sgsx1 
  tb1 = tb1 + sgsy1 
  tc1 = tc1 + sgsz1 
#endif

  if (itrip == 1) then
     call tripping(tb1,td1)
     if (nrank == 0) print *,'TRIPPING KTH STYLE!!'
  endif

end subroutine convdiff
!************************************************************
subroutine scalar(ux1,uy1,uz1,phi1,phis1,phiss1,di1,ta1,tb1,tc1,td1,&
     uy2,uz2,phi2,di2,ta2,tb2,tc2,td2,uz3,phi3,di3,ta3,tb3,epsi,nut1)

  USE param
  USE variables
  USE decomp_2d

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,di1,ta1,tb1,tc1,td1,epsi,nut1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1,phis1,phiss1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: uy2,uz2,di2,ta2,tb2,tc2,td2
  real(mytype),dimension(ysize(1),ysize(2),ysize(3),nphi) :: phi2
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: uz3,di3,ta3,tb3
  real(mytype),dimension(zsize(1),zsize(2),zsize(3),nphi) :: phi3
  integer :: ijk,nvect1,nvect2,nvect3,i,j,k,is

#ifdef ELES
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: kappat1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: sgsphi1
  kappat1 = nut1 / pr_t
  call lesdiff_scalar(phi1, di1, di2, di3, kappat1, sgsphi1)
#endif

  nvect1 = xsize(1)*xsize(2)*xsize(3)
  nvect2 = ysize(1)*ysize(2)*ysize(3)
  nvect3 = zsize(1)*zsize(2)*zsize(3)

  do is = 1, nphi

     !X PENCILS
     call derxS (tb1,phi1(:,:,:,is),di1,sx,ffxpS,fsxpS,fwxpS,xsize(1),xsize(2),xsize(3),1)
     do ijk = 1,nvect1
        tb1(ijk,1,1) = tb1(ijk,1,1)*(ux1(ijk,1,1)+uset(is)*anglex)
     enddo

     call derxxS (ta1,phi1(:,:,:,is),di1,sx,sfxpS,ssxpS,swxpS,xsize(1),xsize(2),xsize(3),1)

     call transpose_x_to_y(phi1(:,:,:,is),phi2(:,:,:,is))

     !Y PENCILS
     call deryS (tb2,phi2(:,:,:,is),di2,sy,ffypS,fsypS,fwypS,ppy,ysize(1),ysize(2),ysize(3),1)
     do ijk = 1,nvect2
       tb2(ijk,1,1) = tb2(ijk,1,1)*(uy2(ijk,1,1)-uset(is)*angley)
     enddo

     if (istret.ne.0) then
        call deryyS (ta2,phi2(:,:,:,is),di2,sy,sfypS,ssypS,swypS,ysize(1),ysize(2),ysize(3),1)
        call deryS (tc2,phi2(:,:,:,is),di2,sy,ffypS,fsypS,fwypS,ppy,ysize(1),ysize(2),ysize(3),1)
        do k = 1,ysize(3)
           do j = 1,ysize(2)
              do i = 1,ysize(1)
                 ta2(i,j,k) = ta2(i,j,k)*pp2y(j)-pp4y(j)*tc2(i,j,k)
              enddo
           enddo
        enddo
     else
        call deryyS (ta2,phi2(:,:,:,is),di2,sy,sfypS,ssypS,swypS,ysize(1),ysize(2),ysize(3),1)
     endif

     call transpose_y_to_z(phi2(:,:,:,is),phi3(:,:,:,is))

     !Z PENCILS
     call derzS (tb3,phi3(:,:,:,is),di3,sz,ffzpS,fszpS,fwzpS,zsize(1),zsize(2),zsize(3),1)
     do ijk = 1,nvect3
        tb3(ijk,1,1) = tb3(ijk,1,1)*uz3(ijk,1,1)
     enddo

     call derzzS (ta3,phi3(:,:,:,is),di3,sz,sfzpS,sszpS,swzpS,zsize(1),zsize(2),zsize(3),1)

     call transpose_z_to_y(ta3,tc2)
     call transpose_z_to_y(tb3,td2)

     !Y PENCILS ADD TERMS
     do ijk = 1,nvect2
        tc2(ijk,1,1) = tc2(ijk,1,1)+ta2(ijk,1,1)
        td2(ijk,1,1) = td2(ijk,1,1)+tb2(ijk,1,1)
     enddo

     call transpose_y_to_x(tc2,tc1)
     call transpose_y_to_x(td2,td1)

     !X PENCILS ADD TERMS
     do ijk = 1,nvect1
        ta1(ijk,1,1) = ta1(ijk,1,1)+tc1(ijk,1,1) !SECOND DERIVATIVE
        tb1(ijk,1,1) = tb1(ijk,1,1)+td1(ijk,1,1) !FIRST DERIVATIVE
     enddo

#ifdef ELES
  if (nrank == 0) print *, "sgsphi",is,"min max= ",minval(sgsphi1(:,:,:,is)),maxval(sgsphi1(:,:,:,is))
     do ijk=1,nvect1
        ta1(ijk,1,1)=(xnu/nsc(is) + kappat1(ijk,1,1) )*ta1(ijk,1,1)-tb1(ijk,1,1)+sgsphi1(ijk,1,1,is)
     enddo
#else
     do ijk=1,nvect1
        ta1(ijk,1,1)=(xnu/nsc(is))*ta1(ijk,1,1)-tb1(ijk,1,1)
     enddo
#endif

     !TIME ADVANCEMENT
     if ((nscheme.eq.1).or.(nscheme.eq.3)) then !AB2 or RK3
        if ((nscheme.eq.1.and.itime.eq.1.and.ilit.eq.0).or.&
             (nscheme.eq.3.and.itr.eq.1)) then
           do ijk = 1,nvect1
              phi1(ijk,1,1,is) = gdt(itr)*ta1(ijk,1,1)+phi1(ijk,1,1,is)
              phis1(ijk,1,1,is) = ta1(ijk,1,1)
           enddo
        else
           do ijk = 1,nvect1
              phi1(ijk,1,1,is) = adt(itr)*ta1(ijk,1,1)+bdt(itr)*phis1(ijk,1,1,is)+phi1(ijk,1,1,is)
              phis1(ijk,1,1,is) = ta1(ijk,1,1)
           enddo
        endif
     endif

     if (nscheme.eq.2) then !AB3
        if ((itime.eq.1).and.(ilit.eq.0)) then
           if (nrank.eq.0) print *,'start with Euler for scalar',itime
           do ijk = 1,nvect1 !start with Euler
              phi1(ijk,1,1,is) = dt*ta1(ijk,1,1)+phi1(ijk,1,1,is)
              phis1(ijk,1,1,is) = ta1(ijk,1,1)
           enddo
        else
           if  ((itime.eq.2).and.(ilit.eq.0)) then
              if (nrank.eq.0) print *,'then with AB2 for scalar',itime
              do ijk = 1,nvect1
                 phi1(ijk,1,1,is) = onepfive*dt*ta1(ijk,1,1)-half*dt*phis1(ijk,1,1,is)+phi1(ijk,1,1,is)
                 phiss1(ijk,1,1,is) = phis1(ijk,1,1,is)
                 phis1(ijk,1,1,is) = ta1(ijk,1,1)
              enddo
           else
              do ijk = 1,nvect1
                 phi1(ijk,1,1,is) = adt(itr)*ta1(ijk,1,1)+bdt(itr)*phis1(ijk,1,1,is)+&
                      cdt(itr)*phiss1(ijk,1,1,is)+phi1(ijk,1,1,is)
                 phiss1(ijk,1,1,is) = phis1(ijk,1,1,is)
                 phis1(ijk,1,1,is) = ta1(ijk,1,1)
              enddo
           endif
        endif
     endif

     if (cont_phi.eq.1) then  !Controle de phi. Ativa no prm.
        do ijk = 1,nvect1
           if (phi1(ijk,1,1,is).lt.zero) phi1(ijk,1,1,is) =  zero
        enddo
     endif

     if (cont_phi.eq.2) then  !Controle de phi. Ativa no prm.
        do ijk = 1,nvect1
           if (phi1(ijk,1,1,is).gt.cp(is)) phi1(ijk,1,1,is) =  cp(is)
           if (phi1(ijk,1,1,is).lt.zero) phi1(ijk,1,1,is) =  zero
        enddo
     endif

  end do !loop nphi

end subroutine scalar
