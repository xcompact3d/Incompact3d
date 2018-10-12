subroutine OVERALL_STATISTIC(ux1,uy1,uz1,phi1,pre1,diss1,ta1,&
     u1sum,v1sum,w1sum,u2sum,v2sum,w2sum,&
     u3sum,v3sum,w3sum,u4sum,v4sum,w4sum,&
     uvsum,uwsum,vwsum,disssum,tsum,&
     psum,ppsum,upsum,vpsum,wpsum)

  USE param
  USE variables
  USE decomp_2d
  USE decomp_2d_io

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,pre1,diss1,ta1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1 

  real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: u1sum,v1sum,w1sum !first
  real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: u2sum,v2sum,w2sum !second
  real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: u3sum,v3sum,w3sum !third
  real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: u4sum,v4sum,w4sum !fourth
  real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: uvsum,uwsum,vwsum !tensors
  real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: disssum,presum,tsum !diss, pressure and auxiliary variable

  real(mytype),dimension(xszS(1),xszS(2),xszS(3),nphi) :: psum,ppsum
  real(mytype),dimension(xszS(1),xszS(2),xszS(3),nphi) :: upsum,vpsum,wpsum
  integer :: is


  !pressure
  call fine_to_coarseS(1,pre1,tsum)
  presum=presum+tsum

  !u1=ux1
  call fine_to_coarseS(1,ux1,tsum)
  u1sum=u1sum+tsum

  !u2=ux1*ux1
  ta1=ux1*ux1
  call fine_to_coarseS(1,ta1,tsum)
  u2sum=u2sum+tsum

  !u3=ux1*ux1*ux1 (THIRD ORDER MOMENTS - SKEWNESS/assimetria)
  ta1=ux1*ux1*ux1
  call fine_to_coarseS(1,ta1,tsum)
  u3sum=u3sum+tsum

  !u4=ux1*ux1*ux1*ux1 (FOURTH ORDER MOMENTS - FLATNESS/achatamento)
  ta1=ux1*ux1*ux1*ux1
  call fine_to_coarseS(1,ta1,tsum)
  u4sum=u4sum+tsum

  !v1=uy1
  call fine_to_coarseS(1,uy1,tsum)
  v1sum=v1sum+tsum

  !v2=uy1*uy1
  ta1=uy1*uy1
  call fine_to_coarseS(1,ta1,tsum)
  v2sum=v2sum+tsum

  !v3sum=uy1*uy1*uy1 (THIRD ORDER MOMENTS - SKEWNESS/assimetria)
  ta1=uy1*uy1*uy1
  call fine_to_coarseS(1,ta1,tsum)
  v3sum=v3sum+tsum

  !v4sum=uy1*uy1*uy1*uy1 (FOURTH ORDER MOMENTS - FLATNESS/achatamento)
  ta1=uy1*uy1*uy1*uy1
  call fine_to_coarseS(1,ta1,tsum)
  v4sum=v4sum+tsum

  !w1=uz1
  call fine_to_coarseS(1,uz1,tsum)
  w1sum=w1sum+tsum

  !w2=uz1*uz1
  ta1=uz1*uz1
  call fine_to_coarseS(1,ta1,tsum)
  w2sum=w2sum+tsum

  !w3=uz1*uz1*uz1 (THIRD ORDER MOMENTS - SKEWNESS/assimetria)
  ta1=uz1*uz1*uz1
  call fine_to_coarseS(1,ta1,tsum)
  w3sum=w3sum+tsum

  !w4=uz1*uz1*uz1*uz1 (FOURTH ORDER MOMENTS - FLATNESS/achatamento)
  ta1=uz1*uz1*uz1*uz1
  call fine_to_coarseS(1,ta1,tsum)
  w4sum=w4sum+tsum

  !uvsum=ux1*uy1
  ta1=ux1*uy1
  call fine_to_coarseS(1,ta1,tsum)
  uvsum=uvsum+tsum

  !uwsum=ux1*uz1
  ta1=ux1*uz1
  call fine_to_coarseS(1,ta1,tsum)
  uwsum=uwsum+tsum

  !vwsum=uy1*uz1
  ta1=uy1*uz1
  call fine_to_coarseS(1,ta1,tsum)
  vwsum=vwsum+tsum

  if (iscalar==1) then
     do is=1, nphi

        !psum=phi1
        call fine_to_coarseS(1,phi1(:,:,:,is),tsum)
        psum(:,:,:,is)=psum(:,:,:,is)+tsum

        !ppsum=phi1*phi1
        ta1=phi1(:,:,:,is)*phi1(:,:,:,is)
        call fine_to_coarseS(1,ta1,tsum)
        ppsum(:,:,:,is)=ppsum(:,:,:,is)+tsum

        !upsum=phi1*ux1
        ta1=phi1(:,:,:,is)*ux1
        call fine_to_coarseS(1,ta1,tsum)
        upsum(:,:,:,is)=upsum(:,:,:,is)+tsum

        !vpsum=phi1*uy1
        ta1=phi1(:,:,:,is)*uy1
        call fine_to_coarseS(1,ta1,tsum)
        vpsum(:,:,:,is)=vpsum(:,:,:,is)+tsum

        !wpsum=phi1*uz1
        ta1=phi1(:,:,:,is)*uz1
        call fine_to_coarseS(1,ta1,tsum)
        wpsum(:,:,:,is)=wpsum(:,:,:,is)+tsum

     enddo
  endif

  if (mod(itime,imodulo)==0) then

     if (save_pre==1) call decomp_2d_write_one(1,presum,'stats/pre.sum',1)

     if (save_ux==1) call decomp_2d_write_one(1,u1sum,'stats/u1.sum',1)
     if (save_ux==1) call decomp_2d_write_one(1,u2sum,'stats/u2.sum',1)
     if (save_ux==1) call decomp_2d_write_one(1,u3sum,'stats/u3.sum',1)
     if (save_ux==1) call decomp_2d_write_one(1,u4sum,'stats/u4.sum',1)

     if (save_uy==1) call decomp_2d_write_one(1,v1sum,'stats/v1.sum',1)
     if (save_uy==1) call decomp_2d_write_one(1,v2sum,'stats/v2.sum',1)
     if (save_uy==1) call decomp_2d_write_one(1,v3sum,'stats/v3.sum',1)
     if (save_uy==1) call decomp_2d_write_one(1,v4sum,'stats/v4.sum',1)

     if (save_uz==1) call decomp_2d_write_one(1,w1sum,'stats/w1.sum',1)
     if (save_uz==1) call decomp_2d_write_one(1,w2sum,'stats/w2.sum',1)
     if (save_uz==1) call decomp_2d_write_one(1,w3sum,'stats/w3.sum',1)
     if (save_uz==1) call decomp_2d_write_one(1,w4sum,'stats/w4.sum',1)

     call decomp_2d_write_one(1,half*(u2sum+v2sum+w2sum),'stats/k.sum',1)

     call decomp_2d_write_one(1,uvsum,'stats/uv.sum',1)
     call decomp_2d_write_one(1,uwsum,'stats/uw.sum',1)
     call decomp_2d_write_one(1,vwsum,'stats/vw.sum',1)

     if (save_phi==1) then
        if (iscalar==1) then
           do is=1, nphi
              call decomp_2d_write_one(1, psum(:,:,:,is), 'stats/c'//char(48+is)//'.sum',1)
              call decomp_2d_write_one(1,ppsum(:,:,:,is),'stats/c'//char(48+is)//'c'//char(48+is)//'.sum',1)
              call decomp_2d_write_one(1,upsum(:,:,:,is),'stats/uc'//char(48+is)//'.sum',1)
              call decomp_2d_write_one(1,vpsum(:,:,:,is),'stats/vc'//char(48+is)//'.sum',1)
              call decomp_2d_write_one(1,wpsum(:,:,:,is),'stats/wc'//char(48+is)//'.sum',1)
           enddo
        endif
     endif
  endif

end subroutine OVERALL_STATISTIC
!############################################################################
subroutine CONVERGENCE_STATISTIC(ux1,ep1,u1sum_tik,u1sum_tak,tsum)

  USE param
  USE variables
  USE decomp_2d
  USE decomp_2d_io

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,ep1
  real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: u1sum_tik,u1sum_tak,tsum
  character(len=30) :: filename

#ifdef IBM
  call fine_to_coarseS(1,(1-ep1)*ux1,tsum)
#else
  call fine_to_coarseS(1,ux1,tsum)
#endif

  u1sum_tik=u1sum_tik+tsum

  if (mod(itime,ntik)==0) then
     if (nrank.eq.0) print *,'Saving uxmean tik =>',(itime/(ntik/2)-1)
     write(filename,"('stats/uxmean',I4.4)") (itime/(ntik/2)-1)
     call decomp_2d_write_one(1,u1sum_tik/ntik,filename,1)
     u1sum_tik = zero
  endif

  if (itime.ge.(ntik/2)) then

     u1sum_tak=u1sum_tak+tsum

     if (itime.gt.(ntik/2).AND.mod(itime+ntik/2,ntik)==0) then
        if (nrank.eq.0) print *,'Saving uxmean tak =>',(itime/(ntik/2)-1)
        write(filename,"('stats/uxmean',I4.4)") (itime/(ntik/2) -1)
        call decomp_2d_write_one(1,u1sum_tak/ntik,filename,1)
        u1sum_tak = zero
     endif

  endif

end subroutine CONVERGENCE_STATISTIC
!############################################################################
subroutine CONVERGENCE_STATISTIC2(ux1,ep1,tik1,tik2,tak1,tak2)

  USE param
  USE variables
  USE decomp_2d
  USE decomp_2d_io

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,ep1
  real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: tik1,tik2,tak1,tak2,tsum
  real(mytype) :: rms1
  character(len=30) :: filename

  rms1 = zero

#ifdef IBM
  call fine_to_coarseS(1,(1-ep1)*ux1,tsum)
#else
  call fine_to_coarseS(1,ux1,tsum)
#endif

  tik1 = tik1 + tsum
  tak1 = tak1 + tsum

  if (mod(itime,ntik)==0) then
     tik2 = tik1/ntik
     tik1 = zero
  end if

  if (itime.eq.(ntik/2.)) tak1 = zero

  if ((itime.gt.(ntik/2.)).AND.(mod(itime+ntik/2,ntik)==0)) then
     tak2 = tak1/ntik
     tak1 = zero
  end if

  if (itime.ge.(int((3./2.)*ntik))) then

     if ((mod(itime,ntik)==0).OR.(mod(itime+ntik/2,ntik)==0)) then

        call RMS(tik2,tak2,rms1)

        if (nrank .eq. 0) then

          print *,'RMS=',rms1

          write(filename,"('stats/rms',I8.8)") itime
          open(67,file=trim(filename),status='unknown',form='formatted')
          write(67,"(2E14.6,I14)") t,rms1,itime
          close(67)

	  rms1=zero

        end if
     endif
  end if

end subroutine CONVERGENCE_STATISTIC2
!############################################################################
subroutine RMS(meanA,meanB,rms1)

  USE param
  USE variables
  USE decomp_2d
  USE MPI

  implicit none

  real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: meanA, meanB
  real(mytype) :: rms0,rms1
  character(len=30) :: filename
  integer :: ijk,code

  rms1 = zero
  do ijk=1,xszS(1)*xszS(2)*xszS(3)
     rms1=rms1+(meanB(ijk,1,1)-meanA(ijk,1,1))**2
  enddo

  rms0 = zero
  rms0 = sqrt(rms1/(nx*ny*nz))
  rms1 = zero

  call MPI_REDUCE(rms0,rms1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
  return

end subroutine RMS
!############################################################################
subroutine EXTRA_STAT (ux1,uy1,uz1,uvisu,tsum,dudxsum,utmapsum)

  USE param
  USE variables
  USE decomp_2d
  USE decomp_2d_io

  implicit none

  real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
  real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu,tsum,dudxsum,utmapsum

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,td1,tf1,di1 !ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,phim1,temp1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tc2,td2,tf2 !ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,phim2,di2,temp2
  !real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: !ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phim3,temp3



  !x-derivatives
  call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
  !call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  !call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  !y-derivatives
  call transpose_x_to_y(ux1,td2)
  !call transpose_x_to_y(uy1,te2)
  call transpose_x_to_y(uz1,tf2)
  call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
  !call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
  call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
  !!z-derivatives
  !call transpose_y_to_z(td2,td3)
  !call transpose_y_to_z(te2,te3)
  !call transpose_y_to_z(tf2,tf3)
  !call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
  !call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
  !call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
  !!all back to x-pencils
  !call transpose_z_to_y(ta3,td2)
  !call transpose_z_to_y(tb3,te2)
  !call transpose_z_to_y(tc3,tf2)
  !call transpose_y_to_x(td2,tg1)
  !call transpose_y_to_x(te2,th1)
  !call transpose_y_to_x(tf2,ti1)
  call transpose_y_to_x(ta2,td1)
  !call transpose_y_to_x(tb2,te1)
  call transpose_y_to_x(tc2,tf1)
  !du/dx=ta1 du/dy=td1 and du/dz=tg1
  !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
  !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1


  !FIRST DERIVATIVE
  call fine_to_coarseS(1,ta1,tsum)
  dudxsum = dudxsum + tsum
  if (mod(itime,imodulo)==0) then
     call decomp_2d_write_one(1,dudx,'stats/dudx.sum',1)
  endif

  !FRICTION VELOCITY OVER ITERATION
  di1 = zero
  di1 = sqrt(sqrt( td1**2 + tf1**2 )*xnu) !du/dy**2 + dw/dy**2

  call fine_to_coarseS(1,ta1,tsum)
  utmapsum=utmapsum+tsum
  if (mod(itime,imodulo)==0) then
     call decomp_2d_write_plane(1,utmapsum,2,1,'stats/utmap.sum')
  endif



end subroutine EXTRA_STAT
