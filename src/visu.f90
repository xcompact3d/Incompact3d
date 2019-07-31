module visu

  implicit none

  private
  public :: write_snapshot, postprocessing

contains

  SUBROUTINE postprocessing(rho1, ux1, uy1, uz1, pp3, phi1, ep1)

    USE decomp_2d, ONLY : mytype, xsize, ph1
    USE case, ONLY : postprocess_case

    USE var, ONLY : nzmsize
    USE var, ONLY : itime
    USE var, ONLY : numscalar, nrhotime, npress
    
    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3)), intent(in) :: ux1, uy1, uz1
    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3),numscalar), intent(in) :: phi1
    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3),nrhotime), intent(in) :: rho1
    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3)), intent(in) :: ep1
    real(mytype), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(in) :: pp3

    CALL write_snapshot(rho1, ux1, uy1, uz1, pp3(:,:,:,1), phi1, ep1, itime)
    CALL postprocess_case(rho1, ux1, uy1, uz1, pp3, phi1, ep1)
    
  END SUBROUTINE postprocessing

  SUBROUTINE write_snapshot(rho1, ux1, uy1, uz1, pp3, phi1, ep1, itime)

    USE decomp_2d, ONLY : transpose_x_to_y, transpose_y_to_z, transpose_z_to_y, transpose_y_to_x
    USE decomp_2d, ONLY : mytype, xsize, ysize, zsize
    USE decomp_2d, ONLY : fine_to_coarseV
    USE decomp_2d_io, ONLY : decomp_2d_write_one

    USE param, ONLY : ivisu, ioutput, nrhotime, ilmn, iscalar, iibm

    USE variables, ONLY : derx, dery, derz 
    USE variables, ONLY : ffx, ffxp, fsx, fsxp, fwx, fwxp
    USE variables, ONLY : ffy, ffyp, fsy, fsyp, fwy, fwyp, ppy
    USE variables, ONLY : ffz, ffzp, fsz, fszp, fwz, fwzp
    USE variables, ONLY : sx, cifip6, cisip6, ciwip6, cifx6, cisx6, ciwx6
    USE variables, ONLY : sy, cifip6y, cisip6y, ciwip6y, cify6, cisy6, ciwy6
    USE variables, ONLY : sz, cifip6z, cisip6z, ciwip6z, cifz6, cisz6, ciwz6
    USE variables, ONLY : numscalar

    USE var, ONLY : one
    USE var, ONLY : uvisu
    USE var, ONLY : pp1, ta1, tb1, tc1, td1, te1, tf1, tg1, th1, ti1, di1, nxmsize
    USE var, ONLY : pp2, ta2, tb2, tc2, td2, te2, tf2, ppi2, di2, dip2, ph2, nymsize
    USE var, ONLY : ppi3, ta3, tb3, tc3, td3, te3, tf3, di3, dip3, ph3, nzmsize

    IMPLICIT NONE

    CHARACTER(len=30) :: filename

    !! Inputs
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ux1, uy1, uz1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ep1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime), INTENT(IN) :: rho1
    REAL(mytype), DIMENSION(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize), INTENT(IN) :: pp3
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), numscalar), INTENT(IN) :: phi1
    INTEGER, INTENT(IN) :: itime

    INTEGER :: is

    IF ((ivisu.NE.0).AND.(MOD(itime, ioutput).EQ.0)) THEN
       !! Write velocity
       uvisu=0.
       if (iibm==2) then
          ta1(:,:,:) = (one - ep1(:,:,:)) * ux1(:,:,:)
       else
          ta1(:,:,:) = ux1(:,:,:)
       endif
       call fine_to_coarseV(1,ta1,uvisu)
990    format('ux',I3.3)
       write(filename, 990) itime/ioutput
       call decomp_2d_write_one(1,uvisu,filename,2)

       uvisu=0.
       if (iibm==2) then
          ta1(:,:,:) = (one - ep1(:,:,:)) * uy1(:,:,:)
       else
          ta1(:,:,:) = uy1(:,:,:)
       endif
       call fine_to_coarseV(1,ta1,uvisu)
991    format('uy',I3.3)
       write(filename, 991) itime/ioutput
       call decomp_2d_write_one(1,uvisu,filename,2)

       uvisu=0.
       if (iibm==2) then
          ta1(:,:,:) = (one - ep1(:,:,:)) * uz1(:,:,:)
       else
          ta1(:,:,:) = uz1(:,:,:)
       endif
       call fine_to_coarseV(1,ta1,uvisu)
992    format('uz',I3.3)
       write(filename, 992) itime/ioutput
       call decomp_2d_write_one(1,uvisu,filename,2)

       !! Write pressure
       !WORK Z-PENCILS
       call interzpv(ppi3,pp3,dip3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
            (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
       !WORK Y-PENCILS
       call transpose_z_to_y(ppi3,pp2,ph3) !nxm nym nz
       call interypv(ppi2,pp2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
            (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
       !WORK X-PENCILS
       call transpose_y_to_x(ppi2,pp1,ph2) !nxm ny nz
       call interxpv(ta1,pp1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
            nxmsize,xsize(1),xsize(2),xsize(3),1)

       uvisu=0._mytype
       if (iibm==2) then
          ta1(:,:,:) = (one - ep1(:,:,:)) * ta1(:,:,:)
       endif
       call fine_to_coarseV(1,ta1,uvisu)
993    format('pp',I3.3)
       write(filename, 993) itime/ioutput
       call decomp_2d_write_one(1,uvisu,filename,2)

       !! Write vorticity
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
994    format('vort',I3.3)
       write(filename, 994) itime/ioutput
       call decomp_2d_write_one(1,uvisu,filename,2)

       !! LMN - write out density
       IF (ilmn) THEN
          uvisu=0.
          call fine_to_coarseV(1,rho1(:,:,:,1),uvisu)
995       format('rho',I3.3)
          write(filename, 995) itime/ioutput
          call decomp_2d_write_one(1,uvisu,filename,2)
       ENDIF

       !! Scalars
       IF (iscalar.NE.0) THEN
996       format('phi',I1.1,I3.3)
          DO is = 1, numscalar
             uvisu=0.
             call fine_to_coarseV(1,phi1(:,:,:,is),uvisu)
             write(filename, 996) is, itime/ioutput
             call decomp_2d_write_one(1,uvisu,filename,2)
          ENDDO
       ENDIF
    ENDIF
  END SUBROUTINE write_snapshot

endmodule visu

!############################################################################
subroutine VISU_INSTA (ux1,uy1,uz1,phi1,ep1,protection)

  USE param
  USE variables
  USE decomp_2d
  USE decomp_2d_io
  USE MPI

  implicit none

  real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
  real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: temp1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: temp2
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: temp3
#ifdef VISUEXTRA
  real(mytype),dimension(ysize(1),ysize(2),ysize(3),numscalar) :: phi2
  real(mytype),dimension(zsize(1),zsize(2),zsize(3),numscalar) :: phi3

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,phim1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,phim2,di2
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phim3
#endif
  real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu

  real(mytype) :: s3df,s2df,ttsize
  real(8) :: t1,tstart,tend
  integer :: n3df,n2df,i,j,k,is
  logical,intent(in) :: protection
  character(len=1) :: a
  character(len=30) :: filename

  !Modificado por Felipe Schuch
  open(10,file='visu.prm',status='unknown',form='formatted')
  read (10,*) a !
  read (10,*) a ! INCOMPACT 3D parameters - Visualization - requires compilation flag -DVISU
  read (10,*) a !
  read (10,*) save_ux
  read (10,*) save_uy
  read (10,*) save_uz
  read (10,*) save_phi
  read (10,*) save_uxm
  read (10,*) save_uym
  read (10,*) save_uzm
  read (10,*) save_phim
  read (10,*) save_pre
  read (10,*) save_prem
  read (10,*) save_ibm
#ifdef VISUEXTRA
  read (10,*) a !
  read (10,*) a ! Extra Visualization - requires compilation flag -DVISUEXTRA
  read (10,*) a !
  read (10,*) save_w
  read (10,*) save_w1
  read (10,*) save_w2
  read (10,*) save_w3
  read (10,*) save_qc
  read (10,*) save_pc
  read (10,*) save_V
  read (10,*) save_dudx
  read (10,*) save_dudy
  read (10,*) save_dudz
  read (10,*) save_dvdx
  read (10,*) save_dvdy
  read (10,*) save_dvdz
  read (10,*) save_dwdx
  read (10,*) save_dwdy
  read (10,*) save_dwdz
  read (10,*) save_dphidx
  read (10,*) save_dphidy
  read (10,*) save_dphidz
  read (10,*) save_dmap
  read (10,*) save_utmap
#endif
  close(10)

  if (nrank .eq. 0) call paraview()

  if (protection) then !Files that must be protected from rewriting
     !when visu runs from post-processing
     save_ux=0; save_uy=0; save_uz=0; save_ibm=0; save_phi=0; save_pre=0
  endif

  !number of 3D fields
  n3df=save_w+save_w1+save_w2+save_w3+save_qc+save_pc+save_ux+save_uy+save_uz+&
       save_phi*numscalar+save_pre+save_dudx+save_dudy+save_dudz+&
       save_dvdx+save_dvdy+save_dvdz+save_dwdx+save_dwdy+save_dwdz+save_V+&
       (save_dphidx+save_dphidy+save_dphidz)*numscalar

  !number of 2D fields
  n2df=save_uxm+save_uym+save_uzm+save_phim*numscalar+save_prem+save_dmap+save_utmap

  s2df=nx*ny*prec   !size of a single 2D field - value in Byte-B
  s3df=s2df*nz      !size of a single 3D field - value in Byte-B
  ttsize=(n3df*s3df+n2df*s2df)   !total size for 2 and 3D fields
  if (save_ibm .eq. 2) ttsize=ttsize+s3df

  tstart=MPI_WTIME()

  if (nrank.eq.0) print *,'==========================================================='
  if (nrank.eq.0) print *,'Writing snapshot =>',itime/ioutput
  if (nrank.eq.0) print *,'This simulation requieres',real((ttsize*int(ilast/ioutput)&
       +(16+3*numscalar)*s3df*int(ilast/icheckpoint))*1e-9,4),'GB' !1e-9 from Byte-B to Gigabyte-GB

  if (save_ux.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,ux1,uvisu)
     write(filename,"('./data/ux',I4.4)") itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_uxm.eq.1) then
     call transpose_x_to_y (ux1,temp2)
     call transpose_y_to_z (temp2,temp3)
     call mean_plane_z(temp3,zsize(1),zsize(2),zsize(3),temp3(:,:,1))
     write(filename,"('./data/uxm',I4.4)") itime/ioutput
     call decomp_2d_write_plane(3,temp3,3,1,filename)
  endif

  if (save_uy.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,uy1,uvisu)
     write(filename,"('./data/uy',I4.4)") itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_uym.eq.1) then
     call transpose_x_to_y (uy1,temp2)
     call transpose_y_to_z (temp2,temp3)
     call mean_plane_z(temp3,zsize(1),zsize(2),zsize(3),temp3(:,:,1))
     write(filename,"('./data/uym',I4.4)") itime/ioutput
     call decomp_2d_write_plane(3,temp3,3,1,filename)
  endif

  if (save_uz.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,uz1,uvisu)
     write(filename,"('./data/uz',I4.4)") itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_uzm.eq.1) then
     call transpose_x_to_y (uz1,temp2)
     call transpose_y_to_z (temp2,temp3)
     call mean_plane_z(temp3,zsize(1),zsize(2),zsize(3),temp3(:,:,1))
     write(filename,"('./data/uzm',I4.4)") itime/ioutput
     call decomp_2d_write_plane(3,temp3,3,1,filename)
  endif

  if (save_V.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,sqrt(ux1**2+uy1**2+uz1**2),uvisu)
     write(filename,"('./data/V',I4.4)") itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (iibm.ne.0) then
     if (save_ibm.ne.0) then
        uvisu=0._mytype
        call fine_to_coarseV(1,ep1,uvisu)
        if (save_ibm.eq.1) write(filename,"('./data/ibm',I4.4)") 0
        if (save_ibm.eq.2) write(filename,"('./data/ibm',I4.4)") itime/ioutput
        call decomp_2d_write_one(1,uvisu,filename,2)
     endif
  endif

  if (iscalar.eq.1) then
     if (save_phi.eq.1) then
        do is=1, numscalar
           uvisu=0._mytype
           call fine_to_coarseV(1,phi1(:,:,:,is),uvisu)
           write(filename,"('./data/phi',I1.1,I4.4)") is, itime/ioutput
           call decomp_2d_write_one(1,uvisu,filename,2)
        enddo
     endif

     if (save_phim.eq.1) then
        do is=1, numscalar
           call transpose_x_to_y (phi1(:,:,:,is),temp2)
           call transpose_y_to_z (temp2,temp3)
           call mean_plane_z(temp3,zsize(1),zsize(2),zsize(3),temp3(:,:,1))
           write(filename,"('./data/phim',I1.1,I4.4)") is, itime/ioutput
           call decomp_2d_write_plane(3,temp3,3,1,filename)
        enddo
     endif
  endif

#ifdef VISUEXTRA
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

  if (save_w.eq.1) then
     di1(:,:,:) = sqrt((tf1(:,:,:) - th1(:,:,:))**2 &
          + (tg1(:,:,:) - tc1(:,:,:))**2 + (tb1(:,:,:) - td1(:,:,:))**2)
     uvisu=0._mytype
     call fine_to_coarseV(1,di1,uvisu)
     write(filename,"('./data/w',I4.4)") itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_w1.eq.1) then
     !dw/dy - dv/dz
     di1(:,:,:) = tf1(:,:,:) - th1(:,:,:)
     uvisu=0.
     call fine_to_coarseV(1,di1,uvisu)
     write(filename,"('./data/w1',I4.4)") itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_w2.eq.1) then
     !du/dz - dw/dx
     di1(:,:,:) = tg1(:,:,:) - tc1(:,:,:)
     uvisu=0._mytype
     call fine_to_coarseV(1,di1,uvisu)
     write(filename,"('./data/w2',I4.4)") itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_w3.eq.1) then
     !dv/dx - du/dy
     di1(:,:,:) = th1(:,:,:) - td1(:,:,:)
     uvisu=0._mytype
     call fine_to_coarseV(1,di1,uvisu)
     write(filename,"('./data/w3',I4.4)") itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_dudx.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,ta1,uvisu)
     write(filename,"('./data/dudx',I4.4)") itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_dudy.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,td1,uvisu)
     write(filename,"('./data/dudy',I4.4)") itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_dudz.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,tg1,uvisu)
     write(filename,"('./data/dudz',I4.4)") itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_dvdx.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,tb1,uvisu)
     write(filename,"('./data/dvdx',I4.4)") itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_dvdy.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,te1,uvisu)
     write(filename,"('./data/dvdy',I4.4)") itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_dvdz.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,th1,uvisu)
     write(filename,"('./data/dvdz',I4.4)") itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_dwdx.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,tc1,uvisu)
     write(filename,"('./data/dwdx',I4.4)") itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_dwdy.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,tf1,uvisu)
     write(filename,"('./data/dwdy',I4.4)") itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_dwdz.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,ti1,uvisu)
     write(filename,"('./data/dwdz',I4.4)") itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_pc.eq.1) then
     !P=-(1./3.)*(ta1**3+te1**3+ti1**3+2*td1*th1*tc1+2*tb1*tf1*tg1+
     !tb1*tg1*tf1+tc1*td1*th1)-ta1*td1*tb1-ta1*tg1*tc1-
     !td1*te1*tb1-tg1*ti1*tc1-te1*th1*tf1-th1*ti1*tf1
     di1(:,:,:)=-(1._mytype/3._mytype)*(ta1(:,:,:)**3+&
          te1(:,:,:)**3+ti1(:,:,:)**3+&
          2.*td1(:,:,:)*th1(:,:,:)*tc1(:,:,:)+&
          2.*tb1(:,:,:)*tf1(:,:,:)*tg1(:,:,:)+&
          tb1(:,:,:)*tg1(:,:,:)*tf1(:,:,:)+&
          tc1(:,:,:)*td1(:,:,:)*th1(:,:,:))-&
          ta1(:,:,:)*td1(:,:,:)*tb1(:,:,:)-&
          ta1(:,:,:)*tg1(:,:,:)*tc1(:,:,:)-&
          td1(:,:,:)*te1(:,:,:)*tb1(:,:,:)-&
          tg1(:,:,:)*ti1(:,:,:)*tc1(:,:,:)-&
          te1(:,:,:)*th1(:,:,:)*tf1(:,:,:)-&
          th1(:,:,:)*ti1(:,:,:)*tf1(:,:,:)
     uvisu=0._mytype
     call fine_to_coarseV(1,di1,uvisu)
     write(filename,"('./data/pc',I4.4)") itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_qc.eq.1) then
     di1=0._mytype
     !Q=-0.5*(ta1**2+te1**2+ti1**2)-td1*tb1-tg1*tc1-th1*tf1
     di1(:,:,:)=-0.5*(ta1(:,:,:)**2+te1(:,:,:)**2+ti1(:,:,:)**2)-&
          td1(:,:,:)*tb1(:,:,:)-tg1(:,:,:)*tc1(:,:,:)-th1(:,:,:)*tf1(:,:,:)
     uvisu=0._mytype
     call fine_to_coarseV(1,di1,uvisu)
     write(filename,"('./data/qc',I4.4)") itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_utmap.eq.1) then
     di1=0._mytype
     write(filename,"('./data/utmap',I4.4)") itime/ioutput
     di1(:,:,:)=sqrt(sqrt((td1(:,:,:)**2)+(tf1(:,:,:)**2))*xnu)
     call decomp_2d_write_plane(1,di1,2,1,filename)
  endif

  if (iscalar.eq.1) then
     if (save_dphidx.eq.1) then
        do is=1, numscalar
           call derxS (tb1,phi1(:,:,:,is),di1,sx,ffxpS,fsxpS,fwxpS,xsize(1),xsize(2),xsize(3),1) !phi even
           write(filename,"('./data/dphi',I1.1,'dx',I4.4)") is, itime/ioutput
           call decomp_2d_write_one(1,tb1,filename,2)
        enddo
     endif

     if (save_dphidy.eq.1) then
        do is=1, numscalar
           call transpose_x_to_y(phi1(:,:,:,is),phi2(:,:,:,is))
           call deryS (tb2,phi2(:,:,:,is),di2,sy,ffypS,fsypS,fwypS,ppy,ysize(1),ysize(2),ysize(3),1) !phi even
           write(filename,"('./data/dphi',I1.1,'dy',I4.4)") is, itime/ioutput
           call decomp_2d_write_one(2,tb2,filename,2)
        enddo
     endif

     if (save_dphidz.eq.1) then
        do is=1, numscalar
           call transpose_x_to_y(phi1(:,:,:,is),phi2(:,:,:,is))
           call transpose_y_to_z(phi2(:,:,:,is),phi3(:,:,:,is))
           call derzS (tb3,phi3(:,:,:,is),di3,sz,ffzpS,fszpS,fwzpS,zsize(1),zsize(2),zsize(3),1) !phi even
           write(filename,"('./data/dphi',I1.1,'dz',I4.4)") is, itime/ioutput
           call decomp_2d_write_one(3,tb3,filename,2)
        enddo
     endif
  endif

#endif

  tend=MPI_WTIME()
  if (nrank.eq.0) print *,'Time in visu (s)', real(tend-tstart,4)
  if (nrank.eq.0) print *,'Writing speed (MB/s)',real((ttsize*1e-6)/(tend-tstart),4)

end subroutine VISU_INSTA
!############################################################################
subroutine VISU_PRE (pp3,ta1,tb1,di1,ta2,tb2,di2,ta3,di3,nxmsize,nymsize,nzmsize,uvisu,pre1)

  USE param
  USE variables
  USE decomp_2d
  USE decomp_2d_io

  implicit none

  integer :: nxmsize,nymsize,nzmsize

  real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu
  real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize) :: pp3
  !Z PENCILS NXM NYM NZM-->NXM NYM NZ
  real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)) :: ta3,di3
  !Y PENCILS NXM NYM NZ -->NXM NY NZ
  real(mytype),dimension(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)) :: ta2
  real(mytype),dimension(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)) :: tb2,di2
  !X PENCILS NXM NY NZ  -->NX NY NZ
  real(mytype),dimension(nxmsize,xsize(2),xsize(3)) :: ta1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tb1,di1,pre1

  character(len=30) filename

  !WORK Z-PENCILS
  call interzpv(ta3,pp3,di3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
       (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
  !WORK Y-PENCILS
  call transpose_z_to_y(ta3,ta2,ph3) !nxm nym nz
  call interypv(tb2,ta2,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
       (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
  !WORK X-PENCILS
  call transpose_y_to_x(tb2,ta1,ph2) !nxm ny nz
  call interxpv(tb1,ta1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
       nxmsize,xsize(1),xsize(2),xsize(3),1)

  pre1=tb1

  if (save_pre.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,pre1,uvisu)
     write(filename,"('./data/pre',I4.4)") itime/ioutput
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_prem.eq.1) then
     tb1=0._mytype
     call mean_plane_z(pre1,xsize(1),xsize(2),xsize(3),tb1(:,:,1))
     write(filename,"('./data/prem',I4.4)") itime/ioutput
     call decomp_2d_write_plane(1,tb1,3,1,filename)
  endif

  return

end subroutine VISU_PRE

!######################################################################################
subroutine mean_plane_x (f1,nx,ny,nz,fm1)

  use param, only : mytype, zero

  implicit none

  integer,intent(in) :: nx, ny, nz
  real(mytype),intent(in),dimension(nx,ny,nz) :: f1
  real(mytype),intent(out),dimension(ny,nz) :: fm1
  integer :: i,j,k

  fm1 = sum(f1,DIM=1)/real(nx,mytype)
  return

end subroutine mean_plane_x
!!######################################################################################
subroutine mean_plane_y (f2,nx,ny,nz,fm2)

  use param, only : mytype, zero

  implicit none

  integer,intent(in) :: nx, ny, nz
  real(mytype),intent(in),dimension(nx,ny,nz) :: f2
  real(mytype),intent(out),dimension(nx,nz) :: fm2
  integer :: i,j,k

  fm2 = sum(f2,DIM=2)/real(ny,mytype)
  return

end subroutine mean_plane_y
!######################################################################################
subroutine mean_plane_z (f3,nx,ny,nz,fm3)

  use param, only : mytype, zero

  implicit none

  integer,intent(in) :: nx, ny, nz
  real(mytype),intent(in),dimension(nx,ny,nz) :: f3
  real(mytype),intent(out),dimension(nx,ny) :: fm3
  integer :: i,j,k

  fm3 = sum(f3,DIM=3)/real(nz,mytype)
  return

end subroutine mean_plane_z
