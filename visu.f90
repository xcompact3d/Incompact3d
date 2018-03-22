!############################################################################
#ifdef VISU
subroutine VISU_INSTA (ux1,uy1,uz1,phi1,ep1,protection)

  USE param
  USE variables
  USE decomp_2d
  USE decomp_2d_io
  USE MPI

  implicit none

  real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
  real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: temp1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: temp2
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: temp3
#ifdef VISUEXTRA
  real(mytype),dimension(ysize(1),ysize(2),ysize(3),nphi) :: phi2
  real(mytype),dimension(zsize(1),zsize(2),zsize(3),nphi) :: phi3

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,phim1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,phim2,di2
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phim3
#endif
  real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu

  real(mytype) :: s3df,s2df,ttsize
  real(8) :: t1,tstart,tend
  integer :: n3df,n2df,i,j,k,ijk,nvect1,is
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
       save_phi*nphi+save_pre+save_dudx+save_dudy+save_dudz+&
       save_dvdx+save_dvdy+save_dvdz+save_dwdx+save_dwdy+save_dwdz+save_V+&
       (save_dphidx+save_dphidy+save_dphidz)*nphi

  !number of 2D fields
  n2df=save_uxm+save_uym+save_uzm+save_phim*nphi+save_prem+save_dmap+save_utmap

  s2df=nx*ny*prec   !size of a single 2D field - value in Byte-B
  s3df=s2df*nz      !size of a single 3D field - value in Byte-B
  ttsize=(n3df*s3df+n2df*s2df)   !total size for 2 and 3D fields
  if (save_ibm .eq. 2) ttsize=ttsize+s3df

  tstart=MPI_WTIME()

  if (nrank.eq.0) print *,'==========================================================='
  if (nrank.eq.0) print *,'Writing snapshot =>',itime/imodulo
  if (nrank.eq.0) print *,'This simulation requieres',real((ttsize*int(ilast/imodulo)&
       +(16+3*nphi)*s3df*int(ilast/isave))*1e-9,4),'GB' !1e-9 from Byte-B to Gigabyte-GB

  nvect1=xsize(1)*xsize(2)*xsize(3)

  if (save_ux.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,ux1,uvisu)
     write(filename,"('./data/ux',I4.4)") itime/imodulo
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_uxm.eq.1) then
     call transpose_x_to_y (ux1,temp2)
     call transpose_y_to_z (temp2,temp3)
     call mean_plane_z(temp3,zsize(1),zsize(2),zsize(3),temp3(:,:,1))
     write(filename,"('./data/uxm',I4.4)") itime/imodulo
     call decomp_2d_write_plane(3,temp3,3,1,filename)
  endif

  if (save_uy.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,uy1,uvisu)
     write(filename,"('./data/uy',I4.4)") itime/imodulo
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_uym.eq.1) then
     call transpose_x_to_y (uy1,temp2)
     call transpose_y_to_z (temp2,temp3)
     call mean_plane_z(temp3,zsize(1),zsize(2),zsize(3),temp3(:,:,1))
     write(filename,"('./data/uym',I4.4)") itime/imodulo
     call decomp_2d_write_plane(3,temp3,3,1,filename)
  endif

  if (save_uz.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,uz1,uvisu)
     write(filename,"('./data/uz',I4.4)") itime/imodulo
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_uzm.eq.1) then
     call transpose_x_to_y (uz1,temp2)
     call transpose_y_to_z (temp2,temp3)
     call mean_plane_z(temp3,zsize(1),zsize(2),zsize(3),temp3(:,:,1))
     write(filename,"('./data/uzm',I4.4)") itime/imodulo
     call decomp_2d_write_plane(3,temp3,3,1,filename)
  endif

  if (save_V.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,sqrt(ux1**2+uy1**2+uz1**2),uvisu)
     write(filename,"('./data/V',I4.4)") itime/imodulo
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

#ifdef IBM
  if (save_ibm.ne.0) then
     uvisu=0._mytype
     call fine_to_coarseV(1,ep1,uvisu)
     if (save_ibm.eq.1) write(filename,"('./data/ibm',I4.4)") 0
     if (save_ibm.eq.2) write(filename,"('./data/ibm',I4.4)") itime/imodulo
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif
#endif

  if (iscalar.eq.1) then
     if (save_phi.eq.1) then
        do is=1, nphi
           uvisu=0._mytype
           call fine_to_coarseV(1,phi1(:,:,:,is),uvisu)
           write(filename,"('./data/phi',I1.1,I4.4)") is, itime/imodulo
           call decomp_2d_write_one(1,uvisu,filename,2)
        enddo
     endif

     if (save_phim.eq.1) then
        do is=1, nphi
           call transpose_x_to_y (phi1(:,:,:,is),temp2)
           call transpose_y_to_z (temp2,temp3)
           call mean_plane_z(temp3,zsize(1),zsize(2),zsize(3),temp3(:,:,1))
           write(filename,"('./data/phim',I1.1,I4.4)") is, itime/imodulo
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
     di1=0._mytype
     do ijk=1,nvect1
        di1(ijk,1,1)=sqrt((tf1(ijk,1,1)-th1(ijk,1,1))**2+&
             (tg1(ijk,1,1)-tc1(ijk,1,1))**2+(tb1(ijk,1,1)-td1(ijk,1,1))**2)
     enddo
     uvisu=0._mytype
     call fine_to_coarseV(1,di1,uvisu)
     write(filename,"('./data/w',I4.4)") itime/imodulo
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_w1.eq.1) then
     di1=0._mytype
     do ijk=1,nvect1 !dw/dy - dv/dz
        di1(ijk,1,1)=tf1(ijk,1,1)-th1(ijk,1,1)
     enddo
     uvisu=0.
     call fine_to_coarseV(1,di1,uvisu)
     write(filename,"('./data/w1',I4.4)") itime/imodulo
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_w2.eq.1) then
     di1=0._mytype
     do ijk=1,nvect1 !du/dz - dw/dx
        di1(ijk,1,1)=tg1(ijk,1,1)-tc1(ijk,1,1)
     enddo
     uvisu=0._mytype
     call fine_to_coarseV(1,di1,uvisu)
     write(filename,"('./data/w2',I4.4)") itime/imodulo
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_w3.eq.1) then
     di1=0._mytype
     do ijk=1,nvect1 !dv/dx - du/dy
        di1(ijk,1,1)=tb1(ijk,1,1)-td1(ijk,1,1)
     enddo
     uvisu=0._mytype
     call fine_to_coarseV(1,di1,uvisu)
     write(filename,"('./data/w3',I4.4)") itime/imodulo
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_dudx.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,ta1,uvisu)
     write(filename,"('./data/dudx',I4.4)") itime/imodulo
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_dudy.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,td1,uvisu)
     write(filename,"('./data/dudy',I4.4)") itime/imodulo
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_dudz.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,tg1,uvisu)
     write(filename,"('./data/dudz',I4.4)") itime/imodulo
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_dvdx.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,tb1,uvisu)
     write(filename,"('./data/dvdx',I4.4)") itime/imodulo
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_dvdy.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,te1,uvisu)
     write(filename,"('./data/dvdy',I4.4)") itime/imodulo
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_dvdz.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,th1,uvisu)
     write(filename,"('./data/dvdz',I4.4)") itime/imodulo
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_dwdx.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,tc1,uvisu)
     write(filename,"('./data/dwdx',I4.4)") itime/imodulo
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_dwdy.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,tf1,uvisu)
     write(filename,"('./data/dwdy',I4.4)") itime/imodulo
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_dwdz.eq.1) then
     uvisu=0._mytype
     call fine_to_coarseV(1,ti1,uvisu)
     write(filename,"('./data/dwdz',I4.4)") itime/imodulo
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_pc.eq.1) then
     di1=0._mytype
     do ijk=1,nvect1
        !P=-(1./3.)*(ta1**3+te1**3+ti1**3+2*td1*th1*tc1+2*tb1*tf1*tg1+
        !tb1*tg1*tf1+tc1*td1*th1)-ta1*td1*tb1-ta1*tg1*tc1-
        !td1*te1*tb1-tg1*ti1*tc1-te1*th1*tf1-th1*ti1*tf1
        di1(ijk,1,1)=-(1._mytype/3._mytype)*(ta1(ijk,1,1)**3+&
             te1(ijk,1,1)**3+ti1(ijk,1,1)**3+&
             2.*td1(ijk,1,1)*th1(ijk,1,1)*tc1(ijk,1,1)+&
             2.*tb1(ijk,1,1)*tf1(ijk,1,1)*tg1(ijk,1,1)+&
             tb1(ijk,1,1)*tg1(ijk,1,1)*tf1(ijk,1,1)+&
             tc1(ijk,1,1)*td1(ijk,1,1)*th1(ijk,1,1))-&
             ta1(ijk,1,1)*td1(ijk,1,1)*tb1(ijk,1,1)-&
             ta1(ijk,1,1)*tg1(ijk,1,1)*tc1(ijk,1,1)-&
             td1(ijk,1,1)*te1(ijk,1,1)*tb1(ijk,1,1)-&
             tg1(ijk,1,1)*ti1(ijk,1,1)*tc1(ijk,1,1)-&
             te1(ijk,1,1)*th1(ijk,1,1)*tf1(ijk,1,1)-&
             th1(ijk,1,1)*ti1(ijk,1,1)*tf1(ijk,1,1)
     enddo
     uvisu=0._mytype
     call fine_to_coarseV(1,di1,uvisu)
     write(filename,"('./data/pc',I4.4)") itime/imodulo
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_qc.eq.1) then
     di1=0._mytype
     do ijk=1,nvect1
        !Q=-0.5*(ta1**2+te1**2+ti1**2)-td1*tb1-tg1*tc1-th1*tf1
        di1(ijk,1,1)=-0.5*(ta1(ijk,1,1)**2+te1(ijk,1,1)**2+ti1(ijk,1,1)**2)-&
             td1(ijk,1,1)*tb1(ijk,1,1)-tg1(ijk,1,1)*tc1(ijk,1,1)-th1(ijk,1,1)*tf1(ijk,1,1)
     enddo
     uvisu=0._mytype
     call fine_to_coarseV(1,di1,uvisu)
     write(filename,"('./data/qc',I4.4)") itime/imodulo
     call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_utmap.eq.1) then
     di1=0._mytype
     write(filename,"('./data/utmap',I4.4)") itime/imodulo
     do ijk=1,nvect1
        di1(ijk,1,1)=sqrt(sqrt((td1(ijk,1,1)**2)+(tf1(ijk,1,1)**2))*xnu)
     enddo
     call decomp_2d_write_plane(1,di1,2,1,filename)
  endif

  if (iscalar.eq.1) then
     if (save_dphidx.eq.1) then
        do is=1, nphi
           call derxS (tb1,phi1(:,:,:,is),di1,sx,ffxpS,fsxpS,fwxpS,xsize(1),xsize(2),xsize(3),1) !phi even
           write(filename,"('./data/dphi',I1.1,'dx',I4.4)") is, itime/imodulo
           call decomp_2d_write_one(1,tb1,filename,2)
        enddo
     endif

     if (save_dphidy.eq.1) then
        do is=1, nphi
           call transpose_x_to_y(phi1(:,:,:,is),phi2(:,:,:,is))
           call deryS (tb2,phi2(:,:,:,is),di2,sy,ffypS,fsypS,fwypS,ppy,ysize(1),ysize(2),ysize(3),1) !phi even
           write(filename,"('./data/dphi',I1.1,'dy',I4.4)") is, itime/imodulo
           call decomp_2d_write_one(2,tb2,filename,2)
        enddo
     endif

     if (save_dphidz.eq.1) then
        do is=1, nphi
           call transpose_x_to_y(phi1(:,:,:,is),phi2(:,:,:,is))
           call transpose_y_to_z(phi2(:,:,:,is),phi3(:,:,:,is))
           call derzS (tb3,phi3(:,:,:,is),di3,sz,ffzpS,fszpS,fwzpS,zsize(1),zsize(2),zsize(3),1) !phi even
           write(filename,"('./data/dphi',I1.1,'dz',I4.4)") is, itime/imodulo
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
subroutine VISU_PRE (pp3,ta1,tb1,di1,ta2,tb2,di2,ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,uvisu,pre1)

  USE param
  USE variables
  USE decomp_2d
  USE decomp_2d_io

  implicit none

  TYPE(DECOMP_INFO) :: phG,ph2,ph3
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

  integer :: nxmsize,nymsize,nzmsize
  character(len=30) filename

  !WORK Z-PENCILS
  call interiz6(ta3,pp3,di3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
  (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
  !WORK Y-PENCILS
  call transpose_z_to_y(ta3,ta2,ph3) !nxm nym nz
  call interiy6(tb2,ta2,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
  (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
  !WORK X-PENCILS
  call transpose_y_to_x(tb2,ta1,ph2) !nxm ny nz
  call interi6(tb1,ta1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
  nxmsize,xsize(1),xsize(2),xsize(3),1)

  pre1=tb1

  if (save_pre.eq.1) then
    uvisu=0._mytype
    call fine_to_coarseV(1,pre1,uvisu)
    write(filename,"('./data/pre',I4.4)") itime/imodulo
    call decomp_2d_write_one(1,uvisu,filename,2)
  endif

  if (save_prem.eq.1) then
    tb1=0._mytype
    call mean_plane_z(pre1,xsize(1),xsize(2),xsize(3),tb1(:,:,1))
    write(filename,"('./data/prem',I4.4)") itime/imodulo
    call decomp_2d_write_plane(1,tb1,3,1,filename)
  endif

  return

end subroutine VISU_PRE
#endif
!######################################################################################
subroutine mean_plane_x (f1,nx,ny,nz,fm1)

  use param, only : mytype, zero

  implicit none

  integer,intent(in) :: nx, ny, nz
  real(mytype),intent(in),dimension(nx,ny,nz) :: f1
  real(mytype),intent(out),dimension(ny,nz) :: fm1
  integer :: i,j,k

  fm1 = zero
  do j=1,ny
    do k=1,nz
      do i=1,nx
        fm1(j,k)=fm1(j,k)+f1(i,j,k)
      enddo
      fm1(j,k)=fm1(j,k)/real(nx,mytype)
    enddo
  enddo
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

  fm2 = zero
  do k=1,nz
    do i=1,nx
      do j=1,ny
        fm2(i,k)=fm2(i,k)+f2(i,j,k)
      enddo
      fm2(i,k)=fm2(i,k)/real(ny,mytype)
    enddo
  enddo
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

  fm3 = zero
  do i=1,nx
    do j=1,ny
      do k=1,nz
        fm3(i,j)=fm3(i,j)+f3(i,j,k)
      enddo
      fm3(i,j)=fm3(i,j)/real(nz,mytype)
    enddo
  enddo
  return

end subroutine mean_plane_z
