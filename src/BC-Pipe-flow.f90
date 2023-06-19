!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module pipe

  use decomp_2d
  use variables
  use param
  !!$use complex_geometry, only: tol

  implicit none

  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: geomcomplex_pipe, init_pipe
  !!$PUBLIC :: geomcomplex_pipe, init_pipe, boundary_conditions_pipe, postprocess_pipe, &
  !!$          momentum_forcing_pipe, axial_averaging

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  subroutine: geomcomplex_pipe
  !!      AUTHOR: Rodrigo Vicente Cruz
  !! DESCRIPTION: Generates epsilon matrix for pipe geometry
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !********************************************************************
  !
  subroutine geomcomplex_pipe(epsi,nxi,nxf,ny,nyi,nyf,nzi,nzf,dx,yp,dz,remp)
  !
  !********************************************************************

    use decomp_2d,only : mytype
    use MPI
    use param,only : zero,one,two,yly,zlz
    use ibm_param

    implicit none

    integer                                         :: nxi,nxf,ny,nyi,nyf,nzi,nzf
    real(mytype),dimension(nxi:nxf,nyi:nyf,nzi:nzf) :: epsi
    real(mytype),dimension(ny)                      :: yp
    real(mytype)                                    :: dx,dz,remp,tol
    !LOCALS
    real(mytype)                                    :: r,ym,zm,yc,zc
    integer                                         :: i,j,k

    epsi(:,:,:) = zero
    yc = yly/two
    zc = zlz/two
    tol=1e-15

    do k=nzi,nzf
        zm=real(k-1,mytype)*dz-zc
        do j=nyi,nyf
            ym=yp(j)-yc
            r=sqrt(ym*ym+zm*zm)
            do i=nxi,nxf
                if (r.gt.ra.and.r.lt.rao) then
                   epsi(i,j,k)=remp
                elseif (abs(r-ra).lt.tol) then
                    epsi(i,j,k)=remp
                elseif (abs(r-rao).lt.tol) then
                    epsi(i,j,k)=remp
                endif
            enddo
        enddo
    enddo
    !
    return

  end subroutine geomcomplex_pipe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  subroutine: init_pipe
  !!      AUTHOR: Rodrigo Vicente Cruz
  !! DESCRIPTION: Initial laminar conditions in pipe
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !********************************************************************
  !
  subroutine init_pipe (ux1,uy1,uz1,ep1,phi1)
  !
  !********************************************************************

    use decomp_2d
    use decomp_2d_io
    use variables
    use param
    use ibm_param
    use MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3))              :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar)    :: phi1

    real(mytype)                    :: r,ym,zm,theta
    real(mytype)                    :: um,yc,zc
    real(mytype)                    :: nu
    integer                         :: k,j,i,fh,ierror,ii,is,it,icht,code
    integer (kind=MPI_OFFSET_KIND)  :: disp
    !
    yc = yly / two
    zc = zlz / two

    if (iscalar.ne.0) then
        !Analytical laminar temperature profile, with nu=4.36
        nu=real(48./11.,8)
        do is=1,numscalar
            do k=1,xsize(3) !MIXED-TYPE BOUNDARY CONDITION : PHI = (Tw-T)/(Tw-Tb)
                zm=dz*real(xstart(3)-1+k-1,mytype)-zc 
                do j=1,xsize(2)
                    if (istret.eq.0) ym=real(j+xstart(2)-1-1,mytype)*dy-yc
                    if (istret.ne.0) ym=yp(j+xstart(2)-1)-yc
                    r=sqrt(ym*ym+zm*zm)
                    do i=1,xsize(1)
                        if (r.le.ra.and.ep1(i,j,k).eq.0) then
                            phi1(i,j,k,is)=two*nu*(three/sixteen + r**four - r**two)
                        else
                            phi1(i,j,k,is)=zero
                        endif
                    enddo
                enddo
            enddo
        enddo
    endif

    ux1=zero;uy1=zero;uz1=zero
    if (iin.ne.0) then
       call system_clock(count=code)
       if (iin.eq.2) code=0
       call random_seed(size = ii)
       call random_seed(put = code+63946*nrank*(/ (i - 1, i = 1, ii) /))
       call random_number(ux1)
       call random_number(uy1)
       call random_number(uz1)
    endif
    !modulation of the random noise + initial velocity profile
    do k=1,xsize(3)
        zm=dz*real(xstart(3)-1+k-1,mytype)-zc 
        do j=1,xsize(2)
            if (istret.eq.0) ym=real(j+xstart(2)-1-1,mytype)*dy-yc
            if (istret.ne.0) ym=yp(j+xstart(2)-1)-yc
            r=sqrt(ym*ym+zm*zm)    
            um=exp(-twenty*r*r)
            !Poiseuille flow
            bxx1(j,k)=two*(one-(ym**two+zm**two)/(ra**two))
            bxy1(j,k)=zero
            bxz1(j,k)=zero
            do i=1,xsize(1)
                if (r.le.ra.and.ep1(i,j,k).eq.0) then
                    ux1(i,j,k)=init_noise*um*(two*ux1(i,j,k)-one)+bxx1(j,k)
                    uy1(i,j,k)=init_noise*um*(two*uy1(i,j,k)-one)+bxy1(j,k)
                    uz1(i,j,k)=init_noise*um*(two*uz1(i,j,k)-one)+bxz1(j,k)
                else
                    ux1(i,j,k)=zero
                    uy1(i,j,k)=zero
                    uz1(i,j,k)=zero
                endif
            enddo
        enddo
    enddo

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init end ok'
#endif
    return

  end subroutine init_pipe
!!$  !********************************************************************
!!$  !
!!$  subroutine boundary_conditions_pipe (ux,uy,uz,phi)
!!$  !
!!$  !********************************************************************
!!$
!!$    use param
!!$    use variables
!!$    use decomp_2d
!!$
!!$    implicit none
!!$
!!$    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
!!$    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
!!$
!!$    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: gx
!!$    real(mytype) :: x, y, z
!!$    integer :: i, j, k, is
!!$!
!!$    return
!!$  end subroutine boundary_conditions_pipe
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !!
!!$  !!  subroutine: pipe_flrt
!!$  !!      AUTHOR: Rodrigo Vicente Cruz
!!$  !! DESCRIPTION: Adjustement of bulk velocity to keep constant 
!!$  !!              flow rate
!!$  !!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !********************************************************************
!!$  !
!!$  subroutine pipe_flrt(ux,uy,uz,ep,constant)
!!$  !
!!$  !********************************************************************
!!$
!!$    use decomp_2d
!!$    use decomp_2d_poisson
!!$    use variables
!!$    use param
!!$    use var
!!$    use ibm, ONLY: ra
!!$    use MPI
!!$
!!$    implicit none
!!$    real(mytype),dimension(xsize(1),xsize(2),xsize(3))  :: ux,uy,uz,ep
!!$    real(mytype)                                        :: constant
!!$    real(mytype)                                        :: qm,qmm
!!$    real(mytype)                                        :: ym,zm,yc,zc,r
!!$    real(mytype)                                        :: ncount,ncountt
!!$    integer                                             :: ivar
!!$    integer                                             :: is,j,i,k,code
!!$    character(len=30)                                   :: filename
!!$
!!$    yc = yly / two
!!$    zc = zlz / two
!!$
!!$    if (itime.eq.ifirst.and.nrank.eq.0) then
!!$       open(96,file='Ub.dat',status='unknown')
!!$    endif
!!$
!!$    !--------------------------- Bulk Velocity ---------------------------
!!$    !Calculate loss of streamwise mean pressure gradient
!!$    qm=zero
!!$    ncount=zero
!!$    do k=1,xsize(3)
!!$        zm=dz*real(xstart(3)-1+k-1,mytype)-zc 
!!$        do j=1,xsize(2)
!!$            if (istret.eq.0) ym=real(j+xstart(2)-1-1,mytype)*dy-yc
!!$            if (istret.ne.0) ym=yp(j+xstart(2)-1)-yc
!!$            r=sqrt(ym*ym+zm*zm)     
!!$            do i=1,xsize(1)
!!$                if (r.le.ra.and.ep(i,j,k).eq.0) then
!!$                    qm=qm+ux(i,j,k)
!!$                    ncount=ncount+one
!!$                endif
!!$            enddo
!!$        enddo
!!$    enddo
!!$    call MPI_ALLREDUCE(qm,qmm,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
!!$    call MPI_ALLREDUCE(ncount,ncountt,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
!!$    qmm=qmm/ncountt
!!$    if (nrank==0) then
!!$       print *,'Velocity:'
!!$       print *,'    Mean velocity before',qmm
!!$       write(96,*) real((itime-1)*dt,mytype), qmm
!!$    endif
!!$
!!$    !Correction
!!$    do j=1,xsize(2)
!!$        if (istret.eq.0) ym=real(j+xstart(2)-1-1,mytype)*dy-yc
!!$        if (istret.ne.0) ym=yp(j+xstart(2)-1)-yc
!!$        do k=1,xsize(3)
!!$            zm=dz*real(xstart(3)-1+k-1,mytype)-zc 
!!$            r=sqrt(ym*ym+zm*zm)     
!!$            do i=1,xsize(1)
!!$                if (r.le.ra.and.ep(i,j,k).eq.0) then
!!$                    ux(i,j,k)=ux(i,j,k)+(constant-qmm)
!!$                endif
!!$            enddo
!!$        enddo
!!$    enddo
!!$
!!$    !Check new bulk velocity
!!$    qmm     = zero
!!$    ncountt = zero
!!$    qm      = zero
!!$    ncount  = zero
!!$    do k=1,xsize(3)
!!$        zm=dz*real(xstart(3)-1+k-1,mytype)-zc 
!!$        do j=1,xsize(2)
!!$            if (istret.eq.0) ym=real(j+xstart(2)-1-1,mytype)*dy-yc
!!$            if (istret.ne.0) ym=yp(j+xstart(2)-1)-yc
!!$            r=sqrt(ym*ym+zm*zm)     
!!$            do i=1,xsize(1)
!!$                if (r.le.ra.and.ep(i,j,k).eq.0) then
!!$                    qm=qm+ux(i,j,k)
!!$                    ncount=ncount+one
!!$                else
!!$                    !Cancel solid zone (ra <= r <= ra+wt) 
!!$                    !and buffer zone (r > ra+wt)
!!$                    ux(i,j,k)=zero
!!$                    uy(i,j,k)=zero
!!$                    uz(i,j,k)=zero
!!$                endif
!!$            enddo
!!$        enddo
!!$    enddo
!!$    call MPI_ALLREDUCE(qm,qmm,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
!!$    call MPI_ALLREDUCE(ncount,ncountt,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
!!$    qmm=qmm/ncountt
!!$    if (nrank==0) print *,'    Mean velocity  after',qmm
!!$        
!!$    return
!!$  end subroutine pipe_flrt
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !!
!!$  !!  subroutine: pipe_blkt
!!$  !!      AUTHOR: Rodrigo Vicente Cruz
!!$  !! DESCRIPTION: Adjustment of bulk temperature according to the 
!!$  !!              thermal boundary condition
!!$  !!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !********************************************************************
!!$  !
!!$  subroutine pipe_blkt(ux,phi,ep,is)
!!$  !
!!$  !********************************************************************
!!$
!!$    use decomp_2d
!!$    use decomp_2d_poisson
!!$    use variables
!!$    use param
!!$    use ibm, ONLY: ra
!!$    use MPI
!!$
!!$    implicit none
!!$    real(mytype),dimension(xsize(1),xsize(2),xsize(3))  :: phi,ux,ep
!!$    real(mytype)                                        :: constant
!!$    real(mytype)                                        :: qv,qvm,qm,qmm
!!$    real(mytype)                                        :: ym,zm,yc,zc,r
!!$    real(mytype)                                        :: ncount,ncountt
!!$    real(mytype)                                        :: phi_out !smoothness for reconstruction
!!$    integer                                             :: ifile,is,j,i,k,code
!!$    character(len=30)                                   :: filename
!!$
!!$    if (iscalar.eq.0) return
!!$
!!$!255 format('Tb_correc_',I2.2,'.dat')
!!$255 format('Tb_',I2.2,'.dat')
!!$256 format(' Scalar:                       #',I2)
!!$
!!$    ifile=50+is
!!$    if (itime.eq.ifirst.and.nrank.eq.0) then
!!$        write(filename,255) is
!!$        open(ifile,file=filename,status='unknown')
!!$    endif
!!$
!!$    if (itbc(is).eq.1) then     !1.MBC
!!$        constant=one
!!$        phi_out =zero
!!$    elseif (itbc(is).eq.2) then !2.IF
!!$        constant=zero
!!$        phi_out =-one/nuw(is)
!!$    elseif (itbc(is).eq.3) then !3.CHT
!!$        constant=zero
!!$        phi_out =-one/nuw(is)
!!$    endif
!!$
!!$    yc = yly / two
!!$    zc = zlz / two
!!$    !--------------------------- Bulk Temperature ---------------------------
!!$    !                  with corrected streamwise velocity
!!$    qm=zero
!!$    qv=zero
!!$    ncount=zero
!!$    do k=1,xsize(3)
!!$        zm=dz*real(xstart(3)-1+k-1,mytype)-zc 
!!$        do j=1,xsize(2)
!!$            if (istret.eq.0) ym=real(j+xstart(2)-1-1,mytype)*dy-yc
!!$            if (istret.ne.0) ym=yp(j+xstart(2)-1)-yc
!!$            r=sqrt(ym*ym+zm*zm)     
!!$            do i=1,xsize(1)
!!$                if (r.le.ra.and.ep(i,j,k).eq.0) then
!!$                    qm=qm+ux(i,j,k)*phi(i,j,k)
!!$                    qv=qv+ux(i,j,k)*ux(i,j,k)
!!$                    ncount=ncount+one
!!$                endif
!!$            enddo
!!$        enddo
!!$    enddo
!!$    call MPI_ALLREDUCE(qm,qmm,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
!!$    call MPI_ALLREDUCE(qv,qvm,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
!!$    call MPI_ALLREDUCE(ncount,ncountt,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
!!$    qmm=qmm/ncountt
!!$    qvm=qvm/ncountt
!!$    if (nrank.eq.0) then
!!$        write(*,256) is
!!$        print *,'         Mean phi before',qmm
!!$        write(ifile,*) real(itime*dt,mytype),(constant-qmm)/(qvm*dt)
!!$    endif
!!$
!!$    !Correction
!!$    do j=1,xsize(2)
!!$        if (istret.eq.0) ym=real(j+xstart(2)-1-1,mytype)*dy-yc
!!$        if (istret.ne.0) ym=yp(j+xstart(2)-1)-yc
!!$        do k=1,xsize(3)
!!$            zm=dz*real(xstart(3)-1+k-1,mytype)-zc 
!!$            r=sqrt(ym*ym+zm*zm)     
!!$            do i=1,xsize(1)
!!$                if (r.le.ra.and.ep(i,j,k).eq.0) then
!!$                    phi(i,j,k)=phi(i,j,k)+ux(i,j,k)*((constant-qmm)/qvm)
!!$                endif
!!$            enddo
!!$        enddo
!!$    enddo
!!$
!!$    !Check new bulk temperature
!!$    qmm     = zero
!!$    ncountt = zero
!!$    qm      = zero
!!$    ncount  = zero
!!$    do k=1,xsize(3)
!!$        zm=dz*real(xstart(3)-1+k-1,mytype)-zc 
!!$        do j=1,xsize(2)
!!$            if (istret.eq.0) ym=real(j+xstart(2)-1-1,mytype)*dy-yc
!!$            if (istret.ne.0) ym=yp(j+xstart(2)-1)-yc
!!$            r=sqrt(ym*ym+zm*zm)     
!!$            do i=1,xsize(1)
!!$                if (r.le.ra.and.ep(i,j,k).eq.0) then
!!$                    qm=qm+ux(i,j,k)*phi(i,j,k)
!!$                    ncount=ncount+one
!!$                else !smoothness for reconstruction
!!$                    phi(i,j,k)=phi_out
!!$                endif
!!$            enddo
!!$        enddo
!!$    enddo
!!$    call MPI_ALLREDUCE(qm,qmm,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
!!$    call MPI_ALLREDUCE(ncount,ncountt,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
!!$    qmm=qmm/ncountt
!!$    if (nrank==0) print *,'          Mean phi after',qmm
!!$
!!$    return
!!$  end subroutine pipe_blkt
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !!
!!$  !!  subroutine: postprocess_pipe
!!$  !!      AUTHOR: Rodrigo Vicente Cruz
!!$  !! DESCRIPTION: Statistics collection for pipe
!!$  !!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !********************************************************************
!!$  !
!!$  subroutine postprocess_pipe(ux1,uy1,uz1,pp3,phi1,ep1) !By Rodrigo Vicente Cruz
!!$  !
!!$  !********************************************************************
!!$
!!$    use MPI
!!$    use decomp_2d
!!$    use decomp_2d_io
!!$    use var, only : umean,vmean,wmean,pmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
!!$    use var, only : phimean, phiphimean,uphimean
!!$    use var, only : phismean, phisphismean
!!$    use var, only : ta1, pp1, di1
!!$    use var, only : ppi3, dip3
!!$    use var, only : pp2, ppi2, dip2
!!$    use var, only : phis1
!!$
!!$    use var, ONLY : nxmsize, nymsize, nzmsize
!!$    use param
!!$    use variables
!!$
!!$    implicit none
!!$    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
!!$    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
!!$    real(mytype), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(in) :: pp3
!!$    character(len=30) :: filename
!!$
!!$    integer :: icht,is
!!$    
!!$    if (itime.lt.initstat) then
!!$       return
!!$    endif
!!$
!!$    !! Mean pressure
!!$    !WORK Z-PENCILS
!!$    call interzpv(ppi3,pp3(:,:,:,1),dip3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
!!$         (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
!!$    !WORK Y-PENCILS
!!$    call transpose_z_to_y(ppi3,pp2,ph3) !nxm nym nz
!!$    call interypv(ppi2,pp2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
!!$         (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
!!$    !WORK X-PENCILS
!!$    call transpose_y_to_x(ppi2,pp1,ph2) !nxm ny nz
!!$    call interxpv(ta1,pp1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
!!$         nxmsize,xsize(1),xsize(2),xsize(3),1)
!!$    if (iibm==2) then
!!$       ta1(:,:,:) = (one - ep1(:,:,:)) * ta1(:,:,:)
!!$    endif
!!$    call axial_averaging(ta1)
!!$    pmean(:,:)=pmean(:,:)+ta1(1,:,:)
!!$
!!$    !! Mean velocity
!!$    !!umean
!!$    if (iibm==2) then
!!$       ta1(:,:,:) = (one - ep1(:,:,:)) * ux1(:,:,:)
!!$    else
!!$       ta1(:,:,:) = ux1(:,:,:)
!!$    endif
!!$    call axial_averaging(ta1)
!!$    umean(:,:)=umean(:,:)+ta1(1,:,:)
!!$    !!vmean
!!$    if (iibm==2) then
!!$       ta1(:,:,:) = (one - ep1(:,:,:)) * uy1(:,:,:)
!!$    else
!!$       ta1(:,:,:) = uy1(:,:,:)
!!$    endif
!!$    call axial_averaging(ta1)
!!$    vmean(:,:)=vmean(:,:)+ta1(1,:,:)
!!$    !!wmean
!!$    if (iibm==2) then
!!$       ta1(:,:,:) = (one - ep1(:,:,:)) * uz1(:,:,:)
!!$    else
!!$       ta1(:,:,:) = uz1(:,:,:)
!!$    endif
!!$    call axial_averaging(ta1)
!!$    wmean(:,:)=wmean(:,:)+ta1(1,:,:)
!!$
!!$
!!$    !! Second-order velocity moments
!!$    !!uumean
!!$    if (iibm==2) then
!!$       ta1(:,:,:) = (one - ep1(:,:,:)) * ux1(:,:,:)*ux1(:,:,:)
!!$    endif
!!$    call axial_averaging(ta1)
!!$    uumean(:,:)=uumean(:,:)+ta1(1,:,:)
!!$    !!vvmean
!!$    if (iibm==2) then
!!$       ta1(:,:,:) = (one - ep1(:,:,:)) * uy1(:,:,:)*uy1(:,:,:)
!!$    endif
!!$    call axial_averaging(ta1)
!!$    vvmean(:,:)=vvmean(:,:)+ta1(1,:,:)
!!$    !!wwmean
!!$    if (iibm==2) then
!!$       ta1(:,:,:) = (one - ep1(:,:,:)) * uz1(:,:,:)*uz1(:,:,:)
!!$    endif
!!$    call axial_averaging(ta1)
!!$    wwmean(:,:)=wwmean(:,:)+ta1(1,:,:)
!!$    !!uvmean
!!$    if (iibm==2) then
!!$       ta1(:,:,:) = (one - ep1(:,:,:)) * ux1(:,:,:)*uy1(:,:,:)
!!$    endif
!!$    call axial_averaging(ta1)
!!$    uvmean(:,:)=uvmean(:,:)+ta1(1,:,:)
!!$    !!uwmean
!!$    if (iibm==2) then
!!$       ta1(:,:,:) = (one - ep1(:,:,:)) * ux1(:,:,:)*uz1(:,:,:)
!!$    endif
!!$    call axial_averaging(ta1)
!!$    uwmean(:,:)=uwmean(:,:)+ta1(1,:,:)
!!$    !!vwmean
!!$    if (iibm==2) then
!!$       ta1(:,:,:) = (one - ep1(:,:,:)) * uy1(:,:,:)*uz1(:,:,:)
!!$    endif
!!$    call axial_averaging(ta1)
!!$    vwmean(:,:)=vwmean(:,:)+ta1(1,:,:)
!!$
!!$    if (iscalar.ne.0) then
!!$        icht=0
!!$        do is=1, numscalar
!!$
!!$            !phimean=phi1
!!$            if (iibm==2) then
!!$               ta1(:,:,:) = (one - ep1(:,:,:)) * phi1(:,:,:,is)
!!$            else
!!$               ta1(:,:,:) = phi1(:,:,:,is)
!!$            endif
!!$            call axial_averaging(ta1)
!!$            phimean(:,:,is)=phimean(:,:,is)+ta1(1,:,:)
!!$
!!$            !phiphimean=phi1*phi1
!!$            if (iibm==2) then
!!$               ta1(:,:,:) = (one - ep1(:,:,:)) * phi1(:,:,:,is)*phi1(:,:,:,is)
!!$            else
!!$               ta1(:,:,:) = phi1(:,:,:,is)*phi1(:,:,:,is)
!!$            endif
!!$            call axial_averaging(ta1)
!!$            phiphimean(:,:,is)=phiphimean(:,:,is)+ta1(1,:,:)
!!$
!!$            !uphimean=ux1*phi1
!!$            if (iibm==2) then
!!$               ta1(:,:,:) = (one - ep1(:,:,:)) * ux1(:,:,:)*phi1(:,:,:,is)
!!$            else
!!$               ta1(:,:,:) = ux1(:,:,:)*phi1(:,:,:,is)
!!$            endif
!!$            call axial_averaging(ta1)
!!$            uphimean(:,:,is)=uphimean(:,:,is)+ta1(1,:,:)
!!$
!!$            if (itbc(is).eq.3) then !Conjugate Heat Transfer, solid statistics
!!$                icht=icht+1
!!$
!!$                !phismean=phis1
!!$                ta1(:,:,:) = ep1(:,:,:) * phis1(:,:,:,icht)
!!$                call axial_averaging(ta1)
!!$                phismean(:,:,icht)=phismean(:,:,icht)+ta1(1,:,:)
!!$
!!$                !phisphismean=phis1*phis1
!!$                ta1(:,:,:) = ep1(:,:,:) * phis1(:,:,:,icht)*phis1(:,:,:,icht)
!!$                call axial_averaging(ta1)
!!$                phisphismean(:,:,icht)=phisphismean(:,:,icht)+ta1(1,:,:)
!!$            endif
!!$
!!$        enddo
!!$    endif
!!$
!!$    if (mod(itime,icheckpoint)==0) then
!!$
!!$        if (nrank==0) then
!!$           print*,'===========================================================<<<<<'
!!$           print *,'Writing stat file',itime
!!$        endif
!!$
!!$        write(filename,"('pmean.dat',I7.7)") itime
!!$        ta1(1,:,:)=pmean(:,:)
!!$        call decomp_2d_write_plane(1,ta1,1,1,filename)
!!$        write(filename,"('umean.dat',I7.7)") itime
!!$        ta1(1,:,:)=umean(:,:)
!!$        call decomp_2d_write_plane(1,ta1,1,1,filename)
!!$        write(filename,"('vmean.dat',I7.7)") itime
!!$        ta1(1,:,:)=vmean(:,:)
!!$        call decomp_2d_write_plane(1,ta1,1,1,filename)
!!$        write(filename,"('wmean.dat',I7.7)") itime
!!$        ta1(1,:,:)=wmean(:,:)
!!$        call decomp_2d_write_plane(1,ta1,1,1,filename)
!!$
!!$        write(filename,"('uumean.dat',I7.7)") itime
!!$        ta1(1,:,:)=uumean(:,:)
!!$        call decomp_2d_write_plane(1,ta1,1,1,filename)
!!$        write(filename,"('vvmean.dat',I7.7)") itime
!!$        ta1(1,:,:)=vvmean(:,:)
!!$        call decomp_2d_write_plane(1,ta1,1,1,filename)
!!$        write(filename,"('wwmean.dat',I7.7)") itime
!!$        ta1(1,:,:)=wwmean(:,:)
!!$        call decomp_2d_write_plane(1,ta1,1,1,filename)
!!$        write(filename,"('uvmean.dat',I7.7)") itime
!!$        ta1(1,:,:)=uvmean(:,:)
!!$        call decomp_2d_write_plane(1,ta1,1,1,filename)
!!$        write(filename,"('uwmean.dat',I7.7)") itime
!!$        ta1(1,:,:)=uwmean(:,:)
!!$        call decomp_2d_write_plane(1,ta1,1,1,filename)
!!$        write(filename,"('vwmean.dat',I7.7)") itime
!!$        ta1(1,:,:)=vwmean(:,:)
!!$        call decomp_2d_write_plane(1,ta1,1,1,filename)
!!$
!!$        write(filename,"('kmean.dat',I7.7)") itime
!!$        ta1(1,:,:)=half*(uumean(:,:)+vvmean(:,:)+wwmean(:,:))
!!$        call decomp_2d_write_plane(1,ta1,1,1,filename)
!!$
!!$        if (iscalar==1) then
!!$            icht=0
!!$            do is=1, numscalar
!!$                write(filename,"('phi',I2.2,'mean.dat',I7.7)") is, itime
!!$                ta1(1,:,:)=phimean(:,:,is)
!!$                call decomp_2d_write_plane(1,ta1,1,1,filename)
!!$                write(filename,"('phiphi',I2.2,'mean.dat',I7.7)") is, itime
!!$                ta1(1,:,:)=phiphimean(:,:,is)
!!$                call decomp_2d_write_plane(1,ta1,1,1,filename)
!!$                write(filename,"('uphi',I2.2,'mean.dat',I7.7)") is, itime
!!$                ta1(1,:,:)=uphimean(:,:,is)
!!$                call decomp_2d_write_plane(1,ta1,1,1,filename)
!!$
!!$                if (itbc(is).eq.3) then !Conjugate Heat Transfer, solid statistics
!!$                    icht=icht+1
!!$                    write(filename,"('phis',I2.2,'mean.dat',I7.7)") is, itime
!!$                    ta1(1,:,:)=phismean(:,:,icht)
!!$                    call decomp_2d_write_plane(1,ta1,1,1,filename)
!!$                    write(filename,"('phisphis',I2.2,'mean.dat',I7.7)") is, itime
!!$                    ta1(1,:,:)=phisphismean(:,:,icht)
!!$                    call decomp_2d_write_plane(1,ta1,1,1,filename)
!!$                endif
!!$           enddo
!!$        endif
!!$        !if (nrank==0) then !! Cleanup old files
!!$        !    if ((itime - icheckpoint).ge.initstat) then
!!$        !        write(filename,"('uphi',I2.2,'mean.dat',I7.7)") is, itime-icheckpoint
!!$        !        call system ("rm " //filename)
!!$        !    endif
!!$        !endif
!!$
!!$    endif
!!$
!!$    return
!!$  end subroutine postprocess_pipe
!!$  !********************************************************************
!!$  !
!!$  subroutine momentum_forcing_pipe(dux1, duy1, ux1, uy1)
!!$  !
!!$  !********************************************************************
!!$
!!$    implicit none
!!$
!!$    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ux1, uy1
!!$    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1
!!$
!!$    return
!!$
!!$  ENDsubroutine momentum_forcing_pipe
!!$
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    !!
!!$    !!  subroutine: set_nbc_coefficients
!!$    !!      AUTHOR: Rodrigo Vicente Cruz
!!$    !! DESCRIPTION: Prepares coefficients of: i)3rd order non-centred derivative
!!$    !!              scheme (*nscy and *nscz); ii) Non-centred extrapo-
!!$    !!              lation scheme (*enscy and *enscz). Useful for Imposed Flux
!!$    !!              (IF) and Conjugate Heat Transfer (CHT) conditions.        
!!$    !!
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !***************************************************************************
!!$  !
!!$  subroutine set_nbc_coefficients()
!!$  !
!!$  !***************************************************************************
!!$  use decomp_2d
!!$  use complex_geometry, ONLY: nobjmax,nobjy,nobjz,yi,yf,zi,zf
!!$  use variables, ONLY: yp,numscalar
!!$  use param, ONLY: one,two,three,dy,dz
!!$  !
!!$  implicit none
!!$  !
!!$  real(mytype)                          :: dyb,dybs,dzb,dzbs
!!$  real(mytype)                          :: r1,r2,r3
!!$  real(mytype)                          :: s0,s1,s2,s3
!!$  integer                               :: nibmax
!!$  integer                               :: jy,kz
!!$  integer                               :: i,j,k,count
!!$  integer                               :: code
!!$   
!!$  nibmax=2 !Max number of immersed boundaries, 2 for each fully immersed object
!!$
!!$  !3rd non-centred first derivative 
!!$  allocate(anscy(nibmax,nobjmax,ysize(3)),bnscy(nibmax,nobjmax,ysize(3)),&
!!$           cnscy(nibmax,nobjmax,ysize(3)),dnscy(nibmax,nobjmax,ysize(3)))
!!$  allocate(anscz(nibmax,nobjmax,zsize(2)),bnscz(nibmax,nobjmax,zsize(2)),&
!!$           cnscz(nibmax,nobjmax,zsize(2)),dnscz(nibmax,nobjmax,zsize(2)))
!!$  !Non-centred boundary extrapolation
!!$  allocate(aenscy(nibmax,nobjmax,ysize(3)),benscy(nibmax,nobjmax,ysize(3)),&
!!$           censcy(nibmax,nobjmax,ysize(3)),denscy(nibmax,nobjmax,ysize(3)))
!!$  allocate(aenscz(nibmax,nobjmax,zsize(2)),benscz(nibmax,nobjmax,zsize(2)),&
!!$           censcz(nibmax,nobjmax,zsize(2)),denscz(nibmax,nobjmax,zsize(2)))
!!$  !Temperature at the wall and tangential heat flux
!!$  allocate(nuw(numscalar))
!!$  allocate(phiw2(nibmax,nobjmax,ysize(1),ysize(3),numscalar),&
!!$           qthetaw2(nibmax,nobjmax,ysize(1),ysize(3),numscalar))
!!$  allocate(phiw3(nibmax,nobjmax,zsize(1),zsize(2),numscalar),&
!!$           qthetaw3(nibmax,nobjmax,zsize(1),zsize(2),numscalar))
!!$
!!$  if (numcht.ne.0) then !If CHT, solid field
!!$      !3rd non-centred first derivative 
!!$      allocate(anscy_s(nibmax,nobjmax,ysize(3)),bnscy_s(nibmax,nobjmax,ysize(3)),&
!!$               cnscy_s(nibmax,nobjmax,ysize(3)),dnscy_s(nibmax,nobjmax,ysize(3)))
!!$      allocate(anscz_s(nibmax,nobjmax,zsize(2)),bnscz_s(nibmax,nobjmax,zsize(2)),&
!!$               cnscz_s(nibmax,nobjmax,zsize(2)),dnscz_s(nibmax,nobjmax,zsize(2)))
!!$      !Non-centred boundary extrapolation
!!$      allocate(aenscy_s(nibmax,nobjmax,ysize(3)),benscy_s(nibmax,nobjmax,ysize(3)),&
!!$               censcy_s(nibmax,nobjmax,ysize(3)),denscy_s(nibmax,nobjmax,ysize(3)))
!!$      allocate(aenscz_s(nibmax,nobjmax,zsize(2)),benscz_s(nibmax,nobjmax,zsize(2)),&
!!$               censcz_s(nibmax,nobjmax,zsize(2)),denscz_s(nibmax,nobjmax,zsize(2)))
!!$      !Temperature at the wall
!!$      !!====DEBUG
!!$      !allocate(phiws2(nibmax,nobjmax,ysize(1),ysize(3),numscalar),&
!!$      !         phiws3(nibmax,nobjmax,zsize(1),zsize(2),numscalar))
!!$      !!=========
!!$      allocate(phiws2(nibmax,nobjmax,ysize(1),ysize(3),numcht),&
!!$               phiws3(nibmax,nobjmax,zsize(1),zsize(2),numcht))
!!$      !!Wall heat flux
!!$      !allocate(wf2(nibmax,nobjmax,ysize(1),ysize(3),numscalar),&
!!$      !         wf3(nibmax,nobjmax,zsize(1),zsize(2),numscalar))
!!$      allocate(qnw2(nibmax,nobjmax,ysize(1),ysize(3),numscalar),&
!!$               qnw3(nibmax,nobjmax,zsize(1),zsize(2),numscalar))
!!$      allocate(qnws2(nibmax,nobjmax,ysize(1),ysize(3),numcht),&
!!$               qnws3(nibmax,nobjmax,zsize(1),zsize(2),numcht))
!!$      allocate(qthetaws2(nibmax,nobjmax,ysize(1),ysize(3),numcht),&
!!$               qthetaws3(nibmax,nobjmax,zsize(1),zsize(2),numcht))
!!$  endif
!!$
!!$  i=1 !arbitrary for simple pipe geometry
!!$  !================ Y-DIRECTION COEFFICIENTS ================
!!$  !Derivative scheme
!!$  anscy(:,:,:)=0.
!!$  bnscy(:,:,:)=0.
!!$  cnscy(:,:,:)=0.
!!$  dnscy(:,:,:)=0.
!!$  !Extrapolation scheme
!!$  aenscy(:,:,:)=0.
!!$  benscy(:,:,:)=0.
!!$  censcy(:,:,:)=0.
!!$  denscy(:,:,:)=0.
!!$  do k=1,ysize(3)
!!$     !if (nobjy(i,k).eq.2) then !2 objects, 4 immersed boundaries
!!$     if (nobjy(i,k).ne.0) then !Immersed objects
!!$
!!$        !FLUID FIELD COEFFICIENTS
!!$        do j=1,nobjy(i,k)
!!$            !1ST IMMERSED BOUNDARY
!!$                jy=1!jy=yi(j,i,k)/dy+1
!!$                do while(yp(jy).lt.yi(j,i,k))
!!$                   jy=jy+1
!!$                enddo
!!$                jy=jy-1
!!$                dyb=abs(yp(jy)-yi(j,i,k))       !distance of closest node without skip
!!$                dybs=abs(yp(jy-1)-yi(j,i,k))    !distance of closest node with skip
!!$                !i) First derivative: with skip
!!$                r1=(dybs+   dy)/(dybs)
!!$                r2=(dybs+two*dy)/(dybs)
!!$                r3=(dybs+three*dy)/(dybs)
!!$                anscy(1,j,k)=-((r1+one)*r2+r1)/(r1*r2*dybs)
!!$                bnscy(1,j,k)=(r1*r2)/(((r1-one)*r2-r1+one)*dybs)
!!$                cnscy(1,j,k)=-(r2)/((((r1**two-r1)*r2)-(r1**three)+(r1**two))*dybs)
!!$                dnscy(1,j,k)=(r1)/(((r2**three)+((-r1-one)*r2**two)+(r1*r2))*dybs) 
!!$                !ii) Extrapolation: with or without skip
!!$                if (iskip_e.eq.0) then !without skip 
!!$                    !====DEBUG
!!$                    !r1=(dyb+   dy)/(dyb)
!!$                    !r2=(dyb+two*dy)/(dyb)
!!$                    !r3=(dyb+three*dy)/(dyb)
!!$                    !=========
!!$                    s0=(dyb)
!!$                    s1=(dyb+dy)
!!$                    s2=(dyb+two*dy)
!!$                    s3=(dyb+three*dy)
!!$                else                   !with skip
!!$                    !====DEBUG
!!$                    !r1=(dybs+   dy)/(dybs)
!!$                    !r2=(dybs+two*dy)/(dybs)
!!$                    !r3=(dybs+three*dy)/(dybs)
!!$                    !=========
!!$                    s0=(dybs)
!!$                    s1=(dybs+dy)
!!$                    s2=(dybs+two*dy)
!!$                    s3=(dybs+three*dy)
!!$                endif
!!$                if (iextp.eq.1) then !1st order
!!$                    aenscy(1,j,k)=one
!!$                    benscy(1,j,k)=zero
!!$                    censcy(1,j,k)=zero
!!$                    denscy(1,j,k)=zero
!!$                elseif (iextp.eq.2) then !2nd order
!!$                    !!====DEBUG
!!$                    !aenscy(1,j,k)=r1/(r1-one)
!!$                    !benscy(1,j,k)=-one/(r1-one)
!!$                    !censcy(1,j,k)=zero
!!$                    !denscy(1,j,k)=zero
!!$                    !!=========
!!$                    aenscy(1,j,k)=s1/(s1-s0)
!!$                    benscy(1,j,k)=-(s0)/(s1-s0)
!!$                    censcy(1,j,k)=zero
!!$                    denscy(1,j,k)=zero
!!$                elseif (iextp.eq.3) then !3rd order
!!$                    !!====DEBUG
!!$                    !aenscy(1,j,k)=r1*r2/((r1-one)*r2-r1+one)
!!$                    !benscy(1,j,k)=-r2/((r1-one)*r2-r1*(r1-one))
!!$                    !censcy(1,j,k)=r1/(r2*r2-(r1+one)*r2+r1)
!!$                    !denscy(1,j,k)=zero
!!$                    !!=========
!!$                    aenscy(1,j,k)=s1*s2/((s1-s0)*s2-s0*s1+s0*s0)
!!$                    benscy(1,j,k)=-(s0*s2)/((s1-s0)*s2-s1*s1+s0*s1)
!!$                    censcy(1,j,k)=(s0*s1)/(s2*s2-s2*(s0+s1)+s0*s1)
!!$                    denscy(1,j,k)=zero
!!$                elseif (iextp.eq.4) then !4th order
!!$                    !!====DEBUG
!!$                    !aenscy(1,j,k)=r1*r2*r3/(((r1-one)*r2-r1+one)*r3+(one-r1)*r2+r1-one)
!!$                    !benscy(1,j,k)=-r2*r3/(((r1-one)*r2-r1*r1+r1)*r3+(one-r1)*r1*r2+r1*r1*(r1-one))
!!$                    !censcy(1,j,k)=r1*r3/((r2*r2-(r1+one)*r2+r1)*r3-r2*r2*r2+(r1+one)*r2*r2-r1*r2)
!!$                    !denscy(1,j,k)=-r1*r2/(r3*r3*r3-(r1+r2+one)*r3*r3+((r1+one)*r2+r1)*r3-r1*r2)         
!!$                    !!=========
!!$                    aenscy(1,j,k)= (s1*s2*s3)/&
!!$                                   ((s2*(s1-s0)-s0*s1+s0*s0)*s3+(s0*s0-s0*s1)*s2+s0*s0*s1-s0*s0*s0)
!!$                    benscy(1,j,k)=-(s0*s2*s3)/&
!!$                                   ((s2*(s1-s0)-s1*s1+s0*s1)*s3+(s0*s1-s1*s1)*s2+s1*s1*s1-s0*s1*s1)
!!$                    censcy(1,j,k)= (s0*s1*s3)/&
!!$                                   ((s2*s2-s2*(s1+s0)+s0*s1)*s3-s2*s2*s2+s2*s2*(s1+s0)-s0*s1*s2)
!!$                    denscy(1,j,k)=-(s0*s1*s2)/&
!!$                                   (s3*s3*s3-s3*s3*(s0+s1+s2)+((s1+s0)*s2+s0*s1)*s3-s0*s1*s2)
!!$                endif
!!$
!!$            !2ND IMMERSED BOUNDARY
!!$                jy=1!jy=(yf(j,i,k)+dy)/dy+1
!!$                !====DEBUG
!!$                !do while(yp(jy).le.yf(j,i,k)) !use if exact_ib_pipe
!!$                do while(yp(jy).lt.yf(j,i,k))
!!$                   jy=jy+1
!!$                enddo
!!$                dyb=abs(yp(jy)-yf(j,i,k))       !distance of closest node without skip
!!$                dybs=abs(yp(jy+1)-yf(j,i,k))    !distance of closest node with skip
!!$                !i) First derivative: with skip
!!$                r1=(dybs+   dy)/(dybs)
!!$                r2=(dybs+two*dy)/(dybs)
!!$                r3=(dybs+three*dy)/(dybs)
!!$                anscy(2,j,k)=-((r1+one)*r2+r1)/(r1*r2*dybs)
!!$                bnscy(2,j,k)=(r1*r2)/(((r1-one)*r2-r1+one)*dybs)
!!$                cnscy(2,j,k)=-(r2)/((((r1**two-r1)*r2)-(r1**three)+(r1**two))*dybs)
!!$                dnscy(2,j,k)=(r1)/(((r2**three)+((-r1-one)*r2**two)+(r1*r2))*dybs) 
!!$                !ii) Extrapolation: with or without skip
!!$                if (iskip_e.eq.0) then !without skip 
!!$                    !====DEBUG
!!$                    !r1=(dyb+   dy)/(dyb)
!!$                    !r2=(dyb+two*dy)/(dyb)
!!$                    !r3=(dyb+three*dy)/(dyb)
!!$                    !=========
!!$                    s0=(dyb)
!!$                    s1=(dyb+dy)
!!$                    s2=(dyb+two*dy)
!!$                    s3=(dyb+three*dy)
!!$                else                   !with skip
!!$                    !====DEBUG
!!$                    !r1=(dybs+   dy)/(dybs)
!!$                    !r2=(dybs+two*dy)/(dybs)
!!$                    !r3=(dybs+three*dy)/(dybs)
!!$                    !=========
!!$                    s0=(dybs)
!!$                    s1=(dybs+dy)
!!$                    s2=(dybs+two*dy)
!!$                    s3=(dybs+three*dy)
!!$                endif
!!$                if (iextp.eq.1) then !1st order
!!$                    aenscy(2,j,k)=one
!!$                    benscy(2,j,k)=zero
!!$                    censcy(2,j,k)=zero
!!$                    denscy(2,j,k)=zero
!!$                elseif (iextp.eq.2) then !2nd order
!!$                    !====DEBUG
!!$                    !aenscy(2,j,k)=r1/(r1-one)
!!$                    !benscy(2,j,k)=-one/(r1-one)
!!$                    !censcy(2,j,k)=zero
!!$                    !denscy(2,j,k)=zero
!!$                    !=========
!!$                    aenscy(2,j,k)=s1/(s1-s0)
!!$                    benscy(2,j,k)=-(s0)/(s1-s0)
!!$                    censcy(2,j,k)=zero
!!$                    denscy(2,j,k)=zero
!!$                elseif (iextp.eq.3) then !3rd order
!!$                    !====DEBUG
!!$                    !aenscy(2,j,k)=r1*r2/((r1-one)*r2-r1+one)
!!$                    !benscy(2,j,k)=-r2/((r1-one)*r2-r1*(r1-one))
!!$                    !censcy(2,j,k)=r1/(r2*r2-(r1+one)*r2+r1)
!!$                    !denscy(2,j,k)=zero
!!$                    !=========
!!$                    aenscy(2,j,k)=s1*s2/((s1-s0)*s2-s0*s1+s0*s0)
!!$                    benscy(2,j,k)=-(s0*s2)/((s1-s0)*s2-s1*s1+s0*s1)
!!$                    censcy(2,j,k)=(s0*s1)/(s2*s2-s2*(s0+s1)+s0*s1)
!!$                    denscy(2,j,k)=zero
!!$                elseif (iextp.eq.4) then !4th order
!!$                    !====DEBUG
!!$                    !aenscy(2,j,k)=r1*r2*r3/(((r1-one)*r2-r1+one)*r3+(one-r1)*r2+r1-one)
!!$                    !benscy(2,j,k)=-r2*r3/(((r1-one)*r2-r1*r1+r1)*r3+(one-r1)*r1*r2+r1*r1*(r1-one))
!!$                    !censcy(2,j,k)=r1*r3/((r2*r2-(r1+one)*r2+r1)*r3-r2*r2*r2+(r1+one)*r2*r2-r1*r2)
!!$                    !denscy(2,j,k)=-r1*r2/(r3*r3*r3-(r1+r2+one)*r3*r3+((r1+one)*r2+r1)*r3-r1*r2)         
!!$                    !=========
!!$                    aenscy(2,j,k)=(s1*s2*s3)/&
!!$                                  ((s2*(s1-s0)-s0*s1+s0*s0)*s3+(s0*s0-s0*s1)*s2+s0*s0*s1-s0*s0*s0)
!!$                    benscy(2,j,k)=-(s0*s2*s3)/&
!!$                                   ((s2*(s1-s0)-s1*s1+s0*s1)*s3+(s0*s1-s1*s1)*s2+s1*s1*s1-s0*s1*s1)
!!$                    censcy(2,j,k)=(s0*s1*s3)/&
!!$                                  ((s2*s2-s2*(s1+s0)+s0*s1)*s3-s2*s2*s2+s2*s2*(s1+s0)-s0*s1*s2)
!!$                    denscy(2,j,k)=-(s0*s1*s2)/&
!!$                                   (s3*s3*s3-s3*s3*(s0+s1+s2)+((s1+s0)*s2+s0*s1)*s3-s0*s1*s2)
!!$                endif
!!$
!!$        enddo
!!$
!!$        !SOLID FIELD COEFFICIENTS (IF CHT)
!!$        if (numcht.ne.0) then 
!!$        do j=1,nobjy(i,k)
!!$            !1ST IMMERSED BOUNDARY
!!$                jy=1!jy=yi(j,i,k)/dy+1
!!$                do while(yp(jy).lt.yi(j,i,k))
!!$                   jy=jy+1
!!$                enddo
!!$                jy=jy-1
!!$                dyb=abs(yp(jy+1)-yi(j,i,k))     !distance of closest node without skip
!!$                dybs=abs(yp(jy+2)-yi(j,i,k))    !distance of closest node with skip
!!$                !i) First derivative: with skip
!!$                r1=(dybs+   dy)/(dybs)
!!$                r2=(dybs+two*dy)/(dybs)
!!$                r3=(dybs+three*dy)/(dybs)
!!$                anscy_s(1,j,k)=-((r1+one)*r2+r1)/(r1*r2*dybs)
!!$                bnscy_s(1,j,k)=(r1*r2)/(((r1-one)*r2-r1+one)*dybs)
!!$                cnscy_s(1,j,k)=-(r2)/((((r1**two-r1)*r2)-(r1**three)+(r1**two))*dybs)
!!$                dnscy_s(1,j,k)=(r1)/(((r2**three)+((-r1-one)*r2**two)+(r1*r2))*dybs) 
!!$                !ii) Extrapolation: with or without skip
!!$                if (iskip_e.eq.0) then !without skip 
!!$                    s0=(dyb)
!!$                    s1=(dyb+dy)
!!$                    s2=(dyb+two*dy)
!!$                    s3=(dyb+three*dy)
!!$                else                   !with skip
!!$                    s0=(dybs)
!!$                    s1=(dybs+dy)
!!$                    s2=(dybs+two*dy)
!!$                    s3=(dybs+three*dy)
!!$                endif
!!$                if (iextp.eq.1) then !1st order
!!$                    aenscy_s(1,j,k)=one
!!$                    benscy_s(1,j,k)=zero
!!$                    censcy_s(1,j,k)=zero
!!$                    denscy_s(1,j,k)=zero
!!$                elseif (iextp.eq.2) then !2nd order
!!$                    aenscy_s(1,j,k)=s1/(s1-s0)
!!$                    benscy_s(1,j,k)=-(s0)/(s1-s0)
!!$                    censcy_s(1,j,k)=zero
!!$                    denscy_s(1,j,k)=zero
!!$                elseif (iextp.eq.3) then !3rd order
!!$                    aenscy_s(1,j,k)=s1*s2/((s1-s0)*s2-s0*s1+s0*s0)
!!$                    benscy_s(1,j,k)=-(s0*s2)/((s1-s0)*s2-s1*s1+s0*s1)
!!$                    censcy_s(1,j,k)=(s0*s1)/(s2*s2-s2*(s0+s1)+s0*s1)
!!$                    denscy_s(1,j,k)=zero
!!$                elseif (iextp.eq.4) then !4th order
!!$                    aenscy_s(1,j,k)=(s1*s2*s3)/&
!!$                                    ((s2*(s1-s0)-s0*s1+s0*s0)*s3+(s0*s0-s0*s1)*s2+s0*s0*s1-s0*s0*s0)
!!$                    benscy_s(1,j,k)=-(s0*s2*s3)/&
!!$                                     ((s2*(s1-s0)-s1*s1+s0*s1)*s3+(s0*s1-s1*s1)*s2+s1*s1*s1-s0*s1*s1)
!!$                    censcy_s(1,j,k)=(s0*s1*s3)/&
!!$                                    ((s2*s2-s2*(s1+s0)+s0*s1)*s3-s2*s2*s2+s2*s2*(s1+s0)-s0*s1*s2)
!!$                    denscy_s(1,j,k)=-(s0*s1*s2)/&
!!$                                     (s3*s3*s3-s3*s3*(s0+s1+s2)+((s1+s0)*s2+s0*s1)*s3-s0*s1*s2)
!!$                endif
!!$
!!$            !2ND IMMERSED BOUNDARY
!!$                jy=1!jy=(yf(j,i,k)+dy)/dy+1
!!$                !====DEBUG
!!$                !do while(yp(jy).le.yf(j,i,k)) !use if exact_ib_pipe
!!$                do while(yp(jy).lt.yf(j,i,k))
!!$                   jy=jy+1
!!$                enddo
!!$                dyb=abs(yp(jy-1)-yf(j,i,k))     !distance of closest node without skip
!!$                dybs=abs(yp(jy-2)-yf(j,i,k))    !distance of closest node with skip
!!$                !i) First derivative: with skip
!!$                r1=(dybs+   dy)/(dybs)
!!$                r2=(dybs+two*dy)/(dybs)
!!$                r3=(dybs+three*dy)/(dybs)
!!$                anscy_s(2,j,k)=-((r1+one)*r2+r1)/(r1*r2*dybs)
!!$                bnscy_s(2,j,k)=(r1*r2)/(((r1-one)*r2-r1+one)*dybs)
!!$                cnscy_s(2,j,k)=-(r2)/((((r1**two-r1)*r2)-(r1**three)+(r1**two))*dybs)
!!$                dnscy_s(2,j,k)=(r1)/(((r2**three)+((-r1-one)*r2**two)+(r1*r2))*dybs) 
!!$                !ii) Extrapolation: with or without skip
!!$                if (iskip_e.eq.0) then !without skip 
!!$                    s0=(dyb)
!!$                    s1=(dyb+dy)
!!$                    s2=(dyb+two*dy)
!!$                    s3=(dyb+three*dy)
!!$                else                   !with skip
!!$                    s0=(dybs)
!!$                    s1=(dybs+dy)
!!$                    s2=(dybs+two*dy)
!!$                    s3=(dybs+three*dy)
!!$                endif
!!$                if (iextp.eq.1) then !1st order
!!$                    aenscy_s(2,j,k)=one
!!$                    benscy_s(2,j,k)=zero
!!$                    censcy_s(2,j,k)=zero
!!$                    denscy_s(2,j,k)=zero
!!$                elseif (iextp.eq.2) then !2nd order
!!$                    aenscy_s(2,j,k)=s1/(s1-s0)
!!$                    benscy_s(2,j,k)=-(s0)/(s1-s0)
!!$                    censcy_s(2,j,k)=zero
!!$                    denscy_s(2,j,k)=zero
!!$                elseif (iextp.eq.3) then !3rd order
!!$                    aenscy_s(2,j,k)=s1*s2/((s1-s0)*s2-s0*s1+s0*s0)
!!$                    benscy_s(2,j,k)=-(s0*s2)/((s1-s0)*s2-s1*s1+s0*s1)
!!$                    censcy_s(2,j,k)=(s0*s1)/(s2*s2-s2*(s0+s1)+s0*s1)
!!$                    denscy_s(2,j,k)=zero
!!$                elseif (iextp.eq.4) then !4th order
!!$                    aenscy_s(2,j,k)=(s1*s2*s3)/&
!!$                                    ((s2*(s1-s0)-s0*s1+s0*s0)*s3+(s0*s0-s0*s1)*s2+s0*s0*s1-s0*s0*s0)
!!$                    benscy_s(2,j,k)=-(s0*s2*s3)/&
!!$                                     ((s2*(s1-s0)-s1*s1+s0*s1)*s3+(s0*s1-s1*s1)*s2+s1*s1*s1-s0*s1*s1)
!!$                    censcy_s(2,j,k)=(s0*s1*s3)/&
!!$                                    ((s2*s2-s2*(s1+s0)+s0*s1)*s3-s2*s2*s2+s2*s2*(s1+s0)-s0*s1*s2)
!!$                    denscy_s(2,j,k)=-(s0*s1*s2)/&
!!$                                     (s3*s3*s3-s3*s3*(s0+s1+s2)+((s1+s0)*s2+s0*s1)*s3-s0*s1*s2)
!!$                endif 
!!$        enddo
!!$        endif
!!$     endif    
!!$  enddo    
!!$
!!$  !================ Z-DIRECTION COEFFICIENTS ================
!!$  !Derivative scheme
!!$  anscz(:,:,:)=0.
!!$  bnscz(:,:,:)=0.
!!$  cnscz(:,:,:)=0.
!!$  dnscz(:,:,:)=0.
!!$  !Extrapolation scheme
!!$  aenscz(:,:,:)=0.
!!$  benscz(:,:,:)=0.
!!$  censcz(:,:,:)=0.
!!$  denscz(:,:,:)=0.
!!$  do j=1,zsize(2)
!!$     !if (nobjz(i,j).eq.2) then !2 objects, 4 immersed boundaries
!!$     if (nobjz(i,j).ne.0) then !Immersed objects
!!$
!!$        !FLUID FIELD COEFFICIENTS
!!$        do k=1,nobjz(i,j)
!!$            !1ST IMMERSED BOUNDARY
!!$                kz=zi(k,i,j)/dz+1
!!$                !====DEBUG
!!$                !if (abs(zp(kz)-zi(k,i,j)).lt.tol) kz=kz-1 !use if exact_ib_pipe
!!$                dzb=abs((kz-1)*dz-zi(k,i,j))    !distance of closest node without skip
!!$                dzbs=abs((kz-1-1)*dz-zi(k,i,j)) !distance of closest node with skip
!!$                !i) First derivative: with skip
!!$                r1=(dzbs+   dz)/(dzbs)
!!$                r2=(dzbs+two*dz)/(dzbs)
!!$                r3=(dzbs+three*dz)/(dzbs)
!!$                anscz(1,k,j)=-((r1+one)*r2+r1)/(r1*r2*dzbs)
!!$                bnscz(1,k,j)=(r1*r2)/(((r1-one)*r2-r1+one)*dzbs)
!!$                cnscz(1,k,j)=-(r2)/((((r1**two-r1)*r2)-(r1**three)+(r1**two))*dzbs)
!!$                dnscz(1,k,j)=(r1)/(((r2**three)+((-r1-one)*r2**two)+(r1*r2))*dzbs) 
!!$                !ii) Extrapolation: with or without skip
!!$                if (iskip_e.eq.0) then !without skip 
!!$                    s0=(dzb)
!!$                    s1=(dzb+dz)
!!$                    s2=(dzb+two*dz)
!!$                    s3=(dzb+three*dz)
!!$                else                   !with skip
!!$                    s0=(dzbs)
!!$                    s1=(dzbs+dz)
!!$                    s2=(dzbs+two*dz)
!!$                    s3=(dzbs+three*dz)
!!$                endif
!!$                if (iextp.eq.1) then !1st order
!!$                    aenscz(1,k,j)=one
!!$                    benscz(1,k,j)=zero
!!$                    censcz(1,k,j)=zero
!!$                    denscz(1,k,j)=zero
!!$                elseif (iextp.eq.2) then !2nd order
!!$                    aenscz(1,k,j)=s1/(s1-s0)
!!$                    benscz(1,k,j)=-(s0)/(s1-s0)
!!$                    censcz(1,k,j)=zero
!!$                    denscz(1,k,j)=zero
!!$                elseif (iextp.eq.3) then !3rd order
!!$                    aenscz(1,k,j)=s1*s2/((s1-s0)*s2-s0*s1+s0*s0)
!!$                    benscz(1,k,j)=-(s0*s2)/((s1-s0)*s2-s1*s1+s0*s1)
!!$                    censcz(1,k,j)=(s0*s1)/(s2*s2-s2*(s0+s1)+s0*s1)
!!$                    denscz(1,k,j)=zero
!!$                elseif (iextp.eq.4) then !4th order
!!$                    aenscz(1,k,j)=(s1*s2*s3)/&
!!$                                  ((s2*(s1-s0)-s0*s1+s0*s0)*s3+(s0*s0-s0*s1)*s2+s0*s0*s1-s0*s0*s0)
!!$                    benscz(1,k,j)=-(s0*s2*s3)/&
!!$                                   ((s2*(s1-s0)-s1*s1+s0*s1)*s3+(s0*s1-s1*s1)*s2+s1*s1*s1-s0*s1*s1)
!!$                    censcz(1,k,j)=(s0*s1*s3)/&
!!$                                  ((s2*s2-s2*(s1+s0)+s0*s1)*s3-s2*s2*s2+s2*s2*(s1+s0)-s0*s1*s2)
!!$                    denscz(1,k,j)=-(s0*s1*s2)/&
!!$                                   (s3*s3*s3-s3*s3*(s0+s1+s2)+((s1+s0)*s2+s0*s1)*s3-s0*s1*s2)
!!$                endif
!!$
!!$            !2ND IMMERSED BOUNDARY
!!$                kz=(zf(k,i,j)+dz)/dz+1
!!$                dzb=abs((kz-1)*dz-zf(k,i,j))    !distance of closest node without skip
!!$                dzbs=abs((kz-1+1)*dz-zf(k,i,j)) !distance of closest node with skip
!!$                !i) First derivative: with skip
!!$                r1=(dzbs+   dz)/(dzbs)
!!$                r2=(dzbs+two*dz)/(dzbs)
!!$                r3=(dzbs+three*dz)/(dzbs)
!!$                anscz(2,k,j)=-((r1+one)*r2+r1)/(r1*r2*dzbs)
!!$                bnscz(2,k,j)=(r1*r2)/(((r1-one)*r2-r1+one)*dzbs)
!!$                cnscz(2,k,j)=-(r2)/((((r1**two-r1)*r2)-(r1**three)+(r1**two))*dzbs)
!!$                dnscz(2,k,j)=(r1)/(((r2**three)+((-r1-one)*r2**two)+(r1*r2))*dzbs) 
!!$                !ii) Extrapolation: with or without skip
!!$                if (iskip_e.eq.0) then !without skip 
!!$                    s0=(dzb)
!!$                    s1=(dzb+dz)
!!$                    s2=(dzb+two*dz)
!!$                    s3=(dzb+three*dz)
!!$                else                   !with skip
!!$                    s0=(dzbs)
!!$                    s1=(dzbs+dz)
!!$                    s2=(dzbs+two*dz)
!!$                    s3=(dzbs+three*dz)
!!$                endif
!!$                if (iextp.eq.1) then !1st order
!!$                    aenscz(2,k,j)=one
!!$                    benscz(2,k,j)=zero
!!$                    censcz(2,k,j)=zero
!!$                    denscz(2,k,j)=zero
!!$                elseif (iextp.eq.2) then !2nd order
!!$                    aenscz(2,k,j)=s1/(s1-s0)
!!$                    benscz(2,k,j)=-(s0)/(s1-s0)
!!$                    censcz(2,k,j)=zero
!!$                    denscz(2,k,j)=zero
!!$                elseif (iextp.eq.3) then !3rd order
!!$                    aenscz(2,k,j)=s1*s2/((s1-s0)*s2-s0*s1+s0*s0)
!!$                    benscz(2,k,j)=-(s0*s2)/((s1-s0)*s2-s1*s1+s0*s1)
!!$                    censcz(2,k,j)=(s0*s1)/(s2*s2-s2*(s0+s1)+s0*s1)
!!$                    denscz(2,k,j)=zero
!!$                elseif (iextp.eq.4) then !4th order
!!$                    aenscz(2,k,j)=(s1*s2*s3)/&
!!$                                  ((s2*(s1-s0)-s0*s1+s0*s0)*s3+(s0*s0-s0*s1)*s2+s0*s0*s1-s0*s0*s0)
!!$                    benscz(2,k,j)=-(s0*s2*s3)/&
!!$                                   ((s2*(s1-s0)-s1*s1+s0*s1)*s3+(s0*s1-s1*s1)*s2+s1*s1*s1-s0*s1*s1)
!!$                    censcz(2,k,j)=(s0*s1*s3)/&
!!$                                  ((s2*s2-s2*(s1+s0)+s0*s1)*s3-s2*s2*s2+s2*s2*(s1+s0)-s0*s1*s2)
!!$                    denscz(2,k,j)=-(s0*s1*s2)/&
!!$                                   (s3*s3*s3-s3*s3*(s0+s1+s2)+((s1+s0)*s2+s0*s1)*s3-s0*s1*s2)
!!$                endif
!!$
!!$        enddo
!!$
!!$        !SOLID FIELD COEFFICIENTS (IF CHT)
!!$        if (numcht.ne.0) then 
!!$        do k=1,nobjz(i,j)
!!$            !1ST IMMERSED BOUNDARY
!!$                kz=zi(k,i,j)/dz+1
!!$                !====DEBUG
!!$                !if (abs(zp(kz)-zi(k,i,j)).lt.tol) kz=kz-1 !use if exact_ib_pipe
!!$                dzb=abs(((kz+1)-1)*dz-zi(k,i,j))    !distance of closest node without skip
!!$                dzbs=abs(((kz+2)-1)*dz-zi(k,i,j))   !distance of closest node with skip
!!$                !i) First derivative: with skip
!!$                r1=(dzbs+   dz)/(dzbs)
!!$                r2=(dzbs+two*dz)/(dzbs)
!!$                r3=(dzbs+three*dz)/(dzbs)
!!$                anscz_s(1,k,j)=-((r1+one)*r2+r1)/(r1*r2*dzbs)
!!$                bnscz_s(1,k,j)=(r1*r2)/(((r1-one)*r2-r1+one)*dzbs)
!!$                cnscz_s(1,k,j)=-(r2)/((((r1**two-r1)*r2)-(r1**three)+(r1**two))*dzbs)
!!$                dnscz_s(1,k,j)=(r1)/(((r2**three)+((-r1-one)*r2**two)+(r1*r2))*dzbs) 
!!$                !ii) Extrapolation: with or without skip
!!$                if (iskip_e.eq.0) then !without skip 
!!$                    s0=(dzb)
!!$                    s1=(dzb+dz)
!!$                    s2=(dzb+two*dz)
!!$                    s3=(dzb+three*dz)
!!$                else                   !with skip
!!$                    s0=(dzbs)
!!$                    s1=(dzbs+dz)
!!$                    s2=(dzbs+two*dz)
!!$                    s3=(dzbs+three*dz)
!!$                endif
!!$                if (iextp.eq.1) then !1st order
!!$                    aenscz_s(1,k,j)=one
!!$                    benscz_s(1,k,j)=zero
!!$                    censcz_s(1,k,j)=zero
!!$                    denscz_s(1,k,j)=zero
!!$                elseif (iextp.eq.2) then !2nd order
!!$                    aenscz_s(1,k,j)=s1/(s1-s0)
!!$                    benscz_s(1,k,j)=-(s0)/(s1-s0)
!!$                    censcz_s(1,k,j)=zero
!!$                    denscz_s(1,k,j)=zero
!!$                elseif (iextp.eq.3) then !3rd order
!!$                    aenscz_s(1,k,j)=s1*s2/((s1-s0)*s2-s0*s1+s0*s0)
!!$                    benscz_s(1,k,j)=-(s0*s2)/((s1-s0)*s2-s1*s1+s0*s1)
!!$                    censcz_s(1,k,j)=(s0*s1)/(s2*s2-s2*(s0+s1)+s0*s1)
!!$                    denscz_s(1,k,j)=zero
!!$                elseif (iextp.eq.4) then !4th order
!!$                    aenscz_s(1,k,j)=(s1*s2*s3)/&
!!$                                    ((s2*(s1-s0)-s0*s1+s0*s0)*s3+(s0*s0-s0*s1)*s2+s0*s0*s1-s0*s0*s0)
!!$                    benscz_s(1,k,j)=-(s0*s2*s3)/&
!!$                                     ((s2*(s1-s0)-s1*s1+s0*s1)*s3+(s0*s1-s1*s1)*s2+s1*s1*s1-s0*s1*s1)
!!$                    censcz_s(1,k,j)=(s0*s1*s3)/&
!!$                                    ((s2*s2-s2*(s1+s0)+s0*s1)*s3-s2*s2*s2+s2*s2*(s1+s0)-s0*s1*s2)
!!$                    denscz_s(1,k,j)=-(s0*s1*s2)/&
!!$                                     (s3*s3*s3-s3*s3*(s0+s1+s2)+((s1+s0)*s2+s0*s1)*s3-s0*s1*s2)
!!$                endif
!!$
!!$            !2ND IMMERSED BOUNDARY
!!$                kz=(zf(k,i,j)+dz)/dz+1
!!$                dzb=abs(((kz-1)-1)*dz-zf(k,i,j))    !distance of closest node without skip
!!$                dzbs=abs(((kz-2)-1)*dz-zf(k,i,j))   !distance of closest node with skip
!!$                !i) First derivative: with skip
!!$                r1=(dzbs+   dz)/(dzbs)
!!$                r2=(dzbs+two*dz)/(dzbs)
!!$                r3=(dzbs+three*dz)/(dzbs)
!!$                anscz_s(2,k,j)=-((r1+one)*r2+r1)/(r1*r2*dzbs)
!!$                bnscz_s(2,k,j)=(r1*r2)/(((r1-one)*r2-r1+one)*dzbs)
!!$                cnscz_s(2,k,j)=-(r2)/((((r1**two-r1)*r2)-(r1**three)+(r1**two))*dzbs)
!!$                dnscz_s(2,k,j)=(r1)/(((r2**three)+((-r1-one)*r2**two)+(r1*r2))*dzbs) 
!!$                !ii) Extrapolation: with or without skip
!!$                if (iskip_e.eq.0) then !without skip 
!!$                    s0=(dzb)
!!$                    s1=(dzb+dz)
!!$                    s2=(dzb+two*dz)
!!$                    s3=(dzb+three*dz)
!!$                else                   !with skip
!!$                    s0=(dzbs)
!!$                    s1=(dzbs+dz)
!!$                    s2=(dzbs+two*dz)
!!$                    s3=(dzbs+three*dz)
!!$                endif
!!$                if (iextp.eq.1) then !1st order
!!$                    aenscz_s(2,k,j)=one
!!$                    benscz_s(2,k,j)=zero
!!$                    censcz_s(2,k,j)=zero
!!$                    denscz_s(2,k,j)=zero
!!$                elseif (iextp.eq.2) then !2nd order
!!$                    aenscz_s(2,k,j)=s1/(s1-s0)
!!$                    benscz_s(2,k,j)=-(s0)/(s1-s0)
!!$                    censcz_s(2,k,j)=zero
!!$                    denscz_s(2,k,j)=zero
!!$                elseif (iextp.eq.3) then !3rd order
!!$                    aenscz_s(2,k,j)=s1*s2/((s1-s0)*s2-s0*s1+s0*s0)
!!$                    benscz_s(2,k,j)=-(s0*s2)/((s1-s0)*s2-s1*s1+s0*s1)
!!$                    censcz_s(2,k,j)=(s0*s1)/(s2*s2-s2*(s0+s1)+s0*s1)
!!$                    denscz_s(2,k,j)=zero
!!$                elseif (iextp.eq.4) then !4th order
!!$                    aenscz_s(2,k,j)=(s1*s2*s3)/&
!!$                                    ((s2*(s1-s0)-s0*s1+s0*s0)*s3+(s0*s0-s0*s1)*s2+s0*s0*s1-s0*s0*s0)
!!$                    benscz_s(2,k,j)=-(s0*s2*s3)/&
!!$                                     ((s2*(s1-s0)-s1*s1+s0*s1)*s3+(s0*s1-s1*s1)*s2+s1*s1*s1-s0*s1*s1)
!!$                    censcz_s(2,k,j)=(s0*s1*s3)/&
!!$                                    ((s2*s2-s2*(s1+s0)+s0*s1)*s3-s2*s2*s2+s2*s2*(s1+s0)-s0*s1*s2)
!!$                    denscz_s(2,k,j)=-(s0*s1*s2)/&
!!$                                     (s3*s3*s3-s3*s3*(s0+s1+s2)+((s1+s0)*s2+s0*s1)*s3-s0*s1*s2)
!!$                endif
!!$        enddo
!!$        endif
!!$
!!$     endif    
!!$  enddo    
!!$
!!$  return
!!$
!!$  end subroutine set_nbc_coefficients
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    !!
!!$    !!  subroutine: phiw_if
!!$    !!      AUTHOR: Rodrigo Vicente Cruz
!!$    !! DESCRIPTION: When Imposed Flux boundary conditions are used:             
!!$    !!              i)Calculation of the tangential heat flux component qtheta  
!!$    !!             ii)Qtheta is then extrapolated to the wall in both transversal
!!$    !!                directions Y and Z. The extrapolated values are then used to
!!$    !!            iii)Determine the phi wall values along Y and Z directions.
!!$    !!                These values enable the estimation of the  Nusselt number and  
!!$    !!                the subsequent reconstruction ensuring normal heat flux equals  
!!$    !!                to unity.
!!$    !!
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !***************************************************************************
!!$  !
!!$  subroutine phiw_if(phi1,ep1,is,iwrite)
!!$  !
!!$  !***************************************************************************
!!$  use decomp_2d
!!$  use decomp_2d_io
!!$  use MPI
!!$  use var, only : ta2,tb2,tc2,td2,ta3,tb3,tc3,td3
!!$  use param, only : zero,one,two,three,four
!!$  use param, only : yly,zlz,dy,dz,dt,itime
!!$  use complex_geometry, only: nobjy,nobjz,yi,yf,zi,zf
!!$  use ibm
!!$  !qn
!!$  !use var, only : qnmean
!!$  use var, only : te1,te2,te3
!!$  use variables, only : ny,nz
!!$  !
!!$  implicit none
!!$  real(mytype),dimension(xsize(1),xsize(2),xsize(3))    :: phi1,ep1
!!$  integer                                               :: is,iwrite
!!$  !LOCALS
!!$  integer                                               :: i,j,k
!!$  real(mytype),dimension(ysize(1),ysize(2),ysize(3))    :: ep2
!!$  real(mytype),dimension(zsize(1),zsize(2),zsize(3))    :: ep3
!!$  integer                                               :: jy,kz   !1st fluid
!!$  integer                                               :: jys,kzs !1st fluid skip
!!$  real(mytype)                                          :: r,yc,zc
!!$  real(mytype)                                          :: ym,zm,theta  !mesh
!!$  real(mytype)                                          :: yw,zw,thetaw !wall
!!$  real(mytype)                                          :: wfy,wfz !wall fluxes
!!$  real(mytype)                                          :: phiwm   !<phi>x,theta
!!$  real(mytype)                                          :: phiwmt
!!$  real(mytype)                                          :: county,countz
!!$  real(mytype)                                          :: countyt,countzt
!!$  integer                                               :: ifile,code
!!$  character(len=40)                                     :: filename
!!$  !
!!$355 format('Nuw_',I2.2,'.dat')
!!$
!!$    ifile=30+is
!!$    if (iwrite.ne.0) then
!!$        if (itime.eq.ifirst.and.nrank.eq.0) then
!!$            write(filename,355) is
!!$            open(ifile,file=filename,status='unknown')
!!$        endif
!!$    endif
!!$
!!$  yc=yly/two
!!$  zc=zlz/two
!!$  !------------------------------------------------------------------------!
!!$  ! The Y and Z derivatives used for the estimation of the tangential flux !
!!$  ! are calculated with a conditional scheme: for the borders of the inner !
!!$  ! fluid zone, a non-centred 2nd order scheme is used                     !
!!$  !------------------------------------------------------------------------!
!!$  call transpose_x_to_y(phi1,ta2)
!!$  call transpose_x_to_y(ep1,ep2)
!!$  !================================= Y-PENCIL ==============================
!!$  !
!!$  !!Pre-reconstruction
!!$  !if (itime.gt.ifirst) then
!!$  !    if (ider.eq.1.and.iibm.eq.2) call nbclagpoly(ta2,is)
!!$  !endif
!!$
!!$  tb2(:,:,:)=zero
!!$  do i=1,ysize(1)
!!$      do k=1,ysize(3)
!!$          zm=real(ystart(3)-1+k-1,mytype)*dz-zc
!!$          do j=1,ysize(2)
!!$              ym=real(ystart(2)-1+j-1,mytype)*dy-yc
!!$              r=sqrt(ym*ym+zm*zm)    
!!$              if (r.le.ra.and.ep2(i,j,k).eq.0) then
!!$
!!$                  !if (ider.eq.1) then !Derivatives with centred second-order scheme
!!$                  !    tb2(i,j,k) = (ta2(i,j+1,k)-ta2(i,j-1,k))/(two*dy)
!!$                  !elseif (ider.eq.2) then !Derivatives with conditional second-order scheme
!!$                      tb2(i,j,k) = (one-ep2(i,j,k))*&
!!$                                   (ep2(i,j-1,k)*(one-ep2(i,j+2,k))*(one-ep2(i,j+1,k))*&
!!$                                   ((-three*ta2(i,j,k)+four*ta2(i,j+1,k)-one*ta2(i,j+2,k))/(two*dy))&
!!$                                  +ep2(i,j-1,k)*ep2(i,j+2,k)*(one-ep2(i,j+1,k))*&
!!$                                   ((-one*ta2(i,j,k)+one*ta2(i,j+1,k))/(dy))&
!!$                                  +ep2(i,j+1,k)*(one-ep2(i,j-2,k))*(one-ep2(i,j-1,k))*&
!!$                                   ((three*ta2(i,j,k)-four*ta2(i,j-1,k)+one*ta2(i,j-2,k))/(two*dy))&
!!$                                  +ep2(i,j+1,k)*ep2(i,j-2,k)*(one-ep2(i,j-1,k))*&
!!$                                   ((one*ta2(i,j,k)-one*ta2(i,j-1,k))/(dy))&
!!$                                  +(one-ep2(i,j+1,k))*(one-ep2(i,j-1,k))*&
!!$                                   ((ta2(i,j+1,k)-ta2(i,j-1,k))/(two*dy)))
!!$                  !endif
!!$              endif
!!$          enddo
!!$      enddo
!!$  enddo
!!$  !
!!$  call transpose_y_to_z(ta2,ta3) !phi
!!$  call transpose_y_to_z(tb2,tc3) !dy(phi)
!!$  call transpose_y_to_z(ep2,ep3)
!!$  !
!!$  !================================= Z-PENCIL ==============================
!!$  !
!!$  !!Pre-reconstruction
!!$  !if (itime.gt.ifirst) then
!!$  !    if (ider.eq.1.and.iibm.eq.2) call nbclagpolz(ta3,is)
!!$  !endif
!!$  !
!!$  tb3(:,:,:)=zero
!!$  td3(:,:,:)=zero
!!$  phiw3(:,:,:,:,is)=zero
!!$  countz=zero
!!$  phiwm=zero
!!$  !
!!$  do i=1,zsize(1)
!!$      do j=1,zsize(2)
!!$          ym=real(zstart(2)-1+j-1,mytype)*dy-yc
!!$  
!!$          !1. Qtheta
!!$          do k=1,zsize(3)
!!$              zm=real(zstart(3)-1+k-1,mytype)*dz-zc
!!$              r=sqrt(ym*ym+zm*zm)    
!!$              if (r.le.ra.and.ep3(i,j,k).eq.0) then
!!$                  !if (ider.eq.1) then !Derivatives with centred second-order scheme
!!$                  !    tb3(i,j,k) = (ta3(i,j,k+1)-ta3(i,j,k-1))/(two*dz)
!!$                  !elseif (ider.eq.2) then !Derivatives with conditional second-order scheme
!!$                      tb3(i,j,k) = (one-ep3(i,j,k))*&
!!$                                   (ep3(i,j,k-1)*(one-ep3(i,j,k+2))*(one-ep3(i,j,k+1))*&
!!$                                   ((-three*ta3(i,j,k)+four*ta3(i,j,k+1)-one*ta3(i,j,k+2))/(two*dz))&
!!$                                  +ep3(i,j,k-1)*ep3(i,j,k+2)*(one-ep3(i,j,k+1))*&
!!$                                   ((-one*ta3(i,j,k)+one*ta3(i,j,k+1))/(dz))&
!!$                                  +ep3(i,j,k+1)*(one-ep3(i,j,k-2))*(one-ep3(i,j,k-1))*&
!!$                                   ((three*ta3(i,j,k)-four*ta3(i,j,k-1)+one*ta3(i,j,k-2))/(two*dz))&
!!$                                  +ep3(i,j,k+1)*ep3(i,j,k-2)*(one-ep3(i,j,k-1))*&
!!$                                   ((one*ta3(i,j,k)-one*ta3(i,j,k-1))/(dz))&
!!$                                  +(one-ep3(i,j,k+1))*(one-ep3(i,j,k-1))*&
!!$                                   ((ta3(i,j,k+1)-ta3(i,j,k-1))/(two*dz)))
!!$                  !endif
!!$  
!!$                  !Heat fluxes in normal and tangential directions (fluid region)
!!$                  theta=atan2(zm,ym)                                         !theta dans [-pi,pi]
!!$                  td3(i,j,k)=-dsin(theta)*tc3(i,j,k)+dcos(theta)*tb3(i,j,k)  !qtheta
!!$                  te3(i,j,k)=-dcos(theta)*tc3(i,j,k)-dsin(theta)*tb3(i,j,k)  !qn
!!$                  !!====DEBUG
!!$                  !if (zstart(2)+j-1.eq.ny/2+1.and.zstart(1)+i-1.eq.nx/2+1) then
!!$                  !    print*, 'qtheta_z=',zm,td3(i,j,k)
!!$                  !endif
!!$                  !!====
!!$              endif
!!$          enddo
!!$  
!!$          !2. Z-direction wall treatment
!!$          if (nobjz(i,j).eq.2) then !Immersed objects
!!$              do k=1,nobjz(i,j)
!!$                  if (k.eq.1) then !1st object
!!$                      kz=(zf(k,i,j)+dz)/dz+1 !without skip
!!$                      !Extrapolation of the wall-tangential heat flux
!!$                        !!$if (iskip_e.eq.0) then !Without skip
!!$                        !!$qthetaw3(2,k,i,j,is)=aenscz(2,k,j)*td3(i,j,kz  )+benscz(2,k,j)*td3(i,j,kz+1)+&
!!$                        !!$                    censcz(2,k,j)*td3(i,j,kz+2)+denscz(2,k,j)*td3(i,j,kz+3)
!!$                        !!$else                   !With skip
!!$                        !!$qthetaw3(2,k,i,j,is)=aenscz(2,k,j)*td3(i,j,kz+1)+benscz(2,k,j)*td3(i,j,kz+2)+&
!!$                        !!$                    censcz(2,k,j)*td3(i,j,kz+3)+denscz(2,k,j)*td3(i,j,kz+4)
!!$                        !!$endif
!!$                        !!!Fourth-order
!!$                        !!qthetaw3(2,k,i,j,is)=aenscz(2,k,j)*td3(i,j,kz  )+benscz(2,k,j)*td3(i,j,kz+1)+&
!!$                        !!                    censcz(2,k,j)*td3(i,j,kz+2)+denscz(2,k,j)*td3(i,j,kz+3)
!!$                        !!!First order with skip
!!$                        !!qthetaw3(2,k,i,j,is)=td3(i,j,kz+1)
!!$                        !First order without skip
!!$                        qthetaw3(2,k,i,j,is)=td3(i,j,kz)
!!$                      
!!$
!!$                      !Temperature at the wall (with O3 non-centred)
!!$                      kzs=kz+1 !with skip
!!$                      zw=zf(k,i,j)-zc
!!$                      thetaw=atan2(zw,ym)
!!$                      wfz=dcos(thetaw)*qthetaw3(2,k,i,j,is)-dsin(thetaw)*one !wall flux in Z
!!$                      phiw3(2,k,i,j,is)=(-bnscz(2,k,j)*ta3(i,j,kzs  )-cnscz(2,k,j)*ta3(i,j,kzs+1)&
!!$                                         -dnscz(2,k,j)*ta3(i,j,kzs+2)+wfz)/anscz(2,k,j)
!!$                      phiwm=phiwm+phiw3(2,k,i,j,is) !averaged value <phi>x,theta for nusselt estimation
!!$                      countz=countz+1.
!!$                  elseif (k.eq.2) then !2nd object
!!$                      kz=zi(k,i,j)/dz+1
!!$                      !====DEBUG
!!$                      !if (abs(zp(kz)-zi(k,i,j)).lt.tol) kz=kz-1 !use if exact_ib_pipe
!!$                      !Extrapolation of the wall-tangential heat flux
!!$                        !!$if (iskip_e.eq.0) then !Without skip
!!$                        !!$qthetaw3(1,k,i,j,is)=aenscz(1,k,j)*td3(i,j,kz  )+benscz(1,k,j)*td3(i,j,kz-1)+&
!!$                        !!$                    censcz(1,k,j)*td3(i,j,kz-2)+denscz(1,k,j)*td3(i,j,kz-3)
!!$                        !!$else                   !With skip
!!$                        !!$qthetaw3(1,k,i,j,is)=aenscz(1,k,j)*td3(i,j,kz-1)+benscz(1,k,j)*td3(i,j,kz-2)+&
!!$                        !!$                    censcz(1,k,j)*td3(i,j,kz-3)+denscz(1,k,j)*td3(i,j,kz-4)
!!$                        !!$endif
!!$                        !!!Fourth-order
!!$                        !!qthetaw3(1,k,i,j,is)=aenscz(1,k,j)*td3(i,j,kz  )+benscz(1,k,j)*td3(i,j,kz-1)+&
!!$                        !!                    censcz(1,k,j)*td3(i,j,kz-2)+denscz(1,k,j)*td3(i,j,kz-3)
!!$                        !!!First-order with skip
!!$                        !!qthetaw3(1,k,i,j,is)=td3(i,j,kz-1)
!!$                        !First order without skip
!!$                        qthetaw3(1,k,i,j,is)=td3(i,j,kz)
!!$  
!!$                      !Temperature at the wall (with O3 non-centred)
!!$                      kzs=kz-1 !with skip
!!$                      zw=zi(k,i,j)-zc
!!$                      thetaw=atan2(zw,ym)
!!$                      wfz=dcos(thetaw)*qthetaw3(1,k,i,j,is)-dsin(thetaw)*one !wall flux in Z
!!$                      phiw3(1,k,i,j,is)=(-bnscz(1,k,j)*ta3(i,j,kzs  )-cnscz(1,k,j)*ta3(i,j,kzs-1)&
!!$                                         -dnscz(1,k,j)*ta3(i,j,kzs-2)-wfz)/anscz(1,k,j)
!!$                      phiwm=phiwm+phiw3(1,k,i,j,is) !averaged value <phi>x,theta for nusselt estimation
!!$                      countz=countz+1.
!!$                  endif
!!$              enddo
!!$          endif
!!$      enddo
!!$  enddo
!!$  call MPI_ALLREDUCE(countz,countzt,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
!!$  !
!!$  call transpose_z_to_y(td3,td2) !qtheta
!!$  !
!!$  !================================= Y-PENCIL ==============================
!!$  !
!!$  phiw2(:,:,:,:,is)=zero
!!$  county=0
!!$  !
!!$  do i=1,ysize(1)
!!$      do k=1,ysize(3)
!!$          zm=real(ystart(3)-1+k-1,mytype)*dz-zc
!!$  
!!$          !1. Y-direction wall treatment
!!$          if (nobjy(i,k).eq.2) then !Immersed objects
!!$              do j=1,nobjy(i,k)
!!$                  if (j.eq.1) then !1st object
!!$                      jy=1!jy=(yf(j,i,k)+dy)/dy+1
!!$                      !====DEBUG
!!$                      !do while(yp(jy).le.yf(j,i,k)) !use if exact_ib_pipe
!!$                      do while(yp(jy).lt.yf(j,i,k))
!!$                         jy=jy+1
!!$                      enddo
!!$                      !Extrapolation of the wall-tangential heat flux
!!$                        !!$if (iskip_e.eq.0) then !Without skip
!!$                        !!$qthetaw2(2,j,i,k,is)=aenscy(2,j,k)*td2(i,jy  ,k)+benscy(2,j,k)*td2(i,jy+1,k)+&
!!$                        !!$                    censcy(2,j,k)*td2(i,jy+2,k)+denscy(2,j,k)*td2(i,jy+3,k)
!!$                        !!$else                   !With skip
!!$                        !!$qthetaw2(2,j,i,k,is)=aenscy(2,j,k)*td2(i,jy+1,k)+benscy(2,j,k)*td2(i,jy+2,k)+&
!!$                        !!$                    censcy(2,j,k)*td2(i,jy+3,k)+denscy(2,j,k)*td2(i,jy+4,k)
!!$                        !!$endif
!!$                        !!!Fourth-order
!!$                        !!qthetaw2(2,j,i,k,is)=aenscy(2,j,k)*td2(i,jy  ,k)+benscy(2,j,k)*td2(i,jy+1,k)+&
!!$                        !!                    censcy(2,j,k)*td2(i,jy+2,k)+denscy(2,j,k)*td2(i,jy+3,k)
!!$                        !!!First-order with skip
!!$                        !!qthetaw2(2,j,i,k,is)=td2(i,jy+1,k)
!!$                        !First-order without skip
!!$                        qthetaw2(2,j,i,k,is)=td2(i,jy,k)
!!$  
!!$                      !Temperature at the wall (with O3 non-centred)
!!$                      jys=jy+1 !with skip
!!$                      yw=yf(j,i,k)-yc
!!$                      thetaw=atan2(zm,yw)
!!$                      wfy=-dsin(thetaw)*qthetaw2(2,j,i,k,is)-dcos(thetaw)*one !wall flux in Y
!!$                      phiw2(2,j,i,k,is)=(-bnscy(2,j,k)*ta2(i,jys  ,k)-cnscy(2,j,k)*ta2(i,jys+1,k)&
!!$                                       -dnscy(2,j,k)*ta2(i,jys+2,k)+wfy)/anscy(2,j,k)
!!$                      phiwm=phiwm+phiw2(2,j,i,k,is) !averaged value <phi>x,theta for nusselt estimation
!!$                      county=county+1.
!!$                  elseif (j.eq.2) then !2nd object
!!$                      jy=1!jy=yi(j,i,k)/dy+1
!!$                      do while(yp(jy).lt.yi(j,i,k))
!!$                         jy=jy+1
!!$                      enddo
!!$                      jy=jy-1
!!$                      !Extrapolation of the wall-tangential heat flux
!!$                        !!$if (iskip_e.eq.0) then !Without skip
!!$                        !!$qthetaw2(1,j,i,k,is)=aenscy(1,j,k)*td2(i,jy  ,k)+benscy(1,j,k)*td2(i,jy-1,k)+&
!!$                        !!$                    censcy(1,j,k)*td2(i,jy-2,k)+denscy(1,j,k)*td2(i,jy-3,k)
!!$                        !!$else                   !With skip
!!$                        !!$qthetaw2(1,j,i,k,is)=aenscy(1,j,k)*td2(i,jy-1,k)+benscy(1,j,k)*td2(i,jy-2,k)+&
!!$                        !!$                    censcy(1,j,k)*td2(i,jy-3,k)+denscy(1,j,k)*td2(i,jy-4,k)
!!$                        !!$endif
!!$                        !!!Fourth-order
!!$                        !!qthetaw2(1,j,i,k,is)=aenscy(1,j,k)*td2(i,jy  ,k)+benscy(1,j,k)*td2(i,jy-1,k)+&
!!$                        !!                    censcy(1,j,k)*td2(i,jy-2,k)+denscy(1,j,k)*td2(i,jy-3,k)
!!$                        !!!First-order with skip
!!$                        !!qthetaw2(1,j,i,k,is)=td2(i,jy-1,k)
!!$                        !First-order without skip
!!$                        qthetaw2(1,j,i,k,is)=td2(i,jy,k)
!!$  
!!$                      !Temperature at the wall (with O3 non-centred)
!!$                      jys=jy-1 !with skip
!!$                      yw=yi(j,i,k)-yc
!!$                      thetaw=atan2(zm,yw)
!!$                      wfy=-dsin(thetaw)*qthetaw2(1,j,i,k,is)-dcos(thetaw)*one !wall flux in Y
!!$                      phiw2(1,j,i,k,is)=(-bnscy(1,j,k)*ta2(i,jys  ,k)-cnscy(1,j,k)*ta2(i,jys-1,k)&
!!$                                       -dnscy(1,j,k)*ta2(i,jys-2,k)-wfy)/anscy(1,j,k)
!!$                      phiwm=phiwm+phiw2(1,j,i,k,is) !averaged value <phi>x,theta for nusselt estimation
!!$                      county=county+1.
!!$                  endif
!!$              enddo
!!$          endif
!!$  
!!$      enddo
!!$  enddo
!!$  call MPI_ALLREDUCE(county,countyt,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
!!$  
!!$  call MPI_ALLREDUCE(phiwm,phiwmt,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
!!$  phiwm=phiwmt/(countyt+countzt)
!!$  nuw(is)=-one/phiwm                !Nu(t)=-1/<phi>x,theta (time average to be performed afterwards
!!$  !if (nrank==0) then
!!$  if (nrank.eq.0.and.iwrite.ne.0) then
!!$     !print *,'Nu(t)       =       ',nuw(is)
!!$     write(ifile,*) real((itime-1)*dt,mytype), nuw(is)
!!$  endif
!!$  !For reconstruction smoothness
!!$  do k=1,xsize(3)
!!$      zm=dz*real(xstart(3)-1+k-1,mytype)-zc 
!!$      do j=1,xsize(2)
!!$          if (istret.eq.0) ym=real(j+xstart(2)-1-1,mytype)*dy-yc
!!$          if (istret.ne.0) ym=yp(j+xstart(2)-1)-yc
!!$          r=sqrt(ym*ym+zm*zm)    
!!$          do i=1,xsize(1)
!!$              if (r.ge.ra.or.ep1(i,j,k).eq.1) then
!!$                  phi1(i,j,k)=-one/nuw(is)
!!$              endif
!!$          enddo
!!$      enddo
!!$  enddo
!!$  
!!$  return
!!$  
!!$  end subroutine phiw_if
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    !!
!!$    !!  subroutine: phiw_cht
!!$    !!      AUTHOR: Rodrigo Vicente Cruz
!!$    !! DESCRIPTION: Predicts wall temperature to be imposed at the wall for
!!$    !!              fluid temperature solution (iflag=1) and solid temperature
!!$    !!              solution. Contains:
!!$    !!              - Prediction of phi^(n+1): Dirichlet (iflag=1)
!!$    !!              - Prediction/extrapolation of tangential heat-flux component
!!$    !!              - Prediction of phis^(n+1): Neumann (iflag=2)
!!$    !!
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !***************************************************************************
!!$  !
!!$  subroutine phiw_cht(phi1,phis1,ep1,is,icht,iflag,iwrite)
!!$  !
!!$  !***************************************************************************
!!$  use decomp_2d
!!$  use decomp_2d_io
!!$  use MPI
!!$  use var, only : ta2,tb2,tc2,td2,tg2,th2,ta3,tb3,tc3,td3,te3,tf3,tg3,th3
!!$  use var, only : di2,di3
!!$  use param, only : zero,one,two,three,four
!!$  use param, only : yly,zlz,dy,dz,dt,itime
!!$  use complex_geometry, only: nobjy,nobjz,yi,yf,zi,zf
!!$  use ibm
!!$  !
!!$  implicit none
!!$  real(mytype),dimension(xsize(1),xsize(2),xsize(3))    :: phi1,phis1,ep1
!!$  integer                                               :: icht,is,iwrite
!!$  integer                                               :: iflag ! 1: fluid solver | 2: solid solver
!!$  !LOCALS
!!$  integer                                               :: i,j,k
!!$  !real(mytype),dimension(ysize(1),ysize(2),ysize(3))    :: ep2
!!$  !real(mytype),dimension(zsize(1),zsize(2),zsize(3))    :: ep3
!!$  integer                                               :: jy,kz   !1st fluid
!!$  real(mytype)                                          :: r,yc,zc
!!$  real(mytype)                                          :: ym,zm       !mesh
!!$  real(mytype)                                          :: yw,zw,theta !wall
!!$  real(mytype)                                          :: wfy,wfz !wall fluxes
!!$  real(mytype)                                          :: phiwm   !<phi>x,theta
!!$  real(mytype)                                          :: phiwmt
!!$  real(mytype)                                          :: county,countz
!!$  real(mytype)                                          :: countyt,countzt
!!$  integer                                               :: ifile,code,irank,iibm_save
!!$  character(len=40)                                     :: filename
!!$  !
!!$
!!$355 format('Nuw_',I2.2,'.dat')
!!$
!!$    ifile=30+is
!!$    if (iwrite.ne.0) then
!!$        if (itime.eq.ifirst.and.iflag.eq.1.and.nrank.eq.0) then
!!$            write(filename,355) is
!!$            open(ifile,file=filename,status='unknown')
!!$        endif
!!$    endif
!!$
!!$    yc=yly/two
!!$    zc=zlz/two
!!$
!!$    call transpose_x_to_y(phi1,ta2)
!!$    call transpose_x_to_y(phis1,tb2)
!!$    !call transpose_x_to_y(ep1,ep2)
!!$
!!$    !Tangential flux treatment
!!$    if (iflag.eq.2) then
!!$        !To prevent Dirichlet reconstruction
!!$        iibm_save=iibm
!!$        iibm=0
!!$
!!$        !Y-PENCILS
!!$        !Fluid: phi^(n+1)
!!$        call chtlagpoly(ta2,is)
!!$        call deryS (tc2,ta2(:,:,:),di2,sy,ffypS,fsypS,fwypS,ppy,ysize(1),ysize(2),ysize(3),1)
!!$        call transpose_y_to_z(tc2,tc3) !dyphi
!!$        !Solid: phis^(n)
!!$        call chtlagpoly_s(tb2,icht)
!!$        call deryS (td2,tb2(:,:,:),di2,sy,ffypS,fsypS,fwypS,ppy,ysize(1),ysize(2),ysize(3),1)
!!$        call transpose_y_to_z(td2,td3) !dyphis
!!$
!!$        !Z-PENCILS
!!$        call transpose_y_to_z(ta2,ta3) !phi
!!$        call transpose_y_to_z(tb2,tb3) !phis
!!$        !Fluid: phi^(n+1)
!!$        call chtlagpolz(ta3,is)
!!$        call derzS (te3,ta3(:,:,:),di3,sz,ffzpS,fszpS,fwzpS,zsize(1),zsize(2),zsize(3),1)
!!$        !Solid: phis^(n)
!!$        call chtlagpolz_s(tb3,icht)
!!$        call derzS (tf3,tb3(:,:,:),di3,sz,ffzpS,fszpS,fwzpS,zsize(1),zsize(2),zsize(3),1)
!!$        iibm=iibm_save
!!$        !!-------------------------
!!$        !!tc3: dyphi
!!$        !!td3: dyphis
!!$        !!te3: dzphi
!!$        !!tf3: dzphis
!!$        !!-------------------------
!!$        !Compute d(phis)/dtheta^(n+1)
!!$        call derphiw_cht(tc3,td3,te3,tf3,is,icht,1)
!!$    endif
!!$
!!$    !=========================================== Y-PENCIL ========================================
!!$    !
!!$    if (iflag.eq.1) then
!!$        !====DEBUG
!!$        !phiw2(:,:,:,:,is)=zero
!!$        !phiws2(:,:,:,:,icht)=zero
!!$        phiwm=zero
!!$        county=zero
!!$    endif
!!$    !
!!$    !
!!$    do i=1,ysize(1)
!!$        do k=1,ysize(3)
!!$            if (nobjy(i,k).ne.0) then !Immersed objects
!!$                do j=1,nobjy(i,k)
!!$                    if (j.eq.1) then !1st object
!!$                        !1ST IMMERSED BOUNDARY
!!$                        jy=1!jy=yi(j,i,k)/dy+1
!!$                        do while(yp(jy).lt.yi(j,i,k))
!!$                           jy=jy+1
!!$                        enddo
!!$                        jy=jy-1
!!$                        if (iflag.eq.1) then !Fluid wall temperature n+1
!!$                            !=====DEBUG
!!$                            !!IF EXACT IB PIPE: FORCE NODES AT IB
!!$                            !if (abs(yp(jy+1)-yi(j,i,k)).lt.tol.and.itime.ne.ifirst)&
!!$                            !    tb2(i,jy+1,k)=phiws2(1,j,i,k,icht)
!!$                            !==========
!!$                            !solid extrapolation
!!$                            phiws2(1,j,i,k,icht)=aenscy_s(1,j,k)*tb2(i,jy+1,k)+benscy_s(1,j,k)*tb2(i,jy+2,k)+&
!!$                                                 censcy_s(1,j,k)*tb2(i,jy+3,k)+denscy_s(1,j,k)*tb2(i,jy+4,k)
!!$                            !fluid (equal)
!!$                            phiw2(1,j,i,k,is)=phiws2(1,j,i,k,icht)
!!$
!!$                        elseif (iflag.eq.2) then !Solid wall temperature n+1
!!$                            !isoflux (IF) qn=1/g2*ra/rao
!!$                            yw=yi(j,i,k)-yc
!!$                            zw=real(ystart(3)-1+k-1,mytype)*dz-zc
!!$                            theta=atan2(zw,yw)
!!$                            !wfy=-dcos(theta)*one
!!$                            wfy=-dcos(theta)*(ra/rao)
!!$                            wfy=wfy/g2(is)
!!$                            !forward scheme
!!$                            phiws2(1,j,i,k,icht)=(-bnscy_s(1,j,k)*tb2(i,jy+2,k)-cnscy_s(1,j,k)*tb2(i,jy+3,k)&
!!$                                                  -dnscy_s(1,j,k)*tb2(i,jy+4,k)+wfy)/anscy_s(1,j,k)
!!$                        endif
!!$
!!$                        !2ND IMMERSED BOUNDARY
!!$                        jy=1!jy=(yf(j,i,k)+dy)/dy+1
!!$                        !====DEBUG
!!$                        !do while(yp(jy).le.yf(j,i,k)) !use if exact_ib_pipe
!!$                        do while(yp(jy).lt.yf(j,i,k))
!!$                           jy=jy+1
!!$                        enddo
!!$                        if (iflag.eq.1) then !Fluid wall temperature n+1
!!$                            !=====DEBUG
!!$                            !!IF EXACT IB PIPE: FORCE NODES AT IB
!!$                            !if (abs(yp(jy-1)-yf(j,i,k)).lt.tol.and.itime.ne.ifirst)&
!!$                            !    tb2(i,jy-1,k)=phiws2(2,j,i,k,icht)
!!$                            !==========
!!$                            !fluid extrapolation
!!$                            phiw2(2,j,i,k,is)=aenscy(2,j,k)*ta2(i,jy  ,k)+benscy(2,j,k)*ta2(i,jy+1,k)+&
!!$                                              censcy(2,j,k)*ta2(i,jy+2,k)+denscy(2,j,k)*ta2(i,jy+3,k)
!!$                            !solid extrapolation
!!$                            phiws2(2,j,i,k,icht)=aenscy_s(2,j,k)*tb2(i,jy-1,k)+benscy_s(2,j,k)*tb2(i,jy-2,k)+&
!!$                                                 censcy_s(2,j,k)*tb2(i,jy-3,k)+denscy_s(2,j,k)*tb2(i,jy-4,k)
!!$                            !fluid: half-sum
!!$                            phiw2(2,j,i,k,is)=half*(phiw2(2,j,i,k,is)+phiws2(2,j,i,k,icht))
!!$                            !!====DEBUG
!!$                            !do irank=0,p_row*p_col-1
!!$                            !  if(nrank.eq.irank) then
!!$                            !    if (ystart(1)+i-1.eq.nx/2+1) print*,'phiw2:' , real((ystart(3)+k-1-1)*dz-zc,4),real(yf(j,i,k),4),phiw2(2,j,i,k,is),is
!!$                            !    if (ystart(1)+i-1.eq.nx/2+1) print*,'phiws2:', real((ystart(3)+k-1-1)*dz-zc,4),real(yf(j,i,k),4),phiws2(2,j,i,k,icht),is
!!$                            !  endif
!!$                            !enddo
!!$                            !!=========
!!$                            if (nobjy(i,k).eq.2) then !average wall temperature 
!!$                                phiwm=phiwm+phiw2(2,j,i,k,is)
!!$                               county=county+1.
!!$                            endif
!!$                        elseif (iflag.eq.2) then !Solid wall temperature n+1
!!$                            if (nobjy(i,k).eq.1) then !only one object
!!$                                !isoflux (IF) qn=1/g2*ra/rao
!!$                                yw=yf(j,i,k)-yc
!!$                                zw=real(ystart(3)-1+k-1,mytype)*dz-zc
!!$                                theta=atan2(zw,yw)
!!$                                !wfy=-dcos(theta)*one
!!$                                wfy=-dcos(theta)*(ra/rao)
!!$                                wfy=wfy/g2(is)
!!$                                !backward scheme
!!$                                phiws2(2,j,i,k,icht)=(-bnscy_s(2,j,k)*tb2(i,jy-2,k)-cnscy_s(2,j,k)*tb2(i,jy-3,k)&
!!$                                                      -dnscy_s(2,j,k)*tb2(i,jy-4,k)-wfy)/anscy_s(2,j,k)
!!$                            else
!!$                                !extrapolate new fluid wall value
!!$                                phiw2(2,j,i,k,is)=aenscy(2,j,k)*ta2(i,jy  ,k)+benscy(2,j,k)*ta2(i,jy+1,k)+&
!!$                                                  censcy(2,j,k)*ta2(i,jy+2,k)+denscy(2,j,k)*ta2(i,jy+3,k)
!!$                                !calculate wall derivative from fluid (forward scheme)
!!$                                wfy=anscy(2,j,k)*phiw2(2,j,i,k,is)+bnscy(2,j,k)*ta2(i,jy+1,k)+&
!!$                                    cnscy(2,j,k)*ta2(i,jy+2,k)    +dnscy(2,j,k)*ta2(i,jy+3,k)
!!$                                wfy=wfy/g2(is)
!!$                                !add tangential flux contribution
!!$                                yw=yf(j,i,k)-yc
!!$                                zw=real(ystart(3)-1+k-1,mytype)*dz-zc
!!$                                theta=atan2(zw,yw)
!!$                                wfy=wfy-dsin(theta)*((g2(is)-one)/g2(is))*qthetaws2(2,j,i,k,icht)
!!$                                !!====DEBUG
!!$                                !if (ystart(1)+i-1.eq.nx/2+1) print*,'wfy', theta, wfy, is
!!$                                !!=========
!!$                                !impose wall derivative on solid (backward scheme)
!!$                                phiws2(2,j,i,k,icht)=(-bnscy_s(2,j,k)*tb2(i,jy-2,k)-cnscy_s(2,j,k)*tb2(i,jy-3,k)&
!!$                                                      -dnscy_s(2,j,k)*tb2(i,jy-4,k)-wfy)/anscy_s(2,j,k)
!!$                                !!====DEBUG
!!$                                !if (ystart(1)+i-1.eq.nx/2+1) print*,'phiws2', theta, phiws2(2,j,i,k,icht), is,itime
!!$                                !!=========
!!$                            endif
!!$                        endif
!!$
!!$                    elseif (j.eq.2) then !2nd object
!!$                        !1ST IMMERSED BOUNDARY
!!$                        jy=1!jy=yi(j,i,k)/dy+1
!!$                        do while(yp(jy).lt.yi(j,i,k))
!!$                           jy=jy+1
!!$                        enddo
!!$                        jy=jy-1
!!$                        if (iflag.eq.1) then !Fluid wall temperature n+1
!!$                            !=====DEBUG
!!$                            !!IF EXACT IB PIPE: FORCE NODES AT IB
!!$                            !if (abs(yp(jy+1)-yi(j,i,k)).lt.tol.and.itime.ne.ifirst)&
!!$                            !    tb2(i,jy+1,k)=phiws2(1,j,i,k,icht)
!!$                            !==========
!!$                            !fluid extrapolation
!!$                            phiw2(1,j,i,k,is)=aenscy(1,j,k)*ta2(i,jy  ,k)+benscy(1,j,k)*ta2(i,jy-1,k)+&
!!$                                              censcy(1,j,k)*ta2(i,jy-2,k)+denscy(1,j,k)*ta2(i,jy-3,k)
!!$                            !solid extrapolation
!!$                            phiws2(1,j,i,k,icht)=aenscy_s(1,j,k)*tb2(i,jy+1,k)+benscy_s(1,j,k)*tb2(i,jy+2,k)+&
!!$                                                 censcy_s(1,j,k)*tb2(i,jy+3,k)+denscy_s(1,j,k)*tb2(i,jy+4,k)
!!$                            !fluid: half-sum
!!$                            phiw2(1,j,i,k,is)=half*(phiw2(1,j,i,k,is)+phiws2(1,j,i,k,icht))
!!$                            !!====DEBUG
!!$                            !do irank=0,p_row*p_col-1
!!$                            !  if(nrank.eq.irank) then
!!$                            !    if (ystart(1)+i-1.eq.nx/2+1) print*,'phiw2:' , real((ystart(3)+k-1-1)*dz-zc,4),real(yi(j,i,k),4),phiw2(1,j,i,k,is),is
!!$                            !    if (ystart(1)+i-1.eq.nx/2+1) print*,'phiws2:', real((ystart(3)+k-1-1)*dz-zc,4),real(yi(j,i,k),4),phiws2(1,j,i,k,icht),is     
!!$                            !  endif
!!$                            !enddo
!!$                            !!=========
!!$                            if (nobjy(i,k).eq.2) then !average wall temperature 
!!$                                phiwm=phiwm+phiw2(1,j,i,k,is)
!!$                                county=county+1.
!!$                            endif
!!$                        elseif (iflag.eq.2) then  !Solid wall temperature n+1
!!$                            !extrapolate new fluid wall value
!!$                            phiw2(1,j,i,k,is)=aenscy(1,j,k)*ta2(i,jy  ,k)+benscy(1,j,k)*ta2(i,jy-1,k)+&
!!$                                              censcy(1,j,k)*ta2(i,jy-2,k)+denscy(1,j,k)*ta2(i,jy-3,k)
!!$                            !calculate wall derivative from fluid (backward scheme)
!!$                            wfy=anscy(1,j,k)*phiw2(1,j,i,k,is)+bnscy(1,j,k)*ta2(i,jy-1,k)+&
!!$                                cnscy(1,j,k)*ta2(i,jy-2,k)    +dnscy(1,j,k)*ta2(i,jy-3,k)
!!$                            wfy=-wfy/g2(is)
!!$                            !add tangential flux contribution
!!$                            yw=yi(j,i,k)-yc
!!$                            zw=real(ystart(3)-1+k-1,mytype)*dz-zc
!!$                            theta=atan2(zw,yw)
!!$                            wfy=wfy-dsin(theta)*((g2(is)-one)/g2(is))*qthetaws2(1,j,i,k,icht)
!!$                            !!====DEBUG
!!$                            !if (ystart(1)+i-1.eq.nx/2+1) print*,'wfy', theta, wfy, is
!!$                            !!=========
!!$                            !impose wall derivative on solid (forward scheme)
!!$                            phiws2(1,j,i,k,icht)=(-bnscy_s(1,j,k)*tb2(i,jy+2,k)-cnscy_s(1,j,k)*tb2(i,jy+3,k)&
!!$                                                  -dnscy_s(1,j,k)*tb2(i,jy+4,k)+wfy)/anscy_s(1,j,k)
!!$                            !!====DEBUG
!!$                            !if (ystart(1)+i-1.eq.nx/2+1) print*,'phiws2', theta, phiws2(1,j,i,k,icht), is,itime
!!$                            !!=========
!!$                        endif
!!$
!!$                        !2ND IMMERSED BOUNDARY
!!$                        jy=1!jy=(yf(j,i,k)+dy)/dy+1
!!$                        !====DEBUG
!!$                        !do while(yp(jy).le.yf(j,i,k)) !use if exact_ib_pipe
!!$                        do while(yp(jy).lt.yf(j,i,k))
!!$                           jy=jy+1
!!$                        enddo
!!$                        if (iflag.eq.1) then !Fluid wall temperature n+1
!!$                            !=====DEBUG
!!$                            !!IF EXACT IB PIPE: FORCE NODES AT IB
!!$                            !if (abs(yp(jy-1)-yf(j,i,k)).lt.tol.and.itime.ne.ifirst)& 
!!$                            !    tb2(i,jy-1,k)=phiws2(2,j,i,k,icht)
!!$                            !==========
!!$                            !solid extrapolation
!!$                            phiws2(2,j,i,k,icht)=aenscy_s(2,j,k)*tb2(i,jy-1,k)+benscy_s(2,j,k)*tb2(i,jy-2,k)+&
!!$                                                 censcy_s(2,j,k)*tb2(i,jy-3,k)+denscy_s(2,j,k)*tb2(i,jy-4,k)
!!$                            !fluid (equal)
!!$                            phiw2(2,j,i,k,is)=phiws2(2,j,i,k,icht)
!!$
!!$                        elseif (iflag.eq.2) then !Solid wall temperature n+1
!!$                            !isoflux (IF) qn=1/g2*ra/rao
!!$                            yw=yf(j,i,k)-yc
!!$                            zw=real(ystart(3)-1+k-1,mytype)*dz-zc
!!$                            theta=atan2(zw,yw)
!!$                            !wfy=-dcos(theta)*one
!!$                            wfy=-dcos(theta)*(ra/rao)
!!$                            wfy=wfy/g2(is)
!!$                            !backward scheme
!!$                            phiws2(2,j,i,k,icht)=(-bnscy_s(2,j,k)*tb2(i,jy-2,k)-cnscy_s(2,j,k)*tb2(i,jy-3,k)&
!!$                                                  -dnscy_s(2,j,k)*tb2(i,jy-4,k)-wfy)/anscy_s(2,j,k)
!!$                        endif
!!$
!!$                    endif
!!$                enddo
!!$            endif
!!$        enddo
!!$    enddo
!!$
!!$    if (iflag.eq.1) call MPI_ALLREDUCE(county,countyt,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
!!$    call transpose_y_to_z(ta2,ta3)     !phi
!!$    call transpose_y_to_z(tb2,tb3)     !phis
!!$    !call transpose_y_to_z(ep2,ep3)
!!$    !=========================================== Z-PENCIL ========================================
!!$    if (iflag.eq.1) then
!!$        !====DEBUG
!!$        !phiw3(:,:,:,:,is)=zero
!!$        !phiws3(:,:,:,:,icht)=zero
!!$        countz=zero
!!$    endif
!!$    do i=1,zsize(1)
!!$        do j=1,zsize(2)
!!$            !if (nobjz(i,j).eq.2) then !Immersed objects
!!$            if (nobjz(i,j).ne.0) then !Immersed objects
!!$                do k=1,nobjz(i,j)
!!$                    if (k.eq.1) then !1st object
!!$                        !1ST IMMERSED BOUNDARY
!!$                        kz=zi(k,i,j)/dz+1
!!$                        !====DEBUG
!!$                        !if (abs(zp(kz)-zi(k,i,j)).lt.tol) kz=kz-1 !use if exact_ib_pipe
!!$                        if (iflag.eq.1) then !Fluid wall temperature n+1
!!$                            !=====DEBUG
!!$                            !!IF EXACT IB PIPE: FORCE NODES AT IB
!!$                            !if (abs(zp(kz+1)-zi(k,i,j)).lt.tol.and.itime.ne.ifirst)& 
!!$                            !    tb3(i,j,kz+1)=phiws3(1,k,i,j,icht)
!!$                            !==========
!!$                            !solid extrapolation                        
!!$                            phiws3(1,k,i,j,icht)=aenscz_s(1,k,j)*tb3(i,j,kz+1)+benscz_s(1,k,j)*tb3(i,j,kz+2)+&
!!$                                               censcz_s(1,k,j)*tb3(i,j,kz+3)+denscz_s(1,k,j)*tb3(i,j,kz+4)
!!$
!!$                            !fluid (equal)
!!$                            phiw3(1,k,i,j,is)=phiws3(1,k,i,j,icht)
!!$                        elseif (iflag.eq.2) then !Solid wall temperature n+1
!!$                            !isoflux (IF) qn=1/g2*ra/rao
!!$                            yw=real(zstart(2)-1+j-1,mytype)*dy-yc
!!$                            zw=zi(k,i,j)-zc
!!$                            theta=atan2(zw,yw)
!!$                            !wfz=-dsin(theta)*one
!!$                            wfz=-dsin(theta)*(ra/rao)
!!$                            wfz=wfz/g2(is)
!!$                            !forward scheme
!!$                            phiws3(1,k,i,j,icht)=(-bnscz_s(1,k,j)*tb3(i,j,kz+2)-cnscz_s(1,k,j)*tb3(i,j,kz+3)&
!!$                                                -dnscz_s(1,k,j)*tb3(i,j,kz+4)+wfz)/anscz_s(1,k,j)
!!$                        endif
!!$
!!$                        !2ND IMMERSED BOUNDARY
!!$                        kz=(zf(k,i,j)+dz)/dz+1 !without skip
!!$                        if (iflag.eq.1) then !Fluid wall temperature n+1
!!$                            !=====DEBUG
!!$                            !!IF EXACT IB PIPE: FORCE NODES AT IB
!!$                            !if (abs(zp(kz-1)-zf(k,i,j)).lt.tol.and.itime.ne.ifirst)& 
!!$                            !    tb3(i,j,kz-1)=phiws3(2,k,i,j,icht)
!!$                            !==========
!!$                            !fluid extrapolation
!!$                            phiw3(2,k,i,j,is)=aenscz(2,k,j)*ta3(i,j,kz  )+benscz(2,k,j)*ta3(i,j,kz+1)+&
!!$                                              censcz(2,k,j)*ta3(i,j,kz+2)+denscz(2,k,j)*ta3(i,j,kz+3)
!!$                            !solid extrapolation
!!$                            phiws3(2,k,i,j,icht)=aenscz_s(2,k,j)*tb3(i,j,kz-1)+benscz_s(2,k,j)*tb3(i,j,kz-2)+&
!!$                                               censcz_s(2,k,j)*tb3(i,j,kz-3)+denscz_s(2,k,j)*tb3(i,j,kz-4)
!!$                            !fluid: half-sum
!!$                            phiw3(2,k,i,j,is)=half*(phiw3(2,k,i,j,is)+phiws3(2,k,i,j,icht))
!!$                            !!====DEBUG
!!$                            !do irank=0,p_row*p_col-1
!!$                            !  if(nrank.eq.irank) then
!!$                            !    if (zstart(1)+i-1.eq.nx/2+1) print*,'phiw3:',  real((zstart(2)+j-1-1)*dy-yc,4),real(zf(k,i,j),4),phiw3(2,k,i,j,is),is           
!!$                            !    if (zstart(1)+i-1.eq.nx/2+1) print*,'phiws3:', real((zstart(2)+j-1-1)*dy-yc,4),real(zf(k,i,j),4),phiws3(2,k,i,j,icht),is           
!!$                            !  endif
!!$                            !enddo
!!$                            !!=========
!!$                            if (nobjz(i,j).eq.2) then !average wall temperature 
!!$                                phiwm=phiwm+phiw3(2,k,i,j,is)
!!$                                countz=countz+1.
!!$                            endif
!!$                        elseif (iflag.eq.2) then !Solid wall temperature n+1
!!$                            if (nobjz(i,j).eq.1) then !only one object
!!$                                !isoflux (IF) qn=1/g2*ra/rao
!!$                                yw=real(zstart(2)-1+j-1,mytype)*dy-yc
!!$                                zw=zf(k,i,j)-zc
!!$                                theta=atan2(zw,yw)
!!$                                !wfz=-dsin(theta)*one
!!$                                wfz=-dsin(theta)*(ra/rao)
!!$                                wfz=wfz/g2(is)
!!$                                !backward scheme
!!$                                phiws3(2,k,i,j,icht)=(-bnscz_s(2,k,j)*tb3(i,j,kz-2)-cnscz_s(2,k,j)*tb3(i,j,kz-3)&
!!$                                                      -dnscz_s(2,k,j)*tb3(i,j,kz-4)-wfz)/anscz_s(2,k,j)
!!$                            else
!!$                                !extrapolate new fluid wall value
!!$                                phiw3(2,k,i,j,is)=aenscz(2,k,j)*ta3(i,j,kz  )+benscz(2,k,j)*ta3(i,j,kz+1)+&
!!$                                                  censcz(2,k,j)*ta3(i,j,kz+2)+denscz(2,k,j)*ta3(i,j,kz+3)
!!$                                !calculate wall derivative from fluid (forward scheme)
!!$                                wfz=anscz(2,k,j)*phiw3(2,k,i,j,is)+bnscz(2,k,j)*ta3(i,j,kz+1)+& 
!!$                                    cnscz(2,k,j)*ta3(i,j,kz+2)    +dnscz(2,k,j)*ta3(i,j,kz+3) 
!!$                                wfz=wfz/g2(is)
!!$                                !add tangential flux contribution
!!$                                yw=real(zstart(2)-1+j-1,mytype)*dy-yc
!!$                                zw=zf(k,i,j)-zc
!!$                                theta=atan2(zw,yw)
!!$                                wfz=wfz+dcos(theta)*((g2(is)-one)/g2(is))*qthetaws3(2,k,i,j,icht)
!!$                                !!====DEBUG
!!$                                !if (zstart(1)+i-1.eq.nx/2+1) print*,'wfz', theta, wfz, is
!!$                                !!=========
!!$                                !impose wall derivative on solid (backward scheme)
!!$                                phiws3(2,k,i,j,icht)=(-bnscz_s(2,k,j)*tb3(i,j,kz-2)-cnscz_s(2,k,j)*tb3(i,j,kz-3)&
!!$                                                      -dnscz_s(2,k,j)*tb3(i,j,kz-4)-wfz)/anscz_s(2,k,j)
!!$                                !!====DEBUG
!!$                                !if (zstart(1)+i-1.eq.nx/2+1) print*,'phiws3', theta, phiws3(2,k,i,j,icht), is,itime
!!$                                !!=========
!!$                            endif
!!$                        endif
!!$
!!$                    elseif (k.eq.2) then !2nd object
!!$                        !1ST IMMERSED BOUNDARY
!!$                        kz=zi(k,i,j)/dz+1
!!$                        !====DEBUG
!!$                        !if (abs(zp(kz)-zi(k,i,j)).lt.tol) kz=kz-1 !use if exact_ib_pipe
!!$                        if (iflag.eq.1) then !Fluid wall temperature n+1
!!$                            !=====DEBUG
!!$                            !!IF EXACT IB PIPE: FORCE NODES AT IB
!!$                            !if (abs(zp(kz+1)-zi(k,i,j)).lt.tol.and.itime.ne.ifirst)& 
!!$                            !    tb3(i,j,kz+1)=phiws3(1,k,i,j,icht)
!!$                            !==========
!!$                            !fluid extrapolation
!!$                            phiw3(1,k,i,j,is)=aenscz(1,k,j)*ta3(i,j,kz  )+benscz(1,k,j)*ta3(i,j,kz-1)+&
!!$                                              censcz(1,k,j)*ta3(i,j,kz-2)+denscz(1,k,j)*ta3(i,j,kz-3)
!!$                            !solid extrapolation
!!$                            phiws3(1,k,i,j,icht)=aenscz_s(1,k,j)*tb3(i,j,kz+1)+benscz_s(1,k,j)*tb3(i,j,kz+2)+&
!!$                                               censcz_s(1,k,j)*tb3(i,j,kz+3)+denscz_s(1,k,j)*tb3(i,j,kz+4)
!!$                            !fluid: half-sum
!!$                            phiw3(1,k,i,j,is)=half*(phiw3(1,k,i,j,is)+phiws3(1,k,i,j,icht))
!!$                            !!====DEBUG
!!$                            !do irank=0,p_row*p_col-1
!!$                            !  if(nrank.eq.irank) then
!!$                            !    if (zstart(1)+i-1.eq.nx/2+1) print*,'phiw3:',  real((zstart(2)+j-1-1)*dy-yc,4),real(zi(k,i,j),4),phiw3(1,k,i,j,is),is
!!$                            !    if (zstart(1)+i-1.eq.nx/2+1) print*,'phiws3:', real((zstart(2)+j-1-1)*dy-yc,4),real(zi(k,i,j),4),phiws3(1,k,i,j,icht),is
!!$                            !  endif
!!$                            !enddo
!!$                            !!=========
!!$                            if (nobjz(i,j).eq.2) then !average wall temperature 
!!$                                phiwm=phiwm+phiw3(1,k,i,j,is)
!!$                                countz=countz+1.
!!$                            endif
!!$                        elseif (iflag.eq.2) then  !Solid wall temperature n+1
!!$                            !extrapolate new fluid wall value
!!$                            phiw3(1,k,i,j,is)=aenscz(1,k,j)*ta3(i,j,kz)+benscz(1,k,j)*ta3(i,j,kz-1)+&
!!$                                              censcz(1,k,j)*ta3(i,j,kz-2)+denscz(1,k,j)*ta3(i,j,kz-3)
!!$                            !calculate wall derivative from fluid (backward scheme)
!!$                            wfz=anscz(1,k,j)*phiw3(1,k,i,j,is)+bnscz(1,k,j)*ta3(i,j,kz-1)+&
!!$                                cnscz(1,k,j)*ta3(i,j,kz-2)+dnscz(1,k,j)*ta3(i,j,kz-3)
!!$                            wfz=-wfz/g2(is)
!!$                            !add tangential flux contribution
!!$                            yw=real(zstart(2)-1+j-1,mytype)*dy-yc
!!$                            zw=zi(k,i,j)-zc
!!$                            theta=atan2(zw,yw)
!!$                            wfz=wfz+dcos(theta)*((g2(is)-one)/g2(is))*qthetaws3(1,k,i,j,icht)
!!$                            !!====DEBUG
!!$                            !if (zstart(1)+i-1.eq.nx/2+1) print*,'wfz', theta, wfz, is
!!$                            !!=========
!!$                            !impose wall derivative on solid (forward scheme)
!!$                            phiws3(1,k,i,j,icht)=(-bnscz_s(1,k,j)*tb3(i,j,kz+2)-cnscz_s(1,k,j)*tb3(i,j,kz+3)&
!!$                                                  -dnscz_s(1,k,j)*tb3(i,j,kz+4)+wfz)/anscz_s(1,k,j)
!!$                            !!====DEBUG
!!$                            !if (zstart(1)+i-1.eq.nx/2+1) print*,'phiws3', theta, phiws3(1,k,i,j,icht), is, itime
!!$                            !!=========
!!$                        endif
!!$
!!$                        !2ND IMMERSED BOUNDARY
!!$                        kz=(zf(k,i,j)+dz)/dz+1 !without skip
!!$                        if (iflag.eq.1) then !Fluid wall temperature n+1
!!$                            !=====DEBUG
!!$                            !!IF EXACT IB PIPE: FORCE NODES AT IB
!!$                            !if (abs(zp(kz-1)-zf(k,i,j)).lt.tol.and.itime.ne.ifirst)& 
!!$                            !    tb3(i,j,kz-1)=phiws3(2,k,i,j,icht)
!!$                            !==========
!!$                            !solid extrapolation
!!$                            phiws3(2,k,i,j,icht)=aenscz_s(2,k,j)*tb3(i,j,kz-1)+benscz_s(2,k,j)*tb3(i,j,kz-2)+&
!!$                                               censcz_s(2,k,j)*tb3(i,j,kz-3)+denscz_s(2,k,j)*tb3(i,j,kz-4)
!!$                            !fluid (equal)
!!$                            phiw3(2,k,i,j,is)=phiws3(2,k,i,j,icht)
!!$                        elseif (iflag.eq.2) then !Solid wall temperature n+1
!!$                            !isoflux (IF) qn=1/g2*ra/rao
!!$                            yw=real(zstart(2)-1+j-1,mytype)*dy-yc
!!$                            zw=zf(k,i,j)-zc
!!$                            theta=atan2(zw,yw)
!!$                            !wfz=-dsin(theta)*one
!!$                            wfz=-dsin(theta)*(ra/rao)
!!$                            wfz=wfz/g2(is)
!!$                            !backward scheme
!!$                            phiws3(2,k,i,j,icht)=(-bnscz_s(2,k,j)*tb3(i,j,kz-2)-cnscz_s(2,k,j)*tb3(i,j,kz-3)&
!!$                                                  -dnscz_s(2,k,j)*tb3(i,j,kz-4)-wfz)/anscz_s(2,k,j)
!!$                        endif
!!$                    endif
!!$                enddo
!!$            endif
!!$        enddo
!!$    enddo
!!$
!!$    !!====DEBUG
!!$    !if (iflag.eq.1) then ! use if exact ib pipe
!!$    !    call transpose_z_to_y(tb3,tb2)
!!$    !    call transpose_y_to_x(tb2,phis1)
!!$    !endif
!!$
!!$    if (iflag.eq.1) then
!!$        call MPI_ALLREDUCE(countz,countzt,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
!!$        call MPI_ALLREDUCE(phiwm,phiwmt,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
!!$        phiwm=phiwmt/(countyt+countzt)
!!$        nuw(is)=-one/phiwm !Nu(t)=-1/<phi>
!!$        !if (nrank==0) then
!!$        if (nrank.eq.0.and.iwrite.ne.0) then
!!$           !print *,'Nu(t)       =       ',nuw(is)
!!$           write(ifile,*) real((itime-1)*dt,mytype), nuw(is)
!!$        endif
!!$        !For reconstruction smoothness
!!$        do k=1,xsize(3)
!!$            zm=dz*real(xstart(3)-1+k-1,mytype)-zc 
!!$            do j=1,xsize(2)
!!$                if (istret.eq.0) ym=real(j+xstart(2)-1-1,mytype)*dy-yc
!!$                if (istret.ne.0) ym=yp(j+xstart(2)-1)-yc
!!$                r=sqrt(ym*ym+zm*zm)    
!!$                do i=1,xsize(1)
!!$                    if (r.ge.ra.or.ep1(i,j,k).eq.1) then
!!$                        phi1(i,j,k)=-one/nuw(is)
!!$                    endif
!!$                enddo
!!$            enddo
!!$        enddo
!!$    endif
!!$    
!!$    return
!!$  
!!$  end subroutine phiw_cht
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !!
!!$  !!  subroutine: derphiw_cht
!!$  !!      AUTHOR: Rodrigo Vicente Cruz
!!$  !! DESCRIPTION: Compute d(phis)/dtheta^(n+1) (tangential flux)
!!$  !!              according to coupling strategy (icoupling)
!!$  !!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !***************************************************************************
!!$  !
!!$  subroutine derphiw_cht(dyphi3,dyphis3,dzphi3,dzphis3,is,icht,icoupling)
!!$  !
!!$  !***************************************************************************
!!$  use decomp_2d
!!$  use decomp_2d_io
!!$  use MPI
!!$  use var, only   : tg2,th2,tg3,th3
!!$  use param, only : zero,one,two,three,four
!!$  use param, only : yly,zlz,dy,dz,dt,itime
!!$  use complex_geometry, only: nobjy,nobjz,yi,yf,zi,zf
!!$  use ibm
!!$  use variables, only : ny,nz
!!$  !
!!$  implicit none
!!$  real(mytype),dimension(zsize(1),zsize(2),zsize(3))    :: dyphi3,dzphi3
!!$  real(mytype),dimension(zsize(1),zsize(2),zsize(3))    :: dyphis3,dzphis3
!!$  integer                                               :: is,icht,icoupling
!!$  !LOCALS
!!$  integer                                               :: i,j,k
!!$  integer                                               :: jy,kz,irank
!!$  real(mytype)                                          :: r,yc,zc
!!$  real(mytype)                                          :: ym,zm,theta  !mesh
!!$  !------------------------------------------------------------------------------!
!!$  ! COUPLING STRATEGIES (icoupling) :                                            !
!!$  !     1 : d(phis)/dtheta^(n+1) = d(phis)/dtheta^(n)                            !
!!$  !     2 : d(phis)/dtheta^(n+1) = d(phi )/dtheta^(n+1)                          !
!!$  !     3 : d(phis)/dtheta^(n+1) = 0.5*(d(phi )/dtheta^(n+1)+d(phis)/dtheta^(n)) !
!!$  !------------------------------------------------------------------------------!
!!$  yc=yly/two
!!$  zc=zlz/two
!!$  !================================= Z-PENCIL ==============================
!!$  th3(:,:,:)=zero
!!$  do i=1,zsize(1)
!!$      do j=1,zsize(2)
!!$          ym=real(zstart(2)-1+j-1,mytype)*dy-yc
!!$
!!$          !Tangential derivative computation
!!$          do k=1,zsize(3)
!!$              zm=real(zstart(3)-1+k-1,mytype)*dz-zc
!!$              r=sqrt(ym*ym+zm*zm)
!!$              theta=atan2(zm,ym)
!!$
!!$              !Fluid: (dphi/dtheta)^(n+1)
!!$              tg3(i,j,k)=-dsin(theta)*dyphi3(i,j,k)+dcos(theta)*dzphi3(i,j,k)
!!$              !Solid: (dphis/dtheta)^(n)
!!$              th3(i,j,k)=-dsin(theta)*dyphis3(i,j,k)+dcos(theta)*dzphis3(i,j,k)
!!$          enddo
!!$        
!!$          !Extrapolation of wall values (z-direction)
!!$          if (nobjz(i,j).eq.2) then !immersed object
!!$              do k=1,nobjz(i,j)
!!$                  if (k.eq.1) then !1st object
!!$                      kz=(zf(k,i,j)+dz)/dz+1 !without skip
!!$                      !wall-tangential derivative (fluid)
!!$                      qthetaw3(2,k,i,j,is)=aenscz(2,k,j)*tg3(i,j,kz  )+benscz(2,k,j)*tg3(i,j,kz+1)+&
!!$                                           censcz(2,k,j)*tg3(i,j,kz+2)+denscz(2,k,j)*tg3(i,j,kz+3)
!!$                      !wall-tangential derivative (solid)
!!$                      qthetaws3(2,k,i,j,icht)=aenscz_s(2,k,j)*th3(i,j,kz-1)+benscz_s(2,k,j)*th3(i,j,kz-2)+&
!!$                                              censcz_s(2,k,j)*th3(i,j,kz-3)+denscz_s(2,k,j)*th3(i,j,kz-4)
!!$                      !coupling strategy
!!$                      if (icoupling.eq.1) then
!!$                          qthetaws3(2,k,i,j,icht)=qthetaws3(2,k,i,j,icht)
!!$                      elseif (icoupling.eq.2) then
!!$                          qthetaws3(2,k,i,j,icht)=qthetaw3(2,k,i,j,is)
!!$                      elseif (icoupling.eq.3) then
!!$                          qthetaws3(2,k,i,j,icht)=half*(qthetaws3(2,k,i,j,icht)+qthetaw3(2,k,i,j,is))
!!$                      endif
!!$                  elseif (k.eq.2) then !2nd object
!!$                      kz=zi(k,i,j)/dz+1
!!$                      !====DEBUG
!!$                      !if (abs(zp(kz)-zi(k,i,j)).lt.tol) kz=kz-1 !use if exact_ib_pipe
!!$                      !wall-tangential derivative (fluid)
!!$                      qthetaw3(1,k,i,j,is)=aenscz(1,k,j)*tg3(i,j,kz  )+benscz(1,k,j)*tg3(i,j,kz-1)+&
!!$                                           censcz(1,k,j)*tg3(i,j,kz-2)+denscz(1,k,j)*tg3(i,j,kz-3)
!!$                      !wall-tangential derivative (solid)
!!$                      qthetaws3(1,k,i,j,icht)=aenscz_s(1,k,j)*th3(i,j,kz+1)+benscz_s(1,k,j)*th3(i,j,kz+2)+&
!!$                                              censcz_s(1,k,j)*th3(i,j,kz+3)+denscz_s(1,k,j)*th3(i,j,kz+4)
!!$                      !coupling strategy
!!$                      if (icoupling.eq.1) then
!!$                          qthetaws3(1,k,i,j,icht)=qthetaws3(1,k,i,j,icht)
!!$                      elseif (icoupling.eq.2) then
!!$                          qthetaws3(1,k,i,j,icht)=qthetaw3(1,k,i,j,is)
!!$                      elseif (icoupling.eq.3) then
!!$                          qthetaws3(1,k,i,j,icht)=half*(qthetaws3(1,k,i,j,icht)+qthetaw3(1,k,i,j,is))
!!$                      endif
!!$                  endif
!!$                  !!====DEBUG
!!$                  !do irank=0,p_row*p_col-1
!!$                  !  if(nrank.eq.irank) then
!!$                  !    if (k.eq.1.and.zstart(1)+i-1.eq.nx/2+1) print*,'qtw3 :', atan2(zf(1,i,j)-zc, ym), qthetaw3 (2,k,i,j,is),is
!!$                  !    if (k.eq.2.and.zstart(1)+i-1.eq.nx/2+1) print*,'qtw3 :', atan2(zi(2,i,j)-zc, ym), qthetaw3 (1,k,i,j,is),is
!!$                  !    if (k.eq.1.and.zstart(1)+i-1.eq.nx/2+1) print*,'qtws3:', atan2(zf(1,i,j)-zc, ym), qthetaws3(2,k,i,j,icht),is
!!$                  !    if (k.eq.2.and.zstart(1)+i-1.eq.nx/2+1) print*,'qtws3:', atan2(zi(2,i,j)-zc, ym), qthetaws3(1,k,i,j,icht),is
!!$                  !  endif
!!$                  !enddo
!!$                  !!=========
!!$              enddo
!!$          endif
!!$      enddo
!!$  enddo
!!$  call transpose_z_to_y(tg3,tg2) !(dphi/dtheta)^(n+1)
!!$  call transpose_z_to_y(th3,th2) !(dphis/dtheta)^(n)
!!$
!!$  !================================= Y-PENCIL ==============================
!!$  do i=1,ysize(1)
!!$      do k=1,ysize(3)
!!$          zm=real(ystart(3)-1+k-1,mytype)*dz-zc
!!$
!!$          !Extrapolation of wall values (y-direction)
!!$          if (nobjy(i,k).eq.2) then !immersed object
!!$              do j=1,nobjy(i,k)
!!$                  if (j.eq.1) then !1st object
!!$                      jy=1!jy=(yf(j,i,k)+dy)/dy+1
!!$                      !====DEBUG
!!$                      !do while(yp(jy).le.yf(j,i,k)) !use if exact_ib_pipe
!!$                      do while(yp(jy).lt.yf(j,i,k))
!!$                         jy=jy+1
!!$                      enddo
!!$                      !wall-tangential derivative (fluid)
!!$                      qthetaw2(2,j,i,k,is)=aenscy(2,j,k)*tg2(i,jy  ,k)+benscy(2,j,k)*tg2(i,jy+1,k)+&
!!$                                           censcy(2,j,k)*tg2(i,jy+2,k)+denscy(2,j,k)*tg2(i,jy+3,k)
!!$                      !wall-tangential derivative (solid)
!!$                      qthetaws2(2,j,i,k,icht)=aenscy_s(2,j,k)*th2(i,jy-1,k)+benscy_s(2,j,k)*th2(i,jy-2,k)+&
!!$                                              censcy_s(2,j,k)*th2(i,jy-3,k)+denscy_s(2,j,k)*th2(i,jy-4,k)
!!$                      !coupling strategy
!!$                      if (icoupling.eq.1) then
!!$                          qthetaws2(2,j,i,k,icht)=qthetaws2(2,j,i,k,icht)
!!$                      elseif (icoupling.eq.2) then
!!$                          qthetaws2(2,j,i,k,icht)=qthetaw2(2,j,i,k,is)
!!$                      elseif (icoupling.eq.3) then
!!$                          qthetaws2(2,j,i,k,icht)=half*(qthetaws2(2,j,i,k,icht)+qthetaw2(2,j,i,k,is))
!!$                      endif
!!$                  elseif (j.eq.2) then !2nd object
!!$                      jy=1!jy=yi(j,i,k)/dy+1
!!$                      do while(yp(jy).lt.yi(j,i,k))
!!$                         jy=jy+1
!!$                      enddo
!!$                      jy=jy-1
!!$                      !wall-tangential derivative (fluid)
!!$                      qthetaw2(1,j,i,k,is)=aenscy(1,j,k)*tg2(i,jy  ,k)+benscy(1,j,k)*tg2(i,jy-1,k)+&
!!$                                           censcy(1,j,k)*tg2(i,jy-2,k)+denscy(1,j,k)*tg2(i,jy-3,k)
!!$                      !wall-tangential derivative (solid)
!!$                      qthetaws2(1,j,i,k,icht)=aenscy_s(1,j,k)*th2(i,jy+1,k)+benscy_s(1,j,k)*th2(i,jy+2,k)+&
!!$                                              censcy_s(1,j,k)*th2(i,jy+3,k)+denscy_s(1,j,k)*th2(i,jy+4,k)
!!$                      !coupling strategy
!!$                      if (icoupling.eq.1) then
!!$                          qthetaws2(1,j,i,k,icht)=qthetaws2(1,j,i,k,icht)
!!$                      elseif (icoupling.eq.2) then
!!$                          qthetaws2(1,j,i,k,icht)=qthetaw2(1,j,i,k,is)
!!$                      elseif (icoupling.eq.3) then
!!$                          qthetaws2(1,j,i,k,icht)=half*(qthetaws2(1,j,i,k,icht)+qthetaw2(1,j,i,k,is))
!!$                      endif
!!$                  endif
!!$                  !!====DEBUG
!!$                  !do irank=0,p_row*p_col-1
!!$                  !  if(nrank.eq.irank) then
!!$                  !    if (j.eq.1.and.ystart(1)+i-1.eq.nx/2+1) print*,'qtw2 :', atan2(zm, yf(1,i,k)-yc), qthetaw2 (2,j,i,k,is),is
!!$                  !    if (j.eq.2.and.ystart(1)+i-1.eq.nx/2+1) print*,'qtw2 :', atan2(zm, yi(2,i,k)-yc), qthetaw2 (1,j,i,k,is),is
!!$                  !    if (j.eq.1.and.ystart(1)+i-1.eq.nx/2+1) print*,'qtws2:', atan2(zm, yf(1,i,k)-yc), qthetaws2(2,j,i,k,icht),is
!!$                  !    if (j.eq.2.and.ystart(1)+i-1.eq.nx/2+1) print*,'qtws2:', atan2(zm, yi(2,i,k)-yc), qthetaws2(1,j,i,k,icht),is
!!$                  !  endif
!!$                  !enddo
!!$                  !!=========
!!$              enddo
!!$          endif
!!$      enddo
!!$  enddo
!!$  !
!!$  return
!!$  !
!!$  end subroutine derphiw_cht
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    !!
!!$    !!  subroutine: phis_condeq
!!$    !!      AUTHOR: Rodrigo Vicente Cruz
!!$    !! DESCRIPTION: Solver for heat diffusion equation in the solid domain
!!$    !!
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !***************************************************************************
!!$  !
!!$  subroutine phis_condeq(phis1,dphis1,is,icht)
!!$  !
!!$  !***************************************************************************
!!$  use param
!!$  use variables
!!$  use decomp_2d
!!$  use MPI
!!$  use ibm
!!$  use time_integrators, only: intt
!!$  use var, ONLY : ep1
!!$
!!$  use var, ONLY : ta1,tb1,di1
!!$  use var, ONLY : ta2,tb2,tc2,di2
!!$  use var, ONLY : ta3,tb3,tc3,di3
!!$
!!$  implicit none
!!$  real(mytype),dimension(xsize(1),xsize(2),xsize(3))       :: phis1
!!$  real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dphis1
!!$  integer                                                  :: is,icht
!!$  !LOCALS
!!$  real(mytype)                                             :: xalpha
!!$  real(mytype)                                             :: zm,ym,zc,yc,r
!!$  integer                                                  :: i,j,k,iibm_save
!!$
!!$    xalpha = xnu/(sc(is)*g1(is))
!!$
!!$    if (ivf.eq.0) then
!!$        !X PENCILS
!!$        iibm_save=iibm
!!$        iibm=0
!!$        call derxxS(ta1,phis1(:,:,:),di1,sx,sfxpS,ssxpS,swxpS,xsize(1),xsize(2),xsize(3),1)
!!$        call transpose_x_to_y(phis1,tb2)
!!$
!!$        !Y PENCILS
!!$        call chtlagpoly_s(tb2,icht)
!!$        call deryyS(ta2,tb2(:,:,:),di2,sy,sfypS,ssypS,swypS,ysize(1),ysize(2),ysize(3),1)        
!!$        !!====DEBUG Y
!!$        !do i=1,ysize(1)
!!$        !    do k=1,ysize(3)
!!$        !        if (ystart(1)+i-1.eq.nx/2+1.and.ystart(3)+k-1.eq.nz/2+1) then
!!$        !            do j=1,ysize(2)
!!$        !                print*,'phis_y:' ,yp(j),tb2(i,j,k),is
!!$        !                print*,'dyy:' ,yp(j),ta2(i,j,k),is
!!$        !            enddo
!!$        !        endif
!!$        !    enddo
!!$        !enddo
!!$        !!====
!!$        call transpose_y_to_z(tb2(:,:,:),tb3(:,:,:))
!!$
!!$        !Z PENCILS
!!$        call chtlagpolz_s(tb3,icht)
!!$        call derzzS(ta3,tb3(:,:,:),di3,sz,sfzpS,sszpS,swzpS,zsize(1),zsize(2),zsize(3),1)
!!$        !!====DEBUG Z
!!$        !do i=1,zsize(1)
!!$        !    do j=1,zsize(2)
!!$        !        if (zstart(1)+i-1.eq.nx/2+1.and.zstart(2)+j-1.eq.ny/2+1) then
!!$        !            do k=1,zsize(3)
!!$        !                print*,'phis_z:', (k-1)*dz,tb3(i,j,k),is
!!$        !                print*,'dzz:' ,  (k-1)*dz,ta3(i,j,k),is
!!$        !            enddo
!!$        !        endif
!!$        !    enddo
!!$        !enddo
!!$        !!====
!!$        iibm=iibm_save
!!$        call transpose_z_to_y(ta3(:,:,:),tc2(:,:,:))
!!$
!!$        !SUM DIFFUSIVE TERMS Y-PENCILS
!!$        ta2(:,:,:)=ta2(:,:,:)+tc2(:,:,:)
!!$        call transpose_y_to_x(ta2(:,:,:),tb1(:,:,:))
!!$
!!$        !SUM DIFFUSIVE TERMS X-PENCILS
!!$        dphis1(:,:,:,1)=xalpha*(ta1(:,:,:)+tb1(:,:,:))
!!$
!!$    else !Viscous filtering
!!$
!!$        dphis1(:,:,:,1)=zero
!!$
!!$    endif
!!$
!!$    !=============== TIME ADVANCEMENT ===============
!!$    CALL intt(phis1(:,:,:), dphis1(:,:,:,:))
!!$    
!!$    if (ivf.ne.0) then !Viscous filtering
!!$        call viscous_filter(phis1(:,:,:),-is,1)
!!$        !call viscous_filter(phis1(:,:,:),-icht,1)
!!$    endif
!!$
!!$    !====DEBUG
!!$    yc=yly/two
!!$    zc=zlz/two
!!$    do k=1,xsize(3)
!!$        zm=dz*real(xstart(3)-1+k-1,mytype)-zc 
!!$        do j=1,xsize(2)
!!$            if (istret.eq.0) ym=real(j+xstart(2)-1-1,mytype)*dy-yc
!!$            if (istret.ne.0) ym=yp(j+xstart(2)-1)-yc
!!$            r=sqrt(ym*ym+zm*zm)    
!!$            do i=1,xsize(1)
!!$                if (ep1(i,j,k).eq.0) then
!!$                    phis1(i,j,k)=zero
!!$                endif
!!$            enddo
!!$        enddo
!!$    enddo
!!$
!!$  end subroutine phis_condeq
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    !!
!!$    !!  subroutine: lagpol*2
!!$    !!      AUTHOR: Rodrigo Vicente Cruz
!!$    !! DESCRIPTION: Reconstruction across transverse-yz domain periodicity
!!$    !!
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !***************************************************************************
!!$  !
!!$  subroutine lagpoly2(u)
!!$  !
!!$  !***************************************************************************
!!$    !
!!$    use param
!!$    use complex_geometry
!!$    use decomp_2d
!!$    use variables
!!$    !
!!$    implicit none
!!$    !
!!$    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: u
!!$    integer                                            :: i,j,k
!!$    real(mytype)                                       :: x,y,z
!!$    integer                                            :: jy              != position du point "zappé"
!!$    integer                                            :: jpif,jpol,nypif
!!$    integer                                            :: jpoli,jpolf     != positions Initiales et Finales du POLynôme considéré
!!$    integer                                            :: nyd,jm
!!$    real(mytype)                                       :: xpol,ypol,dypol !|variables concernant les polynômes
!!$    real(mytype)                                       :: yd
!!$    real(mytype),dimension(10)                         :: xa,ya           !|de Lagrange. A mettre impérativement en 
!!$    integer                                            :: ia,na           !|double précision
!!$    integer                                            :: jj,k_visu
!!$    !
!!$    do k=1,ysize(3)
!!$        do i=1,ysize(1)
!!$            !if(nobjy(i,k).ne.0)then
!!$            if(nobjy(i,k).eq.2)then
!!$                ia=0
!!$                do j=1,nobjy(i,k)          !boucle sur le nombre d'objets par (j,k)
!!$                    !1ère frontière
!!$                    if (j.eq.2) then
!!$                        nypif=npif
!!$                        ia=ia+1
!!$                        xa(ia)=yi(j,i,k)
!!$                        ya(ia)=0.
!!$                        yd=yly
!!$                        nyd=ny
!!$                        xa(ia)=xa(ia)-yd
!!$                        if(yi(j,i,k).gt.0.)then!objet immergé
!!$                            jy=1!jy=yi(j,i,k)/dy+1
!!$                            do while(yp(jy).lt.yi(j,i,k))
!!$                                jy=jy+1
!!$                            enddo
!!$                            jy=jy-1
!!$                            jpoli=jy+1
!!$                            jpoli=jpoli-nyd
!!$                            if(nyipif(j,i,k).lt.npif)nypif=nyipif(j,i,k)
!!$                            do jpif=1,nypif
!!$                                ia=ia+1
!!$                                if(izap.eq.1)then!zapping
!!$                                   !xa(ia)=yp(jy-jpif)!(jy-1)*dy-jpif*dy
!!$                                   xa(ia)=(jy-1)*dy-jpif*dy
!!$                                   ya(ia)=u(i,jy-jpif,k)
!!$                                else             !no zapping
!!$                                   !xa(ia)=yp(jy-jpif+1)!(jy-1)*dy-(jpif-1)*dy
!!$                                   xa(ia)=(jy-1)*dy-(jpif-1)*dy
!!$                                   ya(ia)=u(i,jy-jpif+1,k)
!!$                                endif
!!$                                xa(ia)=xa(ia)-yd
!!$                            enddo
!!$                        else                   !objet semi-immergé
!!$                           jpoli=1
!!$                        endif
!!$                    endif
!!$                    !2ème frontière
!!$                    if (j.eq.1) then
!!$                        nypif=npif
!!$                        ia=ia+1
!!$                        xa(ia)=yf(j,i,k)
!!$                        ya(ia)=0.
!!$                        if(yf(j,i,k).lt.yly)then!objet immergé
!!$                           jy=1!jy=(yf(j,i,k)+dy)/dy+1
!!$                           !====DEBUG
!!$                           !do while(yp(jy).le.yf(j,i,k)) !use if exact_ib_pipe
!!$                           do while(yp(jy).lt.yf(j,i,k))  !there was a bug here yi<-->yf
!!$                              jy=jy+1
!!$                           enddo
!!$                           jpolf=jy-1
!!$                           if(nyfpif(j,i,k).lt.npif)nypif=nyfpif(j,i,k)
!!$                           do jpif=1,nypif
!!$                              ia=ia+1
!!$                              if(izap.eq.1)then!zapping
!!$                                 !xa(ia)=yp(jy+jpif)!(jy-1)*dy+jpif*dy
!!$                                 xa(ia)=(jy-1)*dy+jpif*dy
!!$                                 ya(ia)=u(i,jy+jpif,k)
!!$                              else             !no zapping
!!$                                 !xa(ia)=yp(jy+jpif-1)!(jy-1)*dy+(jpif-1)*dy
!!$                                 xa(ia)=(jy-1)*dy+(jpif-1)*dy
!!$                                 ya(ia)=u(i,jy+jpif-1,k)
!!$                              endif
!!$                           enddo
!!$                        else                   !objet semi-immergé
!!$                           jpolf=ny
!!$                        endif
!!$                    endif
!!$                enddo
!!$                !calcul du polynôme
!!$                na=ia
!!$                do jpol=jpoli,jpolf
!!$                    !xpol=yp(jpol)!dy*(jpol-1)
!!$                    xpol=(jpol-1)*dy
!!$                    call polint(xa,ya,na,xpol,ypol,dypol)
!!$                    !call csplint(xa,ya,na,xpol,ypol)
!!$                    if (jpoli.lt.0.and.jpol.le.0) then
!!$                        jm=jpol+nyd
!!$                        u(i,jm,k)=ypol
!!$                    else
!!$                        u(i,jpol,k)=ypol
!!$                    endif
!!$                enddo
!!$                !!====DEBUG Y
!!$                !if (itime.eq.1) then
!!$                !    if (ystart(1)+i-1.eq.nx/2+1.and.ystart(3)+k-1.eq.nz/2+1) then
!!$                !        do jj=1,ysize(2)
!!$                !            print*,'uy:' ,yp(jj),u(i,jj,k)
!!$                !        enddo
!!$                !        do ia=1,na
!!$                !            print*,'uy:xa,ya:',xa(ia),ya(ia)
!!$                !        enddo
!!$                !    endif
!!$                !    !stop
!!$                !endif
!!$                !!====
!!$                ia=0
!!$            else
!!$                do j=1,ysize(2)
!!$                    u(i,j,k)=zero
!!$                enddo
!!$            endif
!!$        enddo
!!$    enddo
!!$    !
!!$    return
!!$  end subroutine lagpoly2
!!$  !***************************************************************************
!!$  !
!!$  subroutine lagpolz2(u)
!!$  !
!!$  !***************************************************************************
!!$    !
!!$    use param
!!$    use complex_geometry
!!$    use decomp_2d
!!$    use variables
!!$    !
!!$    implicit none
!!$    !
!!$    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: u
!!$    integer                                            :: i,j,k
!!$    real(mytype)                                       :: x,y,z
!!$    integer                                            :: kz              != position du point "zappé"
!!$    integer                                            :: kpif,kpol,nzpif
!!$    integer                                            :: kpoli,kpolf     != positions Initiales et Finales du POLynôme considéré
!!$    integer                                            :: nzd,km
!!$    real(mytype)                                       :: xpol,ypol,dypol !|variables concernant les polynômes
!!$    real(mytype)                                       :: zd
!!$    real(mytype),dimension(10)                         :: xa,ya           !|de Lagrange. A mettre imérativement en 
!!$    integer                                            :: ia,na           !|double précision
!!$    !
!!$    do j=1,zsize(2)
!!$        do i=1,zsize(1)
!!$            !if(nobjz(i,j).ne.0)then
!!$            if(nobjz(i,j).eq.2)then
!!$                ia=0
!!$                do k=1,nobjz(i,j)          !boucle sur le nombre d'objets par couple (i,j)
!!$                    !1ère frontière
!!$                    if (k.eq.2) then
!!$                        nzpif=npif
!!$                        ia=ia+1
!!$                        xa(ia)=zi(k,i,j)
!!$                        ya(ia)=0.
!!$                        zd=zlz
!!$                        nzd=nz
!!$                        xa(ia)=xa(ia)-zd
!!$                        if(zi(k,i,j).gt.0.)then!objet immergé
!!$                            kz=zi(k,i,j)/dz+1
!!$                            !====DEBUG
!!$                            !if (abs(zp(kz)-zi(k,i,j)).lt.tol) kz=kz-1 !use if exact_ib_pipe
!!$                            kpoli=kz+1
!!$                            kpoli=kpoli-nzd
!!$                            if(nzipif(k,i,j).lt.npif)nzpif=nzipif(k,i,j)
!!$                            do kpif=1,nzpif
!!$                                ia=ia+1
!!$                                if(izap.eq.1)then!zapping
!!$                                    xa(ia)=(kz-1)*dz-kpif*dz
!!$                                    ya(ia)=u(i,j,kz-kpif)
!!$                                else             !no zapping
!!$                                    xa(ia)=(kz-1)*dz-(kpif-1)*dz
!!$                                    ya(ia)=u(i,j,kz-kpif+1)
!!$                                endif
!!$                                xa(ia)=xa(ia)-zd
!!$                            enddo
!!$                        else                   !objet semi-immergé
!!$                           kpoli=1
!!$                        endif
!!$                    endif
!!$                    !2ème frontière
!!$                    if (k.eq.1) then
!!$                        nzpif=npif
!!$                        ia=ia+1
!!$                        xa(ia)=zf(k,i,j)
!!$                        ya(ia)=0.
!!$                        if(zf(k,i,j).lt.zlz)then!objet immergé
!!$                           kz=(zf(k,i,j)+dz)/dz+1
!!$                           kpolf=kz-1
!!$                           if(nzfpif(k,i,j).lt.npif)nzpif=nzfpif(k,i,j)
!!$                           do kpif=1,nzpif
!!$                              ia=ia+1
!!$                              if(izap.eq.1)then!zapping
!!$                                 xa(ia)=(kz-1)*dz+kpif*dz
!!$                                 ya(ia)=u(i,j,kz+kpif)
!!$                              else             !no zapping
!!$                                 xa(ia)=(kz-1)*dz+(kpif-1)*dz
!!$                                 ya(ia)=u(i,j,kz+kpif-1)
!!$                              endif
!!$                           enddo
!!$                        else                   !objet semi-immergé
!!$                           kpolf=nz
!!$                        endif
!!$                    endif
!!$                enddo
!!$                !calcul du polynôme
!!$                na=ia
!!$                do kpol=kpoli,kpolf
!!$                    xpol=dz*(kpol-1)
!!$                    call polint(xa,ya,na,xpol,ypol,dypol)
!!$                    if (kpoli.lt.0.and.kpol.le.0) then
!!$                        km=kpol+nzd
!!$                        u(i,j,km)=ypol
!!$                    else
!!$                        u(i,j,kpol)=ypol
!!$                    endif
!!$                enddo
!!$                !!====DEBUG Z
!!$                !if (itime.eq.1) then
!!$                !    if (zstart(1)+i-1.eq.nx/2+1.and.zstart(2)+j-1.eq.nz/2+1) then
!!$                !        !do kk=1,zsize(3)
!!$                !        !    print*,'u_z:' ,(kk-1)*dz,u(i,j,kk)
!!$                !        !enddo
!!$                !        do ia=1,na
!!$                !            print*,'uz:xa,ya:',xa(ia),ya(ia)
!!$                !        enddo
!!$                !    endif
!!$                !    !stop
!!$                !endif
!!$                !!====
!!$                ia=0
!!$            else
!!$                do k=1,zsize(3)
!!$                    u(i,j,k)=zero
!!$                enddo
!!$            endif
!!$        enddo
!!$    enddo
!!$    !
!!$    return
!!$  end subroutine lagpolz2
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    !!
!!$    !!  subroutine: nbclagpol*
!!$    !!      AUTHOR: Rodrigo Vicente Cruz
!!$    !! DESCRIPTION: Performs the reconstruction ensuring a targeted value for 
!!$    !!              the wall derivative (Neumann BC is ensured through a
!!$    !!              targeted Dirichlet)
!!$    !!
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !***************************************************************************
!!$  !
!!$  subroutine nbclagpoly(u,is)
!!$  !
!!$  !***************************************************************************
!!$    !
!!$    use param
!!$    use complex_geometry
!!$    use decomp_2d
!!$    use variables
!!$    !
!!$    implicit none
!!$    !
!!$    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: u
!!$    integer                                            :: is
!!$    integer                                            :: i,j,k
!!$    real(mytype)                                       :: x,y,z
!!$    real(mytype)                                       :: yw,zw,yc,zc,theta
!!$    real(mytype)                                       :: wfy             ! wall flux in Y  
!!$    integer                                            :: jy              != position du point "zappé"
!!$    integer                                            :: jys
!!$    integer                                            :: jpif,jpol,nypif
!!$    integer                                            :: jpoli,jpolf     != positions Initiales et Finales du POLynôme considéré
!!$    real(mytype)                                       :: xpol,ypol,dypol !|variables concernant les polynômes
!!$    real(mytype),dimension(10)                         :: xa,ya           !|de Lagrange. A mettre impérativement en 
!!$    integer                                            :: ia,na           !|double précision
!!$    !
!!$    if (new_rec.eq.1) then
!!$        call nbclagpoly2(u,is)
!!$        return
!!$    endif
!!$    !
!!$    yc=yly/two
!!$    zc=zlz/two
!!$    do k=1,ysize(3)
!!$       do i=1,ysize(1)
!!$          if(nobjy(i,k).ne.0)then
!!$             ia=0
!!$             do j=1,nobjy(i,k)          !boucle sur le nombre d'objets par (j,k)
!!$                !1ère frontière
!!$                nypif=npif
!!$                ia=ia+1
!!$                xa(ia)=yi(j,i,k)
!!$                if (nobjy(i,k).eq.2.and.j.eq.2) then
!!$                    ya(ia)=phiw2(1,j,i,k,is)
!!$                else
!!$                    ya(ia)=-one/nuw(is)
!!$                endif
!!$                if(yi(j,i,k).gt.0.)then!objet immergé
!!$                   jy=1!jy=yi(j,i,k)/dy+1
!!$                   do while(yp(jy).lt.yi(j,i,k))
!!$                      jy=jy+1
!!$                   enddo
!!$                   jy=jy-1
!!$                   jpoli=jy+1
!!$                   if(nyipif(j,i,k).lt.npif)nypif=nyipif(j,i,k)
!!$                   do jpif=1,nypif
!!$                      ia=ia+1
!!$                      if(izap.eq.1)then!zapping
!!$                         xa(ia)=yp(jy-jpif)!(jy-1)*dy-jpif*dy
!!$                         ya(ia)=u(i,jy-jpif,k)
!!$                      else             !no zapping
!!$                         xa(ia)=yp(jy-jpif+1)!(jy-1)*dy-(jpif-1)*dy
!!$                         ya(ia)=u(i,jy-jpif+1,k)
!!$                      endif
!!$                   enddo
!!$                else                   !objet semi-immergé
!!$                   jpoli=1
!!$                endif
!!$                !2ème frontière
!!$                nypif=npif
!!$                ia=ia+1
!!$                xa(ia)=yf(j,i,k)
!!$                if (nobjy(i,k).eq.2.and.j.eq.1) then
!!$                    ya(ia)=phiw2(2,j,i,k,is)
!!$                else
!!$                    ya(ia)=-one/nuw(is)
!!$                endif
!!$                if(yf(j,i,k).lt.yly)then!objet immergé
!!$                   jy=1!jy=(yf(j,i,k)+dy)/dy+1
!!$                   !====DEBUG
!!$                   !do while(yp(jy).le.yf(j,i,k)) !use if exact_ib_pipe
!!$                   do while(yp(jy).lt.yf(j,i,k))
!!$                      jy=jy+1
!!$                   enddo
!!$                   jpolf=jy-1
!!$                   if(nyfpif(j,i,k).lt.npif)nypif=nyfpif(j,i,k)
!!$                   do jpif=1,nypif
!!$                      ia=ia+1
!!$                      if(izap.eq.1)then!zapping
!!$                         xa(ia)=yp(jy+jpif)!(jy-1)*dy+jpif*dy
!!$                         ya(ia)=u(i,jy+jpif,k)
!!$                      else             !no zapping
!!$                         xa(ia)=yp(jy+jpif-1)!(jy-1)*dy+(jpif-1)*dy
!!$                         ya(ia)=u(i,jy+jpif-1,k)
!!$                      endif
!!$                   enddo
!!$                else                   !objet semi-immergé
!!$                   jpolf=ny
!!$                endif
!!$                !calcul du polynôme
!!$                na=ia
!!$                do jpol=jpoli,jpolf
!!$                   xpol=yp(jpol)!dy*(jpol-1)
!!$                   call polint(xa,ya,na,xpol,ypol,dypol)
!!$                   u(i,jpol,k)=ypol
!!$                enddo
!!$                !!====DEBUG Y
!!$                !!if (itime.eq.ilast) then
!!$                !    if (ystart(1)+i-1.eq.nx/2+1.and.ystart(3)+k-1.eq.nz/2+1) then
!!$                !        !do jj=1,ysize(2)
!!$                !        !    print*,'u_y:' ,yp(jj),u(i,jj,k)
!!$                !        !enddo
!!$                !        do ia=1,na
!!$                !            print*,'y:xa,ya:',xa(ia),ya(ia)
!!$                !        enddo
!!$                !    endif
!!$                !    !stop
!!$                !!endif
!!$                !!====
!!$                ia=0
!!$             enddo
!!$          endif
!!$       enddo
!!$    enddo
!!$    !
!!$    return
!!$  end subroutine nbclagpoly
!!$  !***************************************************************************
!!$  !
!!$  subroutine nbclagpolz(u,is)
!!$  !
!!$  !***************************************************************************
!!$    !
!!$    use param
!!$    use complex_geometry
!!$    use decomp_2d
!!$    use variables
!!$    !
!!$    implicit none
!!$    !
!!$    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: u
!!$    integer                                            :: is
!!$    integer                                            :: i,j,k
!!$    real(mytype)                                       :: x,y,z
!!$    real(mytype)                                       :: yw,zw,yc,zc,theta
!!$    real(mytype)                                       :: wfz             ! wall flux in Z  
!!$    integer                                            :: kz              != position du point "zappé"
!!$    integer                                            :: kzs
!!$    integer                                            :: kpif,kpol,nzpif
!!$    integer                                            :: kpoli,kpolf     != positions Initiales et Finales du POLynôme considéré
!!$    real(mytype)                                       :: xpol,ypol,dypol !|variables concernant les polynômes
!!$    real(mytype),dimension(10)                         :: xa,ya           !|de Lagrange. A mettre imérativement en 
!!$    integer                                            :: ia,na           !|double précision
!!$    !
!!$    if (new_rec.eq.1) then
!!$        call nbclagpolz2(u,is)
!!$        return
!!$    endif
!!$    !
!!$    yc=yly/two
!!$    zc=zlz/two
!!$    do j=1,zsize(2)
!!$       do i=1,zsize(1)
!!$          if(nobjz(i,j).ne.0)then
!!$             ia=0
!!$             do k=1,nobjz(i,j)          !boucle sur le nombre d'objets par couple (i,j)
!!$                !1ère frontière
!!$                nzpif=npif
!!$                ia=ia+1
!!$                xa(ia)=zi(k,i,j)
!!$                if (nobjz(i,j).eq.2.and.k.eq.2) then
!!$                    ya(ia)=phiw3(1,k,i,j,is)
!!$                else
!!$                    ya(ia)=-one/nuw(is)
!!$                endif
!!$                if(zi(k,i,j).gt.0.)then!objet immergé
!!$                   kz=zi(k,i,j)/dz+1
!!$                   !====DEBUG
!!$                   !if (abs(zp(kz)-zi(k,i,j)).lt.tol) kz=kz-1 !use if exact_ib_pipe
!!$                   kpoli=kz+1
!!$                   if(nzipif(k,i,j).lt.npif)nzpif=nzipif(k,i,j)
!!$                   do kpif=1,nzpif
!!$                      ia=ia+1
!!$                      if(izap.eq.1)then!zapping
!!$                         xa(ia)=(kz-1)*dz-kpif*dz
!!$                         ya(ia)=u(i,j,kz-kpif)
!!$                      else             !no zapping
!!$                         xa(ia)=(kz-1)*dz-(kpif-1)*dz
!!$                         ya(ia)=u(i,j,kz-kpif+1)
!!$                      endif
!!$                   enddo
!!$                else                   !objet semi-immergé
!!$                   kpoli=1
!!$                endif
!!$                !2ème frontière
!!$                nzpif=npif
!!$                ia=ia+1
!!$                xa(ia)=zf(k,i,j)
!!$                if (nobjz(i,j).eq.2.and.k.eq.1) then
!!$                    ya(ia)=phiw3(2,k,i,j,is)
!!$                else
!!$                    ya(ia)=-one/nuw(is)
!!$                endif
!!$                if(zf(k,i,j).lt.zlz)then!objet immergé
!!$                   kz=(zf(k,i,j)+dz)/dz+1
!!$                   kpolf=kz-1
!!$                   if(nzfpif(k,i,j).lt.npif)nzpif=nzfpif(k,i,j)
!!$                   do kpif=1,nzpif
!!$                      ia=ia+1
!!$                      if(izap.eq.1)then!zapping
!!$                         xa(ia)=(kz-1)*dz+kpif*dz
!!$                         ya(ia)=u(i,j,kz+kpif)
!!$                      else             !no zapping
!!$                         xa(ia)=(kz-1)*dz+(kpif-1)*dz
!!$                         ya(ia)=u(i,j,kz+kpif-1)
!!$                      endif
!!$                   enddo
!!$                else                   !objet semi-immergé
!!$                   kpolf=nz
!!$                endif
!!$                !calcul du polynôme
!!$                na=ia
!!$                do kpol=kpoli,kpolf
!!$                   xpol=dz*(kpol-1)
!!$                   call polint(xa,ya,na,xpol,ypol,dypol)
!!$                   u(i,j,kpol)=ypol
!!$                enddo
!!$                ia=0
!!$             enddo
!!$          endif
!!$       enddo
!!$    enddo
!!$    !
!!$    return
!!$  end subroutine nbclagpolz
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    !!
!!$    !!  subroutine: nbclagpol*2
!!$    !!      AUTHOR: Rodrigo Vicente Cruz
!!$    !! DESCRIPTION: Performs Neumann reconstruction across the transverse-yz
!!$    !!              domain periodicity
!!$    !!
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !***************************************************************************
!!$  !
!!$  subroutine nbclagpoly2(u,is)
!!$  !
!!$  !***************************************************************************
!!$    !
!!$    use param
!!$    use complex_geometry
!!$    use decomp_2d
!!$    use variables
!!$    !
!!$    implicit none
!!$    !
!!$    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: u
!!$    integer                                            :: is
!!$    integer                                            :: i,j,k
!!$    real(mytype)                                       :: x,y,z
!!$    real(mytype)                                       :: yw,zw,yc,zc,theta
!!$    real(mytype)                                       :: wfy             ! wall flux in Y  
!!$    integer                                            :: jy              != position du point "zappé"
!!$    integer                                            :: jys
!!$    integer                                            :: jpif,jpol,nypif
!!$    integer                                            :: jpoli,jpolf     != positions Initiales et Finales du POLynôme considéré
!!$    integer                                            :: nyd,jm
!!$    real(mytype)                                       :: xpol,ypol,dypol !|variables concernant les polynômes
!!$    real(mytype)                                       :: yd
!!$    real(mytype),dimension(10)                         :: xa,ya           !|de Lagrange. A mettre impérativement en 
!!$    integer                                            :: ia,na           !|double précision
!!$    !
!!$    yc=yly/two
!!$    zc=zlz/two
!!$    do k=1,ysize(3)
!!$        do i=1,ysize(1)
!!$            !if(nobjy(i,k).ne.0)then
!!$            if(nobjy(i,k).eq.2)then
!!$                ia=0
!!$                do j=1,nobjy(i,k)          !boucle sur le nombre d'objets par (j,k)
!!$                    !1ère frontière
!!$                    if (j.eq.2) then
!!$                        nypif=npif
!!$                        ia=ia+1
!!$                        xa(ia)=yi(j,i,k)
!!$                        yd=yly
!!$                        nyd=ny
!!$                        xa(ia)=xa(ia)-yd
!!$                        ya(ia)=phiw2(1,j,i,k,is)
!!$                        if(yi(j,i,k).gt.0.)then!objet immergé
!!$                            jy=1!jy=yi(j,i,k)/dy+1
!!$                            do while(yp(jy).lt.yi(j,i,k))
!!$                               jy=jy+1
!!$                            enddo
!!$                            jy=jy-1
!!$                            jpoli=jy+1
!!$                            jpoli=jpoli-nyd
!!$                            if(nyipif(j,i,k).lt.npif)nypif=nyipif(j,i,k)
!!$                            do jpif=1,nypif
!!$                               ia=ia+1
!!$                               if(izap.eq.1)then!zapping
!!$                                    !xa(ia)=yp(jy-jpif)!(jy-1)*dy-jpif*dy
!!$                                    xa(ia)=(jy-1)*dy-jpif*dy
!!$                                    ya(ia)=u(i,jy-jpif,k)
!!$                               else             !no zapping
!!$                                    !xa(ia)=yp(jy-jpif+1)!(jy-1)*dy-(jpif-1)*dy
!!$                                    xa(ia)=(jy-1)*dy-(jpif-1)*dy
!!$                                    ya(ia)=u(i,jy-jpif+1,k)
!!$                               endif
!!$                               xa(ia)=xa(ia)-yd
!!$                            enddo
!!$                        else                   !objet semi-immergé
!!$                           jpoli=1
!!$                        endif
!!$                    endif
!!$                    !2ème frontière
!!$                    if (j.eq.1) then
!!$                        nypif=npif
!!$                        ia=ia+1
!!$                        xa(ia)=yf(j,i,k)
!!$                        ya(ia)=phiw2(2,j,i,k,is)
!!$                        if(yf(j,i,k).lt.yly)then!objet immergé
!!$                           jy=1!jy=(yf(j,i,k)+dy)/dy+1
!!$                           !====DEBUG
!!$                           !do while(yp(jy).le.yf(j,i,k)) !use if exact_ib_pipe
!!$                           do while(yp(jy).lt.yf(j,i,k))
!!$                              jy=jy+1
!!$                           enddo
!!$                           jpolf=jy-1
!!$                           if(nyfpif(j,i,k).lt.npif)nypif=nyfpif(j,i,k)
!!$                           do jpif=1,nypif
!!$                              ia=ia+1
!!$                              if(izap.eq.1)then!zapping
!!$                                 !xa(ia)=yp(jy+jpif)!(jy-1)*dy+jpif*dy
!!$                                 xa(ia)=(jy-1)*dy+jpif*dy
!!$                                 ya(ia)=u(i,jy+jpif,k)
!!$                              else             !no zapping
!!$                                 !xa(ia)=yp(jy+jpif-1)!(jy-1)*dy+(jpif-1)*dy
!!$                                 xa(ia)=(jy-1)*dy+(jpif-1)*dy
!!$                                 ya(ia)=u(i,jy+jpif-1,k)
!!$                              endif
!!$                           enddo
!!$                        else                   !objet semi-immergé
!!$                           jpolf=ny
!!$                        endif
!!$                    endif
!!$                enddo
!!$                !calcul du polynôme
!!$                na=ia
!!$                do jpol=jpoli,jpolf
!!$                    !xpol=yp(jpol)!dy*(jpol-1)
!!$                    xpol=(jpol-1)*dy 
!!$                    call polint(xa,ya,na,xpol,ypol,dypol)
!!$                    if (jpoli.lt.0.and.jpol.le.0) then
!!$                        jm=jpol+nyd
!!$                        u(i,jm,k)=ypol
!!$                    else
!!$                        u(i,jpol,k)=ypol
!!$                    endif
!!$                enddo
!!$                !!====DEBUG Y
!!$                !!if (itime.eq.ilast) then
!!$                !    if (ystart(1)+i-1.eq.nx/2+1.and.ystart(3)+k-1.eq.nz/2+1) then
!!$                !        !do jj=1,ysize(2)
!!$                !        !    print*,'u_y:' ,yp(jj),u(i,jj,k)
!!$                !        !enddo
!!$                !        do ia=1,na
!!$                !            print*,'y:xa,ya:',xa(ia),ya(ia)
!!$                !        enddo
!!$                !    endif
!!$                !    !stop
!!$                !!endif
!!$                !!====
!!$                ia=0
!!$            else
!!$                do j=1,ysize(2)
!!$                    u(i,j,k)=-one/nuw(is)
!!$                enddo
!!$            endif
!!$        enddo
!!$    enddo
!!$    !
!!$    return
!!$  end subroutine nbclagpoly2
!!$  !***************************************************************************
!!$  !
!!$  subroutine nbclagpolz2(u,is)
!!$  !
!!$  !***************************************************************************
!!$    !
!!$    use param
!!$    use complex_geometry
!!$    use decomp_2d
!!$    use variables
!!$    !
!!$    implicit none
!!$    !
!!$    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: u
!!$    integer                                            :: is
!!$    integer                                            :: i,j,k
!!$    real(mytype)                                       :: x,y,z
!!$    real(mytype)                                       :: yw,zw,yc,zc,theta
!!$    real(mytype)                                       :: wfz             ! wall flux in Z  
!!$    integer                                            :: kz              != position du point "zappé"
!!$    integer                                            :: kzs
!!$    integer                                            :: kpif,kpol,nzpif
!!$    integer                                            :: kpoli,kpolf     != positions Initiales et Finales du POLynôme considéré
!!$    integer                                            :: nzd,km
!!$    real(mytype)                                       :: xpol,ypol,dypol !|variables concernant les polynômes
!!$    real(mytype)                                       :: zd
!!$    real(mytype),dimension(10)                         :: xa,ya           !|de Lagrange. A mettre imérativement en 
!!$    integer                                            :: ia,na           !|double précision
!!$    !
!!$    yc=yly/two
!!$    zc=zlz/two
!!$    do j=1,zsize(2)
!!$        do i=1,zsize(1)
!!$            !if(nobjz(i,j).ne.0)then
!!$            if(nobjz(i,j).eq.2)then
!!$                ia=0
!!$                do k=1,nobjz(i,j)          !boucle sur le nombre d'objets par couple (i,j)
!!$                    !1ère frontière
!!$                    if (k.eq.2) then
!!$                        nzpif=npif
!!$                        ia=ia+1
!!$                        xa(ia)=zi(k,i,j)
!!$                        zd=zlz
!!$                        nzd=nz
!!$                        xa(ia)=xa(ia)-zd
!!$                        ya(ia)=phiw3(1,k,i,j,is)
!!$                        if(zi(k,i,j).gt.0.)then!objet immergé
!!$                            kz=zi(k,i,j)/dz+1
!!$                            !====DEBUG
!!$                            !if (abs(zp(kz)-zi(k,i,j)).lt.tol) kz=kz-1 !use if exact_ib_pipe
!!$                            kpoli=kz+1
!!$                            kpoli=kpoli-nzd
!!$                            if(nzipif(k,i,j).lt.npif)nzpif=nzipif(k,i,j)
!!$                            do kpif=1,nzpif
!!$                                ia=ia+1
!!$                                if(izap.eq.1)then!zapping
!!$                                   xa(ia)=(kz-1)*dz-kpif*dz
!!$                                   ya(ia)=u(i,j,kz-kpif)
!!$                                else             !no zapping
!!$                                   xa(ia)=(kz-1)*dz-(kpif-1)*dz
!!$                                   ya(ia)=u(i,j,kz-kpif+1)
!!$                                endif
!!$                                xa(ia)=xa(ia)-zd 
!!$                            enddo
!!$                        else                   !objet semi-immergé
!!$                            kpoli=1
!!$                        endif
!!$                    endif
!!$                    !2ème frontière
!!$                    if (k.eq.1) then
!!$                        nzpif=npif
!!$                        ia=ia+1
!!$                        xa(ia)=zf(k,i,j)
!!$                        ya(ia)=phiw3(2,k,i,j,is)
!!$                        if(zf(k,i,j).lt.zlz)then!objet immergé
!!$                           kz=(zf(k,i,j)+dz)/dz+1
!!$                           kpolf=kz-1
!!$                           if(nzfpif(k,i,j).lt.npif)nzpif=nzfpif(k,i,j)
!!$                           do kpif=1,nzpif
!!$                              ia=ia+1
!!$                              if(izap.eq.1)then!zapping
!!$                                 xa(ia)=(kz-1)*dz+kpif*dz
!!$                                 ya(ia)=u(i,j,kz+kpif)
!!$                              else             !no zapping
!!$                                 xa(ia)=(kz-1)*dz+(kpif-1)*dz
!!$                                 ya(ia)=u(i,j,kz+kpif-1)
!!$                              endif
!!$                           enddo
!!$                        else                   !objet semi-immergé
!!$                           kpolf=nz
!!$                        endif
!!$                    endif
!!$                enddo
!!$                !calcul du polynôme
!!$                na=ia
!!$                do kpol=kpoli,kpolf
!!$                    xpol=(kpol-1)*dz
!!$                    call polint(xa,ya,na,xpol,ypol,dypol)
!!$                    if (kpoli.lt.0.and.kpol.le.0) then
!!$                        km=kpol+nzd
!!$                        u(i,j,km)=ypol
!!$                    else
!!$                        u(i,j,kpol)=ypol
!!$                    endif
!!$                enddo
!!$                !!====DEBUG Z
!!$                !!if (itime.eq.ilast) then
!!$                !    if (zstart(1)+i-1.eq.nx/2+1.and.zstart(2)+j-1.eq.ny/2+1) then
!!$                !        !do kk=1,zsize(3)
!!$                !        !    print*,'u_z:' ,(kk-1)*dz,u(i,j,kk)
!!$                !        !enddo
!!$                !        do ia=1,na
!!$                !            print*,'z:xa,ya:',xa(ia),ya(ia)
!!$                !        enddo
!!$                !    endif
!!$                !    !stop
!!$                !!endif
!!$                !!====
!!$                ia=0
!!$            else
!!$                do k=1,zsize(3)
!!$                    u(i,j,k)=-one/nuw(is)
!!$                enddo
!!$            endif
!!$        enddo
!!$    enddo
!!$    !
!!$    return
!!$  end subroutine nbclagpolz2
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    !!
!!$    !!  subroutine: chtlagpol*
!!$    !!      AUTHOR: Rodrigo Vicente Cruz
!!$    !! DESCRIPTION: When CHT BC are used, performs the reconstruction of fluid
!!$    !!              temperature solution by imposing phiw predicted beforehand 
!!$    !!              (from phiw_cht)
!!$    !!
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !***************************************************************************
!!$  !
!!$  subroutine chtlagpoly(u,is)
!!$  !
!!$  !***************************************************************************
!!$    !
!!$    use param
!!$    use complex_geometry
!!$    use decomp_2d
!!$    use variables
!!$    !use var, ONLY: phiw2
!!$    !
!!$    implicit none
!!$    !
!!$    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: u
!!$    integer                                            :: is
!!$    integer                                            :: i,j,k
!!$    real(mytype)                                       :: x,y,z
!!$    real(mytype)                                       :: yw,zw,yc,zc,theta
!!$    real(mytype)                                       :: wfy             ! wall flux in Y  
!!$    integer                                            :: jy              != position du point "zappé"
!!$    integer                                            :: jys
!!$    integer                                            :: jpif,jpol,nypif
!!$    integer                                            :: jpoli,jpolf     != positions Initiales et Finales du POLynôme considéré
!!$    real(mytype)                                       :: xpol,ypol,dypol !|variables concernant les polynômes
!!$    real(mytype),dimension(10)                         :: xa,ya           !|de Lagrange. A mettre impérativement en 
!!$    integer                                            :: ia,na           !|double précision
!!$    !
!!$    if (new_rec.eq.1) then
!!$        call chtlagpoly2(u,is)
!!$        return
!!$    endif
!!$    !
!!$    yc=yly/two
!!$    zc=zlz/two
!!$    do k=1,ysize(3)
!!$       do i=1,ysize(1)
!!$          if(nobjy(i,k).ne.0)then
!!$             ia=0
!!$             do j=1,nobjy(i,k)          !boucle sur le nombre d'objets par (j,k)
!!$                !1ère frontière
!!$                nypif=npif
!!$                ia=ia+1
!!$                xa(ia)=yi(j,i,k)
!!$                if (nobjy(i,k).eq.2.and.j.eq.2) then
!!$                    ya(ia)=phiw2(1,j,i,k,is)
!!$                else
!!$                    ya(ia)=-one/nuw(is)
!!$                endif
!!$                if(yi(j,i,k).gt.0.)then!objet immergé
!!$                   jy=1!jy=yi(j,i,k)/dy+1
!!$                   do while(yp(jy).lt.yi(j,i,k))
!!$                      jy=jy+1
!!$                   enddo
!!$                   jy=jy-1
!!$                   jpoli=jy+1
!!$                   if(nyipif(j,i,k).lt.npif)nypif=nyipif(j,i,k)
!!$                   do jpif=1,nypif
!!$                      ia=ia+1
!!$                      if(izap.eq.1)then!zapping
!!$                         xa(ia)=yp(jy-jpif)!(jy-1)*dy-jpif*dy
!!$                         ya(ia)=u(i,jy-jpif,k)
!!$                      else             !no zapping
!!$                         xa(ia)=yp(jy-jpif+1)!(jy-1)*dy-(jpif-1)*dy
!!$                         ya(ia)=u(i,jy-jpif+1,k)
!!$                      endif
!!$                   enddo
!!$                else                   !objet semi-immergé
!!$                   jpoli=1
!!$                endif
!!$                !2ème frontière
!!$                nypif=npif
!!$                ia=ia+1
!!$                xa(ia)=yf(j,i,k)
!!$                if (nobjy(i,k).eq.2.and.j.eq.1) then
!!$                    ya(ia)=phiw2(2,j,i,k,is)
!!$                else
!!$                    ya(ia)=-one/nuw(is)
!!$                endif
!!$                if(yf(j,i,k).lt.yly)then!objet immergé
!!$                   jy=1!jy=(yf(j,i,k)+dy)/dy+1
!!$                   !====DEBUG
!!$                   !do while(yp(jy).le.yf(j,i,k)) !use if exact_ib_pipe
!!$                   do while(yp(jy).lt.yf(j,i,k))
!!$                      jy=jy+1
!!$                   enddo
!!$                   jpolf=jy-1
!!$                   if(nyfpif(j,i,k).lt.npif)nypif=nyfpif(j,i,k)
!!$                   do jpif=1,nypif
!!$                      ia=ia+1
!!$                      if(izap.eq.1)then!zapping
!!$                         xa(ia)=yp(jy+jpif)!(jy-1)*dy+jpif*dy
!!$                         ya(ia)=u(i,jy+jpif,k)
!!$                      else             !no zapping
!!$                         xa(ia)=yp(jy+jpif-1)!(jy-1)*dy+(jpif-1)*dy
!!$                         ya(ia)=u(i,jy+jpif-1,k)
!!$                      endif
!!$                   enddo
!!$                else                   !objet semi-immergé
!!$                   jpolf=ny
!!$                endif
!!$                !calcul du polynôme
!!$                na=ia
!!$                do jpol=jpoli,jpolf
!!$                   xpol=yp(jpol)!dy*(jpol-1)
!!$                   call polint(xa,ya,na,xpol,ypol,dypol)
!!$                   u(i,jpol,k)=ypol
!!$                enddo
!!$                !!====DEBUG Y
!!$                !!if (itime.eq.ilast) then
!!$                !    if (ystart(1)+i-1.eq.nx/2+1.and.ystart(3)+k-1.eq.120) then
!!$                !        !do jj=1,ysize(2)
!!$                !        !    print*,'u_y:' ,yp(jj),u(i,jj,k)
!!$                !        !enddo
!!$                !        do ia=1,na
!!$                !            print*,'y:xa,ya:',xa(ia),ya(ia)
!!$                !        enddo
!!$                !    endif
!!$                !    !stop
!!$                !!endif
!!$                !!====
!!$                ia=0
!!$             enddo
!!$          endif
!!$       enddo
!!$    enddo
!!$    !
!!$    return
!!$  end subroutine chtlagpoly
!!$  !***************************************************************************
!!$  !
!!$  subroutine chtlagpolz(u,is)
!!$  !
!!$  !***************************************************************************
!!$    !
!!$    use param
!!$    use complex_geometry
!!$    use decomp_2d
!!$    use variables
!!$    !use var, ONLY: phiw3
!!$    !
!!$    implicit none
!!$    !
!!$    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: u
!!$    integer                                            :: is
!!$    integer                                            :: i,j,k
!!$    real(mytype)                                       :: x,y,z
!!$    real(mytype)                                       :: yw,zw,yc,zc,theta
!!$    real(mytype)                                       :: wfz             ! wall flux in Z  
!!$    integer                                            :: kz              != position du point "zappé"
!!$    integer                                            :: kzs
!!$    integer                                            :: kpif,kpol,nzpif
!!$    integer                                            :: kpoli,kpolf     != positions Initiales et Finales du POLynôme considéré
!!$    real(mytype)                                       :: xpol,ypol,dypol !|variables concernant les polynômes
!!$    real(mytype),dimension(10)                         :: xa,ya           !|de Lagrange. A mettre imérativement en 
!!$    integer                                            :: ia,na           !|double précision
!!$    !
!!$    if (new_rec.eq.1) then
!!$        call chtlagpolz2(u,is)
!!$        return
!!$    endif
!!$    !
!!$    yc=yly/two
!!$    zc=zlz/two
!!$    do j=1,zsize(2)
!!$       do i=1,zsize(1)
!!$          if(nobjz(i,j).ne.0)then
!!$             ia=0
!!$             do k=1,nobjz(i,j)          !boucle sur le nombre d'objets par couple (i,j)
!!$                !1ère frontière
!!$                nzpif=npif
!!$                ia=ia+1
!!$                xa(ia)=zi(k,i,j)
!!$                if (nobjz(i,j).eq.2.and.k.eq.2) then
!!$                    ya(ia)=phiw3(1,k,i,j,is)
!!$                else
!!$                    ya(ia)=-one/nuw(is)
!!$                endif
!!$                if(zi(k,i,j).gt.0.)then!objet immergé
!!$                   kz=zi(k,i,j)/dz+1
!!$                   !====DEBUG
!!$                   !if (abs(zp(kz)-zi(k,i,j)).lt.tol) kz=kz-1 !use if exact_ib_pipe
!!$                   kpoli=kz+1
!!$                   if(nzipif(k,i,j).lt.npif)nzpif=nzipif(k,i,j)
!!$                   do kpif=1,nzpif
!!$                      ia=ia+1
!!$                      if(izap.eq.1)then!zapping
!!$                         xa(ia)=(kz-1)*dz-kpif*dz
!!$                         ya(ia)=u(i,j,kz-kpif)
!!$                      else             !no zapping
!!$                         xa(ia)=(kz-1)*dz-(kpif-1)*dz
!!$                         ya(ia)=u(i,j,kz-kpif+1)
!!$                      endif
!!$                   enddo
!!$                else                   !objet semi-immergé
!!$                   kpoli=1
!!$                endif
!!$                !2ème frontière
!!$                nzpif=npif
!!$                ia=ia+1
!!$                xa(ia)=zf(k,i,j)
!!$                if (nobjz(i,j).eq.2.and.k.eq.1) then
!!$                    ya(ia)=phiw3(2,k,i,j,is)
!!$                else
!!$                    ya(ia)=-one/nuw(is)
!!$                endif
!!$                if(zf(k,i,j).lt.zlz)then!objet immergé
!!$                   kz=(zf(k,i,j)+dz)/dz+1
!!$                   kpolf=kz-1
!!$                   if(nzfpif(k,i,j).lt.npif)nzpif=nzfpif(k,i,j)
!!$                   do kpif=1,nzpif
!!$                      ia=ia+1
!!$                      if(izap.eq.1)then!zapping
!!$                         xa(ia)=(kz-1)*dz+kpif*dz
!!$                         ya(ia)=u(i,j,kz+kpif)
!!$                      else             !no zapping
!!$                         xa(ia)=(kz-1)*dz+(kpif-1)*dz
!!$                         ya(ia)=u(i,j,kz+kpif-1)
!!$                      endif
!!$                   enddo
!!$                else                   !objet semi-immergé
!!$                   kpolf=nz
!!$                endif
!!$                !calcul du polynôme
!!$                na=ia
!!$                do kpol=kpoli,kpolf
!!$                   xpol=dz*(kpol-1)
!!$                   call polint(xa,ya,na,xpol,ypol,dypol)
!!$                   u(i,j,kpol)=ypol
!!$                enddo
!!$                !!====DEBUG Z
!!$                !!if (itime.eq.ilast) then
!!$                !    if (zstart(1)+i-1.eq.nx/2+1.and.zstart(2)+j-1.eq.120) then
!!$                !        !do kk=1,zsize(3)
!!$                !        !    print*,'u_z:' ,(kk-1)*dz,u(i,j,kk)
!!$                !        !enddo
!!$                !        do ia=1,na
!!$                !            print*,'z:xa,ya:',xa(ia),ya(ia)
!!$                !        enddo
!!$                !    endif
!!$                !    !stop
!!$                !!endif
!!$                !!====
!!$                ia=0
!!$             enddo
!!$          endif
!!$       enddo
!!$    enddo
!!$    !
!!$    return
!!$  end subroutine chtlagpolz
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    !!
!!$    !!  subroutine: chtlagpol*2
!!$    !!      AUTHOR: Rodrigo Vicente Cruz
!!$    !! DESCRIPTION: Performs CHT reconstruction across transverse-yz domain
!!$    !!              periodicity
!!$    !!
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !***************************************************************************
!!$  !
!!$  subroutine chtlagpoly2(u,is)
!!$  !
!!$  !***************************************************************************
!!$    !
!!$    use param
!!$    use complex_geometry
!!$    use decomp_2d
!!$    use variables
!!$    !use var, ONLY: phiw2
!!$    !
!!$    implicit none
!!$    !
!!$    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: u
!!$    integer                                            :: is
!!$    integer                                            :: i,j,k
!!$    real(mytype)                                       :: x,y,z
!!$    real(mytype)                                       :: yw,zw,yc,zc,theta
!!$    real(mytype)                                       :: wfy             ! wall flux in Y  
!!$    integer                                            :: jy              != position du point "zappé"
!!$    integer                                            :: jys
!!$    integer                                            :: jpif,jpol,nypif
!!$    integer                                            :: jpoli,jpolf     != positions Initiales et Finales du POLynôme considéré
!!$    integer                                            :: nyd,jm
!!$    real(mytype)                                       :: xpol,ypol,dypol !|variables concernant les polynômes
!!$    real(mytype)                                       :: yd
!!$    real(mytype),dimension(10)                         :: xa,ya           !|de Lagrange. A mettre impérativement en 
!!$    integer                                            :: ia,na           !|double précision
!!$    integer                                            :: jj
!!$    !
!!$    yc=yly/two
!!$    zc=zlz/two
!!$    do k=1,ysize(3)
!!$       do i=1,ysize(1)
!!$          !if(nobjy(i,k).ne.0)then
!!$          if(nobjy(i,k).eq.2)then
!!$             ia=0
!!$             do j=1,nobjy(i,k)          !boucle sur le nombre d'objets par (j,k)
!!$                !1ère frontière
!!$                if (j.eq.2) then !reconstruction par periodicite
!!$                    nypif=npif
!!$                    ia=ia+1
!!$                    xa(ia)=yi(j,i,k)
!!$                    yd=yly
!!$                    nyd=ny
!!$                    xa(ia)=xa(ia)-yd
!!$                    ya(ia)=phiw2(1,j,i,k,is)
!!$                    if(yi(j,i,k).gt.0.)then!objet immergé
!!$                       jy=1!jy=yi(j,i,k)/dy+1
!!$                       do while(yp(jy).lt.yi(j,i,k))
!!$                          jy=jy+1
!!$                       enddo
!!$                       jy=jy-1
!!$                       jpoli=jy+1
!!$                       jpoli=jpoli-nyd
!!$                       if(nyipif(j,i,k).lt.npif)nypif=nyipif(j,i,k)
!!$                       do jpif=1,nypif
!!$                          ia=ia+1
!!$                          if(izap.eq.1)then!zapping
!!$                             xa(ia)=(jy-1)*dy-jpif*dy
!!$                             ya(ia)=u(i,jy-jpif,k)
!!$                          else             !no zapping
!!$                             xa(ia)=(jy-1)*dy-(jpif-1)*dy
!!$                             ya(ia)=u(i,jy-jpif+1,k)
!!$                          endif
!!$                          xa(ia)=xa(ia)-yd
!!$                       enddo
!!$                    else                   !objet semi-immergé
!!$                       jpoli=1
!!$                    endif
!!$                endif
!!$                !2ème frontière
!!$                if (j.eq.1) then
!!$                    nypif=npif
!!$                    ia=ia+1
!!$                    xa(ia)=yf(j,i,k)
!!$                    ya(ia)=phiw2(2,j,i,k,is)
!!$                    if(yf(j,i,k).lt.yly)then!objet immergé
!!$                       jy=1!jy=(yf(j,i,k)+dy)/dy+1
!!$                       !====DEBUG
!!$                       !do while(yp(jy).le.yf(j,i,k)) !use if exact_ib_pipe
!!$                       do while(yp(jy).lt.yf(j,i,k))
!!$                          jy=jy+1
!!$                       enddo
!!$                       jpolf=jy-1
!!$                       if(nyfpif(j,i,k).lt.npif)nypif=nyfpif(j,i,k)
!!$                       do jpif=1,nypif
!!$                          ia=ia+1
!!$                          if(izap.eq.1)then!zapping
!!$                             xa(ia)=(jy-1)*dy+jpif*dy
!!$                             ya(ia)=u(i,jy+jpif,k)
!!$                          else             !no zapping
!!$                             xa(ia)=(jy-1)*dy+(jpif-1)*dy
!!$                             ya(ia)=u(i,jy+jpif-1,k)
!!$                          endif
!!$                       enddo
!!$                    else                   !objet semi-immergé
!!$                       jpolf=ny
!!$                    endif
!!$                endif
!!$            enddo
!!$            !calcul du polynôme
!!$            na=ia
!!$            do jpol=jpoli,jpolf
!!$               !xpol=yp(jpol)!dy*(jpol-1)
!!$               xpol=(jpol-1)*dy
!!$               call polint(xa,ya,na,xpol,ypol,dypol)
!!$               if (jpoli.lt.0.and.jpol.le.0) then
!!$                   jm=jpol+nyd
!!$                   u(i,jm,k)=ypol
!!$               else
!!$                   u(i,jpol,k)=ypol
!!$               endif
!!$            enddo
!!$            !!====DEBUG Y
!!$            !if (is.eq.1.and.itime.eq.30) then
!!$            !    if (ystart(1)+i-1.eq.nx/2+1.and.ystart(3)+k-1.eq.nz/2+1) then
!!$            !        do jj=1,ysize(2)
!!$            !            print*,'phiy:' ,yp(jj),u(i,jj,k),itime
!!$            !        enddo
!!$            !        do ia=1,na
!!$            !            print*,'y:xa,ya:',xa(ia),ya(ia),itime
!!$            !        enddo
!!$            !    endif
!!$            !    !stop
!!$            !endif
!!$            !!====
!!$            ia=0
!!$          else
!!$              do j=1,ysize(2)
!!$                 u(i,j,k)=-one/nuw(is)
!!$              enddo
!!$          endif
!!$       enddo
!!$    enddo
!!$    !
!!$    return
!!$  end subroutine chtlagpoly2
!!$  !***************************************************************************
!!$  !
!!$  subroutine chtlagpolz2(u,is)
!!$  !
!!$  !***************************************************************************
!!$    !
!!$    use param
!!$    use complex_geometry
!!$    use decomp_2d
!!$    use variables
!!$    !use var, ONLY: phiw3
!!$    !
!!$    implicit none
!!$    !
!!$    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: u
!!$    integer                                            :: is
!!$    integer                                            :: i,j,k
!!$    real(mytype)                                       :: x,y,z
!!$    real(mytype)                                       :: yw,zw,yc,zc,theta
!!$    real(mytype)                                       :: wfz             ! wall flux in Z  
!!$    integer                                            :: kz              != position du point "zappé"
!!$    integer                                            :: kzs
!!$    integer                                            :: kpif,kpol,nzpif
!!$    integer                                            :: kpoli,kpolf     != positions Initiales et Finales du POLynôme considéré
!!$    integer                                            :: nzd,km
!!$    real(mytype)                                       :: xpol,ypol,dypol !|variables concernant les polynômes
!!$    real(mytype)                                       :: zd
!!$    real(mytype),dimension(10)                         :: xa,ya           !|de Lagrange. A mettre imérativement en 
!!$    integer                                            :: ia,na           !|double précision
!!$    integer                                            :: kk
!!$    !
!!$    yc=yly/two
!!$    zc=zlz/two
!!$    do j=1,zsize(2)
!!$        do i=1,zsize(1)
!!$            !if(nobjz(i,j).ne.0)then
!!$            if(nobjz(i,j).eq.2)then
!!$                ia=0
!!$                do k=1,nobjz(i,j)          !boucle sur le nombre d'objets par couple (i,j)
!!$                    !1ère frontière
!!$                    if (k.eq.2) then !reconstruction par periodicite
!!$                        nzpif=npif
!!$                        ia=ia+1
!!$                        xa(ia)=zi(k,i,j)
!!$                        zd=zlz
!!$                        nzd=nz
!!$                        xa(ia)=xa(ia)-zd
!!$                        ya(ia)=phiw3(1,k,i,j,is)
!!$                        if(zi(k,i,j).gt.0.)then!objet immergé
!!$                            kz=zi(k,i,j)/dz+1
!!$                            !====DEBUG
!!$                            !if (abs(zp(kz)-zi(k,i,j)).lt.tol) kz=kz-1 !use if exact_ib_pipe
!!$                            kpoli=kz+1
!!$                            kpoli=kpoli-nzd
!!$                            if(nzipif(k,i,j).lt.npif) nzpif=nzipif(k,i,j)
!!$                            do kpif=1,nzpif
!!$                               ia=ia+1
!!$                               if(izap.eq.1)then!zapping
!!$                                  xa(ia)=(kz-1)*dz-kpif*dz
!!$                                  ya(ia)=u(i,j,kz-kpif)
!!$                               else             !no zapping
!!$                                  xa(ia)=(kz-1)*dz-(kpif-1)*dz
!!$                                  ya(ia)=u(i,j,kz-kpif+1)
!!$                               endif
!!$                               xa(ia)=xa(ia)-zd
!!$                            enddo
!!$                        else                   !objet semi-immergé
!!$                            kpoli=1
!!$                        endif
!!$                    endif
!!$                    !2ème frontière
!!$                    if (k.eq.1) then
!!$                        nzpif=npif
!!$                        ia=ia+1
!!$                        xa(ia)=zf(k,i,j)
!!$                        ya(ia)=phiw3(2,k,i,j,is)
!!$                        if(zf(k,i,j).lt.zlz)then !objet immergé
!!$                           kz=(zf(k,i,j)+dz)/dz+1
!!$                           kpolf=kz-1
!!$                           if(nzfpif(k,i,j).lt.npif)nzpif=nzfpif(k,i,j)
!!$                           do kpif=1,nzpif
!!$                              ia=ia+1
!!$                              if(izap.eq.1)then!zapping
!!$                                 xa(ia)=(kz-1)*dz+kpif*dz
!!$                                 ya(ia)=u(i,j,kz+kpif)
!!$                              else             !no zapping
!!$                                 xa(ia)=(kz-1)*dz+(kpif-1)*dz
!!$                                 ya(ia)=u(i,j,kz+kpif-1)
!!$                              endif
!!$                           enddo
!!$                        else                   !objet semi-immergé
!!$                           kpolf=nz
!!$                        endif
!!$                    endif
!!$                enddo
!!$                !calcul du polynôme
!!$                na=ia
!!$                do kpol=kpoli,kpolf
!!$                    xpol=dz*(kpol-1)
!!$                    call polint(xa,ya,na,xpol,ypol,dypol)
!!$                    if (kpoli.lt.0.and.kpol.le.0) then
!!$                        km=kpol+nzd
!!$                        u(i,j,km)=ypol
!!$                    else
!!$                        u(i,j,kpol)=ypol
!!$                    endif
!!$                enddo
!!$                !!====DEBUG Z
!!$                !if (is.eq.1.and.itime.eq.30) then
!!$                !    if (zstart(1)+i-1.eq.nx/2+1.and.zstart(2)+j-1.eq.nz/2+1) then
!!$                !        do kk=1,zsize(3)
!!$                !            print*,'phiz:' ,(kk-1)*dz,u(i,j,kk),itime
!!$                !        enddo
!!$                !        do ia=1,na
!!$                !            print*,'z:xa,ya:',xa(ia),ya(ia),itime
!!$                !        enddo
!!$                !    endif
!!$                !    !stop
!!$                !endif
!!$                !!====
!!$                ia=0
!!$            else
!!$                do k=1,zsize(3)
!!$                    u(i,j,k)=-one/nuw(is)
!!$                enddo
!!$            endif
!!$        enddo
!!$    enddo
!!$    !!====DEBUG Z
!!$    !if (itime.eq.1) stop
!!$    !
!!$    return
!!$  end subroutine chtlagpolz2
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    !!
!!$    !!  subroutine: chtlagpol*_s
!!$    !!      AUTHOR: Rodrigo Vicente Cruz
!!$    !! DESCRIPTION: When CHT BC are used, performs the reconstruction of solid
!!$    !!              temperture solution by imposing phiw predicted beforehand 
!!$    !!              (from phiw_cht)
!!$    !!
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !***************************************************************************
!!$  !
!!$  subroutine chtlagpoly_s(u,icht) !Solid field | By: Rodrigo Vicente Cruz
!!$  !
!!$  !***************************************************************************
!!$    !
!!$    use param
!!$    use complex_geometry
!!$    use decomp_2d
!!$    use variables
!!$    !
!!$    implicit none
!!$    !
!!$    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: u
!!$    integer                                            :: icht
!!$    integer                                            :: i,j,k,jm
!!$    real(mytype)                                       :: x,y,z
!!$    real(mytype)                                       :: yd,yw,zw,yc,zc,theta
!!$    integer                                            :: jy              !skipped solid point
!!$    integer                                            :: jobj            !other object through domain boundary   
!!$    integer                                            :: jpif,jpol,nypif
!!$    integer                                            :: jpoli,jpolf     !initial and final polynomial positions       
!!$    real(mytype)                                       :: xpol,ypol,dypol !|variables for polynomial interpolation.
!!$    real(mytype),dimension(10)                         :: xa,ya           !|use double precision
!!$    integer                                            :: ia,na
!!$    !====DEBUG
!!$    integer                                            :: jj
!!$    !
!!$    yc=yly/two
!!$    zc=zlz/two
!!$    do k=1,ysize(3)
!!$       do i=1,ysize(1)
!!$          if(nobjy(i,k).ne.0)then
!!$             ia=0
!!$             do j=1,nobjy(i,k) !Objects here are the fluid regions
!!$                !1ère frontière
!!$                nypif=npifs
!!$                ia=ia+1
!!$                xa(ia)=yi(j,i,k)
!!$                ya(ia)=phiws2(1,j,i,k,icht)
!!$                jy=1
!!$                do while(yp(jy).lt.yi(j,i,k))
!!$                   jy=jy+1
!!$                enddo
!!$                jpolf=jy-1 !?
!!$                !jpolf=jy
!!$                !if(nyipif(j,i,k).lt.npifs)nypif=nyipif(j,i,k)
!!$                do jpif=1,nypif
!!$                   ia=ia+1
!!$                   if(izap.eq.1)then !skipping
!!$                      !xa(ia)=yp(jy+jpif)
!!$                      xa(ia)=(jy-1)*dy+jpif*dy
!!$                      ya(ia)=u(i,jy+jpif,k)
!!$                   else              !no skipping
!!$                      !xa(ia)=yp(jy+jpif-1)
!!$                      xa(ia)=(jy-1)*dy+(jpif-1)*dy
!!$                      ya(ia)=u(i,jy+jpif-1,k)
!!$                   endif
!!$                enddo
!!$                !2ème frontière
!!$                nypif=npifs
!!$                ia=ia+1
!!$                if (j.eq.1) jobj=2
!!$                if (j.eq.2) jobj=1
!!$                if (nobjy(i,k).eq.1) jobj=j
!!$                jy=1
!!$                !====DEBUG
!!$                !do while(yp(jy).le.yf(jobj,i,k)) !use if exact_ib_pipe
!!$                do while(yp(jy).lt.yf(jobj,i,k))
!!$                   jy=jy+1
!!$                enddo
!!$                jy=jy-1
!!$                jpoli=jy+1 !?
!!$                !jpoli=jy
!!$                if (jpoli.gt.jpolf) then !reconstruction à travers periodicité
!!$                    !yd=yly-dy
!!$                    !jpoli=jpoli-(ny-1)
!!$                    yd=yly
!!$                    jpoli=jpoli-ny
!!$                else
!!$                    yd=0.
!!$                endif
!!$                xa(ia)=yf(jobj,i,k)-yd
!!$                ya(ia)=phiws2(2,jobj,i,k,icht)
!!$                !if(nyfpif(jobj,i,k).lt.npifs)nypif=nyfpif(jobj,i,k)
!!$                do jpif=1,nypif
!!$                   ia=ia+1
!!$                   if(izap.eq.1)then!zapping
!!$                      !xa(ia)=yp(jy-jpif)-yd
!!$                      xa(ia)=(jy-1)*dy-jpif*dy-yd
!!$                      ya(ia)=u(i,jy-jpif,k)
!!$                   else             !no zapping
!!$                      !xa(ia)=yp(jy-jpif+1)-yd
!!$                      xa(ia)=(jy-1)*dy-(jpif+1)*dy-yd
!!$                      ya(ia)=u(i,jy-jpif+1,k)
!!$                   endif
!!$                enddo
!!$                !calcul du polynôme
!!$                na=ia
!!$                do jpol=jpoli,jpolf
!!$                    !xpol=yp(jpol)
!!$                    xpol=(jpol-1)*dy
!!$                    call polint(xa,ya,na,xpol,ypol,dypol)
!!$                    if (jpoli.lt.0.and.jpol.le.0) then
!!$                        jm=jpol+ny
!!$                        u(i,jm,k)=ypol
!!$                    else
!!$                        u(i,jpol,k)=ypol
!!$                    endif
!!$                enddo
!!$                !!====DEBUG Y
!!$                !!if (itime.eq.ifirst) then
!!$                !!if (mod(itime,100).eq.0.and.icht.eq.1) then
!!$                !if (icht.eq.1.and.itime.eq.30) then
!!$                !    if (ystart(1)+i-1.eq.nx/2+1.and.ystart(3)+k-1.eq.nz/2+1) then
!!$                !        if (j.eq.2) then
!!$                !            do jj=1,ysize(2)
!!$                !                print*,'phisy:' ,yp(jj),u(i,jj,k),itime
!!$                !            enddo
!!$                !        endif
!!$                !        do ia=1,na
!!$                !            print*,'y:xa,ya:',xa(ia),ya(ia),itime
!!$                !        enddo
!!$                !    endif
!!$                !    !stop
!!$                !endif
!!$                !!====
!!$                ia=0
!!$             enddo
!!$          endif
!!$       enddo
!!$    enddo
!!$    !
!!$    return
!!$  end subroutine chtlagpoly_s
!!$  !***************************************************************************
!!$  !
!!$  subroutine chtlagpolz_s(u,icht) !Solid field | By: Rodrigo Vicente Cruz
!!$  !
!!$  !***************************************************************************
!!$    !
!!$    use param
!!$    use complex_geometry
!!$    use decomp_2d
!!$    use variables
!!$    !
!!$    implicit none
!!$    !
!!$    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: u
!!$    integer                                            :: icht
!!$    integer                                            :: i,j,k,km
!!$    real(mytype)                                       :: x,y,z
!!$    real(mytype)                                       :: zd,yw,zw,yc,zc,theta
!!$    integer                                            :: kz              !skipped solid point
!!$    integer                                            :: kobj            !other object through domain boundary   
!!$    integer                                            :: kpif,kpol,nzpif
!!$    integer                                            :: kpoli,kpolf     != positions Initiales et Finales du POLynôme considéré
!!$    real(mytype)                                       :: xpol,ypol,dypol !|variables concernant les polynômes
!!$    real(mytype),dimension(10)                         :: xa,ya           !|de Lagrange. A mettre imérativement en 
!!$    integer                                            :: ia,na           !|double précision
!!$    !====DEBUG
!!$    integer                                            :: kk
!!$    !
!!$    yc=yly/two
!!$    zc=zlz/two
!!$    do j=1,zsize(2)
!!$       do i=1,zsize(1)
!!$          if(nobjz(i,j).ne.0)then
!!$             ia=0
!!$             do k=1,nobjz(i,j)   !Objects here are the fluid regions
!!$                !1ère frontière
!!$                nzpif=npifs
!!$                ia=ia+1
!!$                xa(ia)=zi(k,i,j)
!!$                ya(ia)=phiws3(1,k,i,j,icht)
!!$                kz=zi(k,i,j)/dz+1
!!$                !====DEBUG
!!$                !if (abs(zp(kz)-zi(k,i,j)).lt.tol) kz=kz-1 !use if exact_ib_pipe
!!$                kz=kz+1
!!$                kpolf=kz-1 !?
!!$                !kpolf=kz
!!$                !if(nzipif(k,i,j).lt.npif)nzpif=nzipif(k,i,j)
!!$                do kpif=1,nzpif
!!$                   ia=ia+1
!!$                   if(izap.eq.1)then !skipping
!!$                      xa(ia)=(kz-1)*dz+kpif*dz
!!$                      ya(ia)=u(i,j,kz+kpif)
!!$                   else              !no skipping
!!$                      xa(ia)=(kz-1)*dz+(kpif-1)*dz
!!$                      ya(ia)=u(i,j,kz+kpif-1)
!!$                   endif
!!$                enddo
!!$                !2ème frontière
!!$                nzpif=npifs
!!$                ia=ia+1
!!$                if (k.eq.1) kobj=2
!!$                if (k.eq.2) kobj=1
!!$                if (nobjz(i,j).eq.1) kobj=k
!!$                kz=(zf(kobj,i,j)+dz)/dz+1
!!$                kz=kz-1
!!$                kpoli=kz+1 !?
!!$                !kpoli=kz
!!$                if (kpoli.gt.kpolf) then !reconstruction across periodicity
!!$                    zd=zlz
!!$                    kpoli=kpoli-nz
!!$                else
!!$                    zd=0.
!!$                endif
!!$                xa(ia)=zf(kobj,i,j)-zd
!!$                ya(ia)=phiws3(2,kobj,i,j,icht)
!!$                !if(nzipif(kobj,i,j).lt.npif)nzpif=nzipif(kobj,i,j)
!!$                do kpif=1,nzpif
!!$                   ia=ia+1
!!$                   if(izap.eq.1)then !skipping
!!$                      xa(ia)=(kz-1)*dz-kpif*dz-zd
!!$                      ya(ia)=u(i,j,kz-kpif)
!!$                   else              !no skipping
!!$                      xa(ia)=(kz-1)*dz-(kpif-1)*dz-zd
!!$                      ya(ia)=u(i,j,kz-kpif+1)
!!$                   endif
!!$                enddo
!!$                !calcul du polynôme
!!$                na=ia
!!$                do kpol=kpoli,kpolf
!!$                    xpol=(kpol-1)*dz
!!$                    call polint(xa,ya,na,xpol,ypol,dypol)
!!$                    if (kpoli.lt.0.and.kpol.le.0) then
!!$                        km=kpol+nz
!!$                        u(i,j,km)=ypol
!!$                    else
!!$                        u(i,j,kpol)=ypol
!!$                    endif
!!$                enddo
!!$                !!====DEBUG Z
!!$                !!if (itime.eq.ilast) then
!!$                !if (icht.eq.1.and.itime.eq.30) then
!!$                !    if (zstart(1)+i-1.eq.nx/2+1.and.zstart(2)+j-1.eq.ny/2+1) then
!!$                !        if (k.eq.2) then
!!$                !            do kk=1,zsize(3)
!!$                !                print*,'phisz:' ,(kk-1)*dz,u(i,j,kk),itime
!!$                !            enddo
!!$                !        endif
!!$                !        do ia=1,na
!!$                !            print*,'z:xa,ya:',xa(ia),ya(ia),itime
!!$                !        enddo
!!$                !    endif
!!$                !    !stop
!!$                !endif
!!$                !!====
!!$                ia=0
!!$             enddo
!!$          endif
!!$       enddo
!!$    enddo
!!$    !
!!$    return
!!$  end subroutine chtlagpolz_s
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    !!
!!$    !!  subroutine: axial_averaging
!!$    !!      AUTHOR: Rodrigo Vicente Cruz
!!$    !! DESCRIPTION: Performs axial averaging (streamwise-x direction)
!!$    !!
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !***************************************************************************
!!$  !
!!$  subroutine axial_averaging(var)
!!$  !
!!$  !***************************************************************************
!!$  use param
!!$  use variables
!!$  use decomp_2d
!!$  use MPI
!!$  use ibm
!!$  !
!!$  implicit none
!!$  real(mytype),dimension(xsize(1),xsize(2),xsize(3))    :: var
!!$  !LOCALS
!!$  real(mytype),dimension(xsize(2),xsize(3))             :: varm
!!$  integer                                               :: i,j,k
!!$  !
!!$
!!$  varm(:,:)=zero
!!$  do k=1,xsize(3)
!!$      do j=1,xsize(2)
!!$          do i=1,xsize(1)
!!$              varm(j,k) = varm(j,k) + var(i,j,k)
!!$          enddo
!!$      enddo
!!$  enddo
!!$  varm(:,:) = varm(:,:)/dble(xsize(1))
!!$  do k=1,xsize(3)
!!$      do j=1,xsize(2)
!!$          do i=1,xsize(1)
!!$              var(i,j,k) = varm(j,k)
!!$          enddo
!!$      enddo
!!$  enddo
!!$  !
!!$  return
!!$  end subroutine axial_averaging
!!$
end module pipe
