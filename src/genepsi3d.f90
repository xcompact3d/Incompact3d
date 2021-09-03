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
module genepsi

  public

contains

  subroutine epsi_init(ep1)

    USE param, only : zero, one, dx, dz
    USE decomp_2d, only : xstart, xend, xsize, mytype, nrank, decomp_2d_abort
    !USE decomp_2d_io
    USE variables, only : yp, ny

    implicit none
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    logical :: dir_exists

#ifdef DEBG
    if (nrank  ==  0) print *,'# body_init start'
#endif
    !###################################################################
    ! Check if geometry folder exists
    !###################################################################
    if (nrank==0) then
      inquire(file="geometry", exist=dir_exists)
      if (.not.dir_exists) then
        call system("mkdir geometry 2> /dev/null")
      end if
    end if
    !###################################################################
    ep1(:,:,:)=zero
    call geomcomplex(ep1,xstart(1),xend(1),ny,xstart(2),xend(2),xstart(3),xend(3),dx,yp,dz,one)

#ifdef DEBG
    if (nrank  ==  0) print *,'# body_init done'
#endif

    return
  end subroutine epsi_init
!############################################################################
!############################################################################
  subroutine geomcomplex(epsi, nxi, nxf, ny, nyi, nyf, nzi, nzf, dx, yp, dz, remp)

    USE param, ONLY : itype, itype_cyl, itype_hill, itype_channel
    USE decomp_2d, ONLY : mytype
    USE cyl, ONLY : geomcomplex_cyl
    USE hill, ONLY : geomcomplex_hill
    USE channel, ONLY : geomcomplex_channel

    IMPLICIT NONE

    INTEGER :: nxi,nxf,ny,nyi,nyf,nzi,nzf
    REAL(mytype),DIMENSION(nxi:nxf,nyi:nyf,nzi:nzf) :: epsi
    REAL(mytype)               :: dx,dz
    REAL(mytype),DIMENSION(ny) :: yp
    REAL(mytype)               :: remp

    IF (itype.EQ.itype_cyl) THEN

       CALL geomcomplex_cyl(epsi, nxi, nxf, ny, nyi, nyf, nzi, nzf, dx, yp, remp)

    ELSEIF (itype.EQ.itype_hill) THEN

       CALL geomcomplex_hill(epsi,nxi,nxf,ny,nyi,nyf,nzi,nzf,dx,yp,dz,remp)

    ELSEIF (itype.EQ.itype_channel) THEN

       CALL geomcomplex_channel(epsi, nxi, nxf, ny, nyi, nyf, nzi, nzf, yp, remp)

    ENDIF

  end subroutine geomcomplex
!############################################################################
!############################################################################
  subroutine genepsi3d(ep1)

    USE variables, only : nx,ny,nz,nxm,nym,nzm,yp
    USE param, only : xlx,yly,zlz,dx,dy,dz,izap,npif,nclx,ncly,nclz,istret
    USE complex_geometry
    use decomp_2d

    implicit none

    !*****************************************************************!
    ! 0- This program will generate all the files necessary for our
    !    customize IMB based on Lagrange reconstructions
    ! 3- The object is defined in the cylinder subroutine
    ! 4- You can add your own subroutine for your own object
    ! 7- Please cite the following paper if you are using this file:
    ! Gautier R., Laizet S. & Lamballais E., 2014, A DNS study of
    ! jet control with microjets using an alterna ng direc on forcing
    ! strategy, Int. J. of Computa onal Fluid Dynamics, 28, 393--410
    !*****************************************************************!
    !
    logical :: dir_exists
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    !
    if (nrank==0) then
      print *,'==========================================================='
      print *,'Generating the geometry!'
    end if
    !###################################################################
    ! Check if planes folder exists
    !###################################################################
    if (nrank==0) then
      inquire(file="geometry", exist=dir_exists)
      if (.not.dir_exists) then
        call system("mkdir geometry 2> /dev/null")
      end if
    end if
    !###################################################################
    call gene_epsi_3D(ep1,nx,ny,nz,dx,dy,dz,xlx,yly,zlz ,&
         nclx,ncly,nclz,nxraf,nyraf,nzraf   ,&
         xi,xf,yi,yf,zi,zf,nobjx,nobjy,nobjz,&
         nobjmax,yp,nraf)
    call verif_epsi(ep1,npif,izap,nx,ny,nz,nobjmax,&
         nxipif,nxfpif,nyipif,nyfpif,nzipif,nzfpif)
 !   call write_geomcomplex(nx,ny,nz,ep1,nobjx,nobjy,nobjz,xi,xf,yi,yf,zi,zf,&
 !        nxipif,nxfpif,nyipif,nyfpif,nzipif,nzfpif,nobjmax,npif)
    !
  end subroutine genepsi3d
!
!############################################################################
!############################################################################
  subroutine gene_epsi_3D(ep1,nx,ny,nz,dx,dy,dz,xlx,yly,zlz ,&
       nclx,ncly,nclz,nxraf,nyraf,nzraf   ,&
       xi,xf,yi,yf,zi,zf,nobjx,nobjy,nobjz,&
       nobjmax,yp,nraf)
    use param, only : zero,one
    use decomp_2d
    use MPI
    implicit none
    !
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ep1,smoofun,fbcx,fbcy,fbcz
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ep2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ep3
    integer                                            :: nx,ny,nz,nobjmax
    real(mytype)                                       :: dx,dy,dz
    real(mytype)                                       :: xlx,yly,zlz
    logical                                            :: nclx,ncly,nclz
    integer                                            :: nxraf,nyraf,nzraf
    integer                                            :: nraf
    integer,     dimension(xsize(2),xsize(3))          :: nobjx,nobjxraf
    integer,     dimension(ysize(1),ysize(3))          :: nobjy,nobjyraf
    integer,     dimension(zsize(1),zsize(2))          :: nobjz,nobjzraf
    real(mytype),dimension(nobjmax,xsize(2),xsize(3))  :: xi,xf
    real(mytype),dimension(nobjmax,ysize(1),ysize(3))  :: yi,yf
    real(mytype),dimension(nobjmax,zsize(1),zsize(2))  :: zi,zf
    real(mytype),dimension(nxraf,xsize(2),xsize(3))    :: xepsi
    real(mytype),dimension(ysize(1),nyraf,ysize(3))    :: yepsi
    real(mytype),dimension(zsize(1),zsize(2),nzraf)    :: zepsi
    real(mytype),dimension(ny)                         :: yp
    real(mytype),dimension(nyraf)                      :: ypraf
    real(mytype)                     :: dxraf,dyraf,dzraf
    integer                          :: i,j,k
    integer                          :: ii,jj,kk
    real(mytype)                     :: x,y,z
    integer                          :: inum,jnum,knum
    integer                          :: ibug,jbug,kbug
    integer                          :: iobj,jobj,kobj
    integer                          :: iflu,jflu,kflu
    integer                          :: isol,jsol,ksol
    integer                          :: iraf,jraf,kraf
    integer                          :: nobjxmax ,nobjymax ,nobjzmax
    integer                          :: nobjxmaxraf,nobjymaxraf,nobjzmaxraf
    integer                          :: idebraf,jdebraf,kdebraf
    integer                          :: ifinraf,jfinraf,kfinraf
    character(len=4) suffixe
    integer                          :: numvis
    integer                          :: mpi_aux_i, code

    !x-pencil
    ep1=zero
    call geomcomplex(ep1,xstart(1),xend(1),ny,xstart(2),xend(2),xstart(3),xend(3),dx,yp,dz,one)
  !  if (nrank==0) print*,'    step 1'
    if(nclx)then
       dxraf =xlx/nxraf
    else
       dxraf =xlx/(nxraf-1)
    endif
    xepsi=zero
    call geomcomplex(xepsi,1,nxraf,ny,xstart(2),xend(2),xstart(3),xend(3),dxraf,yp,dz,one)
  !  if (nrank==0) print*,'    step 2'
    !y-pencil
    if(ncly)then
       dyraf =yly/nyraf
    else
       dyraf =yly/(nyraf-1)
    endif
    do j=1,ny-1
       do jraf=1,nraf
          ypraf(jraf+nraf*(j-1))=yp(j)+(jraf-1)*(yp(j+1)-yp(j))/nraf
       enddo
    enddo
    if(.not.ncly)ypraf(nyraf)=yp(ny)
    yepsi=zero
    call geomcomplex(yepsi,ystart(1),yend(1),nyraf,1,nyraf,ystart(3),yend(3),dx,ypraf,dz,one)
  !  if (nrank==0) print*,'    step 3'
    !z-pencil
    if(nclz)then
       dzraf=zlz/nzraf
    else
       dzraf=zlz/(nzraf-1)
    endif
    zepsi=zero
    call geomcomplex(zepsi,zstart(1),zend(1),ny,zstart(2),zend(2),1,nzraf,dx,yp,dzraf,one)
  !  if (nrank==0) print*,'    step 4'

    !x-pencil
    nobjx(:,:)=0
    nobjxmax=0
    do k=1,xsize(3)
       do j=1,xsize(2)
          inum=0
          if(ep1(1,j,k) == 1.)then
             inum=1
             nobjx(j,k)=1
          endif
          do i=1,nx-1
             if(ep1(i,j,k) == 0..and.ep1(i+1,j,k) == 1.)then
                inum=inum+1
                nobjx(j,k)=nobjx(j,k)+1
             endif
          enddo
          if(inum > nobjxmax)then
             nobjxmax=inum
          endif
       enddo
    enddo
    call MPI_REDUCE(nobjxmax,mpi_aux_i,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_REDUCE")
  !  if (nrank==0) print*,'        nobjxmax=',mpi_aux_i

    nobjxraf(:,:)=0
    ibug=0
    nobjxmaxraf=0
    inum=0
    do k=1,xsize(3)
       do j=1,xsize(2)
          inum=0
          if(xepsi(1,j,k) == 1.)then
             inum=1
             nobjxraf(j,k)=1
          endif
          do i=1,nxraf-1
             if(xepsi(i,j,k) == 0..and.xepsi(i+1,j,k) == 1.)then
                inum=inum+1
                nobjxraf(j,k)=nobjxraf(j,k)+1
             endif
          enddo
          if(inum > nobjxmaxraf)then
             nobjxmaxraf=inum
          endif
          if(nobjx(j,k) /= nobjxraf(j,k))then
             ibug=ibug+1
          endif
       enddo
    enddo
    call MPI_REDUCE(nobjxmaxraf,mpi_aux_i,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_REDUCE")
  !  if (nrank==0) print*,'        nobjxmaxraf=',mpi_aux_i
    call MPI_REDUCE(ibug,mpi_aux_i,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_REDUCE")
  !  if (nrank==0) print*,'        ibug=',mpi_aux_i
  !  if (nrank==0) print*,'    step 5'

    !y-pencil
    nobjy(:,:)=0
    nobjymax=0
    call transpose_x_to_y(ep1,ep2)
    do k=1,ysize(3)
       do i=1,ysize(1)
          jnum=0
          if(ep2(i,1,k) == 1.)then
             jnum=1
             nobjy(i,k)=1
          endif
          do j=1,ny-1
             if(ep2(i,j,k) == 0..and.ep2(i,j+1,k) == 1.)then
                jnum=jnum+1
                nobjy(i,k)=nobjy(i,k)+1
             endif
          enddo
          if(jnum > nobjymax)then
             nobjymax=jnum
          endif
       enddo
    enddo
    call MPI_REDUCE(nobjymax,mpi_aux_i,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_REDUCE")
  !  if (nrank==0) print*,'        nobjymax=',mpi_aux_i

    nobjyraf(:,:)=0
    jbug=0
    nobjymaxraf=0
    jnum=0
    do k=1,ysize(3)
       do i=1,ysize(1)
          jnum=0
          if(yepsi(i,1,k) == 1.)then
             jnum=1
             nobjyraf(i,k)=1
          endif
          do j=1,nyraf-1
             if(yepsi(i,j,k) == 0..and.yepsi(i,j+1,k) == 1.)then
                jnum=jnum+1
                nobjyraf(i,k)=nobjyraf(i,k)+1
             endif
          enddo
          if(jnum > nobjymaxraf)then
             nobjymaxraf=jnum
          endif
          if(nobjy(i,k) /= nobjyraf(i,k))then
             jbug=jbug+1
          endif
       enddo
    enddo
    call MPI_REDUCE(nobjymaxraf,mpi_aux_i,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_REDUCE")
  !  if (nrank==0) print*,'        nobjymaxraf=',mpi_aux_i
    call MPI_REDUCE(jbug,mpi_aux_i,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_REDUCE")
  !  if (nrank==0) print*,'        jbug=',mpi_aux_i
  !  if (nrank==0) print*,'    step 6'

    !z-pencil
    nobjz(:,:)=0
    nobjzmax=0
    call transpose_y_to_z(ep2,ep3)
    do j=1,zsize(2)
       do i=1,zsize(1)
          knum=0
          if(ep3(i,j,1) == 1.)then
             knum=1
             nobjz(i,j)=1
          endif
          do k=1,nz-1
             if(ep3(i,j,k) == 0..and.ep3(i,j,k+1) == 1.)then
                knum=knum+1
                nobjz(i,j)=nobjz(i,j)+1
             endif
          enddo
          if(knum > nobjzmax)then
             nobjzmax=knum
          endif
       enddo
    enddo
    call MPI_REDUCE(nobjzmax,mpi_aux_i,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_REDUCE")
  !  if (nrank==0) print*,'        nobjzmax=',mpi_aux_i

    nobjzraf(:,:)=0
    kbug=0
    nobjzmaxraf=0
    knum=0
    do j=1,zsize(2)
       do i=1,zsize(1)
          knum=0
          if(zepsi(i,j,1) == 1.)then
             knum=1
             nobjzraf(i,j)=1
          endif
          do k=1,nzraf-1
             if(zepsi(i,j,k) == 0..and.zepsi(i,j,k+1) == 1.)then
                knum=knum+1
                nobjzraf(i,j)=nobjzraf(i,j)+1
             endif
          enddo
          if(knum > nobjzmaxraf)then
             nobjzmaxraf=knum
          endif
          if(nobjz(i,j) /= nobjzraf(i,j))then
             kbug=kbug+1
          endif
       enddo
    enddo
    call MPI_REDUCE(nobjzmaxraf,mpi_aux_i,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_REDUCE")
  !  if (nrank==0) print*,'        nobjzmaxraf=',mpi_aux_i
    call MPI_REDUCE(kbug,mpi_aux_i,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_REDUCE")
  !  if (nrank==0) print*,'        kbug=',mpi_aux_i
  !  if (nrank==0) print*,'    step 7'

    !x-pencil
    do k=1,xsize(3)
       do j=1,xsize(2)
          inum=0
          if(xepsi(1,j,k) == 1.)then
             inum=inum+1
             xi(inum,j,k)=-dx!-xlx
          endif
          do i=1,nxraf-1
             if(xepsi(i,j,k) == 0..and.xepsi(i+1,j,k) == 1.)then
                inum=inum+1
                xi(inum,j,k)=dxraf*(i-1)+dxraf/2.
             elseif(xepsi(i,j,k) == 1..and.xepsi(i+1,j,k) == 0.)then
                xf(inum,j,k)=dxraf*(i-1)+dxraf/2.
             endif
          enddo
          if(xepsi(nxraf,j,k) == 1.)then
             xf(inum,j,k)=xlx+dx!2.*xlx
          endif
       enddo
    enddo

  if(ibug /= 0)then
     do k=1,xsize(3)
        do j=1,xsize(2)
           if(nobjx(j,k) /= nobjxraf(j,k))then
              iobj=0
              if(ep1(1,j,k) == 1.)iobj=iobj+1
              do i=1,nx-1
                 if(ep1(i,j,k) == 0..and.ep1(i+1,j,k) == 1.)iobj=iobj+1
                 if(ep1(i,j,k) == 0..and.ep1(i+1,j,k) == 0.)iflu=1
                 if(ep1(i,j,k) == 1..and.ep1(i+1,j,k) == 1.)isol=1
                 do iraf=1,nraf
                    if(xepsi(iraf+nraf*(i-1)  ,j,k) == 0..and.&
                         xepsi(iraf+nraf*(i-1)+1,j,k) == 1.)idebraf=iraf+nraf*(i-1)+1
                    if(xepsi(iraf+nraf*(i-1)  ,j,k) == 1..and.&
                         xepsi(iraf+nraf*(i-1)+1,j,k) == 0.)ifinraf=iraf+nraf*(i-1)+1
                 enddo
                 if(idebraf /= 0.and.ifinraf /= 0.and.&
                      idebraf < ifinraf.and.iflu == 1)then
                    iobj=iobj+1
                    do ii=iobj,nobjmax-1
                       xi(ii,j,k)=xi(ii+1,j,k)
                       xf(ii,j,k)=xf(ii+1,j,k)
                    enddo
                    iobj=iobj-1
                 endif
                 if(idebraf /= 0.and.ifinraf /= 0.and.&
                      idebraf > ifinraf.and.isol == 1)then
                    iobj=iobj+1
                    do ii=iobj,nobjmax-1
                       xi(ii,j,k)=xi(ii+1,j,k)
                    enddo
                    iobj=iobj-1
                    do ii=iobj,nobjmax-1
                       xf(ii,j,k)=xf(ii+1,j,k)
                    enddo
                 endif
                 idebraf=0
                 ifinraf=0
                 iflu=0
              enddo
           endif
        enddo
     enddo
  endif
!  if (nrank==0) print*,'    step 8'

    !y-pencil
    do k=1,ysize(3)
       do i=1,ysize(1)
          jnum=0
          if(yepsi(i,1,k) == 1.)then
             jnum=jnum+1
             yi(jnum,i,k)=-(yp(2)-yp(1))!-yly
          endif
          do j=1,nyraf-1
             if(yepsi(i,j,k) == 0..and.yepsi(i,j+1,k) == 1.)then
                jnum=jnum+1
                yi(jnum,i,k)=ypraf(j)+(ypraf(j+1)-ypraf(j))/2.!dyraf*(j-1)+dyraf/2.
             elseif(yepsi(i,j,k) == 1..and.yepsi(i,j+1,k) == 0.)then
                yf(jnum,i,k)=ypraf(j)+(ypraf(j+1)-ypraf(j))/2.!dyraf*(j-1)+dyraf/2.
             endif
          enddo
          if(yepsi(i,nyraf,k) == 1.)then
             yf(jnum,i,k)=yly+(yp(ny)-yp(ny-1))/2.!2.*yly
          endif
       enddo
    enddo

  if(jbug /= 0)then
     do k=1,ysize(3)
        do i=1,ysize(1)
           if(nobjy(i,k) /= nobjyraf(i,k))then
              jobj=0
              if(ep2(i,1,k) == 1.)jobj=jobj+1
              do j=1,ny-1
                 if(ep2(i,j,k) == 0..and.ep2(i,j+1,k) == 1.)jobj=jobj+1
                 if(ep2(i,j,k) == 0..and.ep2(i,j+1,k) == 0.)jflu=1
                 if(ep2(i,j,k) == 1..and.ep2(i,j+1,k) == 1.)jsol=1
                 do jraf=1,nraf
                    if(yepsi(i,jraf+nraf*(j-1)  ,k) == 0..and.&
                         yepsi(i,jraf+nraf*(j-1)+1,k) == 1.)jdebraf=jraf+nraf*(j-1)+1
                    if(yepsi(i,jraf+nraf*(j-1)  ,k) == 1..and.&
                         yepsi(i,jraf+nraf*(j-1)+1,k) == 0.)jfinraf=jraf+nraf*(j-1)+1
                 enddo
                 if(jdebraf /= 0.and.jfinraf /= 0.and.&
                      jdebraf < jfinraf.and.jflu == 1)then
                    jobj=jobj+1
                    do jj=jobj,nobjmax-1
                       yi(jj,i,k)=yi(jj+1,i,k)
                       yf(jj,i,k)=yf(jj+1,i,k)
                    enddo
                    jobj=jobj-1
                 endif
                 if(jdebraf /= 0.and.jfinraf /= 0.and.&
                      jdebraf > jfinraf.and.jsol == 1)then
                    jobj=jobj+1
                    do jj=jobj,nobjmax-1
                       yi(jj,i,k)=yi(jj+1,i,k)
                    enddo
                    jobj=jobj-1
                    do jj=jobj,nobjmax-1
                       yf(jj,i,k)=yf(jj+1,i,k)
                    enddo
                 endif
                 jdebraf=0
                 jfinraf=0
                 jflu=0
              enddo
           endif
        enddo
     enddo
  endif
!  if (nrank==0) print*,'    step 9'

    !z-pencil
    do j=1,zsize(2)
       do i=1,zsize(1)
          knum=0
          if(zepsi(i,j,1) == 1.)then
             knum=knum+1
             zi(knum,i,j)=-dz!zlz
          endif
          do k=1,nzraf-1
             if(zepsi(i,j,k) == 0..and.zepsi(i,j,k+1) == 1.)then
                knum=knum+1
                zi(knum,i,j)=dzraf*(k-1)+dzraf/2.
             elseif(zepsi(i,j,k) == 1..and.zepsi(i,j,k+1) == 0.)then
                zf(knum,i,j)=dzraf*(k-1)+dzraf/2.
             endif
          enddo
          if(zepsi(i,j,nzraf) == 1.)then
             zf(knum,i,j)=zlz+dz!2.*zlz
          endif
       enddo
    enddo

kdebraf=0.
  if(kbug /= 0)then
     do j=1,zsize(2)
        do i=1,zsize(1)
           if(nobjz(i,j) /= nobjzraf(i,j))then
              kobj=0
              if(ep3(i,j,1) == 1.)kobj=kobj+1
              do k=1,nz-1
                 if(ep3(i,j,k) == 0..and.ep3(i,j,k+1) == 1.)kobj=kobj+1
                 if(ep3(i,j,k) == 0..and.ep3(i,j,k+1) == 0.)kflu=1
                 if(ep3(i,j,k) == 1..and.ep3(i,j,k+1) == 1.)ksol=1
                 do kraf=1,nraf
                    if(zepsi(i,j,kraf+nraf*(k-1)  ) == 0..and.&
                         zepsi(i,j,kraf+nraf*(k-1)+1) == 1.)kdebraf=kraf+nraf*(k-1)+1
                    if(zepsi(i,j,kraf+nraf*(k-1)  ) == 1..and.&
                         zepsi(i,j,kraf+nraf*(k-1)+1) == 0.)kfinraf=kraf+nraf*(k-1)+1
                 enddo
                 if(kdebraf /= 0.and.kfinraf /= 0.and.&
                      kdebraf < kfinraf.and.kflu == 1)then
                    kobj=kobj+1
                    do kk=kobj,nobjmax-1
                       zi(kk,i,j)=zi(kk+1,i,j)
                       zf(kk,i,j)=zf(kk+1,i,j)
                    enddo
                    kobj=kobj-1
                 endif
                 if(kdebraf /= 0.and.kfinraf /= 0.and.&
                      kdebraf > kfinraf.and.ksol == 1)then
                    kobj=kobj+1
                    do kk=kobj,nobjmax-1
                       zi(kk,i,j)=zi(kk+1,i,j)
                    enddo
                    kobj=kobj-1
                    do kk=kobj,nobjmax-1
                       zf(kk,i,j)=zf(kk+1,i,j)
                    enddo
                 endif
                 kdebraf=0
                 kfinraf=0
                 kflu=0
              enddo
           endif
        enddo
     enddo
  endif
!  if (nrank==0) print*,'    step 10'
  !
  return
end subroutine gene_epsi_3D
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine gene_epsim(epsi,epsim,nx,ny,nz,nobjx,nobjy,nobjz,&
     xi,xf,yi,yf,zi,zf,nobjmax,dx,dz,epm,&
     xlx,yly,zlz,yp)
  implicit none
  !
  real(4),dimension(       nx,ny,nz) :: epsi
  real(4),dimension(       nx,ny,nz) :: epsim
  real(4),dimension(  nobjmax,ny,nz) :: xi,xf
  real(4),dimension(  nobjmax,nx,nz) :: yi,yf
  real(4),dimension(  nobjmax,nx,ny) :: zi,zf
  integer,dimension(          ny,nz) :: nobjx
  integer,dimension(          nx,nz) :: nobjy
  integer,dimension(          nx,ny) :: nobjz
  real(4),dimension(             ny) :: yp
  real(4)                            :: x,y,z
  real(4)                            :: dx,dz
  real(4)                            :: xlx,yly,zlz
  real(4)                            :: xe,ye,ze,epm
  integer                            :: i,j,k
  integer                            :: nx,ny,nz
  integer                            :: ix,jy,kz
  integer                            :: nobjmax
  !
  epsim(:,:,:)=epsi(:,:,:)
  xe=epm*dx
  do k=1,nz
     do j=1,ny
        do i=1,nx
           if(epsi(i,j,k) == 1..and.&
                nobjx (j,k) /= 0)then
              x=dx*(i-1)
              do ix=1,nobjx(j,k)
                 if(x   >= xi(ix,j,k)   .and.&
                      x   <= xi(ix,j,k)+xe.and.&
                      0.  < xi(ix,j,k)   .or. &
                      x   <= xf(ix,j,k)   .and.&
                      x   >= xf(ix,j,k)-xe.and.&
                      xlx > xf(ix,j,k)   )then
                    epsim(i,j,k)=0.
                 endif
              enddo
           endif
        enddo
     enddo
  enddo
  do k=1,nz
     do i=1,nx
        do j=2,ny-1
           if(epsi(i,j,k) == 1..and.&
                nobjy (i,k) > 0)then
              y=yp(j)
              ye=epm*(yp(j+1)-yp(j-1))/2.
              do jy=1,nobjy(i,k)
                 if(y   >= yi(jy,i,k)   .and.&
                      y   <= yi(jy,i,k)+ye.and.&
                      0.  < yi(jy,i,k)   .or. &
                      y   <= yf(jy,i,k)   .and.&
                      y   >= yf(jy,i,k)-ye.and.&
                      yly > yf(jy,i,k)   )then
                    epsim(i,j,k)=0.
                 endif
              enddo
           endif
        enddo
     enddo
  enddo
  ze=epm*dz
  do j=1,ny
     do i=1,nx
        do k=1,nz
           if(epsi(i,j,k) == 1..and.&
                nobjz (i,j) > 0)then
              z=dz*(k-1)
              do kz=1,nobjz(i,j)
                 if(z   >= zi(kz,i,j)   .and.&
                      z   <= zi(kz,i,j)+ze.and.&
                      0.  < zi(kz,i,j)   .or. &
                      z   <= zf(kz,i,j)   .and.&
                      z   >= zf(kz,i,j)-ze.and.&
                      zlz > zf(kz,i,j)   )then
                    epsim(i,j,k)=0.
                 endif
              enddo
           endif
        enddo
     enddo
  enddo
  !
  return
end subroutine gene_epsim
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine verif_epsi(ep1,npif,izap,nx,ny,nz,nobjmax,&
     nxipif,nxfpif,nyipif,nyfpif,nzipif,nzfpif)
  use decomp_2d
  use MPI

    implicit none
    !
    integer                            :: nx,ny,nz,nobjmax
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ep2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ep3
    integer,dimension(0:nobjmax,xsize(2),xsize(3)) :: nxipif,nxfpif
    integer,dimension(0:nobjmax,ysize(1),ysize(3)) :: nyipif,nyfpif
    integer,dimension(0:nobjmax,zsize(1),zsize(2)) :: nzipif,nzfpif
    integer                            :: npif,izap
    integer                            :: i,j,k
    integer                            :: inum ,jnum ,knum
    integer                            :: iflu ,jflu ,kflu
    integer                            :: ising,jsing,ksing,itest
    integer                            :: mpi_aux_i, code

  !x-pencil
  nxipif(:,:,:)=npif
  nxfpif(:,:,:)=npif
  ising=0
  do k=1,xsize(3)
     do j=1,xsize(2)
        inum=0
        iflu=0
        if(ep1(1,j,k) == 1.)inum=inum+1
        if(ep1(1,j,k) == 0.)iflu=iflu+1
        do i=2,nx
           if(ep1(i  ,j,k) == 0.)iflu=iflu+1
           if(ep1(i-1,j,k) == 0..and.&
                ep1(i  ,j,k) == 1.)then
              inum=inum+1
              if(inum == 1)then
                 nxipif(inum  ,j,k)=iflu-izap
                 if(iflu-izap < npif)ising=ising+1
                 if(iflu-izap >= npif)nxipif(inum  ,j,k)=npif
                 iflu=0
              else
                 nxipif(inum  ,j,k)=iflu-izap
                 nxfpif(inum-1,j,k)=iflu-izap
                 if(iflu-izap < npif)ising=ising+1
                 if(iflu-izap >= npif)nxipif(inum  ,j,k)=npif
                 if(iflu-izap >= npif)nxfpif(inum-1,j,k)=npif
                 iflu=0
              endif
           endif
           if(ep1(i,j,k) == 1.)iflu=0
        enddo
        if(ep1(nx,j,k) == 0.)then
           nxfpif(inum,j,k)=iflu-izap
           if(iflu-izap < npif)ising=ising+1
           if(iflu-izap < npif)nxfpif(inum,j,k)=npif
        endif
     enddo
  enddo
  call MPI_REDUCE(ising,mpi_aux_i,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,code)
  if (code /= 0) call decomp_2d_abort(code, "MPI_REDUCE")
!  if (nrank==0) print*,'        number of points with potential problem in X :',mpi_aux_i
!  if (nrank==0) print*,'    step 11'

  !y-pencil
  call transpose_x_to_y(ep1,ep2)
  nyipif(:,:,:)=npif
  nyfpif(:,:,:)=npif
  jsing=0
  do k=1,ysize(3)
     do i=1,ysize(1)
        jnum=0
        jflu=0
        if(ep2(i,1,k) == 1.)jnum=jnum+1
        if(ep2(i,1,k) == 0.)jflu=jflu+1
        do j=2,ny
           if(ep2(i,j  ,k) == 0.)jflu=jflu+1
           if(ep2(i,j-1,k) == 0..and.&
                ep2(i,j  ,k) == 1.)then
              jnum=jnum+1
              if(jnum == 1)then
                 nyipif(jnum  ,i,k)=jflu-izap
                 if(jflu-izap < npif)jsing=jsing+1
                 if(jflu-izap >= npif)nyipif(jnum  ,i,k)=npif
                 jflu=0
              else
                 nyipif(jnum  ,i,k)=jflu-izap
                 nyfpif(jnum-1,i,k)=jflu-izap
                 if(jflu-izap < npif)jsing=jsing+1
                 if(jflu-izap >= npif)nyipif(jnum  ,i,k)=npif
                 if(jflu-izap >= npif)nyfpif(jnum-1,i,k)=npif
                 jflu=0
              endif
           endif
           if(ep2(i,j,k) == 1.)jflu=0
        enddo
        if(ep2(i,ny,k) == 0.)then
           nyfpif(jnum,i,k)=jflu-izap
           if(jflu-izap < npif)jsing=jsing+1
           if(jflu-izap < npif)nyfpif(jnum,i,k)=npif
        endif
     enddo
  enddo
  call MPI_REDUCE(jsing,mpi_aux_i,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,code)
  if (code /= 0) call decomp_2d_abort(code, "MPI_REDUCE")
!  if (nrank==0) print*,'        number of points with potential problem in Y :',mpi_aux_i
!  if (nrank==0) print*,'    step 12'

  !z-pencil
  if(nz > 1)then
     call transpose_y_to_z(ep2,ep3)
     nzipif(:,:,:)=npif
     nzfpif(:,:,:)=npif
     ksing=0
     do j=1,zsize(2)
        do i=1,zsize(1)
           knum=0
           kflu=0
           if(ep3(i,j,1) == 1.)knum=knum+1
           if(ep3(i,j,1) == 0.)kflu=kflu+1
           do k=2,nz
              if(ep3(i,j,k  ) == 0.)kflu=kflu+1
              if(ep3(i,j,k-1) == 0..and.&
                   ep3(i,j,k  ) == 1.)then
                 knum=knum+1
                 if(knum == 1)then
                    nzipif(knum  ,i,j)=kflu-izap
                    if(kflu-izap < npif)ksing=ksing+1
                    if(kflu-izap >= npif)nzipif(knum  ,i,j)=npif
                    kflu=0
                 else
                    nzipif(knum  ,i,j)=kflu-izap
                    nzfpif(knum-1,i,j)=kflu-izap
                    if(kflu-izap < npif)ksing=ksing+1
                    if(kflu-izap >= npif)nzipif(knum  ,i,j)=npif
                    if(kflu-izap >= npif)nzfpif(knum-1,i,j)=npif
                    kflu=0
                 endif
              endif
              if(ep3(i,j,k) == 1.)kflu=0
           enddo
           if(ep3(i,j,nz) == 0.)then
              nzfpif(knum,i,j)=kflu-izap
              if(kflu-izap < npif)ksing=ksing+1
              if(kflu-izap < npif)nzfpif(knum,i,j)=npif
           endif
        enddo
     enddo
     call MPI_REDUCE(ksing,mpi_aux_i,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,code)
     if (code /= 0) call decomp_2d_abort(code, "MPI_REDUCE")
!     if (nrank==0) print*,'        number of points with potential problem in Z :',mpi_aux_i
  endif
!  if (nrank==0) print*,'    step 13'
  !
  return
end subroutine verif_epsi
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
  subroutine write_geomcomplex(nx,ny,nz,ep1,nobjx,nobjy,nobjz,xi,xf,yi,yf,zi,zf,&
       nxipif,nxfpif,nyipif,nyfpif,nzipif,nzfpif,nobjmax,npif)
    use decomp_2d
    USE decomp_2d_io
    implicit none
    !
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    integer                            :: nx,ny,nz,nobjmax
    integer,dimension(xstart(2):xend(2),xstart(3):xend(3)) :: nobjx
    integer,dimension(ystart(1):yend(1),ystart(3):yend(3)) :: nobjy
    integer,dimension(zstart(1):zend(1),zstart(2):zend(2)) :: nobjz
    real(mytype),dimension(nobjmax,xstart(2):xend(2),xstart(3):xend(3)) :: xi,xf
    real(mytype),dimension(nobjmax,ystart(1):yend(1),ystart(3):yend(3)) :: yi,yf
    real(mytype),dimension(nobjmax,zstart(1):zend(1),zstart(2):zend(2)) :: zi,zf
    integer,dimension(0:nobjmax,xstart(2):xend(2),xstart(3):xend(3)) :: nxipif,nxfpif
    integer,dimension(0:nobjmax,ystart(1):yend(1),ystart(3):yend(3)) :: nyipif,nyfpif
    integer,dimension(0:nobjmax,zstart(1):zend(1),zstart(2):zend(2)) :: nzipif,nzfpif
    integer                            :: npif
    integer                            :: i,j,k,count
    !###################################################################
    if (nrank==0) print *,'Writing geometry'
    call decomp_2d_write_one(1,ep1,'geometry/epsilon.dat')
    !###################################################################
    !x-pencil
    open(67,file='geometry/nobjx.dat',form='formatted',access='direct',recl=13)
    do k=xstart(3),xend(3)
       do j=xstart(2),xend(2)
          count = (k-1)*ny+j
          write(67,'(1I12,A)',rec=count) nobjx(j,k),char(10)
       enddo
    enddo
    close(67)
    !y-pencil
    open(67,file='geometry/nobjy.dat',form='formatted',access='direct',recl=13)
    do k=ystart(3),yend(3)
       do i=ystart(1),yend(1)
          count = (k-1)*nx+i
          write(67,'(1I12,A)',rec=count) nobjy(i,k),char(10)
       enddo
    enddo
    close(67)
    !z-pencil
    open(67,file='geometry/nobjz.dat',form='formatted',access='direct',recl=13)
    do j=zstart(2),zend(2)
       do i=zstart(1),zend(1)
          count = (j-1)*nx+i
          write(67,'(1I12,A)',rec=count) nobjz(i,j),char(10)
       enddo
    enddo
    close(67)
    !###################################################################
    !x-pencil
    open(67,file='geometry/nxifpif.dat',form='formatted',access='direct',recl=25)
    do k=xstart(3),xend(3)
       do j=xstart(2),xend(2)
          do i=0,nobjmax
             count = (k-1)*ny*(1+nobjmax)+(j-1)*(1+nobjmax)+i+1
             write(67,'(2I12,A)',rec=count) nxipif(i,j,k),nxfpif(i,j,k),char(10)
          enddo
       enddo
    enddo
    close(67)
    !y-pencil
    open(67,file='geometry/nyifpif.dat',form='formatted',access='direct',recl=25)
    do k=ystart(3),yend(3)
       do i=ystart(1),yend(1)
          do j=0,nobjmax
             count = (k-1)*nx*(1+nobjmax)+(i-1)*(1+nobjmax)+j+1
             write(67,'(2I12,A)',rec=count) nyipif(j,i,k),nyfpif(j,i,k),char(10)
          enddo
       enddo
    enddo
    close(67)
    !z-pencil
    open(67,file='geometry/nzifpif.dat',form='formatted',access='direct',recl=25)
    do j=zstart(2),zend(2)
       do i=zstart(1),zend(1)
          do k=0,nobjmax
             count = (j-1)*nx*(1+nobjmax)+(i-1)*(1+nobjmax)+k+1
             write(67,'(2I12,A)',rec=count) nzipif(k,i,j),nzfpif(k,i,j),char(10)
          enddo
       enddo
    enddo
    close(67)
    !###################################################################
    !x-pencil
    open(67,file='geometry/xixf.dat',form='formatted',access='direct',recl=29)
    do k=xstart(3),xend(3)
       do j=xstart(2),xend(2)
          do i=1,nobjmax
             count = (k-1)*ny*nobjmax+(j-1)*nobjmax+i
             write(67,'(2E14.6,A)',rec=count) xi(i,j,k),xf(i,j,k),char(10)
          enddo
       enddo
    enddo
    close(67)
    !y-pencil
    open(67,file='geometry/yiyf.dat',form='formatted',access='direct',recl=29)
    do k=ystart(3),yend(3)
       do i=ystart(1),yend(1)
          do j=1,nobjmax
             count = (k-1)*nx*nobjmax+(i-1)*nobjmax+j
             write(67,'(2E14.6,A)',rec=count) yi(j,i,k),yf(j,i,k),char(10)
          enddo
       enddo
    enddo
    close(67)
    !z-pencil
    open(67,file='geometry/zizf.dat',form='formatted',access='direct',recl=29)
    do j=zstart(2),zend(2)
       do i=zstart(1),zend(1)
          do k=1,nobjmax
             count = (j-1)*nx*nobjmax+(i-1)*nobjmax+k
             write(67,'(2E14.6,A)',rec=count) zi(k,i,j),zf(k,i,j),char(10)
          enddo
       enddo
    enddo
    close(67)
    !###################################################################
    return
  end subroutine write_geomcomplex
  !############################################################################
  !############################################################################
  subroutine read_geomcomplex()
    !
    USE complex_geometry
    USE decomp_2d
    USE MPI
    !
    implicit none
    !
    integer :: i,j,k
    integer :: code
    !
    if(nrank == 0)then
       open(11,file='nobjx.dat'  ,form='formatted', status='old')
       do k=1,nz
          do j=1,ny
             read(11,*)nobjx(j,k)
          enddo
       enddo
       close(11)
    endif
    call MPI_BCAST(nobjx,ny*nz,MPI_INTEGER,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_BCAST")
    if(nrank == 0)then
       open(12,file='nobjy.dat'  ,form='formatted', status='old')
       do k=1,nz
          do i=1,nx
             read(12,*)nobjy(i,k)
          enddo
       enddo
       close(12)
    endif
    call MPI_BCAST(nobjy,nx*nz,MPI_INTEGER,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_BCAST")
    if(nrank == 0)then
       open(13,file='nobjz.dat'  ,form='formatted', status='old')
       do j=1,ny
          do i=1,nx
             read(13,*)nobjz(i,j)
          enddo
       enddo
       close(13)
    endif
    call MPI_BCAST(nobjz,nx*ny,MPI_INTEGER,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_BCAST")
    if(nrank == 0)then
       open(21,file='nxifpif.dat',form='formatted', status='old')
       do k=1,nz
          do j=1,ny
             do i=0,nobjmax
                read(21,*)nxipif(i,j,k),nxfpif(i,j,k)
             enddo
          enddo
       enddo
       close(21)
    endif
    call MPI_BCAST(nxipif,ny*nz*(nobjmax+1),MPI_INTEGER,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_BCAST")
    call MPI_BCAST(nxfpif,ny*nz*(nobjmax+1),MPI_INTEGER,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_BCAST")
    if(nrank == 0)then
       open(22,file='nyifpif.dat',form='formatted', status='old')
       do k=1,nz
          do i=1,nx
             do j=0,nobjmax
                read(22,*)nyipif(j,i,k),nyfpif(j,i,k)
             enddo
          enddo
       enddo
       close(22)
    endif
    call MPI_BCAST(nyipif,nx*nz*(nobjmax+1),MPI_INTEGER,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_BCAST")
    call MPI_BCAST(nyfpif,nx*nz*(nobjmax+1),MPI_INTEGER,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_BCAST")
    if(nrank == 0)then
       open(23,file='nzifpif.dat',form='formatted', status='old')
       do j=1,ny
          do i=1,nx
             do k=0,nobjmax
                read(23,*)nzipif(k,i,j),nzfpif(k,i,j)
             enddo
          enddo
       enddo
       close(23)
    endif
    call MPI_BCAST(nzipif,nx*ny*(nobjmax+1),MPI_INTEGER,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_BCAST")
    call MPI_BCAST(nzfpif,nx*ny*(nobjmax+1),MPI_INTEGER,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_BCAST")
    if(nrank == 0)then
       open(31,file='xixf.dat'   ,form='formatted', status='old')
       do k=1,nz
          do j=1,ny
             do i=1,nobjmax
                read(31,*)xi(i,j,k),xf(i,j,k)
             enddo
          enddo
       enddo
       close(31)
    endif
    call MPI_BCAST(xi,ny*nz*nobjmax,MPI_REAL,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_BCAST")
    call MPI_BCAST(xf,ny*nz*nobjmax,MPI_REAL,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_BCAST")
    if(nrank == 0)then
       open(32,file='yiyf.dat'   ,form='formatted', status='old')
       do k=1,nz
          do i=1,nx
             do j=1,nobjmax
                read(32,*)yi(j,i,k),yf(j,i,k)
             enddo
          enddo
       enddo
       close(32)
    endif
    call MPI_BCAST(yi,nx*nz*nobjmax,MPI_REAL,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_BCAST")
    call MPI_BCAST(yf,nx*nz*nobjmax,MPI_REAL,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_BCAST")
    if(nrank == 0)then
       open(33,file='zizf.dat'   ,form='formatted', status='old')
       do j=1,ny
          do i=1,nx
             do k=1,nobjmax
                read(33,*)zi(k,i,j),zf(k,i,j)
             enddo
          enddo
       enddo
       close(33)
    endif
    call MPI_BCAST(zi,nx*ny*nobjmax,MPI_REAL,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_BCAST")
    call MPI_BCAST(zf,nx*ny*nobjmax,MPI_REAL,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_BCAST")
    !
    return
  end subroutine read_geomcomplex
!############################################################################
!############################################################################
end module genepsi
