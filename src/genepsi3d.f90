!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module genepsi

  use decomp_2d_constants
  use decomp_2d_mpi
  use mod_stret, only : stretching

  public

contains

  subroutine epsi_init(ep1)

    USE param, only : zero, one, dx, dz
    USE decomp_2d, only : xstart, xend, xsize
    !USE decomp_2d_io
    USE variables, only : yp, ny

    implicit none
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(inout) :: ep1
    logical :: dir_exists

#ifdef DEBG
    if (nrank .eq. 0) write(*,*)'# body_init start'
#endif
    !###################################################################
    ! Check if geometry folder exists
    !###################################################################
    if (nrank==0) then
      inquire(file="geometry", exist=dir_exists)
      if (.not.dir_exists) then
        call execute_command_line("mkdir geometry 2> /dev/null")
      end if
    end if
    !###################################################################
    ep1(:,:,:)=zero
    call geomcomplex(ep1,xstart(1),xend(1),ny,xstart(2),xend(2),xstart(3),xend(3),dx,yp,dz,one)

#ifdef DEBG
    if (nrank .eq. 0) write(*,*)'# body_init done'
#endif

    return
  end subroutine epsi_init
!############################################################################
!############################################################################
  subroutine geomcomplex(epsi, nxi, nxf, ny, nyi, nyf, nzi, nzf, dx, yp, dz, remp)

    USE param, ONLY : itype, itype_cyl, itype_hill, itype_channel,&
                      itype_sandbox, itype_pipe
    USE cyl, ONLY : geomcomplex_cyl
    USE hill, ONLY : geomcomplex_hill
    USE channel, ONLY : geomcomplex_channel
    USE sandbox, ONLY : geomcomplex_sandbox
    USE pipe, ONLY : geomcomplex_pipe

    IMPLICIT NONE

    INTEGER :: nxi,nxf,ny,nyi,nyf,nzi,nzf
    REAL(mytype),DIMENSION(nxi:nxf,nyi:nyf,nzi:nzf),intent(inout) :: epsi
    REAL(mytype)               :: dx,dz
    REAL(mytype),DIMENSION(ny) :: yp
    REAL(mytype)               :: remp

    IF (itype.EQ.itype_cyl) THEN

       CALL geomcomplex_cyl(epsi, nxi, nxf, ny, nyi, nyf, nzi, nzf, dx, yp, remp)

    ELSEIF (itype.EQ.itype_hill) THEN

       CALL geomcomplex_hill(epsi,nxi,nxf,ny,nyi,nyf,nzi,nzf,dx,yp,dz,remp)

    ELSEIF (itype.EQ.itype_channel) THEN

       CALL geomcomplex_channel(epsi, nxi, nxf, ny, nyi, nyf, nzi, nzf, yp, remp)

    ELSEIF (itype.EQ.itype_sandbox) THEN
     
       CALL  geomcomplex_sandbox(epsi, nxi, nxf, ny, nyi, nyf, nzi, nzf, yp, remp)

    ELSEIF (itype.EQ.itype_pipe) THEN

       CALL geomcomplex_pipe(epsi, nxi, nxf, ny, nyi, nyf, nzi, nzf, dx, yp, dz, remp)

    ENDIF

  end subroutine geomcomplex
!############################################################################
!############################################################################
  subroutine genepsi3d(ep1)

    USE variables, only : nx,ny,nz,nxm,nym,nzm,yp, ilist
    USE param, only : xlx,yly,zlz,dx,dy,dz,izap,npif,nclx,ncly,nclz,istret,itype,itype_sandbox
    use param, only : itime
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
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(inout) :: ep1
    !
    if (nrank==0.and.mod(itime,ilist)==0) then
      write(*,*)'==========================================================='
      write(*,*)'Generating the geometry!'
    end if
    !###################################################################
    ! Check if planes folder exists
    !###################################################################
    if (nrank==0) then
      inquire(file="data", exist=dir_exists)
      if (.not.dir_exists) then
        call execute_command_line("mkdir data 2> /dev/null")
      end if
      inquire(file="data/geometry", exist=dir_exists)
      if (.not.dir_exists) then
        call execute_command_line("mkdir data/geometry 2> /dev/null")
      end if
    end if
    !###################################################################
    
    if (itype.eq.itype_sandbox) then !F.Schuch 2020-08-14T11:44:11-03:00
      call geomcomplex_io(nx,ny,nz,ep1,nobjx,nobjy,nobjz,xi,xf,yi,yf,zi,zf,&
           nxipif,nxfpif,nyipif,nyfpif,nzipif,nzfpif,nobjmax,npif,.true.)
    else
      call gene_epsi_3D(ep1,nx,ny,nz,dx,dy,dz,xlx,yly,zlz ,&
           nclx,ncly,nclz,nxraf,nyraf,nzraf   ,&
           xi,xf,yi,yf,zi,zf,nobjx,nobjy,nobjz,&
           nobjmax,yp,nraf)
      call verif_epsi(ep1,npif,izap,nx,ny,nz,nobjmax,&
           nxipif,nxfpif,nyipif,nyfpif,nzipif,nzfpif)
    endif
    ! call geomcomplex_io(nx,ny,nz,ep1,nobjx,nobjy,nobjz,xi,xf,yi,yf,zi,zf,&
    !      nxipif,nxfpif,nyipif,nyfpif,nzipif,nzfpif,nobjmax,npif,.false.)
    !
  end subroutine genepsi3d
!
!############################################################################
!############################################################################
  subroutine gene_epsi_3D(ep1,nx,ny,nz,dx,dy,dz,xlx,yly,zlz ,&
       nclx,ncly,nclz,nxraf,nyraf,nzraf   ,&
       xi,xf,yi,yf,zi,zf,nobjx,nobjy,nobjz,&
       nobjmax,yp,nraf)
    use param, only : istret, zero, half, one, two
    use var, only : ta2, ta3
    use decomp_2d
    use MPI
    use complex_geometry, only: xepsi, yepsi, zepsi
    implicit none
    !
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(inout):: ep1
    !!real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ep1,smoofun,fbcx,fbcy,fbcz
    !!real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ep2
    !!real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ep3
    integer, intent(in)                                  :: nx,ny,nz,nobjmax
    real(mytype), intent(in)                             :: dx,dy,dz
    real(mytype), intent(in)                             :: xlx,yly,zlz
    logical, intent(in)                                  :: nclx,ncly,nclz
    integer, intent(in)                                  :: nxraf,nyraf,nzraf
    real(mytype),dimension(nobjmax,xsize(2),xsize(3)):: xi,xf
    real(mytype),dimension(nobjmax,ysize(1),ysize(3)):: yi,yf
    real(mytype),dimension(nobjmax,zsize(1),zsize(2)):: zi,zf
    integer                                          :: nraf
    integer,     dimension(xsize(2),xsize(3))        :: nobjx
    integer,     dimension(ysize(1),ysize(3))        :: nobjy
    integer,     dimension(zsize(1),zsize(2))        :: nobjz
    real(mytype),dimension(ny)                       :: yp
    
    ! Local variables
    integer,     dimension(xsize(2),xsize(3))          :: nobjxraf
    integer,     dimension(ysize(1),ysize(3))          :: nobjyraf
    integer,     dimension(zsize(1),zsize(2))          :: nobjzraf
    !real(mytype),dimension(nxraf,xsize(2),xsize(3))    :: xepsi
    !real(mytype),dimension(ysize(1),nyraf,ysize(3))    :: yepsi
    !real(mytype),dimension(zsize(1),zsize(2),nzraf)    :: zepsi
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
    integer                          :: nobjxmax,nobjymax ,nobjzmax
    integer                          :: nobjxmaxraf,nobjymaxraf,nobjzmaxraf
    integer                          :: idebraf,jdebraf,kdebraf
    integer                          :: ifinraf,jfinraf,kfinraf
    character(len=4) suffixe
    integer                          :: numvis
    integer                          :: mpi_aux_i, code

    !x-pencil
    ep1=zero
    call geomcomplex(ep1,xstart(1),xend(1),ny,xstart(2),xend(2),xstart(3),xend(3),dx,yp,dz,one)
    ! if (nrank==0) print*,'    step 1'
    if(nclx)then
       dxraf =xlx/real(nxraf, mytype)
    else
       dxraf =xlx/real(nxraf-1, mytype)
    endif
    xepsi=zero
    call geomcomplex(xepsi,1,nxraf,ny,xstart(2),xend(2),xstart(3),xend(3),dxraf,yp,dz,one)
    ! if (nrank==0) print*,'    step 2'
    !y-pencil
    if(ncly)then
       dyraf =yly/real(nyraf, mytype)
    else
       dyraf =yly/real(nyraf-1, mytype)
    endif

    ! Compute ypraf
    if (istret.eq.0) then
       do j = 1, nyraf
          ypraf(j) = (j-1) * dyraf
       end do
    else
       call stretching(nyraf, ypraf, opt_write = .false.)
    end if

    yepsi=zero
    call geomcomplex(yepsi,ystart(1),yend(1),nyraf,1,nyraf,ystart(3),yend(3),dx,ypraf,dz,one)
    ! if (nrank==0) print*,'    step 3'
    !z-pencil
    if(nclz)then
       dzraf=zlz/real(nzraf, mytype)
    else
       dzraf=zlz/real(nzraf-1, mytype)
    endif
    zepsi=zero
    call geomcomplex(zepsi,zstart(1),zend(1),ny,zstart(2),zend(2),1,nzraf,dx,yp,dzraf,one)
    ! if (nrank==0) print*,'    step 4'

    !x-pencil
    nobjx(:,:)=0
    nobjxmax=0
    do k=1,xsize(3)
       do j=1,xsize(2)
          inum=0
          if(ep1(1,j,k).eq.1.)then
             inum=1
             nobjx(j,k)=1
          endif
          do i=1,nx-1
             if(ep1(i,j,k).eq.0..and.ep1(i+1,j,k).eq.1.)then
                inum=inum+1
                nobjx(j,k)=nobjx(j,k)+1
             endif
          enddo
          if(inum.gt.nobjxmax)then
             nobjxmax=inum
          endif
       enddo
    enddo
    call MPI_REDUCE(nobjxmax,mpi_aux_i,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,code)
    ! if (nrank==0) print*,'        nobjxmax=',mpi_aux_i

    nobjxraf(:,:)=0
    ibug=0
    nobjxmaxraf=0
    inum=0
    do k=1,xsize(3)
       do j=1,xsize(2)
          inum=0
          if(xepsi(1,j,k).eq.1.)then
             inum=1
             nobjxraf(j,k)=1
          endif
          do i=1,nxraf-1
             if(xepsi(i,j,k).eq.zero.and.xepsi(i+1,j,k).eq.one)then
                inum=inum+1
                nobjxraf(j,k)=nobjxraf(j,k)+1
             endif
          enddo
          if(inum.gt.nobjxmaxraf)then
             nobjxmaxraf=inum
          endif
          if(nobjx(j,k).ne.nobjxraf(j,k))then
             ibug=ibug+1
          endif
       enddo
    enddo
    call MPI_REDUCE(nobjxmaxraf,mpi_aux_i,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,code)
    ! if (nrank==0) print*,'        nobjxmaxraf=',mpi_aux_i
    call MPI_REDUCE(ibug,mpi_aux_i,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,code)
    ! if (nrank==0) print*,'        ibug=',mpi_aux_i
    ! if (nrank==0) print*,'    step 5'

    !y-pencil
    nobjy(:,:)=0
    nobjymax=0
    call transpose_x_to_y(ep1,ta2)
    do k=1,ysize(3)
       do i=1,ysize(1)
          jnum=0
          if(ta2(i,1,k) == one)then
             jnum=1
             nobjy(i,k)=1
          endif
          do j=1,ny-1
             if(ta2(i,j,k) == zero .and. ta2(i,j+1,k) == one)then
                jnum=jnum+1
                nobjy(i,k)=nobjy(i,k)+1
             endif
          enddo
          if(jnum.gt.nobjymax)then
             nobjymax=jnum
          endif
       enddo
    enddo
    call MPI_REDUCE(nobjymax,mpi_aux_i,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,code)
    ! if (nrank==0) print*,'        nobjymax=',mpi_aux_i

    nobjyraf(:,:)=0
    jbug=0
    nobjymaxraf=0
    jnum=0
    do k=1,ysize(3)
       do i=1,ysize(1)
          jnum=0
          if(yepsi(i,1,k) == one)then
             jnum=1
             nobjyraf(i,k)=1
          endif
          do j=1,nyraf-1
             if(yepsi(i,j,k) == zero .and. yepsi(i,j+1,k) == one)then
                jnum=jnum+1
                nobjyraf(i,k)=nobjyraf(i,k)+1
             endif
          enddo
          if(jnum.gt.nobjymaxraf)then
             nobjymaxraf=jnum
          endif
          if(nobjy(i,k).ne.nobjyraf(i,k))then
             jbug=jbug+1
          endif
       enddo
    enddo
    call MPI_REDUCE(nobjymaxraf,mpi_aux_i,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,code)
    ! if (nrank==0) print*,'        nobjymaxraf=',mpi_aux_i
    call MPI_REDUCE(jbug,mpi_aux_i,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,code)
    ! if (nrank==0) print*,'        jbug=',mpi_aux_i
    ! if (nrank==0) print*,'    step 6'

    !z-pencil
    nobjz(:,:)=0
    nobjzmax=0
    call transpose_y_to_z(ta2,ta3)
    do j=1,zsize(2)
       do i=1,zsize(1)
          knum=0
          if(ta3(i,j,1) == one)then
             knum=1
             nobjz(i,j)=1
          endif
          do k=1,nz-1
             if(ta3(i,j,k) == zero .and. ta3(i,j,k+1) == one)then
                knum=knum+1
                nobjz(i,j)=nobjz(i,j)+1
             endif
          enddo
          if(knum.gt.nobjzmax)then
             nobjzmax=knum
          endif
       enddo
    enddo
    call MPI_REDUCE(nobjzmax,mpi_aux_i,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,code)
    ! if (nrank==0) print*,'        nobjzmax=',mpi_aux_i

    nobjzraf(:,:)=0
    kbug=0
    nobjzmaxraf=0
    knum=0
    do j=1,zsize(2)
       do i=1,zsize(1)
          knum=0
          if(zepsi(i,j,1) == one)then
             knum=1
             nobjzraf(i,j)=1
          endif
          do k=1,nzraf-1
             if(zepsi(i,j,k) == zero .and. zepsi(i,j,k+1) == one)then
                knum=knum+1
                nobjzraf(i,j)=nobjzraf(i,j)+1
             endif
          enddo
          if(knum.gt.nobjzmaxraf)then
             nobjzmaxraf=knum
          endif
          if(nobjz(i,j).ne.nobjzraf(i,j))then
             kbug=kbug+1
          endif
       enddo
    enddo
    call MPI_REDUCE(nobjzmaxraf,mpi_aux_i,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,code)
    ! if (nrank==0) print*,'        nobjzmaxraf=',mpi_aux_i
    call MPI_REDUCE(kbug,mpi_aux_i,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,code)
    ! if (nrank==0) print*,'        kbug=',mpi_aux_i
    ! if (nrank==0) print*,'    step 7'

    !x-pencil
    do k=1,xsize(3)
       do j=1,xsize(2)
          inum=0
          if(xepsi(1,j,k) == one)then
             inum=inum+1
             xi(inum,j,k)=-dx!-xlx
          endif
          do i=1,nxraf-1
             if(xepsi(i,j,k) == zero .and. xepsi(i+1,j,k) == one)then
                inum=inum+1
                xi(inum,j,k)=dxraf*(i-1)+dxraf/2.
             elseif(xepsi(i,j,k) == one .and. xepsi(i+1,j,k)== zero)then
                xf(inum,j,k)=dxraf*(i-1)+dxraf/2.
             endif
          enddo
          if(xepsi(nxraf,j,k)==1.)then
             xf(inum,j,k)=xlx+dx!2.*xlx
          endif
       enddo
    enddo

    if(ibug /= 0)then
       do k=1,xsize(3)
          do j=1,xsize(2)
             if(nobjx(j,k) /= nobjxraf(j,k))then
                iobj=0
                if(ep1(1,j,k) == one)iobj=iobj+1
                do i=1,nx-1
                   if(ep1(i,j,k) == zero .and. ep1(i+1,j,k) ==  one)iobj=iobj+1
                   if(ep1(i,j,k) == zero .and. ep1(i+1,j,k) == zero)iflu=1
                   if(ep1(i,j,k) ==  one .and. ep1(i+1,j,k) ==  one)isol=1
                   do iraf=1,nraf
                      if(xepsi(iraf+nraf*(i-1)  ,j,k) == zero .and.&
                         xepsi(iraf+nraf*(i-1)+1,j,k) ==  one)idebraf=iraf+nraf*(i-1)+1
                      if(xepsi(iraf+nraf*(i-1)  ,j,k) ==  one .and.&
                         xepsi(iraf+nraf*(i-1)+1,j,k) == zero)ifinraf=iraf+nraf*(i-1)+1
                   enddo
                   if(idebraf /= 0 .and. ifinraf /= 0 .and.&
                      idebraf < ifinraf .and. iflu == 1)then
                      iobj=iobj+1
                      do ii=iobj,nobjmax-1
                         xi(ii,j,k)=xi(ii+1,j,k)
                         xf(ii,j,k)=xf(ii+1,j,k)
                      enddo
                      iobj=iobj-1
                   endif
                   if(idebraf /= 0 .and. ifinraf /= 0 .and.&
                      idebraf > ifinraf .and. isol==1)then
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
    !if (nrank==0) write(*,*) '    step 8'

    !y-pencil
    do k=1,ysize(3)
       do i=1,ysize(1)
          jnum=0
          if(yepsi(i,1,k) == one)then
             jnum=jnum+1
             yi(jnum,i,k)=-(yp(2)-yp(1))!-yly
          endif
          do j=1,nyraf-1
             if(yepsi(i,j,k) == zero .and. yepsi(i,j+1,k) == one)then
                jnum=jnum+1
                yi(jnum,i,k)=ypraf(j)+(ypraf(j+1)-ypraf(j))*half!dyraf*(j-1)+dyraf/2.
             elseif(yepsi(i,j,k) == one .and. yepsi(i,j+1,k) == zero)then
                yf(jnum,i,k)=ypraf(j)+(ypraf(j+1)-ypraf(j))*half!dyraf*(j-1)+dyraf/2.
             endif
          enddo
          if(yepsi(i,nyraf,k) == one)then
             yf(jnum,i,k)=yly+(yp(ny)-yp(ny-1))*half!2.*yly
          endif
       enddo
    enddo

    if(jbug /= 0)then
       do k=1,ysize(3)
          do i=1,ysize(1)
             if(nobjy(i,k) /= nobjyraf(i,k))then
                jobj=0
                if(ta2(i,1,k) == one)jobj=jobj+1
                do j=1,ny-1
                   if(ta2(i,j,k) == zero .and. ta2(i,j+1,k) ==  one)jobj=jobj+1
                   if(ta2(i,j,k) == zero .and. ta2(i,j+1,k) == zero)jflu=1
                   if(ta2(i,j,k) ==  one .and. ta2(i,j+1,k) ==  one)jsol=1
                   do jraf=1,nraf
                      if(yepsi(i,jraf+nraf*(j-1)  ,k) == zero .and.&
                         yepsi(i,jraf+nraf*(j-1)+1,k) ==  one)jdebraf=jraf+nraf*(j-1)+1
                      if(yepsi(i,jraf+nraf*(j-1)  ,k) ==  one .and.&
                         yepsi(i,jraf+nraf*(j-1)+1,k) == zero)jfinraf=jraf+nraf*(j-1)+1
                   enddo
                   if(jdebraf /= 0 .and. jfinraf /= 0 .and.&
                      jdebraf < jfinraf.and.jflu == 1)then
                      jobj=jobj+1
                      do jj=jobj,nobjmax-1
                         yi(jj,i,k)=yi(jj+1,i,k)
                         yf(jj,i,k)=yf(jj+1,i,k)
                      enddo
                      jobj=jobj-1
                   endif
                   if(jdebraf /= 0 .and. jfinraf /= 0 .and.&
                      jdebraf > jfinraf .and. jsol == 1)then
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
    !if (nrank==0) write(*,*) '    step 9'

    !z-pencil
    do j=1,zsize(2)
       do i=1,zsize(1)
          knum=0
          if(zepsi(i,j,1) == one)then
             knum=knum+1
             zi(knum,i,j)=-dz!zlz
          endif
          do k=1,nzraf-1
             if(zepsi(i,j,k) == zero .and. zepsi(i,j,k+1) == one)then
                knum=knum+1
                zi(knum,i,j)=dzraf*(k-1)+dzraf*half
             elseif(zepsi(i,j,k) == one .and. zepsi(i,j,k+1) == zero)then
                zf(knum,i,j)=dzraf*(k-1)+dzraf*half
             endif
          enddo
          if(zepsi(i,j,nzraf) == one)then
             zf(knum,i,j)=zlz+dz!2.*zlz
          endif
       enddo
    enddo

    kdebraf=0
    if(kbug.ne.0)then
       do j=1,zsize(2)
          do i=1,zsize(1)
             if(nobjz(i,j) /= nobjzraf(i,j))then
                kobj=0
                if(ta3(i,j,1) == one)kobj=kobj+1
                do k=1,nz-1
                   if(ta3(i,j,k) == zero .and. ta3(i,j,k+1) ==  one)kobj=kobj+1
                   if(ta3(i,j,k) == zero .and. ta3(i,j,k+1) == zero)kflu=1
                   if(ta3(i,j,k) ==  one .and. ta3(i,j,k+1) ==  one)ksol=1
                   do kraf=1,nraf
                      if(zepsi(i,j,kraf+nraf*(k-1)  ) == zero .and.&
                         zepsi(i,j,kraf+nraf*(k-1)+1) ==  one)kdebraf=kraf+nraf*(k-1)+1
                      if(zepsi(i,j,kraf+nraf*(k-1)  ) ==  one .and.&
                         zepsi(i,j,kraf+nraf*(k-1)+1) == zero)kfinraf=kraf+nraf*(k-1)+1
                   enddo
                   if(kdebraf /= 0      .and. kfinraf /= 0 .and.&
                      kdebraf < kfinraf .and. kflu    == 1)then
                      kobj=kobj+1
                      do kk=kobj,nobjmax-1
                         zi(kk,i,j)=zi(kk+1,i,j)
                         zf(kk,i,j)=zf(kk+1,i,j)
                      enddo
                      kobj=kobj-1
                   endif
                   if(kdebraf /= 0      .and. kfinraf /= 0.and.&
                      kdebraf > kfinraf .and.    ksol == 1)then
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
    ! if (nrank==0) print*,'    step 10'
    !
    return
  end subroutine gene_epsi_3D
!############################################################################
!############################################################################
  subroutine verif_epsi(ep1,npif,izap,nx,ny,nz,nobjmax,&
       nxipif,nxfpif,nyipif,nyfpif,nzipif,nzfpif)
    use decomp_2d
    use param, only: zero, one
    use var, only : ta2, ta3
    use MPI

    implicit none
    !
    integer                            :: nx,ny,nz,nobjmax
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    !real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ep2
    !real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ep3
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
          if(ep1(1,j,k) ==  one)inum=inum+1
          if(ep1(1,j,k) == zero)iflu=iflu+1
          do i=2,nx
             if(ep1(i  ,j,k) == zero)iflu=iflu+1
             if(ep1(i-1,j,k) == zero .and.&
                ep1(i  ,j,k) ==  one)then
                inum=inum+1
                if(inum == 1)then
                   nxipif(inum  ,j,k)=iflu-izap
                   if(iflu-izap <  npif)ising=ising+1
                   if(iflu-izap >= npif)nxipif(inum  ,j,k)=npif
                   iflu=0
                else
                   nxipif(inum  ,j,k)=iflu-izap
                   nxfpif(inum-1,j,k)=iflu-izap
                   if(iflu-izap <  npif)ising=ising+1
                   if(iflu-izap >= npif)nxipif(inum  ,j,k)=npif
                   if(iflu-izap >= npif)nxfpif(inum-1,j,k)=npif
                   iflu=0
                endif
             endif
             if(ep1(i,j,k) == one)iflu=0
          enddo
          if(ep1(nx,j,k) == zero)then
             nxfpif(inum,j,k)=iflu-izap
             if(iflu-izap < npif)ising=ising+1
             if(iflu-izap < npif)nxfpif(inum,j,k)=npif
          endif
       enddo
    enddo
    call MPI_REDUCE(ising,mpi_aux_i,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,code)
    !if (nrank==0) write(*,*) '        number of points with potential problem in X :',mpi_aux_i
    !if (nrank==0) write(*,*) '    step 11'

    !y-pencil
    call transpose_x_to_y(ep1,ta2)
    nyipif(:,:,:)=npif
    nyfpif(:,:,:)=npif
    jsing=0
    do k=1,ysize(3)
       do i=1,ysize(1)
          jnum=0
          jflu=0
          if(ta2(i,1,k) ==  one)jnum=jnum+1
          if(ta2(i,1,k) == zero)jflu=jflu+1
          do j=2,ny
             if(ta2(i,j  ,k) == zero)jflu=jflu+1
             if(ta2(i,j-1,k) == zero .and.&
                ta2(i,j  ,k) ==  one)then
                jnum=jnum+1
                if(jnum == 1)then
                   nyipif(jnum  ,i,k)=jflu-izap
                   if(jflu-izap <  npif)jsing=jsing+1
                   if(jflu-izap >= npif)nyipif(jnum  ,i,k)=npif
                   jflu=0
                else
                   nyipif(jnum  ,i,k)=jflu-izap
                   nyfpif(jnum-1,i,k)=jflu-izap
                   if(jflu-izap <  npif)jsing=jsing+1
                   if(jflu-izap >= npif)nyipif(jnum  ,i,k)=npif
                   if(jflu-izap >= npif)nyfpif(jnum-1,i,k)=npif
                   jflu=0
                endif
             endif
             if(ta2(i,j,k) == one)jflu=0
          enddo
          if(ta2(i,ny,k) == zero)then
             nyfpif(jnum,i,k)=jflu-izap
             if(jflu-izap < npif)jsing=jsing+1
             if(jflu-izap < npif)nyfpif(jnum,i,k)=npif
          endif
       enddo
    enddo
    call MPI_REDUCE(jsing,mpi_aux_i,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,code)
    !if (nrank==0) write(*,*) '        number of points with potential problem in Y :',mpi_aux_i
    !if (nrank==0) write(*,*) '    step 12'

    !z-pencil
    if(nz > 1)then
       call transpose_y_to_z(ta2,ta3)
       nzipif(:,:,:)=npif
       nzfpif(:,:,:)=npif
       ksing=0
       do j=1,zsize(2)
          do i=1,zsize(1)
             knum=0
             kflu=0
             if(ta3(i,j,1) ==  one)knum=knum+1
             if(ta3(i,j,1) == zero)kflu=kflu+1
             do k=2,nz
                if(ta3(i,j,k  ) == zero)kflu=kflu+1
                if(ta3(i,j,k-1) == zero .and.&
                   ta3(i,j,k  ) ==  one)then
                   knum=knum+1
                   if(knum == 1)then
                      nzipif(knum  ,i,j)=kflu-izap
                      if(kflu-izap <  npif)ksing=ksing+1
                      if(kflu-izap >= npif)nzipif(knum  ,i,j)=npif
                      kflu=0
                   else
                      nzipif(knum  ,i,j)=kflu-izap
                      nzfpif(knum-1,i,j)=kflu-izap
                      if(kflu-izap <  npif)ksing=ksing+1
                      if(kflu-izap >= npif)nzipif(knum  ,i,j)=npif
                      if(kflu-izap >= npif)nzfpif(knum-1,i,j)=npif
                      kflu=0
                   endif
                endif
                if(ta3(i,j,k) == one )kflu=0
             enddo
             if(ta3(i,j,nz) == zero)then
                nzfpif(knum,i,j)=kflu-izap
                if(kflu-izap < npif)ksing=ksing+1
                if(kflu-izap < npif)nzfpif(knum,i,j)=npif
             endif
          enddo
       enddo
       call MPI_REDUCE(ksing,mpi_aux_i,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,code)
       !if (nrank==0) write(*,*) '        number of points with potential problem in Z :',mpi_aux_i
    endif
    !if (nrank==0) write(*,*) '    step 13'
    !
    return
  end subroutine verif_epsi
!############################################################################
!############################################################################
  subroutine geomcomplex_io(nx,ny,nz,ep1,nobjx,nobjy,nobjz,xi,xf,yi,yf,zi,zf,&
       nxipif,nxfpif,nyipif,nyfpif,nzipif,nzfpif,nobjmax,npif,read_flag)
    use decomp_2d
    USE decomp_2d_io
    implicit none
    !
    logical, intent(in) :: read_flag
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
    character :: tmp_char
    character(len=*), parameter :: io_geom = "io-geom"
    !###################################################################
    if (read_flag) then
      if (nrank==0) print *,'Reading geometry'
      call decomp_2d_read_one(1,ep1,'data/geometry','epsilon.bin',io_geom)   
    else
      if (nrank==0) print *,'Writing geometry'
      call decomp_2d_write_one(1,ep1,'data/geometry','epsilon.bin',0,io_geom)
    endif
    !###################################################################
    !x-pencil
    open(67,file='data/geometry/nobjx.dat',form='formatted',access='direct',recl=13)
    do k=xstart(3),xend(3)
       do j=xstart(2),xend(2)
          count = (k-1)*ny+j
          if (read_flag) then
            read(67,'(1I12,A)',rec=count) nobjx(j,k),tmp_char
          else
            write(67,'(1I12,A)',rec=count) nobjx(j,k),char(10)
          endif
       enddo
    enddo
    close(67)
    !y-pencil
    open(67,file='data/geometry/nobjy.dat',form='formatted',access='direct',recl=13)
    do k=ystart(3),yend(3)
       do i=ystart(1),yend(1)
          count = (k-1)*nx+i
          if (read_flag) then
            read(67,'(1I12,A)',rec=count) nobjy(i,k),tmp_char
          else
            write(67,'(1I12,A)',rec=count) nobjy(i,k),char(10)
         endif
       enddo
    enddo
    close(67)
    !z-pencil
    open(67,file='data/geometry/nobjz.dat',form='formatted',access='direct',recl=13)
    do j=zstart(2),zend(2)
       do i=zstart(1),zend(1)
          count = (j-1)*nx+i
          if (read_flag) then
            read(67,'(1I12,A)',rec=count) nobjz(i,j),tmp_char
          else
            write(67,'(1I12,A)',rec=count) nobjz(i,j),char(10)
         endif
       enddo
    enddo
    close(67)
    !###################################################################
    !x-pencil
    open(67,file='data/geometry/nxifpif.dat',form='formatted',access='direct',recl=25)
    do k=xstart(3),xend(3)
       do j=xstart(2),xend(2)
         do i=0,nobjmax
            count = (k-1)*ny*(1+nobjmax)+(j-1)*(1+nobjmax)+i+1
            if (read_flag) then
              read(67,'(2I12,A)',rec=count) nxipif(i,j,k),nxfpif(i,j,k),tmp_char
            else
              write(67,'(2I12,A)',rec=count) nxipif(i,j,k),nxfpif(i,j,k),char(10)
            endif
          enddo
       enddo
    enddo
    close(67)
    !y-pencil
    open(67,file='data/geometry/nyifpif.dat',form='formatted',access='direct',recl=25)
    do k=ystart(3),yend(3)
       do i=ystart(1),yend(1)
         do j=0,nobjmax
            count = (k-1)*nx*(1+nobjmax)+(i-1)*(1+nobjmax)+j+1
            if (read_flag) then
              read(67,'(2I12,A)',rec=count) nyipif(j,i,k),nyfpif(j,i,k),tmp_char
            else
              write(67,'(2I12,A)',rec=count) nyipif(j,i,k),nyfpif(j,i,k),char(10)
            endif
          enddo
       enddo
    enddo
    close(67)
    !z-pencil
    open(67,file='data/geometry/nzifpif.dat',form='formatted',access='direct',recl=25)
    do j=zstart(2),zend(2)
       do i=zstart(1),zend(1)
         do k=0,nobjmax
            count = (j-1)*nx*(1+nobjmax)+(i-1)*(1+nobjmax)+k+1
            if (read_flag) then
              read(67,'(2I12,A)',rec=count) nzipif(k,i,j),nzfpif(k,i,j),tmp_char
            else
              write(67,'(2I12,A)',rec=count) nzipif(k,i,j),nzfpif(k,i,j),char(10)
            endif
          enddo
       enddo
    enddo
    close(67)
    !###################################################################
    !x-pencil
    open(67,file='data/geometry/xixf.dat',form='formatted',access='direct',recl=49)
    do k=xstart(3),xend(3)
       do j=xstart(2),xend(2)
          do i=1,nobjmax
             count = (k-1)*ny*nobjmax+(j-1)*nobjmax+i
             if (read_flag) then
              read(67,'(2E24.16,A)',rec=count) xi(i,j,k),xf(i,j,k),tmp_char
            else
              write(67,'(2E24.16,A)',rec=count) xi(i,j,k),xf(i,j,k),char(10)
            endif
          enddo
       enddo
    enddo
    close(67)
    !y-pencil
    open(67,file='data/geometry/yiyf.dat',form='formatted',access='direct',recl=49)
    do k=ystart(3),yend(3)
       do i=ystart(1),yend(1)
          do j=1,nobjmax
             count = (k-1)*nx*nobjmax+(i-1)*nobjmax+j
             if (read_flag) then
               read(67,'(2E24.16,A)',rec=count) yi(j,i,k),yf(j,i,k),tmp_char
            else
               write(67,'(2E24.16,A)',rec=count) yi(j,i,k),yf(j,i,k),char(10)
            endif
          enddo
       enddo
    enddo
    close(67)
    !z-pencil
    open(67,file='data/geometry/zizf.dat',form='formatted',access='direct',recl=49)
    do j=zstart(2),zend(2)
       do i=zstart(1),zend(1)
          do k=1,nobjmax
             count = (j-1)*nx*nobjmax+(i-1)*nobjmax+k
             if (read_flag) then
               read(67,'(2E24.16,A)',rec=count) zi(k,i,j),zf(k,i,j),tmp_char
            else
               write(67,'(2E24.16,A)',rec=count) zi(k,i,j),zf(k,i,j),char(10)
            endif
          enddo
       enddo
    enddo
    close(67)
    return
  end subroutine geomcomplex_io
  !############################################################################
end module genepsi
