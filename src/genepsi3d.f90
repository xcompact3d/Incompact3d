#ifdef IBM
subroutine genepsi3d(ep1)

  USE variables, only : nx,ny,nz,nxm,nym,nzm,yp
  USE param, only : xlx,yly,zlz,dx,dy,dz,izap,npif,nclx,ncly,nclz,istret
  USE complex_geometry
  use decomp_2d

  implicit none

  !*****************************************************************!
  ! 0- This program will generate all the files necessary for our
  !    customize IMB based on Lagrange reconstructions
  ! 1- It has to be use with ilag=1 in incompact3d.prm
  ! 3- The object is defined in the cylinder subroutine
  ! 4- You can add your own subroutine for your own object
  ! 7- Please cite the following paper if you are using this file:
  ! Gautier R., Laizet S. & Lamballais E., 2014, A DNS study of
  ! jet control with microjets using an alterna ng direc on forcing
  ! strategy, Int. J. of Computa onal Fluid Dynamics, 28, 393--410
  !*****************************************************************!
  !
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ep1
  !
  if (nrank==0) print *,'Generating the geometry!'
  call gene_epsi_3D(ep1,nx,ny,nz,dx,dy,dz,xlx,yly,zlz ,&
       nclx,ncly,nclz,nxraf,nyraf,nzraf   ,&
       xi,xf,yi,yf,zi,zf,nobjx,nobjy,nobjz,&
       nobjmax,yp,nraf)
  call verif_epsi(ep1,npif,izap,nx,ny,nz,nobjmax,&
       nxipif,nxfpif,nyipif,nyfpif,nzipif,nzfpif)
#ifdef DEBG
  call write_geomcomplex(nx,ny,nz,ep1,nobjx,nobjy,nobjz,xi,xf,yi,yf,zi,zf,&
       nxipif,nxfpif,nyipif,nyfpif,nzipif,nzfpif,nobjmax,npif)
#endif
  !
end subroutine genepsi3d
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine gene_epsi_3D(ep1,nx,ny,nz,dx,dy,dz,xlx,yly,zlz ,&
     nclx,ncly,nclz,nxraf,nyraf,nzraf   ,&
     xi,xf,yi,yf,zi,zf,nobjx,nobjy,nobjz,&
     nobjmax,yp,nraf)
  use param, only : zero,one
  use decomp_2d
  use MPI
  implicit none
  !
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ep1
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
  if (nrank==0) print*,'    step 1'
  if(nclx)then
     dxraf =xlx/nxraf
  else
     dxraf =xlx/(nxraf-1)
  endif
  xepsi=zero
  call geomcomplex(xepsi,1,nxraf,ny,xstart(2),xend(2),xstart(3),xend(3),dxraf,yp,dz,one)
  if (nrank==0) print*,'    step 2'
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
  if (nrank==0) print*,'    step 3'
  !z-pencil
  if(nclz)then
     dzraf=zlz/nzraf
  else
     dzraf=zlz/(nzraf-1)
  endif
  zepsi=zero
  call geomcomplex(zepsi,zstart(1),zend(1),ny,zstart(2),zend(2),1,nzraf,dx,yp,dzraf,one)
  if (nrank==0) print*,'    step 4'

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
  if (nrank==0) print*,'        nobjxmax=',mpi_aux_i

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
           if(xepsi(i,j,k).eq.0..and.xepsi(i+1,j,k).eq.1.)then
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
  if (nrank==0) print*,'        nobjxmaxraf=',mpi_aux_i
  call MPI_REDUCE(ibug,mpi_aux_i,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,code)
  if (nrank==0) print*,'        ibug=',mpi_aux_i
  if (nrank==0) print*,'    step 5'

  !y-pencil
  nobjy(:,:)=0
  nobjymax=0
  call transpose_x_to_y(ep1,ep2)
  do k=1,ysize(3)
     do i=1,ysize(1)
        jnum=0
        if(ep2(i,1,k).eq.1.)then
           jnum=1
           nobjy(i,k)=1
        endif
        do j=1,ny-1
           if(ep2(i,j,k).eq.0..and.ep2(i,j+1,k).eq.1.)then
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
  if (nrank==0) print*,'        nobjymax=',mpi_aux_i

  nobjyraf(:,:)=0
  jbug=0
  nobjymaxraf=0
  jnum=0
  do k=1,ysize(3)
     do i=1,ysize(1)
        jnum=0
        if(yepsi(i,1,k).eq.1.)then
           jnum=1
           nobjyraf(i,k)=1
        endif
        do j=1,nyraf-1
           if(yepsi(i,j,k).eq.0..and.yepsi(i,j+1,k).eq.1.)then
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
  if (nrank==0) print*,'        nobjymaxraf=',mpi_aux_i
  call MPI_REDUCE(jbug,mpi_aux_i,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,code)
  if (nrank==0) print*,'        jbug=',mpi_aux_i
  if (nrank==0) print*,'    step 6'

  !z-pencil
  nobjz(:,:)=0
  nobjzmax=0
  call transpose_y_to_z(ep2,ep3)
  do j=1,zsize(2)
     do i=1,zsize(1)
        knum=0
        if(ep3(i,j,1).eq.1.)then
           knum=1
           nobjz(i,j)=1
        endif
        do k=1,nz-1
           if(ep3(i,j,k).eq.0..and.ep3(i,j,k+1).eq.1.)then
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
  if (nrank==0) print*,'        nobjzmax=',mpi_aux_i

  nobjzraf(:,:)=0
  kbug=0
  nobjzmaxraf=0
  knum=0
  do j=1,zsize(2)
     do i=1,zsize(1)
        knum=0
        if(zepsi(i,j,1).eq.1.)then
           knum=1
           nobjzraf(i,j)=1
        endif
        do k=1,nzraf-1
           if(zepsi(i,j,k).eq.0..and.zepsi(i,j,k+1).eq.1.)then
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
  if (nrank==0) print*,'        nobjzmaxraf=',mpi_aux_i
  call MPI_REDUCE(kbug,mpi_aux_i,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,code)
  if (nrank==0) print*,'        kbug=',mpi_aux_i
  if (nrank==0) print*,'    step 7'

  !x-pencil
  do k=1,xsize(3)
     do j=1,xsize(2)
        inum=0
        if(xepsi(1,j,k).eq.1.)then
           inum=inum+1
           xi(inum,j,k)=-dx!-xlx
        endif
        do i=1,nxraf-1
           if(xepsi(i,j,k).eq.0..and.xepsi(i+1,j,k).eq.1.)then
              inum=inum+1
              xi(inum,j,k)=dxraf*(i-1)+dxraf/2.
           elseif(xepsi(i,j,k).eq.1..and.xepsi(i+1,j,k).eq.0.)then
              xf(inum,j,k)=dxraf*(i-1)+dxraf/2.
           endif
        enddo
        if(xepsi(nxraf,j,k).eq.1.)then
           xf(inum,j,k)=xlx+dx!2.*xlx
        endif
     enddo
  enddo

  if(ibug.ne.0)then
     do k=1,xsize(3)
        do j=1,xsize(2)
           if(nobjx(j,k).ne.nobjxraf(j,k))then
              iobj=0
              if(ep1(1,j,k).eq.1.)iobj=iobj+1
              do i=1,nx-1
                 if(ep1(i,j,k).eq.0..and.ep1(i+1,j,k).eq.1.)iobj=iobj+1
                 if(ep1(i,j,k).eq.0..and.ep1(i+1,j,k).eq.0.)iflu=1
                 if(ep1(i,j,k).eq.1..and.ep1(i+1,j,k).eq.1.)isol=1
                 do iraf=1,nraf
                    if(xepsi(iraf+nraf*(i-1)  ,j,k).eq.0..and.&
                         xepsi(iraf+nraf*(i-1)+1,j,k).eq.1.)idebraf=iraf+nraf*(i-1)+1
                    if(xepsi(iraf+nraf*(i-1)  ,j,k).eq.1..and.&
                         xepsi(iraf+nraf*(i-1)+1,j,k).eq.0.)ifinraf=iraf+nraf*(i-1)+1
                 enddo
                 if(idebraf.ne.0.and.ifinraf.ne.0.and.&
                      idebraf.lt.ifinraf.and.iflu.eq.1)then
                    iobj=iobj+1
                    do ii=iobj,nobjmax-1
                       xi(ii,j,k)=xi(ii+1,j,k)
                       xf(ii,j,k)=xf(ii+1,j,k)
                    enddo
                    iobj=iobj-1
                 endif
                 if(idebraf.ne.0.and.ifinraf.ne.0.and.&
                      idebraf.gt.ifinraf.and.isol.eq.1)then
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
  if (nrank==0) print*,'    step 8'

  !y-pencil
  do k=1,ysize(3)
     do i=1,ysize(1)
        jnum=0
        if(yepsi(i,1,k).eq.1.)then
           jnum=jnum+1
           yi(jnum,i,k)=-(yp(2)-yp(1))!-yly
        endif
        do j=1,nyraf-1
           if(yepsi(i,j,k).eq.0..and.yepsi(i,j+1,k).eq.1.)then
              jnum=jnum+1
              yi(jnum,i,k)=ypraf(j)+(ypraf(j+1)-ypraf(j))/2.!dyraf*(j-1)+dyraf/2.
           elseif(yepsi(i,j,k).eq.1..and.yepsi(i,j+1,k).eq.0.)then
              yf(jnum,i,k)=ypraf(j)+(ypraf(j+1)-ypraf(j))/2.!dyraf*(j-1)+dyraf/2.
           endif
        enddo
        if(yepsi(i,nyraf,k).eq.1.)then
           yf(jnum,i,k)=yly+(yp(ny)-yp(ny-1))/2.!2.*yly
        endif
     enddo
  enddo

  if(jbug.ne.0)then
     do k=1,ysize(3)
        do i=1,ysize(1)
           if(nobjy(i,k).ne.nobjyraf(i,k))then
              jobj=0
              if(ep2(i,1,k).eq.1.)jobj=jobj+1
              do j=1,ny-1
                 if(ep2(i,j,k).eq.0..and.ep2(i,j+1,k).eq.1.)jobj=jobj+1
                 if(ep2(i,j,k).eq.0..and.ep2(i,j+1,k).eq.0.)jflu=1
                 if(ep2(i,j,k).eq.1..and.ep2(i,j+1,k).eq.1.)jsol=1
                 do jraf=1,nraf
                    if(yepsi(i,jraf+nraf*(j-1)  ,k).eq.0..and.&
                         yepsi(i,jraf+nraf*(j-1)+1,k).eq.1.)jdebraf=jraf+nraf*(j-1)+1
                    if(yepsi(i,jraf+nraf*(j-1)  ,k).eq.1..and.&
                         yepsi(i,jraf+nraf*(j-1)+1,k).eq.0.)jfinraf=jraf+nraf*(j-1)+1
                 enddo
                 if(jdebraf.ne.0.and.jfinraf.ne.0.and.&
                      jdebraf.lt.jfinraf.and.jflu.eq.1)then
                    jobj=jobj+1
                    do jj=jobj,nobjmax-1
                       yi(jj,i,k)=yi(jj+1,i,k)
                       yf(jj,i,k)=yf(jj+1,i,k)
                    enddo
                    jobj=jobj-1
                 endif
                 if(jdebraf.ne.0.and.jfinraf.ne.0.and.&
                      jdebraf.gt.jfinraf.and.jsol.eq.1)then
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
  if (nrank==0) print*,'    step 9'

  !z-pencil
  do j=1,zsize(2)
     do i=1,zsize(1)
        knum=0
        if(zepsi(i,j,1).eq.1.)then
           knum=knum+1
           zi(knum,i,j)=-dz!zlz
        endif
        do k=1,nzraf-1
           if(zepsi(i,j,k).eq.0..and.zepsi(i,j,k+1).eq.1.)then
              knum=knum+1
              zi(knum,i,j)=dzraf*(k-1)+dzraf/2.
           elseif(zepsi(i,j,k).eq.1..and.zepsi(i,j,k+1).eq.0.)then
              zf(knum,i,j)=dzraf*(k-1)+dzraf/2.
           endif
        enddo
        if(zepsi(i,j,nzraf).eq.1.)then
           zf(knum,i,j)=zlz+dz!2.*zlz
        endif
     enddo
  enddo

  if(kbug.ne.0)then
     do j=1,zsize(2)
        do i=1,zsize(1)
           if(nobjz(i,j).ne.nobjzraf(i,j))then
              kobj=0
              if(ep3(i,j,1).eq.1.)kobj=kobj+1
              do k=1,nz-1
                 if(ep3(i,j,k).eq.0..and.ep3(i,j,k+1).eq.1.)kobj=kobj+1
                 if(ep3(i,j,k).eq.0..and.ep3(i,j,k+1).eq.0.)kflu=1
                 if(ep3(i,j,k).eq.1..and.ep3(i,j,k+1).eq.1.)ksol=1
                 do kraf=1,nraf
                    if(zepsi(i,j,kraf+nraf*(k-1)  ).eq.0..and.&
                         zepsi(i,j,kraf+nraf*(k-1)+1).eq.1.)kdebraf=kraf+nraf*(k-1)+1
                    if(zepsi(i,j,kraf+nraf*(k-1)  ).eq.1..and.&
                         zepsi(i,j,kraf+nraf*(k-1)+1).eq.0.)kfinraf=kraf+nraf*(k-1)+1
                 enddo
                 if(kdebraf.ne.0.and.kfinraf.ne.0.and.&
                      kdebraf.lt.kfinraf.and.kflu.eq.1)then
                    kobj=kobj+1
                    do kk=kobj,nobjmax-1
                       zi(kk,i,j)=zi(kk+1,i,j)
                       zf(kk,i,j)=zf(kk+1,i,j)
                    enddo
                    kobj=kobj-1
                 endif
                 if(kdebraf.ne.0.and.kfinraf.ne.0.and.&
                      kdebraf.gt.kfinraf.and.ksol.eq.1)then
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
  if (nrank==0) print*,'    step 10'
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
           if(epsi(i,j,k).eq.1..and.&
                nobjx (j,k).ne.0)then
              x=dx*(i-1)
              do ix=1,nobjx(j,k)
                 if(x  .ge.xi(ix,j,k)   .and.&
                      x  .le.xi(ix,j,k)+xe.and.&
                      0. .lt.xi(ix,j,k)   .or. &
                      x  .le.xf(ix,j,k)   .and.&
                      x  .ge.xf(ix,j,k)-xe.and.&
                      xlx.gt.xf(ix,j,k)   )then
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
           if(epsi(i,j,k).eq.1..and.&
                nobjy (i,k).gt.0)then
              y=yp(j)
              ye=epm*(yp(j+1)-yp(j-1))/2.
              do jy=1,nobjy(i,k)
                 if(y  .ge.yi(jy,i,k)   .and.&
                      y  .le.yi(jy,i,k)+ye.and.&
                      0. .lt.yi(jy,i,k)   .or. &
                      y  .le.yf(jy,i,k)   .and.&
                      y  .ge.yf(jy,i,k)-ye.and.&
                      yly.gt.yf(jy,i,k)   )then
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
           if(epsi(i,j,k).eq.1..and.&
                nobjz (i,j).gt.0)then
              z=dz*(k-1)
              do kz=1,nobjz(i,j)
                 if(z  .ge.zi(kz,i,j)   .and.&
                      z  .le.zi(kz,i,j)+ze.and.&
                      0. .lt.zi(kz,i,j)   .or. &
                      z  .le.zf(kz,i,j)   .and.&
                      z  .ge.zf(kz,i,j)-ze.and.&
                      zlz.gt.zf(kz,i,j)   )then
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
        if(ep1(1,j,k).eq.1.)inum=inum+1
        if(ep1(1,j,k).eq.0.)iflu=iflu+1
        do i=2,nx
           if(ep1(i  ,j,k).eq.0.)iflu=iflu+1
           if(ep1(i-1,j,k).eq.0..and.&
                ep1(i  ,j,k).eq.1.)then
              inum=inum+1
              if(inum.eq.1)then
                 nxipif(inum  ,j,k)=iflu-izap
                 if(iflu-izap.lt.npif)ising=ising+1
                 if(iflu-izap.ge.npif)nxipif(inum  ,j,k)=npif
                 iflu=0
              else
                 nxipif(inum  ,j,k)=iflu-izap
                 nxfpif(inum-1,j,k)=iflu-izap
                 if(iflu-izap.lt.npif)ising=ising+1
                 if(iflu-izap.ge.npif)nxipif(inum  ,j,k)=npif
                 if(iflu-izap.ge.npif)nxfpif(inum-1,j,k)=npif
                 iflu=0
              endif
           endif
           if(ep1(i,j,k).eq.1.)iflu=0
        enddo
        if(ep1(nx,j,k).eq.0.)then
           nxfpif(inum,j,k)=iflu-izap
           if(iflu-izap.lt.npif)ising=ising+1
           if(iflu-izap.lt.npif)nxfpif(inum,j,k)=npif
        endif
     enddo
  enddo
  call MPI_REDUCE(ising,mpi_aux_i,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,code)
  if (nrank==0) print*,'        number of points with potential problem in X :',mpi_aux_i
  if (nrank==0) print*,'    step 11'

  !y-pencil
  call transpose_x_to_y(ep1,ep2)
  nyipif(:,:,:)=npif
  nyfpif(:,:,:)=npif
  jsing=0
  do k=1,ysize(3)
     do i=1,ysize(1)
        jnum=0
        jflu=0
        if(ep2(i,1,k).eq.1.)jnum=jnum+1
        if(ep2(i,1,k).eq.0.)jflu=jflu+1
        do j=2,ny
           if(ep2(i,j  ,k).eq.0.)jflu=jflu+1
           if(ep2(i,j-1,k).eq.0..and.&
                ep2(i,j  ,k).eq.1.)then
              jnum=jnum+1
              if(jnum.eq.1)then
                 nyipif(jnum  ,i,k)=jflu-izap
                 if(jflu-izap.lt.npif)jsing=jsing+1
                 if(jflu-izap.ge.npif)nyipif(jnum  ,i,k)=npif
                 jflu=0
              else
                 nyipif(jnum  ,i,k)=jflu-izap
                 nyfpif(jnum-1,i,k)=jflu-izap
                 if(jflu-izap.lt.npif)jsing=jsing+1
                 if(jflu-izap.ge.npif)nyipif(jnum  ,i,k)=npif
                 if(jflu-izap.ge.npif)nyfpif(jnum-1,i,k)=npif
                 jflu=0
              endif
           endif
           if(ep2(i,j,k).eq.1.)jflu=0
        enddo
        if(ep2(i,ny,k).eq.0.)then
           nyfpif(jnum,i,k)=jflu-izap
           if(jflu-izap.lt.npif)jsing=jsing+1
           if(jflu-izap.lt.npif)nyfpif(jnum,i,k)=npif
        endif
     enddo
  enddo
  call MPI_REDUCE(jsing,mpi_aux_i,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,code)
  if (nrank==0) print*,'        number of points with potential problem in Y :',mpi_aux_i
  if (nrank==0) print*,'    step 12'

  !z-pencil
  if(nz.gt.1)then
     call transpose_y_to_z(ep2,ep3)
     nzipif(:,:,:)=npif
     nzfpif(:,:,:)=npif
     ksing=0
     do j=1,zsize(2)
        do i=1,zsize(1)
           knum=0
           kflu=0
           if(ep3(i,j,1).eq.1.)knum=knum+1
           if(ep3(i,j,1).eq.0.)kflu=kflu+1
           do k=2,nz
              if(ep3(i,j,k  ).eq.0.)kflu=kflu+1
              if(ep3(i,j,k-1).eq.0..and.&
                   ep3(i,j,k  ).eq.1.)then
                 knum=knum+1
                 if(knum.eq.1)then
                    nzipif(knum  ,i,j)=kflu-izap
                    if(kflu-izap.lt.npif)ksing=ksing+1
                    if(kflu-izap.ge.npif)nzipif(knum  ,i,j)=npif
                    kflu=0
                 else
                    nzipif(knum  ,i,j)=kflu-izap
                    nzfpif(knum-1,i,j)=kflu-izap
                    if(kflu-izap.lt.npif)ksing=ksing+1
                    if(kflu-izap.ge.npif)nzipif(knum  ,i,j)=npif
                    if(kflu-izap.ge.npif)nzfpif(knum-1,i,j)=npif
                    kflu=0
                 endif
              endif
              if(ep3(i,j,k).eq.1.)kflu=0
           enddo
           if(ep3(i,j,nz).eq.0.)then
              nzfpif(knum,i,j)=kflu-izap
              if(kflu-izap.lt.npif)ksing=ksing+1
              if(kflu-izap.lt.npif)nzfpif(knum,i,j)=npif
           endif
        enddo
     enddo
     call MPI_REDUCE(ksing,mpi_aux_i,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,code)
     if (nrank==0) print*,'        number of points with potential problem in Z :',mpi_aux_i
  endif
  if (nrank==0) print*,'    step 13'
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
  !
  if (nrank==0) print *,'Writing geometry'
  call decomp_2d_write_one(1,ep1,'epsilon.dat')
  !x-pencil
  open(67,file='nobjx.dat',form='formatted',access='direct',recl=13)
  do k=xstart(3),xend(3)
     do j=xstart(2),xend(2)
        count = (k-1)*ny+j
        write(67,'(1I12,A)',rec=count) nobjx(j,k),char(10)
     enddo
  enddo
  close(67)
  !y-pencil
  open(67,file='nobjy.dat',form='formatted',access='direct',recl=13)
  do k=ystart(3),yend(3)
     do i=ystart(1),yend(1)
        count = (k-1)*nx+i
        write(67,'(1I12,A)',rec=count) nobjy(i,k),char(10)
     enddo
  enddo
  close(67)
  !z-pencil
  open(67,file='nobjz.dat',form='formatted',access='direct',recl=13)
  do j=zstart(2),zend(2)
     do i=zstart(1),zend(1)
        count = (j-1)*nx+i
        write(67,'(1I12,A)',rec=count) nobjz(i,j),char(10)
     enddo
  enddo
  close(67)
  !x-pencil
  open(67,file='nxifpif.dat',form='formatted',access='direct',recl=25)
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
  open(67,file='nyifpif.dat',form='formatted',access='direct',recl=25)
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
  open(67,file='nzifpif.dat',form='formatted',access='direct',recl=25)
  do j=zstart(2),zend(2)
     do i=zstart(1),zend(1)
        do k=0,nobjmax
           count = (j-1)*nx*(1+nobjmax)+(i-1)*(1+nobjmax)+k+1
           write(67,'(2I12,A)',rec=count) nzipif(k,i,j),nzfpif(k,i,j),char(10)
        enddo
     enddo
  enddo
  close(67)
  !x-pencil
  open(67,file='xixf.dat',form='formatted',access='direct',recl=29)
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
  open(67,file='yiyf.dat',form='formatted',access='direct',recl=29)
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
  open(67,file='zizf.dat',form='formatted',access='direct',recl=29)
  do j=zstart(2),zend(2)
     do i=zstart(1),zend(1)
        do k=1,nobjmax
           count = (j-1)*nx*nobjmax+(i-1)*nobjmax+k
           write(67,'(2E14.6,A)',rec=count) zi(k,i,j),zf(k,i,j),char(10)
        enddo
     enddo
  enddo
  close(67)
  !
  return
end subroutine write_geomcomplex
#endif
