module flow_type
  use decomp_2d, only : mytype

end module flow_type

subroutine ft_parameter(arg)

  USE param
  USE variables
  USE flow_type
  USE complex_geometry
  USE decomp_2d, only : nrank
  implicit none

  logical,intent(in) :: arg
  character :: a

  open(10,file='BC-TGV.prm',status='unknown',form='formatted')
  read (10,*) a !
  read (10,*) a ! INCOMPACT 3D computational parameters
  read (10,*) a !
  read (10,*) nx
  ny = nx
  nz = nx
  read (10,*) nphi
  read (10,*) p_row
  read (10,*) p_col
  if (arg) then
    close(10)
    return
  endif
  read (10,*) a !
  read (10,*) a ! INCOMPACT 3D Flow parameters
  read (10,*) a !
  read (10,*) xlx
  xlx = xlx*pi
  yly = xlx
  zlz = xlx
  read (10,*) re
  read (10,*) ri
  read (10,*) noise
  read (10,*) dt
  read (10,*) a !
  read (10,*) a ! INCOMPACT3D Flow configuration
  read (10,*) a !
  read (10,*) iscalar
  read (10,*) cont_phi
  read (10,*) jLES
  read (10,*) iin
  read (10,*) ifirst
  read (10,*) ilast
  read (10,*) nscheme
  read (10,*) a !velocity
  read (10,*) nclx1
  read (10,*) nclxn
  read (10,*) ncly1
  read (10,*) nclyn
  read (10,*) nclz1
  read (10,*) nclzn
  read (10,*) a !scalar
  read (10,*) nclxS1
  read (10,*) nclxSn
  read (10,*) nclyS1
  read (10,*) nclySn
  read (10,*) nclzS1
  read (10,*) nclzSn
  read (10,*) a !
  read (10,*) a ! INCOMPACT 3D File parameters
  read (10,*) a !
  read (10,*) ilit
  read (10,*) isave
  read (10,*) imodulo
  read (10,*) imodulo2
  read (10,*) a !
  read (10,*) a ! NUMERICAL DISSIPATION
  read (10,*) a !
  read (10,*) fpi2

  return
end subroutine ft_parameter
!********************************************************************
subroutine init (ux1,uy1,uz1,ep1,phi1,dux1,duy1,duz1,phis1,phiss1)

  USE decomp_2d
  USE decomp_2d_io
  USE variables
  USE param
  USE MPI

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1,phis1,phiss1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1

  real(mytype) :: y,r,um,r3,x,z,h,ct
  real(mytype) :: cx0,cy0,cz0,hg,lg
  integer :: k,j,i,ijk,fh,ierror,is,code
  integer (kind=MPI_OFFSET_KIND) :: disp
  integer, dimension (:), allocatable :: seed
  integer ::  isize

  if (iscalar==1) then

     phi1 = one

  endif

  if (iin.eq.0) then !empty domain

     if (nrank==0) write(*,*) "Empty initial domain!"

     ux1=zero; uy1=zero; uz1=zero

  endif

  if (iin.eq.1) then !generation of a random noise

     if (nrank==0) write(*,*) "Filled initial domain!"

     ux1=zero; uy1=zero; uz1=zero

     !t=zero
     !xv=1._mytype/10zero
     !xxk1=twopi/xlx
     !xxk2=twopi/yly
     do k=1,xsize(3)
        z=real((k+xstart(3)-1-1),mytype)*dz
        do j=1,xsize(2)
           y=real((j+xstart(2)-1-1),mytype)*dy
           do i=1,xsize(1)
              x=real(i-1,mytype)*dx
              !ux1(i,j,k)=sin(2.*pi*x)*cos(2.*pi*y)*cos(2.*pi*z)
              !sin(xxk1*x)*cos(xxk2*y)*exp(-(xxk1*xxk1+xxk2*xxk2)*xv*t)
              !uy1(i,j,k)=sin(2.*pi*y)*cos(2.*pi*x)*cos(2.*pi*z)
              !-xxk1/xxk2*sin(xxk2*y)*cos(xxk1*x)*exp(-(xxk1*xxk1+xxk2*xxk2)*xv*t)
              !uz1(i,j,k)=sin(2.*pi*z)*cos(2.*pi*x)*cos(2.*pi*y)
          ux1(i,j,k)=+sin(x)*cos(y)*cos(z)
          uy1(i,j,k)=-cos(x)*sin(y)*cos(z)
          uz1(i,j,k)=zero
           enddo
        enddo
     enddo

          call random_seed(size=isize)
         allocate (seed(isize))
        seed(:)=67
         call random_seed(put=seed)
     !     call random_number(ux1)
     !     call random_number(uy1)
          call random_number(uz1)

          do k=1,xsize(3)
             do j=1,xsize(2)
                do i=1,xsize(1)
     !              ux1(i,j,k)=noise*(ux1(i,j,k)-half)
     !              uy1(i,j,k)=noise*(uy1(i,j,k)-half)
                   uz1(i,j,k)=0.05*(uz1(i,j,k)-half)
                enddo
             enddo
          enddo

     !     !modulation of the random noise
     !     do k=1,xsize(3)
     !        do j=1,xsize(2)
     !           if (istret.eq.0) y=(j+xstart(2)-1-1)*dy-yly/two
     !           if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/two
     !           um=exp(-0.2*y*y)
     !           do i=1,xsize(1)
     !              ux1(i,j,k)=um*ux1(i,j,k)
     !              uy1(i,j,k)=um*uy1(i,j,k)
     !              uz1(i,j,k)=um*uz1(i,j,k)
     !           enddo
     !        enddo
     !     enddo

  endif

  if (iscalar==1) then
     do is=1,nphi !vectorized for performance Ricardo
        do ijk=1,xsize(1)*xsize(2)*xsize(3)
           phis1(ijk,1,1,is)=phi1(ijk,1,1,is)
           phiss1(ijk,1,1,is)=phis1(ijk,1,1,is)
        enddo
     enddo
  endif

  !  bxx1(j,k)=zero
  !  bxy1(j,k)=zero
  !  bxz1(j,k)=zero

  !INIT FOR G AND U=MEAN FLOW + NOISE
  do k=1,xsize(3)
     do j=1,xsize(2)
        do i=1,xsize(1)
           ux1(i,j,k)=ux1(i,j,k)!+bxx1(j,k)
           uy1(i,j,k)=uy1(i,j,k)!+bxy1(j,k)
           uz1(i,j,k)=uz1(i,j,k)!+bxz1(j,k)
           dux1(i,j,k,1)=ux1(i,j,k)
           duy1(i,j,k,1)=uy1(i,j,k)
           duz1(i,j,k,1)=uz1(i,j,k)
           dux1(i,j,k,2)=dux1(i,j,k,1)
           duy1(i,j,k,2)=duy1(i,j,k,1)
           duz1(i,j,k,2)=duz1(i,j,k,1)
        enddo
     enddo
  enddo

#ifdef DEBG
  if (nrank .eq. 0) print *,'# init end ok'
#endif

  return
end subroutine init
!********************************************************************
subroutine boundary_conditions (ux,uy,uz,phi)

  USE param
  USE variables
  USE decomp_2d

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi


end subroutine boundary_conditions
!********************************************************************
module post_processing

  USE decomp_2d
  USE variables
  USE param

  implicit none

  real(mytype), save, allocatable, dimension(:,:,:) :: vol1,volSimps1
  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  !probes !só vai funcionar se a impressão for em relação ao lapis X!
  integer :: nprobes
  integer, save, allocatable, dimension(:) :: rankprobes, nxprobes, nyprobes, nzprobes

contains

  !############################################################################
  subroutine init_post(ep1)

    USE MPI

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ep1
    real(mytype) :: dxdydz, dxdz, x, xprobes, yprobes, zprobes
    integer :: i,j,k,code
    character :: a

    call alloc_x(vol1, opt_global=.true.)
    call alloc_x(volSimps1, opt_global=.true.)
    
    vol1 = zero
    volSimps1 = zero

    !X PENCILS !Utilizar para integral volumétrica dentro do domínio físico
    dxdydz=dx*dy*dz
    do k=xstart(3),xend(3)
      do j=xstart(2),xend(2)
        do i=xstart(1),xend(1)
          vol1(i,j,k)=dxdydz
          if (i .eq. 1   .or. i .eq. nx) vol1(i,j,k) = vol1(i,j,k)/two
          if (j .eq. 1   .or. j .eq. ny)  vol1(i,j,k) = vol1(i,j,k)/two
          if (k .eq. 1   .or. k .eq. nz)  vol1(i,j,k) = vol1(i,j,k)/two
        end do
      end do
    end do

    !X PENCILS !Utilizar para integral volumétrica dentro do domínio físico (método de Simpson)
    do k=xstart(3),xend(3)
      do j=xstart(2),xend(2)
        do i=xstart(1),xend(1)
             volSimps1(i,j,k)=dxdydz
             if (i .eq. 1   .or. i .eq. nx) volSimps1(i,j,k) = volSimps1(i,j,k) * (five/twelve)
             if (j .eq. 1   .or. j .eq. ny) volSimps1(i,j,k) = volSimps1(i,j,k) * (five/twelve)
             if (k .eq. 1   .or. k .eq. nz) volSimps1(i,j,k) = volSimps1(i,j,k) * (five/twelve)
             if (i .eq. 2   .or. i .eq. nx-1) volSimps1(i,j,k) = volSimps1(i,j,k) * (thirteen/twelve)
             if (j .eq. 2   .or. j .eq. ny-1) volSimps1(i,j,k) = volSimps1(i,j,k) * (thirteen/twelve)
             if (k .eq. 2   .or. k .eq. nz-1) volSimps1(i,j,k) = volSimps1(i,j,k) * (thirteen/twelve)
          end do
       end do
    end do
    
    !probes
    !WORK X-PENCILS
    open(10,file='probes.prm',status='unknown',form='formatted')
    read (10,*) nprobes
    read (10,*) a
    if (nprobes .gt. 0) then
      allocate(nxprobes(nprobes), nyprobes(nprobes), nzprobes(nprobes), rankprobes(nprobes))
      rankprobes(:)=0
      do i=1, nprobes
        read (10,*) xprobes, yprobes, zprobes
        !x
        if (nclx) then
          nxprobes(i)=int(xprobes/dx)
        else
          nxprobes(i)=int(xprobes/dx+1)
        end if
        !y
        if (ncly) then
          nyprobes(i)=int(yprobes/dy)
        else
          nyprobes(i)=int(yprobes/dy+1)
        end if
        !z
        if (nclz) then
          nzprobes(i)=int(zprobes/dz)
        else
          nzprobes(i)=int(zprobes/dz+1)
        end if
        if     (xstart(1) .le. nxprobes(i) .and. nxprobes(i) .le. xend(1)) then
          if   (xstart(2) .le. nyprobes(i) .and. nyprobes(i) .le. xend(2)) then
            if (xstart(3) .le. nzprobes(i) .and. nzprobes(i) .le. xend(3)) then
              rankprobes(i)=1
            endif
          endif
        endif
      enddo
    endif
    close(10)

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init_post ok'
#endif

  end subroutine init_post
  !############################################################################
  subroutine write_probes(ux1,uy1,uz1,phi1) !By Felipe Schuch

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ux1, uy1, uz1
    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),nphi) :: phi1

    integer :: i
    character(len=30) :: filename
    FS = 1+3+nphi !Number of columns
    write(fileformat, '( "(",I4,"(E14.6),A)" )' ) FS
    FS = FS*14+1  !Line width

    do i=1, nprobes
      if (rankprobes(i) .eq. 1) then
        write(filename,"('./out/probe',I4.4)") i
        open(67,file=trim(filename),status='unknown',form='formatted'&
        ,access='direct',recl=FS)
        write(67,fileformat,rec=itime) t,&                         !1
        ux1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !2
        uy1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !3
        uz1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !4
        phi1(nxprobes(i),nyprobes(i),nzprobes(i),:),&         !nphi
        NL                                                    !+1
        close(67)
      endif
    enddo

  end subroutine write_probes
  !############################################################################
  subroutine postprocessing(ux1,uy1,uz1,phi1,ep1) !By Felipe Schuch

    USE decomp_2d
    USE MPI

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1
    
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,ep1,diss1
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
    real(mytype) :: mp(nphi),mps(nphi),vl,es,es1,ek,ek1,ds,ds1

    integer :: i,j,k,is,ijk,code,nvect1
    nvect1=xsize(1)*xsize(2)*xsize(3)

    if (iscalar==1) then

    call suspended(phi1,vol1,mp)
    call suspended(phi1,volSimps1,mps)
  
    if (nrank .eq. 0) write(*,*) 'TRAP SIMPS Mp Mps',mp(1),mps(1)

    if (nrank .eq. 0) then
      FS = 1+1+1 !Number of columns
      write(fileformat, '( "(",I4,"(E14.6),A)" )' ) FS
      FS = FS*14+1  !Line width
      open(67,file='./out/suspended',status='unknown',form='formatted',&
      access='direct',recl=FS)
      write(67,fileformat,rec=itime+1) t,& !1
      mp,&                                 !2
      mps,&                                !3
      NL                                   !+1
      close(67)
    end if

    endif

     if (mod(itime,imodulo2).eq.0) then
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

    es=zero;es1=zero
    !ENSTROPHY
    di1=zero
    do ijk=1,nvect1
      di1(ijk,1,1) = (tf1(ijk,1,1)-th1(ijk,1,1))**2  + (tg1(ijk,1,1)-tc1(ijk,1,1))**2 + (tb1(ijk,1,1)-td1(ijk,1,1))**2
    enddo
!    call MPI_ALLREDUCE(di1,enstrophy,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
!    enstrophy=half*enstrophy/(nx*ny*nz)

    es = sum(di1*volSimps1)
    call MPI_REDUCE(es,es1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

    ek=zero;ek1=zero
    !KINETIC ENERGY
    di1=zero
    do ijk=1,nvect1
      di1(ijk,1,1) = ux1(ijk,1,1)**2 + uy1(ijk,1,1)**2 + uz1(ijk,1,1)**2
    enddo
!    call MPI_ALLREDUCE(di1,eek,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
!    eek=half*eek/(nx*ny*nz)

    ek = sum(di1*volSimps1)
    call MPI_REDUCE(ek,ek1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

    ds=zero;ds1=zero
    !dissipation
    di1=zero
    call dissipation (ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1)
    ds = sum(di1*volSimps1)
    call MPI_REDUCE(ds,ds1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

    
    if (nrank .eq. 0) then
      FS = 1+1+1+1 !Number of columns
      write(fileformat, '( "(",I4,"(E14.6),A)" )' ) FS
      FS = FS*14+1  !Line width
      open(67,file='./out/statistics',status='unknown',form='formatted',&
      access='direct',recl=FS)
      write(67,fileformat,rec=itime/imodulo2+1) t,& !1
      es1,&                               !3
      ek1,&                             !2
      ds1,&
      NL                                !+1
      close(67)
    end if

   endif
   
  end subroutine postprocessing
    !############################################################################
  subroutine suspended(phi1,vol1,mp1)

    USE decomp_2d_io
    USE MPI

    implicit none

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: vol1
    real(mytype),intent(out) :: mp1(1:nphi)

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: temp1
    real(mytype) :: mp(1:nphi)
    integer :: is,code

    mp=zero; mp1=zero

    do is=1, nphi
       temp1 = phi1(:,:,:,is)*vol1(:,:,:)
       mp(is)= sum(temp1)
    end do

    call MPI_REDUCE(mp,mp1,nphi,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

    return
  end subroutine suspended
  !############################################################################
  subroutine dissipation (ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,diss1)

  USE param
  USE variables
  USE decomp_2d

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,diss1
  real(mytype),dimension(3,3,xsize(1),xsize(2),xsize(3)) :: A
  integer :: k,j,i,m,l

  !INSTANTANEOUS DISSIPATION RATE
  diss1=zero
  A(:,:,:,:,:)=zero
  A(1,1,:,:,:)=ta1(:,:,:) !du/dx=ta1
  A(2,1,:,:,:)=tb1(:,:,:) !dv/dx=tb1
  A(3,1,:,:,:)=tc1(:,:,:) !dw/dx=tc1
  A(1,2,:,:,:)=td1(:,:,:) !du/dy=td1
  A(2,2,:,:,:)=te1(:,:,:) !dv/dy=te1
  A(3,2,:,:,:)=tf1(:,:,:) !dw/dy=tf1
  A(1,3,:,:,:)=tg1(:,:,:) !du/dz=tg1
  A(2,3,:,:,:)=th1(:,:,:) !dv/dz=th1
  A(3,3,:,:,:)=ti1(:,:,:) !dw/dz=ti1

  do k=1,xsize(3)
     do j=1,xsize(2)
        do i=1,xsize(1)
           do m=1,3
              do l=1,3
                 diss1(i,j,k)=diss1(i,j,k)+two*xnu*half*half*(A(l,m,i,j,k)+A(m,l,i,j,k))**two
              enddo
           enddo
        enddo
     enddo
  enddo
  return

  end subroutine dissipation
end module post_processing
