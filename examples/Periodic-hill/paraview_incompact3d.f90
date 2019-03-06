program visu_paraview

  implicit none

  integer(4) :: nx,ny,nz
  real(4) :: xlx,yly,zlz,dt,dx,dy,dz
  integer(4) :: nfiles, icrfile, file1, filen, ifile, dig1, dig2, dig3, dig4
  real(4), allocatable :: yp(:),y1(:),y3(:)
  integer(4) :: i, j, k, num, aig, ii, nfil,istret,nclx, ncly, nclz

!IF THE DATA ARE STORED WITH 3 DIGITS, IE UX001,UX002,ETC.
  character(3) :: chits
!IF THE DATA ARE STORED WITH 4 DIGITS, IE UX0001,UX0002,ETC.
!  character(4) :: chits

  write (*,*) 'nx, ny, nz   - Incompact3D'
  read (*,*) nx, ny, nz
  write (*,*) 'xlx, yly, zlz   - Incompact3D'
  read (*,*) xlx, yly, zlz
  write (*,*) 'nclx, ncly, nclz   - Incompact3D'
  read (*,*) nclx, ncly, nclz
  write (*,*) 'n files, first file, last file'
  read (*,*) nfiles,file1, filen
  write (*,*) 'Stretching in the y direction (Y=1/N=0)?'
  read (*,*) istret


  if (nclx==0) dx=xlx/nx
  if (nclx==1 .or. nclx==2) dx=xlx/(nx-1.)
  if (ncly==0) dy=yly/ny
  if (ncly==1.or.ncly==2) dy   =yly/(ny-1.)
  if (nclz==0) dz=zlz/nz
  if (nclz==1.or.nclz==2) dz=zlz/(nz-1.)
  dt=1.

  allocate(y1(nx))
  allocate(yp(ny))
  allocate(y3(nz))
  do i=1,nx
     y1(i)=(i-1)*dx
  enddo
  if (istret==1) then
     print *,'We need to read the yp.dat file'
     open(12,file='yp.dat',form='formatted',status='unknown')
     do j=1,ny
        read(12,*) yp(j)
     enddo
     close(12)
  else
     do j=1,ny
        yp(j)=(j-1)*dy
     enddo
  endif
  do k=1,nz
     y3(k)=(k-1)*dz
  enddo


  nfil=41
  open(nfil,file='visu.xdmf')

  write(nfil,'(A22)')'<?xml version="1.0" ?>'
  write(nfil,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  write(nfil,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
  write(nfil,*)'<Domain>'
  write(nfil,*)'    <Topology name="topo" TopologyType="3DRectMesh"'
  write(nfil,*)'        Dimensions="',nz,ny,nx,'">'
  write(nfil,*)'    </Topology>'
  write(nfil,*)'    <Geometry name="geo" Type="VXVYVZ">'
  write(nfil,*)'    <DataItem Dimensions="',nx,'" NumberType="Float" Precision="4" Format="XML">'
  write(nfil,*)'    ',y1(:) 
  write(nfil,*)'    </DataItem>'
  write(nfil,*)'    <DataItem Dimensions="',ny,'" NumberType="Float" Precision="4" Format="XML">'
  write(nfil,*)'    ',yp(:) 
  write(nfil,*)'    </DataItem>'
  write(nfil,*)'    <DataItem Dimensions="',nz,'" NumberType="Float" Precision="4" Format="XML">'
  write(nfil,*)'    ',y3(:) 
  write(nfil,*)'    </DataItem>'
  write(nfil,*)'    </Geometry>'
  write(nfil,'(/)')
  write(nfil,*)'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
  write(nfil,*)'        <Time TimeType="HyperSlab">'
  write(nfil,*)'            <DataItem Format="XML" NumberType="Float" Dimensions="3">'
  write(nfil,*)'           <!--Start, Stride, Count-->'
  write(nfil,*)'            0.0',dt
  write(nfil,*)'            </DataItem>'
  write(nfil,*)'        </Time>'

  do ifile = file1, filen

!IF THE DATA ARE STORED WITH 4 DIGITS, IE UX0001,UX0002,ETC.
  !   dig1 =   ifile/1000 + 48
  !   dig2 = ( ifile - 1000*( ifile/1000 ) )/100 + 48
  !   dig3 = ( ifile - 100*( ifile/100 ) )/10 + 48
  !   dig4 = ( ifile - 10*( ifile/10 ) )/1 + 48
  !   chits(1:4) = char(dig1)//char(dig2)//char(dig3)//char(dig4)    

!IF THE DATA ARE STORED WITH 3 DIGITS, IE UX001,UX002,ETC.
    dig1 =   ifile/100 + 48
    dig2 = ( ifile - 100*( ifile/100 ) )/10 + 48
    dig3 = ( ifile - 10*( ifile/10 ) )/1 + 48
    chits(1:3) = char(dig1)//char(dig2)//char(dig3)

     write(*,*) ifile, 'file'//chits

     write(nfil,'(/)')
     write(nfil,*)'        <Grid Name="'//chits//'" GridType="Uniform">'
     write(nfil,*)'            <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
     write(nfil,*)'            <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
!SINGLE PRECISION-->Precision=4
!DOUBLE PRECISION-->Precision=8
     write(nfil,*)'            <Attribute Name="ux" Center="Node">'
     write(nfil,*)'               <DataItem Format="Binary" '
     write(nfil,*)'                DataType="Float" Precision="8" Endian="little"'
     write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
     write(nfil,*)'                  ux'//chits
     write(nfil,*)'               </DataItem>'
     write(nfil,*)'            </Attribute>'

!it is possible to add as much field as you want for example uy

    write(nfil,*)'            <Attribute Name="uy" Center="Node">'
    write(nfil,*)'               <DataItem Format="Binary" '
    write(nfil,*)'                DataType="Float" Precision="8" Endian="little"'
    write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
    write(nfil,*)'                  uy'//chits
    write(nfil,*)'               </DataItem>'
    write(nfil,*)'            </Attribute>'

    write(nfil,*)'            <Attribute Name="eps" Center="Node">'
    write(nfil,*)'               <DataItem Format="Binary" '
    write(nfil,*)'                DataType="Float" Precision="8" Endian="little"'
    write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
    write(nfil,*)'                  ep'//chits
    write(nfil,*)'               </DataItem>'
    write(nfil,*)'            </Attribute>'


     write(nfil,*)'        </Grid>'

  enddo
  write(nfil,'(/)')
  write(nfil,*)'    </Grid>'
  write(nfil,*)'</Domain>'
  write(nfil,'(A7)')'</Xdmf>'
  close(nfil)

end program visu_paraview
