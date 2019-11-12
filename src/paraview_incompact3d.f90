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

program visu_paraview

  implicit none

  integer(4) :: nx,ny,nz,ialm,iadm,ivirt,jles,ibuoyancy,isnapshot,iprec,itripping
  real(4) :: xlx,yly,zlz,dt,dx,dy,dz
  integer(4) :: nfiles, icrfile, file1, filen, ifile, isnap, snap1,snapn,dig1, dig2, dig3, dig4
  real(4), allocatable :: yp(:),y1(:),y3(:)
  integer(4) :: i, j, k, num, aig, ii, nfil,istret,nclx, ncly, nclz, meanfil, snapfil
  integer(4) :: imin,imax,jmin,jmax,kmin,kmax,snapx,snapy,snapz
  integer :: ErrFlag, nargin, FNLength, status, DecInd,dynhypvisc
  logical :: back
  character(len=80) :: InputFN, FNBase

  character(3) :: chits
  NAMELIST/PostProcess/nx,ny,nz,xlx,yly,zlz,nclx,ncly,nclz,istret,nfiles,file1,filen,ialm,iadm,ivirt,jles,ibuoyancy,dynhypvisc,&
       itripping,isnapshot,imin,imax,jmin,jmax,kmin,kmax,snap1,snapn,iprec
  !==========================================================================
  ! Handle Input file
  nargin=command_argument_count()
  if (nargin <1) then
     write(6,*) 'Please call the program with the name of the input file on the command line Ex. visualise input.pprc'
     stop
  endif

  call get_command_argument(1,InputFN,FNLength,status)
  back=.true.
  FNBase=inputFN((index(InputFN,'/',back)+1):len(InputFN))
  DecInd=index(FNBase,'.',back)
  if (DecInd >1) then
     FNBase=FNBase(1:(DecInd-1))
  end if
  !===========================================================================

  iprec=8 
  open(10,file=InputFN) 
  read(10,nml=PostProcess)
  close(10) 

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
  if (istret>1) then
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
  write(nfil,*)'    <Topology name="topo" TopologyType="3DCoRectMesh"'
  write(nfil,*)'        Dimensions="',nz,ny,nx,'">'
  write(nfil,*)'        Dimensions="',nz,ny,nx,'">'
  write(nfil,*)'    </Topology>'
  write(nfil,*)'    <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
  write(nfil,*)'        <!-- Origin -->'
  write(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
  write(nfil,*)'        0.0 0.0 0.0'
  write(nfil,*)'        </DataItem>'
  write(nfil,*)'        <!-- DxDyDz -->'
  write(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
  write(nfil,*)'        ',dz,dy,dx
  write(nfil,*)'        </DataItem>'
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
     !     dig1 =  ifile/1000 + 48
     !     dig2 = (ifile - 1000*( ifile/1000 ) )/100 + 48
     !     dig3 = (ifile - 100*( ifile/100 ) )/10 + 48
     !     dig4 = (ifile - 10*( ifile/10 ) )/1 + 48
     !     chits(1:4) = char(dig1)//char(dig2)//char(dig3)//char(dig4)    

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
     write(nfil,*)'                DataType="Float" Precision="',iprec,'" Endian="little"'
     write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
     write(nfil,*)'                  ux'//chits
     write(nfil,*)'               </DataItem>'
     write(nfil,*)'            </Attribute>'

     write(nfil,*)'            <Attribute Name="uy" Center="Node">'
     write(nfil,*)'               <DataItem Format="Binary" '
     write(nfil,*)'                DataType="Float" Precision="',iprec,'" Endian="little"'
     write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
     write(nfil,*)'                  uy'//chits
     write(nfil,*)'               </DataItem>'
     write(nfil,*)'            </Attribute>'

     write(nfil,*)'            <Attribute Name="uz" Center="Node">'
     write(nfil,*)'               <DataItem Format="Binary" '
     write(nfil,*)'                DataType="Float" Precision="',iprec,'" Endian="little"'
     write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
     write(nfil,*)'                  uz'//chits
     write(nfil,*)'               </DataItem>'
     write(nfil,*)'            </Attribute>'

     write(nfil,*)'            <Attribute Name="pp" Center="Node">'
     write(nfil,*)'               <DataItem Format="Binary" '
     write(nfil,*)'                DataType="Float" Precision="',iprec,'" Endian="little"'
     write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
     write(nfil,*)'                  pp'//chits
     write(nfil,*)'               </DataItem>'
     write(nfil,*)'            </Attribute>'

     write(nfil,*)'            <Attribute Name="vort" Center="Node">'
     write(nfil,*)'               <DataItem Format="Binary" '
     write(nfil,*)'                DataType="Float" Precision="',iprec,'" Endian="little"'
     write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
     write(nfil,*)'                  vort'//chits
     write(nfil,*)'               </DataItem>'
     write(nfil,*)'            </Attribute>'

     if(iadm==1) then
        write(nfil,*)'            <Attribute Name="gammadisc" Center="Node">'
        write(nfil,*)'               <DataItem Format="Binary" '
        write(nfil,*)'                DataType="Float" Precision="',iprec,'" Endian="little"'
        write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
        write(nfil,*)'                  gammadisc'//chits
        write(nfil,*)'               </DataItem>'
        write(nfil,*)'            </Attribute>'
     endif

     if(ialm==1) then
        write(nfil,*)'            <Attribute Name="Ftx" Center="Node">'
        write(nfil,*)'               <DataItem Format="Binary" '
        write(nfil,*)'                DataType="Float" Precision="',iprec,'" Endian="little"'
        write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
        write(nfil,*)'                  Ftx'//chits
        write(nfil,*)'               </DataItem>'
        write(nfil,*)'            </Attribute>'

        write(nfil,*)'            <Attribute Name="Fty" Center="Node">'
        write(nfil,*)'               <DataItem Format="Binary" '
        write(nfil,*)'                DataType="Float" Precision="',iprec,'" Endian="little"'
        write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
        write(nfil,*)'                  Fty'//chits
        write(nfil,*)'               </DataItem>'
        write(nfil,*)'            </Attribute>'

        write(nfil,*)'            <Attribute Name="Ftz" Center="Node">'
        write(nfil,*)'               <DataItem Format="Binary" '
        write(nfil,*)'                DataType="Float" Precision="',iprec,'" Endian="little"'
        write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
        write(nfil,*)'                  Ftz'//chits
        write(nfil,*)'               </DataItem>'
        write(nfil,*)'            </Attribute>'
     endif

     if(itripping==1.or.itripping==2) then
        write(nfil,*)'            <Attribute Name="Ftrip" Center="Node">'
        write(nfil,*)'               <DataItem Format="Binary" '
        write(nfil,*)'                DataType="Float" Precision="',iprec,'" Endian="little"'
        write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
        write(nfil,*)'                  Ftrip'//chits
        write(nfil,*)'               </DataItem>'
        write(nfil,*)'            </Attribute>'
     endif

     if(jles.eq.1.and.dynhypvisc.eq.1) then
        write(nfil,*)'            <Attribute Name="dynvisc" Center="Node">'
        write(nfil,*)'               <DataItem Format="Binary" '
        write(nfil,*)'                DataType="Float" Precision="',iprec,'" Endian="little"'
        write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
        write(nfil,*)'                  dynvisc'//chits
        write(nfil,*)'               </DataItem>'
        write(nfil,*)'            </Attribute>'
     endif

     if(ibuoyancy==1) then
        write(nfil,*)'            <Attribute Name="temp" Center="Node">'
        write(nfil,*)'               <DataItem Format="Binary" '
        write(nfil,*)'                DataType="Float" Precision="',iprec,'" Endian="little"'
        write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
        write(nfil,*)'                  temp'//chits
        write(nfil,*)'               </DataItem>'
        write(nfil,*)'            </Attribute>'
     endif

     write(nfil,*)'            <Attribute Name="Q" Center="Node">'
     write(nfil,*)'               <DataItem Format="Binary" '
     write(nfil,*)'                DataType="Float" Precision="',iprec,'" Endian="little"'
     write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
     write(nfil,*)'                  Q'//chits
     write(nfil,*)'               </DataItem>'
     write(nfil,*)'            </Attribute>'

     write(nfil,*)'        </Grid>'

  enddo
  write(nfil,'(/)')
  write(nfil,*)'    </Grid>'
  write(nfil,*)'</Domain>'
  write(nfil,'(A7)')'</Xdmf>'
  close(nfil)



  if (isnapshot==1) then
     snapfil=31
     snapz=kmax-kmin+1
     snapy=jmax-jmin+1
     snapx=imax-imin+1

     open(snapfil,file='visu_snap.xdmf')

     write(snapfil,'(A22)')'<?xml version="1.0" ?>'
     write(snapfil,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
     write(snapfil,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
     write(snapfil,*)'<Domain>'
     write(snapfil,*)'    <Topology name="topo" TopologyType="3DCoRectMesh"'
     write(snapfil,*)'        Dimensions="',snapz,snapy,snapx,'">'
     write(snapfil,*)'        Dimensions="',snapz,snapy,snapx,'">'
     write(snapfil,*)'    </Topology>'
     write(snapfil,*)'    <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
     write(snapfil,*)'        <!-- Origin -->'
     write(snapfil,*)'        <DataItem Format="XML" Dimensions="3">'
     write(snapfil,*)'        0.0 0.0 0.0'
     write(snapfil,*)'        </DataItem>'
     write(snapfil,*)'        <!-- DxDyDz -->'
     write(snapfil,*)'        <DataItem Format="XML" Dimensions="3">'
     write(snapfil,*)'        ',dz,dy,dx
     write(snapfil,*)'        </DataItem>'
     write(snapfil,*)'    </Geometry>'

     write(snapfil,'(/)')
     write(snapfil,*)'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
     write(snapfil,*)'        <Time TimeType="HyperSlab">'
     write(snapfil,*)'            <DataItem Format="XML" NumberType="Float" Dimensions="3">'
     write(snapfil,*)'           <!--Start, Stride, Count-->'
     write(snapfil,*)'            0.0',dt
     write(snapfil,*)'            </DataItem>'
     write(snapfil,*)'        </Time>'

     do isnap = snap1, snapn

        !IF THE DATA ARE STORED WITH 4 DIGITS, IE UX0001,UX0002,ETC.
        !   dig1 =   isnap/1000 + 48
        !   dig2 = ( isnap - 1000*( isnap/1000 ) )/100 + 48
        !   dig3 = ( isnap - 100*( isnap/100 ) )/10 + 48
        !   dig4 = ( isnap - 10*( isnap/10 ) )/1 + 48
        !   chits(1:4) = char(dig1)//char(dig2)//char(dig3)//char(dig4)    

        !IF THE DATA ARE STORED WITH 3 DIGITS, IE UX001,UX002,ETC.
        dig1 =   ifile/100 + 48
        dig2 = ( ifile - 100*( ifile/100 ) )/10 + 48
        dig3 = ( ifile - 10*( ifile/10 ) )/1 + 48
        chits(1:3) = char(dig1)//char(dig2)//char(dig3)

        write(*,*) isnap, 'file'//chits

        write(snapfil,'(/)')
        write(snapfil,*)'        <Grid Name="'//chits//'" GridType="Uniform">'
        write(snapfil,*)'            <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
        write(snapfil,*)'            <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
        !SINGLE PRECISION-->Precision=4
        !DOUBLE PRECISION-->Precision=8
        write(snapfil,*)'            <Attribute Name="snap_ux" Center="Node">'
        write(snapfil,*)'               <DataItem Format="Binary" '
        write(snapfil,*)'                DataType="Float" Precision="8" Endian="little"'
        write(snapfil,*)'                Dimensions="',snapz,snapy,snapx,'">'
        write(snapfil,*)'                  snap_ux'//chits
        write(snapfil,*)'               </DataItem>'
        write(snapfil,*)'            </Attribute>'

        write(snapfil,*)'            <Attribute Name="snap_uy" Center="Node">'
        write(snapfil,*)'               <DataItem Format="Binary" '
        write(snapfil,*)'                DataType="Float" Precision="8" Endian="little"'
        write(snapfil,*)'                Dimensions="',snapz,snapy,snapx,'">'
        write(snapfil,*)'                  snap_uy'//chits
        write(snapfil,*)'               </DataItem>'
        write(snapfil,*)'            </Attribute>'

        write(snapfil,*)'            <Attribute Name="snap_uz" Center="Node">'
        write(snapfil,*)'               <DataItem Format="Binary" '
        write(snapfil,*)'                DataType="Float" Precision="8" Endian="little"'
        write(snapfil,*)'                Dimensions="',snapz,snapy,snapx,'">'
        write(snapfil,*)'                  snap_uz'//chits
        write(snapfil,*)'               </DataItem>'
        write(snapfil,*)'            </Attribute>'

        write(snapfil,*)'        </Grid>'

     enddo
     write(snapfil,'(/)')
     write(snapfil,*)'    </Grid>'
     write(snapfil,*)'</Domain>'
     write(snapfil,'(A7)')'</Xdmf>'
     close(snapfil)

  endif

  meanfil=51 
  open(meanfil,file='visu_mean.xdmf')

  write(meanfil,'(A22)')'<?xml version="1.0" ?>'
  write(meanfil,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  write(meanfil,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
  write(meanfil,*)'<Domain>'
  write(meanfil,*)'    <Topology name="topo" TopologyType="3DRectMesh"'
  write(meanfil,*)'        Dimensions="',nz,ny,nx,'">'
  write(meanfil,*)'    </Topology>'
  write(meanfil,*)'    <Geometry name="geo" Type="VXVYVZ">'
  write(meanfil,*)'    <DataItem Dimensions="',nx,'" NumberType="Float" Precision="8" Format="XML">'
  write(meanfil,*)'    ',y1(:) 
  write(meanfil,*)'    </DataItem>'
  write(meanfil,*)'    <DataItem Dimensions="',ny,'" NumberType="Float" Precision="8" Format="XML">'
  write(meanfil,*)'    ',yp(:) 
  write(meanfil,*)'    </DataItem>'
  write(meanfil,*)'    <DataItem Dimensions="',nz,'" NumberType="Float" Precision="8" Format="XML">'
  write(meanfil,*)'    ',y3(:) 
  write(meanfil,*)'    </DataItem>'
  write(meanfil,*)'    </Geometry>'
  write(meanfil,'(/)')
  write(meanfil,*)'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
  write(meanfil,*)'        <Time TimeType="HyperSlab">'
  write(meanfil,*)'            <DataItem Format="XML" NumberType="Float" Dimensions="3">'
  write(meanfil,*)'           <!--Start, Stride, Count-->'
  write(meanfil,*)'            0.0',dt
  write(meanfil,*)'            </DataItem>'
  write(meanfil,*)'        </Time>'


  write(*,*) 'writing mean variables', meanfil

  write(meanfil,'(/)')
  write(meanfil,*)'        <Grid Name="'//chits//'" GridType="Uniform">'
  write(meanfil,*)'            <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
  write(meanfil,*)'            <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'

  write(meanfil,*)'            <Attribute Name="umean" Center="Node">'
  write(meanfil,*)'               <DataItem Format="Binary" '
  write(meanfil,*)'                DataType="Float" Precision="8" Endian="little"'
  write(meanfil,*)'                Dimensions="',nz,ny,nx,'">'
  write(meanfil,*)'                  umean.dat'
  write(meanfil,*)'               </DataItem>'
  write(meanfil,*)'            </Attribute>'

  write(meanfil,*)'            <Attribute Name="vmean" Center="Node">'
  write(meanfil,*)'               <DataItem Format="Binary" '
  write(meanfil,*)'                DataType="Float" Precision="8" Endian="little"'
  write(meanfil,*)'                Dimensions="',nz,ny,nx,'">'
  write(meanfil,*)'                  vmean.dat'
  write(meanfil,*)'               </DataItem>'
  write(meanfil,*)'            </Attribute>'

  write(meanfil,*)'            <Attribute Name="wmean" Center="Node">'
  write(meanfil,*)'               <DataItem Format="Binary" '
  write(meanfil,*)'                DataType="Float" Precision="8" Endian="little"'
  write(meanfil,*)'                Dimensions="',nz,ny,nx,'">'
  write(meanfil,*)'                  wmean.dat'
  write(meanfil,*)'               </DataItem>'
  write(meanfil,*)'            </Attribute>'

  write(meanfil,*)'            <Attribute Name="uumean" Center="Node">'
  write(meanfil,*)'               <DataItem Format="Binary" '
  write(meanfil,*)'                DataType="Float" Precision="8" Endian="little"'
  write(meanfil,*)'                Dimensions="',nz,ny,nx,'">'
  write(meanfil,*)'                  uumean.dat'
  write(meanfil,*)'               </DataItem>'
  write(meanfil,*)'            </Attribute>'

  write(meanfil,*)'            <Attribute Name="vvmean" Center="Node">'
  write(meanfil,*)'               <DataItem Format="Binary" '
  write(meanfil,*)'                DataType="Float" Precision="8" Endian="little"'
  write(meanfil,*)'                Dimensions="',nz,ny,nx,'">'
  write(meanfil,*)'                  vvmean.dat'
  write(meanfil,*)'               </DataItem>'
  write(meanfil,*)'            </Attribute>'

  write(meanfil,*)'            <Attribute Name="wwmean" Center="Node">'
  write(meanfil,*)'               <DataItem Format="Binary" '
  write(meanfil,*)'                DataType="Float" Precision="8" Endian="little"'
  write(meanfil,*)'                Dimensions="',nz,ny,nx,'">'
  write(meanfil,*)'                  wwmean.dat'
  write(meanfil,*)'               </DataItem>'
  write(meanfil,*)'            </Attribute>'

  write(meanfil,*)'            <Attribute Name="uvmean" Center="Node">'
  write(meanfil,*)'               <DataItem Format="Binary" '
  write(meanfil,*)'                DataType="Float" Precision="8" Endian="little"'
  write(meanfil,*)'                Dimensions="',nz,ny,nx,'">'
  write(meanfil,*)'                  uvmean.dat'
  write(meanfil,*)'               </DataItem>'
  write(meanfil,*)'            </Attribute>'

  write(meanfil,*)'            <Attribute Name="uwmean" Center="Node">'
  write(meanfil,*)'               <DataItem Format="Binary" '
  write(meanfil,*)'                DataType="Float" Precision="8" Endian="little"'
  write(meanfil,*)'                Dimensions="',nz,ny,nx,'">'
  write(meanfil,*)'                  uwmean.dat'
  write(meanfil,*)'               </DataItem>'
  write(meanfil,*)'            </Attribute>'

  write(meanfil,*)'            <Attribute Name="vwmean" Center="Node">'
  write(meanfil,*)'               <DataItem Format="Binary" '
  write(meanfil,*)'                DataType="Float" Precision="8" Endian="little"'
  write(meanfil,*)'                Dimensions="',nz,ny,nx,'">'
  write(meanfil,*)'                  vwmean.dat'
  write(meanfil,*)'               </DataItem>'
  write(meanfil,*)'            </Attribute>'
  if(jLES.ge.2) then 
     write(meanfil,*)'            <Attribute Name="tauxymean" Center="Node">'
     write(meanfil,*)'               <DataItem Format="Binary" '
     write(meanfil,*)'                DataType="Float" Precision="8" Endian="little"'
     write(meanfil,*)'                Dimensions="',nz,ny,nx,'">'
     write(meanfil,*)'                  vwmean.dat'
     write(meanfil,*)'               </DataItem>'
     write(meanfil,*)'            </Attribute>'
     write(meanfil,*)'        </Grid>'
  endif
  write(meanfil,'(/)')
  write(meanfil,*)'    </Grid>'
  write(meanfil,*)'</Domain>'
  write(meanfil,'(A7)')'</Xdmf>'
  close(meanfil)

end program visu_paraview
