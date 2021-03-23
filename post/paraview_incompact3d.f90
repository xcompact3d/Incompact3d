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
!################################################################################
    !!
    !!     PROGRAM: visu_vtk
    !! DESCRIPTION: Create Xdmf file for directly read of Xc3d output in paraview  
    !!
    !!      AUTHOR: 
    !!    MODIFIED: Stefano Rolfo
    !!
!################################################################################

program visu_paraview

  implicit none

  integer(4) :: nx,ny,nz
  real(4) :: xlx,yly,zlz,dt,dx,dy,dz
  integer(4) :: nfiles, icrfile, file1, filen, ifile
  integer(4) :: dig1, dig2, dig3, dig4, dig5, dig6, dig7
  real(4), dimension(:), allocatable :: y2(:),y1(:),y3(:)
  integer(4) :: i, j, k, num, aig, ii, nfil,istret,nclx, ncly, nclz
  integer(4) :: idigits, iprec, ivar, nvar
  character(len=128) :: input_file
  character(len=128), dimension(:), allocatable :: var_name
  !IF THE DATA ARE STORED WITH 3 DIGITS, IE UX001,UX002,ETC.
  character(3) :: chits3
  !IF THE DATA ARE STORED WITH 4 DIGITS, IE UX0001,UX0002,ETC.
  character(4) :: chits4
  !IF THE DATA ARE STORED WITH 5 DIGITS, IE UX00001,UX00002,ETC.
  character(5) :: chits5
  !IF THE DATA ARE STORED WITH 5 DIGITS, IE UX00001,UX00002,ETC.
  character(6) :: chits6
  !IF THE DATA ARE STORED WITH 5 DIGITS, IE UX00001,UX00002,ETC.
  character(7) :: chits7

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
  write (*,*) 'Number of digits in the file name (Allowed 3, 4, 5, 6, 7)?'
  read (*,*) idigits
  write (*,*) 'Writing precision for the XC3D output (Single = 4, Double = 8)?'
  read (*,*) iprec
  write (*,*) 'How many variables you would like to export'
  read (*,*) nvar
  allocate(var_name(nvar))
  do ivar = 1, nvar
     write(*,*) 'XC3D variable name (i.e. ux)'
     read(*,*) var_name(ivar)
  enddo
  
  do ivar = 1, nvar
     write(*,*) 'Var list ', trim(var_name(ivar))
  enddo


  if (idigits < 3 .or. idigits > 7) then
     write(*,*) 'Numbers of digits has to be between 3 and 5 '
     stop
  endif

  if (iprec /= 4) then 
     write(*,*) 'Double precision is assumed'
     iprec = 8
  endif
  
  if (nclx==0) dx=xlx/real(nx,4)
  if (nclx==1 .or. nclx==2) dx=xlx/real(nx-1.,4)
  if (ncly==0) dy=yly/real(ny,4)
  if (ncly==1.or.ncly==2) dy=yly/real(ny-1.,4)
  if (nclz==0) dz=zlz/nz
  if (nclz==1.or.nclz==2) dz=zlz/real(nz-1.,4)
  
  dt=1.

  allocate(y1(nx))
  allocate(y2(ny))
  allocate(y3(nz))
  do i=1,nx
     y1(i)=(i-1)*dx
  enddo
  if (istret==1) then
     write(*,*) 'We need to read the yp.dat file'
     open(12,file='yp.dat',form='formatted',status='unknown')
     do j=1,ny
        read(12,*) y2(j)
     enddo
     close(12)
  else
     do j=1,ny
        y2(j)=(j-1)*dy
     enddo
  endif
  do k=1,nz
     y3(k)=(k-1)*dz
  enddo

  ! Write of the XDMF file with the mesh (rectiliner) and variables
  ! Variables simply point to the XC3D files
  nfil=41
  open(nfil,file='visu.xdmf')
  
  ! write of the mesh (3 arrays)
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
  write(nfil,*)'    ',y2(:) 
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
  
  ! Time loop to write the variables
  do ifile = file1, filen
     ! Get the digits for the file name
     if (idigits == 3) then
        dig1 =   ifile/100 + 48
        dig2 = ( ifile - 100*( ifile/100 ) )/10 + 48
        dig3 = ( ifile - 10*( ifile/10 ) )/1 + 48
        chits3(1:3) = char(dig1)//char(dig2)//char(dig3)
        write(*,*) 'File name has structure var_name'//chits3
     elseif (idigits == 4) then
        dig1 =   ifile/1000 + 48
        dig2 = ( ifile - 1000*( ifile/1000 ) )/100 + 48
        dig3 = ( ifile - 100*( ifile/100 ) )/10 + 48
        dig4 = ( ifile - 10*( ifile/10 ) )/1 + 48
        chits4(1:4) = char(dig1)//char(dig2)//char(dig3)//char(dig4)
        write(*,*) 'File name has structure var_name'//chits4
     elseif (idigits == 5) then 
        dig1 =   ifile/10000 + 48
        dig2 = ( ifile - 10000*( ifile/10000 ) )/1000 + 48
        dig3 = ( ifile - 1000*( ifile/1000 ) )/100 + 48
        dig4 = ( ifile - 100*( ifile/100 ) )/10 + 48
        dig5 = ( ifile - 10*( ifile/10 ) )/1 + 48
        chits5(1:5) = char(dig1)//char(dig2)//char(dig3)//char(dig4)//&
                      char(dig5)
        write(*,*) 'File name has structure var_name'//chits5
     elseif (idigits == 6) then 
        dig1 =   ifile/100000 + 48
        dig2 = ( ifile - 100000*( ifile/100000 ) )/10000 + 48
        dig3 = ( ifile - 10000*( ifile/10000 ) )/1000 + 48
        dig4 = ( ifile - 1000*( ifile/1000 ) )/100 + 48
        dig5 = ( ifile - 100*( ifile/100 ) )/10 + 48
        dig6 = ( ifile - 10*( ifile/10 ) )/1 + 48
        chits6(1:6) = char(dig1)//char(dig2)//char(dig3)//char(dig4)//&
                      char(dig5)//char(dig6)
        write(*,*) 'File name has structure var_name'//chits6
     else
        dig1 =   ifile/1000000 + 48
        dig2 = ( ifile - 1000000*( ifile/1000000 ) )/100000 + 48
        dig3 = ( ifile - 100000*( ifile/100000 ) )/10000 + 48
        dig4 = ( ifile - 10000*( ifile/10000 ) )/1000 + 48
        dig5 = ( ifile - 1000*( ifile/1000 ) )/100 + 48
        dig6 = ( ifile - 100*( ifile/100 ) )/10 + 48
        dig7 = ( ifile - 10*( ifile/10 ) )/1 + 48
        chits7(1:7) = char(dig1)//char(dig2)//char(dig3)//char(dig4)//&
                      char(dig5)//char(dig6)//char(dig7)
        write(*,*) 'File name has structure var_name'//chits7
     endif
    !
     write(nfil,'(/)')
     if (idigits == 3) then
        write(nfil,*)'        <Grid Name="'//chits3//'" GridType="Uniform">'
     elseif (idigits == 4) then 
        write(nfil,*)'        <Grid Name="'//chits4//'" GridType="Uniform">'
     else 
        write(nfil,*)'        <Grid Name="'//chits5//'" GridType="Uniform">'
     endif
     write(nfil,*)'            <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
     write(nfil,*)'            <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
     do ivar=1, nvar
        write(nfil,*)'            <Attribute Name="',trim(var_name(ivar)),'" Center="Node">'
        write(nfil,*)'               <DataItem Format="Binary" '
        if (iprec == 8) then
           write(nfil,*)'                DataType="Float" Precision="8" Endian="little"'
        else
           write(nfil,*)'                DataType="Float" Precision="4" Endian="little"'
        endif
        write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
        if (idigits == 3) then
           input_file = trim(var_name(ivar))//chits3
        elseif (idigits == 3) then
           input_file = trim(var_name(ivar))//chits4
        else
           input_file = trim(var_name(ivar))//chits5
        endif
        write(nfil,*)'                  ',trim(input_file)
        write(nfil,*)'               </DataItem>'
        write(nfil,*)'            </Attribute>'
     enddo
     write(nfil,*)'        </Grid>'
  enddo
  write(nfil,'(/)')
  write(nfil,*)'    </Grid>'
  write(nfil,*)'</Domain>'
  write(nfil,'(A7)')'</Xdmf>'
  close(nfil)

end program visu_paraview
