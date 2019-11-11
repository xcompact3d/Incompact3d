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

SUBROUTINE paraview()

  USE variables
  USE param
  USE decomp_2d

  IMPLICIT NONE

  INTEGER :: i,k,ifile,is,nfil=48
  real(mytype) :: xp(nx),zp(nz)
  CHARACTER(4) :: chits

9999 FORMAT(I4.4)

#ifdef DEBG
  IF (nrank .EQ. 0) PRINT *,'# paraview start'
#endif

  PRINT *,'Writing XDMF files'

  do i=1,nx
     xp(i) = real(i-1,mytype)*dx
  enddo
  do k=1,nz
     zp(k) = real(k-1,mytype)*dz
  enddo

#ifdef DEBG
  IF (nrank .EQ. 0) PRINT *,'# paraview snapshots.xdmf start'
#endif

  OPEN(nfil,file='1_snapshots.xdmf')
  WRITE(nfil,'(A22)')'<?xml version="1.0" ?>'
  WRITE(nfil,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  WRITE(nfil,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
  WRITE(nfil,*)'<Domain>'
  if (istret.ne.0) then
     write(nfil,*)'    <Topology name="topo" TopologyType="3DRectMesh"'
     write(nfil,*)'        Dimensions="',nz,ny,nx,'">'
     write(nfil,*)'    </Topology>'
     write(nfil,*)'    <Geometry name="geo" Type="VXVYVZ">'
     write(nfil,*)'        <DataItem Dimensions="',nx,'" NumberType="Float" Precision="4" Format="XML">'
     write(nfil,*)'        ',xp(:) 
     write(nfil,*)'        </DataItem>'
     write(nfil,*)'        <DataItem Dimensions="',ny,'" NumberType="Float" Precision="4" Format="XML">'
     write(nfil,*)'        ',yp(:) 
     write(nfil,*)'        </DataItem>'
     write(nfil,*)'        <DataItem Dimensions="',nz,'" NumberType="Float" Precision="4" Format="XML">'
     write(nfil,*)'        ',zp(:) 
     write(nfil,*)'        </DataItem>'
     write(nfil,*)'    </Geometry>'
  else
     WRITE(nfil,*)'    <Topology name="topo" TopologyType="3DCoRectMesh"'
     WRITE(nfil,*)'        Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
     WRITE(nfil,*)'    </Topology>'
     WRITE(nfil,*)'    <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
     WRITE(nfil,*)'        <!-- Origin -->'
     WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
     WRITE(nfil,*)'        0.0 0.0 0.0'
     WRITE(nfil,*)'        </DataItem>'
     WRITE(nfil,*)'        <!-- DxDyDz -->'
     WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
     WRITE(nfil,*)'        ',dz*nvisu,dy*nvisu,dx*nvisu
     WRITE(nfil,*)'        </DataItem>'
     WRITE(nfil,*)'    </Geometry>'
  endif
  WRITE(nfil,'(/)')
  WRITE(nfil,*)'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
  WRITE(nfil,*)'        <Time TimeType="HyperSlab">'
  WRITE(nfil,*)'            <DataItem Format="XML" NumberType="Float" Dimensions="3">'
  WRITE(nfil,*)'            <!--Start, Stride, Count-->'
  WRITE(nfil,*)'            0.0',REAL(ioutput*dt,4)
  WRITE(nfil,*)'            </DataItem>'
  WRITE(nfil,*)'        </Time>'

  DO ifile = 0, INT(ilast/ioutput)
     WRITE(chits, 9999) ifile

     WRITE(nfil,'(/)')
     WRITE(nfil,*)'        <Grid Name="'//chits//'" GridType="Uniform">'
     WRITE(nfil,*)'            <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
     WRITE(nfil,*)'            <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
     IF (save_ibm.ne.0) THEN
        WRITE(nfil,*)'            <Attribute Name="ibm" Center="Node">'
        WRITE(nfil,*)'               <DataItem Format="Binary" '
        WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
        WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
        IF (save_ibm.eq.1) WRITE(nfil,*)'                  ./data/ibm0000'
        IF (save_ibm.eq.2) WRITE(nfil,*)'                  ./data/ibm'//chits
        WRITE(nfil,*)'               </DataItem>'
        WRITE(nfil,*)'            </Attribute>'
     END IF
     IF (iscalar.EQ.1.AND.save_phi.EQ.1) THEN
        DO is=1, numscalar
           WRITE(nfil,*)'            <Attribute Name="phi'//CHAR(48+is)//'" Center="Node">'
           WRITE(nfil,*)'               <DataItem Format="Binary" '
           WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
           WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
           WRITE(nfil,*)'                  ./data/phi'//CHAR(48+is)//chits
           WRITE(nfil,*)'               </DataItem>'
           WRITE(nfil,*)'            </Attribute>'
        END DO
     ENDIF
     IF (save_qc.EQ.1) THEN
        WRITE(nfil,*)'            <Attribute Name="q-criterion" Center="Node">'
        WRITE(nfil,*)'               <DataItem Format="Binary" '
        WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
        WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
        WRITE(nfil,*)'                  ./data/qc'//chits
        WRITE(nfil,*)'               </DataItem>'
        WRITE(nfil,*)'            </Attribute>'
     ENDIF
     IF (save_pc.EQ.1) THEN
        WRITE(nfil,*)'            <Attribute Name="p-criterion" Center="Node">'
        WRITE(nfil,*)'               <DataItem Format="Binary" '
        WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
        WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
        WRITE(nfil,*)'                  ./data/pc'//chits
        WRITE(nfil,*)'               </DataItem>'
        WRITE(nfil,*)'            </Attribute>'
     ENDIF
     IF (save_V.EQ.1) THEN
        WRITE(nfil,*)'            <Attribute Name="V" Center="Node">'
        WRITE(nfil,*)'               <DataItem Format="Binary" '
        WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
        WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
        WRITE(nfil,*)'                  ./data/V'//chits
        WRITE(nfil,*)'               </DataItem>'
        WRITE(nfil,*)'            </Attribute>'
     ENDIF
     IF (save_ux.EQ.1) THEN
        WRITE(nfil,*)'            <Attribute Name="ux" Center="Node">'
        WRITE(nfil,*)'               <DataItem Format="Binary" '
        WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
        WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
        WRITE(nfil,*)'                  ./data/ux'//chits
        WRITE(nfil,*)'               </DataItem>'
        WRITE(nfil,*)'            </Attribute>'
     ENDIF
     IF (save_uy.EQ.1) THEN
        WRITE(nfil,*)'            <Attribute Name="uy" Center="Node">'
        WRITE(nfil,*)'               <DataItem Format="Binary" '
        WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
        WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
        WRITE(nfil,*)'                  ./data/uy'//chits
        WRITE(nfil,*)'               </DataItem>'
        WRITE(nfil,*)'            </Attribute>'
     ENDIF
     IF (save_uz.EQ.1) THEN
        WRITE(nfil,*)'            <Attribute Name="uz" Center="Node">'
        WRITE(nfil,*)'               <DataItem Format="Binary" '
        WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
        WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
        WRITE(nfil,*)'                  ./data/uz'//chits
        WRITE(nfil,*)'               </DataItem>'
        WRITE(nfil,*)'            </Attribute>'
     ENDIF
     IF (save_pre.EQ.1) THEN
        WRITE(nfil,*)'            <Attribute Name="pre" Center="Node">'
        WRITE(nfil,*)'               <DataItem Format="Binary" '
        WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
        WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
        WRITE(nfil,*)'                  ./data/pre'//chits
        WRITE(nfil,*)'               </DataItem>'
        WRITE(nfil,*)'            </Attribute>'
     ENDIF
     IF (save_w.EQ.1) THEN
        WRITE(nfil,*)'            <Attribute Name="w" Center="Node">'
        WRITE(nfil,*)'               <DataItem Format="Binary" '
        WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
        WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
        WRITE(nfil,*)'                  ./data/w'//chits
        WRITE(nfil,*)'               </DataItem>'
        WRITE(nfil,*)'            </Attribute>'
     ENDIF
     IF (save_w1.EQ.1) THEN
        WRITE(nfil,*)'            <Attribute Name="w1" Center="Node">'
        WRITE(nfil,*)'               <DataItem Format="Binary" '
        WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
        WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
        WRITE(nfil,*)'                  ./data/w1'//chits
        WRITE(nfil,*)'               </DataItem>'
        WRITE(nfil,*)'            </Attribute>'
     ENDIF
     IF (save_w2.EQ.1) THEN
        WRITE(nfil,*)'            <Attribute Name="w2" Center="Node">'
        WRITE(nfil,*)'               <DataItem Format="Binary" '
        WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
        WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
        WRITE(nfil,*)'                  ./data/w2'//chits
        WRITE(nfil,*)'               </DataItem>'
        WRITE(nfil,*)'            </Attribute>'
     ENDIF
     IF (save_w3.EQ.1) THEN
        WRITE(nfil,*)'            <Attribute Name="w3" Center="Node">'
        WRITE(nfil,*)'               <DataItem Format="Binary" '
        WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
        WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
        WRITE(nfil,*)'                  ./data/w3'//chits
        WRITE(nfil,*)'               </DataItem>'
        WRITE(nfil,*)'            </Attribute>'
     ENDIF
     IF (save_dudx.EQ.1) THEN
        WRITE(nfil,*)'            <Attribute Name="du/dx" Center="Node">'
        WRITE(nfil,*)'               <DataItem Format="Binary" '
        WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
        WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
        WRITE(nfil,*)'                  ./data/dudx'//chits
        WRITE(nfil,*)'               </DataItem>'
        WRITE(nfil,*)'            </Attribute>'
     ENDIF
     IF (save_dudy.EQ.1) THEN
        WRITE(nfil,*)'            <Attribute Name="du/dy" Center="Node">'
        WRITE(nfil,*)'               <DataItem Format="Binary" '
        WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
        WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
        WRITE(nfil,*)'                  ./data/dudy'//chits
        WRITE(nfil,*)'               </DataItem>'
        WRITE(nfil,*)'            </Attribute>'
     ENDIF
     IF (save_dudz.EQ.1) THEN
        WRITE(nfil,*)'            <Attribute Name="du/dz" Center="Node">'
        WRITE(nfil,*)'               <DataItem Format="Binary" '
        WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
        WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
        WRITE(nfil,*)'                  ./data/dudz'//chits
        WRITE(nfil,*)'               </DataItem>'
        WRITE(nfil,*)'            </Attribute>'
     ENDIF
     IF (save_dvdx.EQ.1) THEN
        WRITE(nfil,*)'            <Attribute Name="dv/dx" Center="Node">'
        WRITE(nfil,*)'               <DataItem Format="Binary" '
        WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
        WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
        WRITE(nfil,*)'                  ./data/dvdx'//chits
        WRITE(nfil,*)'               </DataItem>'
        WRITE(nfil,*)'            </Attribute>'
     ENDIF
     IF (save_dvdy.EQ.1) THEN
        WRITE(nfil,*)'            <Attribute Name="dv/dy" Center="Node">'
        WRITE(nfil,*)'               <DataItem Format="Binary" '
        WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
        WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
        WRITE(nfil,*)'                  ./data/dvdy'//chits
        WRITE(nfil,*)'               </DataItem>'
        WRITE(nfil,*)'            </Attribute>'
     ENDIF
     IF (save_dvdz.EQ.1) THEN
        WRITE(nfil,*)'            <Attribute Name="dv/dz" Center="Node">'
        WRITE(nfil,*)'               <DataItem Format="Binary" '
        WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
        WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
        WRITE(nfil,*)'                  ./data/dvdz'//chits
        WRITE(nfil,*)'               </DataItem>'
        WRITE(nfil,*)'            </Attribute>'
     ENDIF
     IF (save_dwdx.EQ.1) THEN
        WRITE(nfil,*)'            <Attribute Name="dw/dx" Center="Node">'
        WRITE(nfil,*)'               <DataItem Format="Binary" '
        WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
        WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
        WRITE(nfil,*)'                  ./data/dwdx'//chits
        WRITE(nfil,*)'               </DataItem>'
        WRITE(nfil,*)'            </Attribute>'
     ENDIF
     IF (save_dwdy.EQ.1) THEN
        WRITE(nfil,*)'            <Attribute Name="dw/dy" Center="Node">'
        WRITE(nfil,*)'               <DataItem Format="Binary" '
        WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
        WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
        WRITE(nfil,*)'                  ./data/dwdy'//chits
        WRITE(nfil,*)'               </DataItem>'
        WRITE(nfil,*)'            </Attribute>'
     ENDIF
     IF (save_dwdz.EQ.1) THEN
        WRITE(nfil,*)'            <Attribute Name="dw/dz" Center="Node">'
        WRITE(nfil,*)'               <DataItem Format="Binary" '
        WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
        WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
        WRITE(nfil,*)'                  ./data/dwdz'//chits
        WRITE(nfil,*)'               </DataItem>'
        WRITE(nfil,*)'            </Attribute>'
     ENDIF

     IF (iscalar.EQ.1.AND.(save_dphidx.EQ.1.OR.save_dphidy.EQ.1.OR.save_dphidz.EQ.1.)) THEN
        DO is=1, numscalar
           IF (save_dphidx.EQ.1) THEN
              WRITE(nfil,*)'            <Attribute Name="dphi'//CHAR(48+is)//'dx" Center="Node">'
              WRITE(nfil,*)'               <DataItem Format="Binary" '
              WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
              WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
              WRITE(nfil,*)'                  ./data/dphi'//CHAR(48+is)//'dx'//chits
              WRITE(nfil,*)'               </DataItem>'
              WRITE(nfil,*)'            </Attribute>'
           ENDIF
           IF (save_dphidy.EQ.1) THEN
              WRITE(nfil,*)'            <Attribute Name="dphi'//CHAR(48+is)//'dy" Center="Node">'
              WRITE(nfil,*)'               <DataItem Format="Binary" '
              WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
              WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
              WRITE(nfil,*)'                  ./data/dphi'//CHAR(48+is)//'dy'//chits
              WRITE(nfil,*)'               </DataItem>'
              WRITE(nfil,*)'            </Attribute>'
           ENDIF
           IF (save_dphidz.EQ.1) THEN
              WRITE(nfil,*)'            <Attribute Name="dphi'//CHAR(48+is)//'dz" Center="Node">'
              WRITE(nfil,*)'               <DataItem Format="Binary" '
              WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
              WRITE(nfil,*)'                Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
              WRITE(nfil,*)'                  ./data/dphi'//CHAR(48+is)//'dz'//chits
              WRITE(nfil,*)'               </DataItem>'
              WRITE(nfil,*)'            </Attribute>'
           ENDIF
        END DO
     ENDIF

     WRITE(nfil,*)'        </Grid>'
  ENDDO
  WRITE(nfil,'(/)')
  WRITE(nfil,*)'    </Grid>'
  WRITE(nfil,*)'</Domain>'
  WRITE(nfil,'(A7)')'</Xdmf>'
  CLOSE(nfil)

#ifdef DEBG
  IF (nrank .EQ. 0) PRINT *,'# paraview snapshots.xdmf done'
  IF (nrank .EQ. 0) PRINT *,'# paraview planes.xdmf start'
#endif

  !###############################################################################
  IF(save_phim.EQ.1.OR.save_uxm.EQ.1.OR.save_uym.EQ.1.OR.save_uzm.EQ.1.OR.save_prem.EQ.1) THEN
     OPEN(nfil,file='2_mean-planes-z.xdmf')
     WRITE(nfil,'(A22)')'<?xml version="1.0" ?>'
     WRITE(nfil,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
     WRITE(nfil,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
     WRITE(nfil,*)'<Domain>'
     if (istret.ne.0) then
        write(nfil,*)'    <Topology name="topo" TopologyType="3DRectMesh"'
        write(nfil,*)'        Dimensions="',1,ny,nx,'">'
        write(nfil,*)'    </Topology>'
        write(nfil,*)'    <Geometry name="geo" Type="VXVYVZ">'
        write(nfil,*)'        <DataItem Dimensions="',nx,'" NumberType="Float" Precision="4" Format="XML">'
        write(nfil,*)'        ',xp(:) 
        write(nfil,*)'        </DataItem>'
        write(nfil,*)'        <DataItem Dimensions="',ny,'" NumberType="Float" Precision="4" Format="XML">'
        write(nfil,*)'        ',yp(:) 
        write(nfil,*)'        </DataItem>'
        write(nfil,*)'        <DataItem Dimensions="1" NumberType="Float" Precision="4" Format="XML">'
        write(nfil,*)'        0.0'
        write(nfil,*)'        </DataItem>'
        write(nfil,*)'    </Geometry>'
     else
        WRITE(nfil,*)'    <Topology name="topo" TopologyType="3DCoRectMesh"'
        WRITE(nfil,*)'        Dimensions="',1,ny,nx,'">'
        WRITE(nfil,*)'    </Topology>'
        WRITE(nfil,*)'    <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
        WRITE(nfil,*)'        <!-- Origin -->'
        WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
        WRITE(nfil,*)'        0.0 0.0 0.0'
        WRITE(nfil,*)'        </DataItem>'
        WRITE(nfil,*)'        <!-- DxDyDz -->'
        WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
        WRITE(nfil,*)'        ',0.0,dy,dx
        WRITE(nfil,*)'        </DataItem>'
        WRITE(nfil,*)'    </Geometry>'
     endif
     WRITE(nfil,'(/)')
     WRITE(nfil,*)'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
     WRITE(nfil,*)'        <Time TimeType="HyperSlab">'
     WRITE(nfil,*)'            <DataItem Format="XML" NumberType="Float" Dimensions="3">'
     WRITE(nfil,*)'           <!--Start, Stride, Count-->'
     WRITE(nfil,*)'            0.0',REAL(ioutput*dt,4)
     WRITE(nfil,*)'            </DataItem>'
     WRITE(nfil,*)'        </Time>'

     DO ifile = 0, INT(ilast/ioutput)
        WRITE(chits, 9999) ifile

        WRITE(nfil,'(/)')
        WRITE(nfil,*)'        <Grid Name="'//chits//'" GridType="Uniform">'
        WRITE(nfil,*)'            <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
        WRITE(nfil,*)'            <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'

        IF (iscalar.EQ.1.AND.save_phim.EQ.1) THEN
           DO is=1, numscalar
              WRITE(nfil,*)'            <Attribute Name="phim'//CHAR(48+is)//'" Center="Node">'
              WRITE(nfil,*)'               <DataItem Format="Binary" '
              WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
              WRITE(nfil,*)'                Dimensions="',1,ny,nx,'">'
              WRITE(nfil,*)'                  ./data/phim'//CHAR(48+is)//chits
              WRITE(nfil,*)'               </DataItem>'
              WRITE(nfil,*)'            </Attribute>'
           END DO
        ENDIF
        IF (save_uxm.EQ.1) THEN
           WRITE(nfil,*)'            <Attribute Name="uxm" Center="Node">'
           WRITE(nfil,*)'               <DataItem Format="Binary" '
           WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
           WRITE(nfil,*)'                Dimensions="',1,ny,nx,'">'
           WRITE(nfil,*)'                  ./data/uxm'//chits
           WRITE(nfil,*)'               </DataItem>'
           WRITE(nfil,*)'            </Attribute>'
        ENDIF
        IF (save_uym.EQ.1) THEN
           WRITE(nfil,*)'            <Attribute Name="uym" Center="Node">'
           WRITE(nfil,*)'               <DataItem Format="Binary" '
           WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
           WRITE(nfil,*)'                Dimensions="',1,ny,nx,'">'
           WRITE(nfil,*)'                  ./data/uym'//chits
           WRITE(nfil,*)'               </DataItem>'
           WRITE(nfil,*)'            </Attribute>'
        ENDIF
        IF (save_uzm.EQ.1) THEN
           WRITE(nfil,*)'            <Attribute Name="uzm" Center="Node">'
           WRITE(nfil,*)'               <DataItem Format="Binary" '
           WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
           WRITE(nfil,*)'                Dimensions="',1,ny,nx,'">'
           WRITE(nfil,*)'                  ./data/uzm'//chits
           WRITE(nfil,*)'               </DataItem>'
           WRITE(nfil,*)'            </Attribute>'
        ENDIF
        IF (save_prem.EQ.1) THEN
           WRITE(nfil,*)'            <Attribute Name="prem" Center="Node">'
           WRITE(nfil,*)'               <DataItem Format="Binary" '
           WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
           WRITE(nfil,*)'                Dimensions="',1,ny,nx,'">'
           WRITE(nfil,*)'                  ./data/prem'//chits
           WRITE(nfil,*)'               </DataItem>'
           WRITE(nfil,*)'            </Attribute>'
        ENDIF
        WRITE(nfil,*)'        </Grid>'
     ENDDO

     WRITE(nfil,'(/)')
     WRITE(nfil,*)'    </Grid>'
     WRITE(nfil,*)'</Domain>'
     WRITE(nfil,'(A7)')'</Xdmf>'
     CLOSE(nfil)

#ifdef DEBG
     IF (nrank .EQ. 0) PRINT *,'# paraview planes.xdmf done'
     IF (nrank .EQ. 0) PRINT *,'# paraview stats_sum.xdmf start'
#endif

  ENDIF

  !###############################################################################
  IF(save_dmap.EQ.1.OR.save_utmap.EQ.1) THEN

     OPEN(nfil,file='3_maps.xdmf')
     WRITE(nfil,'(A22)')'<?xml version="1.0" ?>'
     WRITE(nfil,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
     WRITE(nfil,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
     WRITE(nfil,*)'<Domain>'
     WRITE(nfil,*)'    <Topology name="topo" TopologyType="3DCoRectMesh"'
     WRITE(nfil,*)'        Dimensions="',nz,1,nx,'">'
     WRITE(nfil,*)'    </Topology>'
     WRITE(nfil,*)'    <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
     WRITE(nfil,*)'        <!-- Origin -->'
     WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
     WRITE(nfil,*)'        0.0 0.0 0.0'
     WRITE(nfil,*)'        </DataItem>'
     WRITE(nfil,*)'        <!-- DxDyDz -->'
     WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
     WRITE(nfil,*)'        ',dz,0.0,dx
     WRITE(nfil,*)'        </DataItem>'
     WRITE(nfil,*)'    </Geometry>'
     WRITE(nfil,'(/)')
     WRITE(nfil,*)'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
     WRITE(nfil,*)'        <Time TimeType="HyperSlab">'
     WRITE(nfil,*)'            <DataItem Format="XML" NumberType="Float" Dimensions="3">'
     WRITE(nfil,*)'           <!--Start, Stride, Count-->'
     WRITE(nfil,*)'            0.0',REAL(ioutput*dt,4)
     WRITE(nfil,*)'            </DataItem>'
     WRITE(nfil,*)'        </Time>'

     DO ifile = 0, INT(ilast/ioutput)
        WRITE(chits, 9999) ifile

        WRITE(nfil,'(/)')
        WRITE(nfil,*)'        <Grid Name="'//chits//'" GridType="Uniform">'
        WRITE(nfil,*)'            <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
        WRITE(nfil,*)'            <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'

        IF (iscalar.EQ.1) THEN
           DO is=1, numscalar

              WRITE(nfil,*)'            <Attribute Name="phif'//CHAR(48+is)//'" Center="Node">'
              WRITE(nfil,*)'               <DataItem Format="Binary" '
              WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
              WRITE(nfil,*)'                Dimensions="',nz,1,nx,'">'
              WRITE(nfil,*)'                  ./data/phif'//CHAR(48+is)//chits
              WRITE(nfil,*)'               </DataItem>'
              WRITE(nfil,*)'            </Attribute>'

              IF (save_dmap.EQ.1) THEN
                 WRITE(nfil,*)'            <Attribute Name="dep'//CHAR(48+is)//'" Center="Node">'
                 WRITE(nfil,*)'               <DataItem Format="Binary" '
                 WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
                 WRITE(nfil,*)'                Dimensions="',nz,1,nx,'">'
                 WRITE(nfil,*)'                  ./out/dep'//CHAR(48+is)//chits
                 WRITE(nfil,*)'               </DataItem>'
                 WRITE(nfil,*)'            </Attribute>'
              ENDIF

           ENDDO
        ENDIF

        IF (save_utmap.EQ.1) THEN
           WRITE(nfil,*)'            <Attribute Name="utmap" Center="Node">'
           WRITE(nfil,*)'               <DataItem Format="Binary" '
           WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
           WRITE(nfil,*)'                Dimensions="',nz,1,nx,'">'
           WRITE(nfil,*)'                  ./data/utmap'//chits
           WRITE(nfil,*)'               </DataItem>'
           WRITE(nfil,*)'            </Attribute>'
        ENDIF

        WRITE(nfil,*)'        </Grid>'
     ENDDO
     WRITE(nfil,'(/)')
     WRITE(nfil,*)'    </Grid>'
     WRITE(nfil,*)'</Domain>'
     WRITE(nfil,'(A7)')'</Xdmf>'
     CLOSE(nfil)

  ENDIF

  !###############################################################################
#ifdef STATS

  OPEN(nfil,file='4_stats_sum.xdmf')
  WRITE(nfil,'(A22)')'<?xml version="1.0" ?>'
  WRITE(nfil,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  WRITE(nfil,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
  WRITE(nfil,*)'<Domain>'
  WRITE(nfil,*)'    <Topology name="topo" TopologyType="3DCoRectMesh"'
  WRITE(nfil,*)'        Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'    </Topology>'
  WRITE(nfil,*)'    <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
  WRITE(nfil,*)'        <!-- Origin -->'
  WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
  WRITE(nfil,*)'        0.0 0.0 0.0'
  WRITE(nfil,*)'        </DataItem>'
  WRITE(nfil,*)'        <!-- DxDyDz -->'
  WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
  WRITE(nfil,*)'        ',dz*nstat,dy*nstat,dx*nstat
  WRITE(nfil,*)'        </DataItem>'
  WRITE(nfil,*)'    </Geometry>'
  WRITE(nfil,'(/)')
  WRITE(nfil,*)'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
  WRITE(nfil,*)'        <Time TimeType="HyperSlab">'
  WRITE(nfil,*)'            <DataItem Format="XML" NumberType="Float" Dimensions="3">'
  WRITE(nfil,*)'           <!--Start, Stride, Count-->'
  WRITE(nfil,*)'            0.0',ilast
  WRITE(nfil,*)'            </DataItem>'
  WRITE(nfil,*)'        </Time>'

  WRITE(nfil,'(/)')
  WRITE(nfil,*)'        <Grid Name="" GridType="Uniform">'
  WRITE(nfil,*)'            <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
  WRITE(nfil,*)'            <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'

  IF (save_pre.EQ.1) THEN
     WRITE(nfil,*)'            <Attribute Name="pre.sum" Center="Node">'
     WRITE(nfil,*)'               <DataItem Format="Binary" '
     WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
     WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
     WRITE(nfil,*)'                  ./data/pre.sum'
     WRITE(nfil,*)'               </DataItem>'
     WRITE(nfil,*)'            </Attribute>'
  ENDIF

  WRITE(nfil,*)'            <Attribute Name="c1.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./data/c1.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="c1c1.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./data/c1c1.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="uc1.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./data/uc.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="vc1.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./data/vc1.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="wc1.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./data/wc1.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="u1.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./data/u1.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="v1.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./data/v1.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="w1.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./data/w1.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="u2.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./data/u2.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="v2.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./data/v2.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="w2.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./data/w2.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="uv.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./data/uv.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="uw.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./data/uv.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="vw.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./data/vw.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="k.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./data/k.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  IF (save_dudx.EQ.0) THEN
     WRITE(nfil,*)'            <Attribute Name="dudx.sum" Center="Node">'
     WRITE(nfil,*)'               <DataItem Format="Binary" '
     WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
     WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
     WRITE(nfil,*)'                  ./data/dudx.sum'
     WRITE(nfil,*)'               </DataItem>'
     WRITE(nfil,*)'            </Attribute>'
  ENDIF

  WRITE(nfil,*)'            <Attribute Name="u3.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./data/u3.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="v3.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./data/v3.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="w3.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./data/w3.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="u4.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./data/u4.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="v4.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./data/v4.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="w4.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./data/w4.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'        </Grid>'
  WRITE(nfil,'(/)')
  WRITE(nfil,*)'    </Grid>'
  WRITE(nfil,*)'</Domain>'
  WRITE(nfil,'(A7)')'</Xdmf>'
  CLOSE(nfil)

#ifdef DEBG
  IF (nrank .EQ. 0) PRINT *,'# paraview stats_sum.xdmf done'
  IF (nrank .EQ. 0) PRINT *,'# paraview stats_mean.xdmf start'
#endif

  !###############################################################################

  OPEN(nfil,file='5_stats_mean_convergence.xdmf')
  WRITE(nfil,'(A22)')'<?xml version="1.0" ?>'
  WRITE(nfil,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  WRITE(nfil,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
  WRITE(nfil,*)'<Domain>'
  WRITE(nfil,*)'    <Topology name="topo" TopologyType="3DCoRectMesh"'
  WRITE(nfil,*)'        Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'    </Topology>'
  WRITE(nfil,*)'    <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
  WRITE(nfil,*)'        <!-- Origin -->'
  WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
  WRITE(nfil,*)'        0.0 0.0 0.0'
  WRITE(nfil,*)'        </DataItem>'
  WRITE(nfil,*)'        <!-- DxDyDz -->'
  WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
  WRITE(nfil,*)'        ',dz*nstat,dy*nstat,dx*nstat
  WRITE(nfil,*)'        </DataItem>'
  WRITE(nfil,*)'    </Geometry>'
  WRITE(nfil,'(/)')
  WRITE(nfil,*)'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
  WRITE(nfil,*)'        <Time TimeType="HyperSlab">'
  WRITE(nfil,*)'            <DataItem Format="XML" NumberType="Float" Dimensions="3">'
  WRITE(nfil,*)'           <!--Start, Stride, Count-->'
  WRITE(nfil,*)'            0.0',REAL(ntik*dt,4)
  WRITE(nfil,*)'            </DataItem>'
  WRITE(nfil,*)'        </Time>'

  DO ifile = 0, INT(ilast/ioutput)
     WRITE(chits, 9999) ifile

     WRITE(nfil,'(/)')
     WRITE(nfil,*)'        <Grid Name="'//chits//'" GridType="Uniform">'
     WRITE(nfil,*)'            <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
     WRITE(nfil,*)'            <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'

     WRITE(nfil,*)'            <Attribute Name="ux" Center="Node">'
     WRITE(nfil,*)'               <DataItem Format="Binary" '
     WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
     WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
     WRITE(nfil,*)'                  ./mean/uxmean'//chits
     WRITE(nfil,*)'               </DataItem>'
     WRITE(nfil,*)'            </Attribute>'

     WRITE(nfil,*)'        </Grid>'
  ENDDO
  WRITE(nfil,'(/)')
  WRITE(nfil,*)'    </Grid>'
  WRITE(nfil,*)'</Domain>'
  WRITE(nfil,'(A7)')'</Xdmf>'
  CLOSE(nfil)

#ifdef DEBG
  IF (nrank .EQ. 0) PRINT *,'# paraview stats_mean.xdmf end'
  IF (nrank .EQ. 0) PRINT *,'# paraview end'
#endif

#endif
END SUBROUTINE paraview
