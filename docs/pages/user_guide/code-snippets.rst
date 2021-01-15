Code Snippets
=============

Fortran
-------

- XDMF writer

.. code-block:: fortran

  ! Compile with:
  ! $gfortran -o paraview paraview_incompact3d.f90
  ! Run with:
  ! ./paraview
  !
  ! Observe the numerical precision:
  !     Precision="8" <- If -DDOUBLE_PREC was used at Xcompact3d Makefile
  !     Precision="4" <- Otherwise
  !
  ! Observe if your files are named with 3 or 4 digits, change this file properly
  !
  ! You can add more 3D fields if necessary, see examples below
  !
  !
  program visu_paraview

    implicit none

    integer(4) :: nx,ny,nz
    real(4) :: xlx,yly,zlz,dt,dx,dy,dz
    integer(4) :: nfiles, icrfile, file1, filen, ifile, dig1, dig2, dig3, dig4
    real(4), allocatable :: yp(:),y1(:),y3(:)
    integer(4) :: i, j, k, num, aig, ii, nfil,istret,nclx, ncly, nclz

    character(3) :: chits

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
      open(12,file='yp1.dat',form='formatted',status='unknown')
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
      ! dig1 =   ifile/1000 + 48
      ! dig2 = ( ifile - 1000*( ifile/1000 ) )/100 + 48
      ! dig3 = ( ifile - 100*( ifile/100 ) )/10 + 48
      ! dig4 = ( ifile - 10*( ifile/10 ) )/1 + 48
      ! chits(1:4) = char(dig1)//char(dig2)//char(dig3)//char(dig4)

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
      write(nfil,*)'                DataType="Float" Precision="8" Endian="little"  Seek="0"'
      write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
      write(nfil,*)'                  ux'//chits
      write(nfil,*)'               </DataItem>'
      write(nfil,*)'            </Attribute>'

      !it is possible to add many fields just copying and renaming this block:
      write(nfil,*)'            <Attribute Name="uy" Center="Node">'
      write(nfil,*)'               <DataItem Format="Binary" '
      write(nfil,*)'                DataType="Float" Precision="8" Endian="little"  Seek="0"'
      write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
      write(nfil,*)'                  uy'//chits
      write(nfil,*)'               </DataItem>'
      write(nfil,*)'            </Attribute>'
      !end of block

      write(nfil,*)'            <Attribute Name="uz" Center="Node">'
      write(nfil,*)'               <DataItem Format="Binary" '
      write(nfil,*)'                DataType="Float" Precision="8" Endian="little"  Seek="0"'
      write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
      write(nfil,*)'                  uz'//chits
      write(nfil,*)'               </DataItem>'
      write(nfil,*)'            </Attribute>'

      write(nfil,*)'            <Attribute Name="ibm" Center="Node">'
      write(nfil,*)'               <DataItem Format="Binary" '
      write(nfil,*)'                DataType="Float" Precision="8" Endian="little"  Seek="0"'
      write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
      write(nfil,*)'                  ibm000'!//chits !Just file 000 for ibm
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

Python
------

- Reading the `.i3d` file

.. code-block:: python

  def i3d_to_dict(filename='input.i3d'):
      '''
      This function reads the .i3d file from Xcompact3d and
      returns it into a Python dictionary
      '''

      f = open(filename)

      dict_outer = {}

      for line in f:
          # Remove comments
          line = line.partition('!')[0].replace(' ', '')
          # Remove spaces
          line = " ".join(line.split())

          if line == '':  # Cycle if line is empty
              continue

          # Beginning of a new group
          if line[0] == '&':
              key = line[1:]
              dict_inner = {}
              continue

          # End of the group
          if line.lower() == '/end':
              dict_outer[key] = dict_inner
              continue

          # Get variable's name and value
          param = line.partition('=')[0]
          value = line.partition('=')[-1]

          try:
              # Converting from string according to datatype
              if value[0] == "'" and value[-1] == "'":  # String
                  value = value[1:-1]
              elif value.lower() == '.false.':  # Bool
                  value = False
              elif value.lower() == '.true.':  # Bool
                  value = True
              elif "." in value:  # Float
                  value = float(value)
              else:  # Int
                  value = int(value)
          except:
              print(f"Can't convert {param} : {value}")
              continue

          if "(" in param and ")" == param[-1]:  # Param is a list
              param = param.split('(')[0]
              if param not in dict_inner:
                  dict_inner[param] = []
              dict_inner[param].append(value)
          else:  # Not a list
              dict_inner[param] = value

      f.close()

      return dict_outer

- Reading binary fields

.. code-block:: python

  import numpy as np

  def x3d_readfield(filename, shape, avg=None, dtype=np.float64):
      '''
      This functions reads a binary field and returns a numpy array.
      Args:
          filename: File name.
          shape: Tuple with the array size ordered following nx, ny and nz.
                 It is possible to read 3d fields and planes.
          avg: Compute np.average(axis = avg) if desired.
          dtype: Data type for the Numpy array, use np.float64 if Xcompact
                 was compiled with -DDOUBLE_PREC, np.float32 otherwise.
      Returns:
          A Numpy array with shape (shape)
      Example:
          prm = i3d_to_dict()
          nx = prm['BasicParam']['nx']
          ny = prm['BasicParam']['ny']
          nz = prm['BasicParam']['nz']
          x3d_readfield('some_3d_field', (nx, ny, nz))
          x3d_readfield('some_3d_field', avg = -1, (nx, ny, nz))
          x3d_readfield('some_xy_plane', (nx, ny))
          x3d_readfield('some_xz_plane', (nx, nz))
          x3d_readfield('some_yz_plane', (ny, nz))
      '''
      if avg == None:
          return np.fromfile(filename, dtype=dtype).reshape(shape, order='F')
      else:
          return np.average(np.fromfile(filename, dtype=dtype).reshape(shape,
                                                                       order='F'),
                            axis=avg)


  def x3d_readfield_all(target, shape, avg=None, dtype=np.float64):
      '''
      This functions reads all binary fields that match the pattern and returns
      them in a stacked numpy array.

      Args:
          Target: Pathname pattern.
          shape: Tuple with the array size ordered following nx, ny and nz.
                 It is possible to read 3d fields and planes.
          avg: Compute np.average(axis = avg) if desired.
          dtype: Data type for the Numpy array, use np.float64 if Xcompact
                 was compiled with -DDOUBLE_PREC, np.float32 otherwise.

      Returns:
          A array of shape (shape, nfiles)

      Example:
          x3d_readfield_all('./data/ux*', (nx, ny, nz))
          x3d_readfield_all('./data/uy????', (nx, ny, nz), avg = -1)
      '''
      #
      filenames = sorted(glob.glob(target))
      #
      fields = [
          x3d_readfield(file, shape, avg, dtype=dtype) for file in filenames
      ]
      #
      return np.stack(fields, axis=-1)

- Writing binary fields

.. code-block:: python

  def x3d_writefield(array, filename):
      '''
      This function writes a binary field from a numpy array.
      '''
      # Note that a transposition is necessary if the array's
      # shape is ordered as (nx, ny, nz).
      # It works for planes as well, if array's shape is
      # (nx, ny), (nx, nz) or (ny, nz).
      array.T.tofile(filename)
