!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

!##############################################################################
    !!
    !!     PROGRAM: visu_vtk
    !! DESCRIPTION: Create vtk/vtr files from XC3D output. Output file is txt 
    !!              always with Float32 precision regardless of the precision used
    !!              to write Xc3d output
    !!
    !!      AUTHOR: 
    !!    MODIFIED: Stefano Rolfo
    !!
!##############################################################################


program visu_vtk

  implicit none

  integer(4) :: nx,ny,nz
  real(4) :: xlx,yly,zlz,dt,dx,dy,dz
  real(8),dimension(:,:,:),allocatable :: input_var_64
  real(4),dimension(:,:,:),allocatable :: input_var_32
  real(4),dimension(:,:,:),allocatable :: output_var_32
  integer :: i,j,k,count,nfil, ivar
  integer(4) :: istret,nclx, ncly, nclz, nvar
  integer(4) :: nfiles, icrfile, file1, filen, ifile
  integer(4) :: dig1, dig2, dig3, dig4, dig5, dig6, dig7
  integer(4) :: step
  real(4),dimension(:),allocatable :: y1
  real(4),dimension(:),allocatable :: y2
  real(4),dimension(:),allocatable :: y3
  real(4) :: pi,xa
  character(len=128) :: input_file
  character(len=128) :: output_file
  integer(4) :: idigits, iprec
  integer    :: length_input_file
  integer(4) :: iparaview
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
  if (idigits < 3 .or. idigits > 7) then
     write(*,*) 'Numbers of digits has to be between 3 and 7 '
     stop
  endif
  write (*,*) 'Writing precision for the XC3D postprocessing (Single = 4, Double = 8)?'
  read (*,*) iprec
  write (*,*) 'Paraview version target (i.e. 5.8 => 58)'
  read (*,*) iparaview
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

  ! Allocate of the in/out variables
  allocate(input_var_64(nx,ny,nz))
  allocate(input_var_32(nx,ny,nz))
  allocate(output_var_32(nx,ny,nz))

  if (nclx==0) dx=xlx/real(nx,4)
  if (nclx==1 .or. nclx==2) dx=xlx/real(nx-1.,4)
  if (ncly==0) dy=yly/real(ny,4)
  if (ncly==1.or.ncly==2) dy=yly/real(ny-1.,4)
  if (nclz==0) dz=zlz/nz
  if (nclz==1.or.nclz==2) dz=zlz/real(nz-1.,4)

  allocate(y1(nx))
  allocate(y2(ny))
  allocate(y3(nz))
  do i=1,nx
     y1(i)=(i-1)*dx
  enddo
  if (istret==1) then
     print *,'We need to read the yp.dat file'
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
  
  step = int(filen/nfiles)

  ! Time loop to write the variables
  do ifile = file1, filen, step
     ! Get the digits for the file name
     if (idigits == 3) then
        dig1 =   ifile/100 + 48
        dig2 = ( ifile - 100*( ifile/100 ) )/10 + 48
        dig3 = ( ifile - 10*( ifile/10 ) )/1 + 48
        chits3(1:3) = char(dig1)//char(dig2)//char(dig3)
        write(*,*) 'File name has structure var_name'//chits3
        output_file = 'XC3D'//chits3
     elseif (idigits == 4) then
        dig1 =   ifile/1000 + 48
        dig2 = ( ifile - 1000*( ifile/1000 ) )/100 + 48
        dig3 = ( ifile - 100*( ifile/100 ) )/10 + 48
        dig4 = ( ifile - 10*( ifile/10 ) )/1 + 48
        chits4(1:4) = char(dig1)//char(dig2)//char(dig3)//char(dig4)
        write(*,*) 'File name has structure var_name'//chits4
        output_file = 'XC3D'//chits4
     elseif (idigits == 5) then
        dig1 =   ifile/10000 + 48
        dig2 = ( ifile - 10000*( ifile/10000 ) )/1000 + 48
        dig3 = ( ifile - 1000*( ifile/1000 ) )/100 + 48
        dig4 = ( ifile - 100*( ifile/100 ) )/10 + 48
        dig5 = ( ifile - 10*( ifile/10 ) )/1 + 48
        chits5(1:5) = char(dig1)//char(dig2)//char(dig3)//char(dig4)//&
                      char(dig5)
        write(*,*) 'File name has structure var_name'//chits5
        output_file = 'XC3D'//chits5
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
        output_file = 'XC3D'//chits6
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
        output_file = 'XC3D'//chits7
     endif
     ! Creation of the file output for the time step
     write(*,*) 'output_file ',output_file
     nfil=41
     if (iparaview >= 58 ) then
        output_file=trim(output_file)//'.vtr'
     else
        output_file=trim(output_file)//'.vtk'
     endif
     write(*,*) 'output_file ',output_file
     open(nfil,file=output_file)
     write(nfil,*)'<VTKFile type="RectilinearGrid" version="0.1"',' byte_order="LittleEndian">'
     write(nfil,*)'  <RectilinearGrid WholeExtent=','"1 ',nx,' 1 ',ny,' 1 ',nz,'">'
     write(nfil,*)'    <Piece Extent=','"1 ',nx,' 1 ',ny,' 1 ',nz,'">'
     write(nfil,*)'      <Coordinates>'
     write(nfil,*)'        <DataArray type="Float32"',' Name="X_COORDINATES"',' NumberOfComponents="1">'
     write(nfil,*) (y1(i),i=1,nx)
     write(nfil,*)'        </DataArray>'
     write(nfil,*)'        <DataArray type="Float32"',' Name="Y_COORDINATES"',' NumberOfComponents="1">'
     write(nfil,*) (y2(j),j=1,ny)
     write(nfil,*)'        </DataArray>'
     write(nfil,*)'        <DataArray type="Float32"',' Name="Z_COORDINATES"',' NumberOfComponents="1">'
     write(nfil,*) (y3(k),k=1,nz)
     write(nfil,*)'        </DataArray>'
     write(nfil,*)'      </Coordinates>'
     ! Now we can add the variables
     write(nfil,*)'      <PointData Scalars="scalar">'
     do ivar = 1, nvar
        if (idigits == 3) then
           input_file = trim(var_name(ivar))//'-'//chits3//'.bin'
        elseif (idigits == 3) then
           input_file = trim(var_name(ivar))//'-'//chits4//'.bin'
        else
           input_file = trim(var_name(ivar))//'-'//chits5//'.bin'
        endif
        length_input_file=len_trim(input_file)

        write(*,*) 'length_input_file',length_input_file
        write(*,*) 'File to open ', input_file
        open(10,file=input_file,form='unformatted',&
             access='direct', recl=8, status='old')
        count = 1
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 if (iprec == 8) then
                    read(10,rec=count) input_var_64(i,j,k)
                 else
                    read(10,rec=count) input_var_32(i,j,k)
                 endif
                 count = count + 1
              enddo
           enddo
        enddo
        close(10)

        output_var_32=0.
        if (iprec == 8) then
           output_var_32=real(input_var_64,4)
        else
           output_var_32=input_var_32
        endif

        write(nfil,*)'        <DataArray Name="',trim(var_name(ivar)),'" type="Float32"', &
                              ' NumberOfComponents="1"',' format="ascii">'
        write(nfil,*) (((output_var_32(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        write(nfil,*)'        </DataArray>'
     enddo
     ! End of writing data
     write(nfil,*)'      </PointData>'
     ! Write the end of the file
     write(nfil,*)'    </Piece>'
     write(nfil,*)'  </RectilinearGrid>'
     write(nfil,*)'</VTKFile>'
     close(nfil)
  enddo
  deallocate(y1)
  deallocate(y2)
  deallocate(y3)
  deallocate(input_var_64)
  deallocate(input_var_32)
  deallocate(output_var_32)

end program visu_vtk
