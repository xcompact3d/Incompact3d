Xcompact3d installation documentation
=====================================

## Source Download and Compilation

Xcompact3d sources can be acquired by cloning the git repository: 
```
$ git clone https://github.com/xcompact3d/Incompact3d
```
If you are behind a firewall, you may need to use the `https` protocol instead of the `git` protocol:
```
$ git config --global url."https://".insteadOf git@
```
Be sure to also configure your system to use the appropriate proxy settings, 
e.g. by setting the `https_proxy` and `http_proxy` variables.

## Minimum requirements

**cmake version 3.20 and above**

**a recent modern Fortran compiler (ie, gfortran 9 and above**

### The compiling process

The build system for Xcompact3d is based on CMake. 
It is good practice to directly point to the 
MPI Fortran wrapper that you would like to use to guarantee consistency 
between Fortran compiler and MPI. 
This can be done by setting the default Fortran environmental variable 
```
$ export FC=my_mpif90
```
To generate the build system run 
```
$ cmake -S $path_to_sources -B $path_to_build_directory -DOPTION1 -DOPTION2 ... 
```
for example 
```
$ cmake -S . -B build  
```
By defult the build system will also download 2DECOMP&FFT 
and perform the build install using the
Generic FFT backend. Version 2.0.3 is the default for Xcompact3d building
and all tests are performed against this specific version.
If the directory does not exist it will be generated and it will contain the configuration files.
The configuration can be further
edited by using the `ccmake` utility as
```
$ ccmake $path_to_build_directory
```
To compile the sources 
```
$ cmake --build $path_to_build_directory -j <nproc>
```
appending `-v` will display additional information about the build, such as compiler flags.

After building the library can be tested using `CTest`
for the available options. 
Finally the code can be installed using 
```
$ cmake --install $path_to_build_directory
```
By default the installation directory is located under 
```
$ $path_to_build_directory/opt
```
and it will cointain: 
* The *bin* directory with the execulables **xcompact3d** for the main execution of the code;
* The *example* directory with some examples of input *.i3d* files.

To change the default location the `CMAKE_INSTALL_PREFIX` can be modified using 
```
$ cmake --build $path_to_build_directory -DCMAKE_INSTALL_PREFIX=$path_to_my_opt
```
or via the `ccmake` interface. 

## Examples and Testing
Several input *.i3d* files are available under the [examples](examples) folder.
The input files contain the full set-up, from flow development to statistic collection,
for several canonical test cases.    

### Testing
The testing suite for the **xcompact3d** solver is composed by 13 tests. 
More details are given in [here](tests/README.md)


## Build with an already present 2DECOMP&FFT
If different options from the default 
(i.e. Generic FFT backend and double precision) are necessary, 
2DECOMP&FFT needs to be pre-installed as described [here](https://github.com/2decomp-fft/2decomp-fft/blob/dev/INSTALL.md).
Alternative available options are: 
* FFTW or MKL for the FFT backend engine;
* ADIOS2 instead of MPI-IO for the IO operations;
* SINGLE precision for the build.

2DECOMP&FFT installation provides CMake configuration file that can be used to find the installation directory. 
To allow the `find_package` of `CMake` to work the following variable needs to be set as
```
$ export decomp2d_DIR=/path/to/2decomp/install/opt/lib/decomp2d 
```
Depending on the system *lib* can be *lib64*.

**Note**

Some of the alternative options for FFT and IO backends required additional input.

* For MKL FFT the location of the MKL libraires needs to be passed to the configure as 
for the 2DECOMP&FFT installation with 
```
$ export MKL_DIR=${MKLROOT}/lib/cmake/mkl
```

* For ADIOS the installation directory needs to be passes to the configure as
```
$ cmake -S . -B ./build -DIO_BACKEND=adios2 -Dadios2_DIR=/path/to/adios2/install/lib/cmake/adios2
```

Both steps are necessary for correct linking of the target **xcompact3d** with the libraries 

## Known issues
The tests performed under `CTest` rely on the `CMake` ability to properly find the MPI executable *mpirun*. 
The build system will try to enforce consistency between the MPI Fortran used and the MPI executable, 
for the first iteration of the configure step. 
In case no MPI executable is not found or correct please modify manually the `MPIEXEC_EXECUTABLE` by using 
```
$ cmake -S . -B build -DMPIEXEC_EXECUTABLE=/correct/path/to/mpirun
```
or by using 
```
$ ccmake $path_to_build_directory
```

