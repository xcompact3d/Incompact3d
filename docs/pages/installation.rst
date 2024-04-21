============
Installation
============

The latest version of Xcompact3d only supports `cmake` based builds. The only requirements is a Fortran 90-compatible compiler and a working `MPI` library.


Minimum requirements
====================
**cmake version 3.20 and above**

**a recent modern Fortran compiler (ie, gfortran 9 and above**

Source Download and Compilation
===============================

Xcompact3d sources can be acquired by cloning the git repository: 

``git clone https://github.com/xcompact3d/Incompact3d``

If you are behind a firewall, you may need to use the `https` protocol instead of the `git` protocol:

``git config --global url."https://".insteadOf git@``

Be sure to also configure your system to use the appropriate proxy settings, 
e.g. by setting the `https_proxy` and `http_proxy` variables.


The compiling process
=====================

The build system for Xcompact3d is based on CMake. It is good practice to directly point to the 
MPI Fortran wrapper that you would like to use to guarantee consistency between Fortran compiler and MPI. 
This can be done by setting the default Fortran environmental variable 

``export FC=my_mpif90``

To generate the build system run 

``cmake -S $path_to_sources -B $path_to_build_directory -DOPTION1 -DOPTION2 ... ``

for example 

``cmake -S . -B build``

By defult the build system will also download the 2DECOMP&FFT library  and perform the build install using the
Generic FFT backend. Version 2.0.3 of that library is the default for Xcompact3d building
and all tests are performed against this specific version.
If the directory does not exist it will be generated and it will contain the configuration files.
The configuration can be further
edited by using the `ccmake` utility as

``ccmake $path_to_build_directory``

To compile the sources 

``cmake --build $path_to_build_directory -j <nproc>``

for example, when in the build directory

``cmake --build . -j 8`` will compile the code using 8 CPU cores.

Appending `-v` will display additional information about the build, such as compiler flags.

The executable file **xcompact3d** is located in the build/bin directory.



Testing
=======
The testing suite for the **xcompact3d** solver is composed by 14 tests as follows 

1. Atmospheric Boundary layer (ABL) in neutral conditions (new set-up)

2. Atmospheric Boundary layer (ABL) in neutral conditions (old set-up)

3. Atmospheric Boundary layer (ABL) in convective conditions (old set-up)

4. Atmospheric Boundary layer (ABL) in stable conditions (old set-up)

5. Differentially heated cavity

6. Turbulent Channel Flow with X as streamwise direction

7. Turbulent Channel Flow with Z as streamwise direction

8. Flow around a circular cylinder

9. Flow around a moving circular cylinder

10. Lock exchange

11. Mixing Layer

12. Turbulent Boundary Layer (TBL)

13. Wind Turbine

14. Taylor Green Vortex (TGV)

By default only the  Taylor Green Vortex case is activated, while the full 
testing suite needs to be enable by using the `BUILD_TESTING_FULL` flag as 

``cmake --build $path_to_build_directory -DBUILD_TESTING_FULL=ON``

or by using `ccmake`.

The tests are performed using `CTest` as  

``ctest --test-dir $path_to_build_directory``

Every test is performed in a dedicated working directory that is located under the following path 

``/path/to/build/RunTests``

All standard outputs from all test runs are collated under the file

``/path/to/build/Testing/Temporary/LastTest.log``

together with additional files detailing additional informations such as 
the elapse time for the different tests and the eventual failed cases. 



***Note***
Some of the alternative options for FFT and IO backends required additional input
* For MKL FFT the location of the MKL libraires needs to be passed to the configure as 
for the 2DECOMP&FFT installation with 

``export MKL_DIR=${MKLROOT}/lib/cmake/mkl``

* For ADIOS the installation directory needs to be passes to the configure as

``cmake -S . -B ./build -DIO_BACKEND=adios2 -Dadios2_DIR=/path/to/adios2/install/lib/cmake/adios2``

Both steps are necessary for correct linking of the target **xcompact3d** with the libraries 

Known issues
===============
Some issues with ADIOS2.

The tests performed under `CTest` rely on the `CMake` ability to properly find the MPI executable *mpirun*. 

The build system will try to enforce consistency between the MPI Fortran used and the MPI executable, 
for the first iteration of the configure step. 

In case no MPI executable is not found or correct please modify manually the `MPIEXEC_EXECUTABLE` by using 

``cmake -S . -B build -DMPIEXEC_EXECUTABLE=/correct/path/to/mpirun``

or by using 

``ccmake $path_to_build_directory``

