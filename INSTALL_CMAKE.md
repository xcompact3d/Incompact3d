XCompact3d installation documentation
=====================================

## Source Download and Compilation

Please have a look at [INSTALL.md](INSTALL.md) for the instructions on how to download, build and install 
the code. 

Xcompact3d sources can be acquired by cloning the git repository: 

   git clone https://github.com/xcompact3d/Incompact3d

If you are behind a firewall, you may need to use the `https` protocol instead of the `git` protocol:

   git config --global url."https://".insteadOf git@

Be sure to also configure your system to use the appropriate proxy settings, 
e.g. by setting the `https_proxy` and `http_proxy` variables.

**Note:** The compiling process
The build system for Xcompact3d is based on CMake. 
It is good practice to directly point to the 
MPI Fortran wrapper that you would like to use to guarantee consistency between Fortran compiler and MPI. 
This can be done by setting the default Fortran environmental variable 
```
export FC=my_mpif90
```
To generate the build system run 
```
cmake -S $path_to_sources -B $path_to_build_directory -DOPTION1 -DOPTION2 ... 
```
for example 
```
cmake -S . -B build  
```
By defult the build system will also download 2DECOMP&FFT and perform the build install using the
Generic FFT backend. 
If the directory does not exist it will be generated and it will contain the configuration files.
The configuration can be further
edited by using the `ccmake` utility as
```
ccmake $path_to_build_directory
```
To compile the sources 
```
cmake --build $path_to_build_directory -j <nproc> 
```
By defult the Taylor-Green-Vortex case is also activated and can performed with
```
ctest --test-dir $path_to_build_directory
```
The full test suite, which includes 14 differents tests, can be activated with the variable
`BUILD_TESTING_FULL` as 
```
cmake --build $path_to_build_directory -DBUILD_TESTING_FULL=ON 
```
or by using `ccmake.
`
The installation directory will cointain:
* The *bin* directory with two execulables: **xcompact3d** for the main execution of the code 
* The *example* directory with few example of input *.i3d* files for **Xcompact3d**
* The *lib* directory with the archive for the **decomp2d** library
* The directories *Testing* and *RunTests* with the logs of CTest and the results respectively. 

### CTest
To test your installation you can also type in the terminal from your *build* directory
```
ctest --test-dir $path_to_build_directory 
```
Four tests are performed:
* Taylor Green Vortex (TGV)
* Atmospheric Boundary layer (ABL) in neutral conditions (new set-up)
* Atmospheric Boundary layer (ABL) in neutral conditions (old set-up)
* Atmospheric Boundary layer (ABL) in convective conditions (old set-up)
* Atmospheric Boundary layer (ABL) in stable conditions (old set-up)
* Differentially heated cavity
* Turbulent Channel Flow with X as streamwise direction
* Turbulent Channel Flow with Z as streamwise direction
* Flow around a circular cylinder
* Flow around a moving circular cylinder
* Lock exchange
* Turbulent Boundary Layer (TBL)
* Wind Turbine

The simulations results are located under 
```
$ /path/to/build/RunTests
```
and the standard output from the simulations is in 
```
$ /path/to/build/Testing/Temporary/
```

### Build with an already present 2DECOMP&FFT
If different options from the defualt such as a different backend for the FFT rather than the generic or 
IO with ADIOS2 are necessary, 2DECOMP&FFT needs to be pre-installed as described [here](https://github.com/2decomp-fft/2decomp-fft/blob/dev/INSTALL.md).
2DECOMP&FFT installation provides CMake configuration file that can be used to find the installation directory. 
To allow the `find_package` of CMake to work the following variable needs to be set as
```
$ export decomp2d_DIR=/path/to/2decomp/install/opt/lib/decomp2d 
```
Depending on the system *lib* can be *lib64* 

## Known issues




