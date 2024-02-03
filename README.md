
[![DOI](https://zenodo.org/badge/127266756.svg)](https://zenodo.org/badge/latestdoi/127266756)

## The Xcompact3d code

Xcompact3d is a Fortran-based framework of high-order finite-difference flow solvers 
dedicated to the study of turbulent flows using high fidelity modelling such as 
Direct and Large Eddy Simulations (DNS/LES), for which the largest turbulent scales are simulated.
Xcompact3d can combine the versatility of industrial codes with the accuracy of spectral codes by using 
the Immerse Boundary Method (IBM) to simulate comples geometries, while retaining high order accuracy. 
Its user-friendliness, simplicity, versatility, accuracy, scalability, portability and efficiency 
makes it an attractive tool for the Computational Fluid Dynamics community.

Xcompact3d is currently able to solve the incompressible and low-Mach number variable density 
Navier-Stokes equations up to a sixth-order compact finite-difference schemes 
with a spectral-like accuracy on a monobloc Cartesian mesh.  
It was initially designed in France in the mid-90's for serial processors and later ported to HPC systems. 
It can now be used efficiently on hundreds of thousands CPU cores to investigate turbulence 
and heat transfer problems thanks to the open-source library 2DECOMP&FFT, 
which is a Fortran-based 2D pencil/1D slabs decomposition framework to support building 
large-scale parallel applications on distributed memory systems using MPI. 
The library has a distributed Fast Fourier Transform module as well as I/O capabilities.

Fractional time stepping is using for the time advancement, 
as well as solving a Poisson equation for the incompressible flew. 
The Poisson's equation is fully solved in spectral space via the use of relevant 3D Fast Fourier transforms (FFTs),
allowing the use of any kind of boundary conditions for the velocity field. 
Using the concept of the modified wavenumber (to allow for operations in the spectral space 
to have the same accuracy as if they were performed in the physical space), 
the divergence free condition is ensured up to machine accuracy. 
The pressure field is staggered from the velocity field by half a mesh point 
to avoid spurious oscillations created by the implicit finite-difference schemes. 
The modelling of a fixed or moving solid body inside the computational domain is performed 
with a customised Immersed Boundary Method. 
It is based on a direct forcing term in the Navier-Stokes equations to ensure a no-slip boundary condition 
at the wall of the solid body while imposing non-zero velocities inside the solid body 
to avoid discontinuities on the velocity field. 
This customised IBM, fully compatible with the 2D domain decomposition 
and with a possible mesh refinement at the wall, 
is based on a 1D expansion of the velocity field from fluid regions into solid regions 
using Lagrange polynomials or spline reconstructions. 
In order to reach high Reynolds numbers in a context of LES, 
it is possible to customise the coefficients of the second derivative schemes (used for the viscous term) 
to add extra numerical dissipation in the simulation as a substitute of the missing dissipation 
from the small turbulent scales that are not resolved. 


### External Resources

- [**Twitter**](https://twitter.com/incompact3d)

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
Please have a look at [HOWTO.md](HOWTO.md)

### Releases

Releases are available via git tags, also on the main Github, or via the [Zenodo DOI](https://zenodo.org/badge/latestdoi/127266756) (see top of README page on Github).

