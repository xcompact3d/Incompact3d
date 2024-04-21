
[![DOI](https://zenodo.org/badge/127266756.svg)](https://zenodo.org/badge/latestdoi/127266756)

## The Xcompact3d code

Xcompact3d is a Fortran-based framework of high-order finite-difference flow solvers 
dedicated to the study of turbulent flows using high fidelity modelling such as 
Direct and Large Eddy Simulations (DNS/LES), for which the largest turbulent scales are simulated.
Xcompact3d can combine the versatility of industrial codes with the accuracy of spectral codes by using 
the Immersed Boundary Method (IBM) to simulate complex geometries, while retaining high order accuracy. 
Its user-friendliness, simplicity, versatility, accuracy, scalability, portability and efficiency 
makes it an attractive tool for the Computational Fluid Dynamics community.

Xcompact3d is currently able to solve the incompressible and low-Mach number variable density 
Navier-Stokes equations up to a sixth-order accuracy using compact finite-difference schemes 
with a spectral-like accuracy on a monobloc Cartesian mesh.  
It was initially designed in France in the mid-90's for serial processors and later ported to HPC systems. 
It can now be used efficiently on hundreds of thousands CPU cores to investigate turbulence 
and heat transfer problems thanks to the open-source library 2DECOMP&FFT, 
which is a Fortran-based 2D pencil/1D slabs decomposition framework to support building 
large-scale parallel applications on distributed memory systems using MPI. 
The library has a distributed Fast Fourier Transform module as well as I/O capabilities.

Fractional time stepping is used for the time advancement, solving a Poisson equation to enforce the incompressible condition. 
The Poisson equation is fully solved in spectral space via the use of relevant 3D Fast Fourier transforms (FFTs),
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

Documentation
=============

[Readthedocs](https://xcompact3d.readthedocs.io/en/latest/)



