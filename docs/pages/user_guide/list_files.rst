List of files with description
==============================

SRC directory
-------------
``xcompact3d.f90``: This is the main file where the time advancement of the Navier-Stokes equations is performed, as well as the post-processing and checkpointing procedure. This file also contain the subroutines for the initialisation and finalisation of the code.

``time_integrators.f90``: This file contains the subroutine related to the time advancement of the Navier-Stokes equations via a fractional step method. Typically, you would use an explicit time advancement (Adams Bashforth/ Runge Kutta). For wall-bounded flows, it is possible to use a semi-implicit approach for the viscous terms.

``transeq.f90``: This file contains the subroutines related to the calculcation of the righ hand side terms of the Navier-Stokes equations. Different options are available depending if you are solving the incompressible Navier-Stokes equations or if you are solving the compressible Navier-Stokes equations in the low Mach number limit. It is possible to use a DNS approach, an explicit LES approach or an implicit LES approach.

``poisson.f90``: This file contains the subroutine related to solving the Poisson equation in the spectral space. The pressure field is staggered by half a mesh with respect to the velocity field, hence the need of staggered Fast Fourier Transforms (FFTs). The subroutines ``wave`` and ``abxyz`` are defining the modified wave numbers and transfer functions. Different Poisson solver are available, depending on the type of boundary conditions: (0,0,0), (1/2,0,0), (0,1/2,0),(1/2,1/2,0) and (1/2,1/2,1/2).
If a stretching is used in the Y-direction, then the pressure field can be obtained by inverting a pentadiagonal matrix, hence the need for extra subroutines such as ``matrice_refinement``. Depending on the type of boundary conditions and the stretching, different matrices can be constructed and inverted when a stretching is used.  

``navier.f90``: 

EXAMPLES directory
------------------

2DECOMP directory
------------------
