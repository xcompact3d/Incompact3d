List of files with description
==============================

SRC directory
-------------
This directory contains all the Fortran files related to solving the Navier-Stokes equations via a fractioanl step method on a Cartesian mesh using sixth-order implicit finite-difference schemes.

``xcompact3d.f90``: This is the main file where the time advancement of the Navier-Stokes equations is performed, as well as the post-processing and checkpointing procedure. This file also contain the subroutines for the initialisation and finalisation of the code.

``time_integrators.f90`` and ``implicit.f90``: These files contains the subroutines related to the time advancement of the Navier-Stokes equations via a fractional step method. Typically, you would use an explicit time advancement (Adams Bashforth/ Runge Kutta). For wall-bounded flows, it is possible to use a semi-implicit approach for the viscous terms.

``transeq.f90``: This file contains the subroutines related to the calculcation of the righ hand side terms of the Navier-Stokes equations. Different options are available depending if you are solving the incompressible Navier-Stokes equations or if you are solving the compressible Navier-Stokes equations in the low Mach number limit. It is possible to use a DNS approach, an explicit LES approach or an implicit LES approach.

``poisson.f90``: This file contains the subroutine related to solving the Poisson equation in the spectral space. The pressure field is staggered by half a mesh with respect to the velocity field, hence the need of staggered Fast Fourier Transforms (FFTs). The subroutines ``wave`` and ``abxyz`` are defining the modified wave numbers and transfer functions. Different Poisson solver are available, depending on the type of boundary conditions: (0,0,0), (1/2,0,0), (0,1/2,0),(1/2,1/2,0) and (1/2,1/2,1/2).
If a stretching is used in the Y-direction, then the pressure field can be obtained by inverting a pentadiagonal matrix, hence the need for extra subroutines such as ``matrice_refinement``. Depending on the type of boundary conditions and the stretching, different matrices can be constructed and inverted when a stretching is used.  

``navier.f90``: This files contains a series of subroutines related to solving the Navier-Stokes equations via a fractional step method. Three notable subroutines are "divergence" to compute the divergence of the velocity field on the staggered mesh, ``gradp`` to compute the pressure gradients on the velocity mesh, and "cor_vel" to correct the velocity by the pressure gradient to obtain a divergence free field on the pressure mesh. The subroutine ``pre_correc`` is important as this is where you can impose the boundary conditions on the intermediate velocity field, before computing the Poisson equations. Imposing the boundary conditions after the correction by the pressure gradients (ie, at the end of a time step), would result in losing the divergence free condition. As a result, the boundary conditions are imposed on the intermediate velocity field. Most of the other subroutines are related to the compressible Navier-Stokes equations in the low Mach number limit.

``derive.f90``: this files contains all the subroutine for first and second derivations and interpolations. Different options are available depending on the boundary conditions and depending on the mesh (velocity mesh to pressure mesh [pv], pressure mesh to velocity mesh [vp] or velocity mesh to velocity mesh). A conventional Thomas algorithm is used to inverse the tri-diagonal matrices (as the spatial discretisation schemes are implicit).

``schemes.f90``: This files contains all the coefficients for the derivations and interpolations performed in the ``derive.f90`` file.

``acl_XXX.f90``, ``adm.f90``, ``airfoils.f90``, and ``turbine.f90``: These files only concerns the wind farm simulator. They contains subroutines for the implementation of actuator strategies to model a wind turbine.

``BC_XXX.f90``: These files contains all the subroutines required for a specific flow configuration (initialisation, post-processing, visualisations boundary conditions). Xcompact3d is designed around these files. Each BC file should have a corresponding directory in the examples directory, with a ready-to-run input file. The idea is that when a user wants to simulate a new flow configuration, they will only have to create a new BC file, without having to make changes in the subroutines used to solve the governing equations.

``case.f90``: This file contains pointers to each of the different flow configurations (the BC files) for the initialisation, imposition of the boundary conditions, post-processing and visualisation. It means that the bulk of the code does not need to be changed for each flow configuration.

``les_models.f90``, ``dynstall.f90`` and ``dynstall_legacy.f90``: These files contains the subroutines for various LES  explicit models. Our recommandation is to use an implicit approach for LES but an explicit approach such as the dynamic Smagorinsky model can be used for comparison. Note that these explicit models are not supported and might only work for a limited number of boundary conditions.

``filters.f90``: This file contains subroutines for filtering the data (to be used for an LES approach for example). Different options are possible depending on the boundary conditions.

``forces.f90``: This file contains the subroutines required to compute lift and drag when an immersed object is included in the computational domain. Only tested for the flow around a circular cylinder. Require inputs for the volume control via the input file.

``genepsi3d.f90`` and ``ibm.f90``: These file contain the subroutines related to the various immersed boundary methods implemented in the code. Some of these methods are based on 1D reconstructions of the velocity fields (one reconstruction per velocity component per spatial direction).

``module_param.f90``: This file contains all the modules with 1D and 2D key coefficients/arrays needed to run the code.

``parameters.f90``: This file contains the subroutine ``parameter`` which effectively reads the input files. It also contains the subroutine ``parameter_defaults`` which initialise all the key parameters of the code.

``statistics.f90``, ``visu.f90`` and ``probes.f90``: These files contains some of the subroutines for post-processing the data (to compute mean quantities for example, and to save 2D/3D snapshots of the computational domain).

``variables.f90``: This is where all the variables for a simulation are properly defined (for instance the pressure field, the three components of the velocity field and all the work arrays).

``constant.f90``: Small file where some important constants are defined.

``tools.f90``: This file contains a series of subroutine related to running a simulation, such as the checkpointing procedure (to restart a simulation), subroutines to invert penta, septa and nona diagonal matrices, and subroutines to compute min,max values of an array.

EXAMPLES directory
------------------
This directory contains a series of input files for various ready-to-run simulations.

``Taylor-Green-Vortex``: Folder with input files for the Taylor-Green vortex flow. It also contains some data files for comparison. 

**input_DNS_Re1600.i3d**: input file for a 3D DNS with a Reynolds number equal to 1,600. 

**input_2D.i3d**: input file for a 2D DNS with a Reynolds number equal to 1,600.

**input_ILES_Re5000.i3d**: input file for an 3D implicit LES with a Reynolds number equal to 5,000.

``Cylinder``: Folder with input files for the flow around a circular cylinder. 

**input_DNS300_LR.i3d**: input file for a 3D DNS with a Reynolds number equal to 300 and a fixed cylinder.

**input_DNS300_LR_2D.i3d**: input file for a 2D DNS with a Reynolds number equal to 300 and a fixed cylinder.

**input_DNS300_LR_MOVING.i3d**:input file for a 3D DNS with a Reynolds number equal to 300 and a moving cylinder.

``Channel-Flow``: Folder with input files for the periodic turbulent channel flow. It also contains a Fortran programme called ``stretching_parameter_channel.f90`` to determine the correct ``beta`` parameter for the stretching of the mesh in the wall-normal direction.

**input_DNS_Re180_LR_explicittime.i3d**: input file for a 3D DNS with a Reynolds number of 180.

**input_WALE_LES.i3d**: input file for a 3D explicit LES (WALE model) with a Reynolds number of 180.

``Turbulent-Boundary-Layer``: Folder with input files for a spatially evolving turbulent boundary layer (zero pressure gradient).

**input_LR_lowRe_implicit.i3d**: input fiel for a 3D DNS for low Reynolds numbers.

``2D-hill``: Folder with input files for the turbulent flow around a 2D periodic hill.

**input_LR_Re1000.i3d**: input file for a 3D DNS with a Reynolds number of 1000.


``Cavity``: Folder with input files for a cavity flow.

**input_2D_LR_Re14084.i3d**: input file for a 2D DNS with a Reynolds number of 14084.

``Lock-exchange``: Folder with input files for gravity currents in the lock exchange set-up (finite release of heavy fluid in light ambient).

**input_gravitycurrent_3DRe2236_LR.i3d**: input file for a 3D DNS with a Reynolds number of 2236.

``Mixing-layer``: Folder with input files for a mixing-layer flow.

**input_2D_LR_periodic.i3d**: input file for a 2D DNS with a Reynolds number of 400.

``Wind-Turbine``: Folder with input files for wind turbines simulations.

**two_turbines/input.i3d**: input files for a two-turbine set-up with aligned turbines.

``ABL``: Folder with input files for atmospheric boundary layer simulations.

**input_neutral.i3d**: input file for a 3D LES of a neutrally stable atmospheric boundary layer at 10 m/s.

2DECOMP directory
------------------
This directory contains all the files related to the 2DECOMP&FFT library. More information can be found `here <https://2decomp-fft.github.io/.>`_.
