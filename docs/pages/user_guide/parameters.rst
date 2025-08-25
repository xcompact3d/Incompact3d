List of Parameters
==================

As part of the upgrade from **Incompact3d 2.0** to the now called **Xcompact3d 4.0**, the parameters file was reconfigured to support *NAMELIST I/O*, which produces format-free input for groups of variables. 

For a complete view about how the `.i3d` files are handled, besides to the default value applied to each parameter when not defined by the user, take a look at `parameters.f90 <https://github.com/xcompact3d/Incompact3d/blob/master/src/parameters.f90>`_\ .

BasicParam
----------

* ``p_row`` & ``p_col`` define the domain decomposition for (large-scale) parallel computation. Notice that the product of both should be equal to the number of computational cores where Xcompact3d will run, and ``p_row = p_col = 0`` executes the code in auto-tuning mode. More information can be found at `2DECOMP&FFT <http://www.2decomp.org>`_\ ).

* ``nx``, ``ny`` & ``nz`` are the number of mesh points in each direction; Because of the Fast Fourier Transforms, limitations are in place. Basically you need to pick a combination of power of prime numbers. If the boundary conditions are not periodic (``nclXX`` not equal to zero), then you need to add an extra mesh nodes. For instance 1025 (non periodic boundary conditions) instead of 1024 (periodic boundary conditions).

* ``xlx``, ``yly`` & ``zlz`` are the domain size, normalised with the reference length for the simulation (as an example, for a cylinder case, ``xlx=20D`` where ``D`` is the diameter of the cylinder). You can also decide to run simulations with non-normalised quantities, you just need to be consistent! 

* ``itype`` sets the flow configuration, each one is specified in a different ``Case-<flow-configuration>.f90`` file. They are:

    - 0 - User Custom Configuration;
    - 1 - Turbidity Current in Lock-Release;
    - 2 - Taylor-Green Vortices;
    - 3 - Periodic Turbulent Channel;
    - 4 - Periodic Hill;
    - 5 - Flow over a Cylinder;
    - 7 - Mixing Layer;
    - 9 - Turbulent Boundary Layer;
    - 10 - Atmospheric Boundary Layer;
    - 11 - Uniform flow;
    - 12 - Sandbox configuration;
    - 13 - Differentilly heated cavity;
    - 14 - Pipe flow;
    - 15 - Periodic turbulent boundary layer.
You can modify the ``Case-<flow-configuration>.f90`` file to change the inlet and/or initial conditions. Except for the Atmospheric Boundary Layer, quantities are normalised with a reference velocity, a reference length and a constant density, all equal to 1 (as an example, for the cylinder, the reference velocity is the freestream velocity equal to 1,  the diameter equal to 1 and a constant density equal to 1; as a results the Reynolds number is equal to 1/nu).


* ``istret`` controls mesh refinement in y direction only (it is not possible to refine the mesh in more than one direction):

    - 0 - No refinement (default);
    - 1 - Refinement at the center;
    - 2 - Both sides;
    - 3 - Just near the bottom.
    
    More details about the refinement function (which could be changed, but it would have to be represented in the spectral space by few modes only) can be found in *Laizet, S., & Lamballais, E. (2009)*, **High-order compact schemes for incompressible flows: A simple and efficient method with quasi-spectral accuracy**, *Journal of Computational Physics, 228(16), 5989-6015*. The refinement can be control with the parameter ``beta``.

* ``beta`` is the refinement parameter. Large positive ``beta`` will lead to an almost uniform mesh, small positive beta will lead to very stretched mesh. Best option to find a suitable ``beta`` for your simulation is test and trial errrors!

* ``iin`` defines perturbation at the initial condition:

    - 0 - No random noise (default);
    - 1 - Random noise with amplitude of ``init_noise``;
    - 2 - Random noise with fixed seed (important for development and debugging) and amplitude of ``init_noise``.

    .. note::
      The exactly behavior may be different according to each flow configuration. 

* ``inflow_noise`` Random number amplitude at inflow boundary, expressed as a % of the reference velocity (as an example 0.125 correspond to 12.5% of the reference velocity. It is advise to use small values between 0. and 0.1.

* ``re`` Reynolds number, defined using the reference length scale and reference velocity, hence Re=1/nu.

* ``dt`` Time step, to be inputed manually, depending on the mesh resolution, accuracy of the spatial finite-difference schemes and temporal scheme. By experience, the optimal time step can be found by test and trial errors, as opposed to try to compute the CFL number.

* ``ifirst`` First iteration of the simulation. Do not forget to update if you are using a restart/checkpointing file.

* ``ilast`` Last iteration of the simulation.

* ``numscalar`` Number of scalar in the simulation. When using passive scalars, you can have more than one scalar field at the same time (basically you can solver the same transport equations with different initial conditions.

* ``iscalar`` Enables scalar field(s). It is defined to 1 automatically when ``numscalar > 0``.

* ``iibm`` Flag for Immersed Boundary Method:

    - 0 - Off (default);
    - 1 - On with direct forcing method, i.e., it sets velocity to zero inside the solid body.
    - 2 - On with alternating forcing method, i.e, it uses Lagrangian Interpolators to define the velocity inside the body and imposes no-slip condition at the solid/fluid interface, see *Gautier, R., Laizet, S., & Lamballais, E. (2014)*, **A DNS study of jet control with microjets using an immersed boundary method**, *International Journal of Computational Fluid Dynamics, 28(6-10), 393-410* [no longer supported].
    - 3 - On with alternating forcing method, i.e, it uses Cubic Spline Interpolators to define the velocity inside the body and imposes no-slip condition at the solid/fluid interface. Allows for moving objects, see *Giannenas, A. E., & Laizet, S. (2021)*, **A simple and scalable immersed boundary method for high-fidelity simulations of fixed and moving objects on a Cartesian mesh**, *Applied Mathematical Modelling, 99, 606-627*.

* ``ilmn`` Enables Low Mach Number methodology when set to 1, basically solving the compressible Navier-Stokes equations in the low Mach number limit. If you want to solve the incompressible Navier-Stokes equations, ``ilmn`` should be equal to 0. See *Bartholomew, P., & Laizet, S. (2019)*, **A new highly scalable, high-order accurate framework for variable-density flows: Application to non-Boussinesq gravity currents**, *Computer Physics Communications, 242, 83-94*.

* ``ilesmod`` Enables Large-Eddy simulations (LES) methodologies:

    - 0 - Off.
    - 1 - Smagorinsky 
    - 2 - WALE **- suitable for wall-bounded flows**
    - 3 - dynamic Smagorinsky **- avoid as the filtering procedure is expensive**
    - 4 - ILES **- prefered options for LES**

Please note that we will eventually remove all explicit LES models from the code as our ILES approach is cheaper with the same quality of results, if not better. 
See *Dairay, T., Lamballais, E., Laizet, S., & Vassilicos, J. C. (2017)*, **Numerical dissipation vs. subgrid-scale modelling for large eddy simulation**, *Journal of Computational Physics, 337, 252-274* and *Mahfoze, O. A., & Laizet, S. (2021)*, **Non-explicit large eddy simulations of turbulent channel flows from Reτ= 180 up to Reτ= 5,200**, *Computers & Fluids, 228, 105019*.

* ``nclx1``, ``nclxn``, ``ncly1``, ``nclyn``, ``nclz1`` & ``nclzn`` define the velocity's boundary condition:

    - 0 - Periodic boundary conditions;
    - 1 - Free-slip boundary conditions, with two options depending on the velocity components: symmetry or anti-symmetry via the parameter ``npaire``. Basically, a free-slip periodic boundary condition in the vertical direction corresponds to u=w=constant (similar to a symmetry boundary condition for which the gradient of u and w in the vertical direction is equal to zero), and to v=0 (similar to an anti-symmetry boundary condition);
    - 2 - Dirichlet boundary conditions.
    
    Note that the fractional step method in the code does not need explicit boundary conditions for the pressure field. 

* ``ivisu`` enables I/O for 3D snapshots if equal to 1 (every ``ioutput`` time step);

* ``ipost`` enables online postprocessing if equal to 1.

* ``gravx``, ``gravy`` & ``gravz`` are the three components of the unitary vector pointing in the gravity's direction, only experimented for gravity currents.

* ``cpg`` is a logical parameter for the momentum source term (Turbulent Channel case only). True / False for an imposed pressure gradient / flow rate.

* ``ifilter`` & ``C_filter`` to filter the solution [set to 0 & 0 except for ABL and wind turbines cases. ``ifilter`` activates the filtering (with different direction combinations given by different values, e.g., 1 is for all directions, 2 for x-z filtering etc) and ``C_filter`` is the filter constant (e.g., C_filter=0 filters at 2Delta, C_filter=0.5 does nothing)

* ``iturbine`` case-specific parameter for the wind farm simulator with 1: Actuator line, 2: actuator disk

NumOptions
----------

* ``ifirstder`` Scheme for first order derivative:

    - 1 - 2nd central explicite finite-difference schemes;
    - 2 - 4th central explicite finite-difference schemes **- under development**;
    - 3 - 4th compact finite-difference schemes **- under development**;
    - 4 - 6th compact finite-difference schemes **- prefered option**

* ``isecondder`` Scheme for second derivative:

    - 1 - 2nd central explicite finite-difference schemes;
    - 2 - 4th central explicite finite-difference schemes **- under development**;
    - 3 - 4th compact finite-difference schemes **- under development**;
    - 4 - 6th compact finite-difference schemes **- prefered option for DNS & Explicit LES**
    - 5 - hyperviscous 6th compact finite-difference schemes **- for Implicit LES**
   
   For more details about these different options, and in particular for the customized hyperviscous schemes, please have a look at:
   
   - *Lamballais, E., Fortuné, V., & Laizet, S. (2011).* **Straightforward high-order numerical dissipation via the viscous term for direct and large eddy simulation.** *Journal of Computational Physics, 230(9), 3270-3275.*
   - *Dairay, T., Lamballais, E., Laizet, S., & Vassilicos, J. C. (2017)*, **Numerical dissipation vs. subgrid-scale modelling for large eddy simulation**, *Journal of Computational Physics, 337, 252-274.* 

* ``itimescheme`` Time integration scheme:

    - 1 - Forwards Euler;
    - 2 - Adams-Bashforth 2;
    - 3 - Adams-Bashforth 3 **- recommended option when using immersed boundary methods**;
    - 4 - Adams-Bashforth 4 **- not implemented yet**;
    - 5 - Runge-kutta 3 **- recommended option when not using immersed boundary methods**;
    - 6 - Runge-kutta 4 **- not implemented yet**;
    - 7 - Semi-implict Crank-Nicolson + Adams-Bashforth 3 **- does not work yet for all cases**;
    - 8 - Semi-implict Crank-Nicolson + Runge-kutta 3 **- does not work yet for all cases**.

* ``iimplicit`` Time integration scheme for the Y-diffusive term:

    - 0 - Explicit, default;
    - 1 - Euler implicit **- does not work yet for all cases**;
    - 2 - Crank-Nicolson **- does not work yet for all cases**;

* ``nu0nu`` Ratio between hyperviscosity/viscosity at nu;

* ``cnu`` Ratio between hypervisvosity at :math:`k_m=2/3\pi` and :math:`k_c= \pi`.

For more details about these two parameters, please have a look at:
   
   - *Lamballais, E., Fortuné, V., & Laizet, S. (2011).* **Straightforward high-order numerical dissipation via the viscous term for direct and large eddy simulation.** *Journal of Computational Physics, 230(9), 3270-3275.*
   - *Dairay, T., Lamballais, E., Laizet, S., & Vassilicos, J. C. (2017)*, **Numerical dissipation vs. subgrid-scale modelling for large eddy simulation**, *Journal of Computational Physics, 337, 252-274.* 

Default values are ``cnu=0.44`` and ``nu0nu=4`` which are ideal in a DNS context. 

* ``ipinter`` 
    - 1 - conventional sixth-order interpolation coefficients as described in `Lele 1992 <https://www.sciencedirect.com/science/article/pii/002199919290324R>`_\;
    - 2 - optimal sixth-order interpolation coefficients designed to be as close as possible to spectral interpolators **- recommended option**;
    - 3 - aggressive sixth-order interpolation coefficients designed to add some numerical dissipation at small scales but they could result in spurious oscillations close to a wall.

InOutParam
----------

* ``irestart`` Reads initial flow field if equals to 1. This is when you are using the restarting/checkpointing procedure (have a look at the subroutine ``restart``;
* ``icheckpoint`` Frequency for writing backup file (every ``icheckpoint`` time steps). Note that Xcompact3d will only keep on the disc the latest restart file;
* ``ioutput`` Frequency to generate the 3D snapshots (every ``ioutput`` time steps);
* ``nvisu`` Size for the 3D snapshots to be written on the disc (every ``nvisu`` mesh nodes). By default, use ``nvisu=1`` which corresponds to 3D snapshots of size ``nx x ny x nz``;
* ``iprocessing`` Frequency for online postprocessing **- not supported anymore as the guideline is to use our Python post-processing framework**;
* ``ninflows`` For precursor simulations for atmospheric boundary layers **- will evolve soon!**
* ``ntimesteps`` For precursor simulations for atmospheric boundary layers **- will evolve soon!**
* ``inflowpath`` For precursor simulations for atmospheric boundary layers **- will evolve soon!**
* ``ioutflow`` For precursor simulations for atmospheric boundary layers **- will evolve soon!**
* ``output2D`` **- not supported anymore (will be removed eventually). Keep default value to zero**

    - 0 - 3D, default
    - 1 - 2D, averaged over X
    - 2 - 2D, averaged over Y
    - 3 - 2D, averaged over Z

* ``nprobes`` **- not supported anymore (will be removed eventually). Keep default value to zero**

Statistics
----------

* ``wrotation`` Rotation speed (Rosby number) of a source term in the equations to trigger turbulence (Channel Flow only);
* ``spinup_time`` number of time steps after which the rotation source term is removed (Channel Flow only);
* ``nstat`` Size arrays for statistic collection (every ``nstat`` mesh nodes). By default, use ``nstat=1`` which corresponds to 3D statistic arrays of size ``nx x ny x nz``;
* ``initstat`` Time step when collection of statistics starts.

ProbesParam
-----------

* ``flag_all_digits`` When False (default), 6 digits are recorded. Set to True to record 16 digits **- not supported anymore (will be removed eventually)**;
* ``flag_extra_probes`` Default is False. Set to True to monitor the velocity / pressure / scalars gradients **- not supported anymore (will be removed eventually)**;
* ``xyzprobes`` Array of size (3,nprobes) containing the location of the probes **- not supported anymore (will be removed eventually)**.

ScalarParam
-----------

* ``sc`` Schmidt numbers;
* ``ri`` Richardson numbers;
* ``uset`` Settling velocities;
* ``cp`` Initial concentrations;
* ``nclxS1``, ``nclxSn``, ``nclyS1``, ``nclySn``, ``nclzS1`` & ``nclzSn`` define the scalar's boundary condition:

    - 0 - Periodic;
    - 1 - Odd or Even (default, no-flux);
    - 2 - Dirichlet.

* ``scalar_lbound`` & ``scalar_ubound`` are the Scalar bounds;
* ``sc_even`` True (default) if the scalar is even. False if it is odd;
* ``sc_skew`` Default is False. True to activate the skew-symmetric convection for a scalar;
* ``alpha_sc``, ``beta_sc``, ``g_sc`` are used only when ``iimplicit > 0``. They define the boundary condition for the scalar at the top and bottom walls;
* ``Tref``
* ``iibmS`` Flag for the scalar treatment at the Immersed Boundary Method (pre-release):

    - 0 - Off (default);
    - 1 - On with direct forcing method, i.e., it sets scalar to zero inside the solid body;
    - 2 - On with alternating forcing method, i.e., it uses Lagrangian Interpolators to define the scalar inside the body and imposes zero value at the solid/fluid interface;
    - 3 - On with alternating forcing method, but now the Lagrangian Interpolators are set to impose no-flux for the scalar field at the solid/fluid interface (only recommended if the normal vectors to the object's faces are aligned with one of the coordinate axes).

LESModel
--------

jles, smagcst, walecst, maxdsmagcst, iwall

WallModel
---------

smagwalldamp

Tripping
--------

itrip,A_tr,xs_tr_tbl,ys_tr_tbl,ts_tr_tbl,x0_tr_tbl

ibmstuff
--------

cex,cey,ra,nobjmax,nraf,nvol,iforces

ForceCVs
--------

xld, xrd, yld, yud

LMN
---

dens1, dens2, prandtl, ilmn_bound, ivarcoeff, ilmn_solve_temp, massfrac, mol_weight, imultispecies, primary_species, Fr, ibirman_eos

ABL
---
z_zero, iwallmodel, k_roughness, ustar, dBL, imassconserve, ibuoyancy, iPressureGradient, iCoriolis, CoriolisFreq, istrat, idamping, iheight, TempRate, TempFlux, itherm, gravv, UG, T_wall, T_top

CASE
----

* ``tgv_twod`` Flag used to initialize the Taylor-Green Vortices case with a 2D / 3D initial condition
* ``pfront``

ALMParam
--------
ialmrestart,filealmrestart,iturboutput,NTurbines,TurbinesPath,NActuatorlines,ActuatorlinesPath,eps_factor,rho_air

ADMParam
--------
Ndiscs,ADMcoords,C_T,aind,iturboutput,rho_air
