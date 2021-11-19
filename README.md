## The Xcompact3d code

Xcompact3d is a Fortran-based framework of high-order finite-difference flow solvers dedicated to the study of turbulent flows. Dedicated to Direct and Large Eddy Simulations (DNS/LES) for which the largest turbulent scales are simulated, it can combine the versatility of industrial codes with the accuracy of spectral codes. Its user-friendliness, simplicity, versatility, accuracy, scalability, portability and efficiency makes it an attractive tool for the Computational Fluid Dynamics community.

XCompact3d is currently able to solve the incompressible and low-Mach number variable density Navier-Stokes equations using sixth-order compact finite-difference schemes with a spectral-like accuracy on a monobloc Cartesian mesh.  It was initially designed in France in the mid-90's for serial processors and later converted to HPC systems. It can now be used efficiently on hundreds of thousands CPU cores to investigate turbulence and heat transfer problems thanks to the open-source library 2DECOMP&FFT (a Fortran-based 2D pencil decomposition framework to support building large-scale parallel applications on distributed memory systems using MPI; the library has a Fast Fourier Transform module).
When dealing with incompressible flows, the fractional step method used to advance the simulation in time requires to solve a Poisson equation. This equation is fully solved in spectral space via the use of relevant 3D Fast Fourier transforms (FFTs), allowing the use of any kind of boundary conditions for the velocity field. Using the concept of the modified wavenumber (to allow for operations in the spectral space to have the same accuracy as if they were performed in the physical space), the divergence free condition is ensured up to machine accuracy. The pressure field is staggered from the velocity field by half a mesh to avoid spurious oscillations created by the implicit finite-difference schemes. The modelling of a fixed or moving solid body inside the computational domain is performed with a customised Immersed Boundary Method. It is based on a direct forcing term in the Navier-Stokes equations to ensure a no-slip boundary condition at the wall of the solid body while imposing non-zero velocities inside the solid body to avoid discontinuities on the velocity field. This customised IBM, fully compatible with the 2D domain decomposition and with a possible mesh refinement at the wall, is based on a 1D expansion of the velocity field from fluid regions into solid regions using Lagrange polynomials or spline reconstructions. In order to reach high velocities in a context of LES, it is possible to customise the coefficients of the second derivative schemes (used for the viscous term) to add extra numerical dissipation in the simulation as a substitute of the missing dissipation from the small turbulent scales that are not resolved. 


### External Resources

- [**Twitter**](https://twitter.com/incompact3d)

## Source Download and Compilation

First, make sure you have all the [required dependencies](#required-build-tools-and-external-libraries) installed.
Then, acquire the source code by cloning the git repository:

    git clone git@github.com:xcompact3d/Incompact3d.git

(If you are behind a firewall, you may need to use the `https` protocol instead of the `git` protocol:

    git config --global url."https://".insteadOf git@

Be sure to also configure your system to use the appropriate proxy settings, e.g. by setting the `https_proxy` and `http_proxy` variables.)

By default you will be building the latest unstable version of Incompact3d. However, most users should use the most recent stable version of Incompact3d, which is currently the `3.0` series of releases. You can get this version by changing to the Incompact3d directory and running

    git checkout v3.0

Now run `make` to build the `Incompact3d` executable. To perform a parallel build, use `make -j N` and supply the maximum number of concurrent processes. (See [Platform Specific Build Notes] for details.)
This takes a while, but only has to be done once. If the defaults in the build do not work for you, and you need to set specific make parameters, you can save them in `Make.user`. The build will automatically check for the existence of `Makefile` and use it if it exists.
Building Incompact3d requires very little of disk space and virtual memory.

**Note:** The compiling process

In the Incompact3d, once you have selected the correct options for your Fortran compiler, you just need to do

    make clean
    
to make sure that you will be compiling all the files, and then

    make 
   
Once it is built, you just need to go in one of the examples directories, for instance https://github.com/xcompact3d/Incompact3d/tree/master/examples/Taylor-Green-Vortex and from there use the input.i3d file to configure your simulation. To get to know the code, you can start with a ready-to-run input file, see as an example https://github.com/xcompact3d/Incompact3d/blob/master/examples/Taylor-Green-Vortex/input_DNS_Re1600.i3d which can be use to run the Taylor-Green case in a DNS set-up at Re=1600. Using 16 CPU cores, this simulation should last less than 5 minutes. The command to launch the simulation is
    
    mpirun -np 16 ../../xcompact3d 

or

    nohup mpirun -np 16 ../../xcompact3d > output.out &

### Optional ADIOS2 I/O backend

As part of the ARCHER2 eCSE0302 project, an optional I/O backend using ADIOS2 has been added to the 2DECOMP&FFT library distributed with Xcompact3d.
This is enabled at compile time, by default the original MPIIO backend will be used, to enable the ADIOS2 backend build as

    make clean
    make IO=adios2 ADIOS2DIR=${ADIOS2_DIR}
    
where `${ADIOS2_DIR}` points to the install location of ADIOS2.
ADIOS2 enables configuring the I/O behaviour at runtime using an `xml` configuration file - see the example at `examples/Taylor-Green-Vortex/adios2_config.xml`.
