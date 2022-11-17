===========
Methodology
===========

Xcompact3d is a powerful high-order flow solver for academic research. Dedicated to Direct and Large Eddy Simulations (DNS/LES), it can combine the versatility of industrial codes with the accuracy of spectral codes. It can scale with up to hundreds of thousands computational cores.


A small `video <https://www.youtube.com/watch?v=Uy-IHAjmQ3M>`_ is available with a description of the key ingredients of Xcompact3d!

The simple mesh structure allows for an easy implementation within a 2D domain decomposition strategy, based on the message passing interface (MPI). The computational domain is split into several sub-domains (pencils), which are each assigned an MPI process. The derivatives and interpolations in the *x*, *y*, and *z* directions are performed within the *x*, *y*, and *z* pencils, respectively. The 3D FFTs, required by the Poisson solver, are broken down into a series of 1D FFTs computed in one direction at a time. Global transpositions to switch from one pencil to another are performed with the MPI command MPI_ALLTOALL(V).

The incompressible Navier-Stokes equations are discretized with finite-difference sixth-order schemes on a Cartesian mesh. Explicit or semi-implicit temporal schemes can be used for the time advancement depending on the flow configuration. To treat the incompressibility condition, a fractional step method requires to solve a Poisson equation. This equation is fully solved in spectral space via the use of relevant 3D Fast Fourier transforms (FFTs), allowing any kind of boundary conditions for the velocity field in each spatial direction. Using the concept of the modified wavenumber, the divergence free condition is ensured up to machine accuracy. The pressure field is staggered from the velocity field by half a mesh to avoid spurious oscillations.

The modelling of a solid body inside the computational domain is performed with a customized Immersed Boundary Method. It is based on a direct forcing to ensure a no-slip boundary condition at the wall of the solid body while imposing non-zero velocities inside the solid body to avoid discontinuities on the velocity field. This customised IBM, fully compatible with the 2D domain decomposition and with a possible mesh refinement at the wall, is based on a 1D expansion of the velocity field from fluid regions into solid regions using Lagrange polynomials.

To reach realistic Reynolds numbers, an implicit LES strategy can be implemented to solve the Navier-Stokes equations without any extra explicit modelling. In order to mimic a subgrid-scale model, artificial dissipation can be added via the viscous term thanks to the artificial dissipative features of the high-order compact schemes.

**More information about the numerical methods can be found in:**

* Bartholomew P., Deskos G., Frantz R.A.S., Schuch F.N., Lamballais E. & Laizet S. (2020). **Xcompact3D: An open-source framework for solving turbulence problems on a Cartesian mesh**, SoftwareX, vol 12, 100550.

* Laizet, S., & Lamballais, E. (2009). **High-order compact schemes for incompressible flows: A simple and efficient method with quasi-spectral accuracy**. Journal of Computational Physics, 228(16), 5989-6015.

* Lamballais, E., Fortuné, V., & Laizet, S. (2011). **Straightforward high-order numerical dissipation via the viscous term for direct and large eddy simulation**. Journal of Computational Physics, 230(9), 3270-3275.

**More information about the parallel strategy of the code can be found in:**

* Laizet, S., & Li, N. (2011). **Incompact3d: A powerful tool to tackle turbulence problems with up to O (10\ :sup:`5`\ ) computational cores**. International Journal for Numerical Methods in Fluids, 67(11), 1735-1757.

* Li, N., & Laizet, S. (2010, May). **2DECOMP & FFT-a highly scalable 2d decomposition library and FFT interface**. In Cray User Group 2010 conference (pp. 1-13).

* Laizet, S., Lamballais, E., & Vassilicos, J. C. (2010). **A numerical strategy to combine high-order schemes, complex geometry and parallel computing for high-resolution DNS of fractal generated turbulence**. Computers & Fluids, 39(3), 471-484.

**More information about the immersed boundary methods can be found in:**

* Giannenas, A. E., & Laizet, S. (2021). **A simple and scalable immersed boundary method for high-fidelity simulations of fixed and moving objects on a Cartesian mesh**. Applied Mathematical Modelling, 99, 606-627.

* Gautier, R., Laizet, S., & Lamballais, E. (2014). **A DNS study of jet control with microjets using an immersed boundary method**. International Journal of Computational Fluid Dynamics, 28(6-10), 393-410.

**More information about the high-order numerical dissipation approach can be found in:**

* Mahfoze, O. A., & Laizet, S. (2021). **Non-explicit large eddy simulations of turbulent channel flows from Reτ= 180 up to Reτ= 5,200**. Computers & Fluids, 228, 105019.

* Dairay, T., Lamballais, E., Laizet, S., & Vassilicos, J. C. (2017). **Numerical dissipation vs. subgrid-scale modelling for large eddy simulation**. Journal of Computational Physics, 337, 252-274.

* Lamballais, E., Fortuné, V., & Laizet, S. (2011). **Straightforward high-order numerical dissipation via the viscous term for direct and large eddy simulation**. Journal of Computational Physics, 230(9), 3270-3275.


