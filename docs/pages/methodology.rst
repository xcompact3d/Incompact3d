===========
Methodology
===========

**Xcompact3d 3.0** is based on high-order finite-difference schemes on a Cartesian mesh and conventional time advancement schemes such as Adams–Bashforth and Runge–Kutta methods. 

The main originality of **Xcompact3d 3.0** is that the Poisson equation for the incompressibility of the velocity field is fully solved in spectral space via the use of relevant 3D Fast
Fourier transforms (FFTs). With the help of the concept of modified wavenumber (see `Lele 1992 <https://www.sciencedirect.com/science/article/pii/002199919290324R>`_\), the divergence free condition is ensured up to machine accuracy.  The pressure mesh is staggered from the velocity one by half a mesh to avoid spurious pressure oscillations observed in a fully collocated approach.

The simplicity of the mesh allows an easy implementation of a 2D domain decomposition based on pencils (see the open-source library `2Decomp&FFT <http://www.2decomp.org/>`_\).  The computational domain is split into a number of sub-domains (pencils) which are each assigned to an MPI-process.  The derivatives and interpolations in the x-direction (y-direction, z-direction) are performed in X-pencils (Y-pencils, Z-pencils), respectively. The 3D FFTs required by the Poisson solver are also broken down as series of 1D FFTs computed in one direction at a time. Global transpositions to switch from one pencil to another are performed with the MPI command ``MPI_ALLTOALL(V)``.

Finite-difference schemes
-------------------------

For the first derivatives, it is recommanded to use the classic sixth-order schemes ``ifirstder=4``. This is the set-up by default.

For the interpolations (needed to compute the pressure gradients from the pressure mesh to the velocity mesh, and to compute the right hand side term of the Poisson equation from the velocity mesh to the pressure mesh), three options are available:

*``ipinter=1``: conventional sixth-order interpolation coefficients as described in `Lele 1992 <https://www.sciencedirect.com/science/article/pii/002199919290324R>`_\

*``ipinter=2``: optimal sixth-order interpolation coefficients designed to be as close as possible to spectral interpolators.

*``ipinter=3``: aggressive sixth-order interpolation coefficients designed to add some numerical dissipation at small scales but they could result in spurious oscillations close to a wall.
Ther



More details about the numerical methods can be found in:
* Laizet, S., & Lamballais, E. (2009). High-order compact schemes for incompressible flows: A simple and efficient method with quasi-spectral accuracy. Journal of Computational Physics, 228(16), 5989-6015.

* Laizet, S., & Li, N. (2011). Incompact3d: A powerful tool to tackle turbulence problems with up to O (105) computational cores. International Journal for Numerical Methods in Fluids, 67(11), 1735-1757.

* Bartholomew, P., Deskos, G., Frantz, R. A., Schuch, F. N., Lamballais, E., & Laizet, S. (2020). Xcompact3D: An open-source framework for solving turbulence problems on a Cartesian mesh. SoftwareX, 12, 100550.