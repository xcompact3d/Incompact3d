===========
Methodology
===========

**Xcompact3d 3.0** is based on high-order finite-difference schemes on a Cartesian mesh and conventional time advancement schemes such as Adams–Bashforth and Runge–Kutta methods. 

The main originality of **Xcompact3d 3.0** is that the Poisson equation for the incompressibility of the velocity field is fully solved in spectral space via the use of relevant 3D Fast
Fourier transforms (FFTs). With the help of the concept of modified wavenumber (see `Lele 1992 <https://www.sciencedirect.com/science/article/pii/002199919290324R>`_\), the divergence free condition is ensured up to machine accuracy.  The pressure mesh is staggered from the velocity one by half a mesh to avoid spurious pressure oscillations observed in a fully collocated approach.

