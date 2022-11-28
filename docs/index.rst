======================================
Welcome to Xcompact3d's documentation!
======================================

Xcompact3d is a Fortran-based framework of high-order finite-difference flow solvers dedicated to the study of turbulent flows. Dedicated to Direct and Large Eddy Simulations (DNS/LES) for which the largest turbulent scales are simulated, it can combine the versatility of industrial codes with the accuracy of spectral codes. Its user-friendliness, simplicity, versatility, accuracy, scalability, portability and efficiency makes it an attractive tool for the Computational Fluid Dynamics community.

XCompact3d is currently able to solve the incompressible and low-Mach number variable density Navier-Stokes equations using sixth-order compact finite-difference schemes with a spectral-like accuracy on a monobloc Cartesian mesh.  It was initially designed in France in the mid-90's for serial processors and later converted to HPC systems based on CPU hardware. **A GPU friendly version of Xcompact3d is currently being developed**.


**IMPORTANT:**

* **It is strongly recommended to read the relevant references before starting using Xcompact3d.**

* **We kindly ask you to cite the relevant references (when suitable) in your work based on Xcompact3d.**

.. note::
    This website is still **under development**!
    
.. image:: 3D_VIEW.png
  :width: 1200

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   pages/installation.rst
   pages/benchmark_cases.rst
   pages/user_guide.rst
   pages/methodology.rst
