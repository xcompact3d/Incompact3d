Taylor-Green Vortices
=====================
The Taylor–Green Vortex (TGV) is a well-known benchmark in CFD, modelling the transition from an initially laminar state to a fully turbulent one. It is attractive as a test case due to its simple setup, with the possibility to use a combination of
periodic and free-slip boundary conditions

Have a look at the following `video <https://www.youtube.com/watch?v=yj0njXod7iU>`_ to see how easy it is to run a Taylor Green Vortex simulation with Xcompact3d. 

In the video, it is shown how to perform a Direct Numerical Simulation of the Taylor Green Vortex case with a Reynolds number equal to 1,600.

You will see how to use **Paraview** to visualise the snapshots and **Gnuplot** to visualise the statistics. Examples are provided with free slip and periodic boundary conditions. The key parameters in the input file are also discussed.

For the free slip case, ``input_DNS_Re1600.i3d`` from the ``examples/Taylor-Green-Vortex`` directory is used, and in the video you can see how to modify this input file for the periodic case.

Reference data for Reynolds numbers ranging from 1,250 to 20,000 can be found `here <https://zenodo.org/record/2577239#.YsV6GozMI5k>`_.

.. image:: snap4.png
  :width: 1200
