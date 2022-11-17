How to use the restarting procedure
===================================

Have a look at the following `video <https://youtu.be/O2C9PlfvWko>`_ to see how easy it is to use the restarting procedure in Xcompact3d. 

There are only 3 parameters to deal with:

* ``istart``: first iteration of the simulation

* ``ilast``: last iteration of the simulation

* ``irestart``: switch from starting the simulation from the initial condition defined in the BC files (0) or from starting the simulation from a restart/checkpointing file (1).

Another parameter of interest is ``icheckpoint`` which is the number of time step between writing a restart/checkpointing file.

Please note that Xcompact3d will only keep on the disk the latest restart/checkpointing file (called ``checkpoint``). Details about the restart/checkpointing files can be found in the ``*.info`` files. If you want to save different  restart/checkpointing files (intermediate ones for instance), you will have to do it manually.

The code for the restart/checkpointing files can be found in ``src/tools.f90``.
