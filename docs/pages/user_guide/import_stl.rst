How to import an STL file and run a simulation with an immersed body
====================================================================

**This tutorial is designed for an Ubuntu workstation for the simulation around a sphere**.

Install the xcompact3d toolbox
------------------------------

First you need to make sure to have ``conda`` installed on your machine. See for instance `here <https://linuxconfig.org/installing-anaconda-on-ubuntu-24-04>`_ for an installation on Ubuntu. We only tried to install the xcompact3d toolbox with Python 3.10 and Python 3.12.

.. code-block::

	conda create -n x3d python=3.XX
	conda activate x3d
	pip install xcompact3d-toolbox

Install and compile Xcompact3d
------------------------------

.. code-block::

	git clone https://github.com/xcompact3d/Incompact3d.git
	cd Incompact3d/
	export FC=mpif90
	cmake -S . -B build
	cd build/
	cmake --build . -j 8

The executable file is in the ``build/bin`` directory

Bug fix in the xcompact3d_toolbox package
-----------------------------------------

open ``/home/username/anaconda3/envs/x3d/lib/python3.12/site-packages/xcompact3d_toolbox/sandbox.py`` and change all instances of longdouble to double (lines 438-441)

Generate the initial conditions including the epsilon function
--------------------------------------------------------------
.. code-block::

	cd examples/Sphere
	python3 generate_initial_conditions.py
	
Run the simulation on 8 cores (~7minutes)
-----------------------------------------
.. code-block::
	
	mpirun -np 8 ../../build/bin/xcompact3d


