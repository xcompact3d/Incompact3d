XCompact3d installation documentation
=====================================

## XCompact3d installation documentation


XCompact3d can be installed using the Makefile provided at the top of this repository or alternatively it can 
be build using CMake. 
Instructions on installations using the Makefile provided are given in the README file.
This document provides more detailed instructions on the CMake building and installation together with 
some of the knows issues. 

```
  $ export FC=mpif90
  $ mkdir=build
  $ cd build
  $ cmake ~/path_to_sources
  $ make -j n
  $ make install 
```
where n is the number of tasks that you would like to use and are available on your installation system. 
The default installation will be located under 

    $ ~/path_to_build/opt

The installation directory will cointain:
* The *bin* directory with two execulables: **xcompaxt3d** for the main execution of the code and **xcompact3d_paraview_vtk** to convert the *.bin* files into vtk text format. The converter is useful when the default *xdmf* format is not working with Paraview
* The *example* directory with few example of input *.i3d* files for **Xcompact3d**
* The *lib* directory with the archive for the **decomp2d** library

## CTest
To test your installation you can also type in the terminal from your *build* directory

    $ make test

Four tests are performed:
* Taylor Green Vortex (TGV)
* Turbulent Channel Flow with x as streamwise direction
* Turbulent Channel Flow with z as streamwise direction
* Flow around a circular cylinder

The simulations results are located under 

    $ ~/path_to_build/Test

and the standard output from the simulations is in 

    $ ~/path_to_build/Testing/Temporary/

## Known issues
* Sometimes the *CMake* find MPI function does not properly locate the *mpiexec* for the given compiler. Please look at the output of 

     $ cmake ~/path_to_sources

and make sure that the path to *mpiexec* is the correct one

     $ -- MPI EXEC: ~/my_correct_path_to_mpiexec

If the path is not correct you might have problems in running *CTest*.
To solve the issue do the following 
  * run *ccmake* at the root of the build directory

     $ ccmake . 

  * Toggle the advance mode by hitting *t*
  * Look for the *MPIEXEC_EXECUTABLE* and set it up pointing to the correct *mpiexec*




