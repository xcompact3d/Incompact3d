# Testing 

This folder contains all the necessary files to test Xcompact3d on most of the canonical test cases included into the code. 
Generally the tests run for up to 500 time steps and they can be broadly divided into: 

1. Functional tests to verify that a specific test case reach conclusion without diverging 
1. Verification tests where the solution is compared to a reference dataset

At current stage, all tests belong to category (1) with the exception of the Taylor-Green-Vortex (TGV), where the time evolution of 
some averaged quatities is compared with reference values.
More detials about the TGV comparison is given [here](TGV-Taylor-Green-vortex/README.md)
In all cases the run of the tests, done via ``mpirun``, and the eventual comparison, performed via ``python3``, is instantiated within CTest.
The full lists of the test performed is: 


1. Atmospheric Boundary layer (ABL) in neutral conditions (new set-up)
1. Differentially heated cavity
1. Turbulent Channel Flow with X as streamwise direction
1. Turbulent Channel Flow with Z as streamwise direction
1. Flow around a circular cylinder
1. Flow around a moving circular cylinder
1. Lock exchange
1. Mixing Layer
1. Periodic hill
1. Turbulent Boundary Layer (TBL)
1. Wind Turbine
1. Taylor Green Vortex (TGV)

By default only the  TGV case is activated, while the full 
testing suite needs to be enable by using the `BUILD_TESTING_FULL` flag as 
```
$ cmake --build $path_to_build_directory -DBUILD_TESTING_FULL=ON 
```
or by using `ccmake`.

The tests are performed using `CTest` as  
```
$ ctest --test-dir $path_to_build_directory
```

Every test is performed in a dedicated working directory that is located under the following path 
```
$ /path/to/build/RunTests
```
All standard outputs from all test runs are collated under the file
```
$ /path/to/build/Testing/Temporary/LastTest.log
```
together with additional files detailing additional informations such as 
the elapse time for the different tests and the eventual failed cases. 
