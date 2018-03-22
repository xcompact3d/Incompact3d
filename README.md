## The Incompact3d code

Incompact3d is a Fortran-MPI based, finite difference high-performance code for CFD.

This repository contains a major upgrade in the Incompact3d code.

The new version is faster, more flexible and user friendly.

The main homepage for Incompact3d can be found at [incompact3d.com](https://www.incompact3d.com/).

This is the BitBucket repository of Incompact3d source code, including instructions for running and compiling Incompact3d, below.

## Resources

- **Homepage:** <https://www.incompact3d.com/>
- **Binaries:** <https://www.incompact3d.com/download.html>
- **Documentation:** <https://www.incompact3d.com/docs>
- **Source code:** <https://bitbucket.org/incompact3d/incompact3d/src/v2.0/>
- **Git clone URL:** <git@bitbucket.org:incompact3d/incompact3d.git>
- **Discussion forum:** <https://www.discourse.incompact3d.com/>
- **Mailing lists:** <https://www.incompact3d.com/community.html>

New developers may find the notes in [CONTRIBUTING](https://bitbucket.org/incompact3d/incompact3d/src/v2.0/CONTRIBUTING.md) helpful to start contributing to the Incompact3d codebase.

### External Resources

- [**StackOverflow**](https://stackoverflow.com/questions/tagged/incompact3d)
- [**Youtube**](https://www.youtube.com/channel/incompact3d)
- [**Twitter**](https://twitter.com/incompact3d)


### Main improvments/changes ### 

* There are now more flexibility in terms of boundary conditions. You can choose the combinations 1-2 or 2-1 along with the classic 0-0, 1-1 and 2-2 for any direction for the velocity and the scalar.

* There is now a customized file for each case, i.e. bc_tgv.f90, bc_lock.f90, bc_channel.f90, bc_cylinder.f90, bc_jet.f90... this allows concentrating the case specific in a single file. 
The subroutines init, inflow, outflow, in-situ processing which are case specific are in this file.

* There is no need for recompiling the code if changing the mesh size or the p_row/p_col values.

* The pre_correc was redesigned to work with all boundary conditions. We are still reinfocing the normal velocity equal to zero when using free-slip (this is important when using scalar) 

* All numerical digits were changed to text variables in order to reduce occurances of _mytype and improve the readabilty of the code. Real 0 is now zero, 1 is one... this is not valid for integers! 

* The terms in the scalar subroutine were rewriten to 

* The paraview.f90 is a subroutine that generates the XDMF meta file to open the raw output in Paraview. The .xdmf�s are generated when initiallizing the code and are separated in 3D and 2D fields (i.e. friction velocity map, deposit map, mean planes).

* The restart 

* The visu.prm

* The stretching is defined with a -DSTRETCHING flag on the Makefile. If you do not add this flag and recompile the code, the stretching wont work. The stretching in x, y and z controlled by istretx,istrety and istretz is under development. 

* We removed the deprectaed TWOD flag and the content. If you wish to run a 2D case, you can do a quasi-2D case with nclz1=nclzn=0 and use nz=6 (if DNS) or nz=8 (if ILES). This is the minimum number of points required by the stencil. You can use zlz=0.05, for example. 




### Cannonical cases for code validation ###

1) BC-TGV.f90 ** 3D Taylor-Green problem (TGV)
Reference paper/data: Dairay et al. (2017)

2) BC-Channel-flow.f90 ** Turbulent channel flow
Reference paper/data: Moser & Kim (1999)

3) BC-Cylinder.f90 ** Flow over a cylinder
Reference paper/data: 

4) BC-Lock-exchange.f90 ** Lock-exchange problem
Reference paper/data: Espath et al. (2014)

5) BC-Periodic-hill.f90.f90 ** Periodic Hill
Reference paper/data: 

6) BC-TBL.f90 ** Turbulent Boundary Layer
add the implicit, stretching and tripping...


## Base flow configurations

Incompact3d is built in a flow configuration specific files. You must choose a case and set it on the 'Makefile' and recompile.

| Case             | Filename         |
|:----------------:|:----------------:|
| Turbulent Channel            | Channel-flow     |
| Flow over a Cylinder         | Cylinder         |
| Taylor-Green vortices        | TGV              |
| Periodic Hill                | Periodic-Hill    |
| Density current              | Lock-exchange    |


## Source Download and Compilation

First, make sure you have all the [required dependencies](#required-build-tools-and-external-libraries) installed.
Then, acquire the source code by cloning the git repository:

    git clone git@bitbucket.org:incompact3d/incompact3d.git

(If you are behind a firewall, you may need to use the `https` protocol instead of the `git` protocol:

    git config --global url."https://".insteadOf git@

Be sure to also configure your system to use the appropriate proxy settings, e.g. by setting the `https_proxy` and `http_proxy` variables.)

By default you will be building the latest unstable version of Incompact3d. However, most users should use the most recent stable version of Incompact3d, which is currently the `2.0` series of releases. You can get this version by changing to the Incompact3d directory and running

    git checkout v2.0.0

Now run `make` to build the `Incompact3d` executable. To perform a parallel build, use `make -j N` and supply the maximum number of concurrent processes. (See [Platform Specific Build Notes] for details.)
This takes a while, but only has to be done once. If the defaults in the build do not work for you, and you need to set specific make parameters, you can save them in `Make.user`. The build will automatically check for the existence of `Makefile` and use it if it exists.
Building Incompact3d requires very little of disk space and virtual memory.

**Note:** The build process will fail badly if any of the build directory's parent directories have spaces or other shell meta-characters such as `$` or `:` in their names (this is due to a limitation in GNU make).

Once it is built, you can run the `Incompact3d` executable using its full path in the directory created above (the `Incompact3d` directory).

Now you should be able to run Julia like this:

    mpirun -n 4 ./incompact3d

If everything works correctly, you will see a Incompact3d banner and an interactive prompt into which you can enter expressions for evaluation. (Errors related to libraries might be caused by old, incompatible libraries sitting around in your PATH. In this case, try moving the `Incompact3d` directory earlier in the PATH).

### Updating an existing source tree

If you have previously downloaded `incompact3d` using `git clone`, you can update the
existing source tree using `git pull` rather than starting anew:

    cd incompact3d
    git pull && make

Assuming that you had made no changes to the source tree that will conflict
with upstream updates, these commands will trigger a build to update to the
latest version.

#### General troubleshooting

1. Over time, the base library may accumulate enough changes such that the
   bootstrapping process in building the system image will fail. If this
   happens, the build may fail with an error like

   ```sh
    *** This error is usually fixed by running 'make clean'. If the error persists, try 'make realclean' ***
   ```

   As described, running `make clean && make` is usually sufficient.
   Occasionally, the stronger cleanup done by `make realclean` is needed.

2. New versions of external dependencies may be introduced which may
   occasionally cause conflicts with existing builds of older versions.

   a. To delete existing binaries of `Incompact3d` and all its dependencies,
      delete the `./usr` directory _in the source tree_.

3. If you've updated OS X recently, be sure to run `xcode-select --install` to update the command line tools.

   _To avoid losing work, make sure you know what these commands do before you
   run them. `git` will not be able to undo these changes!_


## Platform-Specific Build Notes

### Linux

#### General

* GCC version 4.7 or later is required to build Julia.
* To use external shared libraries not in the system library search path, set `USE` and `LDFLAGS=-Wl,-rpath,/path/to/dir/contains/libXXX.so` in `Makefile`.
  * Instead of setting `LDFLAGS`, putting the library directory into the environment variable `LD_LIBRARY_PATH` (at both compile and run time) also works.

  - is highly recommended that you remove the limits of the environment (e.g. in `.bash_profile`)

    ulimit -s unlimited
    ulimit -c unlimited

    You must exit and re-login from your terminal for the change to take effect

    ulimit -a


#### Architecture Customization

Incompact3d can be compiled for a non-generic architecture by configuring the `ARCH` Makefile variable. See the appropriate section of `Makefile` for additional customization options, such as `MARCH` and `CPU_TARGET`.

For example, to build for Pentium 4, set `MARCH=pentium4` this will compile the code for the processor specific vectorization level.

You can also set `MARCH=native` for a maximum-performance build customized for the current machine CPU.


#### Ubuntu

In order to compile and execute Incompact3d in the latest Ubuntu version please install the following packages:

    sudo apt install gfortran libopenmpi-dev



#### Fedora/RHEL/CentOS

On RHEL/CentOS 6 systems, the default compiler (`gcc` 4.4) is too old to build Incompact3d.

    sudo dnf install gcc-gfortran

Install or contact your systems administrator to install a more recent version of `gcc`. The [Scientific Linux Developer Toolset](http://linux.web.cern.ch/linux/devtoolset/) works well.

### OS X

You need to have the current Xcode command line utilities installed: run `xcode-select --install` in the terminal.
You will need to rerun this terminal command after each OS X update, otherwise you may run into errors involving missing libraries or headers. You will also need a 64-bit gfortran to compile the code. The gfortran-4.7 (and newer) compilers in Brew, Fink, and MacPorts work for building Incompact3d.

On current systems, we recommend that you install the command line tools as described above. Older systems do not have a separate command line tools package from Apple, and will require a full Xcode install. On these, you will need at least Xcode 4.3.3. In Xcode prior to v5.0, you can alternatively go to Preferences -> Downloads and select the Command Line Utilities.

### External FFT libraries

#### FFTW



### Intel compilers and Math Kernel Library (MKL)

To build Incompact3d using the Intel compilers, and link against
the [MKL] FFTW libraries, first make sure you have a recent version
of the Intel Parallel Studio XE Cluster Edition (version 15 or later).

Qualified students, classroom educators and open-source contributors can obtain a license and download the latest version on https://software.intel.com/en-us/parallel-studio-xe/choose-download/student-linux-fortran

After installation, for a 64-bit architecture, the environment should be set up as follows:

    # bash
    source /path/to/intel/bin/compilervars.sh intel64

The code is compiled with the mpiifort command and the appropriate flags are set in the `Makefile` file.
If you with to compile the code with Intel

    USEICC = 1

## Source Code Organization

The Julia source code is organized as follows:

    contrib/       editor support for Julia source, miscellaneous scripts
    doc/src/manual source for the user manual
    src/           source for Incompact3d code core
    test/          test suites
    ui/            source for various front ends (under development)
