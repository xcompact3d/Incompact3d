============
Installation
============

Xcompact3d supports both traditional Makefile and `cmake` based builds.
In both cases, the only requirements is a Fortran 90-compatible compiler and working `MPI`
installation.

--------
Makefile
--------

After obtaining the Xcompact3d source `cd` into the project root directory and run `make` to build.
By default this will use the `GNU` compiler toolchain, alternative compilers can be used by setting
`CMP`, *e.g.*
``
make CMP=intel
``
permissible values for `CMP` can be found by inspecting the Makefile.
On successful completion you will have an executable `xcompact3d`.

**N.B.** The Makefile currently does not support parallel builds *i.e.* `make -jN`.

^^^^^^^^^^^^^^^^^^^^^^^^
Building on Cray systems
^^^^^^^^^^^^^^^^^^^^^^^^

The Makefile supports building on Cray systems using the Cray compiler with `make CMP=cray`, if
using another environment on a Cray machine (such as `PrgEnv-gnu`) then you should use that
compiler's setting for `CMP` and set `FC=ftn`:
``
# After loading the PrgEnv-gnu environment

make CMP=gnu FC=ftn
``

-----
CMake
-----

To use `cmake` the recommended approach is to create a `build/` directory and configure the build
from there, for example creating `build/` in the Xcompact3d project root directory
``
mkdir build
cd build

cmake ../
``
After which further customisation of the build can be achieved by running `ccmake .` and setting
variables to their desired values (be sure in particular to check the value of
`CMAKE_INSTALL_PREFIX`).

Once the build has been configured run `make` to compile, followed by `make install` which will
install the `xcompact3d` executable to `${CMAKE_INSTALL_PREFIX/bin/}`, you can optionally run tests
on the build by executing `make test`.

---------------
Enabling ADIOS2
---------------

An alternative I/O backend using ADIOS2 has been added to the 2DECOMP&FFT library distributed with
Xcompact3d and can be selected at compile time, following the above build instructions will default
to the original `MPI-IO` backend.

^^^^^^^^^^^^
Dependencies
^^^^^^^^^^^^

Enabling ADIOS2 requires the ADIOS2 library, optionally HDF5 files can be written through the ADIOS2
interface (see runtime configuration).
ADIOS2 (and HDF5) may be available through an HPC machine's module system, if not, and for local
testing and development purposes, the installation process for both is summarised below.

The HDF5 source can be obtained from https://github.com/HDFGroup/hdf5, Xcompact3d with ADIOS2 has
been tested with v1.12.0:
``
git clone git@github.com:HDFGroup/hdf5.git
cd hdf5
git checkout hdf5-1_12_0
``
(release tarfiles are also available).

HDF5 is configured as
``
./configure --prefix=${HDF5_DIR} --enable-parallel --enable-shared --enable-fortran CC=mpicc CXX=mpicxx FC=mpif90
``
where `HDF5_DIR` is the desired install location.
For production use it may be worth exploring the `--enable-build-mode=production` option and other
suggestings in the readmes under `release_docs/`.

After configuring build and install with
``
make
make install
``
this will build and install hdf5 to `${HDF5_DIR}` - check for the presence of `bin/`, `lib/`, etc.
You might also want to add `${HDF5_DIR}` to your path, it contains useful utilities such as `h5dump`
for inspecting hdf5 files.

**N.B.** package manager installations (*e.g.* using `apt-get`) may not be build with
`--enable-parallel` and are therefore unsuitable here.

ADIOS2 can similarly be obtained via git
``
git clone git@github.com:ornladios/ADIOS2.git
cd ADIOS2
git checkout v2.7.1
``
The recommendation is to build with `cmake`:
``
mkdir build
cd build
cmake ../
``
and use `ccmake .` to configure, in particular ensuring the Fortran bindings are enabled, you can
also enable HDF5 and set the path to your HDF5 installation if it is in a non-standard location.
After configuring build and install with
``
make
make install
``
note that you can control the installation location by passing
`-DCMAKE_INSTALL_PREFIX=${ADIOS2_DIR}` to `cmake` or by setting the variable when configuring with
`ccmake`.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Building Xcompact3d with ADIOS2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To build Xcompact3d with ADIOS2 as the I/O backend either use `make`:
``
make clean
make IO=adios2 ADIOS2DIR=${ADIOS2_DIR}
``
or with `cmake` (from the `build/` directory) use `ccmake .` to turn ADIOS2 `ON` and set `ADIOS2DIR`
to `${ADIOS2_DIR}/lib/cmake/adios2` (note some installations use `lib64` in place of `lib`),
followed by `make && make install` as above.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Running Xcompact3d with ADIOS2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Running an ADIOS2-enabled build of Xcompact3d requires an `adios2_config.xml` file to provide the
runtime configuration for ADIOS2, an example can be found in the `Taylor-Green-Vortex` example
directory.
With this it is possible to switch the "engine" for example to change from writing ADIOS2-native
`.bp4` output to HDF5, and various other aspects of the I/O can be controlled at runtime - see the
ADIOS2 documentation for possibilities.

