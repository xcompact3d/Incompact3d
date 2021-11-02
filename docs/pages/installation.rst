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
