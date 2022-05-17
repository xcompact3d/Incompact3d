import os
from numpy.distutils.core import Extension

import decomp2d_options
import decomp2d_io_options

VERSION="0.0.0"
REQUIREMENTS=[
    "numpy",
    "mpi4py"
]

X3DDIR="/home/paul/src/Xcompact3d/Xcompact3d"
D2DDIR=os.path.join(X3DDIR, "decomp2d")
ADIOS2INC="/home/paul/opt/adios2/v2.8.0/include/adios2/fortran"

macros = decomp2d_options.define_macros
macros = macros + decomp2d_io_options.define_macros
libs = []
lib_dirs = []
rlib_dirs = []

with_adios2 = False
for m in macros:
    if "ADIOS2" in m:
        with_adios2 = True

if with_adios2:
    libs = libs + [
        "adios2_fortran_mpi",
        "adios2_fortran"
    ]
    lib_dirs = lib_dirs + [
        "/home/paul/opt/adios2/v2.8.0/lib"
    ]
    rlib_dirs = rlib_dirs + [
        "/home/paul/opt/adios2/v2.8.0/lib"
    ]

d2dext = Extension(
    name = "decomp2d",
    sources = [
        "2decomp4py.f90"
    ],
    include_dirs = [
        os.path.join(D2DDIR, "python"),
        D2DDIR,
        X3DDIR,
        ADIOS2INC
    ],
    define_macros = macros,
    extra_objects = [
        os.path.join(D2DDIR, "decomp_2d.o"),  
        os.path.join(D2DDIR, "io.o")  
    ],
    libraries = libs,
    library_dirs = lib_dirs,
    runtime_library_dirs = rlib_dirs,
    extra_f90_compile_args = [
        "-cpp"
    ]
)

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(
        name = "decomp_2d4py",
        version = VERSION,
        author = "Xcompact3d project",
        description = "decomp2d&fft Python interface",
        license = "Apache 2.0",
        ext_modules = [d2dext],
        install_requires = REQUIREMENTS
    )
