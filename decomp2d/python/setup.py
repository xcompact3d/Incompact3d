import os
from numpy.distutils.core import Extension

import decomp2d_options

VERSION="0.0.0"
REQUIREMENTS=[
    "numpy",
    "mpi4py"
]

X3DDIR="/Users/paulbartholomew/src/fortran/Xcompact3d"
D2DDIR=os.path.join(X3DDIR, "decomp2d")

d2dext = Extension(
    name = "decomp2d",
    sources = [
        "2decomp4py.f90"
    ],
    include_dirs = [
        os.path.join(D2DDIR, "python"),
        D2DDIR,
        X3DDIR
    ],
    define_macros = decomp2d_options.define_macros,
    extra_objects = [
      os.path.join(D2DDIR, "decomp_2d.o")  
    ],
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
