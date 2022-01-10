import os
from numpy.distutils.core import Extension

VERSION="0.0.0"
REQUIREMENTS=[
    "numpy",
    "mpi4py"
]

d2dext = Extension(
    name = "decomp2d",
    sources = [
        "2decomp4py.f90"
    ],
    include_dirs = [
        ".",
        "../",
        "../../"
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
