# Sandbox Flow Configuration

> "A Jupyter sandbox environment coupled into the high-order Navier-Stokes solver Xcompact3d", at JupyterCon 2020 ([@JupyterCon](https://twitter.com/JupyterCon)), by [@fschuch](https://github.com/fschuch), [@momba98](https://github.com/momba98), [@filipi](https://github.com/filipi) & J.H. Silvestrini. More details are available [in this link](https://www.fschuch.com/en/publication/2020-jupytercon/).

> "Sandbox flow configuration: A rapid prototyping tool inside XCompact3d", at XCompact3d 2021 Online Showcase Event, by [@fschuch](https://github.com/fschuch). More details are available [in this link](https://www.fschuch.com/en/talk/sandbox-flow-configuration-a-rapid-prototyping-tool-inside-xcompact3d/).

The `module sandbox` ([BC-Sandbox.f90](../../src/BC-Sandbox.f90)) is a new addition to the code.
With this module, the entire initial set-up for any given flow configuration can be imported from external files, including physical and numerical parameters, initial and boundary conditions, and a solid geometry that can be inserted with Immersed Boundary Method (IBM). There is no need to worry about parallel computation in a distributed-memory system, Message Passing Interface, besides coding, compiling, testing, and debugging in Fortran.

The outcome of the presented framework benefits users from different levels:

- For students in computational fluid dynamics, it provides direct hands-on experience and a safe place for practising and learning;
- For advanced users and code developers, it works as a rapid prototyping tool to test concepts and then compare results to validate any future implementations at the numerical solver;
- Furthermore, it is a useful advance in terms of research reproducibility.

The sandbox flow configuration is activated when `itype=12` at the configuration file. Xcompact3d will then read all the necessary arrays for the initial set-up from the disc, in the raw binary format compatible with [2DECOMP&FFT](http://www.2decomp.org/).
These arrays can be provided from any computational tool, the choice is up to the user: Fortran, Matlab, Python, or any other. In this way, it adds no extra dependencies to your workflow.

## Initial set-up

Try it with the Python package [Xcompact3d-toolbox](https://github.com/fschuch/xcompact3d_toolbox) ([see examples online](https://xcompact3d-toolbox.readthedocs.io/en/latest/tutorial.html#sandbox-examples)). 
Everything needed is obtained in a [xarray.Dataset](http://xarray.pydata.org/en/stable/generated/xarray.Dataset.html), including all variables that should be provided to Xcompact3d and the sandbox flow configuration, according to the computational and physical parameters at the configuration file `input.i3d`. See the example:

```python
>>> import xcompact3d_toolbox as x3d
>>> import xcompact3d_toolbox.sandbox
>>> prm = x3d.Parameters(loadfile='input.i3d')
>>> dataset = x3d.sandbox.init_dataset(prm)
>>> dataset
<xarray.Dataset>
Dimensions:       (n: 1, x: 17, y: 17, z: 17)
Coordinates:
  * x             (x) float64 0.0 0.0625 0.125 0.1875 ... 0.875 0.9375 1.0
  * y             (y) float64 0.0 0.0625 0.125 0.1875 ... 0.875 0.9375 1.0
  * z             (z) float64 0.0 0.0625 0.125 0.1875 ... 0.875 0.9375 1.0
  * n             (n) int32 0
Data variables:
    bxx1          (y, z) float64 0.0 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0 0.0
    bxy1          (y, z) float64 0.0 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0 0.0
    bxz1          (y, z) float64 0.0 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0 0.0
    noise_mod_x1  (y, z) float64 0.0 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0 0.0
    bxphi1        (n, y, z) float64 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0 0.0
    byphi1        (n, x, z) float64 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0 0.0
    byphin        (n, x, z) float64 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0 0.0
    ux            (x, y, z) float64 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0 0.0
    uy            (x, y, z) float64 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0 0.0
    uz            (x, y, z) float64 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0 0.0
    phi           (n, x, y, z) float64 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0
>>> # Sequence of operations to set the value for the variables above
>>> dataset.x3d.write(prm)
```

More details are provided bellow, if you want to code your own solution. The examples are in Python, but it should be compatible with any language, as long as the files are written as raw binaries in the same fashion that [2DECOMP&FFT](http://www.2decomp.org/) would do, and the filenames are correct.

### Initial Condition

The subroutine `init_sandbox` ([BC-Sandbox.f90](../../src/BC-Sandbox.f90)) reads the initial values for velocity from three files:

* `./data/ux.bin`, with dimensions (nx, ny, nz);
* `./data/uy.bin`, with dimensions (nx, ny, nz);
* `./data/uz.bin`, with dimensions (nx, ny, nz).

Each scalar field is loaded (if `numscalar > 0`) from the file:

* `./data/phi<n>.bin`, with dimensions (nx, ny, nz), where 1 <= n <= numscalar.

In Python, any of the arrays can be handled with Numpy, just make sure to use `dtype=np.float64` if Xcompact3d was compiled with the flag `-DDOUBLE_PREC`, use `dtype=np.float32` otherwise. See the example:

```python
import numpy

ux = numpy.zeros(
    shape = (nx, ny, nz),
    dtype = numpy.float64    
)
uy = numpy.zeros_like(ux)
uz = numpy.zeros_like(ux)

phi = numpy.zeros(
    shape = (nx, ny, nz, numscalar),
    dtype = numpy.float64    
)

'''
Sequence of operations to set
the initial velocity and scalar concentration
'''

ux.T.tofile('./data/ux.bin')
uy.T.tofile('./data/uy.bin')
uz.T.tofile('./data/uz.bin')

for n in range(numscalar):
    phi[:,:,:,n].T.tofile('./data/phi{}.bin'.format(n+1))
```

### Geometry

The subroutine `geomcomplex_sandbox` ([BC-Sandbox.f90](../../src/BC-Sandbox.f90)) reads the epsilon field describing the solid geometry that is going to be inserted in the computational domain by IBM:

* `./data/geometry/epsilon.bin`, with dimensions (nx, ny, nz).

It is just necessary if `iibm = 1`. See the example:

```python
epsi = numpy.zeros(
    shape = (nx, ny, nz),
    dtype = numpy.float64    
)

'''
Sequence of operations to set
the geometry as 1 inside the object(s)
and 0 in the fluid region
'''

epsi.T.tofile('./data/geometry/epsilon.bin')
```

More information about the solid object(s) is necessary for the alternating direction forcing strategy (if `iibm = 2`), Xcompact3d expects to read it from three files (`obj-x.csv`, `obj-y.csv`, `obj-z.csv`), consider using [`xcompact3d_toolbox.genepsi.gene_epsi_3D`](xcompact3d_toolbox.genepsi.gene_epsi_3D) in this case ([see example](https://xcompact3d-toolbox.readthedocs.io/en/latest/examples/Cylinder.html)).

### Flow Rate Control

A forcing term is included at `module sandbox` in order to keep a constant flow rate when the flow is periodic in the stream-wise direction (`nclx1=nclxn=0`). To do so, the subroutine `flow_rate_control` ([BC-Sandbox.f90](../../src/BC-Sandbox.f90)) applies a customized operator for the volumetric integration, it can be set in the file:

* `'./data/vol_frc.bin`, with dimensions (nx, ny, nz).

Then, the code will compute de stream-wise flow rate as `int = sum(vol_frc * ux)`, and correct the stream-wise velocity as `ux = ux / int`.

### Boundary Condition

The value of the functions must be provided for Dirichlet boundary condition.
At the moment, `module sandbox` only supports reading the velocity BC where `x=0`, of course, if `nclx1=2`. It can be done using three files:

* `./data/bxx1.bin`, with dimensions (ny,nz), for `ux` where `x=0`;
* `./data/bxy1.bin`, with dimensions (ny,nz), for `uy` where `x=0`;
* `./data/bxz1.bin`, with dimensions (ny,nz), for `uz` where `x=0`.

Additionally, a function can be provided for random noise modulation at the inlet, in order to prevent noise in specific regions of the inflow plane. See `inflow` ([BC-Sandbox.f90](../../src/BC-Sandbox.f90)) for more details. It is done with the file:

* `./data/noise_mod_x1.bin`, with dimensions (ny,nz).

In the case of `nclxn=2`, `module sandbox` will apply convective outflow. Any other boundary set as 2 will be considered no-slip.

For the scalar field(s), `n` inflow profiles (where `x=0`) are loaded when `nclxS1=2`:

* `./data/bxphi1<n>.bin`, with dimensions (ny,nz), for `phi<n>` where `x=0`.

A convective outflow is applied for scalar if `nclxSn=2`. Besides, if the settling velocity is zero (`uset[n]=0`), the value of scalar field(s) at the bottom (`y=0`, when `nclyS1=2`) and top (`y=Ly`, when `nclySx=2`) can be set using the files:

* `./data/byphi1<n>`, with dimensions (nx,nz), for `phi<n>` where `y=0`;
* `./data/byphin<n>`, with dimensions (nx,nz), for `phi<n>` where `y=Ly`.

On the other hand, deposit condition is applied when the settling velocity is greater than zero.

See the example:

```python
bxx1 = numpy.zeros(
    shape = (nx, ny),
    dtype = numpy.float64    
)
bxy1 = numpy.zeros_like(bxx1)
bxy1 = numpy.zeros_like(bxx1)
noise_mod = numpy.ones_like(bxx1)

bxphi1 = numpy.zeros(
    shape = (nx, ny, numscalar),
    dtype = numpy.float64    
)
byphi1 = numpy.zeros(
    shape = (nx, nz, numscalar),
    dtype = numpy.float64    
)
byphin = numpy.ones_like(byphi1)

'''
Sequence of operations to set
the BC for velocity and scalar concentration
'''

bxx1.T.tofile('./data/bxx1.bin')
bxy1.T.tofile('./data/bxy1.bin')
bxz1.T.tofile('./data/bxz1.bin')
noise_mod.T.tofile('./data/noise_mod_x1.bin')

for n in range(numscalar):
    bxphi1[:,:,n].T.tofile('./data/bxphi1{}.bin'.format(n+1))
    byphi1[:,:,n].T.tofile('./data/byphi1{}.bin'.format(n+1))
    byphin[:,:,n].T.tofile('./data/byphiz{}.bin'.format(n+1))
```
