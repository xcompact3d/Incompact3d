# Import necessary libraries
import xcompact3d_toolbox as x3d
import numpy as np
import xarray as xr
import math

# Initialise parameters
x3d.param["mytype"] = np.float64

# Load parameters
prm  = x3d.Parameters(loadfile='input.i3d')
epsi = x3d.init_epsi(prm, dask=True)

# Read the STL file
for key in epsi.keys():
    print(key)
    epsi[key] = epsi[key].geo.from_stl(filename="sphere.stl",
                                       origin=dict(x=2.0,y=2.0,z=2.0), # locates the sphere centre at (2.5, 2.5, 2.5)
                                       scale=1.0/200.0, # Sphere has a diameter of 200, so rescale to a diameter of 1
                                       # rotate=dict(axis=[0,0,1],theta=math.radians(270)), # rotate the geometry
                                       user_tol=1e-3)
dataset = x3d.gene_epsi_3d(epsi, prm)

if prm.iibm >= 2:
    prm.nobjmax = dataset.obj.size
    print(prm.nobjmax) # If using iibm = 2, set nobjmax in input.i3d to the value of this print

ds  = x3d.init_dataset(prm)

# Set boundary conditions
for key in 'bxx1 bxy1 bxz1'.split():
    print(ds[key].attrs['name'])
    ds[key] *= 0.0
    if key=='bxx1':
        ds[key] += 1.0

# Set initial conditions
for key in 'ux uy uz'.split():
    print(ds[key].attrs['name'])
    ds[key] *= 0.0
    if key=='ux':
        ds[key] += 1.0

# Write the fields
prm.dataset.write(ds)
