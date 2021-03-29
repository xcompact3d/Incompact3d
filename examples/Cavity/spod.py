#!/usr/bin/env python3

#
# Spectral POD (pip3 install pyspod)
#

import os
import numpy as np
# Import library specific modules
from pyspod.spod_low_storage import SPOD_low_storage
from pyspod.spod_low_ram     import SPOD_low_ram
from pyspod.spod_streaming   import SPOD_streaming

#
# Read input.i3d file
#
nx = 0
ny = 0
nz = 0
lx = 0.
ly = 0.
lz = 0.
dt = 0.
iout = 0
with open("./input.i3d") as file:
    for line in file.readlines():
        #
        if line[:2] == "nx":
            nx = np.int(line.split("=")[1].split()[0])
        #
        if line[:2] == "ny":
            ny = np.int(line.split("=")[1].split()[0])
        #
        if line[:2] == "nz":
            nz = np.int(line.split("=")[1].split()[0])
        #
        if line[:3] == "xlx":
            lx = np.float(line.split("=")[1].split()[0])
        #
        if line[:3] == "yly":
            ly = np.float(line.split("=")[1].split()[0])
        #
        if line[:3] == "zlz":
            lz = np.float(line.split("=")[1].split()[0])
        #
        if line[:2] == "dt":
            dt = np.float(line.split("=")[1].split()[0])
        #
        if line[:7] == "ioutput":
            iout = np.int(line.split("=")[1].split()[0])

#
# Snapshots from ./data
#
file_t = []
file_ux = []
file_uy = []
for file in os.listdir("data"):
    if file[:5] == "phi1-":
        file_t.append(file)
    if file[:3] == "ux-":
        file_ux.append(file)
    if file[:3] == "uy-":
        file_uy.append(file)

file_t.sort()
file_ux.sort()
file_uy.sort()

#
# Load data
#
files = [file_t, file_ux, file_uy]
variables = ['T', 'ux', 'uy', 'Ec']
dataND = np.zeros((len(file_t), nx, ny, len(variables)))

#
# Read from files first
#
for ivar, flist in enumerate(files):
    for it, ff in enumerate(flist):
        data2D = np.fromfile(open(os.path.join("data",ff), 'rb'), np.float)
    for i in range(nx):
        for j in range(ny):
            dataND[it, i, j, ivar] = data2D[i + j * nx]

#
# Computing after
#    0 => T
#    1 => ux
#    2 => uy
#    3 => Ec
#
dataND[:,:,:,3] = (dataND[:,:,:,1]**2 + dataND[:,:,:,2]**2)/2.

#
# SPOD parameters into a dictionary
#
params = dict()
# -- required parameters
params['dt'          ] = 1                    # data time-sampling
params['nt'          ] = dataND.shape[0]      # number of time snapshots (we consider all data)
params['xdim'        ] = 2                    # number of spatial dimensions (longitude and latitude)
params['nv'          ] = dataND.shape[-1]     # number of variables
params['n_FFT'       ] = 128                  # length of FFT blocks (100 time-snapshots)
params['n_freq'      ] = params['n_FFT'] / 2 + 1             # number of frequencies
params['n_overlap'   ] = np.ceil(params['n_FFT'] * 50 / 100) # dimension block overlap region
params['mean'        ] = 'blockwise'          # type of mean to subtract to the data
params['normalize'   ] = False                # normalization of weights by data variance
params['savedir'     ] = os.path.join(os.getcwd(), 'out') # folder where to save results
# -- optional parameters
params['weights']      = None # if set to None, no weighting (if not specified, Default is None)
params['savefreqs'   ] = np.arange(0,params['n_freq']) # frequencies to be saved
params['n_modes_save'] = 9      # modes to be saved
params['normvar'     ] = False  # normalize data by data variance
params['conf_level'  ] = 0.95   # calculate confidence level
params['savefft'     ] = False  # save FFT blocks to reuse them in the future (saves time)

#
# Run the SPOD
#
spod_ls = SPOD_low_storage(dataND, params=params, data_handler=False, variables=variables)
spod_ls.fit()

#
# Plot some data
#
spod_ls.plot_2D_data(time_idx=[1,2], vars_idx=[0,1,2,3])
spod_ls.plot_data_tracers(coords_list=[(5,2)], time_limits=[0,dataND.shape[0]], vars_idx=[0,1,2,3])

#
# Plot more data
#
T_approx = 3.3 # approximate period
freq = spod_ls.freq
freq_found, freq_idx = spod_ls.find_nearest_freq(freq_required=1/T_approx, freq=freq)
modes_at_freq = spod_ls.get_modes_at_freq(freq_idx=freq_idx)
spod_ls.plot_eigs()
spod_ls.plot_eigs_vs_frequency(freq=freq)
spod_ls.plot_eigs_vs_period(freq=freq)
spod_ls.plot_2D_modes_at_frequency(freq_required=freq_found, freq=freq, x1=np.linspace(0.,lx,nx), x2=np.linspace(0.,ly,ny), modes_idx=[0,1], vars_idx=[0,1,2,3])
