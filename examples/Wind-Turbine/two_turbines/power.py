import numpy as np
file1     = 1       # .perf filenumber at beginning of averaging
nfiles    = 500     # .perf filenumber at end of averaging
itwriting = 200     # iturboutput
dt        = 0.0002  # time step
nturbines = 2       # number of turbines
# Physical time of averaging
time = np.arange(file1, nfiles+1, 1)*dt*itwriting
# Compute thrust and power from .perf files
thrust = np.zeros((nturbines, nfiles+1-file1))
power  = np.zeros((nturbines, nfiles+1-file1))
for iturb in range(1, nturbines+1):
    for ifile in range(file1, nfiles+1):
        data = np.loadtxt("./%d_NTNU_HATT_%d.perf" % (ifile, iturb), delimiter=',', skiprows=2)
        thrust[iturb-1, ifile-file1] = data[11]
        power[iturb-1, ifile-file1]  = data[13]
header = 'Time [s], Power_1 [W], Power_2 [W], Total Power [W]'
np.savetxt("power.txt", np.column_stack((time, power[0,:], power[1,:], power[0,:]+power[1,:])), delimiter=' ')
