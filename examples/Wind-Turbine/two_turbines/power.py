import numpy as np

n_turbines = 2

power = []
for iturb in range(n_turbines):
    data = np.loadtxt(f"./NTNU_HATT_{iturb+1}.perf", delimiter=',', skiprows=2)
    power.append(data[:, 15])
    time = data[:, 1]

total_power = np.zeros_like(power[0])
for iturb in range(n_turbines):
    total_power += power[iturb]

header = 'Time [s], '
for iturb in range(n_turbines):
    header += f'Power_{iturb+1} [W], '
header += 'Total Power [W]'

np.savetxt("power.txt", np.column_stack((time, np.array(power).T, total_power)), delimiter=' ', header=header)
