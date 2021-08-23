#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# Read channel.dat
tmp = np.loadtxt("channel.dat", dtype=str, delimiter="\n", skiprows=1)

# Number of columns (<u>, <v>, <w>, <T>, <u'²>, <v'²>, <w'²>, <T'²>)
ncol = len(tmp[0].split())

# Create numpy array
array=np.zeros((tmp.shape[0], ncol))
for i in range(tmp.shape[0]):
    array[i,:] = np.array([np.float(item) for item in tmp[i].split()])

# Plot
for i in range(ncol):
    plt.plot(array[:,i])
    plt.title("Column " + np.str(i))
    plt.show()
