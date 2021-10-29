#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Plot the L1, L2 and Linf error for the given files.')
parser.add_argument('files', nargs='+', help='list of .dat files')
args = parser.parse_args()

nn = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096]

for filename in args.files:
   #
   tmp=np.loadtxt(filename)
   #
   # Keep log10 of the errors
   errl1=np.zeros(tmp.shape[0])
   errl2=np.zeros(tmp.shape[0])
   errlinf=np.zeros(tmp.shape[0])
   for i in range(tmp.shape[0]):
      errl1[i] = np.log10(np.float(tmp[i][0]))
      errl2[i] = np.log10(np.float(tmp[i][1]))
      errlinf[i] = np.log10(np.float(tmp[i][2]))
   #
   # Plot errors in log scale
   #
   fig, axs = plt.subplots(1)
   fig.suptitle("Convergence for " + filename)
   tmp = 0.1*np.min([10**errl1[0], 10**errl2[0], 10**errlinf[0]])
   axs.plot(np.log2([nn[0], nn[3]]), np.log10([tmp, tmp*(nn[0]/nn[3])**6]), label="Od6")
   axs.plot(np.log2(nn), errl1, label="L1")
   axs.plot(np.log2(nn), errl2, label="L2")
   axs.plot(np.log2(nn), errlinf, label="Linf")
   axs.set_ylabel("Error, log10")
   axs.set_xlabel("Number of nodes, log2")
   axs.legend()
   fig.savefig(filename[:-4] + ".png")

