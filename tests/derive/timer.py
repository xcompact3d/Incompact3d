#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Plot the histogram of the timer for the given files.')
parser.add_argument('files', nargs='+', help='list of .dat files')
args = parser.parse_args()

for filename in args.files:
   #
   tmp=np.loadtxt(filename)
   #
   # Keep log10 of the computation time
   array=np.zeros(tmp.shape[0])
   for i in range(tmp.shape[0]):
      array[i] = np.log10(np.float(tmp[i][-1]))
   #
   # Build bins
   nbins = 15
   bins = [min(array) + i*(max(array)-min(array))/(nbins-1) for i in range(nbins)]
   #
   # Show histogram
   tmpn, tmpbins, tmppatches = plt.hist(array, bins=bins)
   plt.title(filename + " : " + np.str(10**(tmpbins[np.where(tmpn == tmpn.max())][0])) + " seconds")
   plt.ylabel("Number of runs")
   plt.xlabel("log10 of the CPU time")
   plt.savefig(filename[:-4] + ".png")
   plt.show()

