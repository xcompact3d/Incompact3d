import matplotlib.pyplot as plt
import numpy as np

from numpy import *
from scipy.signal import lfilter
from scipy.signal import filtfilt
from scipy.signal import *


# --------------------------- Read .dat File --------------------------- # 


#itime, drag = loadtxt('aerof666',unpack=True, usecols=[0,1])
itime1, drag1 = loadtxt('aerof1',unpack=True, usecols=[0,1])
#itime2, drag2 = loadtxt('aerof2',unpack=True, usecols=[0,1])


fig= plt.figure()
ax2= fig.add_subplot(111)

#ax2.plot(itime,drag,ls='-',color='blue',linewidth=2.0, label='Original')
ax2.plot(itime1,drag1,ls='-',color='red',linewidth=2.0, label='Mov1')


ax2.set(ylim=[1.0,1.25])
#ax2.set(xlim=[195,207.5])
plt.legend(loc='upper right')

plt.savefig('drag2.png', dpi=600)

