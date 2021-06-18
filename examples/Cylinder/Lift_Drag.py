import matplotlib.pyplot as plt
import numpy as np

from numpy import *
from scipy.signal import lfilter
from scipy.signal import filtfilt
from scipy.signal import *


# --------------------------- Read .dat File --------------------------- # 


itime, drag, lift = loadtxt('aerof1',unpack=True, usecols=[0,1,2])
itime2, drag2,lift2 = loadtxt('aerof1-0090000',unpack=True, usecols=[0,1,2])

itime3 = itime2
drag3 = drag2
lift3 = lift2


# Set font to Times New Roman
plt.rcParams['font.family'] = "Times New Roman"


fig= plt.figure()
ax2= fig.add_subplot(111)

ax2.plot(itime2,drag2,ls='-',color='blue',linewidth=2.0, label='Stationary ')

ax2.set(ylim=[0.0,2.5])
#ax2.set(xlim=[0,45000])

plt.legend(loc='upper right')

plt.savefig('drag.png', dpi=600)





#Dt = 0.00032


#fil3_3 = savgol_filter(drag3, 305, 1)
#fil3_3 = savgol_filter(fil3_3, 305, 1)

#print('-------------------------------------------------------------------------')
#print('  Mean Cd for Stationary           :  ', np.mean(drag[150000:250000]))
#print('  Mean Cd for Moving               :  ', np.mean(drag2[10000:20000]))
#print('  Mean Cd for filtered Moving      :  ', np.mean(fil3_3[10000:20000]))
#print('-------------------------------------------------------------------------')
#print('  Mean Cd Amplitude for Stationary :  ', np.amax(drag[150000:250000]) -np.amin(drag[150000:250000]))
#print('  Mean Cd Amplitude for Stationary :  ', np.amax(fil3_3[10000:20000]) -np.amin(fil3_3[10000:20000]))
#print('-------------------------------------------------------------------------')
#fig= plt.figure()
#ax2= fig.add_subplot(111)

#ax2.plot(itime[150000:250000]-150000,drag[150000:250000],ls='-',color='blue',linewidth=2.0, label='Stationary ')
#ax2.plot(itime3-1700,fil3_3,ls='-',color='black',linewidth=2.0, label='Savitzky–Golay Filter k=65')

#ax2.set(ylim=[1.3,1.6])
#ax2.set(xlim=[0,45000])

#plt.legend(loc='upper right')

#plt.savefig('drag2.png', dpi=600)

################################################################################################################
##
##                                       Combine Plots
##
################################################################################################################

#fil3_4 = savgol_filter(lift3, 6005, 1)

#fig= plt.figure()
#ax2= fig.add_subplot(111)

#ax2.plot(itime[150000:250000]-156100,drag[150000:250000],ls='-.',color='black',linewidth=2.0, label='$\mathregular{C_D}$, ADR-S')
#ax2.plot(itime3-500,fil3_3,ls='-',color='tab:orange',linewidth=2.0, label='$\mathregular{C_D}$, ADR-M')

#ax2.plot(itime[150000:250000]-152700-3500,lift[150000:250000],ls='-',color='black',linewidth=2.0, label='$\mathregular{C_L}$, ADR-S')
#ax2.plot(itime3-3500,fil3_4,ls='--',color='tab:orange',linewidth=2.0, label='$\mathregular{C_L}$, ADR-M')

#ax2.set(ylim=[-1.2,2.5])
#ax2.set(xlim=[5000,43750])


#plt.ylabel('$\mathregular{C_D, C_L}$', fontsize=14)
#plt.xlabel('$\mathregular{t}$', fontsize=14)


## Set the ticks and ticklabels for all axes
#labels = ['90', '-1.0', '0.0', '1.0', '2.0']
#plt.setp(ax2, xticks=[5000, 15000, 25000, 35000, 45000], xticklabels=labels)


#plt.legend(loc='upper left')

#plt.savefig('Lift_Drag.pdf', dpi=600, bbox_inches = "tight")
################################################################################################################
##
##                                       Plots for Abstract
##
################################################################################################################
#fig= plt.figure()
#ax2= fig.add_subplot(111)

#ax2.plot(itime[150000:250000]-156100,drag[150000:250000],ls='-.',color='black',linewidth=2.0, label='$\mathregular{C_D}$, ADR-S')
#ax2.plot(itime3-500,fil3_3,ls='-',color='tab:red',linewidth=2.0, label='$\mathregular{C_D}$, ADR-M')

#ax2.plot(itime[150000:250000]-152700-3500,lift[150000:250000],ls='-',color='black',linewidth=2.0, label='$\mathregular{C_L}$, ADR-S')
#ax2.plot(itime3-3500,fil3_4,ls='--',color='tab:red',linewidth=2.0, label='$\mathregular{C_L}$, ADR-M')

#ax2.set(ylim=[-1.2,2.5])
#ax2.set(xlim=[5000,43750])


#plt.ylabel('$\mathregular{C_D, C_L}$', fontsize=14)
#plt.xlabel('$\mathregular{t}$', fontsize=14)


## Set the ticks and ticklabels for all axes
#labels = ['90', '-1.0', '0.0', '1.0', '2.0']
#plt.setp(ax2, xticks=[5000, 15000, 25000, 35000, 45000], xticklabels=labels)


#plt.legend(loc='upper left')

#plt.savefig('Lift_Drag_Abstract.pdf', dpi=600, bbox_inches = "tight")

#fig= plt.figure()
#ax2= fig.add_subplot(111)

#ax2.plot(itime[150000:250000]-156100,drag[150000:250000],ls='-.',color='black',linewidth=2.0, label='$\mathregular{C_D}$, ADR-S')
#ax2.plot(itime3-500,fil3_3,ls='-',color='turquoise',linewidth=2.0, label='$\mathregular{C_D}$, ADR-M')

#ax2.plot(itime[150000:250000]-152700-3500,lift[150000:250000],ls='-',color='black',linewidth=2.0, label='$\mathregular{C_L}$, ADR-S')
#ax2.plot(itime3-3500,fil3_4,ls='--',color='turquoise',linewidth=2.0, label='$\mathregular{C_L}$, ADR-M')

#ax2.set(ylim=[-1.2,2.5])
#ax2.set(xlim=[5000,43750])


#plt.ylabel('$\mathregular{C_D, C_L}$', fontsize=14)
#plt.xlabel('$\mathregular{t}$', fontsize=14)


## Set the ticks and ticklabels for all axes
#labels = ['90', '-1.0', '0.0', '1.0', '2.0']
#plt.setp(ax2, xticks=[5000, 15000, 25000, 35000, 45000], xticklabels=labels)


#plt.legend(loc='upper left')

#plt.savefig('Lift_Drag_Abstract2.pdf', dpi=600, bbox_inches = "tight")





################################################################################################################
##
##                                       Get Strouhal Number
##
################################################################################################################
#import scipy as sp
#from scipy import fftpack


## Get Strouhal Number for the Stationary Case:
#Signal = lift[150000:300000]

#Dt = 0.00032     # Time-Step

#T = len(Signal)*Dt     # Time Unitis (total): itr*Dt
#Fs = 1/Dt              # Sampling frequency
#x = Signal
#N = x.size

## Discrete Fourier Transform
#X = np.fft.fft(x)
#X_db = 20*np.log10(2*np.abs(X)/N)
#f = np.arange(0, N)*Fs/N

## Get size of the array 
#Size = int(len(f)/2)
## Get Index of the max
#index = np.argmax(X_db[0:Size])

## Strouhal Number:
#print('       *    *    Strouhal Number    *    *       ')

#print('------------------------------------------------')
#print('  The Strouhal Number for the Stationary is :  ', f[index])
##print('------------------------------------------------')


## Get Strouhal Number for the Moving Case:
#Signal = lift3[10000:70000]

#Dt = 0.00032     # Time-Step

#T = len(Signal)*Dt     # Time Unitis (total): itr*Dt
#Fs = 1/Dt              # Sampling frequency
#x = Signal
#N = x.size

## Discrete Fourier Transform
#X = np.fft.fft(x)
#X_db = 20*np.log10(2*np.abs(X)/N)
#f = np.arange(0, N)*Fs/N

## Get size of the array 
#Size = int(len(f)/2)
## Get Index of the max
#index = np.argmax(X_db[0:Size])

## Strouhal Number:
##print('------------------------------------------------')
#print('  The Strouhal Number for the Moving is     :  ', f[index])
#print('------------------------------------------------')















## Plot FFT of signal

#fig= plt.figure()
#ax2= fig.add_subplot(111)


#ax2.plot(f, X_db)

##ax2.set(ylim=[0.5,1.75])
##ax2.set(xlim=[-0.2,0.75])

#plt.savefig('FFT.png', dpi=600)


##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################

##temp = lift3[0:70000]

##N = len(temp)

##X = np.fft.fft(x)
##X_db = 20*np.log10(2*np.abs(X)/N)
##f = np.arange(0, N)*Fs/N


##temp_fft = sp.fftpack.fft(temp)
##temp_psd = np.abs(temp_fft) ** 2
##fftfreq = sp.fftpack.fftfreq(len(temp_psd), 1. / (Dt))

##i = fftfreq > 0




##fig= plt.figure()
##ax2= fig.add_subplot(111)

###ax2.plot(itime,lift,ls='-',color='black',linewidth=1.5, label='Moving ')
##ax2.plot(fftfreq[i], 10 * np.log10(temp_psd[i]))

###ax2.set(ylim=[-2,2])
##ax2.set(xlim=[-2,20])

##plt.legend(loc='upper right')

##plt.savefig('Test1.png', dpi=600)




################################################################################################################
##
##                                                 Filter with FFT
##
################################################################################################################


#temp = X_db.copy()
#temp[np.abs(X_db) > 10] = 0

#X1 = (N/2)*(10**(temp/20))

#temp_slow = np.real(np.fft.ifft(X1))



#fig= plt.figure()
#ax2= fig.add_subplot(111)

##ax2.plot(itime,lift,ls='-',color='black',linewidth=1.5, label='Moving ')
#ax2.plot(itime3[0:70000],temp_slow,ls='-.',color='red',linewidth=1.5, label='Moving ')
#ax2.plot(itime[1000000:1040000]-1000000,lift[1000000:1040000],ls='-',color='blue',linewidth=2.0, label='Stationary ')

#ax2.set(ylim=[-2,2])
#ax2.set(xlim=[0,70000])

#plt.legend(loc='upper right')

#plt.savefig('lift2.png', dpi=600)

















