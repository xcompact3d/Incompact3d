# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

# Test 11
A11=np.genfromtxt('/home/gdeskos/xcompact3d_developers/examples/Dbg-schemes/filter_z0210064',delimiter='')
y=A11[:,0]
ftheo=A11[:,1]
fcomp=A11[:,2]
ftheop=A11[:,3]
fcompp=A11[:,4]

plt.plot(y,ftheop,'-bs')
plt.plot(y,fcompp,'-ok')
plt.plot(y,ftheo,'-gs')
plt.plot(y,fcomp,'-or')
plt.show()