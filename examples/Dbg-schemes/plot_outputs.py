# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

A=np.genfromtxt('/home/gdeskos/xcompact3d_developers/examples/Dbg-schemes/filter_y0110064',delimiter='')
y=A[:,0]
ftheo=A[:,1]
fcomp=A[:,2]
ftheop=A[:,3]
fcompp=A[:,4]

plt.plot(y,ftheo,'-b')
plt.plot(y,fcomp,'-ok')
plt.show()