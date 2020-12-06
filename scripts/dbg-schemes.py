import matplotlib.pyplot as plt
import numpy as np
import math
import os
import csv

prm = '../BC-dbg-schemes.prm'

class Error(object):
    def __init__(self, dir,jLes,bc,n):
        self.dir = dir
        self.jLes = jLes
        self.bc = bc
        self.n = n
        filename = "schemes_"+dir+str(jLes)+bc+str(n).zfill(4)
        with open(filename, 'r') as f:
            data = list(csv.reader(f, delimiter=' ', skipinitialspace=True,
                                   quoting=csv.QUOTE_NONNUMERIC))
        data = np.transpose(data)
        self.dfdx = np.amax(data[1])
        self.dfdxp = np.amax(data[2])
        self.dfdxx = np.amax(data[3])
        self.dfdxxp = np.amax(data[4])
        del data

results = {}

for bc in ["00", "11", "22"]:

    for i,dir in enumerate(["x", "y", "z"]):
        os.system("sed -i \'/#ncl"+dir+"1/c\\"+bc[0]+"      #ncl"+dir+"1' "+prm)
        os.system("sed -i \'/#ncl"+dir+"n/c\\"+bc[1]+"      #ncl"+dir+"n' "+prm)
        os.system("sed -i \'/#ncl"+dir+"1/c\\"+bc[0]+"      #ncl"+dir+"1' "+prm)
        os.system("sed -i \'/#ncl"+dir+"n/c\\"+bc[1]+"      #ncl"+dir+"n' "+prm)
        os.system("sed -i \'/#ncl"+dir+"1/c\\"+bc[0]+"      #ncl"+dir+"1' "+prm)
        os.system("sed -i \'/#ncl"+dir+"n/c\\"+bc[1]+"      #ncl"+dir+"n' "+prm)

    if bc == "00":
        points = [16, 32, 64, 128, 256]#, 512, 1024, 2048]
    else:
        points = [17, 33, 65, 129, 257]#, 513, 1025, 2049]
    for n in points:
        for i,dir in enumerate(["x"]):
            os.system("sed -i \'/#nx/c\\"+str(points[0])+"      #nx' "+prm)
            os.system("sed -i \'/#ny/c\\"+str(points[0])+"      #ny' "+prm)
            os.system("sed -i \'/#nz/c\\"+str(points[0])+"      #nz' "+prm)
            os.system("sed -i \'/#n"+dir+"/c\\"+str(n)+"      #n"+dir+"' "+prm)
            for j,jLes in enumerate([0,1]):
                os.system("sed -i \'/#jLES/c\\"+str(jLes)+"      #jLES' "+prm)
                os.system("mpirun -n 1 ../incompact3d")

                filename = "schemes_"+dir+str(jLes)+bc+str(n).zfill(4)
                results[filename] = Error(dir,jLes,bc,n)
                print(filename)


plt.figure(figsize=(16,9), dpi=120)

color = {}
color["00"] = 'C0'
color["11"] = 'C1'
color["22"] = 'C0'

for key, data in results.items():
    plt.scatter(data.n,data.dfdx, color=color[data.bc])
    plt.scatter(data.n,data.dfdxp, color=color[data.bc])
    plt.scatter(data.n,data.dfdxx, color=color[data.bc])
    plt.scatter(data.n,data.dfdxxp, color=color[data.bc])


error_n = np.arange(10, 1000, .5)

error_y = (error_n)**(-6)
plt.plot(error_n,error_y*10000)
error_y = (error_n)**(-6)
plt.plot(error_n,error_y*10000000)
error_y = (error_n)**(-3)
plt.plot(error_n,error_y*1000)
error_y = (error_n)**(-3)
plt.plot(error_n,error_y*400000)

plt.xlabel(r'$n_1$')
plt.ylabel('Error')

plt.xscale('log')
plt.yscale('log')

plt.xlim(10,1000)
plt.ylim(.00000000000001,100)

plt.savefig('max_error.png', format='png')
plt.close('all')
