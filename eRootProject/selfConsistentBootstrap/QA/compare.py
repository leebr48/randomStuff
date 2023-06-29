from scipy.io import netcdf_file
import numpy as np
import matplotlib.pyplot as plt

def makeGrids(ns):
    fullgrid = np.linspace(0,1,num=ns)
    delta = (fullgrid[1] - fullgrid[0]) / 2
    halfgrid = fullgrid - delta
    halfgrid[0] = 0

    return [halfgrid, fullgrid]

oldF = netcdf_file('wout_QA_beta0p025_iota0p42_dreopt_HIGHERRES_2022-04-15.nc', mode='r', mmap=False)
oldIotas = oldF.variables['iotas'][()]
oldIotaf = oldF.variables['iotaf'][()]
oldNs = oldF.variables['ns'][()]
[oldHalf, oldFull] = makeGrids(oldNs)

newF = netcdf_file('wout_yourTry.nc', mode='r', mmap=False)
newIotas = newF.variables['iotas'][()]
newIotaf = newF.variables['iotaf'][()]
newNs = newF.variables['ns'][()]
[newHalf, newFull] = makeGrids(newNs)

plt.plot(oldHalf, oldIotas, label='old half')
plt.plot(newHalf, newIotas, label='new half')

plt.figure()

plt.plot(oldFull, oldIotaf, label='old full')
plt.plot(newFull, newIotaf, label='new full')

plt.show()
