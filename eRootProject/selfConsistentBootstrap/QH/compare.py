from scipy.io import netcdf_file
import numpy as np
import matplotlib.pyplot as plt

def makeGrids(ns):
    fullgrid = np.linspace(0,1,num=ns)
    delta = (fullgrid[1] - fullgrid[0]) / 2
    halfgrid = fullgrid - delta
    halfgrid[0] = 0

    return [halfgrid, fullgrid]

oldF = netcdf_file('wout_20220218-01-021_QH_A6.5_n0_2.2_T0_10_highResVmecForBestFrom020.nc', mode='r', mmap=False)
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
