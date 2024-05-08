# This script plots epseff and DMerc for the initial and final configurations. These are calculated on *ALL SURFACES POSSIBLE* for each configuration.

#Inputs
rhoMin = 0
rhoMax = 1
epseffMax = 5
DMercYMin = -0.001
DMercYMax = 0.01
axisFontSize = 24
legendFontSize = 13
xSizeInches = 8
ySizeInches = 6.0
fileExt = 'pdf'
dpi = 600

########################################################################################

# Import modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import sys
from scipy.io import netcdf_file
sys.path.append('/home/blee/src/stelloptPlusSfincs/src')
from dataProc import createVMECGrids

plt.rc('font', size=axisFontSize)
plt.rc('legend', fontsize=legendFontSize)

# Define handy function
def my_formatter(x, pos):
    if x == 0:
        return str(int(0))
    elif x == 1:
        return str(int(1))
    else:
        return '{0:.1f}'.format(x)

formatter = FuncFormatter(my_formatter)

# Load data
lee1_epseff = np.loadtxt('lee1_epseff.txt')
lee2_epseff = np.loadtxt('lee2_epseff.txt')
lee3_epseff = np.loadtxt('lee3_epseff.txt')
w7xhm_epseff = np.loadtxt('w7xhm_epseff.txt')

lee1_wout = netcdf_file('../DKESstuff/wout_Lee_1.nc', mode='r', mmap=False)
lee2_wout = netcdf_file('../DKESstuff/wout_Lee_2.nc', mode='r', mmap=False)
lee3_wout = netcdf_file('../DKESstuff/wout_Lee_3.nc', mode='r', mmap=False)
w7xhm_wout = netcdf_file('../initialConfigs/wout_W7X_REACTOR_FULL_woptim_forSfincs_correctPres.nc', mode='r', mmap=False)
lee1_ns = lee1_wout.variables['ns'][()]
lee2_ns = lee2_wout.variables['ns'][()]
lee3_ns = lee3_wout.variables['ns'][()]
w7xhm_ns = w7xhm_wout.variables['ns'][()]
if lee1_ns == lee2_ns == lee3_ns == w7xhm_ns:
    ns = lee1_ns
else:
    raise IOError('"ns" for the initial and optimized configurations differ - something is wrong.')
lee1_DMerc = lee1_wout.variables['DMerc'][()]
lee2_DMerc = lee2_wout.variables['DMerc'][()]
lee3_DMerc = lee3_wout.variables['DMerc'][()]
w7xhm_DMerc = w7xhm_wout.variables['DMerc'][()]

# Process data
_, fulls = createVMECGrids(ns)
sgrid = fulls[1:] # Ignore the magnetic axis because eps_eff cannot be evaluated there by NEO.
rhogrid = np.sqrt(sgrid)
data = np.column_stack((rhogrid, lee1_epseff, lee2_epseff, lee3_epseff, w7xhm_epseff, lee1_DMerc[1:], lee2_DMerc[1:], lee3_DMerc[1:], w7xhm_DMerc[1:])) # Need to ignore axis for magnetic well
filteredData = data[(data[:, 0] >= rhoMin) & (data[:, 0] <= rhoMax)]

# Plot data
plt.subplots(figsize=(xSizeInches, ySizeInches))
plt.plot(filteredData[:,0], filteredData[:,1])
plt.plot(filteredData[:,0], filteredData[:,2])
plt.plot(filteredData[:,0], filteredData[:,3])
plt.plot(filteredData[:,0], filteredData[:,4])
plt.yticks(np.arange(0,np.ceil(np.max(filteredData[:,1]))+1, 1))
plt.xlim(xmin=0, xmax=1)
plt.ylim(ymax=epseffMax)
plt.gca().xaxis.set_major_formatter(formatter)
plt.xlabel(r'$\rho$')
plt.ylabel(r'$\epsilon_\mathrm{eff}$ (%)')
plt.legend(['Configuration 1', 'Configuration 2', 'Configuration 3', 'W7-X High-Mirror'])
plt.savefig('epseff'+'.'+fileExt, bbox_inches='tight', dpi=dpi)

plt.subplots(figsize=(xSizeInches, ySizeInches))
plt.plot(filteredData[:,0], filteredData[:,5])
plt.plot(filteredData[:,0], filteredData[:,6])
plt.plot(filteredData[:,0], filteredData[:,7])
plt.plot(filteredData[:,0], filteredData[:,8])
plt.xlim(xmin=0, xmax=1)
plt.ylim(ymin=DMercYMin, ymax=DMercYMax)
plt.gca().xaxis.set_major_formatter(formatter)
plt.xlabel(r'$\rho$')
plt.ylabel(r'$D_\mathrm{Merc}$')
plt.legend(['Configuration 1', 'Configuration 2', 'Configuration 3', 'W7-X High-Mirror'])
plt.savefig('DMerc'+'.'+fileExt, bbox_inches='tight', dpi=dpi)

plt.show()
