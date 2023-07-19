# This script plots epseff and DMerc for the initial and final configurations. These are calculated on *ALL SURFACES POSSIBLE* for each configuration.

#Inputs
rhoMin = 0
rhoMax = 0.85
DMercYMin = -0.001
DMercYMax = 0.01
fileExt = 'png'

########################################################################################

# Import modules
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.io import netcdf_file
sys.path.append('/home/blee/src/stelloptPlusSfincs/src')
from dataProc import createVMECGrids

# Load data
epseff_opt = np.loadtxt('epseff_opt.txt')
epseff_w7xhm = np.loadtxt('epseff_w7xhm.txt')

wout_opt = netcdf_file('../DKESstuff/wout_Lee_1.nc', mode='r', mmap=False)
wout_init = netcdf_file('../initialConfigs/wout_W7X_REACTOR_FULL_woptim_forSfincs_correctPres.nc', mode='r', mmap=False)
ns_opt = wout_opt.variables['ns'][()]
ns_init = wout_init.variables['ns'][()]
if ns_opt == ns_init:
    ns = ns_opt
else:
    raise IOError('"ns" for the initial and optimized configurations differ - something is wrong.')
DMerc_opt = wout_opt.variables['DMerc'][()]
DMerc_init = wout_init.variables['DMerc'][()]

# Process data
_, fulls = createVMECGrids(ns)
sgrid = fulls[1:] # Ignore the magnetic axis because eps_eff cannot be evaluated there by NEO.
rhogrid = np.sqrt(sgrid)
data = np.column_stack((rhogrid, epseff_opt, epseff_w7xhm, DMerc_opt[1:], DMerc_init[1:])) # Need to ignore axis for magnetic well
filteredData = data[(data[:, 0] >= rhoMin) & (data[:, 0] <= rhoMax)]

# Plot data
plt.figure()
plt.plot(filteredData[:,0], filteredData[:,1])
plt.plot(filteredData[:,0], filteredData[:,2])
plt.yticks(np.arange(0,np.ceil(np.max(filteredData[:,1]))+1, 1))
plt.xlim(xmin=0)
plt.xlabel(r'$\rho$')
plt.ylabel(r'$\epsilon_\mathrm{eff}$ (%)')
plt.legend(['Opt. Config.', 'W7-X High-Mirror'])
plt.savefig('epseff'+'.'+fileExt, bbox_inches='tight', dpi=400)

plt.figure()
plt.plot(filteredData[:,0], filteredData[:,3])
plt.plot(filteredData[:,0], filteredData[:,4])
plt.xlim(xmin=0)
plt.ylim(ymin=DMercYMin, ymax=DMercYMax)
plt.xlabel(r'$\rho$')
plt.ylabel(r'$D_\mathrm{Merc}$')
plt.legend(['Opt. Config.', 'Init. Config.'])

plt.savefig('DMerc'+'.'+fileExt, bbox_inches='tight', dpi=400)

plt.show()
