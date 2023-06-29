# This script plots the results from epseff.txt and magwell.txt, which were produced using NEO and STELLOPT for *ALL SURFACES POSSIBLE* in the configuration.

#Inputs
rhoMin = 0
rhoMax = 0.85
fileExt = 'png'

########################################################################################

# Import modules
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/blee/src/stelloptPlusSfincs/src')
from dataProc import createVMECGrids

# Load data
epseff = np.loadtxt('epseff.txt')
epseff_w7xhm = np.loadtxt('epseff_w7xhm.txt')
magwell = np.loadtxt('magwell.txt')
magwell_w7xhm = np.loadtxt('magwell_w7xhm.txt')

# Process data
ns = len(epseff) + 1 # The magnetic axis is ignored by NEO, so we must include the + 1
_, fulls = createVMECGrids(ns)
sgrid = fulls[1:] # Again, ignore the magnetic axis
rhogrid = np.sqrt(sgrid)
data = np.column_stack((rhogrid, epseff, epseff_w7xhm, magwell[1:], magwell_w7xhm[1:])) # Need to ignore axis for magnetic well
filteredData = data[(data[:, 0] >= rhoMin) & (data[:, 0] <= rhoMax)]

# Plot data
plt.figure()
plt.plot(filteredData[:,0], filteredData[:,1])
plt.plot(filteredData[:,0], filteredData[:,2])
plt.yticks(np.arange(0,np.ceil(np.max(filteredData[:,1]))+1, 1))
plt.xlabel(r'$\rho$')
plt.ylabel(r'$\epsilon_\mathrm{eff}$ (%)')
plt.legend(['Opt. Config.','Init. Config.'])
plt.savefig('epseff'+'.'+fileExt, bbox_inches='tight', dpi=400)

plt.figure()
plt.plot(filteredData[:,0], filteredData[:,3])
plt.plot(filteredData[:,0], filteredData[:,4])
plt.xlabel(r'$\rho$')
plt.ylabel(r'Magnetic well')
plt.legend(['Opt. Config.','Init. Config.'])

plt.savefig('magwell'+'.'+fileExt, bbox_inches='tight', dpi=400)

plt.show()
