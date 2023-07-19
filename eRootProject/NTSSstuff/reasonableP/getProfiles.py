# Load modules
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('/home/blee/src/stelloptPlusSfincs/src')
from IO import makeStringForStellopt

# Important constant
eVToJ = 1.602176634e-19 # In STELLOPT, temperatures are written in eV but pressures are written in Pa

# Load data
data = np.loadtxt('profs_reasonableP', skiprows=1)

# Load minor radius
r = data[:,0]
a = np.max(r)

# Handy function
def newProf(ind, dtype):
    prof = data[:,ind]
    s = (r ** 2) / (a ** 2) # Profiles are arranged such that the first column is r (m), and the second column is the value
    if dtype == 't':
        profUse = prof * 1000 # use eV instead of keV
    elif dtype == 'n':
        profUse = prof * 1e20 # count by 1's, not 10^20's
    else:
        raise IOError('Invalid <dtype> choice.')
    return np.column_stack((s, profUse))

# Load profile information
ne = newProf(1, 'n')
nD = newProf(2, 'n')
nT = newProf(3, 'n')
nHe = newProf(4, 'n')

te = newProf(5, 't')
tD = newProf(6, 't')
tT = newProf(7, 't')
tHe = tT

# Output profile information
# Note that all the profiles are evaluated on the same r grid, so we don't need to worry about fitting splines or anything
out = '&BEAMS3D_INPUT' + '\n'
out += makeStringForStellopt('NE_AUX_S', ne[:,0])
out += makeStringForStellopt('NE_AUX_F', ne[:,1])
out += makeStringForStellopt('TE_AUX_S', te[:,0])
out += makeStringForStellopt('TE_AUX_F', te[:,1])
out += makeStringForStellopt('NI_AUX_S', nD[:,0])
out += makeStringForStellopt('NI_AUX_F(1,:)', nD[:,1])
out += makeStringForStellopt('NI_AUX_F(2,:)', nT[:,1])
out += makeStringForStellopt('NI_AUX_F(3,:)', nHe[:,1])
out += makeStringForStellopt('TI_AUX_S', tD[:,0])
out += makeStringForStellopt('TI_AUX_F', tD[:,1])
out += makeStringForStellopt('NI_AUX_M', [3.343583746e-27, 5.008267660e-27, 6.646479070e-27])
out += makeStringForStellopt('NI_AUX_Z', [1, 1, 2])
out += '/\n'

with open('input.reactorProfiles', 'w') as f:
    f.write(out)

# Create one file with profiles that's easily Numpy-readable
out = np.column_stack((ne[:,0], ne[:,1], nD[:,1], nT[:,1], nHe[:,1], te[:,1], tD[:,1], tD[:,1], tD[:,1]))
hdr = 'VMEC S, ne (m^-3), nD (m^-3), nT (m^-3), nHe (m^-3), Te (eV), TD (eV), TT (eV), THe (eV)'
np.savetxt('reactorProfiles.dat', out, header=hdr)

# Calculate and output pressure profile
p = eVToJ * (ne[:,1]*te[:,1] + nD[:,1]*tD[:,1] + nT[:,1]*tT[:,1] + nHe[:,1]*tHe[:,1]) # Pascals

presString = makeStringForStellopt('GAMMA', [0])
presString += makeStringForStellopt('PMASS_TYPE', "'akima_spline'")
presString += makeStringForStellopt('PRES_SCALE', [1])
presString += makeStringForStellopt('AM_AUX_S', ne[:,0])
presString += makeStringForStellopt('AM_AUX_F', p)

with open('pressureProfile.txt', 'w') as f:
    f.write(presString)

# Make some plots
plt.plot(ne[:,0], ne[:,1] / 1e20, label='e')
plt.plot(nD[:,0], nD[:,1] / 1e20, label='D')
plt.plot(nT[:,0], nT[:,1] / 1e20, label='T')
plt.plot(nHe[:,0], nHe[:,1] / 1e20, label='He')
plt.legend(loc='best')
plt.xlabel('Normalized toroidal flux $s$')
plt.ylabel('Density ($\mathrm{10^{20}/m^{3}}$)')
plt.savefig('n.pdf', bbox_inches='tight', dpi=400)

plt.figure()
plt.plot(te[:,0], te[:,1], label='e')
plt.plot(tD[:,0], tD[:,1], label='ions')
plt.legend(loc='best')
plt.xlabel('Normalized toroidal flux $s$')
plt.ylabel('Temperature (keV)')
plt.savefig('T.pdf', bbox_inches='tight', dpi=400)

plt.figure()
plt.plot(ne[:,0], p)
plt.xlabel('Normalized toroidal flux $s$')
plt.ylabel('Pressure (Pa)')
plt.savefig('p.pdf', bbox_inches='tight', dpi=400)

plt.figure()
plt.plot(np.sqrt(ne[:,0]), p)
plt.xlabel(r'Normalized radius $\rho$')
plt.ylabel('Pressure (Pa)')
plt.savefig('p_rho.pdf', bbox_inches='tight', dpi=400)
