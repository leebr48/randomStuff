import numpy as np
from scipy.io import netcdf_file
import matplotlib.pyplot as plt
import sys
sys.path.insert(1,'/home/blee/src/stelloptPlusSfincs/src')
from IO import makeStringForStellopt

# Fudge factor
neScale = 1.0243902439127575

# Inputs from Michael's input file
AM =  [8.2e+05, 0.0, -8.2e+05, -8.2e+05, 4.1e+05, 4.1e+05]
pres_scale = 1.004
ns = 64

# Calculate your approximation of the pressure vector
eVToJ = 1.602176634e-19 # In STELLOPT, temperatures are written in eV but pressures are written in Pa

def polynomial(coeffs):
    return lambda x: sum(a*x**i for i, a in enumerate(coeffs))

presPoly = polynomial(AM)

sVals = np.linspace(0, 1, num=21)
rawPVals = presPoly(sVals)
scaledPVals = pres_scale * np.array(rawPVals) # Pascals

# Get his version of the pressure vector
fullgrid = np.linspace(0, 1, num=ns)
diff = (fullgrid[1] - fullgrid[0]) / 2 
halfgrid = fullgrid - diff
halfgrid[0] = 0

f = netcdf_file('wout_DrevlakERoot.nc', mode='r', mmap=False)
presf = f.variables['presf'][()] # Pascals
pres = f.variables['pres'][()] # Pascals

# Plot to compare
plt.plot(fullgrid, presf, label='presf')
#plt.plot(halfgrid, pres, label='pres')
plt.plot(sVals, scaledPVals, label='new')
plt.legend(loc='best')

# Get NE and TE curves. We can assume there is one ion species for now and use addIons.py to fix that later.

# Assume the TE curve is flat:
teFunc = lambda s: 19000 * (1-s) + 1000 # eV

te = teFunc(sVals) # eV
ti = te

ne = scaledPVals / (2 * eVToJ * te)
ne = ne * neScale # NOTE: adjustment so that addIons.py will give use the right pressure
ni = ne

plt.figure()
plt.plot(sVals, ne)
#plt.show()

# Print out the vectors you have for now
print(makeStringForStellopt('AM_AUX_S', sVals).strip())
print(makeStringForStellopt('AM_AUX_F', scaledPVals).strip())
print('****************************')
print(makeStringForStellopt('NE_AUX_S', sVals).strip())
print(makeStringForStellopt('NE_AUX_F', ne).strip())
print(makeStringForStellopt('TE_AUX_S', sVals).strip())
print(makeStringForStellopt('TE_AUX_F', te).strip())
print(makeStringForStellopt('NI_AUX_S', sVals).strip())
print(makeStringForStellopt('NI_AUX_F', ni).strip())
print(makeStringForStellopt('TI_AUX_S', sVals).strip())
print(makeStringForStellopt('TI_AUX_F', ti).strip())
print(makeStringForStellopt('NI_AUX_M', [3.343583746e-27]).strip())
print(makeStringForStellopt('NI_AUX_Z', [1]).strip())
print('Note that changing the number of ions with addIons.py will change P. You will need to change the scale factor in this script to account for that. (Just divide the first element of the original P vector by the first element of the addIons P vector, and you will get the fudge factor.)')
