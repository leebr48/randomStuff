dataFile = 'stellopt.W7X_REACTOR_woptim_forSfincs'
fname = 'lee1_epseff'
ftype = 'txt'

import numpy as np

data = np.loadtxt(dataFile, skiprows=7, max_rows=127, usecols=(2))

epseff = data ** (2/3) * 100

np.savetxt(fname + '.' + ftype, epseff)
