# User options
file = 'wout_twentySecondObj4_cont13.nc' # Relative path to a VMEC wout file, or 'MICHAEL' for Michael Drevlak's robust electron root configuration
rnorm = 0.998
calcEpsEff = True
calcQIres = True # Note that this cannot be used with 'MICHAEL'
numQIres = 10
fileType = 'pdf' # Can be 'pdf' or 'png'
showPlots = True

#############################################################################################################################

# Module import
import sys
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

geomtools_lib_path = '/home/blee/src/geometrytools/python/lib'
sys.path.insert(0, geomtools_lib_path)

from geomlib import bcgeom, vmecgeom
from fluxcoorddiscr import fluxcoorddiscr
from eps_eff import eps_eff

if calcQIres:
    QI_lib_path = '/home/blee/src/QI_Target_Updates'
    sys.path.insert(0, QI_lib_path)
    from qi_functions import QuasiIsodynamicResidual as QIres # NOTE: you are using a custom version of this code
    from simsopt.mhd.vmec import Vmec
    import numpy as np

# Process configuration
if file == 'MICHAEL':
    file = '/home/blee/Downloads/hydra_Np04_20190108.bc'
    Geom = bcgeom(file, verbose=0)
else:
    Geom = bcgeom(vmecgeom(file), parallelize=True)

Booz = fluxcoorddiscr(Geom, rnorm=rnorm, Npol=101, Ntor=101, name='Boozer')

if calcEpsEff:
    epseff = eps_eff(Geom, Npol=199, Ntor=199, Nper=49)

if calcQIres:
    vmec = Vmec(file)
    vmec.keep_all_files = False
    svals = np.linspace(0, 1, num=numQIres)
    QIresVals = QIres(vmec, svals)

# Plot
fileSaveName = file.split('/')[-1]
if calcEpsEff:
    plt.plot(Geom.s, epseff)
    plt.xlabel('s')
    plt.ylabel('eps_eff')
    plt.savefig(fileSaveName+'_epseff.'+fileType, bbox_inches='tight', dpi=400)
    plt.figure()
if calcQIres:
    plt.plot(svals, QIresVals)
    plt.xlabel('s')
    plt.ylabel('QIres')
    plt.savefig(fileSaveName+'_QIres.'+fileType, bbox_inches='tight', dpi=400)
    plt.figure()
plt.plot(Geom.s, Geom.iota)
plt.xlabel('s')
plt.ylabel('iota')
plt.savefig(fileSaveName+'_iota.'+fileType, bbox_inches='tight', dpi=400)
fig1, ax1 = Booz.plot('B', cmap='jet')
fig1.savefig(fileSaveName+'_contours'+'_rnorm'+str(rnorm)+'.'+fileType, bbox_inches='tight', dpi=400)
fig2, ax2 = Booz.plot3d('B', title='Magnetic field, rnorm='+str(rnorm), torstride=2, polstride=2, cmap=None)
fig2.savefig(fileSaveName+'_3d'+'_rnorm'+str(rnorm)+'.'+fileType, bbox_inches='tight', dpi=400)
if showPlots:
    plt.show()
