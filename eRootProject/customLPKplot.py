fname = './DKESstuff/wout_Lee_2.nc'
figName = 'lee2_lpk'
xmin = 17 # either 17 or 35
xmax = 25 # either 25 or 43
ymin = -5.25
ymax = 5.25
numticks = 5
axisFontSize = 24
xSizeInches = xmax - xmin
ySizeInches = ymax - ymin
dpi = 600
fileExt = 'pdf'

import matplotlib as mpl
from matplotlib.ticker import LinearLocator
import matplotlib.pyplot as plt
import numpy as np
from simsopt.mhd.vmec import Vmec
import math
import sys, os

plt.rc('font', size=axisFontSize)

vmec = Vmec(fname)

nfp = vmec.wout.nfp
lasym = vmec.wout.lasym
xn = vmec.wout.xn
xm = vmec.wout.xm
xn_nyq = vmec.wout.xn_nyq
xm_nyq = vmec.wout.xm_nyq
rmnc = vmec.wout.rmnc.T
zmns = vmec.wout.zmns.T
raxis_cc = vmec.wout.raxis_cc
zaxis_cs = vmec.wout.zaxis_cs
nmodes = len(xn)
ns = vmec.wout.ns
if lasym == 1:
    rmns = vmec.wout.rmns
    zmnc = vmec.wout.zmnc
    raxis_cs = vmec.wout.raxis_cs
    zaxis_cc = vmec.wout.zaxis_cc
else:
    rmns = 0*rmnc
    zmnc = 0*rmnc
    raxis_cs = 0*raxis_cc
    zaxis_cc = 0*raxis_cc

fig = plt.figure(figsize=(xSizeInches, ySizeInches))
fig.patch.set_facecolor('white')

ntheta = 200
nzeta = 5
theta = np.linspace(0,2*np.pi,num=ntheta)
zeta = np.linspace(0,np.pi/nfp,num=nzeta,endpoint=True)
iradius = ns-1

R = np.zeros((ntheta,nzeta))
Z = np.zeros((ntheta,nzeta))
for itheta in range(ntheta):
    for izeta in range(nzeta):
        for imode in range(nmodes):
            angle = xm[imode]*theta[itheta] - xn[imode]*zeta[izeta]
            R[itheta,izeta] = R[itheta,izeta] + rmnc[iradius,imode]*math.cos(angle) + rmns[iradius,imode]*math.sin(angle)
            Z[itheta,izeta] = Z[itheta,izeta] + zmns[iradius,imode]*math.sin(angle) + zmnc[iradius,imode]*math.cos(angle)

for ind in range(nzeta):
    plt.plot(R[:,ind], Z[:,ind], '-')
    plt.plot(R[:,ind], Z[:,ind], '-')
    plt.plot(R[:,ind], Z[:,ind], '-')
    plt.plot(R[:,ind], Z[:,ind], '-')

plt.axis((xmin, xmax, ymin, ymax))
plt.gca().xaxis.set_major_locator(LinearLocator(numticks=numticks))
plt.xlabel('R (m)')
plt.ylabel('Z (m)')

plt.savefig(figName+'.'+fileExt, bbox_inches='tight', dpi=dpi)
plt.show()
plt.close()
