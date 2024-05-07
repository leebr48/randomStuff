# This script was written by Alan Goodman.

import matplotlib.pyplot as plt
import numpy as np
from simsopt.mhd.vmec import Vmec
import math
import sys, os

fname = './DKESstuff/wout_Lee_1.nc'
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

fig = plt.figure("Poincare Plots",figsize=(14,7))
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

numCols = 3
numRows = 2
plotNum = 1

plt.subplot(numRows,numCols,plotNum)
#plt.subplot(1,1,1)
plotNum += 1
for ind in range(nzeta):
    plt.plot(R[:,ind], Z[:,ind], '-')
    plt.plot(R[:,ind], Z[:,ind], '-')
    plt.plot(R[:,ind], Z[:,ind], '-')
    plt.plot(R[:,ind], Z[:,ind], '-')
plt.gca().set_aspect('equal',adjustable='box')
#plt.legend(fontsize='x-small')
plt.xlabel('R')
plt.ylabel('Z')
    
ntheta = 500
nzeta = 5
nradius = 10
radii = np.linspace(1,ns-1,nradius)
radii = np.floor(radii)
theta = np.linspace(0,2*np.pi,num=ntheta)
zeta = np.linspace(0,2*np.pi/nfp/2,num=nzeta,endpoint=True)

from fractions import Fraction
denoms = np.linspace(0, 1, num=nzeta)
titles = []
for idenom in range(len(denoms)):
    denom = denoms[idenom]
    anglestr = str(Fraction(denom).limit_denominator().numerator) + r'$\pi/$' +  str(Fraction(denom).limit_denominator().denominator)
    if anglestr[:7] == '1$\pi/$':
        anglestr = '$\pi/$' + anglestr[7:]
    titles.append( r'$\phi=$' + anglestr )
titles[0] = r'$\phi=0$'
titles[-1] = r'$\phi=\pi$'

def FindBoundary(theta,phi,iradius):
    angle = xm*theta + xn*phi
    rb = np.sum(rmnc[iradius,:] * np.cos(angle))
    zb = np.sum(zmns[iradius,:] * np.sin(angle))
    if lasym == True:
        rb += np.sum(rmns[iradius,:] * np.sin(angle))
        zb += np.sum(zmnc[iradius,:] * np.cos(angle))
    return rb,zb

iradii = np.linspace(0,ns-1,num=nradius).round()
iradii = [int(i) for i in iradii]
R = np.zeros((ntheta,nzeta,nradius))
Z = np.zeros((ntheta,nzeta,nradius))
for itheta in range(ntheta):
    for izeta in range(nzeta):
        for iradius in range(len(radii)):
            rad = int(radii[iradius])
            R[itheta,izeta,iradius], Z[itheta,izeta,iradius] = FindBoundary(theta[itheta],zeta[izeta],rad)

Raxis = np.zeros(nzeta)
Zaxis = np.zeros(nzeta)
for jn in range(len(raxis_cc)):
    n = jn * nfp
    sinangle = np.sin(n * zeta)
    cosangle = np.cos(n * zeta)

    Raxis += raxis_cc[jn] * cosangle
    Zaxis += zaxis_cs[jn] * sinangle

for izeta in range(nzeta):
    plt.subplot(numRows,numCols,plotNum)
    plotNum += 1
    for iradius in range(nradius):
        plt.plot(R[:,izeta,iradius], Z[:,izeta,iradius], '-')
    plt.plot(Raxis[izeta],Zaxis[izeta],'xr')
    plt.gca().set_aspect('equal',adjustable='box')
    plt.xlabel('R')
    plt.ylabel('Z')
    plt.title((titles[izeta]))

plt.subplots_adjust(wspace=0.39, hspace=0.444)
plt.figtext(0.5,0.99,os.path.abspath(fname),ha='center',va='top',fontsize=6)

plt.show()
plt.close()
