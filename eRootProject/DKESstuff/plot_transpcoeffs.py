#!/usr/bin/env python

# This script was originally written by Hakan Smith.

# Inputs
lowerXLim = 1.0e-4
upperXLim = 1.0e0
lowerYLim = 1.5e-3
upperYLim = 5.0e1
axisFontSize = 24
legendFontSize = 14
xSizeInches = 7.9
ySizeInches = 7.9
fileExt = 'png'

# Code
import sys, os
import numpy as np
lib_path=os.getenv('NEOTRANSP_PYTHON_LIB')
#  Then you can put a copy of this file wherever you want. 
sys.path.insert(0,lib_path)
from neolib import DKESdata, profile_data, transp_data, wait_for_user
import matplotlib.pyplot as plt

plt.rc('font', size=axisFontSize)
plt.rc('legend', fontsize=legendFontSize)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% LOAD DATA FILES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dkfilehm1='w7x-hm1.dk' # hm='high mirror'
bcfilehm1='w7x-hm1.bc'
dkfileDrevlak='hydra_Np04_20190108.dk'
bcfileDrevlak='hydra_Np04_20190108.bc'
dkfileLee='Lee_1b.dk'
bcfileLee='Lee_1.bc'

dkhm1=DKESdata(dkfilehm1,bcfilehm1,reduced=False,minorradiusdefinition='W7AS')
dkDrevlak=DKESdata(dkfileDrevlak,bcfileDrevlak,reduced=False,minorradiusdefinition='W7AS')
dkLee=DKESdata(dkfileLee,bcfileLee,reduced=False,minorradiusdefinition='W7AS')

rind=2 #Chosen radius index
rhm1=dkhm1.r[rind]
rDrevlak=dkDrevlak.r[rind]
rLee=dkLee.r[rind]

fighm1, axhm1 =dkhm1.plotD11star_vs_nustar(rhm1, errorbars=True, title='')
axhm1.set_xlim(left=lowerXLim, right=upperXLim)
axhm1.set_ylim(bottom=lowerYLim, top=upperYLim)
fighm1.set_size_inches(xSizeInches, ySizeInches)
fighm1.savefig('w7xhm_r' + str(rhm1) + '.' + fileExt, bbox_inches='tight', dpi=400)
fighm1.show()

figDrevlak, axDrevlak =dkDrevlak.plotD11star_vs_nustar(rDrevlak, errorbars=True, title='')
axDrevlak.set_xlim(left=lowerXLim, right=upperXLim)
axDrevlak.set_ylim(bottom=lowerYLim, top=upperYLim)
figDrevlak.set_size_inches(xSizeInches, ySizeInches)
figDrevlak.savefig('drevlak_r' + str(rDrevlak) + '.' + fileExt, bbox_inches='tight', dpi=400)
figDrevlak.show()

figLee, axLee =dkLee.plotD11star_vs_nustar(rLee, errorbars=True, title='')
axLee.set_xlim(left=lowerXLim, right=upperXLim)
axLee.set_ylim(bottom=lowerYLim, top=upperYLim)
figLee.set_size_inches(xSizeInches, ySizeInches)
figLee.savefig('lee_r' + str(rLee) + '.' + fileExt, bbox_inches='tight', dpi=400)
figLee.show()

wait_for_user()
