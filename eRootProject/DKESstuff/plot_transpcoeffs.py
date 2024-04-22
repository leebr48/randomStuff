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
fileExt = 'pdf'

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
dkfileLee1='Lee_1b.dk'
bcfileLee1='Lee_1.bc'
dkfileLee2='Lee_2.dk'
bcfileLee2='Lee_2.bc'
dkfileLee3='Lee_3.dk'
bcfileLee3='Lee_3.bc'

dkhm1=DKESdata(dkfilehm1,bcfilehm1,reduced=False,minorradiusdefinition='W7AS')
dkDrevlak=DKESdata(dkfileDrevlak,bcfileDrevlak,reduced=False,minorradiusdefinition='W7AS')
dkLee1=DKESdata(dkfileLee1,bcfileLee1,reduced=False,minorradiusdefinition='W7AS')
dkLee2=DKESdata(dkfileLee2,bcfileLee2,reduced=False,minorradiusdefinition='W7AS')
dkLee3=DKESdata(dkfileLee3,bcfileLee3,reduced=False,minorradiusdefinition='W7AS')

rind=2 #Chosen radius index
rhm1=dkhm1.r[rind]
rDrevlak=dkDrevlak.r[rind]
rLee1=dkLee1.r[rind]
rLee2=dkLee2.r[rind]
rLee3=dkLee3.r[rind]

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

figLee1, axLee1 =dkLee1.plotD11star_vs_nustar(rLee1, errorbars=True, title='')
axLee1.set_xlim(left=lowerXLim, right=upperXLim)
axLee1.set_ylim(bottom=lowerYLim, top=upperYLim)
figLee1.set_size_inches(xSizeInches, ySizeInches)
figLee1.savefig('lee1_r' + str(rLee1) + '.' + fileExt, bbox_inches='tight', dpi=400)
figLee1.show()

figLee2, axLee2 =dkLee2.plotD11star_vs_nustar(rLee2, errorbars=True, title='')
axLee2.set_xlim(left=lowerXLim, right=upperXLim)
axLee2.set_ylim(bottom=lowerYLim, top=upperYLim)
figLee2.set_size_inches(xSizeInches, ySizeInches)
figLee2.savefig('lee2_r' + str(rLee2) + '.' + fileExt, bbox_inches='tight', dpi=400)
figLee2.show()

figLee3, axLee3 =dkLee3.plotD11star_vs_nustar(rLee3, errorbars=True, title='')
axLee3.set_xlim(left=lowerXLim, right=upperXLim)
axLee3.set_ylim(bottom=lowerYLim, top=upperYLim)
figLee3.set_size_inches(xSizeInches, ySizeInches)
figLee3.savefig('lee3_r' + str(rLee3) + '.' + fileExt, bbox_inches='tight', dpi=400)
figLee3.show()

wait_for_user()
