#!/usr/bin/env python
import sys, os
import numpy as np
lib_path=os.getenv('NEOTRANSP_PYTHON_LIB')
#  Then you can put a copy of this file wherever you want. 
sys.path.insert(0,lib_path)
from neolib import DKESdata, profile_data, transp_data, wait_for_user

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% LOAD DATA FILES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dkfilehm1='./DKESstuff/w7x-hm1.dk' # hm='high mirror'
bcfilehm1='./DKESstuff/w7x-hm1.bc'
dkfileDrevlak='./DKESstuff/hydra_Np04_20190108.dk'
bcfileDrevlak='DKESstuff/hydra_Np04_20190108.bc'
#dkfileLee='./DKESstuff/Lee_1b.dk'
#bcfileLee='./DKESstuff/Lee_1.bc'
dkfileLee='./DKESstuff/Lee_2.dk'
bcfileLee='./DKESstuff/Lee_2.bc'

dkhm1=DKESdata(dkfilehm1,bcfilehm1,reduced=False,minorradiusdefinition='W7AS')
dkDrevlak=DKESdata(dkfileDrevlak,bcfileDrevlak,reduced=False,minorradiusdefinition='W7AS')
dkLee=DKESdata(dkfileLee,bcfileLee,reduced=False,minorradiusdefinition='W7AS')

rind=3 #Chosen radius index
rhm1=dkhm1.r[rind]
rDrevlak=dkDrevlak.r[rind]
rLee=dkLee.r[rind]

fighm1, axhm1 =dkhm1.plotD11star_vs_nustar(rhm1,errorbars=True)
axhm1.set_xlim(left=2e-5, right=4e2)
axhm1.set_ylim(bottom=1e-3, top=2e2)
fighm1.show()

figDrevlak, axDrevlak =dkDrevlak.plotD11star_vs_nustar(rDrevlak,errorbars=True)
axDrevlak.set_xlim(left=2e-5, right=4e2)
axDrevlak.set_ylim(bottom=1e-3, top=2e2)
figDrevlak.show()

figLee, axLee =dkLee.plotD11star_vs_nustar(rLee,errorbars=True)
axLee.set_xlim(left=2e-5, right=4e2)
axLee.set_ylim(bottom=1e-3, top=2e2)
figLee.show()

wait_for_user()
