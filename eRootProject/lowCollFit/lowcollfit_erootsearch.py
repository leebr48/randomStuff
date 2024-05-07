#!/usr/bin/env python

# This script was written by Hakan Smith.

import sys, os
import numpy as np
lib_path=os.getenv('NEOTRANSP_PYTHON_LIB')
sys.path.insert(0,lib_path)
#sys.path.insert(0, os.getenv('GEOMETRYTOOLS_HOME')+'/python/lib')
from neolib import DKESdata, profile_data, physcnst, transp_data, wait_for_user
import matplotlib.pyplot as plt
#from W7Xdatabase import neoclass_database

#Reference
#equil='w7x-hm1-pb4.4'
equil='hydra_Np04_20190108'

if equil=='hydra_Np04_20190108':
  B00=5.1 #T
  #scaling=0.9138 (Craig's value)
  scaling=1 #I simplify here
  #Input profiles
  TekeV=15.1
  TikeV=14.0
  ne20=1.64
  dlnTedr=-0.81
  dlnTidr=-0.83
  dlnndr=-0.47

elif equil=='w7x-hm1-pb4.4':
  B00=5.4 #T
  scaling=3.6477
  #Input profiles
  TekeV=13.6
  TikeV=12.6
  ne20=1.95
  dlnTedr=-0.54
  dlnTidr=-0.53
  dlnndr=-0.31
  

#basepath='/afs/ipp-garching.mpg.de/home/s/smithh/Forskning/CraigBeidler/Reactor_calc2023/'
dk=DKESdata(equil+'.dk.nc',minorradiusdefinition='W7AS')

dk.set_scaling(scaling)
#print(scaling)
#print(dk.get_scaling())
#print(dk.radius_scaling) #VERY STRANGE, RETURNS 1
#print('----------')


physconst=physcnst()

#constants
B=3
dkrind=3 #3:halfradius
#r=dk.r_interp(rho,'r',interpvar='rho')
#r=dk.r[dkrind]*dk.get_scaling()
r=dk.get_r()[dkrind]
rho=r/dk.minorradius()
#r=dk.minorradius()*rho
R=dk.r_interp(r,'R')
absiota=np.abs(dk.r_interp(r,'iota'))
drdrHM=1/np.abs(dk.GeomLH.psi_a*2/dk.r_interp(r,'GeomLH.B00')/dk.GeomLH.minorradiusW7AS**2)
PSfactor=np.abs(dk.r_interp(r,'PSfactor'))*drdrHM**2

#Temporary corr:
#PSfactor/=scaling**2 #might have to chang here if I change in neolib
#print('rho='+str(rho))
#print('dk.r[dkrind]='+str(dk.r[dkrind]))
#print('r='+str(r))
#print('R='+str(R))
#print('PSfactor='+str(PSfactor))
#sys.exit()



prof=profile_data(['e','H'],r=r)
prof.set_T_keV('e',TekeV)
prof.set_T_keV('H',TikeV)
prof.set_dlnTdr('e',dlnTedr)
prof.set_dlnTdr('H',dlnTidr)
prof.set_n_20m3('e',ne20)
prof.set_dlnndr('e',dlnndr)
prof.make_quasineutral('H')
vTe=prof.vT('e')
vTi=prof.vT('H')

qe=-physconst.e
qi=physconst.e
C_e=np.pi/16*physconst.me**2*vTe**3/(absiota*R*B**2*qe**2)*4/np.sqrt(np.pi)
C_i=np.pi/16*physconst.mp**2*vTi**3/(absiota*R*B**2*qi**2)*4/np.sqrt(np.pi)

xmin=1e-2
xmax=1e3
Nx=50
x=np.logspace(np.log10(xmin),np.log10(xmax),Nx)
nuoverv_e=prof.nuoverv('e',x)[0]
nuoverv_i=prof.nuoverv('H',x)[0]
nustar_e=nuoverv_e*(R/absiota)
nustar_i=nuoverv_i*(R/absiota)
ve=vTe*x
vi=vTi*x

#D11star at Er=0
(g110e,dum1,dum2)=dk.transpcoeff(r,nuoverv_e,0.0,Er_isowndimension=False)
(g110i,dum1,dum2)=dk.transpcoeff(r,nuoverv_i,0.0,Er_isowndimension=False)
D11star0e=np.abs(8.0/np.pi*R*absiota*(-g110e)*drdrHM**2)
D11star0i=np.abs(8.0/np.pi*R*absiota*(-g110i)*drdrHM**2)

EkVovBmin=0.0 #kV must be 0
EkVovBmax=10.5   #kV
NE=148
#EkVovB=np.insert(np.logspace(np.log10(EkVovBmin),np.log10(EkVovBmax),NE-1),0,0.0)
#EkVovB=np.insert(np.linspace(EkVovBmin,EkVovBmax,NE-1),0,0.0)
EkVovB=np.linspace(EkVovBmin,EkVovBmax,NE)
EovB=EkVovB*1e3

Estar_e=EovB[:,None]/ve[None,:]
Estar_i=EovB[:,None]/vi[None,:]
y_e=Estar_e/nustar_e[None,:]
y_i=Estar_i/nustar_i[None,:]
qEovT_e=-1*EkVovB/TekeV*B
qEovT_i= 1*EkVovB/TikeV*B
A1e=dlnndr-qEovT_e-1.5*dlnTedr
A1i=dlnndr-qEovT_i-1.5*dlnTidr
A1e0=dlnndr-1.5*dlnTedr
A1i0=dlnndr-1.5*dlnTidr
A1eE=-qEovT_e
A1iE=-qEovT_i
A2e=dlnTedr
A2i=dlnTidr

#print(qEovT_e)

#parameters
#g11_0=2.2e-5
#Npl=8*R*absiota/np.pi*g11_0*300
#nustarPS=8*absiota**2/np.pi*g11_0*30000
#nustar0=8*R**2/np.pi*g11_0
EovervBres=absiota*r/R
print('EovervBres='+str(EovervBres))
yr_e=Estar_e/EovervBres
yr_i=Estar_i/EovervBres

EovervB=[0.0,1e-5,3e-5,1e-4,3e-4,1e-3,3e-3]

#Example comparison with DKES data

Transp=transp_data(prof,dk,roots='i&e',B00=B,Ermin_kVm=-19,Ermax_kVm=28,momcorr=False,
                   prespline=False,parallelize=False)#,fignr=7)
ErkVmoverB_eroot=Transp.get_ErkVm(root='e')[0]/B
ErkVmoverB_iroot=Transp.get_ErkVm(root='i')[0]/B
print('eroot Er='+str(ErkVmoverB_eroot))
print('iroot Er='+str(ErkVmoverB_iroot))
ErkVmoverB_myroot=ErkVmoverB_eroot
if np.isnan(ErkVmoverB_eroot):
  ErkVmoverB_myroot=ErkVmoverB_iroot
figdk11p, axdk11p =dk.plotD11star_vs_nustar(r,prof=prof,spec=['e','H'],ErkVmoverB=ErkVmoverB_myroot,
                                            EovervB=EovervB,showGaussLaguerreterms=True,showequivtok=False)
#figdk11p.show()

L=Transp.transport_matrix_allradii(EkVovB[None,:]*B)[:,0,:,:,:]
#print(L.shape)
Ge0=-(L[0,:,0,0]*A1e0+L[0,:,0,1]*A2e)
Gi0=-(L[1,:,0,0]*A1i0+L[1,:,0,1]*A2i)
GeE=-(L[0,:,0,0]*A1eE)
GiE=-(L[1,:,0,0]*A1iE)

figA,axA=plt.subplots(1,1,sharex=True)
#axA[0].plot(EkVovB,A1e,'r',
#            EkVovB,A1i,'b')
#axA[0].set_ylabel(r'$A_1$')
axA.plot(EkVovB,Ge0,'r:',EkVovB,GeE,'r--',EkVovB,GeE+Ge0,'r',
            -EkVovB,Ge0,'r:',-EkVovB,-GeE,'r--',-EkVovB,-GeE+Ge0,'r',
            EkVovB,Gi0,'b:',EkVovB,GiE,'b--',EkVovB,GiE+Gi0,'b',
            -EkVovB,Gi0,'b:',-EkVovB,-GiE,'b--',-EkVovB,-GiE+Gi0,'b',
         EkVovB,GiE+Gi0-(GeE+Ge0),'k',-EkVovB,-GiE+Gi0-(-GeE+Ge0),'k')
axA.set_ylabel(r'$\Gamma$')
#figA.show()  
#wait_for_user()

if False:#True:#False:#True:
  profEscan=prof.repeat(len(EkVovB))
  TranspEscan=transp_data(profEscan,dk,roots='givenErkVm',givenErkVm=EkVovB*B,
                     B00=B,Ermin_kVm=-19,Ermax_kVm=28,momcorr=False,
                     prespline=False,parallelize=False,fignr=7)
  figtrE,axtrE=TranspEscan.plot(xlabel='Er')
  figtrE.show()
 
  wait_for_user()

dk_alpha=1.5
dk_y0=dk.fit.g11_sq[dkrind]/dk.fit.g11_0[dkrind]*(absiota/R)**dk_alpha
dk_gamma=dk.fit.ex_er[dkrind]
dk_nustar0=8*(R/scaling)**2/np.pi*dk.fit.g11_0[dkrind]#/scaling**2
dk_nustarPS=1/(16*absiota**2/3/np.pi*PSfactor)
#g11_0=2.2e-5
#Npl=8*R*absiota/np.pi*dk.fit.g11_0[dkrind]*300
#nustarPS=8*absiota**2/np.pi*g11_0*30000
#nustar0=8*R**2/np.pi*g11_0
#print('Npl='+str(Npl))
print('dk_nustarPS='+str(dk_nustarPS))
#print('dk_nustarPS='+str(1/(16*absiota**2/3/np.pi*dk.PSfactor[dkrind])))

#if equil=='w7x-hm1'
#  my_gamma=1.8
#  my_alpha=1.6
#  my_y0=dk_y0*2*1.5
if equil=='hydra_Np04_20190108':
  my_gamma=1.0
  my_alpha=1.5
  my_y0=0.18 
  #my_y0=0.05 
  my_nustar0=5.4e-4#/scaling**2
  my_nustarPS=dk_nustarPS
  Npl=0.09
  my_yr0=0.1
  my_yfexp=2
  my_PSplexp=1
if equil=='w7x-hm1-pb4.4':
  my_gamma=1.2
  my_alpha=1.5
  my_y0=0.33
  my_nustar0=dk_nustar0*1.05
  my_nustarPS=dk_nustarPS
  my_yr0=0.1
  my_yfexp=2
  if False:
    Npl=0.04
    my_PSplexp=0.5
  else:
    Npl=0.09
    my_PSplexp=1
print('my_y0='+str(my_y0))
print('my_nustar0='+str(my_nustar0))
#print('='+str())

y0str='_y0'+str(my_y0)

figdk11, axdk11 =dk.plotD11star_vs_nustar(r,EovervB=EovervB,showequivtok=False)
xlim=axdk11.get_xlim()
ylim=axdk11.get_ylim()

def Hfun(y,yr,y0,yr0,yfexp,gamma,alpha):
  y0fun=y0*(1+(yr/yr0)**yfexp)
  return (1+(y/y0fun)**(alpha*gamma))**(-1/gamma)
  
def Nfun1(nustar,nustar0,Npl,nustarPS,PSplexp=1):
  return nustar0/nustar+Npl*(1+(nustar/nustarPS/Npl)**PSplexp)**(1/PSplexp)

def Nfun2(nustar,nustar0,Npl,nustarPS,PSplexp=1):
  nuoverv=nustar/R*absiota
  (g110,dum1,dum2)=dk.transpcoeff(r,nuoverv,0.0,Er_isowndimension=False)
  D11starDKES=np.abs(8.0/np.pi*R*absiota*(-g110)*drdrHM**2)
  goodi=np.argmin((nustar-1e-7)**2)
  #print('D11starDKES[goodi]*nustar[goodi]='+str(D11starDKES[goodi]*nustar[goodi]))
  #print('nustar0='+str(nustar0))
  nustar0DKES=D11starDKES[goodi]*nustar[goodi]
  return D11starDKES+(nustar0-nustar0DKES)/nustar

#Nfun=Nfun1
Nfun=Nfun2

colors=['b','g','y','c','r','k']*5
for Eii in range(len(EovervB)):
  cmpEstar_e=EovervB[Eii]
  cmpy_e=cmpEstar_e/nustar_e
  cmpyr_e=cmpEstar_e/EovervBres
  #cmpy0fun_e=my_y0*(1+(cmpyr_e/my_yr0)**my_yfexp)
  #cmpHe=(1+(cmpy_e/cmpy0fun_e)**(my_alpha*my_gamma))**(-1/my_gamma)
  cmpHe=Hfun(cmpy_e,cmpyr_e,my_y0,my_yr0,my_yfexp,my_gamma,my_alpha)
  #cmpNe=my_nustar0/nustar_e+Npl*(1+(nustar_e/my_nustarPS/Npl)**my_PSplexp)**(1/my_PSplexp)
  cmpNe=Nfun(nustar_e,my_nustar0,Npl,my_nustarPS,PSplexp=my_PSplexp)
  cmpD11star_e=cmpNe*cmpHe
  axdk11.loglog(nustar_e,cmpD11star_e,colors[Eii]+'-.')

axdk11.set_xlim(xlim)
axdk11.set_ylim(ylim)
figdk11.show()
figdk11.savefig('lcf_'+equil+'_D11stvsnust'+y0str+'.eps')

#figdk113d,(axdk113d0,axdk113d1,axdk113d2) = dk.plot(dkrind,mode='lin')#,prof=prof,spec='H',ErkVmoverB=ErkVmoverB_myroot)
#figdk113d.show()


#wait_for_user()

#Calculation for the dk case
He=Hfun(y_e,yr_e,my_y0,my_yr0,my_yfexp,my_gamma,my_alpha)
Hi=Hfun(y_i,yr_i,my_y0,my_yr0,my_yfexp,my_gamma,my_alpha)
#y0fun_e=my_y0*(1+(yr_e/my_yr0)**my_yfexp)
#y0fun_i=my_y0*(1+(yr_i/my_yr0)**my_yfexp)
#He=(1+(y_e/y0fun_e)**(my_alpha*my_gamma))**(-1/my_gamma)
#Hi=(1+(y_i/y0fun_i)**(my_alpha*my_gamma))**(-1/my_gamma)
Ne=Nfun(nustar_e,my_nustar0,Npl,my_nustarPS,PSplexp=my_PSplexp)
Ni=Nfun(nustar_i,my_nustar0,Npl,my_nustarPS,PSplexp=my_PSplexp)
#Ne=dk_nustar0/nustar_e+Npl*(1+(nustar_e/my_nustarPS/Npl)**my_PSplexp)**(1/my_PSplexp)
#Ni=dk_nustar0/nustar_i+Npl*(1+(nustar_i/my_nustarPS/Npl)**my_PSplexp)**(1/my_PSplexp)
D11star_e=Ne[None,:]*He
D11star_i=Ni[None,:]*Hi

L11e=C_e*np.trapz((x**5*np.exp(-x**2))[None,:]*D11star_e,x=x,axis=1)
L12e=C_e*np.trapz((x**7*np.exp(-x**2))[None,:]*D11star_e,x=x,axis=1)
L11i=C_i*np.trapz((x**5*np.exp(-x**2))[None,:]*D11star_i,x=x,axis=1)
L12i=C_i*np.trapz((x**7*np.exp(-x**2))[None,:]*D11star_i,x=x,axis=1)

Govn_e0=-L11e*A1e0-L12e*A2e
Govn_i0=-L11i*A1i0-L12i*A2i
Govn_eE=-L11e*A1eE
Govn_iE=-L11i*A1iE
Govn_ep=Govn_e0+Govn_eE
Govn_ip=Govn_i0+Govn_iE
Govn_em=Govn_e0-Govn_eE
Govn_im=Govn_i0-Govn_iE
DGonv_p=Govn_ip-Govn_ep
DGonv_m=Govn_im-Govn_em

def calcMaxw(EkVovB,DGonv_p,DGonv_m):
  Ipos=np.nan
  Ineg=np.nan
  tmp=np.where(np.diff(np.sign(DGonv_p))==2)[0]
  erootexists=(np.size(tmp)>0)
  if erootexists:
    i0=tmp[0]
    Ipos=np.trapz(DGonv_p[:i0],x=EkVovB[:i0])
    Ipos+=-DGonv_p[i0]**2/2*(EkVovB[i0+1]-EkVovB[i0])/(DGonv_p[i0+1]-DGonv_p[i0])
  tmp=np.where(np.diff(np.sign(DGonv_m))==-2)[0]
  irootexists=(np.size(tmp)>0)
  if irootexists:
    i0=tmp[0]
    Ineg=np.trapz(DGonv_m[:i0],x=EkVovB[:i0])
    Ineg+=-DGonv_m[i0]**2/2*(EkVovB[i0+1]-EkVovB[i0])/(DGonv_m[i0+1]-DGonv_m[i0])
  MaxwI=Ipos+Ineg

  tmp=np.where(np.diff(np.sign(DGonv_p))==-2)[0]
  urootexists=(np.size(tmp)>0)
  if urootexists:
    i0=tmp[0]
    Iposu=np.trapz(DGonv_p[:i0],x=EkVovB[:i0])
    Iposu+=-DGonv_p[i0]**2/2*(EkVovB[i0+1]-EkVovB[i0])/(DGonv_p[i0+1]-DGonv_p[i0])
  MaxwQ=np.nan
  if erootexists and urootexists and irootexists:
    MaxwQ=(Iposu+Ineg)/(Iposu+Ineg-(Ipos-Iposu))
  return MaxwQ#,MaxwI

MaxwQ=calcMaxw(EkVovB,DGonv_p,DGonv_m)
print('MaxwQ='+str(MaxwQ))
  
axA.plot( EkVovB,Govn_e0,'m:', EkVovB,Govn_eE,'m--',  EkVovB,Govn_ep,'m',
         -EkVovB,Govn_e0,'m:',-EkVovB,-Govn_eE,'m--',-EkVovB,Govn_em,'m',
          EkVovB,Govn_i0,'c:', EkVovB,Govn_iE,'c--',  EkVovB,Govn_ip,'c',
          -EkVovB,Govn_i0,'c:',-EkVovB,-Govn_iE,'c--',-EkVovB,Govn_im,'c',
         EkVovB,Govn_iE+Govn_i0-(Govn_eE+Govn_e0),'y',-EkVovB,-Govn_iE+Govn_i0-(-Govn_eE+Govn_e0),'y')
axA.grid(True)
ylim=axA.get_ylim()
EkvovB_res_i=EovervBres*vTi/1e3
EkvovB_0_i=EkvovB_res_i*my_yr0
axA.plot([EkvovB_res_i*0.30]*2,ylim,'b:')
#figA.show()

#simpler plot
figA2,axA2=plt.subplots(1,1)
#axA[0].plot(EkVovB,A1e,'r',
#            EkVovB,A1i,'b')
#axA[0].set_ylabel(r'$A_1$')
axA2.plot(EkVovB,GiE+Gi0-(GeE+Ge0),'k',
          EkVovB,Govn_iE+Govn_i0-(Govn_eE+Govn_e0),'k-.',
          [EkvovB_0_i]*2,ylim,'k:',
          -EkVovB,-GiE+Gi0-(-GeE+Ge0),'k',
          -EkVovB,-Govn_iE+Govn_i0-(-Govn_eE+Govn_e0),'k-.')
axA2.set_ylabel(r'$\Delta\Gamma/n$ [m/s]')
axA2.set_xlabel(r'$E_r/B$ [kV/T]')
axA2.set_ylim((-0.5,0.9))
axA2.legend(('DKES-Neotransp','model',r'$\tilde{E}=\tilde{E}_0$ for thermal ions'))
axA2.grid(True)
axA2.set_title(r'$y_0=$'+str(my_y0)+', $I_\mathrm{pos}/(I_\mathrm{pos}+|I_\mathrm{neg}|)=$'+'{:.4f}'.format(MaxwQ))
figA2.show()
figA2.savefig('lcf_'+equil+'_GammavsEr'+y0str+'.eps')

wait_for_user('press ENTER to continue to parameter scan')

#LOOP over gamma and y
#Ngamma=7
#gammaA=my_gamma*0.7
#gammaB=my_gamma*1.3
#dgamma=(gammaB-gammaA)/(Ngamma-1)
Nnustar0=109
nustar0A=my_nustar0*0.3
nustar0B=my_nustar0*6
Ny0=119
y0A=my_y0*0.2
y0B=my_y0*2.0
#dy0=(y0B-y0A)/(Ny0-1)
#(gammaM,y0M)=np.mgrid[slice(gammaA,gammaB + dgamma,dgamma),slice(y0A,y0B + dy0,dy0)]
#(gammaM,y0M)=np.mgrid[gammaA:gammaB:1j*Ngamma,y0A:y0B:1j*Ny0]
(nustar0M,y0M)=np.mgrid[nustar0A:nustar0B:1j*Nnustar0,y0A:y0B:1j*Ny0]
#print(gridy0)
#sys.exit()
#gammas=gammaM[:,0]
gamma=my_gamma
nustar0s=nustar0M[:,0]
y0s=y0M[0]

#EminM=np.nan*np.zeros((Ngamma,Ny0))
DeltaminM=np.nan*np.zeros((Nnustar0,Ny0))
ErootM=np.nan*np.zeros((Nnustar0,Ny0))
MaxwQM=np.nan*np.zeros((Nnustar0,Ny0))

for y0i in range(Ny0):
  for nustar0i in range(Nnustar0):
    #gamma=gammas[gammai]
    nustar0=nustar0s[nustar0i]
    y0=y0s[y0i]
    

    He=(1+(y_e/y0)**(my_alpha*gamma))**(-1/gamma)
    Hi=(1+(y_i/y0)**(my_alpha*gamma))**(-1/gamma)

    Ne=nustar0/nustar_e+Npl*(1+(nustar_e/my_nustarPS/Npl)**my_PSplexp)**(1/my_PSplexp)
    Ni=nustar0/nustar_i+Npl*(1+(nustar_i/my_nustarPS/Npl)**my_PSplexp)**(1/my_PSplexp)

    D11star_e=Ne[None,:]*He
    D11star_i=Ni[None,:]*Hi

    L11e=C_e*np.trapz((x**5*np.exp(-x**2))[None,:]*D11star_e,x=x,axis=1)
    L12e=C_e*np.trapz((x**7*np.exp(-x**2))[None,:]*D11star_e,x=x,axis=1)
    L11i=C_i*np.trapz((x**5*np.exp(-x**2))[None,:]*D11star_i,x=x,axis=1)
    L12i=C_i*np.trapz((x**7*np.exp(-x**2))[None,:]*D11star_i,x=x,axis=1)

    Govn_e=-L11e*A1e-L12e*A2e
    Govn_i=-L11i*A1i-L12i*A2i
    Govn_e0=-L11e*A1e0-L12e*A2e
    Govn_eE=-L11e*A1eE
    Govn_i0=-L11i*A1i0-L12i*A2i
    Govn_iE=-L11i*A1iE
    Govn_ep=Govn_e0+Govn_eE
    Govn_ip=Govn_i0+Govn_iE
    Govn_em=Govn_e0-Govn_eE
    Govn_im=Govn_i0-Govn_iE
    DGonv_p=Govn_ip-Govn_ep
    DGonv_m=Govn_im-Govn_em

    MaxwQM[nustar0i,y0i]=calcMaxw(EkVovB,DGonv_p,DGonv_m)
  
    DeltaGovn=Govn_i-Govn_e
    Eimin=np.argmin(DeltaGovn)
    if Eimin!=0 and Eimin!=len(EkVovB)-1:
      DeltaminM[nustar0i,y0i]=np.min(DeltaGovn)
    else:
      DeltaminM[nustar0i,y0i]=np.nan

    tmp=np.where(np.diff(np.sign(DeltaGovn))==2)[0]
    if len(tmp)>0:
      ErootM[nustar0i,y0i]=EkVovB[tmp[0]]
      
    if False:#True:#Ny0==1 and Ngamma==1:#False:#True:#False:#True:#False:  
      fig11,ax11=plt.subplots(1,1)
      for Ei in range(NE):
        ax11.loglog(nustar_e,D11star_e[Ei],'r',
                    nustar_i,D11star_i[Ei],'b')
      ax11.set_title('nustar0='+str(nustar0)+', y0='+str(y0))  
      fig11.show()

      figG,axG=plt.subplots(1,1)
      for Ei in range(NE):
        axG.semilogx(EkVovB,Govn_e,'r',
                 EkVovB,Govn_i,'b',
                 EkVovB,DeltaGovn,'k')
      axG.set_title('nustar0='+str(nustar0)+', y0='+str(y0))  
      figG.show()

      #print(tmp)
      wait_for_user()

epseffM=(nustar0M*(3*np.pi/4)**2)**(2/3)/2
my_epseff=(my_nustar0*(3*np.pi/4)**2)**(2/3)/2
lblchoice='epseff'
if lblchoice=='epseff':
  lblM=epseffM
  my_lblval=my_epseff
  lblstr=r'$\epsilon_\mathrm{eff}$'
else:
  lblM=nustar0M
  lblstr=r'$\nu^\star_0$'
  my_lblval=my_nustar0
if Nnustar0>2 and Ny0>2:      
  figop,axop=plt.subplots(1,1)
  if False:
    im=axop.pcolormesh(y0M,lblM,DeltaminM)
    figop.colorbar(im,ax=axop)
  else:
    Nlevels=25
    mx=np.max(np.abs(np.where(np.isnan(DeltaminM),0.0,DeltaminM)))
    print('mx='+str(mx))
    im=axop.contour(y0M,lblM,DeltaminM,Nlevels,cmap='jet',vmin=-mx,vmax=mx)
    #figop.colorbar(im,ax=axop)
    plt.clabel(im)
  axop.set_ylabel(lblstr)
  axop.set_xlabel(r'$y_0$')
  axop.set_title(r'min($Z\Gamma_i-\Gamma_e$)')
  axop.plot(my_y0,my_lblval,'r+')
  figop.show()
  figop.savefig('lcf_'+equil+'_search_DeltaGamma.eps')

  figEr,axEr=plt.subplots(1,1)
  im=axEr.pcolormesh(y0M,lblM,ErootM)
  figEr.colorbar(im,ax=axEr)
  axEr.set_ylabel(lblstr)
  axEr.set_xlabel(r'$y_0$')
  axEr.set_title(r'Electron root $E_r/B$ [kV/T]')
  axEr.plot(my_y0,my_lblval,'r+')
  figEr.show()
  figEr.savefig('lcf_'+equil+'_search_Er.eps')

  figMxw,axMxw=plt.subplots(1,1)
  im=axMxw.pcolormesh(y0M,lblM,MaxwQM)
  figMxw.colorbar(im,ax=axMxw)
  axMxw.set_ylabel(lblstr)
  axMxw.set_xlabel(r'$y_0$')
  axMxw.set_title(r'Maxwell quotient')
  axMxw.plot(my_y0,my_lblval,'r+')
  figMxw.show()
  figMxw.savefig('lcf_'+equil+'_search_MaxwQ.eps')

wait_for_user('press ENTER to end')
