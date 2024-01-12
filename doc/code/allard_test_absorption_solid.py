# -*- coding: utf-8 -*-
"""
Allard test example for test of full porous application

Created on Tue Jan 31 23:18:27 2023

@author: alexander
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyva.models as mds
import pyva.systems.infiniteLayers as iL
import pyva.properties.structuralPropertyClasses as stPC
import pyva.data.matrixClasses as mC
import pyva.properties.materialClasses as matC

# my packages

import pyva.properties.materialClasses as matC

#plt.close('all')



wb = pd.read_excel('data/elasto_absorber.xlsx')
omega_VA1 = 2*np.pi*wb.Frequency.values
abs_5_VA1 = wb['Foam55 treatment 5'].values
abs_DAF_VA1 = wb['Foam55 treatment 5'].values

# default fluid
air    = matC.Fluid(eta = 0.0)

# constants

foam     = matC.IsoMat(E=55000.,rho0=70,nu=0.4,eta=0.01)
foam20mm = stPC.PlateProp(0.02, foam)
foam10mm = stPC.PlateProp(0.01, foam)

# frequency and wavenumber
omega  = omega_VA1 # 2*np.pi*np.logspace(1,4,100)
omega0 = 2*np.pi*1000
angles  = np.linspace(0,np.pi/2)
angle45 = np.pi/4 
k1      = air.wavenumber(omega0)*np.sin(angles)

kx45    = air.wavenumber(omega)*np.sin(angle45)




z0     = air.z0

solid_foam = iL.SolidLayer(foam20mm)
solid_foam_1cm = iL.SolidLayer(foam10mm)


TMM_solid  = mds.TMmodel((solid_foam,))
#TMM_solid  = mds.TMmodel((solid_foam_1cm,solid_foam_1cm))


# Allard version
alpha_solid_DAF = TMM_solid.absorption_diffuse(omega,theta_max=np.pi/2,allard=True,signal = False,boundary_condition = 'fixed',)
alpha_solid_0   = TMM_solid.absorption(omega,kx=0,allard=True,boundary_condition = 'fixed',signal = False)
Z0_alla         = TMM_solid.impedance_allard(omega,kx=kx45,ID=0,boundary_condition = 'fixed',signal = False)





#%% solve


#%% plot impedance of solid layer
plt.figure(1)
plt.plot(omega,np.imag(Z0_alla),label='Im 0deg')
plt.plot(omega,np.real(Z0_alla),':',label='Re 0deg')


plt.xscale('log')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('z/Pa m s$^{-1}$')
plt.legend()

#%% plot abs for angle 0 deg
plt.figure(2)
plt.plot(omega,alpha_solid_0,label='0 degalla')
plt.plot(omega_VA1,abs_5_VA1,label='0-5 deg VA1')

plt.xscale('log')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('absorption')
plt.legend()

#%% plot abs for DAF
plt.figure(3)
plt.plot(omega,alpha_solid_DAF,label='alla')
plt.plot(omega_VA1,abs_DAF_VA1,label='DAF VA1')

plt.xscale('log')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('absorption')
plt.legend()
