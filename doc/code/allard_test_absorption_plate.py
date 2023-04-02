# -*- coding: utf-8 -*-
"""
Allard test example for test of full porous application

Created on Tue Jan 31 23:18:27 2023

@author: alexander
"""

import numpy as np
import matplotlib.pyplot as plt
import pyva.models as mds
import pyva.systems.infiniteLayers as iL
import pyva.properties.structuralPropertyClasses as stPC
import pyva.data.matrixClasses as mC
import pyva.properties.materialClasses as matC

# my packages

import pyva.properties.materialClasses as matC

plt.close('all')

omega = 2*np.pi*np.logspace(1,4,100)
omega0 = 2*np.pi*1000

# default fluid
air    = matC.Fluid(eta = 0.0)
fibre1 = matC.EquivalentFluid(porosity = 0.98, \
                               flow_res = 25000.,\
                               tortuosity = 1.02, \
                               length_visc = 90.e-6, \
                               length_therm = 180.e-6,\
                               rho_bulk = 0.98*1.20 + 30. , \
                               rho0 = 1.208, \
                               dynamic_viscosity = 1.81e-5 )
# Thickness
h1 = 0.1
h2 = 0.2

alu = matC.IsoMat()
alu2mm = stPC.PlateProp(0.002, alu)

angles  = np.linspace(0,np.pi/2)
angle45 = np.pi/4 
#kx45    = air.wavenumber(omega)*np.sin(angle45)
#kx      = air.wavenumber(omega)*np.sin(angle45)

z0     = air.z0

fibre_10cm  = iL.FluidLayer(h1,fibre1)
perforate   = iL.PerforatedLayer(0.005, 0.001, distance = 0.02)
plate_alu2mm = iL.PlateLayer(alu2mm)

TMM_fibre_10      = mds.TMmodel((fibre_10cm,))
TMM_plate_fibre_10 = mds.TMmodel((plate_alu2mm, fibre_10cm,))


# Old version 
alpha_fibre_10  = TMM_fibre_10.absorption_diffuse(omega,theta_max=np.pi/2,signal = False)
z_fibre_10      = TMM_fibre_10.impedance(omega,kx=0.,ID=0,boundary_condition = 'fixed')
alpha_plate_fibre_10  = TMM_plate_fibre_10.absorption_diffuse(omega,theta_max=np.pi/2,signal = False)
z_plate_fibre_10      = TMM_plate_fibre_10.impedance(omega,kx=0.,ID=0,boundary_condition = 'fixed')

# Allard version
alpha_fibre_10_alla  = TMM_fibre_10.absorption_diffuse(omega,theta_max=np.pi/2,allard=True,signal = False)
Z0_fibre_10 = TMM_fibre_10.impedance_allard(omega,kx=0.,ID=0,boundary_condition = 'fixed',signal = False)
alpha_plate_fibre_10_alla  = TMM_plate_fibre_10.absorption_diffuse(omega,theta_max=np.pi/2,allard=True,signal = False)
Z0_plate_fibre_10 = TMM_plate_fibre_10.impedance_allard(omega,kx=0.,ID=0,boundary_condition = 'fixed',signal = False)




# default fluid
theta  = 50/180*np.pi
kx     = air.wavenumber(omega)*np.sin(theta)
omega1 = 6000
angles = np.linspace(0,np.pi/2)
k1     = air.wavenumber(omega1)*np.sin(angles)



#%% solve


#%% plot impedance of fibre results for publishing
plt.figure(1)
#plt.plot(omega,alpha_fibre_20,label='20cm fibre')
plt.plot(omega,np.real(z_fibre_10.ydata.flatten()),label='Re 10cm fibre')
plt.plot(omega,np.imag(z_fibre_10.ydata.flatten()),label='Im 10cm fibre')
plt.plot(omega,np.real(Z0_fibre_10),label='Re 10cm fibre Alla')
plt.plot(omega,np.imag(Z0_fibre_10),label='Im 10cm fibre Alla')
plt.xscale('log')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('z/Pa m s$^{-1}$')
plt.legend()


plt.figure(2)
plt.plot(omega,np.real(z_plate_fibre_10.ydata.flatten()),label='Re 10cm fibre')
plt.plot(omega,np.imag(z_plate_fibre_10.ydata.flatten()),label='Im 10cm fibre')
plt.plot(omega,np.real(Z0_plate_fibre_10),label='Re 10cm plate fibre Alla')
plt.plot(omega,np.imag(Z0_plate_fibre_10),label='Im 10cm plate fibre Alla')

plt.xscale('log')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('z/Pa m s$^{-1}$')
plt.legend()

#%% plot impedance of fibre results for publishing
plt.figure(3)
plt.plot(omega,alpha_fibre_10,label='10cm fibre')
plt.plot(omega,alpha_fibre_10_alla,label='10cm fibre alla')
plt.plot(omega,alpha_plate_fibre_10,label='10cm plate fibre')
plt.plot(omega,alpha_plate_fibre_10_alla,label='10cm plate fibre alla')
#plt.plot(omega,alpha_plate_fibre_10,label='perforate + 10cm fibre' )

plt.xscale('log')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('absorption')
plt.legend()

