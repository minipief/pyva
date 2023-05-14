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

plt.close('all')

wb = pd.read_excel('data/elasto_absorber.xlsx')
omega_VA1 = 2*np.pi*wb.Frequency.values
abs_5_VA1 = wb['Foam55 + Alu 3mm 5'].values
abs_DAF_VA1 = wb['Foam55 + Alu 3mm 78'].values

omega = 2*np.pi*np.logspace(2,4,200)
omega0 = 2*np.pi*1000

# default fluid
air    = matC.Fluid(eta = 0.0)
foam     = matC.IsoMat(E=55000.,rho0=70,nu=0.4,eta=0.01)
foam50mm = stPC.PlateProp(0.05, foam)
alu = matC.IsoMat()
alu3mm = stPC.PlateProp(0.003, alu)

angles  = np.linspace(0,np.pi/2)
angle45 = np.pi/4
kx45    = air.wavenumber(omega)*np.sin(angle45)
#kx      = air.wavenumber(omega)*np.sin(angle45)
# default fluid
theta  = 50/180*np.pi
kx     = air.wavenumber(omega)*np.sin(theta)
omega1 = 6000
k1     = air.wavenumber(omega1)*np.sin(angles)



z0     = air.z0

plate_alu3mm = iL.PlateLayer(alu3mm)
solid_alu3mm = iL.SolidLayer(alu3mm)
foam_50mm    = iL.SolidLayer(foam50mm)

TMM_plate_foam = mds.TMmodel((plate_alu3mm, foam_50mm))
TMM_solid_foam = mds.TMmodel((solid_alu3mm, foam_50mm))


alpha_plate_foam_78  = TMM_plate_foam.absorption_diffuse(omega,theta_max=np.pi*78/180,theta_step = np.pi/180,signal = False,allard=True)
alpha_plate_foam_0   = TMM_plate_foam.absorption(omega,kx=0.0,signal = False,allard=True)
alpha_solid_foam_78  = TMM_solid_foam.absorption_diffuse(omega,theta_max=np.pi*78/180,theta_step = np.pi/180,signal = False,allard=True)
alpha_solid_foam_0   = TMM_solid_foam.absorption(omega,kx=0.0,signal = False,allard=True)

z_plate_foam      = TMM_plate_foam.impedance_allard(omega,kx=0.,ID=0,boundary_condition = 'fixed',signal = False)
z_solid_foam      = TMM_solid_foam.impedance_allard(omega,kx=0.,ID=0,boundary_condition = 'fixed',signal = False)
z_plate_foam_45   = TMM_plate_foam.impedance_allard(omega,kx=kx45,ID=0,boundary_condition = 'fixed',signal = False)
z_solid_foam_45   = TMM_solid_foam.impedance_allard(omega,kx=kx45,ID=0,boundary_condition = 'fixed',signal = False)

#D1,F = TMM_solid_fibre_10.allard_matrix(omega1, kx = k1,reduced = True,boundary_condition = 'fixed')





#%% solve


#%% plot impedance of fibre results for publishing
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(omega,np.real(z_plate_foam),label='Re plate foam')
plt.plot(omega,np.real(z_solid_foam),':',label='Re solid foam')
plt.xscale('log')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('z/Pa m s$^{-1}$')
plt.legend()

plt.subplot(2,1,2)
plt.plot(omega,np.imag(z_plate_foam),label='Im plate foam')
plt.plot(omega,np.imag(z_solid_foam),':',label='Im solid foam')
#plt.plot(omega,np.real(Z0_mass_fibre_10),'--',label='Re 10cm mass fibre Alla')
#plt.plot(omega,np.imag(Z0_mass_fibre_10),'--',label='Im 10cm mass fibre Alla')


plt.xscale('log')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('z/Pa m s$^{-1}$')
plt.legend()

#%% plot impedance of fibre results for publishing
plt.figure(2)
plt.title('45 degrees')
plt.subplot(2,1,1)
plt.plot(omega,np.real(z_plate_foam_45),label='Re plate foam')
plt.plot(omega,np.real(z_solid_foam_45),':',label='Re solid foam')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('z/Pa m s$^{-1}$')
plt.legend()

plt.subplot(2,1,2)
plt.plot(omega,np.imag(z_plate_foam_45),label='Im plate foam')
plt.plot(omega,np.imag(z_solid_foam_45),':',label='Im solid foam')
#plt.plot(omega,np.real(Z0_mass_fibre_10),'--',label='Re 10cm mass fibre Alla')
#plt.plot(omega,np.imag(Z0_mass_fibre_10),'--',label='Im 10cm mass fibre Alla')


plt.xscale('log')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('z/Pa m s$^{-1}$')
plt.legend()

#%% plot impedance of fibre results for publishing
plt.figure(3)
plt.plot(omega,alpha_plate_foam_0,label='plate foam')
plt.plot(omega,alpha_solid_foam_0,'.:',label='solid foam')
plt.plot(omega_VA1,abs_5_VA1,label='VA1')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('absorption')
plt.legend()

#%% plot impedance of fibre results for publishing
plt.figure(4)
plt.plot(omega,alpha_plate_foam_78,label='plate foam')
plt.plot(omega,alpha_solid_foam_78,'.-',label='solid foam')
plt.plot(omega_VA1,abs_DAF_VA1,label='VA1')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('absorption')
plt.legend()
