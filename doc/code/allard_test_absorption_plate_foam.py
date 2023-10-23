# -*- coding: utf-8 -*-
"""
Allard example of section 11.7.2 of [All2009]_ 
Example with lay-up given in figure 11.16 and table 11.7. 

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

plt.close('all')

freq  = np.linspace(200,1400,50)
omega = 2*np.pi*freq

# default fluid
air    = matC.Fluid(eta = 0.0)
carpet_solid = matC.IsoMat(E=20000.,rho0=60,nu=0.,eta=0.5)
carpet = matC.PoroElasticMat(carpet_solid, \
                            flow_res = 5000., \
                            porosity = 0.99, \
                            tortuosity = 1., \
                            length_visc = 23.E-6, length_therm = 28.E-6)
    
screen_mat = matC.IsoMat(E=30000.,rho0 = 2000, nu=0.49)
screen     = mat

fibre_solid = matC.IsoMat(E=100000.,rho0=60,nu=0.,eta=0.88)
fibre = matC.PoroElasticMat(fibre_solid, \
                            flow_res = 33000., \
                            porosity = 0.98, \
                            tortuosity = 1.1, \
                            length_visc = 50.E-6, length_therm = 110.E-6)

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
