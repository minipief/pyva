# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 09:04:17 2022

"""

import numpy as np
import matplotlib.pyplot as plt

import pyva.properties.materialClasses as matC

plt.close('all')

#Define the fluid
air   = matC.Fluid()
water = matC.Fluid(c0=1500, rho0=1000,dynamic_viscosity=1.0087,eta = 0.) 

omega = np.geomspace(100,10000,5)

air.c_freq(omega)
air.c_freq()

z1 = air.impedance()*1.5

reflection_factor = air.reflection_factor(omega, z1, theta = 0)
absorption_coefficent = air.absorption(omega, z1, theta = 0) 
absorption_diffuse    = air.absorption_diffuse(omega, z1)

# Equivalent fluid
air    = matC.Fluid(eta = 0.0)
fibre_limp = matC.EquivalentFluid(porosity = 0.98, \
                               flow_res = 25000.,\
                               tortuosity = 1.02, \
                               length_visc = 90.e-6, \
                               length_therm = 180.e-6,\
                               rho_bulk = 31.176 , \
                               rho0 = 1.208, \
                               dynamic_viscosity = 1.81e-5 )
    
fibre_rigid = matC.EquivalentFluid(porosity = 0.98, \
                               flow_res = 25000.,\
                               tortuosity = 1.02, \
                               length_visc = 90.e-6, \
                               length_therm = 180.e-6,\
                               rho_bulk = 31.176 , \
                               limp = False, \
                               rho0 = 1.208, \
                               dynamic_viscosity = 1.81e-5 )

omega = np.geomspace(100,10000,100)

c_limp  = fibre_limp.c_freq(omega)
c_rigid = fibre_rigid.c_freq(omega)

rho_limp  = fibre_limp.rho_freq(omega)
rho_rigid = fibre_rigid.rho_freq(omega)

#%% plot sound speed

plt.figure(1)
plt.plot(omega,np.real(c_limp),'-',label='Re limp')
plt.plot(omega,np.imag(c_limp),':',label='Im limp')
plt.plot(omega,np.real(c_rigid),'.-',label='Re rigid')
plt.plot(omega,np.imag(c_rigid),'--',label='Im rigid')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('$c/$m s$^{-1}$')
plt.legend()

plt.savefig('../source/images/equiv_sound_speed.png')

#%% plot density

plt.figure(2)
plt.plot(omega,np.real(rho_limp),'-',label='Re limp')
plt.plot(omega,np.imag(rho_limp),':',label='Im limp')
plt.plot(omega,np.real(rho_rigid),'.-',label='Re rigid')
plt.plot(omega,np.imag(rho_rigid),'--',label='Im rigid')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('$\\rho/$kg m$^{-3}$')
plt.ylim(-50,50 )
plt.legend()

plt.savefig('../source/images/equiv_density.png')



    
    