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

air_20C_980 = matC.Fluid.air(294.15,0.980)
air_21C_1013 = matC.Fluid.air(295.15,1.023)
air_10C_1000 = matC.Fluid.air(284.15,1.0)


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

fibre_delany = matC.DelanyBazley(flow_res = 25000.,\
                               rho0 = 1.208, \
                               dynamic_viscosity = 1.81e-5 )

fibre_miki = matC.DelanyBazley(flow_res = 25000.,\
                               rho0 = 1.208, \
                               dynamic_viscosity = 1.81e-5,miki = True )

    
omega = np.geomspace(100,10000,100)

c_limp  = fibre_limp.c_freq(omega)
c_rigid = fibre_rigid.c_freq(omega)
c_delany = fibre_delany.c_freq(omega)
c_miki  = fibre_miki.c_freq(omega)

rho_limp  = fibre_limp.rho_freq(omega)
rho_rigid = fibre_rigid.rho_freq(omega)
rho_delany = fibre_delany.rho_freq(omega)
rho_miki  = fibre_miki.rho_freq(omega)


#%% plot sound speed
plt.close(1)
plt.figure(1)
plt.plot(omega,np.real(c_limp),'-',label='Re limp',color = 'C0')
plt.plot(omega,np.imag(c_limp),'--',label='Im limp',color = 'C0')
plt.plot(omega,np.real(c_rigid),'-',label='Re rigid',color = 'C1')
plt.plot(omega,np.imag(c_rigid),'--',label='Im rigid',color = 'C1')
plt.plot(omega,np.real(c_delany),'-',label='Re Delany',color = 'C2')
plt.plot(omega,np.imag(c_delany),'--',label='Im Delany',color = 'C2')
plt.plot(omega,np.real(c_miki),'-',label='Re Miki',color = 'C3')
plt.plot(omega,np.imag(c_miki),'--',label='Im Miki',color = 'C3')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('$c/$m s$^{-1}$')
plt.legend()

plt.savefig('../source/images/equiv_sound_speed.png')

#%% plot density
plt.close(2)
plt.figure(2)
plt.plot(omega,np.real(rho_limp),'-',label='Re limp',color = 'C0')
plt.plot(omega,np.imag(rho_limp),'--',label='Im limp',color = 'C0')
plt.plot(omega,np.real(rho_rigid),'-',label='Re rigid',color = 'C1')
plt.plot(omega,np.imag(rho_rigid),'--',label='Im rigid',color = 'C1')
plt.plot(omega,np.real(rho_delany),'-',label='Re Delany',color = 'C2')
plt.plot(omega,np.imag(rho_delany),'--',label='Im Delany',color = 'C2')
plt.plot(omega,np.real(rho_miki),'-',label='Re Miki',color = 'C3')
plt.plot(omega,np.imag(rho_miki),'--',label='Im Miki',color = 'C3')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('$\\rho/$kg m$^{-3}$')
plt.ylim(-50,50 )
plt.legend()

plt.savefig('../source/images/equiv_density.png')



    
    