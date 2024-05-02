# -*- coding: utf-8 -*-
"""
Biot inversion test

Step1: Check with numerical test case
Created on Fri Apr  5 09:36:11 2024

@author: FYR04BX
"""

# System
import copy

# Numerics
import numpy as np
import matplotlib.pyplot as plt

# pyva packages
import pyva.properties.materialClasses as matC
import pyva.models as mds
import pyva.systems.infiniteLayers as iL

# Initialize
plt.close('all')

omega_limits = 2*np.pi*np.array((100,2000))
omega = np.geomspace(omega_limits[0],omega_limits[-1],50)

# default fluid
air   = matC.Fluid(eta = 0.0)
fibre = matC.EquivalentFluid(porosity = 0.98, \
                               flow_res = 25000.,\
                               tortuosity = 1.02, \
                               length_visc = 90.e-6, \
                               length_therm = 180.e-6,\
                               rho_bulk = 0.98*1.20 + 30. , \
                               rho0 = air.rho0, \
                               dynamic_viscosity = air.dynamic_viscosity,
                               c0 = air.c0, kappa=air.kappa, heat_conductivity=air.heat_conductivity,\
                               limp = False)
    
# create similar material with  
fibre_limp = copy.deepcopy(fibre)
fibre_limp.limp = True
    
# Numerical test data
rho_eq = fibre.rho_freq(omega)
K_eq   = fibre.bulk_modulus(omega)
rho_eq_limp = fibre_limp.rho_freq(omega)
K_eq_limp   = fibre_limp.bulk_modulus(omega)

# recreate from this
fibre_inv = matC.EquivalentFluid.inverse_parameter_derivation(25000., 0.98, 0.98*1.20 + 30., air, rho_eq, K_eq, omega,limp=False)
fibre_inv_limp = matC.EquivalentFluid.inverse_parameter_derivation(25000., 0.98, 0.98*1.20 + 30., air, rho_eq_limp, K_eq_limp, omega,limp=True)
print(fibre_inv)

# Check and plot

#%% plot1
plt.close(1)
plt.figure(1)
plt.plot(omega,fibre_inv.tortuosity,label='inv')
plt.plot(omega,fibre_inv_limp.tortuosity,'--',label='inv limp')
plt.plot(omega_limits,[fibre.tortuosity]*2,'.',label='inv')
plt.xlabel('$\omega/s^{-1}$')
plt.ylabel('tortuosity')
delta = 0.1 #np.abs(fibre.tortuosity-fibre_inv.tortuosity[0])
plt.ylim((fibre.tortuosity-3*delta,fibre.tortuosity+3*delta))
plt.legend()

#%% plot2
plt.close(2)
plt.figure(2)
plt.plot(omega,fibre_inv.length_visc*1E6,label='inv visc')
plt.plot(omega,fibre_inv.length_therm*1E6,label='inv therm')
plt.plot(omega,fibre_inv_limp.length_visc*1E6,label='inv visc limp')
plt.plot(omega,fibre_inv_limp.length_therm*1E6,label='inv therm limp')
plt.plot(omega_limits,[fibre.length_visc*1E6]*2,'.',label='visc ref.')
plt.plot(omega_limits,[fibre.length_therm*1E6]*2,'d',label='therm ref.')
plt.xlabel('$\omega/s^{-1}$')
plt.ylabel('thermal lentgths/microns')
#plt.ylim((fibre.length_visc*0.8E6,fibre.length_visc*1.2E6))
plt.show()
plt.legend()


