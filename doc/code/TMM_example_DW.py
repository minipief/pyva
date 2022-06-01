# -*- coding: utf-8 -*-
#
# Example with 3 tubes showing the capabilities of DynamicMatrix
import numpy as np
import matplotlib.pyplot as plt
import pyva.models as mds
import pyva.systems.infiniteLayers as iL
import pyva.properties.structuralPropertyClasses as stPC

# my packages

import pyva.properties.materialClasses as matC

plt.close('all')

omega = 2*np.pi*np.logspace(1.5,4,1000)


# default fluid
air    = matC.Fluid(eta = 0.01)
rho_bulk = 0.98*1.20 + 30.
fibre1 = matC.EquivalentFluid(porosity = 0.98, \
                               flow_res = 25000.,\
                               tortuosity = 1.02, \
                               length_visc = 90.e-6, \
                               length_therm = 180.e-6,\
                               rho_bulk = rho_bulk , \
                               rho0 = 1.208, \
                               dynamic_viscosity = 1.81e-5 )

# Create props
alu = matC.IsoMat()
alu1mm = stPC.PlateProp(0.001,alu)
    
    
# Thickness
omega1 = 6000
k1     = np.real(air.wavenumber(omega1))
angles = np.linspace(0,np.pi/2)
#kx     = k1*np.sin(angles)
#kx     = mC.DataAxis(kx,typestr='wavenumber')

theta  = 50/180*np.pi
kx     =  air.wavenumber(omega)*np.sin(theta)

z0     = air.z0

# Fluid layer
air_5cm  = iL.FluidLayer(0.05,air)
fibre_5cm  = iL.FluidLayer(0.05,fibre1)
# Limp mass layer
heavy_1kg = iL.MassLayer(0.001, 1000)
heavy_2kg7 = iL.MassLayer(0.001, 2700) 
# Plate layer
iL_alu1mm = iL.PlateLayer(alu1mm,)

#heavy_1kg = iL.PlateLayer(alu03mm)
T_alu            = mds.TMmodel((iL_alu1mm,))
T_alu_air_mass   = mds.TMmodel((iL_alu1mm,air_5cm,heavy_1kg))
T_alu_fibre_mass10 = mds.TMmodel((iL_alu1mm,fibre_5cm,heavy_1kg))
T_alu_fibre_mass27 = mds.TMmodel((iL_alu1mm,fibre_5cm,heavy_2kg7))

tau_alu = T_alu.transmission_diffuse(omega,signal = False)
tau_alu_air_mass = T_alu_air_mass.transmission_diffuse(omega,theta_step=np.pi/1000,signal = False)
tau_alu_fibre_mass10 = T_alu_fibre_mass10.transmission_diffuse(omega,signal = False)
tau_alu_fibre_mass27 = T_alu_fibre_mass27.transmission_diffuse(omega,signal = False)


#%% fig1
plt.figure(1)

plt.plot(omega,-10*np.log10(tau_alu),label='alu 1mm')
plt.plot(omega,-10*np.log10(tau_alu_air_mass),label='alu-air5cm-1kg')
plt.plot(omega,-10*np.log10(tau_alu_fibre_mass10),label='alu-fibre5cm-1kg')
plt.plot(omega,-10*np.log10(tau_alu_fibre_mass27),label='alu-fibre5cm-2.7kg')
plt.xscale('log')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('$TL$')
plt.legend()

plt.savefig('../source/images/TMM_DW_transmission.png')

#%%fig 2



