# -*- coding: utf-8 -*-
#
# Example with 3 tubes showing the capabilities of DynamicMatrix
import numpy as np
import matplotlib.pyplot as plt
import pyva.systems.infiniteLayers as iL
import pyva.properties.structuralPropertyClasses as stPC

# my packages

import pyva.properties.materialClasses as matC

plt.close('all')

omega = 2*np.pi*np.logspace(3,5,100)


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

tau_mass0  = heavy_2kg7.transmission_coefficient(omega,0.)
tau_mass30 = heavy_2kg7.transmission_coefficient(omega,30.*np.pi/180)

tau_plate0  = iL_alu1mm.transmission_coefficient(omega,0.)
tau_plate30 = iL_alu1mm.transmission_coefficient(omega,30.*np.pi/180)




#%% fig1
plt.figure(1)

plt.plot(omega,-10*np.log10(tau_mass0),label='mass 2.7kg m$^2$ 0 deg')
plt.plot(omega,-10*np.log10(tau_mass30),label='mass 2.7kg m$^2$ 30 deg')
plt.plot(omega,-10*np.log10(tau_plate0),':',label='alu 1mm  0 deg')
plt.plot(omega,-10*np.log10(tau_plate30),':',label='alu 1mm 30 deg')
plt.xscale('log')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('$TL$')
plt.legend()

plt.savefig('../source/images/infinite_layer_TL.png')

#%%fig 2



