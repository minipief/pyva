# -*- coding: utf-8 -*-
#
# Example with 3 tubes showing the capabilities of DynamicMatrix
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyva.models as mds
import pyva.systems.acousticRadiators as aR
import pyva.systems.infiniteLayers as iL
import pyva.properties.structuralPropertyClasses as stPC

# my packages

import pyva.properties.materialClasses as matC

plt.close('all')

wb1 = pd.read_excel('data/fibre_test.xlsx')
omega_VA1 = 2*np.pi*wb1.frequency.values
TL_alu1mm_VA1 = wb1['alu1mm'].values
#TL_alu1mm_foam_1kg_VA1 = wb1['alu1mm+5cm_fibre-1kg/m^2'].values

omega0 = 2*np.pi*1000

omega = 2*np.pi*np.logspace(2,4.5,50)


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
#alu = matC.IsoMat() 
alu = matC.IsoMat(E=7.1E10,rho0=2700.,nu=0.34,eta = 0.)
alu1mm    = stPC.PlateProp(0.001,alu)
       
# Thickness
omega1 = 6000
k1     = np.real(air.wavenumber(omega1))
angles = np.linspace(0,np.pi/2)
#kx     = k1*np.sin(angles)
#kx     = mC.DataAxis(kx,typestr='wavenumber')

# Angle definition
theta  = 50/180*np.pi
kx     =  air.wavenumber(omega)*np.sin(theta)

z0     = air.z0
SIF_out = aR.HalfSpace(air)

# Fluid layer
fibre_5cm  = iL.FluidLayer(0.05,fibre1)
# Limp mass layer
heavy_1kg = iL.MassLayer(0.001, 1000)
# Plate layer
iL_alu1mm = iL.PlateLayer(alu1mm)
iL_alu1mm_solid   = iL.SolidLayer(alu1mm)

#heavy_1kg = iL.PlateLayer(alu03mm)
T_alu              = mds.TMmodel((iL_alu1mm,))
T_alu_solid        = mds.TMmodel((iL_alu1mm_solid,))

print('Starting single wall calculation')

z_p   = alu1mm.transfer_impedance(omega, kx)
z_out = SIF_out.radiation_impedance_wavenumber(omega, kx)
z_ref = z_p+z_out
p_ref = z_out/(z_ref)
z_in_plate = T_alu.impedance(omega,kx,boundary_condition='equivalent_fluid',signal = False)
z_in_plate_allard = T_alu.impedance_allard(omega,kx,boundary_condition='equivalent_fluid',signal = False)
z_in_solid_allard = T_alu_solid.impedance_allard(omega,kx,boundary_condition='equivalent_fluid',signal = False)

# Plane wave
tau_alu_50 = T_alu.transmission_coefficient(omega,kx,signal = False)
tau_alu_allard_50 = T_alu.transmission_allard(omega,kx,signal = False)
tau_alu_solid_50 = T_alu_solid.transmission_allard(omega,kx,signal = False)

# Diffuse field
tau_alu = T_alu.transmission_diffuse(omega,theta_step=np.pi/360,signal = False)
tau_alu_allard = T_alu.transmission_diffuse(omega,theta_step=np.pi/360,signal = False,allard = True)
tau_alu_solid = T_alu_solid.transmission_diffuse(omega,theta_step=np.pi/360,signal = False,allard = True)

#%% fig1
plt.figure(1)

plt.plot(omega,-10*np.log10(tau_alu_50),label='alu 1mm')
plt.plot(omega,-10*np.log10(tau_alu_allard_50),label='alu 1mm allard')
plt.plot(omega,-10*np.log10(tau_alu_solid_50),':',label='alu 1mm solid')

plt.xscale('log')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('$TL$')
plt.legend()

#%% fig2
plt.figure(2)

plt.plot(omega,-10*np.log10(tau_alu),label='alu 1mm')
plt.plot(omega,-10*np.log10(tau_alu_allard),label='alu 1mm allard')
plt.plot(omega,-10*np.log10(tau_alu_solid),':',label='alu 1mm solid')
#plt.plot(omega_VA1,TL_alu1mm_VA1,'d',label = 'VAOne')
plt.xscale('log')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('$TL$')
plt.legend()


#%%fig 3
plt.figure(3)
plt.subplot(2,1,1)
plt.plot(omega,np.real(z_ref),'d',label='Re ref')
plt.plot(omega,np.real(z_in_plate),'-',label='Re plate')
plt.plot(omega,np.real(z_in_plate_allard),'-.',label='Re plate allard')
plt.plot(omega,np.real(z_in_solid_allard),':',label='Re solid allard')
# plt.plot(omega2_VA1,abs_alu1mm_foam_1kg_VA1,'.',label = 'VAOne source')
# plt.plot(omega2_VA1,abs_1kg_foam_alu1mm_VA1,'.',label = 'VAOne receiver')
plt.xscale('log')
#plt.yscale('symlog')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel(r'$z$')
plt.legend()

plt.subplot(2,1,2)
plt.plot(omega,np.imag(z_ref),'d',label='Im ref')
plt.plot(omega,np.imag(z_in_plate),'-',label='Im plate')
plt.plot(omega,np.imag(z_in_plate_allard),'-.',label='Im plate allard')
plt.plot(omega,np.imag(z_in_solid_allard),':',label='Im solid allard')
# plt.plot(omega2_VA1,abs_alu1mm_foam_1kg_VA1,'.',label = 'VAOne source')
# plt.plot(omega2_VA1,abs_1kg_foam_alu1mm_VA1,'.',label = 'VAOne receiver')
plt.xscale('log')
#plt.yscale('symlog')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel(r'$z$')
plt.legend()


#%%fig 3
#plt.figure(3)

# plt.plot(omega,abs_alu1mm_foam_1kg,':',label='abs alu side')
# plt.plot(omega,abs_1kg_foam_alu1mm,':',label='abs 1kg side')
# plt.plot(omega2_VA1,abs_alu1mm_foam_1kg_VA1,'.',label = 'VAOne source')
# plt.plot(omega2_VA1,abs_1kg_foam_alu1mm_VA1,'.',label = 'VAOne receiver')
# plt.xscale('log')
# #plt.yscale('log')
# plt.xlabel('$\omega/$s$^{-1}$')
# plt.ylabel(r'$\alpha$')
# plt.legend()

