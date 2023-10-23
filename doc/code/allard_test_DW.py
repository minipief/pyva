# -*- coding: utf-8 -*-
#
# Example with 3 tubes showing the capabilities of DynamicMatrix
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyva.models as mds
import pyva.systems.infiniteLayers as iL
import pyva.properties.structuralPropertyClasses as stPC
import pyva.properties.materialClasses as matC

plt.close('all')

wb1 = pd.read_excel('data/fibre_test.xlsx')
omega_VA1 = 2*np.pi*wb1.frequency.values
TL_alu1mm_VA1 = wb1['alu1mm'].values
#TL_alu1mm_foam_1kg_VA1 = wb1['alu1mm+5cm_fibre-1kg/m^2'].values

wb2 = pd.read_excel('data/aluminium_foam_1kg.xlsx')
omega2_VA1 = 2*np.pi*wb2.Frequency.values
abs_alu1mm_foam_1kg_VA1 = wb2['Absorption Source'].values
abs_1kg_foam_alu1mm_VA1 = wb2['Absorption Receiver'].values
TL_alu1mm_foam_1kg_VA1 = wb2['TL'].values



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
fibre_5cm  = iL.FluidLayer(0.05,fibre1)
# Limp mass layer
heavy_1kg = iL.MassLayer(0.001, 1000)
# Plate layer
iL_alu1mm = iL.PlateLayer(alu1mm)
iL_alu1mm_solid = iL.SolidLayer(alu1mm)

#heavy_1kg = iL.PlateLayer(alu03mm)
T_alu              = mds.TMmodel((iL_alu1mm,))
T_alu_solid        = mds.TMmodel((iL_alu1mm_solid,))
T_alu_plate_fibre_mass = mds.TMmodel((iL_alu1mm,fibre_5cm,heavy_1kg))
T_alu_solid_fibre_mass = mds.TMmodel((iL_alu1mm_solid,fibre_5cm,heavy_1kg))
T_mass_solid_fibre_alu = mds.TMmodel((heavy_1kg,fibre_5cm,iL_alu1mm_solid))

print('Starting single wall calculation')
tau_alu = T_alu.transmission_diffuse(omega,signal = False)
tau_alu_allard = T_alu.transmission_diffuse(omega,signal = False,allard = True)
tau_alu_solid = T_alu_solid.transmission_diffuse(omega,signal = False,allard = True)
print('Starting double wall calculation')
tau_alu_plate_fibre_mass = T_alu_plate_fibre_mass.transmission_diffuse(omega,theta_max = 78/180*np.pi,signal = False)
tau_alu_plate_fibre_mass_alla = T_alu_plate_fibre_mass.transmission_diffuse(omega,theta_max = 78/180*np.pi,signal = False,allard = True)
tau_alu_solid_fibre_mass = T_alu_solid_fibre_mass.transmission_diffuse(omega,theta_max = 78/180*np.pi,signal = False,allard = True)
tau_mass_solid_fibre_alu = T_mass_solid_fibre_alu.transmission_diffuse(omega,theta_max = 78/180*np.pi,signal = False,allard = True)
abs_alu1mm_foam_1kg = T_alu_solid_fibre_mass.absorption_diffuse(omega,theta_max = 78/180*np.pi,signal = False,allard = True)
abs_1kg_foam_alu1mm = T_mass_solid_fibre_alu.absorption_diffuse(omega,theta_max = 78/180*np.pi,signal = False,allard = True)

#%% fig1
plt.figure(1)

plt.plot(omega,-10*np.log10(tau_alu),label='alu 1mm')
plt.plot(omega,-10*np.log10(tau_alu_allard),label='alu 1mm allard')
plt.plot(omega,-10*np.log10(tau_alu_solid),':',label='alu 1mm solid')
plt.plot(omega_VA1,TL_alu1mm_VA1,'d',label = 'VAOne')
plt.xscale('log')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('$TL$')
plt.legend()


#%%fig 2
plt.figure(2)

plt.plot(omega,-10*np.log10(tau_alu_plate_fibre_mass),label='plate')
plt.plot(omega,-10*np.log10(tau_alu_plate_fibre_mass_alla),label='plate Allard')
plt.plot(omega,-10*np.log10(tau_alu_solid_fibre_mass),'.-',label='alu-fibre5cm-1kg solid')
plt.plot(omega,-10*np.log10(tau_mass_solid_fibre_alu),'<:',label='1kg-fibre5cm-alu solid')
plt.plot(omega2_VA1,TL_alu1mm_foam_1kg_VA1,'d',label = 'VAOne')
plt.xscale('log')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('$TL$')
plt.legend()

#%%fig 3
plt.figure(3)

plt.plot(omega,abs_alu1mm_foam_1kg,':',label='abs alu side')
plt.plot(omega,abs_1kg_foam_alu1mm,':',label='abs 1kg side')
plt.plot(omega2_VA1,abs_alu1mm_foam_1kg_VA1,'.',label = 'VAOne source')
plt.plot(omega2_VA1,abs_1kg_foam_alu1mm_VA1,'.',label = 'VAOne receiver')
plt.xscale('log')
#plt.yscale('log')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel(r'$\alpha$')
plt.legend()

