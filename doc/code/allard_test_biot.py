# -*- coding: utf-8 -*-
#
# Example with 
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

# wb1 = pd.read_excel('data/fibre_test.xlsx')
# omega_VA1 = 2*np.pi*wb1.frequency.values
# TL_alu1mm_VA1 = wb1['alu1mm'].values
#TL_alu1mm_foam_1kg_VA1 = wb1['alu1mm+5cm_fibre-1kg/m^2'].values

# default fluid
air    = matC.Fluid(eta = 0.01)
z0     = air.z0

freq  = np.linspace(200, 1500,100)
omega = np.pi*2*freq
theta = 45/180*np.pi
kx45  =  air.wavenumber(omega)*np.sin(theta)


rho_bulk = 0.98*1.20 + 30.
# fibre1 = matC.EquivalentFluid(porosity = 0.98, \
#                                flow_res = 25000.,\
#                                tortuosity = 1.02, \
#                                length_visc = 90.e-6, \
#                                length_therm = 180.e-6,\
#                                rho_bulk = rho_bulk , \
#                                rho0 = 1.208, \
#                                dynamic_viscosity = 1.81e-5 )


# Domisol Coffrage nu = 0 ??? page 124 Allard [All2009]_
E = 2*2200000.*(1+0.) # Calculate E from G in reference 
ela_vac = matC.IsoMat(E,130., 0., 0.1) # Frame in vaccuum    
poroela = matC.PoroElasticMat(ela_vac, \
                            flow_res = 40000., \
                            porosity = 0.94, \
                            tortuosity = 1.06, \
                            length_visc = 56.E-6, length_therm = 110.E-6)
    
poro_rigid = matC.EquivalentFluid(rho_bulk = poroela.rho_bulk, \
                            flow_res = 40000., \
                            porosity = 0.94, \
                            tortuosity = 1.06, \
                            length_visc = 56.E-6, length_therm = 110.E-6,limp=False)

poro_limp = matC.EquivalentFluid(rho_bulk = poroela.rho_bulk, \
                                flow_res = 40000., \
                                porosity = 0.94, \
                                tortuosity = 1.06, \
                                length_visc = 56.E-6, length_therm = 110.E-6,limp=True)

    
# Create Layers    
poro_10cm = iL.PoroElasticLayer(poroela, 0.1) # Change order of paramters
equiv_limp_10cm = iL.FluidLayer(0.1,fluid=poro_limp) 
equiv_rigid_10cm = iL.FluidLayer(0.1,fluid=poro_rigid) 
    
il_poro = mds.TMmodel((poro_10cm,))
il_equiv_limp = mds.TMmodel((equiv_limp_10cm,))
il_equiv_rigid = mds.TMmodel((equiv_rigid_10cm,))

abs0 = il_poro.absorption(omega, kx = 0., allard=True)

Z0 = il_poro.impedance_allard(omega, kx = 0.)
Z0_ref = poroela.surface_impedances(omega, 0.1)
Z0_equiv_limp = il_equiv_limp.impedance(omega, kx = 0.)
Z0_equiv_rigid = il_equiv_rigid.impedance(omega, kx = 0.)

Z45 = il_poro.impedance_allard(omega, kx = kx45)

    
# Thickness

kl2,ka2,ks2,mu1,mu2,mu3 = poroela.wavenumbers(omega)
Zf1,Zf2,Zs1,Zs2 = poroela.impedances(omega)

ka = np.sqrt(ka2)
kl = np.sqrt(kl2)

#%% Plot1
plt.figure(1)
plt.plot(omega/2/np.pi,np.real(ka),'d',label='Re 1st')
plt.plot(omega/2/np.pi,np.imag(ka),'+',label='Im 1st')
plt.plot(omega/2/np.pi,np.real(kl),'d',label='Re 2nd')
plt.plot(omega/2/np.pi,np.imag(kl),'+',label='Im 2nd')
plt.legend()
plt.title('compressional wave-numbers')
plt.xlabel('f/Hz')
plt.ylabel('k_a/m^-1')

#%% Plot2
plt.figure(2)
plt.plot(freq,np.real(Zf1/z0),label='Re Z_f1')
plt.plot(freq,np.imag(Zf1/z0),label='Im Z_f1')
plt.plot(freq,np.real(Zf2/z0),label='Re Z_f2')
plt.plot(freq,np.imag(Zf2/z0),label='Im Z_f2')
plt.legend()
plt.xlabel('f/Hz')
plt.ylabel('Z_f/z0')

#%% Plot3
plt.figure(3)
plt.plot(freq,np.real(Zs1/z0),label='Re Zs1')
plt.plot(freq,np.imag(Zs1/z0),label='Im Zs1')
plt.plot(freq,np.real(Zs2/z0),label='Re Zs2')
plt.plot(freq,np.imag(Zs2/z0),label='Im Zs2')
plt.legend()
plt.xlabel('f/Hz')
plt.ylabel('Z_s/z0')


#%% Plot3
plt.figure(4)
plt.plot(freq,np.real(Z0_ref/z0),label='Re')
plt.plot(freq,np.imag(Z0_ref/z0),label='Im')
#plt.plot(freq,np.real(Z0.ydata[0]/z0),'+:',label='Re ref')
#plt.plot(freq,np.imag(Z0.ydata[0]/z0),'+:',label='Im ref')
plt.plot(freq,np.real(Z0_equiv_limp.ydata[0]/z0),':',label='Re porous limp')
plt.plot(freq,np.imag(Z0_equiv_limp.ydata[0]/z0),':',label='Im porous limp')
plt.plot(freq,np.real(Z0_equiv_rigid.ydata[0]/z0),'-.',label='Re porous rigid')
plt.plot(freq,np.imag(Z0_equiv_rigid.ydata[0]/z0),'-.',label='Im porous rigid')

plt.legend()


#%% Plot4
plt.figure(5)
plt.plot(freq,np.real(Z45.ydata[0]/z0),'+:',label='Re')
plt.plot(freq,np.imag(Z45.ydata[0]/z0),'+:',label='Im')
plt.xlim((200,1420))
plt.ylim((-3,3))
plt.legend()



