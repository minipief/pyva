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


# default fluid
air    = matC.Fluid(eta = 0.01)
z0     = air.z0

freq  = np.linspace(250, 10000,100)
omega = np.pi*2*freq
theta = 45/180*np.pi
kx45  =  air.wavenumber(omega)*np.sin(theta)

alu     = matC.IsoMat(E=7.2e10,rho0=2800,nu=0.3,eta=0.007)
ela_vac = matC.IsoMat(E=2.93e5,rho0=11.2,nu=0.2,eta=0.06) # Frame in vaccuum    E6
poroela = matC.PoroElasticMat(ela_vac, \
                            flow_res = 6600., \
                            porosity = 0.98, \
                            tortuosity = 1.03, \
                            length_visc = 200.E-6, length_therm = 380.E-6)

poroequiv = matC.EquivalentFluid(flow_res = 6600.,rho_bulk= 0.98*air.rho0+11.2, \
                            porosity = 0.98, \
                            tortuosity = 1.03, \
                            length_visc = 200.E-6, length_therm = 380.E-6)

    
foam_25_4mm = iL.PoroElasticLayer(poroela, 0.0254)
foam_25_4mm_equiv = iL.FluidLayer(0.0254,poroequiv)

alu_1_6mm = stPC.PlateProp(0.0016, alu)
#iL_alu_1_6mm   = iL.SolidLayer(alu_1_6mm)
iL_alu_1_6mm   = iL.PlateLayer(alu_1_6mm)

#layer_set = mds.TMmodel((iL_alu_1_6mm,foam_25_4mm))
#layer_set = mds.TMmodel((foam_25_4mm_equiv,iL_alu_1_6mm))
layer_set = mds.TMmodel((foam_25_4mm,iL_alu_1_6mm))

tau0 = layer_set.transmission_allard(omega,kx=0.)
tau = layer_set.transmission_diffuse(omega,theta_max=78/180*np.pi,allard=True)

    
# Thickness

kl2,ka2,ks2,mu1,mu2,mu3 = poroela.wavenumbers(omega)
ka = np.sqrt(ka2)
kl = np.sqrt(kl2)

#%% Plot1
plt.figure()
plt.plot(omega/2/np.pi,np.real(ka),'d',label='Re')
plt.plot(omega/2/np.pi,np.imag(ka),'+',label='Im')
plt.plot(omega/2/np.pi,np.real(kl),'d',label='Re')
plt.plot(omega/2/np.pi,np.imag(kl),'+',label='Im')
plt.legend()
plt.xlabel('f/Hz')
plt.ylabel('k/m^-1')

#%% Plot2
plt.figure(2)
tau.plot(2,res='dB')
# plt.plot(freq,np.real(Z0_ref/z0),label='Re')
# plt.plot(freq,np.imag(Z0_ref/z0),label='Im')
# plt.plot(freq,np.real(Z0.ydata[0]/z0),'+:',label='Re')
# plt.plot(freq,np.imag(Z0.ydata[0]/z0),'+:',label='Im')
# plt.legend()

