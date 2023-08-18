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

wb1 = pd.read_excel('data/melamin_study.xlsx',sheet_name = "TL 78",usecols = "A:C",header = 0)
omega_VA1 = 2*np.pi*wb1.Frequency.values
freq_VA1 = wb1.Frequency.values
TL_melamin3cm_VA1 = wb1['30mm Melamin'].values
TL_melamin3cm_2_4kg_VA1 = wb1['30mm Melamin + 2.4kg'].values

wb1 = pd.read_excel('data/melamin_study.xlsx',sheet_name = "TL 78 (decoup)",usecols = "A:C",header = 0)
TL_melamin3cm_dec_VA1 = wb1['30mm Melamin'].values
TL_melamin3cm_2_4kg_dec_VA1 = wb1['30mm Melamin + 2.4kg'].values


wb2 = pd.read_excel('data/melamin_study.xlsx',sheet_name = "abs 78",usecols = "A:E",header = 0)
abs_melamin3cm_VA1 = wb2['melamin'].values
abs_melamin3cm_dec_VA1 = wb2['melamin decoupled'].values
abs_melamin3cm_2_4kg_VA1 = wb2['melamin-mass'].values
abs_melamin3cm_2_4kg_dec_VA1 = wb2['melamin-mass decoupled'].values



# default fluid
air    = matC.Fluid(rho0=1.21,eta = 0.0)
z0     = air.z0

freq  = np.geomspace(100, 16000,100)
omega = np.pi*2*freq
theta = 45/180*np.pi
kx45  =  air.wavenumber(omega)*np.sin(theta)

steel   = matC.IsoMat(E=2.0e11,rho0=7700,nu=0.27,eta=0.00)
alu     = matC.IsoMat(eta = 0.1)
rubber  = matC.IsoMat(E=2.6e6,rho0=1200,nu=0.49,eta=0.00)

ela_vac = matC.IsoMat(E=3.0e5,rho0=12.0,nu=0.4,eta=0.1) # Frame in vaccuum    E6
poroela = matC.PoroElasticMat(ela_vac, \
                            flow_res = 30000., \
                            porosity = 0.99, \
                            tortuosity = 1.01, \
                            length_visc = 250.E-6, length_therm = 550.E-6)


# Foams
foam_3cm = iL.PoroElasticLayer(poroela, 0.03)
foam_1_5cm = iL.PoroElasticLayer(poroela, 0.015)
foam_1cm = iL.PoroElasticLayer(poroela, 0.01)
foam_1um = iL.PoroElasticLayer(poroela, 0.00001)
foam_equiv = iL.MassLayer(0.03, 12.)




# metal
steel_1mm    = stPC.PlateProp(0.001, steel)
alu_1um    = stPC.PlateProp(1.E-6, alu) #
rubber_2mm    = stPC.PlateProp(0.002, rubber)

iL_steel_1mm = iL.PlateLayer(steel_1mm)
#iL_alu_1um = iL.PlateLayer(alu_1um)
iL_alu_1um = iL.SolidLayer(alu_1um)
iL_rubber_2mm = iL.SolidLayer(rubber_2mm)
heavy2_4kg   = iL.MassLayer(0.002, 1200)
iL_nothing = iL.MassLayer(1e-6,1.)

#layer_set1 = mds.TMmodel((iL_alu_1um,foam_3cm))
layer_set1 = mds.TMmodel((foam_3cm,iL_alu_1um))
layer_set1a = mds.TMmodel((foam_3cm,))
#layer_set1 = mds.TMmodel((foam_3cm,))
#layer_set2  = mds.TMmodel((heavy2_4kg,foam_3cm,iL_alu_1um)) # 
layer_set2  = mds.TMmodel((iL_rubber_2mm,foam_3cm,iL_alu_1um))
layer_set2a  = mds.TMmodel((iL_rubber_2mm,foam_3cm))


steel_set = mds.TMmodel((iL_steel_1mm,))
mass_set  = mds.TMmodel((foam_equiv,))

tau1 = layer_set1.transmission_diffuse(omega,theta_max=78/180*np.pi,allard=True)
tau1a = layer_set1a.transmission_diffuse(omega,theta_max=78/180*np.pi,allard=True)
tau2 = layer_set2.transmission_diffuse(omega,theta_max=78/180*np.pi,allard=True)
tau2a = layer_set2a.transmission_diffuse(omega,theta_max=78/180*np.pi,allard=True)

abs1 = layer_set1.absorption_diffuse(omega,theta_max=78/180*np.pi,allard=True)
abs2 = layer_set2.absorption_diffuse(omega,theta_max=78/180*np.pi,allard=True)


    
# Thickness

kl2,ka2,ks2,mu1,mu2,mu3 = poroela.wavenumbers(omega)
ka = np.sqrt(ka2)
kl = np.sqrt(kl2)

#%% Plot1
plt.figure(1)
plt.plot(freq,-10*np.log10(tau1.ydata[0]),'.:',label='TL_this alu + foam')
plt.plot(freq,-10*np.log10(tau1.ydata[0]),'.:',label='TL_foam')
plt.plot(freq_VA1,TL_melamin3cm_VA1,label='VA1 mela')
plt.plot(freq_VA1,TL_melamin3cm_dec_VA1,label='VA1 mela dec')
plt.xlabel('f/Hz')
plt.ylabel('TL/dB')
plt.xscale("log")
plt.legend()


#%% Plot2
plt.figure(2)
plt.plot(freq,-10*np.log10(tau2.ydata[0]),label='TL_mass_foam_thin alu')
plt.plot(freq,-10*np.log10(tau2a.ydata[0]),label='TL_mass_foam')
plt.plot(freq_VA1,TL_melamin3cm_2_4kg_VA1,label='VA1 mela + mass')
plt.plot(freq_VA1,TL_melamin3cm_2_4kg_dec_VA1,label='VA1 mela + mass dec')
plt.xlabel('f/Hz')
plt.ylabel('TL/dB')
plt.xscale("log")
plt.legend()

#%% Plot3
plt.figure(3)
plt.plot(freq,abs1.ydata[0],'.-',label='foam')
plt.plot(freq_VA1,abs_melamin3cm_VA1,label='VA1 mela')
plt.plot(freq_VA1,abs_melamin3cm_dec_VA1,label='VA1 mela dec')
plt.xlabel('f/Hz')
plt.ylabel('abs')
plt.xscale("log")
plt.yscale("log")
plt.legend()

#%% Plot4
plt.figure(4)
plt.plot(freq,abs2.ydata[0],label='mass_foam')
plt.plot(freq_VA1,abs_melamin3cm_2_4kg_VA1,label='VA1 mela + mass')
plt.plot(freq_VA1,abs_melamin3cm_2_4kg_dec_VA1,label='VA1 mela + mass dec')
plt.xlabel('f/Hz')
plt.ylabel('abs')
plt.xscale("log")
plt.yscale("log")
plt.legend()

