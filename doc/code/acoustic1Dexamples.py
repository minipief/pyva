# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 09:04:17 2022

"""

import numpy as np
import matplotlib.pyplot as plt

import pyva.systems.acoustic1Dsystems as ac1Dsys
import pyva.properties.materialClasses as matC

#Define the fluid
air = matC.Fluid()
L = 1.
S = 0.01

# Tube
tube = ac1Dsys.AcousticTube(L,air,S)

omega = np.linspace(100.,200.,11)
tube_elem = tube.acoustic_FE(omega,ID =[1,2])

# Mass stiffness
mass_stiffness_layer = ac1Dsys.MassStiffness(0.001, 2.)

rubber = matC.Fluid(1500,1.1)
T = 200    # tension
h = 0.001  # thickness
my_mem = ac1Dsys.Membrane(rubber,h,T,S)

# Perforate

# Perforate parameter
thickness = 0.002
holeR     = 0.0001
dist      = 0.01

perf = ac1Dsys.PerforatedLayer(thickness,holeR,distance = dist)
perf.porosity

omega = 2*np.pi*np.linspace(1.,1000.)
perf.plot(omega,res='real')
#plt.xscale('log')
#plt.yscale('log')
plt.figure()
plt.ylim((90000,130000))

#plt.savefig('../source/images/perforate_resistance.png')


# Helmholtz paramters
omega = 2*np.pi*np.linspace(100.,5000.,200)

# Perforate for HR
thickness = 0.0002 
holeR     = 0.0001
porosity  = 0.0072


# Helmholtz parameter
V         = 0.000001 
L         = 0.005
R         = 0.002
Ac        = np.pi*R**2

# Create perforate for top
perf_HR   = ac1Dsys.PerforatedLayer(thickness,holeR,Ac,porosity = dist)

# Create pure and covered HR
resPure  = ac1Dsys.HelmholtzResonator(V,L,R,air)
resPerf  = ac1Dsys.HelmholtzResonator(V,L,R,air,0.85,end_impedance=perf_HR.radiation_impedance)

Za_pure = resPure.radiation_impedance(omega)
Za_perf = resPerf.radiation_impedance(omega)

#%% HR plot

plt.figure()
plt.plot(omega,np.real(Za_pure),label = 'Re pure')
plt.plot(omega,np.imag(Za_pure),label = 'Im pure')
plt.plot(omega,np.real(Za_perf),label = 'Re perf')
plt.plot(omega,np.imag(Za_perf),label = 'Im perf')
plt.xlabel('$\omega/s^{-1}$')
plt.ylabel('radiation impedance')

plt.xscale('log')
plt.legend(loc=4)

#plt.savefig('../source/images/HR_acoustic_impedance.png')

# Quarterwave Parameters


quarter_perf  = ac1Dsys.QuarterWaveResonator(4*L,R,air,0,end_impedance=perf.radiation_impedance)
Za_perf = quarter_perf.radiation_impedance(omega)

#%% QWR plot

plt.figure()
plt.plot(omega,np.real(Za_perf),label = 'Re perf')
plt.plot(omega,np.imag(Za_perf),label = 'Im perf')
plt.xlabel('$\omega/s^{-1}$')
plt.ylabel('radiation impedance')

plt.xscale('log')
plt.legend(loc=4)

#plt.savefig('../source/images/QWR_acoustic_impedance.png')















