# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 14:45:39 2022

@author: Johannes Seidel
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

# pyva packages
import pyva.coupling.junctions as con

import pyva.systems.structure2Dsystems as st2Dsys

import pyva.properties.structuralPropertyClasses as stPC
import pyva.properties.materialClasses as matC

import pyva.data.matrixClasses as mC

# plate dimensions
lx = 1
ly = 1
t  = 0.027

# material
steel = matC.IsoMat(E=210e9, nu=0.3, rho0=7890.)

# property
steel_27mm = stPC.PlateProp(t, steel, 'Steel t=27mm')
print(steel_27mm)

# Create plate subsystems
plate1 = st2Dsys.RectangularPlate(1, lx, ly, prop=steel_27mm, eta = 0.03)
plate2 = st2Dsys.RectangularPlate(2, lx, ly, prop=steel_27mm, eta = 0.03)

# Frequency range
omega = mC.DataAxis.octave_band(f_max=2*np.pi*10000)
om    = omega.data
freq  = om/2/np.pi 

steel_27mm.plot_wavenumbers(om)

k = np.zeros(len(om))
for idx, o in enumerate(om):
    k[idx] = steel_27mm.wavenumber_L(o) / steel_27mm.wavenumber_B(o)
       
plt.figure()
plt.plot(freq, k, label = r'$k_L/k_B$')
plt.xlabel(r'$f / [Hz]$')
plt.ylabel(r'$k_L/k_B$')
plt.legend()
plt.show()
    
    
#%% Figure 2
o = 2*np.pi*10000 # 10000Hz will give k_l/k_b = 0.3 for the 27mm thick steel plate

nSamples = 100

tau_BB = np.zeros(nSamples)
tau_BL = np.zeros(nSamples)
tau_BS = np.zeros(nSamples)

phi = np.linspace(0, np.pi, nSamples)

for idx, angle in enumerate(phi):

    lj = con.LineJunction( [plate1, plate2], 1, [0, angle])

    taus = lj.transmission_wavenumber_diffuse(omega      = o,
                                              i_sys      = (0,1), 
                                              i_in_wave  = (3,3,3),
                                              i_out_wave = (1,2,3),
                                              N_step     = 100,
                                              method     = 'langley',
                                              Signal     = False)
    
    tau_BL[idx] = taus[0,:]
    tau_BS[idx] = taus[1,:]
    tau_BB[idx] = taus[2,:]

phi = phi * 180/np.pi
 
plt.figure()
plt.title(r'Figure 2. Diffuse wave transmission coefficients for a "V" plate junction with subtended angle $\theta$ and $k_L$/$k_B$ = 0.3')
plt.plot(phi, tau_BB, label = r'$\tau^{12}_{BB}$', linestyle = '-')
plt.plot(phi, tau_BS, label = r'$\tau^{12}_{BS}$', linestyle = '--')
plt.plot(phi, tau_BL, label = r'$\tau^{12}_{BL}$', linestyle = '-.')
plt.xlabel(r'$\theta(degrees)$')
plt.ylabel(r'$\tau(\omega)$')
plt.xticks((0,45,90,135,180))
plt.yticks((0,0.25,0.5,0.75,1.0))
plt.xlim(0, 180)
plt.ylim(0, 1)
plt.legend()
plt.show()