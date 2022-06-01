# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 09:04:17 2022

"""

import numpy as np
import matplotlib.pyplot as plt
import pyva.properties.materialClasses as matC
import pyva.properties.geometricalPropertyClasses as geoPC
import pyva.properties.structuralPropertyClasses as stPC

#Define the solids
alu   = matC.IsoMat()
steel = matC.IsoMat(E=2.1e11,rho0=7850, nu = 0.3)
print(alu)

alu.G
steel.G

alu.G_complex

# Beam constants
rho  = 2700.
E    = 71.e9
nu   = 0.34

h    = 0.02
b    = 0.03

A    = h*b
Ix   = b**3*h/12
Iy   = b*h**3/12
Ixy  = 0.

beam_sec1 = geoPC.CrossSection(Ix, Iy, Ixy, A)
beam_sec2 = geoPC.RectBeam(h, b)


print(beam_sec1)
print(beam_sec2)

omega = np.geomspace(100,10000,5)

beam_prop = stPC.BeamProp(beam_sec2,alu)
print(beam_prop)

beam_prop.Bx
beam_prop.By

alu   = matC.IsoMat()
alu4mm = stPC.PlateProp(0.004, alu)

alu4mm.c_L()
alu4mm.c_S()
alu4mm.c_B_phase(omega)
alu4mm.c_B_group(omega)


#%% plot sound speed

# plt.figure(1)
# plt.plot(omega,np.real(c_limp),label='Re limp')
# plt.plot(omega,np.imag(c_limp),label='Im limp')
# plt.plot(omega,np.real(c_rigid),label='Re rigid')
# plt.plot(omega,np.imag(c_rigid),label='Im rigid')
# plt.xlabel('$\omega/$s$^{-1}$')
# plt.ylabel('$c/$m s$^{-1}$')
# plt.legend()

#plt.savefig('../source/images/equiv_sound_speed.png')

#%% plot density

# plt.figure(2)
# plt.plot(omega,np.real(rho_limp),label='Re limp')
# plt.plot(omega,np.imag(rho_limp),label='Im limp')
# plt.plot(omega,np.real(rho_rigid),label='Re rigid')
# plt.plot(omega,np.imag(rho_rigid),label='Im rigid')
# plt.xlabel('$\omega/$s$^{-1}$')
# plt.ylabel('$\\rho/$kg m$^{-3}$')
# plt.ylim(-50,50 )
# plt.legend()

#plt.savefig('../source/images/equiv_density.png')



    
    