# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 09:04:17 2022

"""

import numpy as np
import matplotlib.pyplot as plt

import pyva.systems.structure2Dsystems as st2Dsys
import pyva.properties.materialClasses as matC
import pyva.properties.structuralPropertyClasses as stPC

import pyva.systems.acousticRadiators as aR


# Define the properties
alu = matC.IsoMat() # Alu is default

# Plate constants
h    = 0.02
Lx   = 2
Ly   = 3
area = Lx*Ly
perimeter = 2*(Lx+Ly)

# Plate property
plate_prop    = stPC.PlateProp(h,alu)

plate = st2Dsys.Structure2DSystem(1, area, plate_prop, perimeter = perimeter)
plate_square = st2Dsys.Structure2DSystem(1, area, plate_prop, perimeter = 0)

plate_square.Lx

omega = np.geomspace(100,10000,200)*2*np.pi
freq  = omega/2/np.pi

rad_eff = plate.radiation_efficiency(omega)
rad_eff_simple = plate.radiation_efficiency_simple(omega)

# RectangularPlate
rec_plate = st2Dsys.RectangularPlate(2, Lx,Ly, prop=plate_prop)

# Excitation position 
x0 = 0.71
y0 = 1.22

_,Ns = rec_plate.get_modes_index(omega[-1])
N_max = Ns[-1]

# Modal displacement
w0 = rec_plate.w_modal_force(omega, N_max, 1., x0, y0, x0, y0)
# Infinite plate displacement
w0_inf = rec_plate.prop.w_inf(omega,0., 1.)





#%% Plot Radiation efficiency
plt.figure(1)
plt.plot(freq,rad_eff,label = 'Leppington')
plt.plot(freq,rad_eff_simple,label = 'ISO EN 12354-1')
plt.xscale('log')
plt.xlabel('$f/$Hz')
plt.ylabel('radiation efiiciency')
plt.legend()

#plt.savefig('../source/images/plate_radiation_efficiency.png')

#%% Plot Diplacement response
plt.figure(2)
plt.plot(freq,np.real(w0),label = 'Re w modal')
plt.plot(freq,np.imag(w0),label = 'Im w modal')
#plt.plot(freq,np.real(w0_inf),label = 'Re w inf')
plt.plot(freq,np.imag(w0_inf),label = 'Im w inf')
plt.xscale('log')
plt.ylim(-1e-7,1e-7)
plt.xlabel('$f/$Hz')
plt.ylabel('displacement / m')
plt.legend()

plt.savefig('../source/images/plate_point_displacement.png')

#%% Discrete use of rectangular plates

Lx = 0.8
Ly = 0.5
h  = 0.004

# Plate property 6mm
plate_prop  = stPC.PlateProp(h,alu)
plate_prop.coincidence_frequency()
rec_plate   = st2Dsys.RectangularPlate(3, Lx,Ly, prop=plate_prop, eta = 0.05)

half_air = aR.HalfSpace()
rad_mesh = half_air.get_mesh(omega[-1], Lx, Ly)

tau_inf   = rec_plate.transmission_coefficient_discrete(omega, (half_air,))
tau_modal = rec_plate.modal_transmission_coefficient_discrete(omega,(half_air,),)

mode_shapes,mesh = rec_plate.normal_modes(5000.)



#%% Plot TL
plt.figure(3)
plt.plot(freq,-10*np.log10(tau_modal),label = 'modal')
plt.plot(freq,-10*np.log10(tau_inf),label = 'inf')
plt.xscale('log')
plt.xlabel('$f/$Hz')
plt.ylabel('$TL$ / dB')
plt.legend()

plt.savefig('../source/images/plate_transmission_discrete.png')

#%% PLot mode

mode_shapes.plot3d(4,1)
plt.savefig('../source/images/plate_mode.png')