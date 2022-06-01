# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 09:04:17 2022

"""

import numpy as np
import matplotlib.pyplot as plt

import pyva.systems.acoustic3Dsystems as ac3Dsys
import pyva.properties.materialClasses as matC

#Define the fluid
air = matC.Fluid()

# Cavity Parameters
Lx = 6.
Ly = 4.
Lz = 3.

room = ac3Dsys.RectangularRoom(1, Lx, Ly, Lz, air)
print(room)

#%% Modal density
omega = np.geomspace(100,10000,num=64)

mod_dens        = room.modal_density(omega)
mod_dens_precise,om_c = room.modal_density_precise(omega)
 
#%% HR plot

plt.figure()
plt.plot(omega,mod_dens,label = 'estimated')
plt.plot(om_c,mod_dens_precise,label = 'precise')
plt.xlabel('$\omega/s^{-1}$')
plt.ylabel('modal density')

plt.xscale('log')
plt.yscale('log')
plt.legend(loc=2)

#plt.savefig('../source/images/room_modal_density.png')
















