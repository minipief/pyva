# -*- coding: utf-8 -*-
"""
Allard test example for test of full porous application

Created on Tue Jan 31 23:18:27 2023

@author: alexander
"""

import numpy as np
import matplotlib.pyplot as plt
import pyva.models as mds
import pyva.systems.infiniteLayers as iL
import pyva.properties.structuralPropertyClasses as stPC

# my packages

import pyva.properties.materialClasses as matC

plt.close('all')

omega = 2*np.pi*np.logspace(3,5,100)


# default fluid
air    = matC.Fluid(eta = 0.01)
theta  = 50/180*np.pi
kx     = air.wavenumber(omega)*np.sin(theta)
omega1 = 6000
angles = np.linspace(0,np.pi/2)
k1     = air.wavenumber(omega1)*np.sin(angles)



# Create props
alu = matC.IsoMat()

# thicknesses
h1 = 0.001
alu1mm = stPC.PlateProp(0.001,alu)
il_alu_solid_1mm = iL.SolidLayer(alu1mm,)
il_alu_solid_2mm = iL.SolidLayer(alu1mm,)
il_alu_plate_1mm = iL.PlateLayer(alu1mm,)


T1 = il_alu_solid_1mm.transfer_impedance(omega1,kx)

T1M_solid = mds.TMmodel((il_alu_solid_1mm,))
T1M_plate = mds.TMmodel((il_alu_plate_1mm,))
T1M_solid_solid = mds.TMmodel((il_alu_solid_1mm,il_alu_solid_1mm))


V0,V1 = T1M_solid.V0()
v0,v1 = T1M_solid.V0(boundary_condition = 'fixed')
V01,V02 = T1M_solid_solid.V0(boundary_condition = 'fixed')

#D0 = T1M_solid.allard_matrix(air, omega1, kx = k1)
#D01 = T1M_solid_solid.allard_matrix(air, omega1, kx = k1, boundary_condition = 'fixed')
D0   = T1M_solid.allard_matrix(air, omega1, kx = k1,reduced = False)
D1,F = T1M_solid.allard_matrix(air, omega1, kx = k1,reduced = True)

V000 = D1.solve(F)

