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






# Create props and materials
# Fluids
air = matC.Fluid()
poro_limp = matC.EquivalentFluid(rho_bulk = 30., \
                                flow_res = 40000., \
                                porosity = 0.94, \
                                tortuosity = 1.06, \
                                length_visc = 56.E-6, length_therm = 110.E-6,limp=True)
# Poroelastic
ela_vac = matC.IsoMat(4400000.0,130., 0., 0.1) # Frame in vaccuum    
poroela = matC.PoroElasticMat(ela_vac, \
                            flow_res = 40000., \
                            porosity = 0.94, \
                            tortuosity = 1.06, \
                            length_visc = 56.E-6, length_therm = 110.E-6)

# Solid
alu = matC.IsoMat()
alu1mm = stPC.PlateProp(0.001,alu)

# Layers
il_alu_solid_1mm = iL.SolidLayer(alu1mm,)
il_air_10mm      = iL.FluidLayer(0.01,fluid=air)
il_fibre_20mm    = iL.FluidLayer(0.02,fluid=poro_limp)
il_poro_ela      = iL.PoroElasticLayer(poroela, 0.05) 

# default fluid
theta  = 50/180*np.pi
kx     = air.wavenumber(omega)*np.sin(theta)
omega1 = 6000
angles = np.linspace(0,np.pi/2)
k1     = air.wavenumber(omega1)*np.sin(angles)


TM_all = mds.TMmodel((il_poro_ela,il_air_10mm,il_fibre_20mm,il_alu_solid_1mm))

V0,V1 = TM_all.V0()

v0,v1 = TM_all.V0(boundary_condition = 'fixed')



D0   = TM_all.allard_matrix(6000., kx = 0.,reduced = False)
D1,F = TM_all.allard_matrix(6000., kx = 0.,reduced = True)

Vs = D1.solve(F)

