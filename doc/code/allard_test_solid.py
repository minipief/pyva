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

#plt.close('all')

wb1 = pd.read_excel('data/fibre_test.xlsx')
omega_VA1 = 2*np.pi*wb1.frequency.values
TL_alu1mm_VA1 = wb1['alu1mm'].values
#TL_alu1mm_foam_1kg_VA1 = wb1['alu1mm+5cm_fibre-1kg/m^2'].values

omega0 = 2*np.pi*1000

omega = 2*np.pi*np.logspace(2,4.5,50)


# default fluid
air    = matC.Fluid(eta = 0.01)

# Create props
#alu = matC.IsoMat() 
alu = matC.IsoMat(E=7.1E10,rho0=2700.,nu=0.34,eta = 0.)
alu1mm    = stPC.PlateProp(0.001,alu)
alu0_5mm  = stPC.PlateProp(0.0005,alu)
    
    
# Thickness
omega1 = 6000
k1     = np.real(air.wavenumber(omega1))
angles = np.linspace(0,np.pi/2)
#kx     = k1*np.sin(angles)
#kx     = mC.DataAxis(kx,typestr='wavenumber')

# Angle definition
theta  = 50/180*np.pi
kx     =  air.wavenumber(omega)*np.sin(theta)

z0     = air.z0
SIF_out = aR.HalfSpace(air)

# Limp mass layer
heavy_1kg = iL.MassLayer(0.001, 1000)
# Plate layer
iL_alu1mm = iL.PlateLayer(alu1mm)
iL_alu1mm_solid   = iL.SolidLayer(alu1mm)
iL_alu05mm_solid  = iL.SolidLayer(alu0_5mm)
# Layup
T_alu_solid        = mds.TMmodel((iL_alu1mm_solid,))
T_alu_solid3       = mds.TMmodel((iL_alu05mm_solid,iL_alu05mm_solid))


# Matrix check
T1 = iL_alu1mm_solid.transfer_impedance(omega,kx)
T2 = iL_alu05mm_solid.transfer_impedance(omega,kx,ID = [1,2])
T3 = iL_alu05mm_solid.transfer_impedance(omega,kx,ID = [2,3])

T1A = T2.dot(T3)

diffT = T1.Dindex(20)-T1A.Dindex(20)
rel_T = np.abs(diffT)/np.abs(T1.Dindex(20))

A1,F1 = T_alu_solid.allard_matrix(omega,kx,boundary_condition='equivalent_fluid',reduced=True)
A3,F3 = T_alu_solid3.allard_matrix(omega,kx,boundary_condition='equivalent_fluid',reduced=True)

diffA = A1.Dindex(20)-A3.Dindex(20)
rel_A = np.abs(diffA)/np.abs(A1.Dindex(20))






