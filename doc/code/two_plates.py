# -*- coding: utf-8 -*-
"""
Created on Sat Oct 10 16:55:57 2020

two plate script for book presentation 

"""
import numpy as np
import matplotlib.pyplot as plt

import pyva.coupling.junctions as jun
import pyva.properties.structuralPropertyClasses as stPC
import pyva.systems.structure2Dsystems as st2Dsys
import pyva.data.matrixClasses as mC
import pyva.properties.materialClasses as matC

import pyva.useful as uf

# my packages



# set print
printsw = False
#printsw = True

plt.close('all')
 

# x-axis tics
fc,fclabels = uf.get_3rd_oct_axis_labels()
fc,fclabels = fc[1:],fclabels[1:]

# Plate dimensions
Lx1 = 2.5
Lx2 = 1.7
Ly = 1.7
area1 = Lx1*Ly
area2 = Lx2*Ly

# Plate thickness
t1 = 2.0E-3;
t2 = 3.0E-3;

# junction properties
angle1 = 0
angle2 = 90*np.pi/180

# Create materials
alu = matC.IsoMat(nu=0.3,eta = 0.0)

# Create props
alu1mm = stPC.PlateProp(0.001,alu)
alu2mm = stPC.PlateProp(t1,alu)
alu3mm = stPC.PlateProp(t2,alu)

# Create plate subsystems
plate1 = st2Dsys.Structure2DSystem(1,area1,alu1mm)
plate2 = st2Dsys.Structure2DSystem(2,area2,alu2mm)

J12 = jun.LineJunction((plate1,plate2),Ly,(angle1,angle2))
#om = 2*np.pi*np.logspace(2,np.log10(10000),3)

print (J12)

dofs = J12.wave_DOF
J12.res_DOF
J12.exc_DOF
JM = J12.junction_matrix(np.array([1000,2000]))





# Frequency range
omega = mC.DataAxis.octave_band()
omega0 = 5000*2*np.pi
om    = omega.data
freq  = om/2/np.pi 


# show the complicated shape over angle
max_k = alu1mm.wavenumber_B(omega0)
min_k = alu2mm.wavenumber_L(omega0)

kx = np.linspace(0.,max_k,200)

tau5000 = J12.transmission_wavenumber(omega0,kx,(0,1), i_in_wave = (3,3,5,5) , i_out_wave= (5,3,5,3))
tau5000.plot(1)

plt.savefig('../source/images/line_junction_tau.png')

# check tau first
taus = J12.transmission_wavenumber_diffuse(omega.angular_frequency, (0,1), i_in_wave = (3,3) , i_out_wave= (5,3))
taus.plot(2)
plt.savefig('../source/images/line_junction_tau_diffuse.png')

