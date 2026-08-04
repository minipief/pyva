# -*- coding: utf-8 -*-
"""
Pyva box example without isolation
"""

import numpy as np
import matplotlib.pyplot as plt

import pyva.models as mds
import pyva.coupling.junctions as jun
import pyva.properties.structuralPropertyClasses as stPC
import pyva.systems.structure2Dsystems as st2Dsys
import pyva.systems.acoustic3Dsystems as ac3Dsys
import pyva.loads.loadCase as lC

import pyva.useful as uf

# my packages

import pyva.data.dof as dof
import pyva.data.matrixClasses as mC
import pyva.properties.materialClasses as matC

plt.close('all')

# x-axis tics for better readability
fc,fclabels = uf.get_3rd_oct_axis_labels()
fc,fclabels = fc[1:],fclabels[1:]

# Frequency range
omega = mC.DataAxis.octave_band(f_max=2*np.pi*10000)

# Plate dimensions
Lx = 1.2
Ly = 1
Lz = 1.1

# Box dimensions
V = Lx*Ly*Lz
A = 2*(Lx*Ly+Ly*Lz+Lx*Lz)
P = 4*(Lx+Ly+Lz)

# Plates thickness
t = 2.0E-3;

# Create materials
steel = matC.IsoMat(E=210e9,nu=0.3,rho0=7800, eta = 0.0)
air   = matC.Fluid() 

# Create props
steel2mm = stPC.PlateProp(t,steel)

area_dof = dof.DOF(0,0,dof.DOFtype(typestr='area'))


# %% create models
f_c = steel2mm.coincidence_frequency(air.c0)/2/np.pi
print('coincidence frequency = {0}'.format(f_c))

# Create plate subsystems
plate1 = st2Dsys.RectangularPlate(1,Lx,Lz,prop = steel2mm)

#room     = ac3Dsys.Acoustic3DSystem(6, V , A, P, air)

# create semi infinite fluids
sif1 = jun.SemiInfiniteFluid((plate1,), air)

#om = 2*np.pi*np.logspace(2,np.log10(10000),3)




# define load
# room power
#power1Watt = lC.Load(omega, np.ones(omega.shape), dof.DOF(1,0,dof.DOFtype(typestr = 'power')), name = '1Watt')
# plate power
power1mWatt = lC.Load(omega, 0.001*np.ones(omega.shape), dof.DOF(1,3,dof.DOFtype(typestr = 'power')), name = '1mWatt')

#create SEA model
plate = mds.HybridModel((plate1,),xdata=omega)

#connect both
plate.add_SIF({'sif1' : sif1})

plate.add_load('1Watt',power1mWatt) # add 1mWatt per band 


#%% solving
plate.create_SEA_matrix(sym = 1)
plate.solve()

#%% plotting 1
plt.close(1)
plate.result.plot(1,ID = [1],xscale = 'log',yscale = 'log')
plt.figure(1)
plt.yscale('log')
#plt.ylim(1e-5,0.01)
plt.xticks(2*np.pi*fc,fclabels)
plt.xlabel('$f_c/$Hz')
#plt.ylabel('$TL/$dB')
plt.legend()
plt.show()


#%% more info
sif1_in = plate.power_input('sif1')

sif_all = sif1_in.sum()

#%% plot 3
plt.close(3)
sif1_in.plot(3,xscale='log',yscale='log')
#plt.plot(om_VA1,power_in_1_nonres,':',label = 'VA1 non-res')
#plt.yscale('log')

#plt.xscale('log')
plt.xlabel('$f_c/$Hz')
plt.ylabel('$\Pi_{in}/$W')

plt.xticks(2*np.pi*fc,fclabels)

plt.legend()
plt.tight_layout()
plt.show()


