# -*- coding: utf-8 -*-
#
# Example with 3 tubes arranged as expansion chamber
import numpy as np
import matplotlib.pyplot as plt

# my packages
import pyva.data.dof as dof
import pyva.systems.acoustic1Dsystems as ac1D
import pyva.systems.acousticRadiators as acR
 
import pyva.properties.materialClasses as matC
import pyva.data.matrixClasses as mC
import pyva.models as mds
import pyva.loads.loadCase as lC
import pyva.useful as uf


plt.close('all')
 
# Define frequency axis
deltaF = 5
f0     = 10
f1     = 4000

xdata  = mC.DataAxis(2*np.pi*np.arange(f0,f1/np.pi,deltaF),typestr='angular frequency')

# Tube parameter
R1 = 0.05
R2 = R1*np.sqrt(10)
R3 = 0.05

A1 = np.pi*R1**2
A2 = np.pi*R2**2
A3 = np.pi*R3**2

L1 = 0.2
L2 = 0.3
L3 = 0.2

# The fluid
air   = matC.Fluid(eta=0.01)

# radiation in to free field pie Z_a = rho c / S
Pow_free = np.real(0.5*air.impedance(0)/A3*1**2)
print('Power in free field pipe P_in = {0}'.format(Pow_free))
print('Radiation mobility Y_a = {0}'.format(A1/air.impedance(0)))

# S1/s2 ratios h
hs = (2,5,10)

# Plot testcase
om = xdata.data # 2*np.pi*np.arange(f0,f1/np.pi,deltaF)
k  = air.wavenumber(om)

# Define single tubes
tube1 = ac1D.AcousticTube(L1,air,A1)
tube2 = ac1D.AcousticTube(L2,air,A2)
tube3 = ac1D.AcousticTube(L3,air,A3)
end4  = acR.CircularPiston(R3,air)
# Reference tube of similar length
tube_ref = ac1D.AcousticTube(L1+L2+L3,air,A1)

# Create finite elements
elem1 = tube1.acoustic_FE(xdata,ID=[1,2])
elem2 = tube2.acoustic_FE(xdata,ID=[2,3])
elem3 = tube3.acoustic_FE(xdata,ID=[3,4])
rad4  =  end4.acoustic_FE(xdata,ID=[4])
rad4free  =  air.acoustic_FE(xdata,A3,ID=[4])
entry4free  =  air.acoustic_FE(xdata,A1,ID=[1])

elem_ref = tube_ref.acoustic_FE(xdata,ID=[1,4])

# Creare the 1D FE models

# Define required DOFtype
Qdof = dof.DOFtype(typestr=('volume flow'))
Pdof = dof.DOFtype(typestr=('pressure'))
# Nodes
NIDs   = [1,2,3,4]
# Response and excitation DOFs  
excdof = dof.DOF(NIDs,[0],Pdof,repetition = True )
resdof = dof.DOF(NIDs,[1],Qdof,repetition = True )

# Empty matrices
#FE2Dshape = (len(NIDs)*(len(NIDs)+1)//2,len(xdata))
#FEshape   = (len(NIDs),len(NIDs),len(xdata))

# Creat VAmodel
tube_network = mds.VAmodel(None,xdata, excdof, resdof, sym=1, dtype=complex) 
tube_ref_network = mds.VAmodel(None,xdata, excdof[0:4:3], resdof[0:4:3], sym=1, dtype=complex) 


# Add the elements
tube_network += elem1
tube_network += elem2
tube_network += elem3
# free of real pelase select
tube_network += rad4
tube_network += entry4free

# Same for reference
tube_ref_network += elem_ref
tube_ref_network += rad4
tube_ref_network += entry4free

# Define loadcases and solve
volume_source = lC.Load(xdata, 0.001*np.ones(len(xdata)), dof.DOF([1],[1],Qdof), name = 'VolumeFlow')
tube_network.add_load({1:volume_source})
tube_ref_network.add_load({1:volume_source})

tube_network.solve(loadresponse=True)
tube_ref_network.solve(loadresponse=True)

#%% fig2


#%% fig3


# Deterimine Transmission loss
print('Calculate in- and output power')

#%% fig10


pow_in  = tube_network.power(1,1)
pow_in.plot(10)
pow_out = tube_network.power(1,4,boundary = rad4) #free
pow_out.plot(10,cs='r')

plt.savefig('../source/images/power_expansion_chamber.png')

#%% fig12

pow_in  = tube_ref_network.power(1,1)
pow_in.plot(12)
pow_ref  = tube_ref_network.power(1,4,boundary = rad4)
pow_ref.plot(12,cs='r')

plt.savefig('../source/images/power_expansion_chamber_reference.png')

#%% Check transmission loss ( works with open end)
IL = pow_out.transfer(pow_ref,IDs=[4,4])
IL.plot(13,res = 'dB')

#plt.savefig('../source/images/power_expansion_chamber_IL.png')








