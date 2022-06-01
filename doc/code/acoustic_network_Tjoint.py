# -*- coding: utf-8 -*-
#
# Example with 2 tubes and Helmholts resonator
import numpy as np
import matplotlib.pyplot as plt
import copy

# my packages
import pyva.data.dof as dof
import pyva.systems.acoustic1Dsystems as ac1D
import pyva.systems.acousticRadiators as acR
 
import pyva.properties.materialClasses as matC
import pyva.data.matrixClasses as mC
import pyva.models as mds
import pyva.loads.loadCase as lC


plt.close('all')


# Define frequency axis
deltaF = 5
f0     = 10
f1     = 6000/2/np.pi
xdata  = mC.DataAxis(2*np.pi*np.arange(f0,f1,deltaF),typestr='angular frequency')

# Tube parameter
R1 = 0.01
A1 = np.pi*R1**2

L1 = 0.20
L3 = 0.20

# The fluid
air   = matC.Fluid(eta=0.0001)

# Perforate parameter
thickness = 0.0002 
holeR     = 0.0001
porosity  = 0.05

# Helmholtz parameter
V0        = 0.0002 
LH        = 0.02
R         = 0.01
Ac        = np.pi*R**2

# Reference acoustic impedance
Za0 = air.z0/Ac

# Define the perforate and HR
myPerf    = ac1D.PerforatedLayer(thickness,holeR,Ac,porosity = porosity)
myResPerf = ac1D.HelmholtzResonator(V0,LH,R,air,0.85,end_impedance=myPerf.radiation_impedance)

# Resonancefrequncies
om01         = myResPerf.omega_resonance()
print("HR resonance Frequency: {0:.1f}".format(om01/2/np.pi))

# Helmholtz
Za       = myResPerf.radiation_impedance(xdata.data)

#%% fig1
plt.figure(1)
plt.plot(xdata.data,np.real(Za),label='resistance',lw=2)
plt.plot(xdata.data,np.imag(Za),label='reactance',lw=2)
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('$Z_a/$Pa m$^2$')
plt.legend(loc = 0)

plt.tight_layout
plt.show
plt.savefig('../source/images/tjoint_HR_impeance.png')

#%% Model setup

# Define single tubes
tube1 = ac1D.AcousticTube(L1,air,A1)
tube3 = ac1D.AcousticTube(L3,air,A1)
end3  = acR.CircularPiston(R1,air)
tube_ref = ac1D.AcousticTube(L1+L3,air,A1)

# Create elements
elem1 = tube1.acoustic_FE(xdata,ID=[1,2])
elem3 = tube3.acoustic_FE(xdata,ID=[2,3])
rad3  =  end3.acoustic_FE(xdata,ID=[3])
entry4free  =  air.acoustic_FE(xdata,A1,ID=[1])

helmPerf = myResPerf.acoustic_FE(xdata,[2])

elem_ref = tube_ref.acoustic_FE(xdata,ID=[1,3])






# Creare the 1D FE models

# Define DOF the DOFs
Qdof = dof.DOFtype(typestr=('volume flow'))

# Nodes
NIDs   = [1,2,3]
# Response and excitation DOFs  
excdof = dof.DOF(NIDs,[0],dof.DOFtype(typestr=('pressure')),repetition = True )
resdof = dof.DOF(NIDs,[1],Qdof,repetition = True )


# Creat VAmodel
tube_network = mds.VAmodel(None,xdata,excdof,resdof,sym=1,dtype=complex) 

# Add the elements
tube_network += elem1
tube_network += elem3
# free of real pelase select
tube_network += rad3
tube_network += entry4free

tube_network += helmPerf

tube_ref  = mds.VAmodel(None,xdata,excdof[0:3:2],resdof[0:3:2],sym=1,dtype=complex) 
tube_ref  += entry4free
tube_ref  += elem_ref
tube_ref  += rad3


# Define loadcase
volume_source = lC.Load(xdata, 0.001*np.ones(len(xdata)), dof.DOF([1],[1],Qdof), name = 'VolumeFlow')

tube_network.add_load({1:volume_source}) # double pressure because of radiation into both directions
tube_network.solve(loadresponse=True)

tube_ref.add_load({1:volume_source})
tube_ref.solve(loadresponse=True)


#%% fig2
# Do the insertion loss
pow_out = tube_network.power(1,3,boundary = rad3) #free
pow_ref = tube_ref.power(1,3,boundary = rad3)

IL = pow_out.transfer(pow_ref,IDs=[3,3])
IL.plot(2,res='dB')

plt.savefig('../source/images/tjoint_IL.png')




