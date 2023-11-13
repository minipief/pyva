# pyva double wall exampple
# using Allards transfer matrix method
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyva.models as mds
import pyva.systems.infiniteLayers as iL
import pyva.properties.structuralPropertyClasses as stPC
import pyva.properties.materialClasses as matC

plt.close('all')

# default fluid
air    = matC.Fluid(rho0=1.21,eta = 0.0)
z0     = air.z0

# Angles and frequncies
freq  = np.geomspace(100, 16000,100)
omega = np.pi*2*freq


# Isotropic materials
alu     = matC.IsoMat(eta = 0.1)
rubber  = matC.IsoMat(E=2.6e6,rho0=1200,nu=0.49,eta=0.00)

# Melamin foam paramters
melamin_vac = matC.IsoMat(E=3.0e5,rho0=12.0,nu=0.4,eta=0.1) # Frame in vaccuum    E6
melamin = matC.PoroElasticMat(melamin_vac, \
                            flow_res = 30000., \
                            porosity = 0.99, \
                            tortuosity = 1.01, \
                            length_visc = 250.E-6, length_therm = 550.E-6)


# plate properties
alu_1mm        = stPC.PlateProp(0.001, alu)
rubber_2mm     = stPC.PlateProp(0.002, rubber)
# test the foam as solid
foam_3cm_solid = stPC.PlateProp(0.03, melamin_vac) 
        
# Foam Layers
iL_foam_3cm = iL.PoroElasticLayer(melamin, 0.03)
iL_foam_3cm_solid = iL.SolidLayer(foam_3cm_solid)

# rubber and alu as solid- and screen layer
iL_rubber_solid_2mm = iL.SolidLayer(rubber_2mm)
iL_rubber_imper_2mm = iL.ImperviousScreenLayer(rubber_2mm)
iL_alu_solid_1mm    = iL.SolidLayer(alu_1mm)
iL_alu_imper_1mm    = iL.ImperviousScreenLayer(alu_1mm)

# Mass of Fluid as gap
# iL_nothing    = iL.MassLayer(1e-6,1.)
iL_nothing    = iL.FluidLayer(0.00001,air)

# TMmpodel using the solid or impervious screen formulation 
alu_melamin_rubber_solid = mds.TMmodel((iL_alu_solid_1mm,iL_foam_3cm,iL_rubber_solid_2mm))
alu_melamin_rubber_imper = mds.TMmodel((iL_alu_imper_1mm,iL_foam_3cm,iL_rubber_imper_2mm))
# Try other variations
alu_melamin_rubber_decoup = mds.TMmodel((iL_alu_imper_1mm,iL_nothing,iL_foam_3cm,iL_rubber_imper_2mm))
alu_melamin_rubber_foam_as_solid = mds.TMmodel((iL_alu_imper_1mm,iL_foam_3cm_solid,iL_rubber_imper_2mm))

tau_imper = alu_melamin_rubber_imper.transmission_diffuse(omega,theta_max=78/180*np.pi,allard=True,signal=False)
tau_solid = alu_melamin_rubber_solid.transmission_diffuse(omega,theta_max=78/180*np.pi,allard=True,signal=False)
tau_decoup = alu_melamin_rubber_decoup.transmission_diffuse(omega,theta_max=78/180*np.pi,allard=True,signal=False)
tau_foam_as_solid = alu_melamin_rubber_foam_as_solid.transmission_diffuse(omega,theta_max=78/180*np.pi,allard=True,signal=False)

    
#%% Plot1
plt.figure(1)
plt.plot(freq,-10*np.log10(tau_imper),label='impervious screen')
plt.plot(freq,-10*np.log10(tau_solid),label='solid')
plt.plot(freq,-10*np.log10(tau_decoup),label='decoup')
plt.plot(freq,-10*np.log10(tau_foam_as_solid),label='foam as solid')
plt.xlabel('f/Hz')
plt.ylabel('TL/dB')
plt.xscale("log")
plt.legend()

plt.savefig('../source/images/allard_DW_TL.png')


