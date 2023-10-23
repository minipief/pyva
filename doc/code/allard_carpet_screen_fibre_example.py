# -*- coding: utf-8 -*-
"""
Allard example of section 11.7.2 of [All2009]_ 
Example with lay-up given in figure 11.16 and table 11.7.
Table 11.7 is not correct - the fibrous layer thickness is 12.5mm noz 1.25mm.
Figure 11.17 shows the resuls but in Hz not in kHz as in [All2009]_

Alexander Peiffer
"""

import numpy as np
import matplotlib.pyplot as plt
import pyva.models as mds
import pyva.systems.infiniteLayers as iL
import pyva.properties.structuralPropertyClasses as stPC
import pyva.properties.materialClasses as matC

plt.close('all')

# Frequency Range
freq  = np.linspace(200,1400,100) # kHz
omega = 2*np.pi*freq

# Define Materials
air    = matC.Fluid(eta = 0.0)
# Carpet as Poroelastic
carpet_solid = matC.IsoMat(E=20000.,rho0=60,nu=0.,eta=0.5)
carpet_mat   = matC.PoroElasticMat(carpet_solid, \
                            flow_res = 5000., \
                            porosity = 0.99, \
                            tortuosity = 1., \
                            length_visc = 23.E-6, length_therm = 28.E-6)

# Impervious scree as PlateProp
screen_mat  = matC.IsoMat(E=30000.,rho0 = 2000, nu=0.49)
screen_prop = stPC.PlateProp(0.003, screen_mat)

# Fibre as Poroelastic
fibre_solid = matC.IsoMat(E=100000.,rho0=60,nu=0.,eta=0.88)
fibre_mat   = matC.PoroElasticMat(fibre_solid, \
                            flow_res = 33000., \
                            porosity = 0.98, \
                            tortuosity = 1.1, \
                            length_visc = 50.E-6, length_therm = 110.E-6)

# Air impedance for normalisation 
z0     = air.z0

# Define infiniteLayers
carpet1  = iL.PoroElasticLayer(carpet_mat, 0.0035)
carpet2  = iL.PoroElasticLayer(carpet_mat, 0.0035)
screen   = iL.ImperviousScreenLayer(screen_prop)
fibre    = iL.PoroElasticLayer(fibre_mat, 0.0125)

# Create lay-up as TMmodel
TMM_layup = mds.TMmodel((carpet1,carpet2,screen,fibre))

# Calculate normal surface impedance
Z          = TMM_layup.impedance_allard(omega,kx=0,signal = False)
# Calculate normal absorption
alpha0     = TMM_layup.absorption(omega,kx=0.0,signal = False,allard=True)
# Calculate diffure absorption
alpha_diff = TMM_layup.absorption_diffuse(omega,theta_max=np.pi*78/180,theta_step = np.pi/180,signal = False,allard=True)


#%% plot impedance of according to fig 11.17 of [All2009]_
# The frequency is in Hz not in kHz as given in fig 11.17
plt.figure(1)
plt.plot(freq,np.real(Z)/z0 ,label='Re')
plt.plot(freq,np.imag(Z)/z0 ,label='Im')
plt.xlabel('$f$/Hz')
plt.ylabel('$z/z0$')
plt.legend()

plt.savefig('../source/images/carpet_fibre_impedance.png')

#%% plot absorption
plt.figure(2)
plt.plot(freq,alpha0 ,label='normal')
plt.plot(freq,alpha_diff, label='diffuse')
plt.xlabel('$f$/Hz')
plt.ylabel( r'$\alpha$')
plt.legend()

plt.savefig('../source/images/carpet_fibre_absorption.png')


