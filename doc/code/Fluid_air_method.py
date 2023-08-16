# Example with 3 layups for sound absorption
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pyva.models as mds
import pyva.systems.infiniteLayers as iL

# my packages
import pyva.properties.materialClasses as matC

plt.close('all')

P0 = 1.01325 # standart atm

Hs = np.linspace(0, 60,4)
Ps = (1.+np.linspace(-0.010,0.010,5))*P0
Ts = np.linspace(0,30,7)
T0 = 273.15


matC.Fluid.air(296.15, 1.012,25.)


#%%dry air
#%%speed of sound
plt.figure()
for i,P in enumerate(Ps):
    plt.plot(Ts,matC.Fluid.air(T0+Ts,P).c0,label = 'Ps={0:.4f} Bar'.format(P))
plt.xlabel('T/Celsius')
plt.ylabel('c0/(m/s)')
plt.legend()

#%%density
plt.figure()
for i,P in enumerate(Ps):
    plt.plot(Ts,matC.Fluid.air(T0+Ts,P).rho0,label = 'Ps={0:.4f} Bar'.format(P))
plt.xlabel('T/Celsius')
plt.ylabel('density/(kg/m^3)')
plt.legend()

#%%density
plt.figure()
for i,P in enumerate(Ps):
    plt.plot(Ts,matC.Fluid.air(T0+Ts,P).dynamic_viscosity,label = 'Ps={0:.4f} Bar'.format(P))
plt.xlabel('T/Celsius')
plt.ylabel('dyamic visc/(Pa s)')
plt.legend()

#%%pressure fixed
#%%speed of sound
plt.figure()
for i,H in enumerate(Hs):
    plt.plot(Ts,matC.Fluid.air(T0+Ts,P0,H).c0,label = 'h_rel={0:.0f}%'.format(H))
plt.xlabel('T/Celsius')
plt.ylabel('c0/(m/s)')
plt.legend()

#%%density
plt.figure()
for i,H in enumerate(Hs):
    plt.plot(Ts,matC.Fluid.air(T0+Ts,P0,H).rho0,label = 'h_rel={0:.0f}%'.format(H))
plt.xlabel('T/Celsius')
plt.ylabel('density/(kg/m^3)')
plt.legend()

#%%density
plt.figure()
for i,H in enumerate(Hs):
    plt.plot(Ts,matC.Fluid.air(T0+Ts,P0,H).dynamic_viscosity,label = 'h_rel={0:.0f}%'.format(H))
plt.xlabel('T/Celsius')
plt.ylabel('dyamic visc/(Pa s)')
plt.legend()



