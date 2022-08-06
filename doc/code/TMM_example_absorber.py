# Example with 3 layups for sound absorption
import numpy as np
import matplotlib.pyplot as plt
import pyva.models as mds
import pyva.systems.infiniteLayers as iL

# my packages
import pyva.properties.materialClasses as matC

plt.close('all')

omega = 2*np.pi*np.logspace(1,4,100)

# default fluid
air    = matC.Fluid(eta = 0.0)
fibre1 = matC.EquivalentFluid(porosity = 0.98, \
                               flow_res = 25000.,\
                               tortuosity = 1.02, \
                               length_visc = 90.e-6, \
                               length_therm = 180.e-6,\
                               rho_bulk = 0.98*1.20 + 30. , \
                               rho0 = 1.208, \
                               dynamic_viscosity = 1.81e-5 )
# Thickness
h1 = 0.1
h2 = 0.2

angles  = np.linspace(0,np.pi/2)
angle45 = np.pi/4 
kx45    = air.wavenumber(omega)*np.sin(angle45)
#kx     = mC.DataAxis(kx,typestr='wavenumber')

z0     = air.z0

fibre_10cm  = iL.FluidLayer(h1,fibre1)
fibre_20cm  = iL.FluidLayer(h2,fibre1)
perforate   = iL.PerforatedLayer(0.005, 0.001, distance = 0.02)

TMM_fibre_10      = mds.TMmodel((fibre_10cm,))
TMM_fibre_20      = mds.TMmodel((fibre_20cm,))
TMM_perf_fibre_20 = mds.TMmodel((perforate, fibre_10cm,))

alpha_fibre_10      = TMM_fibre_10.absorption_diffuse(omega,theta_max=np.pi/2,signal = False)
alpha_fibre_20      = TMM_fibre_20.absorption_diffuse(omega,theta_max=np.pi/2,signal = False)
alpha_perf_fibre_10 = TMM_perf_fibre_20.absorption_diffuse(omega,theta_max=np.pi/2,signal = False)


# %% plot1

# plot impedance of fibre results for publishing
plt.figure(1)
plt.plot(omega,alpha_fibre_20,label='20cm fibre')
plt.plot(omega,alpha_fibre_10,label='10cm fibre')
plt.plot(omega,alpha_perf_fibre_10,label='perforate + 10cm fibre' )

plt.xscale('log')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('absorption')
plt.legend()

plt.savefig('../source/images/TMM_absorber_abs_diffuse.png')

