# Example with 3 layups for sound absorption
import numpy as np
import matplotlib.pyplot as plt
import pyva.models as mds
import pyva.systems.infiniteLayers as iL
import pyva.properties.structuralPropertyClasses as stPC


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

steel = matC.IsoMat(E=15.E10,rho0=7800.,nu=0.27)
steel5mm = stPC.PlateProp(0.005, steel)
nothing  = stPC.PlateProp(0.000001, steel)
    
angles  = np.linspace(0,np.pi/2)
angle45 = np.pi/4 
kx45    = air.wavenumber(omega)*np.sin(angle45)
#kx     = mC.DataAxis(kx,typestr='wavenumber')

z0     = air.z0

il_fibre_10cm  = iL.FluidLayer(0.1,fibre1)
il_perforate   = iL.PerforatedLayer(0.005, 0.001, distance = 0.02)
il_steel_5mm   = iL.PlateLayer(steel5mm)
il_nothing_perf  = iL.PlateLayer(nothing,perforation = il_perforate)
il_steel_5mm_perf = iL.PlateLayer(steel5mm,perforation = il_perforate)
il_steel_5mm_s_perf = iL.ImperviousScreenLayer(steel5mm,perforation = il_perforate)

 



TMM_fibre_10         = mds.TMmodel((il_fibre_10cm,))
TMM_perf_fibre_10    = mds.TMmodel((il_perforate, il_fibre_10cm,))
TMM_steel_fibre_10   = mds.TMmodel((il_steel_5mm, il_fibre_10cm,))
TMM_steel_perf_fibre_10   = mds.TMmodel((il_steel_5mm_perf, il_fibre_10cm,))
TMM_steel_s_perf_fibre_10   = mds.TMmodel((il_steel_5mm_s_perf, il_fibre_10cm,))
TMM_nothing_perf_fibre_10   = mds.TMmodel((il_nothing_perf, il_fibre_10cm,))


alpha_steel_perf_fibre_10 = TMM_steel_perf_fibre_10.absorption_diffuse(omega,theta_max=np.pi/2,signal = False)
alpha_steel_s_perf_fibre_10 = TMM_steel_s_perf_fibre_10.absorption_diffuse(omega,theta_max=np.pi/2,signal = False,allard=True)
alpha_nothing_perf_fibre_10 = TMM_nothing_perf_fibre_10.absorption_diffuse(omega,theta_max=np.pi/2,signal = False)

alpha_fibre_10      = TMM_fibre_10.absorption_diffuse(omega,theta_max=np.pi/2,signal = False)
alpha_perf_fibre_10 = TMM_perf_fibre_10.absorption_diffuse(omega,theta_max=np.pi/2,signal = False)
alpha_steel_fibre_10 = TMM_steel_fibre_10.absorption_diffuse(omega,theta_max=np.pi/2,signal = False)




# %% plot1

# plot impedance of fibre results for publishing
plt.figure(1)
plt.plot(omega,alpha_fibre_10,label='10cm fibre')
plt.plot(omega,alpha_perf_fibre_10,label='perforate + 10cm fibre' )
plt.plot(omega,alpha_steel_fibre_10,label='5mm steel + 10cm fibre' )
plt.plot(omega,alpha_steel_perf_fibre_10,':',label='5mm perforated steel  + 10cm fibre' )
plt.plot(omega,alpha_nothing_perf_fibre_10,':',label='5mm perforated nothing :) + 10cm fibre' )
plt.plot(omega,alpha_steel_s_perf_fibre_10,'-.',label='5mm perforated solid steel  + ...' )

plt.xscale('log')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('absorption')
plt.legend()

plt.savefig('../source/images/TMM_perforate_test_abs_diffuse.png')

