# Example with 3 layups for sound absorption
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pyva.models as mds
import pyva.systems.infiniteLayers as iL

# my packages
import pyva.properties.materialClasses as matC

plt.close('all')

omega = 2*np.pi*np.logspace(2,4,100)
freq  = omega/2/np.pi

# default fluid
air    = matC.Fluid(eta = 0.0)
fibre1 = matC.EquivalentFluid(porosity = 0.999, \
                              flow_res = 380519.,\
                              tortuosity = 1.0, \
                              length_visc = 56.5e-6, \
                              length_therm = 113.0e-6,\
                              rho_bulk = 71.67 , \
                              rho0 = 1.23, \
                              dynamic_viscosity = 1.81e-5,
                              limp = True )

fibre2 = matC.EquivalentFluid(porosity = 0.999, \
                              flow_res = 438610.,\
                              tortuosity = 1.0, \
                              length_visc = 66.1e-6, \
                              length_therm = 66.1e-6,\
                              rho_bulk = 82.03 , \
                              rho0 = 1.23, \
                              dynamic_viscosity = 1.81e-5,
                              limp = True )

    
# Thickness
h = 0.012

# Create function from model
def fibre_fit(por,sigma,tor,Lam_visc,Lam_term):
#flow_res,porosity,tortuosity,rho_bulk    
    return matC.EquivalentFluid(sigma,por,\
                              tor, 71.67 ,\
                              Lam_visc, \
                              Lam_term,\
                              rho0 = 1.23, \
                              dynamic_viscosity = 1.81e-5,
                              limp = True )
        
def layer_fit(por,sigma,tor,Lam_visc,Lam_term):
    
    return mds.TMmodel((iL.FluidLayer(h,fibre_fit(por,sigma,tor,Lam_visc,Lam_term)),))


def alpha_fit(f,por,sigma,tor,Lam_visc,Lam_term):
    
    return layer_fit(por,sigma,tor,Lam_visc,Lam_term).absorption(2*np.pi*f,0.).ydata.flatten()





z0     = air.z0
fibre_12mm  = iL.FluidLayer(h,fibre1)
fibre_2_12mm  = iL.FluidLayer(h,fibre2)
TMM_fibre_12      = mds.TMmodel((fibre_12mm,))
TMM_fibre_2_12      = mds.TMmodel((fibre_2_12mm,))
alpha_fibre_12      = TMM_fibre_12.absorption(omega,0).ydata.flatten()
alpha_fibre_2_12      = TMM_fibre_2_12.absorption(omega,0).ydata.flatten()

# Import test data
f_test,alpha_test1,alpha_test2 = np.loadtxt ('.//data//impedance_test_results.csv',
                    unpack = True,
                    usecols = (0,1,2), skiprows = 1,
                    delimiter = ',')

# Import test data
f_test_Feb,alpha_test_Feb = np.loadtxt ('.//data//impedance_test_02_22.csv',
                    unpack = True,
                    usecols = (0,1), skiprows = 1,
                    delimiter = ',')

bounds = ([0.,20000,1.,1.e-6,1.e-6],[1.,600000,4.,5.e-4,5.e-4])

popt,pcov = curve_fit(alpha_fit, f_test_Feb[0:], alpha_test_Feb[0:],p0 = [0.999,438610., 1., 66.1e-6, 66.1e-6],bounds = bounds)



# %% plot1

# plot impedance of fibre results for publishing
plt.figure(1)
#alpha_fibre_12.plot(1)
#plt.plot(freq,alpha_fibre_12,label='12mm fibre')
plt.plot(freq,alpha_fibre_2_12,label='12mm fibre Feb')
#plt.plot(f_test,alpha_test1,'d:',label='test1')
plt.plot(f_test_Feb,alpha_test_Feb,'d:',label='test_02')

#plt.plot(f_test,alpha_test2,'o:',label='test2')

# check fit function f,por,sigma,Lam_visc,Lam_term
#plt.plot(freq,alpha_fit(freq, 0.999,438610.,  66.1e-6, 66.1e-6), label = 'fit_test') 
plt.plot(freq,alpha_fit(freq, *popt), label = 'fit_test') 


x0 = 4000
plt.text(x0, 0.6, r"$\sigma = {0:.2f}$".format(popt[1]), fontsize=12)
plt.text(x0, 0.55, r"$\Phi = {0:.2f}\mu$m".format(popt[0]), fontsize=12)
plt.text(x0, 0.5, r"$\alpha_\infty = {0:.2f}$".format(popt[2]), fontsize=12)
plt.text(x0, 0.45, r"$\Lambda = {0:.2f}\mu$m".format(popt[3]*1e6), fontsize=12)
plt.text(x0, 0.4, r"$\Lambda' = {0:.2f}\mu$m".format(popt[4]*1e6), fontsize=12)


#plt.xscale('log')

plt.xlim((250,6250))
plt.ylim((0,1))
plt.xlabel('$f/$Hz')
plt.ylabel('absorption')
plt.legend()

#plt.savefig('../source/images/TMM_absorber_abs_diffuse.png')

