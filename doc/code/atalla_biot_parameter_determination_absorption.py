# Example of JCA-parameter identification
#
# See Atalla Y, Panneton, R.; Inverse Acoustical Characterizatin of Open Cell Porous Media Using Impedance Tube Measurements
# 
# Taken from figure 7 - not 6

# Numerics
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import dual_annealing

# pyva packages
import pyva.properties.materialClasses as matC
import pyva.models as mds
import pyva.systems.infiniteLayers as iL


plt.close('all')

omega = 2*np.pi*np.geomspace(200,2000,100)
freq  = omega/2/np.pi

# default fluid
air    = matC.Fluid(eta = 0.0)

# Frequency range
f_min = 500

# select test
test = input('Enter test selector 1-4!')
test = int(test)

# switch if fibres are calculated based on limp model
fibre_sw = True # True # False means rigid

# name test
test_str = ('foam1','foam2','fibrous1','fibrous2')

# Given test data
if test == 1:
    # Foam 1
    flow_res = 4971.
    porosity = 0.97
    rho_bulk = 21.6
    h = 49.92e-3
    #h = 98.7e-3
    limp = False
    # reference values
    tor_ref = 1.31
    Lam_visc_ref = 123.19e-6
    Lam_term_ref = 289.54e-6
    # Environment
    T = 21.5+273.15
    P0 = 0.992
    
elif test == 2:
    # Foam 2
    flow_res = 8197.
    porosity = 0.95
    rho_bulk = 23.9
    h = 34.35e-3
    #h = 68.77e-3
    limp = False
    # reference values
    tor_ref = 1.42
    Lam_visc_ref = 133.10e-6
    Lam_term_ref = 212.65e-6
    # Environment
    T = 24.6+273.15
    P0 = 0.997

elif test == 3:
    # Fibrous 1
    flow_res = 21235.
    porosity = 0.94
    rho_bulk = 89.6
    h = 23.37e-3
    #h = 46.83e-3
    limp = fibre_sw
    # reference values
    tor_ref = 1.00
    Lam_visc_ref = 48.62e-6
    Lam_term_ref = 114.39e-6
        # Environment
    T = 22.6+273.15
    P0 = 1.002

elif test == 4:
    # Fibrous 2
    flow_res = 50470.
    porosity = 0.89
    rho_bulk = 150.0
    h = 37.38e-3
    #h = 18.98e-3
    limp = fibre_sw
    # reference values
    tor_ref = 1.00
    Lam_visc_ref = 41.51e-6
    Lam_term_ref = 114.27e-6
        # Environment
    T = 24.6+273.15
    P0 = 1.002


# Set air under test conditions
air_test = matC.Fluid.air(T,P0)


# Create function from model
def fibre_fit(tor,Lam_visc,Lam_term):
    """
    Function of fibre material with 3 parameters for optimisation

    Parameters
    ----------
    tor : float
        tortuosity.
    Lam_visc : float
        characteristic viscous length.
    Lam_term : float
        characteristic termal length.

    Returns
    -------
    EquivalentFluid object with 3 paramters
        
    """
    
    return matC.EquivalentFluid(flow_res, porosity , tor, rho_bulk, \
                              Lam_visc, Lam_term,\
                              rho0 = air_test.rho0, \
                              c0 = air_test.c0,\
                              dynamic_viscosity = air_test.dynamic_viscosity,\
                              Cp = air_test.Cp, heat_conductivity=air_test.heat_conductivity, \
                              limp = limp )
        
def layer_fit(tor,Lam_visc,Lam_term):
    """
    Function of fibre material layer with 3 parameters for optimisation.

    Parameters
    ----------
    tor : float
        tortuosity.
    Lam_visc : float
        characteristic viscous length.
    Lam_term : float
        characteristic termal length.

    Returns
    -------
    TMmodel object with 3 paramters

    """
    return mds.TMmodel((iL.FluidLayer(h,fibre_fit(tor,Lam_visc,Lam_term)),))


def impedance_fit(f,tor,Lam_visc,Lam_term):
    """
    Function of fibre material layer with 3 parameters for optimisation

    Parameters
    ----------
    f : ndarray of float
        frequnecy.
    tor : float
        tortuosity.
    Lam_visc : float
        characteristic viscous length.
    Lam_term : float
        characteristic termal length.

    Returns
    -------
    Signal
        surface impedance.

    """
    
    return layer_fit(tor,Lam_visc,Lam_term).impedance(2*np.pi*f,0.).ydata.flatten()





z0     = air_test.z0

# Import test data
f_test,Zs_re,Zs_im = np.loadtxt ('.//data//'+test_str[test-1]+'.csv',
                    unpack = True,
                    usecols = (0,1,2), skiprows = 1,
                    delimiter = ',')

# create absolute values
Zs = (Zs_re+1j*Zs_im)*z0
reflection = (Zs-z0)/(Zs+z0)
alpha_test = 1-abs(reflection**2)

# index for frequency range selection
i_freq = f_test >= f_min

# Create cost function with given test data
def cost_function(x):
    """
    Cost function of mean square residual 

    Parameters
    ----------
    x : array of 3 parameters 
        (tortuosity, visc_length, term_length).

    Returns
    -------
    float
        sum of squared difference between test and model impedance.

    """
    return np.sum(np.abs(impedance_fit(f_test[i_freq],x[0],x[1],x[2])-Zs[i_freq])**2)



# set bounds
lw = [1.,1.e-6,1.e-6] # lower bounds
up = [4.,4.e-4,4.e-4] # upper bounds
bounds=list(zip(lw, up))
# set start value
x0     = [1., 30.1e-6, 60.1e-6]



#res = least_squares(cost_function, x0,bounds = bounds)
res = dual_annealing(cost_function, bounds = bounds)


#%% plot1

# plot impedance of fibre results for publishing
plt.figure(1)
plt.plot(f_test,Zs_re,'d:',label='Re')
plt.plot(f_test,Zs_im,'d:',label='Im')

# check fit function f,por,sigma,Lam_visc,Lam_term
plt.plot(freq,np.real(impedance_fit(freq, *res.x)/z0), label = 'pyva_Re') 
plt.plot(freq,np.imag(impedance_fit(freq, *res.x)/z0), label = 'pyva_Im') 

# check Atalla solution ,por,sigma,Lam_visc,Lam_term
plt.plot(freq,np.real(impedance_fit(freq, tor_ref, Lam_visc_ref,Lam_term_ref)/z0), label = 'Atalla_Re') 
plt.plot(freq,np.imag(impedance_fit(freq, tor_ref, Lam_visc_ref,Lam_term_ref)/z0), label = 'Attala_Im') 

X0 = 550; X1 = 900; Y0 = -3.6; DY = -0.6
plt.text(X0, Y0, "Fit pyva:", fontsize=12)
plt.text(X1, Y0, r"$\alpha_\infty = {0:.2f}$".format(res.x[0]), fontsize=12)
plt.text(X1, Y0+DY, r"$\Lambda = {0:.2f}\mu$m".format(res.x[1]*1e6), fontsize=12)
plt.text(X1, Y0+2*DY, r"$\Lambda' = {0:.2f}\mu$m".format(res.x[2]*1e6), fontsize=12)
plt.text(X0, Y0+3*DY, "Atalla:", fontsize=12)
plt.text(X1, Y0+3*DY, r"$\alpha_\infty = {0:.2f}$".format(tor_ref), fontsize=12)
plt.text(X1, Y0+4*DY, r"$\Lambda = {0:.2f}\mu$m".format(Lam_visc_ref*1e6), fontsize=12)
plt.text(X1, Y0+5*DY, r"$\Lambda' = {0:.2f}\mu$m".format(Lam_term_ref*1e6), fontsize=12)


#plt.xscale('log')

plt.xlim((200,2000))
plt.ylim((-7,3))
plt.xlabel('$f/$Hz')
plt.ylabel('$Z_s$')

plt.title(test_str[test-1])
plt.legend()

plt.savefig('../source/images/atalla_JCA_parameter_'+test_str[test-1]+'.png')

