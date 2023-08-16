# -*- coding: utf-8 -*-
"""
Module for dynamic material descriptions

The materialClasses package handles the material properties of fluids and
structures.

"""

import numpy as np
import scipy.integrate as integrate
import pyva.data.matrixClasses as mC
import pyva.data.dof as dof
import pyva.useful as uf

import matplotlib.pyplot as plt

def isscalar(data):
    return np.isscalar(data) and not(isinstance(data,str))

class Fluid:
    """
    The fluid class deals with fluid properties including damping
    
    This class implements some methods that may not seem reasonable
    for standart fluids, for example the frequency dependent methods
    for speed of sound and impedance.
    They are defined in such a way so that every daugther class can implement
    those methods especially when they are frequency dependent as for example 
    in case of the equivalent fluid model for fibre materials
    
    
    Attributes
    ----------
    c0: complex
        Speed of sound
    rho0: complex
        Density
    eta: float
        Damping loss
    dampingModell: str
        Identifier for damping model 
    """
    
    def __init__(self,c0=343.,rho0=1.23,eta=0.01,dynamic_viscosity=1.84e-5,kappa = 1.4, \
                      Cp=1005.1, heat_conductivity = 0.0257673, **kwargs):
        """
        Constructor for Fluid

        Parameters
        ----------
        c0 : complex, optional
            Speed of sound. The default is 343..
        rho0 : complex, optional
            density. The default is 1.23.
        eta : float, optional
            Damping loss. The default is 0.01.
        dynamic_viscosity : float, optional
            Dynamic viscosity. The default is 1.84e-5.
        kappa : float, optional
            ratio of specific heat capacities. The default is 1.4.
        Cp : float, optional
            Heat capacity at constant pressure. The default is 1.0051.
        heat_conductivity : float, optional
            Heat conductivity. The default is 0.0257673.
        **kwargs : dict
            Arbitrary keyword argument list.
        dampingModell: str
            Identifier for damping model 

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        None.
                
        Examples
        --------
            import materialClasses as mc
            myAir = mc.Fluid()
                        
        """
        
        #self.name = name
        self.c0 = c0;
        self.rho0 = rho0
        self.eta  = eta
        self.dynamic_viscosity = dynamic_viscosity
        self.kappa = kappa 
        self.Cp  = Cp
        self.heat_conductivity = heat_conductivity

        
        for kw in kwargs:
            if kw == 'damping_model':
                if kwargs[kw] in ['constant','viscous']:
                    self.damping_model = kwargs[kw]
                else:
                    raise ValueError('Unkown value {1} for keyword {0}'.format(kw,kwargs[kw]))
            else:
                raise ValueError('Unkown keyword {0}'.format(kw))

 
        
        
        
        
    def __str__(self):
        """
        str for fluids

        Returns
        -------
        _str : str
            Fluid description.

        """
        #_str  = "name            : {0}\n".format(self.name)
        _str  = "c0                : {0}\n".format(self.c0)
        _str += "rho0              : {0}\n".format(self.rho0)
        _str += "nu0               : {0}\n".format(self.nu0)
        _str += "eta               : {0}\n".format(self.eta)
        _str += "dynamic_visc      : {0}\n".format(self.dynamic_viscosity)
        _str += "Cp                : {0}\n".format(self.Cp)
        _str += "heat_conductivity : {0}\n".format(self.heat_conductivity)
        _str += "Pr                : {0}\n".format(self.Pr)
        _str += "kappa             : {0}\n".format(self.kappa)
        
        return _str

    def __repr__(self):
        _str = 'Fluid(c0={0},rho0={1},eta={2})'.format(self.c0, self.rho0,self.eta)
        return _str
    
    @staticmethod
    def air(temperature, pressure,h_rel=0.):
        """
        Determine precise properties of air from ambient conditions.
        
        Implementation is based on [1] that collects numerous papers dealing with various
        properties of air. Rasmussen condenses all papers into one formula collection.
        
        [1] Rasmussen, K. (1997). Calculation methods for the physical properties of air 
        used in the calibration of microphones.

        Parameters
        ----------
        temperature : float
            temperature in Kelvin
        pressure : float
            atmospheric pressure in bar
        h_rel : float
            relative humidity in percent. The default is 0 (dry air).

        Returns
        -------
        Fluid 
            air with properties according to environmental conditions.

        """
        T0 = 273.15 # temperature in Kelvin at 0C
        R_l = 287.05 # gas constant at dry air
        R_d = 461.   # gas constant of water steam
        M_d = 29     # molecular weight of dry air
        M_w = 18     # molecular weight of water
        p_sr = 101325 # The reference static pressure in Pa
        x_c = 0.0004 # Mole fraction of CO2 - recommended value under lab conditions
        
        # Table A.1 coefficients
        a_sv = [1.2378847E-5, -1.9121316E-02, 33.93711047,-6.3431645E+3]
        a_f  = [1.00062, 3.14E-8, 5.6E-7]
        a_Z  = [1.58123E-6, -2.9331E-8, 1.1043E-10, 5.707E-6, -2.051E-8,\
                1.9898E-4,  -2.376E-6,  1.83E-11 , -7.65E-9 ]
        a_c  = [331.5024, 0.603055, -5.28E-4, 51.471935,	0.1495874, \
                -7.82E-4, -1.82E-7,  3.73E-8, -2.93E-10,-85.20931, \
                -0.228525, 5.91E-5, -2.835149,-2.15E-13, 29.179762, 4.86E-4 ]
        a_k  = [ 1.400822, -1.75E-5, -1.73E-7, -0.0873629, -1.665E-4, \
                -3.26E-6,  2.047E-8, -1.26E-10, 5.939E-14, -0.1199717,\
                -8.693E-4, 1.979E-6,	 -1.104E-2, -3.478E-16, 4.50616E-2, 1.82E-6 ]
        a_e  = [84.986, 7.0, 113.157, -1., -3.7501E-3,-100.015 ]
        a_ka = [60.054, 1.846, 2.06E-6, 40.,	-1.775E-4]
        a_cp = [0.251625, -9.2525E-5, 2.1334E-7,	-1.0043E-10, 0.12477, \
                -2.283E-5, 1.267E-7, 0.01116, 4.61E-6, 1.74E-8 ]
        
        # unit changes
        T  = temperature
        T2 = T*T
        t = T-T0 # tmeperature in Celsius
        t2 = t*t # square
        p_s = pressure*100000
        p_s2 = p_s*p_s 
        H = h_rel
        kcal = 4186.8 # kilo calories in J 
        
        # speed of sound in dry air
        #c_d = 331.45*np.sqrt(1+TC/T0) # [1] Eq.(10.1)
        # Vapor pressure of water in Pa [2] eq. 22
        p_sv = np.exp(a_sv[0]*T2+a_sv[1]*T+a_sv[2]+a_sv[3]/T) # Buck Equation
        # Vapor pressure of water in Pa [2] eq. 22
        #C = 4.6151 - 6.8346*(T0+0.01/T)**1.261
        #p_sv = p_sr*10**C
        # Enhancement factor
        f = a_f[0]+a_f[1]*p_s+a_f[2]*t2
        # Mole fraction of water vapor in air
        x_w = H/100*p_sv/p_s*f
        # Compressibility factor
        ps_T = p_s/T
        x_w2 = x_w*x_w
        Z = 1-ps_T*(a_Z[0]+a_Z[1]*t+a_Z[2]*t2+(a_Z[3]+a_Z[4]*t)*x_w +\
                     (a_Z[5]+a_Z[6]*t)*x_w2) + ps_T*ps_T*(a_Z[7]+a_Z[8]*x_w2)
                
        # density of air
        rho = (3.48349+1.44*(x_c-0.0004))/1000*p_s/Z/T*(1-0.378*x_w)

        # Zero frequency speed of sound
        c_0 = a_c[0]+a_c[1]*t+a_c[2]*t2+(a_c[3]+a_c[4]*t+a_c[5]*t2)*x_w + \
              (a_c[6]+a_c[7]*t+a_c[8]*t2)*p_s +(a_c[9]+a_c[10]*t+a_c[11]*t2)*x_c + \
              a_c[12]*x_w2+a_c[13]*p_s2+a_c[14]*x_c*x_c+a_c[15]*x_w*p_s*x_c

        # Ratio of specific heats
        kappa = a_k[0]+a_k[1]*t+a_k[2]*t2+(a_k[3]+a_k[4]*t+a_k[5]*t2)*x_w + \
              (a_k[6]+a_k[7]*t+a_k[8]*t2)*p_s +(a_k[9]+a_k[10]*t+a_k[11]*t2)*x_c + \
              a_k[12]*x_w2+a_k[13]*p_s2+a_k[14]*x_c*x_c+a_k[15]*x_w*p_s*x_c
        # Viscosity of air
        eta = (a_e[0]+a_e[1]*T+(a_e[2]+a_e[3]*T)*x_w+a_e[4]*T2+a_e[5]*x_w2)*1.E-8
        # Thermal Conductivity
        k_a = (a_ka[0] + a_ka[1]*T +a_ka[2]*T2 +(a_ka[3] +a_ka[4]*T)*x_w)*1.E-8
        # Specific heat at constant pressure
        C_p = a_cp[0]+a_cp[1]*T+a_cp[2]*T2+a_cp[3]*T*T2+ \
              (a_cp[4]+a_cp[5]*T+a_cp[6]*T2)*x_w + \
              (a_cp[7]+a_cp[8]*T+a_cp[9]*T2)*x_w2
        
        return Fluid(c0=c_0,rho0=rho,eta=0,dynamic_viscosity=eta,\
                     Cp=C_p*kcal, heat_conductivity=k_a*kcal,kappa = kappa)

        
    @property    
    def z0(self):
        """
        Real characteristic impedance without damping

        Returns
        -------
        float
            Characteristic impedance.

        """
        return self.c0*self.rho0
    
    @property    
    def Pr(self):
        """
        Prandtl number.

        Returns
        -------
        float
            Prandtl number.

        """
        return self.Cp/self.heat_conductivity*self.dynamic_viscosity
    
    @property    
    def kinematic_viscosity(self):
        """
        Kinematic viscosity.

        Returns
        -------
        float
            Kinematic viscosity.

        """
        return self.nu0

    @property    
    def diffusivity(self):
        """
        Diffusivity.

        Returns
        -------
        float
            Diffusivity.

        """
        return self.heat_conductivity/self.Cp/self.rho0

    
    @property    
    def nu0(self):
        """
        Kinematic viscosity.
        

        Returns
        -------
        float
            Kinematic viscosity
            
        """
        return self.dynamic_viscosity/self.rho0


    def impedance(self,omega=0.):
        """
        Complex characteristic impedance including damping
        

        Parameters
        ----------
        omega : float or ndarray, optional
            angular frequency. The default is 0..

        Returns
        -------
        complex
            characteristic impedance.

        """
           
        return self.c_freq(omega)*self.rho_freq(omega)
        
    def wavenumber(self,omega):
        """
        Wavenumber including damping
        
        Parameters
        ----------
        omega : TYPE
            DESCRIPTION.

        Returns
        -------
        complex
            Complex wavenumber
            
        """
        return omega/self.c_freq(omega)

    def wavelength(self,omega):
        """
        wavelength
        
        Parameters
        ----------
        omega : TYPE
            DESCRIPTION.

        Returns
        -------
        float
            wavelength
            
        """
        return 2*np.pi/np.real(self.wavenumber(omega))

    def shear_wavenumber(self,omega):
        """
        Shear wavenumber according to Maa's theory
        
        Parameters
        ----------
        omega : float or ndarray
            angular frequency.

        Returns
        -------
        complex
            wavenumber
            
        """
        return np.sqrt(omega*self.rho0/self.dynamic_viscosity)
        
    def c_freq(self,omega = 0.):
        """
        Complex, frequency dependent speed of sound

        Parameters
        ----------
        omega : float
            Angular frequency.

        Returns
        -------
        complex
            speed of sound.

        """
        
        return self.c0/self.damping(omega)
        
    def rho_freq(self,omega):
        """
        density of fluid
        
        For simple fluid the output is the attribute rho0 but for later 
        implementation of fibre material this method is created to initiate the frequency
        dependence in the mother class
        
        Parameters
        ----------
        omega : float
            Angular frequency.         
        
        Returns
        -------
        complex
            density
            
        """

        if isscalar(omega):
            _rho = self.rho0
        else:
            _rho = np.tile(self.rho0,np.size(omega))
        return _rho


    def damping(self,omega):
        """
        Damping loss

        Parameters
        ----------
        omega : float or ndarray
            DESCRIPTION.

        Returns
        -------
        float or ndarray
            damping loss.

        """
        
        if isscalar(omega):
            eta = self.eta
        else:
            eta = np.ones(np.size(omega))*self.eta
        return 1-0.5j*eta
        
    def reflection_factor(self,omega,impedance,theta=0,area=1.):
        """
        Refection factor for interface to other fluids

        Parameters
        ----------
        omega : float
            angular frequency.
        impedance : complex
            impedance of interfacing fluid.
        theta : float, optional
            angle of incidence. The default is 0.
        area : float, optional
            ?????? . The default is 1..

        Returns
        -------
        complex
            reflection coefficient.

        """
        
        if callable(impedance):
            _z  = impedance(omega)
        else:
            _z  = impedance
        _z1 = self.impedance(omega)/area 
        return (_z*np.cos(theta)-_z1)/(_z*np.cos(theta)+_z1)

    def absorption(self,omega,impedance,theta=0):
        """
        Absorption coefficient for interface to other fluids

        Parameters
        ----------
        omega : float
            angular frequency.
        impedance : complex
            impedance of interfacing fluid.
        theta : float, optional
            angle of incidence. The default is 0.

        Returns
        -------
        complex
            Absorption coefficient.

        """
  
        return 1-np.abs(self.reflection_factor(omega,impedance,theta)**2)
        
    def absorption_diffuse(self,omega,z,theta_max=np.pi/2*0.99,theta_step = np.pi/200):
        """
        Diffuse absorption coefficient for interface to other fluids

        Parameters
        ----------
        omega : float
            angular frequency.
        z : complex
            impedance.
        theta_max : float, optional
            Maximum integration angle. The default is np.pi/2*0.99.
        theta_step : float, optional
            Angle integration step. The default is np.pi/200.

        Returns
        -------
        alpha : float
            diffuse field absorption coefficient.

        """
        
        theta_ = np.linspace(0,theta_max,int(np.floor((theta_max+theta_step)/theta_step)))
        denom = np.sin(theta_max)**2
        
        if uf.isscalar(omega):
                abs_theta = self.absorption(omega,z,theta_)
                #remove nans
                abs_theta[np.isnan(abs_theta)] = 0.
                alpha = 2*integrate.simps(abs_theta*np.sin(theta_)*np.cos(theta_), theta_)/denom
        else:
            alpha = np.zeros(omega.shape)
            for i_,om_ in enumerate(omega):
                abs_theta = self.absorption(om_,z,theta_)
                #remove nans
                abs_theta[np.isnan(abs_theta)] = 0.
                alpha[i_] = 2*integrate.simps(abs_theta*np.sin(theta_)*np.cos(theta_), theta_)/denom
        
        return alpha
    
    def acoustic_FE(self,omega,S,ID=[1],**kwargs):
        """
        Acoustic Finite Element radiator/end condition of plane wave fluid
        
        See also
        --------
        piston
        
        Parameters
        ----------
        omega : float
            angular frequency.
        S: float
            tube cross section
        ID : list of int, optional
            list of input ID. The default is [1].
        **kwargs : TYPE
            DESCRIPTION.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        DynamicMatrix
              0D aoustic radiation mobility S/Z * p = Q    Za = p/Q = p/(v*S) = Z/S     
        
        """
        
#        for kw in kwargs:
#            if kw == 'velocity':
#                if kwargs[kw]=='q':
#                    volumeflow = True
#                else:
#                    volumeflow = False
#            else:
#                raise ValueError('Unkown keyword {0}'.format(kw))
    
        # Check if omega is a xdata object
        if isinstance(omega, mC.DataAxis):
            if omega.type.type == 18:
                omega = omega.data*2*np.pi
            elif omega.type.type == 21:
                omega = omega.data
            else:
                raise ValueError('Omega must be an instance of DataAxis of type frequency')
                
                
        excdof = dof.DOF(ID,[0],dof.DOFtype(typestr=('pressure')) )
        resdof = dof.DOF(ID,[1],dof.DOFtype(typestr=('volume flow')) )
        xdata  = mC.DataAxis(omega,typestr='angular frequency')
        
        data   = np.zeros((1,len(omega)),dtype=complex)
        
        for iomega in range(len(omega)):
            data[0,iomega] = S/self.impedance(omega[iomega])
              
        return mC.DynamicMatrix(data,xdata,excdof,resdof,sym=1,shape=(1,1,len(omega)))
    
    
    def infinite_layer_TM(self,omega,wavenumber,thickness,ID=[1,2],**kwargs):
        """
        Calculated transfermatrix of infinite fluid layer 
        
        Deprecated: This method is part of the FluidLayer class

        Parameters
        ----------
        omega : float
            angular frequency.
        wavenumber : float
            wavenumbner in plane direction.
        thickness : float
            thickness of fluid layer.
        ID : list if int, optional
            list of input and output ID. The default is [1,2].
        **kwargs : TYPE
            DESCRIPTION.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        DynamicMatrix
            transfer matrix of infinite fluid layer.

        """
        
        # Check if wavenumber is an xdata object
        if isinstance(wavenumber, mC.DataAxis):
            if wavenumber.type.type == 26:
                wavenumber = wavenumber.data
                Nx = len(wavenumber)
            else:
                raise ValueError('Wavenumber must be an instance of DataAxis of type wavenumber')
        elif isinstance(wavenumber,np.ndarray):
            Nx = wavenumber.size
        elif isscalar(wavenumber):
            Nx = 1
        else:
            raise ValueError('Wavenumber must be a numpy vector or a scalar')
                
                
        
        k  = self.wavenumber(omega)
        z  = self.impedance(omega)
        kz = np.sqrt(k*k-wavenumber*wavenumber)
        
        tmm = np.zeros((2,2,Nx),dtype = np.complex128)
        tmm[0,0,:] = np.cos(kz*thickness)
        tmm[0,1,:] = 1j*z*np.sin(kz*thickness)
        tmm[1,0,:] = 1j*np.sin(kz*thickness)/z
        tmm[1,1,:] = tmm[0,0,:]

        
        dof1 = dof.DOF([ID[0],ID[0]],[0,3],('pressure','velocity'))
        dof2 = dof.DOF([ID[1],ID[1]],[0,3],('pressure','velocity'))

        xdata  = mC.DataAxis(np.array(wavenumber).flatten(),typestr='wavenumber')

        return mC.DynamicMatrix(tmm,xdata,dof2,dof1) # excdof are first argument     

        
class IsoMat:
    """
    Class for isotropic, solid materials
    
    Attributes
    ----------
        E: float
            young modulus
        rho0: float
            Density
        nu : float, optional
            Poisson numbeer.
        eta: float
            Damping loss
        dampingModell: str
            Identifier for damping model 
    """
    
    def __init__(self,E=7.1E10,rho0=2700.,nu=0.34,eta=0.01,**params):
        """
        Constructor of IsoMat
        
        The default meterial is aluminium

        Parameters
        ----------
        E : float, optional
            young modulus. The default is 7.1E10.
        rho0 : float, optional
            density. The default is 2700..
        nu : float, optional
            Poisson numbeer. The default is 0.34.
        eta : float, optional
            Damping loss. The default is 0.01.
        **params : dict
            arbitrary list of keyword arguments.

        Returns
        -------
        None.

        """
        
        #self.name = name
        self.E = E
        self.rho0 = rho0
        self.nu   = nu
        self.eta  = eta
        
    def __str__(self):
        #_str  = "name           : {0}\n".format(self.name)
        _str = "E              : {0}\n".format(self.E)
        _str += "rho0           : {0}\n".format(self.rho0)
        _str += "nu             : {0}\n".format(self.nu)
        _str += "eta            : {0}\n".format(self.eta)
            

        return _str

    def __repr__(self):
        return "IsoMat(E={0},rho0={1},nu={2},eta={3})".format(self.E,self.rho0,self.nu,self.eta) 
      
    @property    
    def c_T(self):
        """
        Complex shear wave speed
        
        This property includes the complex component from damping.

        Returns
        -------
        complex
            shear wave speed.

        """
        return np.sqrt(self.G_complex/self.rho0)

    @property    
    def c_S(self):
        """
        Complex shear wave speed
        
        This property includes the complex component from damping.

        Returns
        -------
        complex
            shear wave speed.

        """
        return np.sqrt(self.G_complex/self.rho0)
     
    @property
    def c_L(self):
        """
        Longitudinal wave speed
        
        This property includes the complex component from damping.

        Returns
        -------
        complex
            Longitudinal wave speed.

        """
        return np.sqrt(self.E_complex*(1-self.nu)/self.rho0/(1+self.nu)/(1-2*self.nu))

    @property
    def G(self):
        """
        Shear modulus

        Returns
        -------
        float
            shear modulus.

        """
        return self.E/2/(1+self.nu)
        
    @property
    def G_complex(self):
        """
        Complex shear modulus

        This property includes the complex component from damping.

        Returns
        -------
        complex
            complex shear modulus.

        """
        return self.G*(1+1j*self.eta)

    @property
    def E_complex(self):
        """
        Complex young modulus

        This property includes the complex component from damping.

        Returns
        -------
        complex
            complex young modulus.

        """
        return self.E*(1+1j*self.eta)
    
    def lambda_lame(self):
        """
        1st Laméconstant lambda.

        Returns
        -------
        float
            1st Lame constant.

        """
        G = self.G
        nu = self.nu
        return 2*G*nu/(1-2*nu)
    
    def wavenumber_L(self,omega):
        """
        Calculate longitudinal wavenumber in bulk material.

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        longitudinal wavenmuber.

        """
        return omega/self.c_L
    
    def wavenumber_S(self,omega):
        """
        Calculate shear wavenumber in bulk material.

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        shear wavenmuber.

        """
        return omega/self.c_S
    
    @property
    def bulk_modulus(self):
        """
        

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        
        return self.E_complex/(3-6*self.nu) 
    
class EquivalentFluid(Fluid):
    """
    The EquivalentFluid class deals with fibre models based on Champoux Allard with 
    limp/rigid frame theory
    
    This class is a doughter class of fluid.
    
    Attributes
    ----------
    flow_res: flow resistivity
    porosity: volume porosity
    tortuosity: tortuosity (alpha_inf)
    rho_bulk: apearant density of fluid and matrix
    length_visc: viscous characteristic length
    length_therm:  thermal characteristic length
    limp: Switch for use of limp model, false means rigid frame
    
    See [All2005]_ for details of the theory

    """
    
    def __init__(self,flow_res,porosity,tortuosity,rho_bulk,length_visc,length_therm,limp = True,\
                     c0=343.,rho0=1.23,eta=0.0,dynamic_viscosity=1.84e-5,kappa = 1.4, \
                     Cp=1005.1, heat_conductivity = 0.0257673):
        """
        Construcutor of equivalent fluid class
        
        The equivalent fluid model deals with the propagation of sound waves through
        a geometrical porous construction. It models the acoustic motion inside the 
        skeleon and the viscous and thermal interactions between fluid and skeleton.
        
    

        Parameters
        ----------
        flow_res : float
            flow resistivity
        porosity : float
            volume porosity.
        tortuosity : float
            tortuosity (alpha_inf)N.
        rho_bulk : float
            apearant density of fluid and matrix.
        length_visc : float
            viscous characteristic lengt.
        length_therm : float
            thermal characteristic length.
        limp : bool, optional
            Switch for use of limp model. The default is True.
        c0 : complex, optional
            Speed of sound. The default is 343..
        rho0 : complex, optional
            density. The default is 1.23.
        eta : float, optional
            Damping loss. The default is 0.01.
        dynamic_viscosity : float, optional
            Dynamic viscosity. The default is 1.84e-5.
        kappa : float, optional
            ratio of specific heat capacities. The default is 1.4.
        Pr : float, optional
            Prandtl number. The default is 0.71.

        Returns
        -------
        None.

        """
    
        super().__init__(c0,rho0,eta,dynamic_viscosity,kappa,Cp,heat_conductivity)
    
        self.flow_res   = flow_res
        self.porosity   = porosity
        self.tortuosity = tortuosity
        self.rho_bulk   = rho_bulk
        self.length_visc = length_visc
        self.length_therm = length_therm
        self.limp       = limp
        
    def __repr__(self):
        _str = 'EquivalentFluid(flow_res={0},porosity={1},tortuosity={2},rho_bulk={3},length_visc={4},length_visc={5},limp={6},'.\
            format(self.flow_res, self.porosity,self.tortuosity,self.rho_bulk,self.length_visc,self.length_visc,self.limp)
        _str += 'c0={0},rho0={1},eta={2})'.format(self.c0, self.rho0,self.eta)
        return _str
   
    def __str__(self):
        """
        str for EquivalentFluids

        Returns
        -------
        _str : str
            Fluid description.

        """
        _str  = super().__str__()
        _str += "flow_res          : {0}\n".format(self.flow_res)
        _str += "porosity          : {0}\n".format(self.porosity)
        _str += "tortuosity        : {0}\n".format(self.tortuosity)
        _str += "rho_bulk          : {0}\n".format(self.rho_bulk)
        _str += "lentgh_visc       : {0}\n".format(self.length_visc)
        _str += "lentgh_therm      : {0}\n".format(self.length_therm)
        
        return _str

    def G(self,omega):    
        """
        Characteristic function of the equivalent fluid
        
        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        complex
            Characteristic function.

        """
        tor2  = self.tortuosity*self.tortuosity
        rho0  = self.rho0
        
        return np.sqrt(1+4j*tor2*self.nu0*rho0*omega/ \
                    (self.flow_res*self.length_visc*self.porosity)/ \
                    (self.flow_res*self.length_visc*self.porosity));

    def rho_rigid(self,omega):
        """
        Dynamic complex density of the equivalent fluid with rigid matrix
        
        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        complex
            Dynamic complex density.

        """
        return self.tortuosity*self.rho0/self.porosity* \
               (1+self.flow_res*self.porosity*self.G(omega)/(1j*omega*self.rho0*self.tortuosity))

    def rho_freq(self,omega):
        """
        Complex density of the equivalent fluid
        
        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        complex
            Complex density.

        """
        rho02 = self.rho0**2
        rho   = self.rho_rigid(omega)

        if self.limp: # Limp correction
            return (self.rho_bulk*rho-rho02)/(rho+self.rho_bulk-2*self.rho0)
        else:
            return rho

    def bulk_modulus(self,omega):
        """
        Dynamic complex bulk modulus of the equivalent fluid
        
        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        complex
            Dynamic complex bulk modulus.
        """
        
        
        eta0 = self.dynamic_viscosity
        
        # /porosity added
        return self.c0*self.z0/self.porosity/(self.kappa-(self.kappa-1)/ \
            (1+8*eta0*np.sqrt(1+(1j*self.rho0*omega*self.Pr*self.length_therm**2)/(16*eta0))/\
            (1j*self.rho0*omega*self.Pr*self.length_therm**2)))

    def impedance(self,omega):
        """
        Complex characteristic impedance of the equivalent fluid

        Parameters
        ----------
        omega : float or ndarray
            angular frequency. 

        Returns
        -------
        complex
            characteristic impedance.

        """
         
        return np.sqrt(self.bulk_modulus(omega)*self.rho_freq(omega));

    def propagation_constant(self,omega):
        """
        propagation_constant

        Parameters
        ----------
        omega : float or ndarray
            angular frequency

        Returns
        -------
        complex
            propagation constant.
        """
         
        return 1j*omega*np.sqrt(self.rho_freq(omega)/self.bulk_modulus(omega))


    def c_freq(self,omega):
        """
        Complex, frequency dependent speed of sound of the equivalent fluid

        Parameters
        ----------
        omega : float
            Angular frequency.

        Returns
        -------
        complex
            speed of sound.

        """
         
        return 1j*omega/self.propagation_constant(omega)

class PoroElasticMat(EquivalentFluid):
    """
    The PoroElasticMat class deals porous and elastic material

    The material model is implemented according to Allard [All2009]_ and requires a geometry information 
    of the frame given by the EquivalentFlluid attributes and the bulk properties elastic frame matrial.
    In addition the bulk modulus of the material the frame is made of can be given, but most cases this can 
    be considered as much stiffer as the bulk material and is therefor considered as inf.
    
    Attributes
    ----------
    frame: flow resistivity
    porosity: volume porosity
    tortuosity: tortuosity (alpha_inf)
    rho_bulk: apearant density of fluid and matrix
    length_visc: viscous characteristic length
    length_therm:  thermal characteristic length
    E
    
    
    See [All2005]_ for details of the theory

    """
    
    def __init__(self,solid_mat,\
                      flow_res,porosity,tortuosity,length_visc,length_therm,limp = False,\
                      c0=343.,rho0=1.23,dynamic_viscosity=1.84e-5,kappa = 1.4, \
                      Cp=1005.1, heat_conductivity = 0.0257673 , Ks = np.Inf):
        """
        Construcutor of equivalent fluid class.
        
        The equivalent fluid model deals with the propagation of sound waves through
        a geometrical porous construction. It models the acoustic motion inside the 
        skeleon and the viscous and thermal interactions between fluid and skeleton.
        
        Parameters
        ----------
        solid_mat : IsoMat
            Bulk frame material (in vaccuum). 
        flow_res : float
            flow resistivity
        porosity : float
            volume porosity.
        tortuosity : float
            tortuosity (alpha_inf)N.
        rho_bulk : float
            apearant density of fluid and matrix.
        length_visc : float
            viscous characteristic lengt.
        length_therm : float
            thermal characteristic length.
        limp : bool, optional
            Switch for use of limp model. The default is True.
        c0 : complex, optional
            Speed of sound. The default is 343..
        rho0 : complex, optional
            density. The default is 1.23.
        eta : float, optional
            Damping loss. The default is 0.01.
        dynamic_viscosity : float, optional
            Dynamic viscosity. The default is 1.84e-5.
        kappa : float, optional
            ratio of specific heat capacities. The default is 1.4.
        Pr : float, optional
            Prandtl number. The default is 0.71.

        Returns
        -------
        None.

        
        """
        # Calculate rho_bulk - not used in Bio theory
        rho_bulk = solid_mat.rho0 + porosity * rho0
        super().__init__(flow_res,porosity,tortuosity,rho_bulk,length_visc,length_therm,limp,c0,rho0,0.,dynamic_viscosity,kappa,Cp,heat_conductivity) 
        self.solid_mat = solid_mat     
    
    def __repr__(self):
        """
        reps for PoroElasticMat

        Returns
        -------
        str_ : string
            PoroElasticMat expression.

        """
        _str = 'PoroElasticMat({0},flow_res={1},porosity={2},tortuosity={3},rho_bulk={4},length_visc={5},length_visc={6},'.\
            format(repr(self.solid_mat),self.flow_res, self.porosity,self.tortuosity,self.rho_bulk,self.length_visc,self.length_visc)
        _str += 'c0={0},rho0={1},eta={2})'.format(self.c0, self.rho0,self.eta)
        return _str
    
    def P(self,omega):
        """
        P function of poroelastic material 
        
        Eq. (6.28) of [All2009]_ assuming Ks >> Kb

        Parameters
        ----------
        omega : float
            Angular frequency.

        Returns
        -------
        complex
            P.

        """
        
        return 4/3*self.solid_mat.G_complex+self.solid_mat.bulk_modulus+(1-self.porosity)**2*self.bulk_modulus(omega) # /self.porosity current version is K/phi
    
    def Q(self,omega):
        """
        Q function of poroelastic material 
        
        Eq. (6.27) of [All2009]_ assuming Ks >> Kb

        Parameters
        ----------
        omega : float
            Angular frequency.

        Returns
        -------
        complex
            Q.

        """
        
        return (1-self.porosity)*self.bulk_modulus(omega)*self.porosity
    
    def R(self,omega):
        """
        R function of poroelastic material 
        
        Eq. (6.26) of [All2009]_ assuming Ks >> Kb

        Parameters
        ----------
        omega : float
            Angular frequency.

        Returns
        -------
        complex
            R.

        """
        
        return self.porosity**2*self.bulk_modulus(omega)

    def N(self,omega):
        """
        N function of poroelastic material, equal to shear modulus
        
        Just for usage of similar symbols

        Parameters
        ----------
        omega : float
            Angular frequency.

        Returns
        -------
        complex
            shear modulus.

        """
        
        return self.solid_mat.G_complex

    @property  
    def rho(self):
        """
        Joint static density of air and frame.
        
        Returns
        -------
        complex
            Dynamic coupled density.

        """
        
        return self.porosity*np.real(self.fluid.rho0)+self.solid_mat.rho0
        
    def rho12(self,omega):
        """
        Dynamic coupled density rho_12.
        
        Eq. (6.56)

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        complex
            Dynamic coupled density.

        """
        
        rho12_ = -self.porosity*self.rho0*(self.tortuosity-1)
        return rho12_+1j*self.flow_res*self.porosity**2*self.G(omega)/omega
    
    def rho11(self,omega):
        """
        Dynamic frame density rho_11.
        
        Eq. (6.56)

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        complex
            Dynamic frame density.

        """        
        return self.solid_mat.rho0 - self.rho12(omega)        
    
    def rho22(self,omega):
        """
        Dynamic frame density rho_11.
        
        Eq. (6.56)

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        complex
            Dynamic frame density.

        """        
        return self.porosity*self.rho0 - self.rho12(omega)
    
    def Delta(self,omega):
        """
        Delta helper function for Biot Theory
        
        Eq(6.69)

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        complex
            Delta.

        """
        
        P     = self.P(omega)
        Q     = self.Q(omega)
        R     = self.R(omega)
        rho11 = self.rho11(omega)
        rho12 = self.rho12(omega)
        rho22 = self.rho22(omega)
        
        return (P*rho22+R*rho11-2*Q*rho12)**2-4*(P*R-Q*Q)*(rho11*rho22-rho12)
    
    def wavenumbers(self,omega):
        """
        Provide wavenumber coefficients of poroelastic materials

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        delta1 : complex
            squared first compressional wavenmuber.
        delta2 : complex
            squared first compressional wavenmuber.
        delta3 : complex
            squared shear wavenmuber.
        mu1 : complex
            ratio of the velocity of air to the velocity of the frame of 1st compr. wave.
        mu2 : complex
            ratio of the velocity of air to the velocity of the frame of 2nd compr. wave.
        mu3 : complex
            ratio of the velocity of air to the velocity if the frame of shear wave.

        """
        

        P     = self.P(omega)
        Q     = self.Q(omega)
        R     = self.R(omega)
        rho11 = self.rho11(omega)
        rho12 = self.rho12(omega)
        rho22 = self.rho22(omega)
        rho_  = rho11*rho22-rho12*rho12

        omega2 = omega*omega
        
        PR_Q2 = P*R-Q*Q
        
        #ka    = omega*np.sqrt(rho22/R)
        
        PR2Q  = P*rho22+R*rho11-2*Q*rho12
        Delta = PR2Q*PR2Q-4*PR_Q2*rho_

        # delta are the squared values, because they are mainly used
        # compressional waves 
        delta1 = omega2/2/PR_Q2*(PR2Q-np.sqrt(Delta)) # Eq (6.67)
        delta2 = omega2/2/PR_Q2*(PR2Q+np.sqrt(Delta)) # Eq (6.68)
        mu1    = (P*delta1-omega2*rho11)/(omega2*rho12-Q*delta1) # Eq (6.71)
        mu2    = (P*delta2-omega2*rho11)/(omega2*rho12-Q*delta2)

        # shear wave
        delta3 = omega*omega/self.solid_mat.G_complex*rho_/rho22 # Eq (6.83)
        mu3    = -rho12/rho22 # Eq (6.84)
        
        return (delta1,delta2,delta3,mu1,mu2,mu3)
    
    def impedances(self,omega):
        """
        Characteristic impedances of poroelastic material.

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        Zf1 : complex
            Characterisitc impedance of the air for the first compressinal wave.
        Zf2 : complex
            Characterisitc impedance of the air for the second compressinal wave.
        Zs1 : complex
            Characterisitc impedance of the frame for the first compressinal wave.
        Zs2 : complex
            Characterisitc impedance of the frame for the second compressinal wave.


        """
        
        
        delta1_2,delta2_2,delta3_2,mu1,mu2,mu3 = self.wavenumbers(omega)
        delta1 = np.sqrt(delta1_2)
        delta2 = np.sqrt(delta2_2)
        
        R = self.R(omega)
        Q = self.Q(omega)
        P = self.P(omega)
        
        Zf1 = (R+Q/mu1)*delta1/self.porosity/omega # Eq (6.74)
        Zf2 = (R+Q/mu2)*delta2/self.porosity/omega # Eq (6.75)
        
        Zs1 = (P+Q*mu1)*delta1/omega # Eq (6.77)
        Zs2 = (P+Q*mu2)*delta2/omega
        
        return (Zf1,Zf2,Zs1,Zs2)
    
    def surface_impedances(self,omega,thickness):
        """
        Surface impedance of poroelastic layer with hard wall backing

        Parameters
        ----------
        omega : float
            angular frequency.
        thickness : float
            thickness of layer.

        Returns
        -------
        Z : ndarray
            surface impedance.

        """
        h = thickness
        phi = self.porosity
        
        delta1_2,delta2_2,delta3_2,mu1,mu2,mu3 = self.wavenumbers(omega)
        delta1 = np.sqrt(delta1_2)
        delta2 = np.sqrt(delta2_2)
        
        Zf1,Zf2,Zs1,Zs2 = self.impedances(omega)
        
        # Eq. (6.108)
        D = (1-phi+phi*mu2)*(Zs1-(1-phi)*Zf1*mu1)*np.tan(delta2*h) + \
            (1-phi+phi*mu1)*(Zf2*mu2*(1-phi)-Zs2)*np.tan(delta1*h) 
        # Eq. (6.107)
        Z = -1j*(Zs1*Zf2*mu2-Zs2*Zf1*mu1)/D
        
        return Z

    