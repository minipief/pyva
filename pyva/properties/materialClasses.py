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
    def air(temperature, pressure):
        """
        Determines precise properties from ambient conditions

        Parameters
        ----------
        temperature : float
            temperature in Kelvin
        pressure : float
            atmospheric pressure

        Returns
        -------
        Fluid 
            air with properties according to environmental conditions.

        """
        
        T0 = 273.15
        TC = temperature-T0
        rho_ = 1.293*T0/temperature*pressure/1.013
        eta_ = 1.e-5*(1.723+0.0047*TC)
        lambda_ = (0.02427+7.130e-5*TC)
        c_p_ = 1003.+0.1*TC
        c0_ = 331.5+0.6*TC
        
        return Fluid(c0=c0_,rho0=rho_,eta=0,dynamic_viscosity=eta_,Cp=c_p_, heat_conductivity=lambda_)
        
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
        Prandtl number

        Returns
        -------
        float
            Prandtl number.

        """
        return self.Cp/self.heat_conductivity*self.dynamic_viscosity
    
    @property    
    def nu0(self):
        """
        Kinematic viscosity 
        

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
        return "IsoMat(E={0},,rho0={1},nu={2},eta={3})".format(self.E,self.rho0,self.nu,self.eta) 
      
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
    
class EquivalentFluid(Fluid):
    """
    The EquivalentFluid class deals fibre models based on Champoux Allard with 
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
        flow_res : TYPE
            flow resistivity
        porosity : TYPE
            volume porosity.
        tortuosity : TYPE
            tortuosity (alpha_inf)N.
        rho_bulk : TYPE
            apearant density of fluid and matrix.
        length_visc : TYPE
            viscous characteristic lengt.
        length_therm : TYPE
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
        
        return np.sqrt(1+4*1j*tor2*self.nu0*rho0*omega/ \
                       ( self.flow_res*self.length_visc*self.porosity)/ \
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

    