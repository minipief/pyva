# -*- coding: utf-8 -*-
"""
Module for acoustic one-dimensional systems

The acoustic1Dsystems module defines acoustic one dimensional systems and
related tools and methods to calculate their dynamic bahaviour.
Acoustic systems means 'fluid' systems.
"""

import pyva.properties.materialClasses as mc
import numpy as np
import pyva.data.matrixClasses as mC
import pyva.data.dof as dof
import pyva.useful as uF
import scipy.special as spl


class AcousticTube:
    """
    The acoustic tube class deals with one dimensional tubes filled with fluid
    
    Many system desciptions are implemented because of the examples in 
    [Pei2022]_ and therefor usually not required for system modelling but only
    for presentation purpose.
    
    Attributes
    ----------
    L:      length of the tube
    fluid:  fluid in the tube
    area:      cross section
    """
    
    def __init__(self,L,fluid=mc.Fluid,area=1.):
        """
        Class contructor for acoustic tube
        

        Parameters
        ----------
        L : float
            Length.
        fluid : fluid, optional
            fluid in the tube. The default is mc.Fluid.
        area : float, optional
            area of tube cross section. The default is 1..

        Returns
        -------
        None.
                
        Examples
        --------
        import acoustic1Dsystems as ac1Dsys
        myTube = ac1Dsys.AcousticTube(2.)
                        
        """
        
        
        self.L = L #< length of the tube
        self.fluid = fluid #< fluid in the tube
        self.area = area #< cross section of tube
        
    def __str__(self):
        """
        String output of acoustic tube

        Returns
        -------
        str

        """
        _str  = "AcousticTube: \n"
        _str += "L               : {0}\n".format(self.L)
        _str += str(self.fluid)
        _str += "area            : {0}\n".format(self.area)
        
        return _str
            
    def __repr__(self):
        """
        Reps of acoutic tube

        Returns
        -------
        None.

        """
        return "AcousticTube({0},fluid={1},S={2})".format(self.L,self.fluid.__repr__(),self.area)
        

        
    def p_simpleV1(self,omega,v):
        """
        Pressure at excited port 1 of the tube with vibrating surface
        
        Parameters
        ----------
        omega : float
            angular frequency.
        v : complex
            velocity at port 1.

        Returns
        -------
        complex
            pressure response.

        """
        
        c = self.fluid.cFreq(omega)
        z = self.fluid.impedance(omega) 
        return z*v/(1j*np.tan(omega/c*self.L))
    
    def p_amp(self,n,v):
        """
        Pressure amplitude for n-th resonance and velocity v
        
        Culculates the pressure at the resonance frequencies

        Parameters
        ----------
        n : int
            number of resonance.
        v : float
            velocit.

        Returns
        -------
        float
            pressure at resonance.

        """
        z0 = self.fluid.z0;
        return z0*v/np.tanh(n*np.pi*self.fluid.eta/2)

    def pressure(self,omega,x,v):
        """
        Pressure along the tube with vibrating surface at port 1
        
        take simple 1D transmission line model based on AcousticTube properties

        Parameters
        ----------
        omega : float
            angular frequency
        x : float
            position in tube
        v : complex
            velocity
                
                
        Returns
        -------
            pressure amplitude due to velocity         
        
        """
        
        z = self.fluid.impedance(omega)
        k = self.fluid.wavenumber(omega)
        L = self.L
        return z/(1j*np.sin(k*L))*(v[0]*np.cos(k*(L-x))-v[1]*np.cos(k*x))

    def velocity(self,omega,v,x):
        """
        Velocity along the tube with vibrating surface at port 1
        
        take simple 1D transmission line model based on AcousticTube properties

        Parameters
        ----------
        omega : float
            angular frequency
        x : float
            position in tube
        v : complex
            velocity
                
        Returns
        -------
            velocity amplitude due to velocity at port 1        
        
        """
        
        k = self.fluid.wavenumber(omega)
        L = self.L
        return (v[0]*np.sin(k*(L-x))+v[1]*np.sin(k*x))/(1j * np.sin(k*L))
        
    def p_mode(self,x,n):
        """
        Pressure mode shape for fixed ends
        
        Takes simple 1D transmission line model based on AcousticTube properties

        Parameters
        ---------
        x : position in tube
        n : number of mode
                
        Returns
        -------
            pressure shape of mode         
        
        """

        return np.where(n<1,np.sqrt(1./self.L)*np.cos(n*x), \
                            np.sqrt(2./self.L)*np.cos(n*np.pi/self.L*x))
            
    def k_mode(self,n):
        """
        Wavenumnber for mode n
        
        take simple 1D transmission line model based on AcousticTube properties

        Parameters
        ----------
        n : mode number
                
        Returns
        -------
        float
            wavenumber  
        """

        return n*np.pi/self.L
        
    def k2_mode(self,n):
        """
        Squared wavenumnber for mode n
         
        take simple 1D transmission line model based on AcousticTube properties

        Parameters
        ----------
        n : mode number
                
        Returns
        -------
        float
            squared wavenumber
        """
        
        return self.k_mode(n)*self.k_mode(n)    

    def f_mode(self,n):
        """
        Modal frequency

        Parameters
        ----------
        n : mode number
                
        Returns
        -------
        float
            modal frequency
        """
        return n*self.fluid.c0/(2*self.L)
        
    def omega_mode(self,n):
        """
        Modal angular frequency
        
        Parameters
        ----------
        n : mode number
                
        Returns
        -------
        float
            angular modal frequency
        """

        return n*self.fluid.c0*np.pi/self.L
        
    def p_mode_free(self,x,n):
        """
        Pressure mode shape for free ends
        
        Parameters
        ----------
        x : position in tube
        n : number of mode
                
        Returns
        -------
            pressure shape of mode         
        
        """

        if n<1:
            return np.sqrt(1./self.L)
        else:
            return np.sqrt(2./self.L)*np.sin(n*np.pi/self.L*x)
        
    def p_N(self,omega,q,x0,n):
        """
        Modal coordinates for acoustic source q at x0
        
        Takes simple 1D transmission line model based on AcousticTube properties

        Parameters
        ----------
        omega: float
            angular frequency
        q : complex
            volume source strengh, volume flow rate
        x : float
            position in tube
        n : int
            number of mode
                
                
        Returns
        -------
        modal coordinate for n         
        
        """
        
        rho0 = self.fluid.rho0
        c0   = self.fluid.c0
        eta  = self.fluid.eta
        norm = np.sqrt(2./self.L)

        #detect zero index
        norm=np.where(n<1,np.sqrt(1./self.L),np.sqrt(2./self.L))
        return 1j*norm*omega*rho0*self.p_mode(x0,n)*q/(self.k2_mode(n)*(1+1j*eta)-(omega/c0)**2)
            
        
    def p_modal(self,omega,x,N,v1 = 1.):
        """
        Modal pressure reponse for velocity at port 1
        
        Takes simple 1D transmission line model based on AcousticTube properties

        Parameters
        ----------
        omega : float
            angular frequency
        v1: complex
            velocity at port 1
        x: float
            position in tube
        N: int
            maximum number of mode
                                
        Returns
        -------
        pressure along tube         
        
        """
        
        # modal pressure with given velocity    
        p = self.p_N(omega,v1,0.,0)*self.p_mode(x,0)
        #print('Solution for: ',omega)
        for im in range(1,int(N+1)):
            p  += self.p_N(omega,v1,0.,im)*self.p_mode(x,im)
        return p
                   
    def power(self,omega,v):
        """
        Power input into tube at port 1
        
        Takes simple 1D transmission line model based on AcousticTube properties

        Parameters
        ----------
        omega: float
            angular frequency
        v : complex
            velocity at node 1
                
        Returns
        -------
            power input         
        
        """

        k0  = omega/self.fluid.c0
        z0  = self.fluid.z0
        eta = self.fluid.eta
        L   = self.L
        return 0.5*z0*np.abs(v*v)*(np.tanh(-k0*eta*L))/(np.cos(2*k0*L)/np.cosh(k0*eta*L)-1.)
            
    def acoustic_impedance(self,omega):
        """
        Acoustic impedance at port 1 with fixed end at port 2
        
        Parameters
        ----------
        omega: float
            angular frequency
                
        Returns
        -------
        complex
            impedance       
        
        """
        
        kL  = self.fluid.wavenumber(omega)*self.L
        z   = self.fluid.impedance(omega)
        return z/(1j*np.tan(kL))
        
    def transfer_impedance(self,omega,ID=[1,2],velocity = 'v'):
        """
        Transferimpedance of acoustic tubes
        
        Takes simple 1D transmission line model based on AcousticTube properties
        to model the so called transfer impedance, given by
        
        .. math:: 
            \\begin{Bmatrix} p_1 \\\\ v_1  \\end{Bmatrix} =
            \\begin{bmatrix} 
            T_{11} & T_{12} \\\\
            T_{21} & T_{22} 
            \\end{bmatrix} 
            \\begin{Bmatrix} p_2 \\\\ v_2 \\end{Bmatrix}
        
        or 
        
        .. math:: 
            \\begin{Bmatrix} p_1 \\\\ q_1  \\end{Bmatrix} =
            \\begin{bmatrix} 
            T_{a,11} & T_{a,12} \\\\
            T_{a,21} & T_{a,22} 
            \\end{bmatrix} 
            \\begin{Bmatrix} p_2 \\\\ q_2 \\end{Bmatrix}
        
        depending on the velocity paramter
        
        Parameters
        ----------
        omega : float
            angular frequency
        ID : list of int
            IDs of input and output port
        Kwargs:
            velocity = 'q':  
                
        Returns
        -------
        DynamicMatrix
            transfermatrix         
        
        """

        volumeflow = False
        
        if velocity == 'q':
            volumeflow = True
        elif velocity == 'v':
            volumeflow = False
        else:
            raise ValueError("velocity must be 'q' or 'v'")
        
        
        kz = lambda omega: self.fluid.wavenumber(omega)
        z  = lambda omega: self.fluid.impedance(omega)
        
        
        T11 = lambda omega: np.cos(kz(omega))
        if volumeflow:
            T12 = lambda omega: 1j*np.sin(kz(omega))*z(omega)/self.area
            T21 = lambda omega: 1j*np.sin(kz(omega))/z(omega)*self.area
            Tdof = dof.DOFtype(typestr=('pressure','volume flow'))
        else:
            T12 = lambda omega: 1j*np.sin(kz(omega))*z(omega)
            T21 = lambda omega: 1j*np.sin(kz(omega))/z(omega)
            Tdof = dof.DOFtype(typestr=('pressure','velocity'))
        
        excdof = dof.DOF([ID[1], ID[1]],[0,1],Tdof)
        resdof = dof.DOF([ID[0], ID[0]],[0,1],Tdof)
        xdata  = mC.DataAxis(omega,dof.DOFtype(typestr='angular frequency'))
        
        data   = np.zeros((2,2,len(omega)))
        for iomega in range(len(omega)):
            data[0,0,iomega] = T11(omega)
            data[1,1,iomega] = T11(omega)
            data[0,1,iomega] = T12(omega)
            data[1,0,iomega] = T21(omega)
              
        return mC.dynamicmatrix(mC.LinearMatrix(data),xdata,excdof,resdof)

    def acoustic_FE(self,omega,ID=[1,2]):
        """
        Acoustic Finite Element of acoustic tubes
        
        Takes simple 1D transmission line model based on AcousticTube properties
        to model the so mobility matrix, given by
        
        .. math:: 
            \\begin{bmatrix} 
            Y_{a,11} & Y_{a,12} \\\\
            Y_{a,21} & Y_{a,22} 
            \\end{bmatrix} 
            \\begin{Bmatrix} Q_1 \\\\ Q_2  \\end{Bmatrix} = 
            \\begin{Bmatrix} p_1 \\\\ p_2 \\end{Bmatrix}

        Parameters
        ----------
        omega : float
            angular frequency.
        ID : list of int, optional
            IDs of element ports. The default is [1,2].

        Returns
        -------
        DynamicMatrix
            Mobility matrix of tube         
        
        """
        

    
        # Check if omega is a xdata object
        if isinstance(omega, mC.DataAxis):
            if omega.type.type == 18:
                omega = omega.data*2*np.pi
            elif omega.type.type == 21:
                omega = omega.data
            else:
                raise ValueError('Omega must be an instance of DataAxis of type frequency')
                
        
        
        
        kz = lambda omega: self.fluid.wavenumber(omega)
        z  = lambda omega: self.fluid.impedance(omega)
                
        K11 = lambda omega: self.area/(1j*np.tan(kz(omega)*self.L)*z(omega) )
        K12 = lambda omega: self.area/(1j*np.sin(kz(omega)*self.L)*z(omega) )

        
        excdof = dof.DOF(ID,[0,0],dof.DOFtype(typestr=('pressure')) )
        resdof = dof.DOF(ID,[1,1],dof.DOFtype(typestr=('volume flow')) )
        xdata  = mC.DataAxis(omega,typestr='angular frequency')
        
        data   = np.zeros((3,len(omega)),dtype=complex)
        
        for iomega in range(len(omega)):
            data[0,iomega] = K11(omega[iomega])
#            data[1,iomega] = K11(omega[iomega])
            data[1,iomega] = K12(omega[iomega])
            data[2,iomega] = K11(omega[iomega])
              
        return mC.DynamicMatrix(data,xdata,excdof,resdof,sym=1,shape=(2,2,len(omega)))
    
class LumpedAcoustic:
    """
    The LumpedAcoustic class represents the lumped elements with no wave motion along the duct
    
    The aim of this class is mainly to be a mother class to all following limp elements.
    Thus, there are only two attributes with the key attribute imnpedance that defines the dynamics and 
    is a function.

    
    Attributes
    ----------
    impedance: constant or  frequency depenent impedance of the layer in the 1D element

    """
    
    def __init__(self,impedance,area=1):
        """
        Class contructor for LumpedAcoustic objects        

        Parameters
        ----------
        impedance : function 
            transfer impedance function with omega as argument.
        area : float, optional
            cross section area. The default is 1.

        Returns
        -------
        None.

        """
        
        self.impedance = impedance
        self.area      = area
        
    def plot(self,omega,fig=1,**kwargs):
        """
        Plots the transfer impedance for the related part

        Parameters
        ----------
        omega : float
            angular frequency.
        fig : int, optional
            figure ID. The default is 1.
        **kwargs : dict
            Arbitrary keyword parameter list passed to Signal.plot.

        Returns
        -------
        None.

        """
        
        _Z = self.Signal(omega)
        _Z.plot(fig, **kwargs)
 

    def Signal(self,omega):
        """
        generates a Signal from the requested frequency range

        Parameters
        ----------
        omega : float
            angular frequency.
            
        Returns
        -------
        Signal
            slope of transfer impedance.

        """
        
        _omega_val       = uF.get_omega_values(omega) 
        _omega_DataAxis = uF.get_omega_DataAxis(omega)
    
        return mC.Signal(_omega_DataAxis,self.impedance(_omega_val),dof = dof.DOF([0],[0],dof.DOFtype(typestr='pressure',xtypestr='velocity')))        
    
    def radiation_impedance(self,omega):
        return self.impedance(omega)/self.area

        
    def transfer_impedance(self,omega,ID=[1,2],**kwargs):
        """
        Transferimpedance of LumpedAcoustic
        
        Takes simple 1D transmission line model based on AcousticTube properties
        to model the so called transfer impedance, given by :math:`Z = (p_2-p_1)/v` 
        
        This class is an abstrat class that has no implementation of the specific
        transfer impedance. This must be done  by the daughter classes 
        
        
        .. math:: 
            \\begin{Bmatrix} p_1 \\\\ v_1  \\end{Bmatrix} =
            \\begin{bmatrix} 
            1 & Z \\\\
            0 & 1 
            \\end{bmatrix} 
            \\begin{Bmatrix} p_2 \\\\ v_2 \\end{Bmatrix}
        
        

        Parameters
        ----------
        omega : float
            angular frequency
        ID: list of int
            IDs of the two ports
                                
        Returns
        -------
        DynamicMatrix
            [2 x 2] array of transferimpedance         
        
        """
        

        # init
        lDOF = [0,1] # default pressure and velocity DOF 
        
        for kw in kwargs:
            if kw == 'velocity':
                if kwargs[kw]=='q':
                    volumeflow = True
                else:
                    volumeflow = False
            elif kw == 'DOF':
                lDOF = kwargs[kw] 
            else:
                raise ValueError('Unkown keyword {0}'.format(kw))
        
        S  = self.area
        omega = uF.get_omega_values(omega)
        
        if volumeflow:
            T12 = lambda omega: self.impedance(omega)/S
            Tdof = ('pressure','volume flow')
        else:
            T12 = lambda omega: self.impedance(omega)
            Tdof = ('pressure','velocity')
        
        excdof = dof.DOF([ID[1],ID[1]],lDOF,Tdof)
        resdof = dof.DOF([ID[0],ID[0]],lDOF,Tdof)
        xdata  = mC.DataAxis(omega,typestr='angular frequency')
        
        data   = np.zeros((2,2,len(omega)),dtype = np.complex128)
        for iomega in range(len(omega)):
            data[0,0,iomega] = 1
            data[1,1,iomega] = 1
            data[0,1,iomega] = T12(omega[iomega])
              
        return mC.DynamicMatrix(mC.LinearMatrix(data),xdata,excdof,resdof)

    def acoustic_FE(self,omega,ID=[1,2],**kwargs):
        """
        Acoustic Finite Element of MassLayers
        
        Takes simple 1D transmission line model based on AcousticTube properties
        to model the so mobility matrix, given by
        
        .. math:: 
            \\frac{S}{Z(\\omega)} 
            \\begin{bmatrix} 
            1 & 1 \\\\
            1 & 1 
            \\end{bmatrix} 
            \\begin{Bmatrix} p_1 \\\\ p_2  \\end{Bmatrix} = 
            \\begin{Bmatrix} Q_1 \\\\ Q_2 \\end{Bmatrix}

        Parameters
        ----------
        omega : float
            angular frequency.
        ID : list of int, optional
            IDs of element ports. The default is [1,2].

        Returns
        -------
        DynamicMatrix
            Mobility matrix of limped acoustic         
        
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
        omega = uF.get_omega_values(omega)
                                
        fac = lambda omega: self.area/self.impedance(omega)

        excdof = dof.DOF(ID,[0,0],dof.DOFtype(typestr=('pressure')) )
        resdof = dof.DOF(ID,[1,1],dof.DOFtype(typestr=('volume flow')) )
        xdata  = mC.DataAxis(omega,typestr='angular frequency')
        
        data   = np.ones((3,len(omega)),dtype=complex)
        
        for iomega in range(len(omega)):
            data[0,iomega] = fac(omega[iomega])
            data[1,iomega] = fac(omega[iomega])
            data[2,iomega] = fac(omega[iomega])
              
        return mC.DynamicMatrix(data,xdata,excdof,resdof,sym=1,shape=(2,2,len(omega)))


class MassStiffness(LumpedAcoustic):
    """
    The MassStiffness class provides a generic mass, spring and damping element in an acoustic 1D network
    
    Attributes
    ----------
    mass : of the layer in the 1D element
    stiffness : of the layer in the 1D element
    area : the mass/stiffness is extended over 
    vdamping : vicous damping
    damping : damping loss
    """
    
    def __init__(self,mass,stiffness,area=1.,damping=0,vdamping=0):
        """
        Class contructor for MassStiffness objects
        
        Parameters
        ----------
        mass : float
            of the layer in the 1D element
        stiffness : float
            of the layer in the 1D element
        area : float
            area of the mass/stiffness is extended over 
        vdamping: float
            vicous damping
        damping: float
            damping loss of the tube
                
        Examples
        --------
        import acoustic1Dsystems as ac1Dsys
        myTube = ac1Dsys.1Dlumped(1,0)
                        
        """
        
        self.mass = mass
        self.stiffness = stiffness
        self.damping = damping
        self.vdamping = vdamping

        # Implementation of impedance function of mother class
        impedance = lambda omega: 1j*omega*mass/area + stiffness*(1+1j*damping)/(1j*omega*area) + vdamping/area 
        
        super().__init__(impedance,area)
        
    def __str__(self):
        """
        str for MassStiffness

        Returns
        -------
        str

        """
        _str  = "MassStiffness: \n"
        _str += "m              : {0}\n".format(self.mass)
        _str += "k_s            : {0}\n".format(self.stiffness)
        _str += "c_v            : {0}\n".format(self.vdamping)
        _str += "eta            : {0}\n".format(self.damping)
        _str += "area           : {0}\n".format(self.area)
        
        return _str
        
        
        return "AcousticTube of length {0}, cross section {1}".format(self.L,self.area)
    
    def __repr__(self):
        """
        Reps of acoutic tube

        Returns
        -------
        None.

        """
        return str(self)

    @property        
    def mass_per_area(self):
        """
        Mass per area 

        Returns
        -------
        float
            specific mass per area.

        """
        
        return self.mass/self.area

    @property        
    def stiffness_per_area(self):
        """
        Stiffness per area

        Returns
        -------
        float 
            ftiffness per area.

        """
        
        return self.stiffness/self.area
        
class Membrane(MassStiffness):
    """
    The Membrane class deals with membranes
    
    The behaviour is calculated assuming a circular membrane with similar area
    
    Attributes
    ----------
    thickness : of the membrane
    material : of the membrane
    tension : of the membrane
    """
    
    def __init__(self,material,thickness,tension,area,**kwargs):
        """
        Class contructor for MassLayer
        
        Parameters
        ----------
        thickness : float
            thickness of the membrane 
        material : Fluid or isoMat
            material of the membrane 
        tension : float
            tension of the membrane 
        area : float
            area of the membrane
                
        Examples
        --------
        import acoustic1Dsystems as ac1Dsys
        my_mem = ac1Dsys.membrane(rubber,0.001,200,0.01)
                        
        """
        

        mass = 4/3*material.rho0*thickness*area
        stiffness = 8*tension
        
        super().__init__(mass,stiffness,area,**kwargs)
        
        self.thickness = thickness
        self.material = material        
        self.tension = tension
        
        
        
class PerforatedLayer(LumpedAcoustic):
    """
    The PerforatedLayer class deals perforated layers in acoustic networks
    
    The behaviour is calculated using flow throught the holes and radiation be the disk radiator pattern
    
    Attributes
    ----------
    thickness : of the perforted layer
    pattern : of the holes (quadradic, triangular, )
    hole_radius : of the hole sin the membrane
    porosity : area ration of holes and closed surface
    fluid : fluid in the holes
    alpha: correction constant for resistivity correction approx 4.-2. for
           sharp to round edges
    """
    

    
    
    def __init__(self,thickness,hole_radius,area=1.,fluid=mc.Fluid(damping_model='viscous'),pattern='square',alpha=2.,**kwargs):
        """
        Class contructor for PerforatedLayer
        
        Parameters
        ----------
        thickness : float
            thickness of perforate
        hole_radius : float
            radius of perforate holes
        area : float
            area of perforate plate
        fluid : fluid
            fluid in the hole
        pattern : str
            'square','triangular','hexagonal' identifier for the hole pttern
        alpha : float
            factor for edge sharpness 
        kwargs : dict
            Arbitrary keyword parameter list
        distance : float
            distance between holes
        porosity: float
            surface porosity of perforate
                
        Examples
        --------
        import acoustic1Dsystems as ac1Dsys
        myTube = ac1Dsys.material(rubber,0.001)
                        
        """
        
        self.thickness   = thickness
        self.hole_radius = hole_radius        
        self.fluid       = fluid
        self.alpha       = alpha
        
        if pattern in ('square','triangular','hexagonal'):
            self.pattern = pattern
        else:
            raise ValueError('Unknown value for pattern argument!')
            
        for kw in kwargs:
            if kw == 'distance':
                self.distance = kwargs[kw]
                if pattern == 'square':
                    self.porosity = np.pi*self.hole_radius**2/self.distance**2
                elif pattern == 'triangular':
                    self.porosity = 2*np.pi*self.hole_radius**2/np.sqrt(3)/self.distance**2
                elif pattern == 'hexagonal':
                    self.porosity = 4*np.pi*self.hole_radius**2/3/np.sqrt(3)/self.distance**2
            elif kw == 'porosity':
                self.porosity = kwargs[kw]
                if pattern == 'square':
                    self.distance = np.sqrt(np.pi/self.porosity)*self.hole_radius
                elif pattern == 'triangular':
                    self.distance = np.sqrt(2*np.pi/np.sqrt(3)/self.porosity)*self.hole_radius
                elif pattern == 'hexagonal':
                    self.distance = np.sqrt(4*np.pi/3/np.sqrt(3)/self.porosity)*self.hole_radius
                
            else:
                raise ValueError('Unkown keyword {0}'.format(kw))
        
        # Define impedance funtion for attribue
        def impedance_(omega):  
            _val = self.transfer_resistivity(omega)+ \
                   1j*self.transfer_reactance(omega)
        
            return _val
        
        # assign attributes of mother classes
        super().__init__(impedance_,area)
        
        
    def transfer_resistivity(self,omega):
        """
        Calculates the scalar transfer impedance of the flow in one perforate
        
        Parameters
        ----------
        omega : float
            angular frequency
                
        Returns
        -------
        float
            transfer resistivity
                        
        """        
        
        ks  = self.fluid.shear_wavenumber(omega)
        _R  = self.hole_radius
        
        _f1 = 8*self.thickness*self.fluid.dynamic_viscosity/_R**2/self.porosity
        _f2 = np.sqrt(1+((ks*_R)**2)/8) + (self.alpha*np.sqrt(2)*ks*_R**2)/8/self.thickness
                
        return _f1*_f2
    
    def transfer_reactance(self,omega):
        """
        Calculates the resistive end correction of both siges of one perforate
        
        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        float
            transfer reactance.

                                
        """        

        ks  = self.fluid.shear_wavenumber(omega)
        _R  = self.hole_radius

        _f1 = omega*self.thickness*self.fluid.rho0/self.porosity
        _f2 = 1 + 1/np.sqrt(9 + ((ks*_R)**2)/2) + self.reactive_end_correction_perforate(omega)
        
     
        return _f1*_f2
        
    def reactive_end_correction_perforate(self,omega):
        """
        Calculates the reactive end correctoin of both sides of one perforate
        
        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        float
            reactive end correction.
        """        
        
        if self.pattern == 'square':
            _eps = 0.47
        elif self.pattern == 'triangular':
            _eps = 0.31
        else: # self.pattern == 'hexagonal':
            _eps = 0.62    
        
        _fint = (1+_eps)*np.sqrt(self.porosity)-_eps*np.sqrt(self.porosity**3)
        _Xcorr = 16/3/np.pi*self.hole_radius/self.thickness*(1-_fint)
        
        return _Xcorr
    
class HelmholtzResonator:
    """
    The HelmholtzResonator class provides methods for modelling of HR
        
    Attributes
    ----------
    volume : volume of the resonace volume
    length : length of the neck
    radius : radius of neck 
    length_correction : sum of upper and lower correction part
    Fluid : fluid in the system
    end_impedance : function handle for end transfer impedance
    """
        
    
    def __init__(self,volume,length,radius,fluid = mc.Fluid('air'),length_correction_factor = 1.7,**kwargs):
        """
        Class contructor for helholtzResonator
        
        Parameters
        ----------
        volume : float
            volume of the resonace volume
        length : float
            length of the neck
        radius : float
            radius of neck 
        length_correction_factor: float
            sum of upper and lower end correction as factor of radius
        fluid: fluid
            fluid in the system
        KwArgs: dict
            Arbitrary keyword parameter list
        end_impedance: function
             radiation impedance of end correction
                                        
        """
        
        self.volume      = volume
        self.radius      = radius     
        self.length      = length
        self.fluid       = fluid
        self.length_correction = length_correction_factor*radius
        self.end_impedance   = lambda omega: 0*omega
                    
        for kw in kwargs:
            if kw == 'end_impedance':
                self.end_impedance = kwargs[kw]
            else:
                raise ValueError('Unkown keyword {0}'.format(kw))
                
    @property
    def area(self):
        """
        Area property

        Returns
        -------
        float
            area.

        """
        
        return np.pi*self.radius**2
                
    def radiation_impedance(self,omega):
        """
        Radiation impedance of a HelmholtzResonator Za = p/Q
        
        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        Z : complex
            radiation impedance.

        """
        
        _L = self.length+self.length_correction

        Z = self.end_impedance(omega) + \
            1j*(omega*self.fluid.rho0*_L/self.area - self.fluid.rho0*self.fluid.c0**2/(omega*self.volume)) 

        return Z
        
    def omega_resonance(self):
        """
        Angular resonance frequency of HR
        
        Assumes the resistance free resonance frequency
                
        Returns
        -------
        float
            resonance angular frequency
        """
        
        return self.fluid.c0*np.sqrt(self.area/(self.length+self.length_correction)/self.volume)

    def acoustic_FE(self,omega,ID=[1],**kwargs):
        """
        Acoustic Point Element of HR
        
        A p = q with A = 1/Z_rad
        
        Uses the HR radiation impedance to determine the properties of this point element.
             
        Parameters
        ----------
        omega : float
            angular frequency.
        ID : list of int
            node ID.

        Returns
        -------
        Dynamic matrix
            single DOF FE representation of HR
        """
           
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
            data[0,iomega] = 1./self.radiation_impedance(omega[iomega])
              
        return mC.DynamicMatrix(data,xdata,excdof,resdof,sym=1,shape=(1,1,len(omega)))
    
class QuarterWaveResonator:
    """
    The QuarterWaveResonator class provides methods for modelling of this
        
    Attributes
    ----------
    length : length of the neck
    radius : radius of resonator 
    length_correction: sum of upper and lower correction part
    fluid : fluid in the system
    end_impedance : function handle for end transfer impedance
    """
    
    air = mc.Fluid('air') #,damping_model='viscous')
    
    
    def __init__(self,length,radius,fluid = air,length_correction_factor = 0.85,**kwargs):
        """
        Class contructor for helholtzResonator

        Parameters
        ----------
        length : float
            length of the neck
        radius : float
            radius of neck 
        length_correction_factor: float
            sum of upper and lower end correction as factor of radius
        fluid: fluid
            fluid in the system
        KwArgs: dict
            Arbitrary keyword parameter list
        end_impedance: function
             radiation impedance of end correction
                                        
        """
        
        
        self.radius      = radius     
        self.length      = length
        self.fluid       = fluid
        self.length_correction = length_correction_factor*radius
        self.end_impedance   = lambda omega: 0*omega
                    
        for kw in kwargs:
            if kw == 'end_impedance':
                self.end_impedance = kwargs[kw]
            else:
                raise ValueError('Unkown keyword {0}'.format(kw))
                
    @property
    def area(self):
        return np.pi*self.radius**2
                
    def radiation_impedance(self,omega):
        """
        Radiation impedance of a QuarterWaveResonator Za = p/Q
        
        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        Z : complex
            radiation impedance.

        """

        _L = self.length_correction

        Z = self.end_impedance(omega) + \
            (1j*omega*self.fluid.rho0*_L + self.fluid.impedance(omega)/(1j*np.tan(self.fluid.wavenumber(omega)*self.length) ))/self.area 

        return Z
        
    def omega_resonance(self,n=0):
        """
        Angular resonance frequency
        
        Assumes the resistance free resonance frequency
                
        Returns
        -------
        float
            resonance angular frequency
        """
        
        return self.fluid.c0*(2*n+1)*np.pi/self.length

     
        