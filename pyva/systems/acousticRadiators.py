# -*- coding: utf-8 -*-
"""
Module for acoustic sources

The acousticRadiators module defines several classes of radiotors into homogeneous
fluids, sometimes infinity but often into semi-infinity half spaces.

"""

import numpy as np
import types

import pyva.properties.materialClasses as mc
import pyva.data.matrixClasses as mC
import pyva.geometry.meshClasses as meshC
import matplotlib.pyplot as plt
import scipy.special as spl
import scipy.integrate as integrate
import pyva.data.dof as dof
import pyva.useful as uf


class CircularPiston:
    """
    The CircularPiston class defines the dynamics and acoustics of a piston 
    in a rigid wall radiating into a semi-infinite fluid
    
    Attributes
    ----------
    radius : float
        Radius of piston
    fluid : fluid
        fluid of half space 
    """
    
    def __init__(self,radius,fluid = mc.Fluid()):
        """
        Class contructor for CircularPiston
        
        Parameters
        ----------
        radius : float
            radius of the piston
        fluid : fluid
            fluid in half space
                                        
        """

        self.radius = radius
        self.fluid  = fluid
        
    def __str__(self):
        """
        str for CircularPiston

        Returns
        -------
        _str : str
            CircularPiston description.

        """
        #_str  = "name            : {0}\n".format(self.name)
        _str  = "CircularPiston:\n"
        _str += "radius:{0}\n".format(self.radius)
        _str += "fluid:\n{0}".format(self.fluid)
        
        return _str

    def __repr__(self):
        _str = 'CircularPiston({0},fluid={1})'.format(self.radius,self.fluid.__repr__())
        return _str
    
    @property
    def area(self):
        """
        Area of Piston
        
        Returns
        -------
        Surface area of piston         
        """

        return np.pi*self.radius**2
        
    def acousticImpedance(self,omega):
        """
        Radiation impedance of circular piston Za = p/v
        
        The pressure is averaged over the piston surface
        
        Parameters
        ----------
        omega : float
            angular frequency
                               
        Returns
        -------
        acoustic radiation impedance of piston         
        
        """

        kr = omega/self.fluid.c0*self.radius
        return self.fluid.z0*(1-spl.j1(2*kr)/(kr)+1j*spl.struve(1,2*kr)/(kr))
 
    def radiationImpedance(self,omega):
        """
        Radiation impedance of circular piston Za = p/Q
        
        The pressure is averaged over the piston surface
        
        Parameters
        ----------
        omega : float
            angular frequency
                               
        Returns
        -------
        volume radiation impedance of piston         
        
        """

        return self.acousticImpedance(omega)/self.area

        
    def mechanicalImpedance(self,omega):
        """
        Mechanical impedance of cisular piston Z = F/v
        
        The pressure is average over the piston surface
        
        Parameters
        ----------
        omega : float
            angular frequency
                               
        Returns
        -------
        mechanical impedance of radiating piston         
        
        """
        return self.area*self.acousticImpedance(omega)
        
    def pressure(self,omega,dist,theta=0):
        """
        Field pressure of cisular piston
        
        Parameters
        ----------
        omega : float
            angular frequency
        dist : float
            distance to piston center
        theta : float
            angle to piston normal
                               
        Returns
        -------
        field pressure of acoustic field in fluid         
        
        """

        c0 = self.fluid.c0; rho0 = self.fluid.rho0
        k  = self.fluid.wavenumber(omega)
        r  = dist
        ka = k*self.radius
        return 1j*k*c0*rho0/(2*np.pi*r)*spl.jn(1,ka*np.sin(theta))/(ka*np.sin(theta))
    
    def acoustic_FE(self,omega,ID=[1],**kwargs):
        """
        Acoustic Point Element of piston radiators
        
            A p = q with A = 1/Zrad
        
        Uses the priston radiation function to simulate a circular open end
        in a baffle
        
        Parameters
        ----------
        omega : float
            angular frequency
                
        Returns
        -------
        Dynamic matrix
            [1 x 1 x Nfreq]          
        
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
            data[0,iomega] = 1./self.radiationImpedance(omega[iomega])
              
        return mC.DynamicMatrix(data,xdata,excdof,resdof,sym=1,shape=(1,1,len(omega)))
          
class Sphere:
    """
    The Sphere class defines an acoustic sphere 
    that may operate as breathing or oscillatinf sphere
    
    Attributes
    ----------
    radius : float
        Radius of sphere
    fluid : fluid
        fluid in free space 
    """

    def __init__(self,radius,fluid=mc.Fluid()):
        """
        Class contructor for sphere
        
        Parameters
        ----------
        radius : float
            radius of the sphere
        fluid : fluid
            fluid in free space
                
        """

        self.radius = radius
        self.fluid  = fluid 
        
    def radiation_impedance(self,omega):
        """
        Radiation impedance of breathing sphere Zrad = p/Q
        
        The pressure is average over the sphere surface
        
        Parameters
        ----------
        omega : float
            angular frequency
                               
        Returns:
            radiation impedance of sphere         
        
        """
        area = 4*np.pi*self.radius**2
        return self.acoustic_impedance(omega,self.radius)/area
    
    def velocity_potential(self,omega,dist,Q=1.):
        """
        Velocity potential of breathing sphere

        Parameters
        ----------
        omega : float
            angular frequency.
        dist : float
            distance to sphere center.
        Q : complex or nd.array of complex, optional
            Source streng in volume flow rate, volume per time. The default is 1. .

        Returns
        -------
        complex or nd.array of complex
            velocity potential.

        """

        k = self.fluid.wavenumber(omega)
        c_1 = -Q/4/np.pi/dist / (1 + 1j*k*self.radius)
        
        return c_1*np.exp(-1j*k*(dist-self.radius)) 
    
    
    def pressure(self,omega,dist,Q=1.):
        """
        Pressure of breathing sphere

        Parameters
        ----------
        omega : float
            angular frequency.
        dist : float
            distance to sphere center.
        Q : complex or nd.array of complex, optional
            Source streng in volume flow rate, volume per time. The default is 1. .

        Returns
        -------
        complex or nd.array of complex
            pressure of breathing sphere.

        """
        
        k = self.fluid.wavenumber(omega)
        z0 = self.fluid.impedance(omega)
        return -1j*k*z0*self.velocity_potential(omega,dist)
    
    def velocity(self,omega,dist,Q=1.):
        """
        Velocity of breathing sphere

        Parameters
        ----------
        omega : float
            angular frequency.
        dist : float
            distance to sphere center.
        Q : complex or nd.array of complex, optional
            Source streng in volume flow rate, volume per time. The default is 1. .

        Returns
        -------
        complex or nd.array of complex
            velocity of breathing sphere.

        """
        
        k   = self.fluid.wavenumber(omega)
        c_1 = -(1+1j*k*dist)/dist
        return c_1*self.velocity_potential(omega,dist)
    
    def acoustic_impedance(self,omega,dist):
        """
        Acoustic impedance of breathing sphere

        Parameters
        ----------
        omega : float
            angular frequency.
        dist : float
            distance to sphere center.

        Returns
        -------
        complex or nd.array of complex
            acoustic impedance of sphere soud field.

        """
        kr = self.fluid.wavenumber(omega)*dist
        z0 = self.fluid.impedance(omega)
        
        return 1j*z0*kr/(1+1j*kr)
        
        
                
class Monopole:
    """
    The Monopole class defines an acoustic momopole
    
    Attributes:
            Q:   Source strength in volume flow rate amplitude Q = dV/dt
        fluid:   fluid of free space 
    """
    def __init__(self,Q,fluid=mc.Fluid):
        """
        Class contructor for Monopole
        
        Args:
            Q:     Source strength in volume flow rate amplitude Q = dV/dt
            fluid: fluid in free space
                
        Examples:
            import acousticRadiators as ar
            myMono = ar.Monopole(1.0)
                        
        """
        self.Q = Q
        self.fluid = fluid
        
    def radiation_impedance(self,omega):
        """
        Acoustic radiation impedance of Monopole Zrad = p/Q
        
        The pressure is average over the sphere surface
        
        Parameters
        ----------
        omega : float
            angular frequency
                               
        Returns:
            radiation impedance of sphere         
        
        """
        Z0 = np.abs(self.fluid.impedance())
        Z = np.real(self.fluid.wavenumber(omega))**2*Z0/4./np.pi
        return Z
        
    def pAmp(self,omega,r,Q=1):
        return self.fluid.wavenumber(omega)*self.fluid.impedance(omega)*Q/4/np.pi/r
        
class HalfSpace:
    """
    Class for HalfSpace operations

    Attributes:
        fluid:      fluid in free space
        treatment:  treatment for radiating area
    
    """
    def __init__(self,fluid=mc.Fluid(),treatment = 'none'):
        """
        Class contructor for HalfSpace
        
        Parameters
        ----------
        fluid : Fluid
            of half space
        treatment: str, Signal 
            treatment between vibrating surface and half space
                
        Examples
        --------
            import acousticRadiators as aR
            myHalfSpace = aR.HalfSpace()
                        
        """

        self.fluid = fluid
        self.treatment = treatment
        
    def __str__(self):
        """
        str for HalfSpace

        Returns
        -------
        _str : str
            Fluid description.

        """
        #_str  = "name            : {0}\n".format(self.name)
        _str  = "HalfSpace:\n"
        _str += "fluid:\n{0}".format(self.fluid)
        _str += "treatment: \n {0}\n".format(self.treatment)
        
        return _str

    def __repr__(self):
        _str = 'HalfSpace({0},treatment={1})'.format(self.fluid.__repr__(), self.treatment.__repr__())
        return _str
        
    def radiation_stiffness_wavenumber(self,omega,kx,ky=0):
        """
        Acoustic radiation stiffness of half space in wavenumber domain D = F/x
                
        Parameters
        ----------
        omega : float
            angular frequency
        kx : float
            x-wavenumber of exciting wave
        ky : float
            y-wavenumber of exciting wave
                               
        Returns
        -------
        radiation stiffness of semi infinity fluid half space        
        """
        
        return -1j*omega*self.radiation_impedance_wavenumber(omega,kx,ky)

    def radiation_impedance_wavenumber(self,omega,kx,ky=0):
        """
        Acoustic radiation impedance of half space in wavenumber domain Z = p/v
                
        Parameters
        ----------
        omega : float
            angular frequency
        kx : float
            x-wavenumber of exciting wave
        ky : float
            y-wavenumber of exciting wave
                               
        Returns
        -------
        radiation stiffness of semi infinity fluid half space        
        """
        
        ka = np.real(self.fluid.wavenumber(omega))
        ka2 = ka*np.conj(ka)
        kz = np.lib.scimath.sqrt(ka2-kx*kx-ky*ky)
        
        
        #om2 = omega*omega
        return omega*self.fluid.rho0/kz
        

        
    def radiation_stiffness_piston(self,omega,dist,dA1=1.0,dA2=1.0):
        """
        Acoustic radiation stiffness of separated vibrating source D = F(x1)/w(X2)
        
        This method returns the force generated at x1 by a vibrating element at x2.
                
        Parameters
        ----------
        omega : ndarray or float
            angular frequency
        dist : ndarray or float
            distance between x1 and x2 dist = abs(x1-x2)
        dA1 : float
             Element area at x1
        dA2 : float
            Element area at x2
        method : str
            Identifier for calculation method 'piston' or 'wavelett'
            
        Returns
        -------
        res : ndarray of either len(omega) or len(dist)
            radiation stiffness        
        """
        
        if uf.isscalar(omega):
    
            dist = np.array(dist).flatten()
            # radius according to element area
            r1 = np.sqrt(dA2/np.pi) 
            ka = np.real(self.fluid.wavenumber(omega))
            ix = dist != 0.
            kr = ka*dist[ix]
            rhoF = self.fluid.rho0 #  Freq(omega)
    
            res = np.zeros(dist.shape,dtype=np.complex128)
            
            res[np.logical_not(ix)] = omega*omega*self.fluid.rho0*dA1*dA2/np.pi*(1j*omega/2/self.fluid.c0-8/3/np.pi/r1) 
            res[ix] = omega**3*dA1*dA2*rhoF/(2*np.pi*self.fluid.c0)*(-np.cos(kr)/kr + 1j*np.sinc(kr/np.pi))
        elif uf.isscalar(dist):
            # radius according to element area
            r1 = np.sqrt(dA2/np.pi) 
            ka = np.real(self.fluid.wavenumber(omega))
            kr = ka*dist
            rhoF = self.fluid.rho0 #  Freq(omega)
    
            #res = np.zeros(omega.shape,dtype=np.complex128)
            
            if dist == 0:
                res = omega*omega*self.fluid.rho0*dA1*dA2/np.pi*(1j*omega/2/self.fluid.c0-8/3/np.pi/r1)
            else:
                res = omega**3*dA1*dA2*rhoF/(2*np.pi*self.fluid.c0)*(-np.cos(kr)/kr + 1j*np.sinc(kr/np.pi))
            
        return res
            


    def radiation_stiffness_wavelet(self,omega,dist,ks):
        """
        Acoustic radiation stiffness of separated vibrating source D = F(x1)/w(X2)
        
        This method implements the wavelet method from Langley [Lan2007]_

        Parameters
        ----------
        omega : ndarray or float
            angular frequency
        dist : ndarray or float
            distance between radiating elments.
        ks : float
            wavenumber parameter.

        Returns
        -------
        res : ndarray of either len(omega) or len(dist)
            radiation stiffness        

        """

        # Workaround for itj0y0 using Section 6.511 Eq. 6 of I.S. Gradshteyn and I.M. Ryzhik, Tables of intgrals
        def itj0(x):
            return x*spl.j0(x)+np.pi*x/2*(spl.j1(x)*spl.struve(0,x)-spl.j0(x)*spl.struve(1,x))

        if uf.isscalar(omega):
    
            dist = np.array(dist).flatten()

            #ks = 2*np.pi/np.sqrt(dx*dx+dy*dy)
            ka = np.real(self.fluid.wavenumber(omega))
            ix = dist != 0.
            kr = ka*dist[ix]
            
            fak = 2*np.pi**3*omega**3*self.fluid.rho0/(ks**4*self.fluid.c0)
            
            # @todo The current itj0y0 does not work - reuse when issue soves
            #bes,_ = spl.itj0y0(dist[ix]*ks)
            bes = itj0(dist[ix]*ks)
            res = np.zeros(dist.shape,dtype=np.complex128)
            
            res[ix]    = fak*(1j*np.sinc(kr/np.pi)-(np.cos(kr)-1+bes)/kr)
            res[np.logical_not(ix)] = fak*(1j-ks/ka)
        elif uf.isscalar(dist):
            
            #ks = 2*np.pi/np.sqrt(dx*dx+dy*dy)
            ka = np.real(self.fluid.wavenumber(omega))
            kr = ka*dist
            
            fak = 2*np.pi**3*omega**3*self.fluid.rho0/(ks**4*self.fluid.c0)
            # @todo The current itj0y0 does not work - reuse when issue soves
            #bes,_ = spl.itj0y0(dist[ix]*ks)
            bes = itj0(dist[ix]*ks)
            res = np.zeros(omega.shape,dtype=np.complex128)
            
            if dist != 0:
                res    = fak*(1j*np.sinc(kr/np.pi)-(np.cos(kr)-1+bes)/kr)
            else:
                res = fak*(1j-ks/ka)
                
        return res         
    
    def get_mesh(self,omega_max,Lx,Ly,N=4):
        """
        Creates mesh for HalfSpace
        
        The element size is calculated from the wavelength of the fluid.

        Parameters
        ----------
        omega_max : double
            Maximum angular frequency for determination of smallest wavelength
        Lx : float
            length in x direction
        Ly : float
            length in y direction
        N : interger, optional
            Minimal number of nodes per wavelength

        Returns
        -------
        mesh : regmesh2D
            The mesh covering the given conditions

        """
        
        # determine wavelength of fluid
        lam = self.fluid.wavelength(omega_max)
        Dmin = lam/N
        NX = int(np.ceil(Lx/Dmin+1))
        NY = int(np.ceil(Ly/Dmin+1))
    
        mesh = meshC.RegMesh2D(0., 0., Lx, Ly, NX, NY)
        
        return mesh
        
    def radiation_stiffness_mesh(self,omega,mesh,method='wavelet'):
        """
        Acoustic radiation stiffness matrix of a regular mesh
                        
        Parameters
        ----------
        omega : float
            angular frequency
        mesh : mesh
            mesh of surface points
        method : str
            Identifier for calculation method 'piston' or 'wavelet'
            
        Returns
        -------
        LinearMatrix
            radiation stiffness matrix of semi infinite fluid half space        
        
        """
        
        
        distances, index = mesh.distance()
        
        Ndist  = len(distances)
        Nx     = len(omega)
        Nmesh  = mesh.Nmesh
        dA     = mesh.dA
        ks     = mesh.ks
        
        _res    = np.zeros(((Nmesh*(Nmesh+1))//2,Nx),dtype=np.complex128)
        _radbuf = np.zeros(Ndist,dtype=np.complex128) 
        
        if method == 'piston':
           for i in range(Nx):
                #for id in range(Ndist):
                _radbuf = self.radiation_stiffness_piston(omega[i],distances,dA,dA)
                    
                _res[:,i] = _radbuf[index]
        elif method == 'wavelet':
            for i in range(Nx):
                #for id in range(Ndist):
                _radbuf = self.radiation_stiffness_wavelet(omega[i],distances,ks)
                    
                _res[:,i] = _radbuf[index]
        else: 
            raise ValueError('method: {0} not known'.format(method))
                
            
        return mC.LinearMatrix(_res,sym = 1,shape=(Nmesh,Nmesh,Nx) )

    def radiation_stiffness_mesh_single(self,omega,mesh,method='piston'):
        """
        Acoustic radiation stiffness matrix of a regular mesh
                        
        Parameters
        ----------
        omega : float
            angular frequency
            mesh:  of surface points
            method: Identifier for calculation method 'piston' or 'wavelet'
            
        Returns:
            radiation stiffness matrix of semi infinite fluid half space @type DynamicMatrix        
        
        """
        
        
        distances, index = mesh.distance()
        
        Ndist  = len(distances)
        Nx     = len(omega)
        Nmesh  = mesh.Nmesh
        dA     = mesh.dA
        ks     = mesh.ks
        
        _res    = np.zeros(((Nmesh*(Nmesh+1))//2,Nx),dtype=np.complex128)
        _radbuf = np.zeros(Ndist,dtype=np.complex128) 
        
        if method == 'piston':
           for i in range(Nx):
                #for id in range(Ndist):
                _radbuf = self.radiation_stiffness_piston(omega[i],distances,dA,dA)
                    
                _res[:,i] = _radbuf[index]
        elif method == 'wavelet':
            for i in range(Nx):
                #for id in range(Ndist):
                _radbuf = self.radiation_stiffness_wavelet(omega[i],distances,ks)
                    
                _res[:,i] = _radbuf[index]
        else: 
            raise ValueError('method: {0} not known'.format(method))
                
            
        return mC.LinearMatrix(_res,sym = 1,shape=(Nmesh,Nmesh,Nx) )

            
            
    def shape_radiation_power(self,omega,shape,method='piston',normalise=True):
        
        D = self.radiationStiffnessMesh(omega,shape)
        
        Nx     = len(omega)
        #dA     = shape.dA
        
        
        
        ydata  = np.zeros(Nx,dtype=np.float64)
        
        for iw in range(len(omega)):

            Drad      = D.Dindex(iw)
       
            ydata[iw] = 0.5*omega[iw]*np.imag(shape.Zvec(False).conj().transpose() @ Drad @ shape.Zvec(False))
            
            # deal with shapes
         
        xdata = mC.DataAxis(omega,typestr='angular frequency')
            
        return mC.Signal(xdata,ydata,dof.DOF(0,0,dof.DOFtype(typestr='sound power')))

            
            
    def shape_radiation_efficiency(self,omega,shape,method='piston',normalise='True', signal = True):
        """
        Shape radiation stiffness.
        
        The shape radiation stiffness is the radiation stiffness defined on the
        generalised shape coordinates of the surface, for example the mode shape 
        of a plate

        Parameters
        ----------
        omega : float
            angular frequency
        shape : RegShape2D
            radiating shape.
        method : str, optional
            method identifier. The default is 'piston'.
        normalise : bool, optional
            switch for shape nrormalisation. The default is 'True'.
        signal : bool, optional
            switch for signal output. The default is 'True'.
            

        Returns
        -------
        Signal
            radiation stiffness.

        """
        
        D = self.radiation_stiffness_mesh(omega,shape,method)
        
        Nx     = len(omega)
        dA     = shape.dA
        
        c2 = shape.mean_square(normalise)
                
        ydata = np.zeros(Nx,dtype=np.float64)
        
        _shape = shape.Zvec(normalise)
        
        for iw in range(len(omega)):
            Drad      = D.Dindex(iw)
            const     = dA*self.fluid.z0*omega[iw]*c2*2
                        
            ydata[iw] = np.imag(_shape.conj().transpose() @ Drad @ _shape/const)
            
            # deal with shapes
         
            
        if signal:
            xdata = mC.DataAxis(omega,typestr='angular frequency')
            return mC.Signal(xdata,ydata,dof.DOF(0,0,dof.DOFtype(typestr='general')))
        else:
            return ydata
            
    
    def radiation_efficiency_leppington(self,omega,kx,ky,Lx,Ly,simple_muGT1=False):
        """
        Acoustic radiation stiffness of half space in wavenumber domain  using Leppingtons
        Approximation.
        
        See [Pei2022]_ and [Lep1982]_ for details
                
        Parameters
        ----------
        omega : float
            angular frequency
        kx : float
            x-wavenumber of shape
        ky : float
            y-wavenumber of shape
        Lx : float
            length in x-direction
        Ly : float
            length in y-direcition
                                           
        Returns
        -------
        radiation efficiency        
        
        """
        
        if omega.size ==1:
            omega_vec = False # omega is skalar and kx and ky scalar
            sigma = np.zeros(np.size(kx),float)
        else:
            omega_vec = True # omega is vector and kx and ky scalar
            sigma = np.zeros(np.size(omega),float)
            
        ka  = np.real(self.fluid.wavenumber(omega))
        

        alpha = kx/ka
        beta  = ky/ka
        # this should be a vector in any case
        mu    = np.sqrt(alpha**2+beta**2)
        if Lx < Ly:
            eps1  = np.pi/ka/Lx*0.9
        else:
            eps1  = np.pi/ka/Ly*0.9
            
        
        #if Lx < Ly:
        # #    a1 = Lx
        # #    a2 = Ly
        # #else:
        #    Lx,Ly = Ly,Lx

        def A(alpha,beta,omega_vec): 
            """
            eq. (4.29)
            """
            a2b2 = alpha**2+beta**2
            # switch formula depending if ky or ka is scalar
            if omega_vec:
                return 0.5/(Ly*ky*np.sqrt(a2b2-1))*(1+beta**2/(a2b2-1)) 
            else:
                return 0.5/(Ly*ka*beta*np.sqrt(a2b2-1))*(1+beta**2/(a2b2-1)) 

                
        def B(alpha,beta,omega_vec):
            """
            eq. (4.30)
            """
            a2b2 = alpha**2+beta**2
            # switch formula depending if ky or ka is scalar
            if omega_vec:
                return 0.5/(Lx*kx*np.sqrt(a2b2-1))*(1+alpha**2/(a2b2-1))
            else:
                return 0.5/(Lx*ka*alpha*np.sqrt(a2b2-1))*(1+alpha**2/(a2b2-1))
                
        def N1(alpha,beta,ka,omega_vec):
            """
            N1 from eq (6.23)

            """
            out = np.zeros(len(alpha))
            ix  = Lx*beta<Ly*alpha
            if omega_vec:
                out[ix] = 0.5*ka[ix]*Lx/alpha[ix]*np.abs(eps(alpha[ix],beta[ix]))
                ix = np.logical_not(ix)
                out[ix] = 0.5*ka[ix]*Ly/beta[ix]*np.abs(eps(alpha[ix],beta[ix]))
            else:
                out[ix] = 0.5*ka*Lx/alpha[ix]*np.abs(eps(alpha[ix],beta[ix]))
                ix = np.logical_not(ix)
                out[ix] = 0.5*ka*Ly/beta[ix]*np.abs(eps(alpha[ix],beta[ix]))
                
            return out

        def T(alpha,beta):
            out     = np.zeros(len(alpha))
            ix      = Lx*beta<Ly*alpha
            out[ix] = Lx*beta[ix]/Ly/alpha[ix]
            ix = np.logical_not(ix)
            out[ix] = Ly*alpha[ix]/Lx/beta[ix]
            return out
            
            
        def eps(alpha,beta):
            return 2*(np.sqrt(alpha**2+beta**2)-1)
        
        def epssign(eps):
            return -np.sign(eps)
        
        def F(z):
            _s,_c = spl.fresnel(np.sqrt(2/np.pi)*z) 
            return (_c+1j*_s)*np.sqrt(np.pi/2)
        
        def Si(x):
            out = np.zeros(np.shape(x))
            for ii,x_ in enumerate(x):
               out[ii] = integrate.quad(lambda x: np.sinc(x/np.pi), 0, x_)[0]
               
            return out
        
        
        # important parameter for approximation aroud coincidence
        muTol = 0.1
        
        def Hfun(x):
            return 0.5-0.15*x
        

        
        # defining indexes for cases
        iMuLT1 = mu < 1 # condition for above coincidence
        iMuGT1 = mu > 1 
        # simple self defined tolerance
        #iMunear   = np.abs(mu-1) < muTol # this defineds the indexes near 1
        # eps1 use at unit circle
        iMunear  = np.abs(mu-1) < eps1
        # combination with all zones
        #
        # define specific area for points of intersection
        #iMuleft   = np.logical_and(iMunear,iMuGT1)
        #iMuright  = np.logical_and(iMunear,iMuLT1)
        
        # buffer for smooth takeover
        # sigma_near = np.zeros(np.shape(iMunear))
        iMu1 = mu == 1
        
        # Eq. (3.1) in [1] mu < 1
        sigma[iMuLT1] = 1/np.sqrt(1-mu[iMuLT1]**2)
        # simple version Eq. (7.7) mu > 1 
        if simple_muGT1:
            if omega_vec:
                sigma[iMuGT1] = (Lx+Ly)/(np.pi*mu[iMuGT1]*ka[iMuGT1]*Lx*Ly*np.sqrt(mu[iMuGT1]**2-1))\
                       *(np.log((mu[iMuGT1]+1)/(mu[iMuGT1]-1))+2*mu[iMuGT1]/(mu[iMuGT1]**2-1))
            else:
                sigma[iMuGT1] = (Lx+Ly)/(np.pi*mu[iMuGT1]*ka*Lx*Ly*np.sqrt(mu[iMuGT1]**2-1))\
                       *(np.log((mu[iMuGT1]+1)/(mu[iMuGT1]-1))+2*mu[iMuGT1]/(mu[iMuGT1]**2-1))
        else:
            # more indexes
            iregI   = np.logical_and(mu > 1, alpha>1, beta<1)
            iregII  = np.logical_and(mu > 1, alpha<1, beta>1)
            iregIII = np.logical_and(mu > 1, alpha<1, beta<1)
            iregIV  = np.logical_and(alpha>1, beta>1)
            sigma[iregI]   = 2*B(alpha[iregI],beta[iregI],omega_vec)
            sigma[iregII]  = 2*A(alpha[iregII],beta[iregII],omega_vec)
            sigma[iregIII] = 2*(A(alpha[iregIII],beta[iregIII],omega_vec)\
                                +B(alpha[iregIII],beta[iregIII],omega_vec))
                
            aa = alpha[iregIV]
            bb = beta[iregIV]
            pp = np.sqrt(aa**2+bb**2-1)
       
            buf1 = 2*pp**4+pp**2-3*(aa*bb)**2
            buf2 = 2*pp**4+pp**2+3*(aa*bb)**2

            if omega_vec:
                ka_ = ka[iregIV]
            else:
                ka_ = ka
            
            #eq. (4.41)
            #sigma[iregIV] = 1/(np.pi*ka_**2*Lx*Ly*aa*bb)*(2*aa*bb*buf1/(pp**4*(aa**2-1)*(bb**2-1)) + \
            #                                buf2/pp**5*np.log(np.abs((aa*bb+pp)/(aa*bb-pp))))
    
            # iC18  = np.abs(alpha-1)< eps1
            # aa = alpha[iC18] 
            # if omega_vec:
            #     ka_ = ka[iC18]
            # else:
            #     ka_ = ka
                
            # AA = 1/(ka_*Ly*beta[iC18]**2)
                
            # sigma[iC18] += 2*AA/np.pi*\
            #     (Si(ka_*Lx*(1-aa))-np.pi/2*np.sign(1-aa)-(1-np.cos(ka_*Lx*(1-aa)))/(ka_*Lx*(1-aa)))
        
            # iC23  = np.abs(beta-1) < eps1
            # bb = beta[iC23] 
            # if omega_vec:
            #     ka_ = ka[iC23]
            # else:
            #     ka_ = ka
                
            # BB = 1/(ka_*Lx*alpha[iC23]**2)
                
            # sigma[iC23] += 2*BB/np.pi*\
            #     (Si(ka_*Ly*(1-bb))-np.pi/2*np.sign(1-bb)-(1-np.cos(ka_*Ly*(1-bb)))/(ka_*Ly*(1-bb)))
            
            
              
                
                       
        # Eq. (6.21)
        _alpha= alpha[iMunear]
        _beta = beta[iMunear]
        _eps  = eps(_alpha,_beta)
        _sign = epssign(_eps)
        _exp4 = np.exp(_sign/4*1j*np.pi)
        if omega_vec:
            _N1   = N1(_alpha,_beta,ka[iMunear],omega_vec)
        else:
            _N1   = N1(_alpha,_beta,ka,omega_vec)
            
        _N1sqr= np.sqrt(_N1)
        _F    = F(_N1sqr)
        
        sigma_buf = _sign/np.sqrt(np.pi*np.abs(_eps))*\
            np.imag(_exp4*((2.-1j/_N1)*_F+1j/_N1sqr*np.exp(1j*_N1))+\
            T(_alpha,_beta)*_exp4*(1.5/(_N1*_N1sqr)*np.exp(1j*_N1)-(1j/_N1+3/2/_N1**2)*_F))
        

        #plt.figure(10)
        #plt.semilogx(omega,sigma)
        #plt.semilogx(omega[iMunear],sigma_buf)
        #plt.ylim((0,3))

        
        sigma[iMunear] = sigma_buf
        #ix_L, = np.nonzero(iMuleft)
        #ix_R, = np.nonzero(iMuright)
        #NL = len(ix_L)
        #NR = len(ix_R)
        
        #minL = np.amin(np.abs(sigma_buf[:NL] - sigma[iMuleft]))
        #minR = np.amin(np.abs(sigma_buf[-NR:] - sigma[iMuright]))
                                              
        #i_L, = np.where(minL == np.abs(sigma_buf[:NL] - sigma[iMuleft]))
        #i_R, = np.where(minR == np.abs(sigma_buf[-NR:] - sigma[iMuright]))
                
        #i_intersect_L = np.argwhere(sigma_buf <= sigma[iMunear]).flatten()
        #iMuL = ix_L[0]  # get 1st position of iMunear
        #iMuR = ix_R[0] # get 1st position of iMunear
        
        #ix_near = np.arange(i_L[-1],i_R[0]+NL,dtype = int)
        
        #sigma[iMuL+ix_near] = sigma_buf[ix_near] 
        #sigma[iMunear] = sigma_buf 

        #plt.semilogx(omega[iMuL+ix_near],sigma_buf[ix_near],'r-')
        #plt.semilogx(omega[iMunear],sigma_buf,'r-')
        
        #plt.show     
        
        return sigma
              
                
