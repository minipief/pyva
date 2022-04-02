# -*- coding: utf-8 -*-
"""
Module for acoustic thee-dimensional systems or cavities

The acoustic3Dsystems module deals with the acoustic three dimensional systems.
This includes classes and methods for the simulation of rectangular rooms, but also 
for more generic SEA cavities 

All classes in this module extend the SEA_system class.
"""

import numpy as np
import pyva.properties.materialClasses as mc
import pyva.systems.SEA_system as SEAsys
import pyva.coupling.junctions as con

import pyva.useful as uf


class Acoustic3DSystem(SEAsys.SEA_system):
    """
    class for acoustical 3D SEA system or cavities
        
    Attributes
    ----------
    volume : float
        volume of the cavity
    surface : float
        surface area of the cavity
    perimeter : float
        perimeter of the cavitiy
    fluid : fluid
        fluid of the cavity
    absorption_area : float
        cumulative absorption area
    damping_type : str
        identifier for damping type
    """

    def __init__(self,ID,volume,surface,perimeter,fluid,absorption_area = 0.,damping_type = ['eta',]): 
        """
        Constructor for Acoustic3DSystem

        Parameters
        ----------
        ID : int
            ID of SEA system.
        volume : float
            volume of SEA cavity.
        surface : float
            surface of SEA cavitiy.
        perimeter : float
            perimeter of SEA cavity.
        fluid : fluid
            fluid of SEA cavity.
        absorption_area : float of Signal, optional
            absorption area of cavity. The default is 0..
        damping_type : list of str, optional
            identifyer for dampping method. The default is ['eta',].

        Returns
        -------
        None.

        """

               
        super().__init__(ID,[0],typestr='pressure',eta=fluid.eta)


        self.volume = volume
        self.surface = surface
        self.perimeter = perimeter
        self.fluid = fluid
        self.damping_type = damping_type
        # Potentialy frequency dependent parameters 
        self._absorption_area = absorption_area
        
    def __repr__(self):
        _str = 'Acoutic3DSystem({0},{1},{2},{3},{4})'.format(self.ID,self.volume,self.surface,self.perimeter,repr(self.fluid))
        return _str

    def __str__(self):
        _str = 'SEA cavity system with ID:{0} \n'.format(self.ID)
        _str += 'volume          : {}\n'.format(self.volume)
        _str += 'surface         : {}\n'.format(self.surface)
        _str += 'perimeter       : {}\n'.format(self.perimeter)
        _str += 'fluid:\n------\n{}------\n'.format(self.fluid)
        _str += 'damping_type    : {}\n'.format(self.damping_type)
        return _str
        

        
    def modal_density(self, omega, wave_DOF = 0):
        """
        Modal density estimation
                
        Parameters
        ----------
        omega : float 
            angular frequency
        wave_DOF : int
            wave degree of freedom. Default is 0.
                
        Returns
        -------
            Modal density
        """
        
        if wave_DOF != 0:
            raise ValueError('modal density of acoustic systme require wavedof = 0')
        
        V = self.volume
        S = self.surface
        P = self.perimeter
        c = np.real(self.fluid.c_freq(omega))

        return V*omega**2/(2*np.pi**2*c**3)+S*omega/(8*np.pi*c**2)+P/(16*np.pi*c)
    
    def absorption_area(self,omega):
        """
        Absorption area of cavity

        In case of constant absorption area the frequency paramter
        is ignored.

        Parameters
        ----------
        omega : float 
            angular frequency
        

        Returns
        -------
        absorption area

        """
        
        if uf.isscalar(self._absorption_area):
            return self._absorption_area
        else:
            # requires inter1 function for Signals
            return self._absorption_area.interp(omega)
        
        

    def damping_loss(self,omega,wave_DOF=0):
        """
        Damping loss of cavity SEA systemns
        
        Parameters
        ----------
        omega : float
            frequency 
        wave_DOF: int
            dummy argument (to full fill SEA class requirements)
            
        """
        
        eta = np.zeros(np.shape(omega))
        
        for it,dtype in enumerate(self.damping_type):
            if dtype == 'eta':
                eta[:] += self.eta
            elif dtype == 'surface':
                eta[:] += self.absorption_area(omega)*np.real(self.fluid.c_freq(omega))/4/self.volume/omega
        return eta
    
    def physical_unit(self,omega,energy):
        """
        Physical unit pressure calculated from system energy

        Parameters
        ----------
        omega : ndarray
            angular frequency

        Returns
        -------
        rms pressure

        """

        # real values considered because they only contribute to energy!
        c   = np.real(self.fluid.c_freq(omega))
        rho = np.real(self.fluid.rho_freq(omega))                
        
        return c*np.sqrt(rho*energy/self.volume)    
    
    
    def iscavity(self):
        """
        Check if SEA system is a cavity 

        Returns
        -------
        bool
            True.

        """
        
        return True


class RectangularRoom(Acoustic3DSystem):
    """
    This class deals with three dimensional, rectangular rooms filled with fluid
    
    Practical cavities are rarely. This class is applied in 
    chapters 4 and 6 [Pei2022]_ for demonstration of complex and random systems.
    
    Attributes
    ----------
    Lx : float
        Room length in x-direction
    Ly : float
        Room width  in y-direction
    Lz : float
        Room height in z-direction
    fluid: fluid
        fluid in the room
    """
    
    def __init__(self,ID,Lx,Ly,Lz,fluid,absorption_area = 0.,damping_type = ['eta',]):
        """
        Class contructor for acoustic tube
        
        Parameters
        ----------
        Lx : float
            Room length in x-direction
        Ly : float
            Room width  in y-direction
        Lz : float
            Room height in z-direction
        fluid : Fluid
            fluid in the room
        absorption_area : float or function, optional
            cumulative absorption area. The default is 0..
        damping_type : list of str, optional
            damping type identifier. The default is ['eta',].        
        """
        
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        
        super().__init__(ID,Lx*Ly*Lz,2*(Lx*Ly+Ly*Lz+Lx*Lz),4*(Lx+Ly+Lz),fluid,\
                         absorption_area=absorption_area,damping_type=damping_type)
                
    def p_mode(self,r,n):
        """
        Modal pressure shape for room with rigid walls
        
        Implementation of simple recangular room model 

        Parameters
        ----------
        r : list or tuple of float
            (x,y,z) room coodinates
        n : int
            number of mode
                
        Returns
        -------
            pressure distribution of mode         
        
        """

        aleph = 8
        lx = self.Lx
        ly = self.Ly
        lz = self.Lz
        V  = self.volume
        for i in range(3): # loop over room dimensions
            if n[i]<1:
                aleph/=2
        return np.sqrt(aleph/V)*np.cos(n[0]*np.pi/lx*r[0])*np.cos(n[1]*np.pi/ly*r[1])*np.cos(n[2]*np.pi/lz*r[2])
        
        
        
    def pms_modal_random(self,omega,Q=1,As=0):
        """
        Mean square pressure for diffuse wave field 
        
        Simple 1-system SEA cavity model for demonstration purpose
        mainly.

        Parameters
        ----------
            omega: angular frequency
                Q: source strength
               As: absorption area
                
        Returns
        -------
        mean square pressure
        """
        
        k  = np.real(self.fluid.wavenumber(omega))
        z0 = np.real(self.fluid.impedance())
        eta0= self.fluid.eta + As*self.fluid.c0/4/self.volume/omega
        return k*z0*z0/(4*np.pi*eta0*self.volume*Q**2)
        


    def k_mode_square(self,n):
        """
        Squared modal wavenumber of mode n 
        
        Parameters
        ----------
        n : tuple of int
            (nx,ny,nz) tupel of mode number
                
        Returns
        -------
            squared modal wavenumber
        """
        
        return np.pi*np.pi*((n[0]/self.Lx)**2+(n[1]/self.Ly)**2+(n[2]/self.Lz)**2)
        
    def f_mode(self, n):
        """
        modal frequency of mode n 
                
        Parameters
        ----------
        n : tuple of int
            (nx,ny,nz) tupel of mode number
                
        Returns
        -------
        modal frequency
        """

        return np.sqrt(self.k_mode_square(n))*self.fluid.c0/np.pi/2
        
    def omega_mode(self, n):
        """
        modal angular frequency of mode n 
        
        Parameters
        ----------
        n : tuple of int
            (nx,ny,nz) tupel of mode number
                
        Returns
        -------
        modal angular frequency
        """
        return np.sqrt(self.k_mode_square(n))*self.fluid.c0
            
    def green_modal(self,omega,r,r0,N=(50,50,50),source=1,Bnn=lambda n: 0):
        """
        Modal Green function of room 
        
        Parameters
        ----------
        r : tuple of float
            receiver location r = (x,y,z)
        r0 : tuple of float
            source location  r0 = (x0,y0,z0)
        N : tuple of int
            (Nx,Ny,Nz) tupel of maximal mode number
        Bnn : function of mode index tuple
            modal damping loss factor
                
        Returns:
            Modal approximation of Green function
        """
        p  = 0.
        c0   = self.fluid.c0 
        rho0 = self.fluid.rho0
        k    = omega/c0
        k2   = k**2
        for nx in range(N[0]):
            for ny in range(N[1]):
                for nz in range(N[2]):
                     n =  (nx,ny,nz)
                     p -= 1j*omega*rho0*source*self.p_mode(r,n)*self.p_mode(r0,n)/ \
                          (k2-self.k_mode_square(n)*(1+1j*self.fluid.eta)-1j*k*Bnn(n))
        return p
    
    def mechanical_point_impedance(self,omega,r0,N,S):
        """
        Modal mechanical point impedance
        
        This method requires a suface S to create a force from the acoustic impedance
        following from the Green function
        
        Parameters
        ----------
        omega : float
            angular frequency
        r0 : tuple of float
            source location  r0 = (x0,y0,z0)
        N : tuple of int
            (Nx,Ny,Nz) tuple of maximal mode number
          
        Returns
        -------
        Mechancal point impedance from modal approximation
        """
        
        return self.green_modal(omega,r0,r0,N,S)*S # siehe S als quellstärke

    def radiation_point_impedance(self,omega,r0,N,Bnn=lambda: 0):
        """
        Radiation point impedance
                
        Parameters
        ----------
        omega: float
            angular frequency
        r0 : tuple of float
            source location (x0,y0,z0)
        N : tuple of int
            (Nx,Ny,Nz) tuple of maximal mode number
        Bnn: function
            modal damping loss factor
             
        Returns
        -------
        acoustical point impedance from modal approximation
        """
        return self.green_modal(omega,r0,r0,N,1)
               
    def estimated_mode_count(self,omega):
        """
        Number of mode estimation
                
        Parameters
        ----------
        omega : float
            angular frequency
                
        Returns
        -------
        Number of modes below omega
        """

        fc = omega/(2*np.pi*np.real(self.fluid.c_freq(omega)))
        V = self.volume
        S = self.surface
        P = self.perimeter

        return 4*np.pi*V/3*(fc)**3+np.pi*S/4*(fc*fc)+P/8*fc 
        
    def mode_count(self,omegaVec):
        """
        Exact number count
        
        This method counts modes until angular limit frequencies
                
        Parameters
        ----------
        omegaVec: ndarray
            angular frequency vector
                
        Returns
        -------
        Mode count vector
        """
        
        iN=0; nx=0; ny=0; nz=0
        maxOmega = np.max(omegaVec);
        klim = np.real(self.fluid.wavenumber(maxOmega))**2
        
        Nest   = round(1.3*self.estimated_mode_count(maxOmega))
        omegas = np.tile(maxOmega+1.,Nest)
        N      = np.zeros(len(omegaVec))
        while self.k_mode_square((nx,ny,nz))<=klim:
            while self.k_mode_square((nx,ny,nz))<=klim:
                while self.k_mode_square((nx,ny,nz))<=klim:
                    omegas[iN] = self.omega_mode((nx,ny,nz))
                    nz+=1        
                    iN+=1
                    print(iN,omegas[iN-1])
                ny+=1
                nz =0
            nx +=1
            ny =0
                    
        for iom in range(len(omegaVec)):
            indexes = np.where(omegas<omegaVec[iom],1.,0. )
            N[iom]  = np.sum(indexes)
            
        return N

    def modal_density_precise(self,omega):
        """
        Exact modal density
        
        This method derives the modal density by counting modes in given bands
                
        Parameters
        ----------
        omega : ndarray
            angular frequency vector
                
        Returns
        -------
        Modal density vector
        """
        
        delta = omega[1:]-omega[:-1]
        Ns    = self.mode_count(omega)
        Ns    = Ns[1:]-Ns[:-1]
        return Ns/delta, omega[0:-1]+delta/2
            
    def surfint(self,y0,y1,z0,z1,n):
        """
        Surface integral required for rectangluar piston
        
        Integral implementation of equation (4.92) from [Pei2022]_
        
        Parameters
        ----------
        y0 : float
            minimum y-coodinate of piston
        y1 : float
            maximum y-coodinate of piston
        z0 : float
            minimum z-coodinate of piston
        z1 : float
            maximum z-coodinate of piston
        omega : float
            angular frequency vector
        n : tuple of int
            mode number n = (nx,ny,nz)
                
        Returns
        -------
        integral value
        """
        
        ly = self.Ly
        lz = self.Lz
        lx = self.Lx
        
        aleph = 8
        if n[0]<1:
            aleph/=2
    
        if n[1] > 0:
            cy = n[1]*np.pi/ly
            dy = (np.sin(cy*y1)-np.sin(cy*y0))/cy
        else:
            dy = np.abs(y1-y0)
            aleph/=2
            
        if n[2] > 0:
            cz = n[2]*np.pi/lz
            dz = (np.sin(cz*z1)-np.sin(cz*z0))/cz
        else:
            dz = np.abs(z1-z0)
            aleph/=2
            
        return dy*dz/np.sqrt(self.volume)*np.sqrt(aleph)
        #return dy*dz*aleph/self.volume
           
    def modal_impedance_rectangular(self,omega,y0,y1,z0,z1,n):
        """
        Mechanical radiation impedance of rectangular piston
        
        Parameters
        ----------
        y0 : float
            minimum y-coodinate of piston
        y1 : float
            maximum y-coodinate of piston
        z0 : float
            minimum z-coodinate of piston
        z1 : float
            maximum z-coodinate of piston
        omega : float
            angular frequency vector
        n : tuple of int
            mode number n = (nx,ny,nz)
                
        Returns
        -------
        Modal impedance for mode n
        """

        k2 = (omega/self.fluid.c0)**2
        c1 = -1j*omega*self.fluid.rho0/(k2-self.k_mode_square(n)*(1+1j*self.fluid.eta))

        return c1*self.surfint(y0,y1,z0,z1,n)
    
    def Z_rect(self,omega,y0,y1,z0,z1,N):
        """
        Mechanical radiation impedance of rectangular piston
        
        Parameters
        ----------
        y0 : float
            minimum y-coodinate of piston
        y1 : float
            maximum y-coodinate of piston
        z0 : float
            minimum z-coodinate of piston
        z1 : float
            maximum z-coodinate of piston
        omega : float
            angular frequency vector
        N : tuple of int
            maximal mode number N = (Nx,Ny,Nz)
                
        Returns
        -------
        Mechanical impedance
        """
        
        Z = 0.
        A = np.abs(y1-y0)*np.abs(z1-z0)
        for nx in range(N[0]):
            for ny in range(N[1]):
                for nz in range(N[2]):
                    n =  (nx,ny,nz)
                    Z += self.modal_impedance_rectangular(omega,y0,y1,z0,z1,n)*self.surfint(y0,y1,z0,z1,n)/A
        
        return Z
