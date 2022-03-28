# -*- coding: utf-8 -*-
"""
Module for structural two-dimensional systems

Generally aimed at two-dimensional systems this module includes exclusively 
flat plates. In the future more realistic system will be included
"""

import numpy as np
import scipy.integrate as integrate
import scipy.linalg as linalg
import scipy.special as special
import pyva.properties.structuralPropertyClasses as stPC
import pyva.properties.materialClasses as mc
import pyva.data.matrixClasses as mC
import pyva.data.dof as dof
import pyva.systems.acousticRadiators as aR
import pyva.useful as uf
import pyva.systems.SEA_system as SEAsys
import pyva.models as mds
#import pyva.TMmodel as TMM
import pyva.systems.infiniteLayers as iL
import pyva.geometry.meshClasses as meshC


import random

#@todo define abstract mother class for 2D system for global method implementation 

class Structure2DSystem(SEAsys.SEA_system):
    """
    Class for structural 2D SEA system
    
    This class provide methods for the random description of
    plate SEA systems. Plate systems are described by the wave fields 
    that occur in plates. The wave field is identified by the wave_DOF
    attribute. In classical SEA [Lan1990]_ there are 3 propagating wave types.
    
    1. Longitudinal waves
    2. Shear waves
    3. Bending waves
    
    The wave_DOF correspons to the number in the enumeration.
    As the bending wave equation has two solutions a wave_DOF=4 corresponds
    also to the bending wave.
    
    In [Pei2022]_ it is shown the diffuse field reciprocity is valid for bending
    or a combination of shear- and longitudinal wave. Thus, wave_DOF=5
    identifies in-plane waves. If all waves are not separated and considered 
    in total, this is denoted by wave_DOF = 0.
        
    Attributes
    ----------
    area : float
        area of the plate
    perimeter : float
        perimeter of the plate
    curvature : str
        topology identifier 'flat', 'singlecurved', 'doublycurved' 
    """

    def __init__(self,ID,area,prop, 
                 curvature='flat', 
                 wave_DOF = [3,5], # bending must be first
                 eta = 0.01, 
                 perimeter=0.,
                 trim = (None,None)):
        """
        Constructor for two dimensional structural SEA systems

        Parameters
        ----------
        ID : int
            System identifier.
        area : float
            area of the plate.
        prop : PlateProp
            Property of plate.
        curvature : str, optional
            identifier for curvature. The default is 'flat'.
        wave_DOF : in, optional
            wave degrees of freedom. The default is [3,5].
        # bending must be first
        eta : float or Signal, optional
            damping loss factor. The default is 0.01.
        perimeter : float, optional
            perimeter of plate. The default is 0..
        trim : tuple of TMmodel, optional
            noise control treatement lay-up on front and back side. The default is ('none','none').

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        if wave_DOF[0] == 3 or wave_DOF[0] == 0:
            super().__init__(ID,wave_DOF,eta=eta)
        else:
            raise ValueError('First wavefield must be bending or all')

        self.area = area
        self.prop = prop
        if perimeter == 0.:
            self.perimeter = 4*np.sqrt(area)
        else:
            self.perimeter = perimeter
        self.curvature = curvature
        self.trim_sw = [False,False]
        self.trim = trim
        if len(trim) == 2:
            for i,it in enumerate(trim):
                if isinstance(it,mds.TMmodel):
                    self.trim_sw[i] = True
        
    def __repr__(self):
        _str = 'SEA shell system with ID:{0}\t'.format(self.ID)
        _str += 'reverberant wave_DOF(s):{0}'.format(str(self.wave_DOF.dof)) 
        return _str
    
    @property
    def Lx(self):
        """
        Lx estimator 
        
        This method calculates the best fit of area and perimeter of a rectangular plate
        for the calculation of the radiation efficiency based on the assuption
        that it can be estimated from a rectangular plate of similar shape.
        
        Returns
        -------
        float
            Lx.

        """
        
        if self.perimeter <= 4*np.sqrt(self.area):
            Lx = np.sqrt(self.area)
        else:
            Pq = self.perimeter/4 # P/4
            A  = self.area
            Lx = Pq + np.sqrt(Pq*Pq-A)
        
        return Lx
    
    @property
    def Ly(self):
        """
        Ly estimator 
        
        This method calculates the best fit of area and perimeter of a rectangular plate
        for the calculation of the radiation efficiency based on the assuption
        that it can be estimated from a rectangular plate of similar shape.
        
        Returns
        -------
        float
            Ly.

        """
        
        A  = self.area
        
        Ly = A/self.Lx
        
        return Ly
    
    
        
    def modal_density(self,omega,wave_DOF=3):
        """
        Modal density of plate system
        
        Parameters
        ----------
        omega : float
            angular frequency (not used for flat plates)
        wave_DOF : int
            wavetype 1:longitudinal 2:shear 3/4: bending 5:in_plane 0:all 
        
        Returns
        -------
        float
            modal_density
        
        """
        
        if wave_DOF == 1:
            k    = np.real(self.prop.wavenumber_L(omega))
            c_gr = np.real(self.prop.c_L())
        elif wave_DOF == 2:
            k    = np.real(self.prop.wavenumber_T(omega))
            c_gr = np.real(self.prop.c_T())
        elif wave_DOF in (3,4):
            k    = self.prop.wavenumber_B(omega)
            c_gr = self.prop.c_B_group(omega)
        elif wave_DOF == 5:
            return self.modal_density(omega,1)+self.modal_density(omega,2)
        elif wave_DOF == 0:
            return self.modal_density(omega,5)+self.modal_density(omega,3)
        else:
            raise(ValueError('Unknown wave_DOF argument {0}'.format(wave_DOF)))
        
        return self.area*k/2/np.pi/np.real(c_gr)
        
    def mode_count(self,omega,wave_DOF=3):
        """
        Number of modes until omega 

        Parameters
        ----------
        omega : double
            maximum angular frequency.
        wave_DOF : integer, optional
            wave degree of freedom. The default is 3.

        Returns
        -------
        integer
            number of modes.

        """
        
        if wave_DOF == 1:
            k    = np.real(self.prop.wavenumber_L(omega))
        elif wave_DOF == 2:
            k    = np.real(self.prop.wavenumber_T(omega))
        elif wave_DOF in (3,4):
            k    = self.prop.wavenumber_B(omega)
        elif wave_DOF == 5:
            return self.modal_density(omega,1)+self.modal_density(omega,2)
        elif wave_DOF == 0:
            return self.modal_density(omega,5)+self.modal_density(omega,3)
        else:
            raise(ValueError('Unknown wave_DOF argument {0}'.format(wave_DOF)))
        
        return self.area*k**2/4/np.pi
        
    def damping_loss(self,omega,wave_DOF):
        """
        Damping loss of SEA systemns
        
        Differntiation between wavetype not yet implemented
        
        Parameters
        ----------
        omega : double
            maximum angular frequency.
        wave_DOF : integer, optional
            wave degree of freedom. The default is 3.

        Returns
        -------
        float
            damping_loss of plate
            
        """
        return np.array([self.eta]*len(omega))
    
    @property
    def mass(self):
        """
        Total mass of plate

        Returns
        -------
        float
            plate mass.

        """
        
        return self.area*self.prop.mass_per_area
    
    def isplate(self):
        """
        Confirms that SEA system is a plate

        Returns
        -------
        bool
            True.

        """
        return True
    
    def iscavity(self):
        """
        Confirms that SEA system is not a cavity

        Returns
        -------
        bool
            False.

        """
        return False

        
    def force_excitation_power(self,omega,force = 1):
        """
        Power input due to normal force excitation
        
        This method assumes free field bending wave radation in the plate.
        This is a very strict assumption only fulfilled by large plates and
        high frequency. 

        Parameters
        ----------
        omega : float
            angular frequency.
        force : complex, optional
            force rms value. The default is 1.

        Returns
        -------
        float
            input power.

        """
        return self.prop.force_excitation_power(omega,force)        
               
            
    def physical_unit(self,omega,energy,restype = 'velocity'):
        """
        Provides physical unit / velolcity of plate for energy

        Parameters
        ----------
        omega : ndarray
            angular frequency
        energy : float
            energy of wave field
  
        Returns
        -------
        ndarray 
            physical unit

        """
        
        if restype == "velocity":
            denom = self.prop.mass_per_area*self.area
        elif restype == "displacement":
            denom = self.prop.mass_per_area*self.area*omega**2
        
        return np.sqrt(energy/denom)
        
        
        
    def non_resonant_TMM(self):
        """
        Provides the mass formulation of the 2Dsystem

        Returns
        -------
        TMmodel
            Transfer matrix object of non-resonant lay-up

        """
        
        TM_non_res = [iL.MassLayer(self.prop.thickness,self.prop.material.rho0)]
        if self.trim_sw[0]:
            #TM_non_res.insert(self.trim[0].layers[::-1])
            TM_non_res = self.trim[0].layers[::-1] + TM_non_res
        if self.trim_sw[1]:
            #TM_non_res.append(self.trim[1].layers[:])
            TM_non_res += self.trim[1].layers[:]
        
        return mds.TMmodel(TM_non_res)
        
    def resonant_TMM(self,trim = True):
        """
        Provides the resonant formulation of the 2Dsystem

        Resonant means here under consideration of bending stiffness

        Returns
        -------
        TMmodel
            Transfer matrix object of resonant layer

        """
        
        TM_res = [iL.PlateLayer(self.prop)]
        if trim:
            if self.trim_sw[0]:
                #TM_non_res.insert(self.trim[0].layers[::-1])
                TM_res = self.trim[0].layers[::-1] + TM_res
            if self.trim_sw[1]:
                #TM_non_res.append(self.trim[1].layers[:])
                TM_res += self.trim[1].layers[:]
           
        return mds.TMmodel(TM_res)

        

    def edge_radiation_stiffness_wavenumber(self,omega,wavenumber,wtype=0):
        """
        radiation_stiffness_wavenumber taken from property method
        
        See also
        --------
        :meth:`pyva.properties.structuralPropertyClasses.PlateProp.edge_radiation_stiffness_wavenumber`
        
        """
        return self.prop.edge_radiation_stiffness_wavenumber(omega,wavenumber,wtype)

    def edge_imaginary_radiation_stiffness_wavenumber(self,omega,wavenumber,wtype=1):
        """
        imaginary_radiation_stiffness_wavenumber taken from property method

        See also
        --------
        :meth:`pyva.properties.structuralPropertyClasses.PlateProp.edge_imaginary_radiation_stiffness_wavenumber`

        """
        return self.prop.edge_imaginary_radiation_stiffness_wavenumber(omega,wavenumber,wtype)

        
    def plate_wavenumber(self,omega,i_wave):
        """
        plate_wavenumber taken from property method

        See also
        --------
        :meth:`pyva.properties.structuralPropertyClasses.PlateProp.plate_wavenumber`
        """
        return self.prop.plate_wavenumber(omega,i_wave)

    def wave_transformation_matrix(self,omega,wavenumber,inv=False):
        """
        wave_transform taken from property method

        See also
        --------
        :meth:`pyva.properties.structuralPropertyClasses.PlateProp.wave_transformation_matrix`
        """
        return self.prop.wave_transformation_matrix(omega,wavenumber,inv)           
    
    def edge_wave_amplitude_radiation_stiffness(self,omega,wavenumber,i_wave):
        """
        edge_wave_amplitude_radiation_stiffness taken from property method

        See also
        --------
        :meth:`pyva.properties.structuralPropertyClasses.PlateProp.edge_wave_amplitude_radiation_stiffness`
        """
        
        return self.prop.edge_wave_amplitude_radiation_stiffness(omega,wavenumber,i_wave)     

    def edge_wave_excitation_force(self,omega,wavenumber,i_wave,matrix):
        """
        wave_excitation_force taken from property method
        
        See also
        --------
        :meth:`pyva.properties.structuralPropertyClasses.PlateProp.edge_wave_excitation_force`
        """
        return self.prop.edge_wave_excitation_force(omega,wavenumber,i_wave,matrix)

    def edge_wave_excitation_displacement(self,omega,wavenumber,i_wave):
        """
        wave_excitation_force taken from property method

        See also
        --------
        :meth:`pyva.properties.structuralPropertyClasses.PlateProp.edge_wave_excitation_displacement`
        """
        return self.prop.edge_wave_excitation_displacement(omega,wavenumber,i_wave)
    
    def edge_wave_excitation_force_cross_correlation(self,omega,wavenumber,i_wave,matrix = False):
        """
        wave_excitation_force_cross_correlation taken from property method

        See also
        --------
        :meth:`pyva.properties.structuralPropertyClasses.PlateProp.edge_wave_excitation_force_cross_correlation`
        """
        return self.prop.edge_wave_excitation_force_cross_correlation(omega,wavenumber,i_wave,matrix)


    def edge_wave_amplitude_radiated_power(self,Psi,omega,wavenumber,i_wave):
        """
        wave_amplitude radiated_power taken from property method

        See also
        --------
        :meth:`pyva.properties.structuralPropertyClasses.PlateProp.edge_wave_amplitude_radiated_power`
        """
        return self.prop.edge_wave_amplitude_radiated_power(Psi,omega,wavenumber,i_wave)          

    def radiation_efficiency(self,omega,fluid = mc.Fluid(),Nstep = 90,simple_muGT1=True):
        """
        Radiation efficieny of rectangular plates using leppingtons theory
        
        This method applies half circle integration over all avaiable mode shapes
        at one frequency.

        Parameters
        ----------
        omega : ndarray
            angular freuquency
        fluid : fluid
            fluid material for half space. The default is mc.Fluid().
        Nstep : int, optional
            Number of integration steps for half circle integration. The default is 90.
        simple_muGT1 : bool, optional
            Switch for simplified Leppington apporach. The default is True.

        Returns
        -------
        sigma : float
            Radiation efficiency of plate.

        """

        ar = aR.HalfSpace(fluid)
            
        phi   = np.linspace(0,np.pi/2,Nstep)  
        phi   = phi[1:-1]
        sigma = np.zeros(np.shape(omega)) 
        k_B = self.plate_wavenumber(omega,3)

        sigs = np.zeros(np.shape(phi))
        
        for ifreq,om in enumerate(omega):
            k_y  = k_B[ifreq]*np.cos(phi)
            k_x  = k_B[ifreq]*np.sin(phi)
            sigs = ar.radiation_efficiency_leppington(om, k_x, k_y, self.Lx, self.Ly,simple_muGT1=simple_muGT1) 
                
            sigma[ifreq] = 2/np.pi*integrate.trapz(sigs,phi)
            
        return sigma
            
    def radiation_efficiency_simple(self,omega,fluid = mc.Fluid()):
        """
        Radiation efficieny of rectangular plates using ISO EN 12354-1
        
        Equations in the code are from [1] 
        
        [1] D. Johansson, P. Comnell, Statistical Energy Analysis Software, Master Thesis, 
        Chalmers University, 2010
        
        Parameters
        ----------
        omega : ndarray
            angular freuquency
        fluid : fluid
            fluid material for half space

        Returns
        -------
        float
            radiation efficiency for reververant wave field of plate

        """



        c0 = fluid.c0
        c2 = c0**2
        Lx = self.Lx
        Ly = self.Ly
        
        f   = omega/2/np.pi
        fc  = self.prop.coincidence_frequency(c0)/2/np.pi
        f11 = c2/4/fc*(1/Lx**2+1/Ly**2)
        sigma = np.zeros(np.shape(omega)) 

        sigma1 = 1/np.lib.scimath.sqrt(1-fc/f)
        sigma2 = 4*Lx*Ly*(f/c0)**2
        lam2   = f/fc
        lam    = np.sqrt(lam2)

        if f11<=fc/2:
            # See frequency and index diagram
            #
            #-----f11----fc/2------------fc---------------------------
            # ix1  !   ix2 !   -ix3-     !        ix4
            #     ix12     !
            #      !    ix23             !
            #            ix0             !
            
            ix4  = f >= fc
            #ix1  = f < f11
            ix12 = f <= fc/2
            ix0  = f < fc
            #ix3  = np.logical_and(f < fc,f>=f11)
            #ix2  = np.logical_and(f<=fc/2,f>=f11) 
            
            sigma[ix4] = np.real(sigma1[ix4])
            
            delta1 = ((1-lam2[ix0])*np.log((1+lam[ix0])/(1-lam[ix0]))+2*lam[ix0])/(4*np.pi**2*(1-lam2[ix0])**1.5)
            sigma[ix0] = 2*(Lx+Ly)/(Lx*Ly)*c0/fc*delta1
            # calculate delta2 for f<=fc/2
            delta2     = 8*c2*(1-2*lam2[ix12])/(fc**2*np.pi**4*Lx*Ly*lam[ix12]*np.sqrt(1-lam2[ix12]))
            sigma[ix12] += delta2
            ix_eq4_54 = np.logical_and(f<f11,sigma>sigma2,sigma<=2)
            sigma[ix_eq4_54] = sigma2[ix_eq4_54]
        else:
            #------   ----fc/2-----f11---fc---------------------------
            # ix1  !   ix2 !    ix3      !        ix4
            #      !    ix23             !
            #            ix0             !
            
            sigma3 = np.sqrt(omega*(Lx+Ly)/16/c0)
            
            i_eq_4_56 = np.logical_and(f<fc,sigma2<sigma3)
            i_eq_4_57 = np.logical_and(f>fc,sigma1<sigma3)
            i_else    = np.ones(np.shape(omega),dtype = bool)
            i_else    = np.logical_and(i_else,np.logical_not(np.logical_or(i_eq_4_56,i_eq_4_57)))
            
            sigma[i_eq_4_56] = sigma2[i_eq_4_56] 
            sigma[i_eq_4_57] = sigma1[i_eq_4_57]
            sigma[i_else]    = sigma3[i_else]
            
        sigma[sigma > 2.] = 2 # Eq 4.55 and 4.59
            
        return sigma
    
    def w_random(self,omega,F):
        """
        Single SEA plate subsystem response to rms force

        Parameters
        ----------
        omega : nd.arreay
            angular frequency
        F : float
            rms force

        Returns
        -------
        float
            random rms displacement.

        """
        return F**2/(8.*self.area*np.sqrt(self.prop.B*self.prop.mass_per_area**3)*omega**3*self.prop.material.eta)
    




class RectangularPlate(Structure2DSystem):
    """
    Class for simple rectangular plate methods
    
    The rectangular plate extends the SEA plate description by deterministic 
    properties that follow from the simple rectangular geometry and that allow 
    for modal simulation methods.
    
    Attributes
    ----------
    Lx : float
        length in x-direction
    Ly : float
        length in y-direction
    
    """
    def __init__(self,ID,Lx,Ly, curvature='flat',
                 prop=stPC.PlateProp,
                 wave_DOF = [3,5],
                 eta=0.01,
                 trim = ('none','none')):
        """
        Construtor of rectangular plate

        Parameters
        ----------
        ID : int
            System identifier.
        Lx : float
            length of plate in x-direction.
        Ly : float
            length of plate in y-direction.
        prop : PlateProp
            Property of plate.
        curvature : str, optional
            identifier for curvature. The default is 'flat'.
        wave_DOF : list of int, optional
            wave degrees of freedom. The default is [3,5] bending must be first
        eta : float or Signal, optional
            damping loss factor. The default is 0.01.
        trim : tuple of TMmodel, optional
            noise control treatement lay-up on front and back side. The default is ('none','none').

        Returns
        -------
        None.

        """
        super().__init__(ID,Lx*Ly,prop,wave_DOF = wave_DOF, eta = eta, trim = trim,perimeter = 2*(Lx+Ly))
        self._Lx = Lx
        self._Ly = Ly
        
    @property
    def Lx(self):
        return self._Lx
    
    @property
    def Ly(self):
        return self._Ly
    
        
    def __str__(self):
        _str = 'SEA rectangular plate system with ID:{0}\t'.format(self.ID)
        _str += 'reverberant wave_DOF(s):{0}'.format(str(self.wave_DOF.dof)) 
        return _str
    
    # Methods for application of analytical theory in book, less use for practical simulation
    def w_mode(self,nx,ny,x,y):
        """
        Displacement mode shape for bending waves

        Parameters
        ----------
        nx : int
            mode index in x.
        ny : int
            mode index in y.
        x : float
            x-coordinate.
        y : float
            y-coordinate.

        Returns
        -------
        complex
            displacement.

        """
        
        return 2/np.sqrt(self.area)*np.sin(nx*np.pi*x/self.Lx)*np.sin(ny*np.pi*y/self.Ly)

    def w_mode_map(self,nx,ny,x,y,dX2,dY2):
        """
        Mapped displacement mode shape for bending waves
        
        The mapped mode shape simulates a numeric mode shape by sampling
        the analytical solution be a regular mesh

        Parameters
        ----------
        nx : int
            mode index in x.
        ny : int
            mode index in y.
        x : float
            x-coordinate.
        y : float
            y-coordinate.
        dX2 : float
            half element edge length in x.
        dY2 :  float
            half element edge length in y.

        Returns
        -------
        float or ndarray
            mode shape.

        """
        Lx = self.Lx
        Ly = self.Ly
        dA = 4*dX2*dY2
        nxpi = nx*np.pi/Lx
        nypi = ny*np.pi/Ly
        return 2*np.sqrt(self.area)/np.pi**2/nx/ny/dA\
                *(np.cos(nxpi*(x+dX2))-np.cos(nxpi*(x-dX2)))\
                *(np.cos(nypi*(y+dY2))-np.cos(nypi*(y-dY2)))
        
    def w_modal_force(self,omega,N,Fz,x,y,x0,y0):
        """
        Model frequency response due to force exciation

        Parameters
        ----------
        omega : float
            angular frequency.
        N : tuple of int
            Maximum mode number.
        Fz : complex
            Force amplitude in z.
        x : float
            response x-coordinate.
        y : float
            response y-coordinate.
        x0 : float
            force x-coordinate.
        y0 : float
            force y-coordinate.

        Returns
        -------
        w : complex
            response displacement.

        """
        
        w   = 0.
        mpA = self.prop.mass_per_area
        for nx in range(1,N[0]+1):
            for ny in range(1,N[1]+1):
                omn2 = self.omega_mode(nx,ny)**2
                A    = np.conj(self.w_mode(nx,ny,x0,y0))*Fz/mpA
                w    +=self.w_mode(nx,ny,x,y)*A/(omn2*(1+1j*self.prop.material.eta)-omega*omega)
        return w
        
    def w_modal_force_random(self,omega,N,Fz,x,y,x0,y0,var,seed):
        """
        Model frequency response due to force exciation with random variation
        
        The var parameter varies mass and modal frequency

        Parameters
        ----------
        omega : float
            angular frequency.
        N : int
            Maximum mode number.
        Fz : complex
            Force amplitude in z.
        x : float
            response x-coordinate.
        y : float
            response y-coordinate.
        x0 : float
            force x-coordinate.
        y0 : float
            force y-coordinate.
        var : float
            relative variation of mass and modal frequeny
            The variation interval is [1-var/2 , 1+var/2].
        seed : float
            seed for random functioun.

        Returns
        -------
        w : complex
            mode shape.

        """
        
        random.seed(seed)
        w = 0.

        for nx in range(1,N[0]+1):
            for ny in range(1,N[1]+1):
                omn2 = (self.omega_mode(nx,ny)*random.uniform(1-var/2,1+var/2))**2
                Mn   = self.prop.mass_per_area()*random.uniform(1-var/2,1+var/2)
                A    = np.conj(self.w_mode(nx,ny,x0,y0))*Fz/Mn
                w    +=self.w_mode(nx,ny,x,y)*A/(omn2*(1+1j*self.prop.material.eta)-omega*omega)
        return w
        
    
    
    def omega_mode(self,nx,ny):
        """
        Angular modal frequency of plate mode with index nx, ny

        Parameters
        ----------
        nx : integer > 0 
            mode index in x-direction.
        ny : integer > 0 
            mode index in y-direction.

        Returns
        -------
        float
            Angular modal frequency

        """
    
        return np.sqrt(self.prop.B_per_M)*((nx*np.pi/self.Lx)**2+(ny*np.pi/self.Ly)**2)
    
    def zPlate(self,omega,N,x0,y0):
        """
        Normal impedance 
        
        Calculated with the modal frequency response

        Parameters
        ----------
        omega : float
            angular frequency.
        N : int
            maximum mode number.
        x0 : float
            x-position.
        y0 : float
            y-position.

        Returns
        -------
        complex
            normal point impedance of plate.

        """
        return 1./(1j*omega*self.w_modal_force(omega,N,1.,x0,y0,x0,y0))

    def z_plate_random(self,omega,N,x0,y0,var,seed):
        """
        Normal impedance with variation 
        
        Calculated with the modal frequency response

        Parameters
        ----------
        omega : float
            angular frequency.
        N : int
            maximum mode number.
        x0 : float
            x-position.
        y0 : float
            y-position.
        var : float
            relative variation of mass and modal frequeny
            The variation interval is [1-var/2 , 1+var/2].
        seed : float
            seed for random functioun.

        Returns
        -------
        complex
            normal point impedance of plate.

        """
        
        return 1./(1j*omega*self.w_modal_force_random(omega,N,1.,x0,y0,x0,y0,var,seed))
    
    def get_mesh(self,omega_max,N=6):
        """
        Creates mesh for plate
        
        The element size is calculated from the wavelength of the plate.

        Parameters
        ----------
        omega_max : double
            Maximum angular frequency for determination of smallest wavelength
        N : interger, optional
            Minimal number of nodes per wavelength

        Returns
        -------
        mesh : regmesh2D
            The mesh covering the given conditions

        """
        
        # determine wavelength plate
        lam = self.prop.wavelength_B(omega_max)
        Dmin = lam/N
        NX = int(np.ceil(self.Lx/Dmin+1))
        NY = int(np.ceil(self.Ly/Dmin+1))
    
        mesh = meshC.RegMesh2D(0., 0., self.Lx, self.Ly, NX, NY)
        
        return mesh
    
    def get_modes_index(self,omega_max):
        """
        gets sorted modal frequencies and indexs

        Parameters
        ----------
        omega_max : float
            Maximum angular frequency.

        Returns
        -------
        oms : ndarray float
            modal angular frequencies.
        ns : ndarray in int
            [Nmode x 2] array with sorted mode index pairs.

        """
        
        
        # determine modes up to omega_max
        bpm = self.prop.B_per_M
        Lx = self.Lx
        Ly = self.Ly

        # modal wavenumber steps
        D_kx = np.pi/Lx
        #D_ky = np.pi/Ly
        
        # find max wavenumber corresponding to omega max
        k_max = self.prop.wavenumber_B(omega_max)
        
        # initialize node index array
        Nmax = int(np.round(self.mode_count(omega_max))) # two for some safety
        ns  = np.zeros((Nmax,2),dtype=np.uint32)
        oms = np.zeros((Nmax,))

        # start with max Nx
        #Nx1 = int(np.floor(k_max/D_kx))
        Nx = int(np.floor(np.sqrt(omega_max*(Lx/np.pi)**2/np.sqrt(bpm)-(Lx/Ly)**2)))
        
        
        n0 = 0
        
        for ix in range(Nx):
            # current omega y
            nx = ix+1
            if k_max  > nx*D_kx:
                Ny = int(np.floor(np.sqrt(omega_max*(Ly/np.pi)**2/np.sqrt(bpm)-(nx*Ly/Lx)**2)))
                #Ny       = int(np.ceil(ky_max/D_ky)) 
                iy = np.arange(Ny)
            else:
                iy = 0
            
            
            oms[iy+n0] = self.omega_mode(nx,iy+1)
            ns[iy+n0,0] = nx 
            ns[iy+n0,1] = iy+1
            n0 += Ny
            
        # remove remaining zeros
        ix = oms != 0
        oms = oms[ix]
        ns  = ns[ix,:]
        # sort modes in frequency orde
        ix = np.argsort(oms)
        oms = oms[ix]
        ns  = ns[ix,:]
                    
        return oms,ns

    def normal_modes(self,omega_max,norm = 'mass',mapping='element',N=4):
        """
        Normal mode shapes mapped from analytic solution

        Parameters
        ----------
        omega_max : float
            Maximum modal angular frequency.
        norm : str, optional
            Identifier for normalisation. The default is 'mass'.
        mapping : str, optional
            Identifier for mapping. The default is 'element'.
        N : int, optional
            Number of samplex per wavlength. The default is 4.

        Returns
        -------
        modes : shape2D
            modes.
        mesh : mesh2D
            mesh for node location.

        """
        
        # get frequuencies and indexes
        omega_n,ns = self.get_modes_index(omega_max)
        xdata      = mC.DataAxis(omega_n,typestr = 'angular frequency')
        
        # create sppropriate mesh
        mesh       = self.get_mesh(omega_max,N=N)
        
        # .. and according meshSignal
        ydata = np.zeros((mesh.Nmesh,len(omega_n)))
        x,y = mesh.nodes()
        
        modedof = dof.DOF(np.arange(1,mesh.Nmesh+1,dtype=np.int64),[3],dof.DOFtype(typestr = 'mass', exponent = -0.5),repetition = True)

        X,Y = mesh.nodes()
        dX2 = mesh.dX/2
        dY2 = mesh.dY/2


                
        for ix,nn in enumerate(ns):
            if mapping == 'element':
                ydata[:,ix] = np.sqrt(self.area/self.mass)*self.w_mode_map(nn[0],nn[1],X,Y,dX2,dY2).flatten()
            elif mapping == 'mesh':
                ydata[:,ix] = np.sqrt(self.area/self.mass)*self.w_mode(nn[0],nn[1],X,Y).flatten()
            
        modes = mC.ShapeSignal(mesh,xdata,ydata,modedof)
        
        return (modes,mesh)
        
    def dynamic_stiffness_mesh(self,omega,mesh,method='infinite'):
        """
        Dynamic stiffness or rectangular plate assuming infinite dimentions
        
        Parameters
        ----------
        omega : float
            angular frequency.
        mesh : mesh2D
            radiating mesh
        method : TYPE, optional
            DESCRIPTION. The default is 'infinite'.

        Returns
        -------
        values : DynamicMatrix
            Dynamic stiffness matrix of infinite plate
        """
        
                
        
        # Calculate distances and index to repetitive distances
        distances, index = mesh.distance()
        
        Ndist  = len(distances)
        Nx     = len(omega)
        Nmesh  = mesh.Nmesh

        def Rx(x,omega):
            B = self.prop.B_complex
            k = self.prop.wavenumber_B(omega,comp = True)
            
            if x > 0:
                values = 1/(8j*B*k*k)*(special.hankel2(0,k*x)-special.hankel2(0,-1j*k*x))
            else:
                values = -1j/(8*B*k*k)
            return values
        

        #_res  = np.zeros(((Nmesh*(Nmesh+1))//2,Nx),dtype=np.complex128)
        _res  = np.zeros((Nmesh,Nmesh,Nx),dtype=np.complex128)
        _Rbuf = np.zeros(Ndist,dtype=np.complex128) 
        _A    = mC.LinearMatrix(np.zeros(((Nmesh*(Nmesh+1))//2,1),dtype=np.complex128),sym=1,shape=(Nmesh,Nmesh,1))        

        for i in range(Nx):
            for id in range(Ndist):
                _Rbuf[id] = Rx(distances[id],omega[i])
                
            _A._data[:,0] = _Rbuf[index]
            _res[:,:,i] = linalg.inv(_A.Dindex(0)) # triu indexing required!!!!

               
        return mC.LinearMatrix(_res)#,sym = 1,shape=(Nmesh,Nmesh,Nx) )
        
        
    def transmission_coefficient_discrete(self,omega,half_spaces,mesh=None,method='piston'):
        """
        method for transmission coeffient calculation based on discrete radiation stiffness
        
        This method is based om the the dicrete plate stiffness of infinite plates
        and not on the plate modes

        Parameters
        ----------
        omega : float 
            angular frequencies
        half_spaces : tuple of HalfSpace
            connectet HalfSpaces, two spaces when fluid differs 
        mesh : mesh2D
            mesh for discrete simulation, default is None for automatic calculation
        method : str
            Identifier for calculation method 'piston' or 'wavelett'

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        tau : float
            transmission coefficient.        
        """
        
        # Determine number of half spaces
        
        if isinstance(half_spaces,(tuple,list)) and len(half_spaces) in [1,2]:
            if len(half_spaces)==1:
                    
                hs1 = half_spaces[0]
                hs2 = hs1

                if mesh is None:
                    mesh = hs1.get_mesh(np.max(omega),self.Lx,self.Ly,N=2)

                # D1 = hs1.radiation_stiffness_mesh(omega,mesh,method)
                # D2 = D1
                two_fluids = False
            else:
                hs1 = half_spaces[0]
                hs2 = half_spaces[1]

                if mesh is None:
                    # determine smallest wavelength
                    lam1 =  hs1.fluid.wavelength(np.max(omega),self.Lx,self.Ly)
                    lam2 =  hs2.fluid.wavelength(np.max(omega),self.Lx,self.Ly)
                    
                    if lam1 < lam2:
                        mesh = hs1.get_mesh(np.max(omega),N=2)
                    else:
                        mesh = hs2.get_mesh(np.max(omega),N=2)

                two_fluids = True

            print('Method used mesh with {0} elements\n'.format(mesh.Nmesh))
                # D1 = hs1.radiation_stiffness_mesh(omega,mesh,method)
                # D2 = hs2.radiation_stiffness_mesh(omega,mesh,method)
        else:
            raise ValueError('half_spaces must be tuple of length 1 or 2 of HalfSpace objects')
                     
        tau = np.zeros(omega.shape)
       
        
        for ifreq,om_ in enumerate(omega):
            # Create Dtot
            D1 = hs1.radiation_stiffness_mesh([om_],mesh,method)
            if two_fluids:
                D2 = hs2.radiation_stiffness_mesh([om_],mesh,method)
            else:
                D2 = D1
                
            Dtot = (D1 + D2 + self.dynamic_stiffness_mesh([om_],mesh)).Dindex(0)

            print('Dealing with omega = {0:.1f}'.format(omega[ifreq]))
            D1_  = np.imag(D1.Dindex(0))
            Dtot_= linalg.inv(Dtot)
            D2_  = np.imag(D2.Dindex(0))
            
            Dtot__ = np.dot(Dtot_,D2_)
            Dtot__ = np.dot(Dtot__,np.conjugate(np.transpose(Dtot_)))
            tau[ifreq] = 16*np.pi/(hs2.fluid.wavenumber(omega[ifreq])**2*self.area)*(D1_*Dtot__).sum() 
            
        return tau
    
    def modal_transmission_coefficient_discrete(self,omega,half_spaces,trim=(None,None),method='piston',modal_factor = 1.2,N=4,mapping='mesh'):
        """
        modal transmission coeffient based on discrete radiation stiffness
                        
        Parameters
        ----------
        omega : float 
            angular frequencies
        half_spaces : tuple of HalfSpace
            connectet HalfSpaces, two spaces when fluid differs 
        trim : TYPE, optional
            DESCRIPTION. The default is (None,None).
        method : str
            Identifier for calculation method 'piston' or 'wavelett'
        modal_factor : float, optional
            factor for upper frequency of used modes. The default is 1.2.
        N : int, optional
            nodes per wavelength. The default is 4.
        mapping : str, optional
            Identifier for mapping method. The default is 'mesh'.
            
        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        tau : ndarray of float
            transmission coefficient of the MassLayer
        """
        
        # Get modal base 
        modes,mesh = self.normal_modes(omega.max()*modal_factor,N=N,mapping=mapping)

        
        # Determine number of half spaces
        if isinstance(half_spaces,(tuple,list)) and len(half_spaces) in [1,2]:
            N_hs = len(half_spaces)
            if N_hs==1:
                hs1 = half_spaces[0]
                hs2 = hs1
            else:
                hs1 = half_spaces[0]
                hs2 = half_spaces[1]
        else:
            raise ValueError('half_spaces must be tuple of length 1 or 2 of HalfSpace objects')
        
        
        # arrange modes on modal matrix
        MM    = modes.ydata
        om_n  = modes.xdata.angular_frequency
        om_n2 = om_n**2*(1+1j*self.damping_loss(om_n,wave_DOF=3))
        
        tau  = np.zeros(np.shape(omega))
        Ddir = 2*[None]

        trim_list_sw = [False]*2

        # Deal with trim
        if any([x != None for x in trim]):
            D11 = 2*[None]
            D12 = 2*[None]
            D22 = 2*[None]
            Dred = 2*[None]
            D_F2_F1= 2*[None]

            
        # frequency loop
        print('Dealing with omega =           ',end='')
        for i_om,om in enumerate(omega):
            print('\b\b\b\b\b\b\b\b\b\b\b\b\b {0:>8.1f} Hz'.format(om),end='')

            #convert to modal space
            D1 = hs1.radiation_stiffness_mesh([om],mesh,method)            
            D1_ = D1.Dindex(0)
            Ddir[0] = np.conj(MM.T).dot(D1_).dot(MM)
            # from now D1_ and D2_ are the inverse
 
            if N_hs == 2:
                D2 = hs2.radiation_stiffness_mesh([om],mesh,method)
                D2_ = D2.Dindex(0)
                Ddir[1] = np.conj(MM.T).dot(D2_).dot(MM)
            else:
                Ddir[1] = Ddir[0]
                
            HH = np.diag(om_n2-om**2)
                
            for ix in range(2):
                if isinstance(trim[ix],mds.TMmodel):
                    trim_list_sw[ix] = True
                    D11_,D12_,D22_ = trim[ix].stiffness_matrix_mesh([om],mesh)
                    D11[ix] = np.conj(MM.T).dot(D11_.Dindex(0)).dot(MM)
                    D12[ix] = np.conj(MM.T).dot(D12_.Dindex(0)).dot(MM)
                    D22[ix] = np.conj(MM.T).dot(D22_.Dindex(0)).dot(MM)
                    # 
                    D_F2_F1[ix] = D12[ix].dot(np.linalg.inv(D22[ix]+Ddir[ix]))
                    Dred[ix] = D_F2_F1[ix] @ Ddir[ix] @ D_F2_F1[ix].conj().T 
                    HH+= D_F2_F1[ix]

                else:
                    HH += Ddir[ix]
            
            D_tot_inv = np.linalg.inv(HH)
 
            if all([x == None for x in trim]):
                t_ = (np.imag(Ddir[1])*(D_tot_inv.dot(np.imag(Ddir[0])).dot(np.conj(D_tot_inv.T))) ).sum()
            elif trim[0]==None and trim_list_sw[1]: 
                t_ = (np.imag(Dred[1])*(D_tot_inv.dot(np.imag(Ddir[0])).dot(np.conj(D_tot_inv.T))) ).sum()
            elif trim[1]==None and trim_list_sw[0]: 
                t_ = (np.imag(Ddir[1])*(D_tot_inv.dot(np.imag(Dred[0])).dot(np.conj(D_tot_inv.T))) ).sum()
            else:
                t_ = (np.imag(Dred[1])*(D_tot_inv.dot(np.imag(Dred[0])).dot(np.conj(D_tot_inv.T))) ).sum()
                
        
            ka = np.real(hs1.fluid.wavenumber(om))**2
            tau[i_om] = 16*np.pi/ka/self.area*np.real(t_)
        print('\nfinished')    
        return tau
    



class MassLayer:
    """
    The mass Layer class provides methods in the system context
    
    The implementation here focusses on the 2D or plate like properties
    """
    
    def __init__(self,Lx,Ly,mass_per_area):
        """
        Constructor of mass layer

        Parameters
        ----------
        Lx : float
            length in x-direction.
        Ly : float
            length in y-direction.
        mass_per_area : float
            area weight.

        Returns
        -------
        None.

        """
        self.Lx = Lx
        self.Ly = Ly
        self.mass_per_area = mass_per_area

    @property    
    def area(self):
        """
        property method for area

        Returns
        -------
        float
            area.

        """
        
        return self.Lx*self.Ly
        
    def dynamic_stiffness_mesh(self,omega,mesh):
        """
        Dynamic matrix of mass layer for regular mesh
        

        Parameters
        ----------
        omega : float
            angular frequency.
        mesh : mesh2D
            mesh for matrix generation.

        Returns
        -------
        DynamicMatrix
            Dynammic stiffness matrix of mass layer.

        
        """
                
        Nx     = len(omega)
        Nmesh  = mesh.Nmesh
        
        _res    = np.zeros((Nmesh,Nx),dtype=np.complex128)
        
        for i in range(Nx):
            _res[:,i] = -self.mass_per_area*mesh.dA*omega[i]**2
            
        return mC.LinearMatrix(_res,sym = 3,shape=(Nmesh,Nmesh,Nx) )
        
    def transmission_coefficient_discrete(self,omega,half_spaces,mesh,method='piston'):
        """
        transmission coeffient based on discrete radiation stiffness
        
        Parameters
        ----------
        omega : float
            angular frequency.
        half_spaces : tuple of HalfSpace
            Fluid half spaces.
        mesh : mesh2D
            Mesh for discretisation.
        method : str, optional
            Identifier for radiation method. The default is 'piston'.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        tau : float
            transmission coefficient of the MassLayer
        """
        
        # Determine number of half spaces
        
        if isinstance(half_spaces,(tuple,list)) and len(half_spaces) in [1,2]:
            if len(half_spaces)==1:
                hs1 = half_spaces[0]
                hs2 = hs1
                D1 = hs1.radiation_stiffness_mesh(omega,mesh,method)
                D2 = D1
            else:
                hs1 = half_spaces[0]
                hs2 = half_spaces[1]
                D1 = hs1.radiation_stiffness_mesh(omega,mesh,method)
                D2 = hs2.radiation_stiffness_mesh(omega,mesh,method)
        else:
            raise ValueError('half_spaces must be tuple of length 1 or 2 of HalfSpace objects')
             
        # Create Dtot
        Dtot = D1 + D2 + self.dynamic_stiffness_mesh(omega,mesh)
        #Dtot.inv
        
        tau = np.zeros(omega.shape)
       
        
        for ifreq in range(len(omega)):
            D1_  = np.imag(D1.Dindex(ifreq))
            Dtot_= linalg.inv(Dtot.Dindex(ifreq))
            D2_  = np.imag(D2.Dindex(ifreq))
            
            Dtot__ = np.dot(Dtot_,D2_)
            Dtot__ = np.dot(Dtot__,np.conjugate(np.transpose(Dtot_)))
            tau[ifreq] = 16*np.pi/(hs2.fluid.wavenumber(omega[ifreq])**2*self.area)*(D1_*Dtot__).sum() 
            
        return tau
        
    def transmission_coefficient_angular(self,omega,fluid1,theta=0,fluid2='none'):
        """
        angular transmission coeffient based on wave transmission of infinite plate
        
        Parameters
        ----------
        omega : float
            angular frequency.
        fluid1 : fluid
            Fluid on excitation and radiation side.
        theta : float, optional
            Angle of incidence. The default is 0.
        fluid2 : fluid, optional
            Radiatinh fluid. The default is 'none'.

        Returns
        -------
        float
            transmission coefficient.
        """
        
        # Determine number of half spaces
        
        return 1/(1 + (self.mass_per_area*omega/2/fluid1.z0*np.cos(theta))**2)
    
    def transmission_coefficient_diffuse(self,omega,fluid1,fluid2='none',theta_max=np.pi/2*0.99):
        """
        diffuse transmission coeffient calculation based on discrete radiation stiffness
                        
        Parameters
        ----------
        omega : float
            angular frequency.
        fluid1 : fluid
            Fluid on excitation and radiation side.
        theta_max : float, optional
            Angle of incidence. The default is <pi/2.
        fluid2 : fluid, optional
            Radiating fluid. The default is 'none'.

        Returns
        -------
            angular transmission coefficient of the MassLayer
        """
        
        # Determine number of half spaces
        
        tau = np.zeros(omega.shape)
        
        # Get average denominator
        den = 0.5*np.sin(theta_max)**2
        
        
        for ifreq in range(len(omega)):
            f = lambda x: self.transmission_coefficient_angular(omega[ifreq],fluid1,x)*np.cos(x)*np.sin(x)/den
            tau[ifreq],prec = integrate.quad(f, 0, theta_max)

        return tau    
        
