# -*- coding: utf-8 -*-
"""
Module with abstract SEA base class

The SEA_system module contains the main and abstract SEA class that
forms the base for all datailed SEA system implementatin as plates, shell and cavities
"""

import numpy as np
import pyva.useful as uF
import pyva.data.dof as dof
from abc import ABC

#@todo define abstract mother class for 2D system for global method implementation 

class SEA_system: #(ABC):
    """
    class for SEA systems
    
    The wave_DOF are the correspondant of the local degrees of freedom 
    for deterministic quantities. The have a different nomentclature.
    
    0. Scalar DOF e.g. pressure or all waves
    1. Longitudinal wave
    2. Shear wave
    3. Bending wave
    4. also bending wave (equal to 3)
    5. in-plane wave (1.+2.)
        
    Attributes
    ----------
    wave_DOF: DOF
        wave degrees of freedom of the system
    eta: float or Signal
        Damping loss        

    """

    def __init__(self,ID,wave_DOF,typestr = 'velocity', eta=0.01):
        """
        Class constructor of SEA_system

        Parameters
        ----------
        ID : int
            Syste ID.
        wave_DOF : int
            local wave degree of freedom.
        typestr : str, optional
            Identifier for phyical quantity. The default is 'velocity'.
        eta : float, optional
            Damping loss. The default is 0.01.

        Returns
        -------
        None.

        """
        
        self.wave_DOF = dof.DOF(ID,wave_DOF,dof.DOFtype(typestr=typestr),repetition = True)
        self.eta = eta
        
    def __str__(self):
        _str = 'SEA_system with ID:{0}\t'.format(self.ID)
        _str += 'reverberant wave_DOF(s):{0}'.format(str(self.wave_DOF.dof)) 
        return _str

        
    #@abstractmethod 
    def modal_density(self,omega,wave_DOF=0):
        """
        Modal density of SEA systemns
        
        Parameters
        ----------
        omega : float
            angular frequency 
        wave_DOF: int
            wave degree of freedom 
        
        """
        pass
        
    def modes_in_band(self,omega,wave_DOF = 0, btype = 'oct'):
        """
        Modes in band of cavity system
        
        Parameters
        ----------
        omega : float
            angular frequency
        btype: str
            type of band 'oct' for factored steps and 'lin' for linear steps
        
        Returns
        -------
            modes in band
        
        """
        
        if uF.isscalar(omega):
            raise ValueError('omega must be an array of at least 2 elements')
    
        
        #freq = omega/2/np.pi
    
        if btype == 'oct':
            d  = np.sqrt(omega[1]/omega[0])
            return omega*(d-1/d)*self.modal_density(omega,wave_DOF = wave_DOF)
        elif btype == 'lin':
            d  = omega[1]-omega[0] 
            return d*self.modal_density(omega,wave_DOF = wave_DOF)
        else:
            raise ValueError('Unknown option for btype')
            
    def modal_overlap(self,omega,wave_DOF = 0):
        """
        Modal overlap 
        
        Absorption area and fluid damping is considered
        
        Parameters
        ----------
        omega : ndarray
            angular frequency vector
                
        Returns
        -------
        Modal overlap
        """
        
        return omega*self.damping_loss(omega,wave_DOF=wave_DOF)*self.modal_density(omega, wave_DOF=wave_DOF)
        
        
    #@abstractmethod
    def damping_loss(self,omega,wave_DOF = 0):
        """
        Damping loss of SEA systemns
        
        Parameters
        ----------
        omega : float
            angular frequency 
        wave_DOF: int
            wavetype ID. The default is 0.
        
        """

        
    #@abstractmethod
    def physical_unit(self,omega,energy,restype = 'velocity'):
        """
        Physical unit rms value
        
        Parameters
        ----------
        omega : float
            angular frequency 
        energy : float
            energy of system
        restype : str
            identifier for physical unit
        
        """
    
    @property
    def N_wave_fields(self):
        """
        Number of wave field in physical SEA system
        
        Some physical SEA systems are constituted by several wave fields that
        can be considered as single SEA systems or reverberant fields.
        This method shall provide the number of wave field that are used

        Returns
        -------
        int
            Number of wave fields.

        """
        
        return len(self.wave_DOF)
        
    @property
    def ID(self):
        """
        unique ID of SEA system
        
        Returns
        -------
        int
            ID of SEA system
        """
        
        return self.wave_DOF._ID[0] # must be single DOF
    
    # set of default meshod 
    def isplate(self):
        """
        Checks if system is a plate

        Sets default to False, so only one method must be implemented in 
        daugther class

        Returns
        -------
        bool
            True if system is a plate (must be overloaded).

        """
        return False

    def iscavity(self):
        """
        Checks if system is a plate
        
        Sets default to False, so only one method must be implemented in 
        daugther class
                
        Returns
        -------
        bool
            True if system is a cavity (must be overloaded).

        """
        return False
    
    def isSIF(self):
        """
        Checks if system is a semi infinite fluid
        
        Sets default to False, so only one method must be implemented in 
        daugther class
                
        Returns
        -------
        bool
            True if system is a SIF (must be overloaded).

        """
        return False

    
    
        
        
        



