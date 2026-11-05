"""
Module for one-dimensional structural subsystems

In most cases this concerns beams.
"""

import numpy as np


class beam:
    """
    The beam class deals with dynamic motion of beams.
    
    Attributes
    ----------
    L : float
        length of beam
    beam_prop : BeamProp
        property of beam
    
    """
    def __init__(self,L,beam_prop):
        """
        Constructor for beam

        Parameters
        ----------
        L : float
            length of beam
        beam_prop : BeamProp
            property of beam.

        Returns
        -------
        None.

        """
        
        self.L = L
        self.beam_prop = beam_prop
    
    def __str__(self):
        """
        String output of beam

        Returns
        -------
        str

        """
        _str  = "beam: \n"
        _str += "L               : {0}\n".format(self.L)
        _str += "beam_prop:\n{0}\n".format(self.beam_prop)
        
        return _str
        
    def __repr__(self):
        """
        Repr of beam

        Returns
        -------
        str.

        """
        return "beam({0},{1})".format(self.L,self.beam_prop.__repr__())
    
 
    # Simple properties from combination of CrossSection and IsoMat

    def w_mode(self,x,n):
        """
        Displacement mode shape

        Parameters
        ----------
        x : float
            x-coordinate along length.
        n : int
            mode number.

        Returns
        -------
        float
            displacement of mode shape.

        """
        
        L = self.L
        return np.sqrt(2/L)*np.sin(n*np.pi/L*x)

    def Rx_mode(self,x,n):
        """ 
        Rotation of mode shape

        Parameters
        ----------
        x : float
            x-coordinate along length.
        n : int
            mode number.

        Returns
        -------
        float
            rotation of mode shape.

        """
        
        L = self.L
        return np.sqrt(2/L)*np.cos(n*np.pi/L*x)

    def omega_mode(self,n):
        """
        Angular modal frequency

        Parameters
        ----------
        n : int
            mode number.

        Returns
        -------
        float
            angular modal frequency.

        """
        
        BperM = self.beam_prop.Bx/self.beam_prop.mass_per_length
        return (n*np.pi/self.L)**2*np.sqrt(BperM)

    def w_modal_force(self,omega,x,N,F,x0):
        """
        Modal displacement response due to force excitation

        Parameters
        ----------
        omega : float
            angular frequency.
        x : float
            x-coordinate along length.
        N : int
            Maximum mode number.
        F : complex
            amplitude of force exciation.
        x0 : float
            x-coordinat of exciation.

        Returns
        -------
        w : complex
            displacement response.

        """
        w = 0.
        eta = self.beam_prop.iso_mat.eta
        for im in range(1,N+1):
            #print(im)
            omn2 = self.omega_mode(im)**2
            A=np.conj(self.w_mode(x0,im))*F/self.beam_prop.mass_per_length
            w += self.w_mode(x,im)*A/(omn2*(1+1j*eta)-omega**2)
        return w

    def w_modal_moment(self,omega,x,N,M,x0):
        """
        Modal displacement response due to moment excitation

        Parameters
        ----------
        omega : float
            angular frequency.
        x : float
            x-coordinate along length.
        N : int
            Maximum mode number.
        F : complex
            amplitude of force exciation.
        x0 : float
            x-coordinat of exciation.

        Returns
        -------
        w : complex
            displacement response.

        """        
        
        w = 0.
        eta = self.beam_prop.iso_mat.eta
        for im in range(1,N+1):
            #print(im)
            omn2 = self.omega_mode(im)**2
            A = np.conj(self.Rx_mode(x0,im))*M/self.beam_prop.mass_per_length
            w += self.Rx_mode(x,im)/(omn2**2*(1+1j*eta)-omega**2)*A
        return w
    
    def z_beam(self,omega,x0,N):
        """
        Mechanical point impedance of beam

        Parameters
        ----------
        omega : float
            angular frequency.
        x0 : float
            x-coordinat of exciation.
        N : int
            Maximum mode number.

        Returns
        -------
        complex
            point impedance.

        """
        return 1./(1j*omega*self.w_modal_force(omega,x0,N,1.,x0))    
    
    


