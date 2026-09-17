"""
Module for condensed / lumped systens

The lumpedSystems package defines the simple mechanical systems as masses, springs or 
harmonic oscillators
"""
import numpy as np

class HarmonicOscillator:
    """
    The harmonocOscillator class defines the dynamics of a simple damped or undamped 
    harmonic oscillator.
    
    The main purpose of this class is the presentation of the oscillator dynamics
    in [Pei2022]_
    
    Attributes
    ----------
    mass : float
        oscillating mass
    stiffness: float
        stiffness of spring
    damping: float
        viscous damping
    """
    
    def __init__(self,mass,stiffness,damping=0.):
        """
        Class contructor for HarmonicOscillator
        
        Parameters
        ----------
        mass : float
            oscillating mass
        stiffness : float
            stiffness of spring
        damping : float
            viscous damping
                
        """
        
        self.mass      = mass
        self.stiffness = stiffness
        self.damping   = damping
    
    # Simple properties from combination of CrossSection and IsoMat

    @property
    def omega_mode(self):
        """
        Resonance frequency or harmonic oscillator
        
        Returns
        -------
        float
            angular resonance frequency without damping
        
        """
        
        return np.sqrt(self.stiffness/self.mass)

    @property
    def omega_resonance(self):
        """
        Resonance frequency or harmonic oscillator
        
        Returns: float
            angular resonance frequency *with* damping
        
        """
        
        zeta = self.critical_damping_ratio
        
        if zeta < 1.:
            return self.omega_mode*np.sqrt(1.-self.critical_damping_ratio**2)
        else:
            raise ValueError("Oscillator is not under-damped")

    @property
    def critical_damping_ratio(self):
        """
        Ratio of viscous- to critical-viscous damping 
        
        Returns
        -------
        float
            The ratio of damping to critical viscous damping
        
        """
        return self.damping/self.critical_viscous_damping

    @property
    def critical_viscous_damping(self):
        """
        The critical viscous damping 
        
        Returns
        -------
        float
            Critical viscous damping
        
        """
        return 2*np.sqrt(self.mass*self.stiffness)
    
    @property
    def damping_loss(self):
        """
        The damping loss 
        
        Returns
        -------
        float
            damping loss
        
        """
        return 2*self.damping/self.critical_viscous_damping
    

    def displacement_amplitude(self,time,x0,v0=0):
        """
        Complex amplitude due to initial conditions
        
        This method separates between the tree cases. It provides the complex 
        oscillation amplitude X(t) for underdamped case. And the real displacement
        over time for the damped cases
        
        Parameters
        ----------
        time : float
            time > 0
        x0 : float
            initial displacement at time=0
        v0 : float
            initial velocity     at time=0
        
        Returns
        -------
        complex 
            displacement amplitude over time
        
        """
        
        zeta = self.critical_damping_ratio
        
        if zeta < 1.: # under damped
            om0 = self.omega_mode
            omD = self.omega_resonance
            Xamp = np.sqrt((x0*omD)**2+(v0+zeta*om0*x0)**2)/omD
            phi  = np.arctan((v0+zeta*om0*x0)/(x0*omD))
            return Xamp*np.exp(-1j*phi-zeta*om0*time)
        elif zeta > 1.: # overdamped
            om0 = self.omega_mode
            sqr   = np.sqrt(zeta*zeta-1)
            
            B1 =  (x0*om0*(zeta+sqr) + v0)/2/om0/sqr
            B2 = -(x0*om0*(zeta-sqr) + v0)/2/om0/sqr
            return B1*np.exp((-zeta+sqr)*om0*time)+B2*np.exp(-(zeta+sqr)*om0*time)
        else: # critically damped 
            om0 = self.omega_mode
            B3 = x0
            B4 = v0 + om0*x0
            return (B3+B4*time)*np.exp(-om0*time)

    def velocity_amplitude(self,time,x0,v0=0):
        """
        Complex velocity amplitude due to initial conditions
        
        This method separates between the tree cases. It provides the complex 
        oscillation velocity V(t) for underdamped case. And the real velocity
        over time for the damped cases
        
        Parameters
        ----------
        time : float
            time > 0
        x0 : float
            initial displacement at time=0
        v0 : float
            initial velocity     at time=0

        Returns
        -------
        complex
            velocity amplitude over time
        
        """
        
        zeta = self.critical_damping_ratio
        
        if zeta < 1.:
            om0 = self.omega_mode
            omD = self.omega_resonance
            Xamp = np.sqrt((x0*omD)**2+(v0+zeta*om0*x0)**2)/omD
            phi  = np.arctan((v0+zeta*om0*x0)/(x0*omD))
            return (1j*omD-zeta*om0)*Xamp*np.exp(-1j*phi-zeta*om0*time)
        elif zeta > 1.:         
            om0 = self.omega_mode
            omD = self.omega_resonance
            sqr   = np.sqrt(zeta*zeta-1)
            
            B1 =  (x0*om0*(zeta+sqr) + v0)/2/om0/sqr*(-zeta+sqr)*om0 # incl derivative
            B2 = -(x0*om0*(zeta-sqr) - v0)/2/om0/sqr*(-zeta-sqr)*om0
            return B1*np.exp((-zeta+sqr)*om0*time)+B2*np.exp((-zeta-sqr)*om0*time)
        else:
            om0 = self.omega_mode
            B3 = x0
            B4 = v0 + om0*x0
            return (-om0*(B3+B4*time)+B4)*np.exp(-om0*time)

    def displacement(self,time,x0,v0=0):
        """
        Displacement of harmonic oscillator over time
        
        Parameters
        ----------
        time : float
            time > 0
        x0 : float
            initial displacement at time=0
        v0 : float
            initial velocity     at time=0

        Returns
        -------
        float
            displacement over time
        
        """
        
        zeta = self.critical_damping_ratio
        X0   = self.displacement_amplitude(time,x0,v0)
        
        
        if  zeta < 1.:
            omD  = self.omega_resonance
            return np.real(X0*np.exp(1j*omD*time))            
        else:
            return X0
            
    def velocity(self,time,x0,v0=0):
        """
        Velocity of harmonic oscillator over time
        
        Parameters
        ----------
        time : float
            time > 0
        x0 : float
            initial displacement at time=0
        v0 : float
            initial velocity     at time=0

        Returns
        -------
        float
            displacement over time
        
        """
        zeta = self.critical_damping_ratio
        V0   = self.velocity_amplitude(time,x0,v0)
        
        if  zeta < 1.:
            omD  = self.omega_resonance
            return np.real(V0*np.exp(1j*omD*time))            
        else:
            return V0
        
            
    def kinetic_energy(self,time,x0,v0=0):
        """
        kinetic energy of harmonic oscillator over time
        
        Parameters
        ----------
        time : float
            time > 0
        x0 : float
            initial displacement at time=0
        v0 : float
            initial velocity     at time=0

        Returns
        -------
        float
            kinetic energy over time
        
        """
    
        return 0.5*self.mass*self.velocity(time,x0,v0)**2
        
    def potential_energy(self,time,x0,v0=0):
        """
        potential energy of harmonic oscillator over time
        
        Parameters
        ----------
        time : float
            time > 0
        x0 : float
            initial displacement at time=0
        v0 : float
            initial velocity     at time=0

        Returns
        -------
        float
            potential energy over time
        
        """

        return 0.5*self.stiffness*self.displacement(time,x0,v0)**2
        
    def energy(self,time,x0,v0=0):
        """
        total energy of harmonic oscillator over time
        
        Parameters
        ----------
        time : float
            time > 0
        x0 : float
            initial displacement at time=0
        v0 : float
            initial velocity     at time=0

        Returns
        -------
        float
            energy over time
        
        """
        
        return self.kinetic_energy(time,x0,v0)+self.potential_energy(time,x0,v0)
    
    def u_force(self,omega,F):
        """
        Harmonic forced response of harmonic oscillator
        
        Parameters
        ----------
        omega : float
            angular frequency
        F : complex
            amplitude of force excitation
            
        Returns
        -------
        complex
            displacement amplitude
        """
        
        return F/(self.stiffness-self.mass*omega**2+1j*self.damping*omega)
        

        
            
                    