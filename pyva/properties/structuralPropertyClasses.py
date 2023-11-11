"""
Module for structural physical properties.
"""

import numpy as np
import pyva.properties.geometricalPropertyClasses as geoPC
import pyva.properties.materialClasses as mc
import pyva.data.matrixClasses as mC
import scipy.special as scs
import scipy.integrate as integrate

class BeamProp:
    """
    Class for beams.
    
    The BeamProp class deals with the bending and longitudinal wave properties 
    of beams. In future implementations all three motions of beams should be 
    implemented.
    """
    
    def __init__(self,cross_section=geoPC.CrossSection,iso_mat=mc.IsoMat):
        """
        Constructor of BeamProp.

        Parameters
        ----------
        cross_section : CrossSection, optional
            Cross section of beam. The default is geoPC.CrossSection.
        iso_mat : IsoMat, optional
            Isotropic material of beam. The default is mc.IsoMat.

        Returns
        -------
        None.

        """
        self.cross_section = cross_section
        self.iso_mat = iso_mat

    def __str__(self):
        """
        String output of BeamProp.

        Returns
        -------
        str

        """
        _str  = "BeamProp: \n"
        _str += "cross_ection:\n{0}".format(self.cross_section)
        _str += "iso_mat:\n{0}".format(self.iso_mat)
        
        return _str
            
    def __repr__(self):
        """
        Reps of acoutic tube

        Returns
        -------
        None.

        """
        
        return "BeamProp(cross_section={0},iso_mat={1})".format(self.cross_section.__repr__(),self.iso_mat.__repr__())


    @property
    def Bx(self):
        """
        Bending stiffness around x-axis.

        Returns
        -------
        float
            bending stiffness.

        """
        return self.iso_mat.E*self.cross_section.Ix
    
    @property
    def By(self):
        """
        Bending stiffness in around x

        Returns
        -------
        float
            bending stiffness.

        """
        return self.iso_mat.E*self.cross_section.Iy
    
    @property
    def mass_per_length(self):
        """
        mass per beam length

        Returns
        -------
        float
            mass per length.

        """
        return self.cross_section.area*self.iso_mat.rho0
    
    def kx(self,omega):
        """
        wavenumber for bending around x

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        float
            bending wavenumber.

        """
        BperM = self.Bx/self.mass_per_length
        return (omega*omega/BperM)**0.25

    def ky(self,omega):
        """
        wavenumber for bending around y

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        float
            bending wavenumber.

        """
        BperM = self.By/self.mass_per_length
        return (omega*omega/BperM)**0.25

    def wavelength_x(self,omega):
        """
        wavelength for bending around x

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        float
            bending wavelength.

        """
        return np.pi*2/self.kx(omega)

    def wavelength_y(self,omega):
        """
        wavelength for bending around y

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        float
            bending wavelength.

        """
        return np.pi*2/self.ky(omega)
            
    def c_phase(self,omega):
        """
        phase wave speed for bending around x

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        float
            phase wave speed.

        """
        BperM = self.Bx/self.mass_per_length
        return np.sqrt(omega)*(BperM)**0.25
    
    def c_group(self,omega):
        """
        group wave speed for bending around x

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        float
            phase wave speed.

        """
        return 2*self.c_phase(omega)
    
    def z_beam_inf_x(self,omega):
        """
        Point impedance of half beam for bending in x

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        complex
            half beam impedance.

        """
        
        return (2+2j)*self.Bx*self.kx(omega)**3/omega
    
    def z_beam_inf_y(self,omega):
        """
        Point impedance of half beam for bending in y

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        complex
            half beam impedance.

        """
        return (2+2j)*self.By*self.ky(omega)**3/omega

    def stiffness_bending_x(self,omega):
        """
        Point stiffness of half beam for bending in x

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        complex
            half beam stiffness.

        """        
        
        return (self.Bx*self.kx(omega)**2)*(0.5j-0.5)
    
    def stiffness_bending_y(self,omega):
        """
        Point stiffness of half beam for bending in y

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        complex
            half beam stiffness.

        """        
        
        return (self.By*self.ky(omega)**2)*(0.5j-0.5)

    def stiffness_longitudinal_z(self,omega):
        """
        Point stiffness of half beam for longitudinal in z.

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        complex
            half beam stiffness.

        """
        return 1j*omega*self.cross_section.area*np.sqrt(self.iso_mat.E*self.iso_mat.rho0)        

        
class PlateProp:
    """
    Class for thin plates.
    
    The plate property class deals with the dynamics of thin plates.
    This comprises all waves and dynamics of isotropic plates.
    
    Due to this simplicity there are only two attributes
    
    Attributes
    ----------
    thickness : float
        plate thickness.
    material : IsoMat
        isotropic plate material.
    
    
    """
    
    def __init__(self,thickness,material):
        """
        Class constructor of plate property

        Parameters
        ----------
        thickness : float
            plate thickness.
        material : IsoMat
            isotropic plate material.

        Returns
        -------
        None.

        """
        
        self.thickness = thickness
        self.material  = material
    
    def __str__(self):
        _str  = "PlateProp object with properties:\n"
        _str += "thickness      : {0}\n".format(self.thickness)
        _str += "material\n--------\n{0}".format(self.material)
        return _str

    def __repr__(self):
        str_ = "PlateProp({0},{1})".format(self.thickness,self.material)
        return str_
    
    @property
    def B(self):
        """
        real bending stiffness
        
        Returns
        -------
        float
            bending stiffness.

        """
        
        nu = self.material.nu
        return self.material.E*self.thickness**3/(12*(1-nu*nu))
    
    @property
    def B_complex(self):
        """
        complex bending stiffness with consideration of damping loss.
        
        Returns
        -------
        complex
            bending stiffness.

        """
        
        eta = self.material.eta
        return self.B*(1+1j*eta) 

    @property      
    def S(self):
        """
        transversal stiffness
        
        Returns
        -------
        float
            transversal stiffness.

        """
        return self.material.E*self.thickness/2/(1+self.material.nu)
    
    @property      
    def S_complex(self):
        """
        complex transversal stiffness.
        
        Returns
        -------
        complex
            complex transversal stiffness.

        """
        eta = self.material.eta
        return self.material.E*(1+1j*eta)*self.thickness/2/(1+self.material.nu)


    @property      
    def C(self):
        """
        longitudinal stiffness

        Returns
        -------
        float
            transversal stiffness.
        """

        return self.material.E*self.thickness/(1-self.material.nu**2)

        
    @property
    def B_per_M(self):
        """
        Bending stiffness per area mass

        Returns
        -------
        float
            Bending stiffness per area specific mass.

        """
        
        return self.B/self.mass_per_area

    @property
    def B_per_M_complex(self):
        """
        Complex bending stiffness per area mass

        Returns
        -------
        complex
            Bending stiffness per area specific mass.

        """ 
        
        return self.B_complex/self.mass_per_area
        
    @property
    def mass_per_area(self):
        """
        Mass per area

        Returns
        -------
        float
            mass per area.

        """
        
        return self.thickness*self.material.rho0    
        
    def c_B_phase(self,omega):
        """
        Bending phase wave speed

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        float
            phase wave speed.

        """
        
        BperA = self.B/self.mass_per_area
        return np.sqrt(omega)*(BperA)**0.25

    def c_L(self):
        """
        Longitudinal wave speed

        Returns
        -------
        float
            wave speed.

        """        
        
        return np.sqrt(self.material.E_complex/(1-self.material.nu**2)/self.material.rho0)

    def c_T(self):
        """
        Shear wave speed

        Returns
        -------
        float
            wave speed.

        """        
        return self.material.c_S

    def c_S(self):
        """
        Shear wave speed

        Returns
        -------
        float
            wave speed.

        """        
        return self.material.c_S

    
    def c_B_group(self,omega):
        """
        Bending group wave speed

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        float
            group wave speed.

        """        
        
        return 2*self.c_B_phase(omega)

    def wavenumber_B(self,omega,comp = False):
        """
        Bending wavenumber

        Parameters
        ----------
        omega : float
            angular frequency.
        comp : bool, optional
            Switch for complex consideration of damping. The default is False.

        Returns
        -------
        complex
            wavenumber.

        """

        if comp:
            BperM = self.B_per_M_complex
        else:
            BperM = self.B_per_M
                
        return np.sqrt(omega)/((BperM)**0.25)

    def wavenumber_constants(self,omega):
        """
        Bending wavenumber

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        tuple
            (kL,LS,kB).

        """

        BperM = self.B_per_M
                
        return (omega/self.c_L(),\
                omega/self.c_S(),\
                np.sqrt(omega)/((BperM)**0.25) )

    def wavenumber_B_4(self,omega):
        BperM = self.B_per_M
        return omega**2/BperM

    def wavenumber_T(self,omega=0.):
        return omega/self.material.c_S

    def wavenumber_S(self,omega=0.):
        """
        Shear wavenumber

        Parameters
        ----------
        omega : float, optional
            angular frequency. The default is 0.

        Returns
        -------
        float
            wavenumber.

        """
        return omega/self.material.c_S


    def wavenumber_L(self,omega=0.):
        """
        Longitudinal wavenumber

        Parameters
        ----------
        omega : float, optional
            angular frequency. The default is 0.

        Returns
        -------
        float
            wavenumber.

        """
        return omega/self.c_L()
      
    def wavelength_B(self,omega):
        """
        Bending wavelength

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        float
            wavelength.

        """
        return 2*np.pi*self.c_B_phase(omega)/omega

    def wavelength_T(self,omega=0.):
        return 2*np.pi*self.c_S()/omega

    def wavelength_S(self,omega=0.):
        """
        Shear wavelength

        Parameters
        ----------
        omega : float, optional
            angular frequency. The default is 0.

        Returns
        -------
        float
            wavelength.

        """
        return 2*np.pi*self.c_S()/omega

    def wavelength_L(self,omega):
        """
        Longitudial wavelength

        Parameters
        ----------
        omega : float, optional
            angular frequency. The default is 0.

        Returns
        -------
        float
            wavelength.

        """
        return 2*np.pi*self.c_L()/omega

    def muL(self,omega,kx):
        """
        Longitudiunal wave propagation constant in y-direction 

        Parameters
        ----------
        omega : TYPE
            DESCRIPTION.
        kx : float
            edge wavenumber.

        Returns
        -------
        uL : complex
            propagation constant.

        """
        kL = self.wavenumber_L(omega)
        uL = -np.sqrt(kx*kx-np.complex128(kL)**2)
        return uL


    def muS(self,omega,kx):
        """
        Shear wave propagation constant in y-direction 

        Parameters
        ----------
        omega : TYPE
            DESCRIPTION.
        kx : float
            edge wavenumber.

        Returns
        -------
        uS : complex
            propagation constant.

        """
        kS = self.wavenumber_T(omega)
        uS = -np.sqrt(kx*kx-np.complex128(kS)**2)
        return uS

    def muT(self,omega,kx):
        kT = self.wavenumber_T(omega)
        uT = -np.sqrt(kx*kx-np.complex128(kT)**2)
        return uT

             
    def w_inf(self,omega,r,Fz):
        """
        Displacement due to normal point force excitatoin 

        Parameters
        ----------
        omega : float
            angular frequency.
        r : float
            distance source - receiver.
        Fz : complex
            Normal force amplitude.

        Returns
        -------
        complex
            Displacement.

        """
        
        k = self.wavenumber_B(omega)
        B = self.B
        return np.where(r==0,\
                     Fz/(8j*B*k*k),
                     Fz/(8j*B*k*k)*(scs.hankel2(0,k*np.abs(r))- \
                        scs.hankel2(0,-1j*k*np.abs(r))))
    @property    
    def point_impedance(self):
        """
        Normal force point impedance

        Returns
        -------
        float
            mechanical impedance.

        """
        
        return 8*np.sqrt(self.B*self.mass_per_area)

    def force_excitation_power(self,omega,force=1.):
        """
        Power input due to normal force input

        Parameters
        ----------
        omega : float
            angular frequency.
        force : complex, optional
            normal force. The default is 1.

        Returns
        -------
        float
            power.

        """
        return force**2/self.point_impedance
        
    def point_stiffness(self,omega):
        """
        Normal force point stiffness

        Returns
        -------
        float
            mechanical stiffness.

        """

        return 1j*omega*self.point_impedance
            
    @property    
    def point_impedance_edge(self):
        """
        Normal force point edge impedance

        Returns
        -------
        float
            mechanical impedance.

        """
        return 3.5*np.sqrt(self.B*self.mass_per_area)
        
    def point_stiffness_edge(self,omega):
        """
        Normal force point edge stiffness

        Returns
        -------
        float
            mechanical stiffness.

        """
        return 1j*omega*self.point_impedance_edge
    
    def coincidence_frequency(self,c0 = 343.):
        """
        coincidence frequency of flat plate 

        Parameters
        ----------
        c0 : float
            speed of sound of fluid

        Returns
        -------
        float
            angular coincidence frequency in s^(-1)

        """
                
        return np.sqrt(1/np.real(self.B_per_M))*c0**2
    
    def transfer_impedance(self,omega,kx):
        """
        Transfer impedance of infinite plates
        
        The transfer impedanc is used in the transfer matrix applications
        using infinite plate theory
        
        Parameters
        ----------
        omega : float
            angular frequency
        kx : float
            wavenumer of incomiing wave
               
        Returns
        -------
        complex
            The transfer impedance of infinite plates
        
        """
        kB4 = omega**2/self.B_per_M_complex
        return 1j*omega*self.mass_per_area*(1-kx**4/kB4)
    
    def transmission_coefficient_angular(self,omega,fluid1,theta=0,fluid2='none'):
        """
        Transmission coefficient of plane wave transmission through infinite plates
                        
        Parameters
        ----------
        omega : float
            angular frequency
        theta : float
            angle of incidence in radiants 0 <= theta <= pi/2
        fluid1 : fluid
            fluid on irradiation side or both if fluid 2 is not given
        fluid2 : fluid
            fluid on transmission side 
            
        Returns
        -------
        float
            angular transmission coefficient of plateprop
        """
        
        # Determine number of half spaces
        
        ka  = np.real(fluid1.wavenumber(omega))
        m   = self.mass_per_area
        kB4 = omega**2/self.B_per_M_complex
        z0  = fluid1.z0
        
        return 1/(1 + (m*omega/2/z0*np.cos(theta)*np.abs((ka*np.sin(theta))**4/kB4-1))**2)
    
    def transmission_coefficient_diffuse(self,omega,fluid1,fluid2='none',theta_max=np.pi*78/180,theta_step = np.pi/180):
        """
        Diffuse transmission coeffient of infinite plates
        
        The diffuse field is implemented by angle intergration
                        
        Parameters
        ----------
        omega : float
            angular frequency
        fluid1 : fluid
            fluid on irradiation side or both if fluid 2 is not given
        fluid2 : fluid
            fluid on transmission side 
        theta_max : float, optional
            max angle of incidence in radiants 0 <= theta <= pi/2. Default value is np.pi*78/180
        theta_step : float, optional
            max angle of incidence in radiants 0 <= theta <= pi/2. Default value is np.pi/180
            
            
        Returns
        -------
        float
            Diffuse field transmission coefficient
        """
        
        # Determine number of half spaces
        
        tau = np.zeros(omega.shape)
        
        # Get average denominator
        denom = 0.5*np.sin(theta_max)**2
        
        theta_ = np.linspace(0,theta_max,int(np.floor((theta_max+theta_step)/theta_step)))
        
        for ifreq in range(len(omega)):
            #f = lambda x: self.transmission_coefficient_angular(omega[ifreq],fluid1,x)*np.cos(x)*np.sin(x)/denom
            tau_kx = self.transmission_coefficient_angular(omega[ifreq],fluid1,theta_)
            tau_kx[np.isnan(tau_kx)] = 0.
            tau_kx[tau_kx<0] = 0.

            #tau[ifreq],prec = integrate.quad(f, 0, theta_max)
            tau[ifreq] = integrate.trapz(tau_kx*np.sin(theta_)*np.cos(theta_), theta_)/denom

        return tau  
        
        
      
    def edge_radiation_stiffness_wavenumber_LM(self,omega,wavenumber,wave_DOF = 0):
        """
        Structural radiation stiffness of straight plate edges

        This method calculates the structural dynamic radiation
        stiffness matrix of a semi-infinite plate in the wave number domain. The
        force-displacement relations in wave number domain are derived from the
        harmonic solution of the vibrating plate equation, as explained in [Pei2022]_.
        Given a cartesian system with
        the x-axis along the plate edge and the y-axix looking inside the plate, 
        four degrees of freedom are considered at the edge of the plate: three
        displacement (u,v,w) and the rotation (theta) in the edge direction (x),
        so that the force-displacement function is a 4 X 4 matrix, function of
        wavenumber and time frequency, where the out-of-plane behavior is 
        decoupled from the in-plane behavior. The matrix is given in terms of its
        single coefficients a_ij.
              
        Parameters
        ----------
        omega : float
            angular frequency
        wavenumber : float
            wavenumber kx 
        wave_DOF : int
            identifier of radiation matrix 0 all waves, 1,2,3 correponds to longitudinal, shear and bending wave 

        Returns
        -------
        LinearMatrix
            matrix [4 x 4] of the stiffness element
            
        """
         
        Kx=wavenumber
        # Npoint=len(Kx);
 
 
        #eta = self.material.eta
        kB = self.wavenumber_B(omega) #bending wavenumber
        kL = self.wavenumber_L(omega) #longitudinal wavenumber
        kS = self.wavenumber_T(omega) #shear wavenumber
        B = self.B # bending stiffness with damping
        S = self.S #shear stiffness with damping
        nu= self.material.nu #poisson
          
 
        # in-plane motion, roots
        
        if wave_DOF ==  5:
            uL= self.muL(omega,wavenumber)
            uS= self.muT(omega,wavenumber)
        else:
            uL= -np.sqrt(Kx**2-np.complex128(kL)**2)
            uS= -np.sqrt(Kx**2-np.complex128(kS)**2)

        #kyL = np.sqrt(np.complex128(kL**2-Kx**2))
        kyS = np.sqrt(np.complex128(kS**2-Kx**2))

        data_ = np.zeros((4,4,Kx.size),dtype=np.complex128)

        Sfac  = S/(Kx**2-uL*uS)
        #Denominator for in-plane waves with kx < kL
        #Sdenom = (Kx**2+kyL*kyS)
        #Sfac   = S/Sdenom
        #Denominator for in-plane waves with kx > kL
        
        S2fac  = S/(Kx**2-1j*kyS*uL) 
    
        # out-of-plane motion, roots
        uB1 = -np.sqrt(Kx**2+kB**2)
        uB2 = -np.sqrt(Kx**2-np.complex128(kB)**2)

        #uB1[np.isreal(uB1)] = 1.E100

        if wave_DOF == 0:
            #uL[np.isreal(uL)] = 0.
            #uS[np.isreal(uS)] = 0.
            
            # in-plane motion, matrix function elements langs_SL = -1 gives Langleys convention
            data_[0,0,:] = -Sfac*(uL*kS**2) 
            data_[0,1,:] = -Sfac*(1j*Kx*(2*uS*uL+kS**2-2*Kx**2)) 
            data_[1,0,:] = -data_[0,1,:]#*Sfac*(1j*Kx*(2*uS*uL+kS**2-2*Kx**2))#asym_SL*data_[0,1,:]
            data_[1,1,:] = -Sfac*(uS*kS**2)

     
            # out-of-plane motion, matrix function elements
            data_[2,2,:] = -B*(uB1**2*uB2+uB2**2*uB1) 
            data_[2,3,:] = B*(uB2*uB1+nu*Kx**2)
            data_[3,2,:] = data_[2,3,:]
            data_[3,3,:] = -B*(uB2+uB1)
            
            
        elif wave_DOF == 1: #in [1]: NOt VALID
            data_[0,0,:] = -Sfac*2*(uL*Kx**2) 
            data_[0,1,:] = -Sfac*2j*Kx*uS*uL 
            data_[1,0,:] = -Sfac*1j*Kx*(2*Kx**2-kS**2)
            data_[1,1,:] = Sfac*uS*(2*Kx**2-kS**2)

        elif wave_DOF == 2: # not valid
            data_[0,0,:] = Sfac*uL*(2*Kx**2-kS**2)
            data_[0,1,:] = Sfac*1j*Kx*(2*Kx**2-kS**2) 
            data_[1,0,:] = Sfac*2j*Kx*uS*uL
            data_[1,1,:] = -Sfac*(uS*Kx**2)
            
        elif wave_DOF in (3,4):
            # out-of-plane motion, matrix function elements
            data_[2,2,:] = -B*(uB1**2*uB2+uB2**2*uB1) 
            data_[2,3,:] = B*(uB2*uB1+nu*Kx**2)
            data_[3,2,:] = data_[2,3,:]
            data_[3,3,:] = -B*(uB2+uB1)


            
        elif wave_DOF == 5: #in plane Eq. (8.121)
            # in-plane motion, matrix function elements
            data_[0,0,:] = -Sfac*(uL*kS**2) 
            data_[0,1,:] = -Sfac*(1j*Kx*(2*uS*uL+kS**2-2*Kx**2)) 
            data_[1,0,:] = -data_[0,1,:]#*Sfac*(1j*Kx*(2*uS*uL+kS**2-2*Kx**2))#asym_SL*data_[0,1,:]
            data_[1,1,:] = -Sfac*(uS*kS**2)
            
            # seperate version for kS > kX > kL 
            # Did not work
            #ix = np.logical_and(Kx < kS,Kx >= kL)
            #data_[0,0,ix] =  -S2fac[ix]*(uL[ix]*kS**2)  
            #data_[0,1,ix] =  -S2fac[ix]*Kx[ix]*(2*kyS[ix]*uL[ix]+1j*(kS**2-2*Kx[ix]**2)) 
            #data_[1,0,ix] =  -data_[0,1,ix]
            #data_[1,1,ix] =  -S2fac[ix]*1j*kyS[ix]*kS**2       
            
            
        return mC.LinearMatrix(data_)

    def edge_radiation_stiffness_wavenumber(self,omega,wavenumber,wave_DOF = 0):
        """
        Structural radiation stiffness of straight plate edges

        This method calculates the structural dynamic radiation
        stiffness matrix of a semi-infinite plate in the wave number domain. The
        force-displacement relations in wave number domain are derived from the
        harmonic solution of the vibrating plate equation, as explained in [Pei2022]_.
        Given a cartesian system with
        the x-axis along the plate edge and the y-axix looking inside the plate, 
        four degrees of freedom are considered at the edge of the plate: three
        displacement (u,v,w) and the rotation (theta) in the edge direction (x),
        so that the force-displacement function is a 4 X 4 matrix, function of
        wavenumber and time frequency, where the out-of-plane behavior is 
        decoupled from the in-plane behavior. The matrix is given in terms of its
        single coefficients a_ij.
              
        Parameters
        ----------
        omega : float
            angular frequency
        wavenumber : float
            wavenumber kx 
        wave_DOF : int
            identifier of radiation matrix 0 all waves, 1,2,3 correponds to longitudinal, shear and bending wave 

        Returns
        -------
        nd.array
            matrix [Nkx x 4 x 4] of the stiffness element
            
        """
         
        Kx=wavenumber
        # Npoint=len(Kx);
 
 
        #eta = self.material.eta
        # kB = self.wavenumber_B(omega) #bending wavenumber
        # kL = self.wavenumber_L(omega) #longitudinal wavenumber
        # kS = self.wavenumber_T(omega) #shear wavenumber
        B = self.B # bending stiffness with damping
        S = self.S #shear stiffness with damping
        nu= self.material.nu #poisson
          
        #use more efficient method that provide all results in one step
        (kL,kS,kB) = self.wavenumber_constants(omega)
        
    
 
        # in-plane motion, roots
        uL= -np.sqrt(Kx**2-np.complex128(kL)**2)
        uS= -np.sqrt(Kx**2-np.complex128(kS)**2)

        data_ = np.zeros((Kx.size,4,4),dtype=np.complex128)

        Sfac  = S/(Kx**2-uL*uS)
    
        # out-of-plane motion, roots
        uB1 = -np.sqrt(Kx**2+kB**2)
        uB2 = -np.sqrt(Kx**2-np.complex128(kB)**2)

        #uB1[np.isreal(uB1)] = 1.E100

        if wave_DOF in (3,4):
            # out-of-plane motion, matrix function elements
            data_[:,2,2] = -B*(uB1**2*uB2+uB2**2*uB1) 
            data_[:,2,3] = B*(uB2*uB1+nu*Kx**2)
            data_[:,3,2] = data_[:,2,3]
            data_[:,3,3] = -B*(uB2+uB1)
            
        elif wave_DOF == 5: #in plane
            # seperate version for kX > kL
            # in-plane motion, matrix function elements
            data_[:,0,0] = -Sfac*(uL*kS**2) 
            data_[:,0,1] = -Sfac*(1j*Kx*(2*uS*uL+kS**2-2*Kx**2)) 
            data_[:,1,0] = -data_[:,0,1]#*Sfac*(1j*Kx*(2*uS*uL+kS**2-2*Kx**2))#asym_SL*data_[0,1,:]
            data_[:,1,1] = -Sfac*(uS*kS**2)
        elif wave_DOF == 0:
            
            # in-plane motion, matrix function elements langs_SL = -1 gives Langleys convention
            data_[:,0,0] = -Sfac*(uL*kS**2) 
            data_[:,0,1] = -Sfac*(1j*Kx*(2*uS*uL+kS**2-2*Kx**2)) 
            data_[:,1,0] = -data_[:,0,1]#*Sfac*(1j*Kx*(2*uS*uL+kS**2-2*Kx**2))#asym_SL*data_[0,1,:]
            data_[:,1,1] = -Sfac*(uS*kS**2)

     
            # out-of-plane motion, matrix function elements
            data_[:,2,2] = -B*(uB1**2*uB2+uB2**2*uB1) 
            data_[:,2,3] = B*(uB2*uB1+nu*Kx**2)
            data_[:,3,2] = data_[:,2,3]
            data_[:,3,3] = -B*(uB2+uB1)
            
            
        elif wave_DOF == 1: #in [1]: NOt VALID
            data_[:,0,0] = -Sfac*2*(uL*Kx**2) 
            data_[:,0,1] = -Sfac*2j*Kx*uS*uL 
            data_[:,1,0] = -Sfac*1j*Kx*(2*Kx**2-kS**2)
            data_[:,1,1] = Sfac*uS*(2*Kx**2-kS**2)

        elif wave_DOF == 2: # not valid
            data_[:,0,0] = Sfac*uL*(2*Kx**2-kS**2)
            data_[:,0,1] = Sfac*1j*Kx*(2*Kx**2-kS**2) 
            data_[:,1,0] = Sfac*2j*Kx*uS*uL
            data_[:,1,1] = -Sfac*(uS*Kx**2)
            

            
        return data_

    
    def edge_imaginary_radiation_stiffness_wavenumber(self,omega,wavenumber,wave_DOF=0):
        """
        Structural edge radiation stiffness for specific wave types
        
        Due to special relationships there are different formulations of the
        imaginary radiation stiffness as following directly from the 
        radation_stiffness_wavenumber except for bendig waves where both
        expressions are similar. See [Pei2022]_ for details
        
        Will be deprecated in the new future, because it is equivalent to the skew hermitian matrix.
        
        Parameters
        ----------
        omega : float
            angular frequency
        wavenumber : float
            wavenumber kx 
        wave_DOF : int
            identifier of radiation matrix 0 all waves, 1,2,(3 or 4) correponds 
            to longitudinal, shear and bending wave
            5 stands for inplane (1+2)

        Returns
        -------
        LinearMatrix
            matrix [4 x 4] of the stiffness element
            
        """

        if wave_DOF == 0: # Return standard values
            return self.edge_radiation_stiffness_wavenumber(omega,wavenumber,wave_DOF=wave_DOF).imag()

         
        Kx=wavenumber
        # Npoint=len(Kx);
 
 
        #eta = self.material.eta
        kB = self.wavenumber_B(omega) #bending wavenumber
        kL = self.wavenumber_L(omega) #longitudinal wavenumber
        kS = self.wavenumber_T(omega) #shear wavenumber
        B = self.B # bending stiffness with damping
        S = self.S #shear stiffness with damping
        nu= self.material.nu #poisson
          
 
        # in-plane motion, roots
        
        uL= self.muL(omega,wavenumber)
        uS= self.muT(omega,wavenumber)
        
        kyL = np.sqrt(np.complex128(kL**2-Kx**2))
        kyS = np.sqrt(np.complex128(kS**2-Kx**2))
                 
        # Try removing not allowed values
        #uL[np.isreal(uL)] = 0
        #uS[np.isreal(uS)] = 0

        data_ = np.zeros((4,4,Kx.size),dtype=np.complex128)

        #Denominator for in-plane waves with kx < kL
        Sdenom = (Kx**2+kyL*kyS)
        #Denominator for in-plane waves with kx > kL
        S2denom = (Kx**4+uL**2*kyS**2)

        
        Sfac   = S/Sdenom
        S2fac  = S/S2denom*kS**2*kyS    

        # out-of-plane motion, roots
        uB1 = -np.sqrt(Kx**2+kB**2)
        uB2 = -np.sqrt(Kx**2-np.complex128(kB)**2)

        #uB1[np.isreal(uB1)] = 1.E100
        
            
        if wave_DOF == 1: #in [1] Only for demonstration purpose, see Appendix B
            fak1 = Sfac/Sdenom*kS**2*kyL
            ix = Kx < kL
            data_[0,0,ix] =  fak1[ix]*Kx[ix]**2 
            data_[0,1,ix] = -fak1[ix]*kyS[ix]*Kx[ix] 
            data_[1,0,ix] =  data_[0,1,ix]
            data_[1,1,ix] =  fak1[ix]*kyS[ix]**2

        elif wave_DOF == 2: # Only for demonstration purpose, see Appendix B
            fak1 = Sfac/Sdenom*kS**2*kyS # Sdenom is squared for single wave
            ix = Kx <= kS
            data_[0,0,ix] =  fak1[ix]*kyL[ix]**2 
            data_[0,1,ix] =  fak1[ix]*kyL[ix]*Kx[ix] 
            data_[1,0,ix] =  data_[0,1,ix]
            data_[1,1,ix] =  fak1[ix]*Kx[ix]**2
            ix = np.logical_and(Kx < kS,Kx >= kL)
            fak1 = S/S2denom*kS**2*kyS
            data_[0,0,ix] =  fak1[ix]*uL[ix]**2 
            data_[0,1,ix] =  1j*fak1[ix]*uL[ix]*Kx[ix] 
            data_[1,0,ix] =  -data_[0,1,ix]
            data_[1,1,ix] =  fak1[ix]*Kx[ix]**2
            
        elif wave_DOF in (3,4): # Eq before (B.53)
            # out-of-plane motion, matrix function elements
            data_[2,2,:] = np.imag(-B*(uB1**2*uB2+uB2**2*uB1)) 
            data_[2,3,:] = np.imag(B*(uB2*uB1+nu*Kx**2))
            data_[3,2,:] = data_[2,3,:]
            data_[3,3,:] = np.imag(-B*(uB2+uB1))
        elif wave_DOF == 5: #in plane
            # for Kx < kL
            ix = Kx < kL
            data_[0,0,ix] = np.imag(-Sfac[ix]*(uL[ix]*kS**2)) 
            data_[0,1,ix] = np.imag(-Sfac[ix]*(1j*Kx[ix]*(2*uS[ix]*uL[ix]+kS**2-2*Kx[ix]**2))) 
            data_[1,0,ix] = -data_[0,1,ix]#*Sfac*(1j*Kx*(2*uS*uL+kS**2-2*Kx**2))#asym_SL*data_[0,1,:]
            data_[1,1,ix] = np.imag(-Sfac[ix]*(uS[ix]*kS**2))                
            # for kL < Kx <= kS Eq (B.76) must be used
            
            ix = np.logical_and(Kx < kS,Kx >= kL)
            data_[0,0,ix] =  S2fac[ix]*uL[ix]**2 
            data_[0,1,ix] =  1j*S2fac[ix]*uL[ix]*Kx[ix] 
            data_[1,0,ix] =  -data_[0,1,ix]
            data_[1,1,ix] =  S2fac[ix]*Kx[ix]**2            
            
            
                                                                          
        return mC.LinearMatrix(data_)
    
    def edge_skew_radiation_stiffness_wavenumber_LM(self,omega,wavenumber,wave_DOF=0):
        """
        Structural edge radiation stiffness for specific wave types
        
        Due to special relationships there are different formulations of the
        imaginary radiation stiffness as following directly from the 
        radation_stiffness_wavenumber except for bendig waves where both
        expressions are similar. See [Pei2022]_ for details
                
        Parameters
        ----------
        omega : float
            angular frequency
        wavenumber : float
            wavenumber kx 
        wave_DOF : int
            identifier of radiation matrix 0 all waves, 1,2,(3 or 4) correponds 
            to longitudinal, shear and bending wave
            5 stands for inplane (1+2)

        Returns
        -------
        LinearMatrix
            matrix [4 x 4] of the stiffness element
            
        """

        M = self.edge_radiation_stiffness_wavenumber_LM(omega,wavenumber,wave_DOF=wave_DOF)

        return -0.5j*(M-M.H()) 

    def edge_skew_radiation_stiffness_wavenumber(self,omega,wavenumber,wave_DOF=0):
        """
        Structural edge radiation stiffness for specific wave types.
        
        Due to special relationships there are different formulations of the
        imaginary radiation stiffness as following directly from the 
        radation_stiffness_wavenumber except for bendig waves where both
        expressions are similar. See [Pei2022]_ for details
                
        Parameters
        ----------
        omega : float
            angular frequency
        wavenumber : float
            wavenumber kx 
        wave_DOF : int
            identifier of radiation matrix 0 all waves, 1,2,(3 or 4) correponds 
            to longitudinal, shear and bending wave
            5 stands for inplane (1+2)

        Returns
        -------
        np.ndarray
            matrix [Nx x 4 x 4] of the stiffness element
            
        """
        M = self.edge_radiation_stiffness_wavenumber(omega,wavenumber,wave_DOF=wave_DOF)

        return -0.5j*(M-mC.hermitian(M)) 

    
    def wave_transformation_matrix_LM(self,omega,wavenumber,inv=False):
        """
        Transformation matrix from wave amplitude coordinates into edge harmonic displacement.
    
        LinearMatrix version
    
        Parameters
        ----------
        omega : float
            angular frequency
        wavenumber : float
            wavenumber
        inv : bool
            switch for inverse version
            
        Returns
        -------
        LinearMmatrix
            [4 x 4] transformation matrix
        """         
        Kx=np.real(wavenumber);
        # Npoint=len(Kx);
 
 
        # use methods from property2structure to get properties
        #eta = self.material.eta
        kB = self.wavenumber_B(omega) #bending wavenumber
        #kL = self.wavenumber_L(omega) #longitudinal wavenumber
        #kS = self.wavenumber_T(omega) #shear wavenumber
        #nu = self.material.nu #poisson
          
 
        # in-plane motion, roots
        #uL= -np.sqrt(Kx**2-np.complex128(kL)**2)
        #uS= -np.sqrt(Kx**2-np.complex128(kS)**2)
        uL= self.muL(omega,wavenumber)
        uS= self.muT(omega,wavenumber)        

        # out-of-plane motion, roots
        uB1 = -np.sqrt(Kx**2+kB**2)
        uB2 = -np.sqrt(Kx**2-np.complex128(kB)**2)

        # Try removing not allowed values
        #uL[np.isreal(uL)] = 0.
        #uS[np.isreal(uS)] = 0.
        
        data_ = np.zeros((4,4,len(wavenumber)),dtype=np.complex128)

        if inv:
            facSL = 1/(Kx**2-uS*uL)
            # in-plane motion, matrix function elements
            data_[0,0,:] = facSL*Kx                        # Fx - u
            data_[0,1,:] = 1j*facSL*uS                     # Fx - v 
            data_[1,0,:] = 1j*facSL*uL                     # Fy - u
            data_[1,1,:] = -facSL*Kx                       # Fy - v
            # out-of-plane motion, matrix function elements
            facB = -1/(uB1-uB2)
            data_[2,2,:] = facB*uB2                     # Fx - w
            data_[2,3,:] = -facB                     # Fx - beta
            data_[3,2,:] = -facB*uB1                   # Mx - w
            data_[3,3,:] = facB                  # Mx - beta
            
        else:
            
            # in-plane motion, matrix function elements
            data_[0,0,:] = Kx                        # Fx - u
            data_[0,1,:] = 1j*uS                     # Fx - v 
            data_[1,0,:] = 1j*uL                     # Fy - u
            data_[1,1,:] = -Kx                       # Fy - v
            # out-of-plane motion, matrix function elements
            data_[2,2,:] = 1                     # Fx - w
            data_[2,3,:] = 1                     # Fx - beta
            data_[3,2,:] = uB1                   # Mx - w
            data_[3,3,:] = uB2                   # Mx - beta
                
        return mC.LinearMatrix(data_)
    
    def wave_transformation_matrix(self,omega,wavenumber,inv=False):
        """
        Transformation matrix from wave amplitude coordinates into edge harmonic displacement.
    
        Parameters
        ----------
        omega : float
            angular frequency
        wavenumber : float
            wavenumber
        inv : bool
            switch for inverse version

        Returns
        -------
        np.ndarray
            [Nkx x 4 x 4] transformation matrix
        """         
        Kx=np.real(wavenumber);
        # Npoint=len(Kx);
 
 
        # use methods from property2structure to get properties
        kB = self.wavenumber_B(omega) #bending wavenumber
          
 
        # in-plane motion, roots
        uL= self.muL(omega,wavenumber)
        uS= self.muT(omega,wavenumber)        

        # out-of-plane motion, roots
        uB1 = -np.sqrt(Kx**2+kB**2)
        uB2 = -np.sqrt(Kx**2-np.complex128(kB)**2)
        
        data_ = np.zeros((len(wavenumber),4,4),dtype=np.complex128)

        if inv:
            facSL = 1/(Kx**2-uS*uL)
            # in-plane motion, matrix function elements
            data_[:,0,0] = facSL*Kx                        # Fx - u
            data_[:,0,1] = 1j*facSL*uS                     # Fx - v 
            data_[:,1,0] = 1j*facSL*uL                     # Fy - u
            data_[:,1,1] = -facSL*Kx                       # Fy - v
            # out-of-plane motion, matrix function elements
            facB = -1/(uB1-uB2)
            data_[:,2,2] = facB*uB2                     # Fx - w
            data_[:,2,3] = -facB                     # Fx - beta
            data_[:,3,2] = -facB*uB1                   # Mx - w
            data_[:,3,3] = facB                  # Mx - beta
            
        else:
            
            # in-plane motion, matrix function elements
            data_[:,0,0] = Kx                        # Fx - u
            data_[:,0,1] = 1j*uS                     # Fx - v 
            data_[:,1,0] = 1j*uL                     # Fy - u
            data_[:,1,1] = -Kx                       # Fy - v
            # out-of-plane motion, matrix function elements
            data_[:,2,2] = 1                     # Fx - w
            data_[:,2,3] = 1                     # Fx - beta
            data_[:,3,2] = uB1                   # Mx - w
            data_[:,3,3] = uB2                   # Mx - beta
                
        return data_
    
    def plate_wavenumber(self,omega,wave_DOF):
        """
        Wavenumber of plates

        Parameters
        ----------
        omega : float
            angular frequency.
        wave_DOF : int
            wave degree of freedom.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        k_plate : TYPE
            DESCRIPTION.

        """
        
        
        if wave_DOF in (1,):
            k_plate = np.real(self.wavenumber_L(omega))
        elif wave_DOF in (2,5): # 5 for common LS treatment takes k_T because k_T > k_L
            k_plate = np.real(self.wavenumber_S(omega))
        elif wave_DOF in (3,4):
            k_plate = np.real(self.wavenumber_B(omega))
        else:
            raise ValueError('wave_DOF argument must be in [1,2,3,4,5] but is {0}'.format(wave_DOF))
            
        return k_plate
    
    def edge_wave_amplitude_radiation_stiffness(self,omega,wavenumber,wave_DOF):
        """
        Efficient radiation stiffness related to wave amplitude
        
        This method uses the power expressions to derive efficient stiffness 
        values

        Parameters
        ----------
        omega : float
            angular frequency.
        wavenumber : float
            edge wavenumber.
        wave_DOF : int
            wave degree of freedom.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        float
            efficient wave radiation stiffeness.

        """
        
        
        
        k_plate = self.plate_wavenumber(omega,wave_DOF) 

        phi = np.arccos(wavenumber/k_plate)
        phi[np.isnan(phi)]=0.
        
        sin_phi = np.sin(phi)
        
        if wave_DOF == 1:
            return self.mass_per_area*omega*omega*k_plate*sin_phi        
        elif wave_DOF == 2:
            return self.mass_per_area*omega*omega*k_plate*sin_phi        
        elif wave_DOF in (3,4):
            return 2*self.mass_per_area*omega*omega/k_plate*sin_phi
        else:
            raise ValueError('wave_DOF argument must be in [1,2,3,4]')
        
    def edge_wave_excitation_force(self,omega,wavenumber,wave_DOF,matrix = False):
        """
        wave_excitation_force calculates the blocked force wave amplitudes into edge harmonic diplacement
    
        Parameters
        ----------
        omega : float
            angular frequency.
        wavenumber : float
            edge wavenumber.
        wave_DOF : int
            wave degree of freedom.

        Returns
        -------
        LinearMatrix
            vector [Fx,Fz,Fz,Mx]^T of force
        """
         
        Kx=np.real(wavenumber);
        # Npoint=len(Kx);
 
        # all following wavenumbers are tuned to the positive propagation i.e. 
        # mu is ALLWAYS negative
 
        # use methods from property2structure to get properties
        #eta = self.material.eta
        kB = self.wavenumber_B(omega) #bending wavenumber
        kL = self.wavenumber_L(omega) #longitudinal wavenumber
        kS = self.wavenumber_T(omega) #shear wavenumber
        B = self.B # bending stiffness with damping
        S = self.S #shear stiffness with damping
        nu= self.material.nu #poisson

        f_ = np.zeros((4,1,len(wavenumber)),dtype=np.complex128)

        if wave_DOF == 5: # 0 + 1 together
            f  = self.wave_excitation_force(omega,wavenumber,1,matrix)
            f += self.wave_excitation_force(omega,wavenumber,2,matrix)
            return f
        
        if matrix: # Blocked force is calculated from force and displacment exlicitely
                        
            # In this branch the Langley version is used:
            # f = D.q - F and not analytically as in my version
            
            # Thus, we need the matrix for the force calulation
            q_ = np.zeros((4,1,len(wavenumber)),dtype=np.complex128)
            D_edge = self.edge_radiation_stiffness_wavenumber_LM(omega,wavenumber)
            
            if wave_DOF == 1: # Longitudiunal
                #uL= -np.sqrt(Kx**2-np.complex128(kL)**2)
                uL= self.muL(omega,wavenumber)
                # Prepare vectors and matrices
                f_[0,0,:] = S*2*Kx*uL
                f_[1,0,:] = -S*1j*(2*Kx**2 - kS**2)
                q_[0,0,:] = Kx
                q_[1,0,:] = -1j*uL

            elif wave_DOF == 2: # shear
                #uS= -np.sqrt(Kx**2-np.complex128(kS)**2)
                uS= self.muT(omega,wavenumber)
                f_[0,0,:] = -S*1j*(2*Kx**2 - kS**2)
                f_[1,0,:] = -S*2*Kx*uS
                q_[0,0,:] = -1j*uS
                q_[1,0,:] = -Kx

            elif wave_DOF in (3,4): # bending
                uB2 = -np.sqrt(Kx**2-np.complex128(kB)**2)
                # Prepare vectors and matrices
                f_[2,0,:] = -B*(uB2**3-(2-nu)*Kx**2*uB2)
                f_[3,0,:] = -B*(uB2**2-nu*Kx**2)
                q_[2,0,:] = 1
                q_[3,0,:] = -uB2

            q = mC.LinearMatrix(q_)
            f = mC.LinearMatrix(f_)
            f = D_edge.dot(q)-f_ # claculate force excplicitely
        
        else: # use the analytical solution of the above equation 
          
            if wave_DOF in (1,2):
            
                # in-plane motion, roots
                
                # Try removing not allowed values
                #uL[np.isreal(uL)] = 0.
                #uS[np.isreal(uS)] = 0.
                
                if wave_DOF == 1: # Longitudiunal
                    uL= -np.sqrt(Kx**2-np.complex128(kL)**2)
                    uS= -np.sqrt(Kx**2-np.complex128(kS)**2)
                    # in-plane motion, matrix function elements
                    denom = Kx**2-uL*uS
                    f_[0,0,:] = -2*S*Kx*uL*kS**2/denom # 4*Kx*uL         -2*uL*Kx**2*kS/denom               # Fx - u
                    f_[1,0,:] = 2j*S*uS*uL*kS**2/denom
                elif wave_DOF == 2:
                    uL= -np.sqrt(Kx**2-np.complex128(kL)**2)
                    uS= -np.sqrt(Kx**2-np.complex128(kS)**2)
                    denom = Kx**2-uL*uS
                    # in-plane motion, matrix function elements
                    f_[0,0,:] = 2j*S*(kS**2*uL*uS)/denom                        # Fx - u
                    f_[1,0,:] = 2*S*(kS**2*Kx*uS)/denom
            elif wave_DOF in (3,4):
                # out-of-plane motion, roots
                uB1 = -np.sqrt(Kx**2+kB**2)
                uB2 = -np.sqrt(Kx**2-np.complex128(kB)**2)
                
                #uB2[np.isreal(uB2)] = 0.
    
                # Old version not separating between uB2 + and uB2-
                #data_[2,0,:] = -2*B*(uB2**3+(2-nu)*Kx**2*uB2)                        # Fx - u
                #data_[3,0,:] = -2*B*(uB2**2-nu*Kx**2)
                
                # New version
                #data_[2,0,:] = 2*B*(uB2**2*(uB2+uB1))                        # Fx - u
                f_[2,0,:] = -2*B*(uB1*uB2**2+uB2*uB1**2)                        # Fx - u
                f_[3,0,:] = 2*B*(uB2**2+uB1*uB2)
            
            f = mC.LinearMatrix(f_)
                
        return f

    def edge_wave_excitation_force_cross_correlation(self,omega,wavenumber,wave_DOF,matrix = False):
        """
        Cross correlation of excitation forces due to wave theory
  
        The main purpose of this method is to check if the blocked force
        assumption leads to the same csd as the diffuse field reciprocity
  
        Parameters
        ----------
        omega : float
            angular frequency.
        wavenumber : float
            edge wavenumber.
        wave_DOF : int
            wave degree of freedom.
        matrix : bool
            switch wether the analytical or the direct solution is used- Default value is True-        

        Returns
        -------
        LinearMatrix
            [4 x 4] csd matrix
        """
        
        if wave_DOF == 5:
            Sff = self.wave_excitation_force_cross_correlation(omega,wavenumber,1,matrix) \
                + self.wave_excitation_force_cross_correlation(omega,wavenumber,2,matrix)
        else:
            f = self.wave_excitation_force(omega,wavenumber,wave_DOF,matrix)
            fH = f.H()
            Sff = f.dot(fH) # this is surprisingly simple
                         
        return Sff


    def edge_wave_excitation_displacement(self,omega,wavenumber,wave_DOF):
        """
        Calculates the edge displacement due to outgoing wave
    
        Parameters
        ----------
        omega : float
            angular frequency.
        wavenumber : float
            edge wavenumber.
        wave_DOF : int
            wave degree of freedom.

        Returs
        ------
            vector [4 x 1] of the displacement [u,v,w,\beta] 
        """
         
        Kx=np.real(wavenumber);
        # Npoint=len(Kx);
 
    
 
        # use methods from property2structure to get properties
        #eta = self.material.eta
        kB = self.wavenumber_B(omega) #bending wavenumber
        kL = self.wavenumber_L(omega) #longitudinal wavenumber
        kS = self.wavenumber_T(omega) #shear wavenumber
    
        # Thus, we need the matrix for the force calulation
        q_ = np.zeros((4,1,len(wavenumber)),dtype=np.complex128)
        
        if wave_DOF == 5:
            q =  self.wave_excitation_displacement(omega,wavenumber,1)
            q += self.wave_excitation_displacement(omega,wavenumber,2)
            return q
        
        if wave_DOF == 1: # Longitudiunal
            uL= np.sqrt(Kx**2-np.complex128(kL)**2)
            # Prepare vectors and matrices
            q_[0,0,:] = Kx
            q_[1,0,:] = 1j*uL

        elif wave_DOF == 2:
            uS= np.sqrt(Kx**2-np.complex128(kS)**2)
            q_[0,0,:] = 1j*uS
            q_[1,0,:] = -Kx

        elif wave_DOF in (3,4):
            uB2 = np.sqrt(Kx**2-np.complex128(kB)**2)
            # Prepare vectors and matrices
            q_[2,0,:] = 1
            q_[3,0,:] = uB2

        q = mC.LinearMatrix(q_)
                
        return q

        
    def edge_wave_amplitude_radiated_power(self,Psi,omega,wavenumber,wave_DOF):
        """
        Calculates the radiated power from wave amplidute

        Parameters
        ----------
        Psi : complex
            wave amplitude.
        omega : float
            angular frequency.
        wavenumber : float
            edge wavenuber.
        wave_DOF : int
            wave degrees of freedom.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        float
            radiated power.

        """
        
        
        
        if wave_DOF == 1:
            k_plate = np.real(self.wavenumber_L(omega))
        elif wave_DOF == 2:
            k_plate = np.real(self.wavenumber_T(omega))
        elif wave_DOF in (3,4):
            k_plate = np.real(self.wavenumber_B(omega))
        else:
            raise ValueError('wave_DOF argument must be in [1,2,3 or 4] but is {0}'.format(wave_DOF))

        phi = np.zeros(wavenumber.shape)
        ix = wavenumber <= k_plate
        phi[ix] = np.arccos(wavenumber[ix]/k_plate)
        #phi[np.isnan(phi)]=0.
        
        sin_phi = np.sin(np.real(phi))
        
        Psi2 = np.abs(Psi)**2
        
        if wave_DOF == 1:
            return 0.5*self.mass_per_area*omega**3*k_plate*Psi2*sin_phi        
        elif wave_DOF == 2:
            return 0.5*self.mass_per_area*omega**3*k_plate*Psi2*sin_phi        
        elif wave_DOF in (3,4):
            return self.mass_per_area*omega**3*Psi2/k_plate*sin_phi
        
       




        
    
    
    
    
    
