# -*- coding: utf-8 -*-
"""
Module for infinite and flat layers

The infiniteLayers module provides classes for the simulation of 
infinitely extended layer of sequential flat systems
"""

import pyva.properties.materialClasses as mc
import pyva.systems.acoustic1Dsystems as ac1D
import numpy as np
import pyva.data.matrixClasses as mC
import pyva.data.dof as dof
import scipy.special as spl
import pyva.useful as uf
import copy

class AcousticLayer:
    """
    Base class for all acoustic infinite layer
    
    This class constitutes the system properties, of one layer or element 
    of acoustic networks or layer models
    
    The XXXXdof attributes follow from the properties of the daughter classes
    
    Attributes
    ----------
    left_dof : DOF
        degrees of freedom of left connection 
    right_dof : DOF
        degrees of freedom of left connection 
        
    """
    
    def __init__(self,thickness,left_dof,right_dof):
        """
        Constructor of AcousticLayer
       

        Parameters
        ----------
        thickness : float
            Thickness of layer.
        left_dof : DOF
            degrees of freedom of left connection 
        right_dof : DOF
            degrees of freedom of left connection 

        Returns
        -------
        None.
        """        


        
        if isinstance(left_dof,dof.DOF):
            self.left_dof = left_dof
        else:
            raise(ValueError,"leftdof argument must be an instance of the dof class")

        if isinstance(right_dof,dof.DOF):
            self.right_dof = right_dof
        else:
            raise(ValueError,"rightdof argument must be an instance of the dof class")
            
        self.thickness = thickness
        
    def get_xdata(self,omega,kx):
        """
        Determines the appropriate xdata from omega and kx.
        
        Infinite Layer theory often involves integration over wavenumber.
        Thus, the kx is the integration variable.
        
        If omega is scalar xdata will be wavenumber in x (kx)
        if kx is scalar xdata will be angular frequency
        if kx and omega have same dimension it is assumed that kx belongs to a constant angle 
        and is given by kx = omega/c0*sin(theta) and xdata will be angular frequency
        
        Parameters
        ----------
        omega : float or ndarray
            angular frequency
        kx :float or ndarray
             wavenumber
            
        Returns
        -------
        DataAxis
            wavenumber or omega as DataAxis object 
            
        """
        
        if uf.isscalar(omega):
            xdata  = mC.DataAxis(np.array(kx).flatten(),typestr = 'wavenumber')
        elif uf.isscalar(kx):
            xdata  = mC.DataAxis(np.array(omega).flatten(),typestr='angular frequency')
        elif omega.shape == kx.shape:
            xdata  = mC.DataAxis(np.array(omega).flatten(),typestr='angular frequency')
        else:
            raise ValueError('One argument from kx or omega must be scalar or both are of same shape')
            
        return xdata


class FluidLayer(AcousticLayer):
    """
    The AcousticLayer class deals with the infinite layers of fluid
    
    Attributes
    ----------
    thickness : float
        thickness of the layer
    fluid : fluid, optional
        Fluid of the layer. The default is mc.Fluid.    
    ID : list of int
        [left ID, right ID] of fluid layer
    """
    
    def __init__(self,thickness,fluid=mc.Fluid,ID = [1,2]):
        """
        Constructor of FluidLayer

        Parameters
        ----------
        thickness : float
            thickness of the layer
        fluid : fluid, optional
            Fluid of the layer. The default is mc.Fluid.
        ID : list of int
            [left ID, right ID] of fluid layer

        Returns
        -------
        None.

        """
        
        
        # set DOF according to ID and natural DOF of the layer
        Tdof = ('pressure','velocity')
        _left_dof = dof.DOF([ID[0], ID[0]],[0,3],Tdof)
        _right_dof = dof.DOF([ID[1], ID[1]],[0,3],Tdof)
                       
        super().__init__(thickness,_left_dof,_right_dof)
        self.fluid = fluid       # fluid in the layer
        

    def __str__(self):
        str_ = 'fluid layer of thickness {0}'.format(self.thickness)
        return str_
        
    def transfer_impedance(self,omega,kx=0, ID = None):
        """
        Transferimpedance of fluid layer

        Parameters
        ----------
        omega : float or ndarray
            scalar angular frequency
        kx : float or ndarray, optional
            In-plane wavenumber. The default is 0.
        ID : list of int, optional
            Left and right when overwritten, None takes object ID. The default is None.
                      
        Returns
        -------
        DynamicMatrix
            [2 x 2] array of transferimpedance         
        
        """
        
        if ID is None:
            ID = self.ID

        kz = lambda omega_,kx_: np.sqrt(self.fluid.wavenumber(omega_)**2-kx_**2)
        rho_om  = lambda omega_: self.fluid.rho_freq(omega_)*omega_
        h  = self.thickness
        
        # Old - not correctr version
        #T11 = lambda omega_,kx_: np.cos(kz(omega_,kx_)*h)
        #T12 = lambda omega_,kx_: 1j*np.sin(kz(omega_,kx_)*h)*z(omega_)
        #T21 = lambda omega_,kx_: 1j*np.sin(kz(omega,kx_)*h)/z(omega_)
        # new
        T11 = lambda omega_,kx_: np.cos(kz(omega_,kx_)*h)
        T12 = lambda omega_,kx_: 1j*np.sin(kz(omega_,kx_)*h)*rho_om(omega_)/kz(omega_,kx_)
        T21 = lambda omega_,kx_: 1j*np.sin(kz(omega,kx_)*h)*kz(omega_,kx_)/rho_om(omega_)
        
        xdata = self.get_xdata(omega, kx )   
                       
        data   = np.zeros((2,2,len(xdata)),dtype = np.complex128)

        data[0,0,:] = T11(omega,kx)
        data[1,1,:] = T11(omega,kx) # T11 = T22
        data[0,1,:] = T12(omega,kx)
        data[1,0,:] = T21(omega,kx)
        
        # Update ID in DOF attribute
        ldof = dof.DOF([ID[0],ID[0]], self.left_dof.dof,self.left_dof.type)
        rdof = dof.DOF([ID[1],ID[1]], self.right_dof.dof,self.right_dof.type)
        
        return mC.DynamicMatrix(data,xdata,rdof,ldof)

class FluidLayerHoneyComb(FluidLayer):
    """
    The FluidLayer class deals with the infinite honcomb core layers in the absorber context
    
    Implementation is equal to fluid layer, except that the wavenumbe in in-plane direction is alway zero.
    
    Attributes
    ----------
    thickness : float
        thickness of the layer
    fluid : fluid, optional
        Fluid of the layer. The default is mc.Fluid.    
    ID : list of int
        [left ID, right ID] of fluid layer
    
    """
    
    def __init__(self,thickness,fluid=mc.Fluid,ID=[1,2]):
        """
        Class contructor for acoustic tube
        
        Parameters
        ----------
        thickness : float
            thickness of the layer
        fluid : fluid, optional
            Fluid of the layer. The default is mc.Fluid.
        ID : list of int
            [left ID, right ID] of fluid layer

        Returns
        -------
        None.
                        
        """
        
        super().__init__(thickness,fluid)
        

    def __str__(self):
        str_ = 'fluid layer Honeycomb of thickness {0}'.format(self.thickness)
        return str_
        
    def transfer_impedance(self,omega,kx = 0, ID = None):
        """
        Transferimpedance of honey comb layers
        
        Takes simple 1D transmission line model based on AcousticTube properties.
        kx Parameter is ignored in the transfermatrix but used in the xdata generation
        
        Parameters
        ----------
        omega : float or ndarray
            scalar angular frequency
        kx : float or ndarray, optional
            In-plane wavenumber. The default is 0.
        ID : list of int, optional
            Left and right when overwritten, None takes object ID. The default is None.
                      
        Returns
        -------
        DynamicMatrix
            [2 x 2] array of transferimpedance         
        
        """
        
        if ID is None:
            ID = self.ID

        kz = lambda omega_: np.sqrt(self.fluid.wavenumber(omega_)**2)
        rho_om  = lambda omega_: self.fluid.rho_freq(omega_)*omega_
        h  = self.thickness
        
        
        # Old - not correctr version
        #T11 = lambda omega_,kx_: np.cos(kz(omega_,kx_)*h)
        #T12 = lambda omega_,kx_: 1j*np.sin(kz(omega_,kx_)*h)*z(omega_)
        #T21 = lambda omega_,kx_: 1j*np.sin(kz(omega,kx_)*h)/z(omega_)
        # new
        T11 = lambda omega_: np.cos(kz(omega_)*h)
        T12 = lambda omega_: 1j*np.sin(kz(omega_)*h)*rho_om(omega_)/kz(omega_)
        T21 = lambda omega_: 1j*np.sin(kz(omega_)*h)*kz(omega_)/rho_om(omega_)
        
        xdata = self.get_xdata(omega, kx)
        
                       
        data   = np.zeros((2,2,len(xdata)),dtype = np.complex128)

        # if xdata is frequency this works, if omega is scalar, value is written in array
        data[0,0,:] = T11(omega)
        data[1,1,:] = T11(omega) # T11 = T22
        data[0,1,:] = T12(omega)
        data[1,0,:] = T21(omega)
        
        # Update ID in DOF attribute
        ldof = dof.DOF([ID[0],ID[0]], self.left_dof.dof,self.left_dof.type)
        rdof = dof.DOF([ID[1],ID[1]], self.right_dof.dof,self.right_dof.type)
        
        return mC.DynamicMatrix(data,xdata,rdof,ldof)

    
class MassLayer(AcousticLayer):
    """
    The massAcoustic class represents the mass layer 
    
    Attributes
    ----------
    thickness : float
        thickness of the layer
    rho : float
        Density.    
    ID : list of int
        [left ID, right ID] of MassLayer
    """
    
    def __init__(self,thickness,rho,ID=[1,2]):
        """
        Constructror of mass layer

        Parameters
        ----------
        thickness : float
            thickness of the layer.
        rho : float
            Density.
        ID : list of int
            [left ID, right ID] of MassLayer

        Returns
        -------
        None.

        """
                
        # set DOF according to ID and natural DOF of the layer
        Tdof = ('pressure','velocity')
        _left_dof = dof.DOF([ID[0], ID[0]],[0,3],Tdof)
        _right_dof = dof.DOF([ID[1], ID[1]],[0,3],Tdof)
                       
        super().__init__(thickness,_left_dof,_right_dof)
        
        self.rho = rho
        
        # create lumped element component
        self._lumped = ac1D.MassStiffness(rho*thickness,0)
        
    def __str__(self):
        str_ = 'mass layer of thickness  {0}'.format(self.thickness)
        return str_

    def transfer_impedance(self,omega,kx = 0,ID=None):
        """
        Transferimpedance of honey mass layer
        
        Takes simple 1D transmission line model based on AcousticTube properties.
        kx Parameter is ignored in the transfermatrix but used in the xdata generation
        
        Parameters
        ----------
        omega : float or ndarray
            scalar angular frequency
        kx : float or ndarray, optional
            In-plane wavenumber. The default is 0.
        ID : list of int, optional
            Left and right when overwritten, None takes object ID. The default is None.
                      
        Returns
        -------
        DynamicMatrix
            [2 x 2] array of transferimpedance    |
        
        """

        if ID is None:
            ID = self.ID


        # Update ID in DOF attribute
        ldof = dof.DOF([ID[0],ID[0]], self.left_dof.dof,self.left_dof.type)
        rdof = dof.DOF([ID[1],ID[1]], self.right_dof.dof,self.right_dof.type)
        
        # convert parameters to np.array in any case        
        kx    = np.array(kx).flatten()
        omega = np.array(omega).flatten()
        
        xdata = self.get_xdata(omega, kx)
        
        if xdata.type.typestr == 'wavenumber':
            # xdata ,ust be replaced for original mehod. 
            _T = self._lumped.transfer_impedance(omega,ID,velocity = 'v')
            # just matrix repetition, there is no wavenumber dependence
            T_kx = np.tile(_T.data,(xdata.shape))            
            return mC.DynamicMatrix(T_kx,xdata,rdof,ldof)
            
        elif xdata.type.typestr == 'angular frequency':
            # method from 1Dsystem can be used as is
            return self._lumped.transfer_impedance(omega,ID,velocity = 'v',DOF=[0,3])

    @property
    def mass_per_area(self):
        """
        Mass per area of MassLayer
        
        Returns
        -------
        float
            mass per area.

        """
        return self.rho*self.thickness
    
    def transmission_coefficient(self,omega,theta = 0.,fluid = mc.Fluid()):
        """
        Mass law transmission coefficient 
        
        Parameters
        ----------
        omega : float
            angular frequency.
        theta : float, optional
            angle of incidence related to normal. The default is 0..
        fluid : fluid, optional
            surrounding fluid. The default is mc.Fluid().

        Returns
        -------
        float
            transmission coefficient.

        """
                     
        return 1/(1+(self.mass_per_area*omega/2/fluid.z0*np.cos(theta))**2)
        
                      

class PlateLayer(AcousticLayer):
    """
    The PlateLayer class represents the infinite plate layer 
    
    Attributes:
      plate_prop: of plate
    """
    
    def __init__(self,plate_prop,ID=[1,2]):
        """
        Class contructor for PlateLayer objects
        
        Args:
             impedance: transfer impedance as complex function of omega @type function
                  area: cross sectio of the element   
                        
        """
                
        # set DOF according to ID and natural DOF of the layer
        Tdof = ('pressure','velocity')
        _left_dof = dof.DOF([ID[0], ID[1]],[0,3],Tdof)
        _right_dof = dof.DOF([ID[1], ID[1]],[0,3],Tdof)
                       
        super().__init__(plate_prop.thickness,_left_dof,_right_dof)
        
        self.prop = plate_prop
    
    def __repr__(self):
        str_ = 'plate layer of thickness {0}'.format(self.thickness)
        return str_

    @property
    def mass_per_area(self):
        """
        Mass per area of PlateLayer
        
        Returns
        -------
        float
            mass per area.

        """
        
        return self.prop.mass_per_area

    def transfer_impedance(self,omega,kx=0,ID = None):
        """
        Transferimpedance of plates
        
        Takes simple infinite plate model to calculate the transfer impedance
        
        Parameters
        ----------
        omega : float or ndarray
            scalar angular frequency
        kx : float or ndarray, optional
            In-plane wavenumber. The default is 0.
        ID : list of int, optional
            Left and right when overwritten, None takes object ID. The default is None.
                      
        Returns
        -------
        DynamicMatrix
            [2 x 2] array of transferimpedance
        
        """

        if ID is None:
            ID = self.ID
                
        xdata = self.get_xdata(omega, kx)
                       
        data   = np.zeros((2,2,len(xdata)),dtype = np.complex128)

        data[0,0,:] = 1.
        data[1,1,:] = 1.
        data[0,1,:] = self.prop.transfer_impedance(omega,kx)
        #data[1,0,:] = 0.
        
        # Update ID in DOF attribute
        ldof = dof.DOF([ID[0],ID[0]], self.left_dof.dof,self.left_dof.type)
        rdof = dof.DOF([ID[1],ID[1]], self.right_dof.dof,self.right_dof.type)
        
        return mC.DynamicMatrix(data,xdata,rdof,ldof)
    
    def transmission_coefficient(self,omega,theta = 0.,fluid = mc.Fluid()):
        """
        Mass law transmission coefficient

        Parameters
        ----------
        omega : float
            angular frequency.
        theta : float, optional
            angle of incidence related to normal. The default is 0..
        fluid : fluid, optional
            surrounding fluid. The default is mc.Fluid().

        Returns
        -------
        float
            transmission coefficient.
        
        """
        
        ka = fluid.wavenumber(omega)
        kB = self.prop.wavenumber_B(omega)
        
        return 1/(1+(self.mass_per_area*omega/2/fluid.z0*((ka*np.sin(theta)/kB)**4-1)*np.cos(theta))**2)
    
class PerforatedLayer(AcousticLayer):
    """
    The PerforatedLayer class deals perforated layers in acoustic networks
    
    The behaviour is calculated using flow throught the holes and radiation be the disk radiator pattern
    
    Attributes
    ----------
    thickness: float
        thickness of the perforted layer
    pattern: str
        identifier for hole pattern (quadradic, triangular, )
    hole_radius: float
        radius of the holes in the membrane
    porosity: float 
        area ratio of holes to totalsurface
    fluid: fluid
        fluid in the hole
    alpha: float
        correction constant for resistivity correction approx 4.-2. for
        sharp to round edges
    """
   
    
    # Model is based on fluid viscosity
    visc_air = mc.Fluid(damping_model='viscous')
        
    def __init__(self,thickness,hole_radius,fluid=visc_air,ID=[1,2],pattern='square',alpha=2.,**kwargs):
        """
        Class contructor for PerforatedLayer
        

        Parameters
        ----------
        thickness: float
            thickness of the perforted layer
        hole_radius: float
            radius of the holes in the membrane
        fluid: fluid
            fluid in the hole. The default is visc_air.
        ID : TYPE, optional
            DESCRIPTION. The default is [1,2].
        pattern: str
            identifier for hole pattern ('square',triangular'',). The default is 'square'.
        alpha: float
            correction constant for resistivity correction approx 4.-2. for
            sharp to round edges: The default is 2.
        **kwargs : dict
            arbitrary keyword parameter list
        distance: float
            distance between holes
        porosity:  float 
            surface porosity of perforate

        Returns
        -------
        None.
                        
        """
 
        # set DOF according to ID and natural DOF of the layer
        Tdof = ('pressure','velocity')
        _left_dof = dof.DOF([ID[0], ID[1]],[0,3],Tdof)
        _right_dof = dof.DOF([ID[1], ID[1]],[0,3],Tdof)
         
        # assign attributes of mother classes
        super().__init__(thickness,_left_dof,_right_dof)
        self._lumped = ac1D.PerforatedLayer(thickness,hole_radius,1.,fluid=fluid, \
                                            pattern = pattern,alpha = alpha,**kwargs)
        
    def transfer_impedance(self,omega,kx = 0,ID=[None]):
        """
        Transferimpedance of perforateLayer
        
        The transferimpedance is caluclated based on the lumped_acoustic
        model, because it is independent from the kx
        
        Parameters
        ----------
        omega : float or ndarray
            scalar angular frequency
        kx : float or ndarray, optional
            In-plane wavenumber. The default is 0.
        ID : list of int, optional
            Left and right when overwritten, None takes object ID. The default is None.
                      
        Returns
        -------
        DynamicMatrix
            [2 x 2] array of transferimpedance
        """

        # Update ID in DOF attribute
        
        if ID[0] is None:
            _ID = [self.left_dof.ID[0],self.right_dof.ID[1]]
        else:
            _ID = ID
            
        ldof = dof.DOF([_ID[0],_ID[0]], self.left_dof.dof,self.left_dof.type)
        rdof = dof.DOF([_ID[1],_ID[1]], self.right_dof.dof,self.right_dof.type)
        
        # convert parameters to np.array in any case        
        kx    = np.array(kx).flatten()
        omega = np.array(omega).flatten()
        
        xdata = self.get_xdata(omega, kx)
        
        if xdata.type.typestr == 'wavenumber':
            # xdata ,ust be replaced for original mehod. 
            _T = self._lumped.transfer_impedance(omega,ID,velocity = 'v')
            # just matrix repetition, there is no wavenumber dependence
            T_kx = np.tile(_T.data,(xdata.shape))            
            return mC.DynamicMatrix(T_kx,xdata,rdof,ldof)
            
        elif xdata.type.typestr == 'angular frequency':
            # method from 1Dsystem can be used as is
            return self._lumped.transfer_impedance(omega,_ID,velocity = 'v',DOF=[0,3])
        
    def plot(self,omega,**kwargs):
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
                
        self._lumped.plot(omega,**kwargs)
        
    @property
    def porosity(self):
        """
        porosity parameter

        Returns
        -------
        None.

        """
        
        return self._lumped.porosity
        
                    
