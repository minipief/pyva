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

# contants
PRESSURE = dof.DOFtype(typestr='pressure')
VELOCITY = dof.DOFtype(typestr='velocity')
STRESS   = dof.DOFtype(typestr='stress')


def I_solid_fluid(xdata,left_ID,right_ID):
    """
    Provide contact matrix I of solid - fluid intefaces.
    
    The matrix corresponds to Eq. (11.71) in [All2008_]
    
    When solid and fluid are exchanged J and I are interchanged. 
    The response DOFs stay the same, because they defined the row index in the
    global allard matrix.

    Parameters
    ----------
    xdata : DataAxis
        xdata for contact matrix definition, wavenumber or frequency.
    left_ID : integer
        left ID of connection.
    right_ID : integer
        right ID of connection.

    Returns
    -------
    DynamicMatrix
        I matrix of solid fluid connection.

    """
    #
    res_dof = solid_fluid_res_dof(left_ID)
    exc_dof = solid_exc_dof(right_ID)
    
    #data = np.zeros((3,4,xdata.Nxdata))
    data = np.array([[0,1,0,0],
                     [0,0,1,0],
                     [0,0,0,1]])
    # stack data along depth dimension
    data = np.dstack([data for _ in range(len(xdata))]) 
    
    return mC.DynamicMatrix(data, xdata, exc_dof, res_dof)

def J_solid_fluid(xdata,left_ID,right_ID):
    """
    Provide contact matrix JU of solid - fluid intefaces.
    
    The matrix corresponds to Eq. (11.71) in [All2008_]
    When solid and fluid are exchanged J and I are interchanged. 
    The response DOFs stay the same, because they defined the row index in the
    global allard matrix.

    Parameters
    ----------
    xdata : DataAxis
        xdata for contact matrix definition, wavenumber or frequency.
    left_ID : integer
        left ID of connection.
    right_ID : integer
        right ID of connection.

    Returns
    -------
    DynamicMatrix
        J matrix of solid fluid connection.

    """
    #
    res_dof = solid_fluid_res_dof(left_ID)
    exc_dof = fluid_exc_dof(right_ID)
    
    #data = np.zeros((3,4,xdata.Nxdata))
    data = np.array([[0,-1],
                     [1, 0],
                     [0, 0]])
    # stack data along depth dimension
    data = np.dstack([data for _ in range(len(xdata)) ] ) 
    
    return mC.DynamicMatrix(data, xdata, exc_dof, res_dof)

def Y_fluid_fixed(xdata,left_ID):
    
    #
    res_dof = dof.DOF([left_ID+1],[3],[VELOCITY])
    exc_dof = fluid_exc_dof(left_ID)
    
    #data = np.zeros((3,4,xdata.Nxdata))
    data = np.array([[0,1]])
    # stack data along depth dimension
    data = np.dstack([data for _ in range(len(xdata)) ] ) 
    
    return mC.DynamicMatrix(data, xdata, exc_dof, res_dof)

def Y_solid_fixed(xdata,left_ID):
    
    #
    res_dof = dof.DOF([left_ID+1,left_ID+1],[1,3],[VELOCITY,VELOCITY])
    exc_dof = solid_exc_dof(left_ID)
    
    #data = np.zeros((3,4,xdata.Nxdata))
    data = np.array([[1,0,0,0],
                     [0,1,0,0]])
    # stack data along depth dimension
    data = np.dstack([data for _ in range(len(xdata)) ] ) 
    
    return mC.DynamicMatrix(data, xdata, exc_dof, res_dof)

def I_fluid_fluid(xdata,left_ID):
    """
    Provide contact matrix I of fluid - fluid intefaces.
    
    Parameters
    ----------
    xdata : DataAxis
        xdata for contact matrix definition, wavenumber or frequency.
    left_ID : integer
        left ID of connection.
    right_ID : integer
        right ID of connection.

    Returns
    -------
    DynamicMatrix
        I matrix of solid fluid connection.

    """
    #
    res_dof = fluid_fluid_res_dof(left_ID)
    exc_dof = fluid_fluid_res_dof(left_ID)
    
    #data = np.zeros((3,4,xdata.Nxdata))
    data = np.eye(2)
    # stack data along depth dimension
    data = np.dstack([data for _ in range(len(xdata))]) 
    
    return mC.DynamicMatrix(data, xdata, exc_dof, res_dof)

def J_fluid_fluid(xdata,left_ID,right_ID):
    """
    Provide contact matrix J of fluid - fluid intefaces.
    
    Parameters
    ----------
    xdata : DataAxis
        xdata for contact matrix definition, wavenumber or frequency.
    left_ID : integer
        left ID of connection.
    right_ID : integer
        right ID of connection.

    Returns
    -------
    DynamicMatrix
        I matrix of solid fluid connection.

    """
    #
    res_dof = fluid_fluid_res_dof(left_ID)
    exc_dof = fluid_fluid_res_dof(right_ID)
    
    #data = np.zeros((3,4,xdata.Nxdata))
    data = -np.eye(2)
    # stack data along depth dimension
    data = np.dstack([data for _ in range(len(xdata))]) 
    
    return mC.DynamicMatrix(data, xdata, exc_dof, res_dof)

def I_porous_fluid(xdata,porosity,ID):
    
    #res_dof = dof()
    phi = porosity
    return np.array([[0,1-phi,phi,0,0,0],
                     [0,  0  , 0 ,1,0,0],
                     [0,  0  , 0 ,0,1,0],
                     [0,  0  , 0 ,0,0,1]])

def J_porous_fluid(porosity):
    
    phi = porosity
    return np.array([[0    ,-1],
                     [1-phi, 0],
                     [0    , 0],
                     [0    , 0]])

# functions for DOF generation

def solid_fluid_res_dof(left_ID):
    
    ID_   = [left_ID]*3
    dof_  = [3,3,1]
    type_ = [ VELOCITY , STRESS, STRESS ] 
    return dof.DOF(ID_,dof_,type_)

def solid_exc_dof(left_ID):
    
    ID_   = [left_ID]*4 
    dof_  = [1,3,3,1]
    type_ = [VELOCITY, VELOCITY , STRESS, STRESS] 
        
    return dof.DOF(ID_,dof_,type_)

# def fluid_solid_dof(left_ID):
    
#     ID_   = [left_ID]*2
#     dof_  = [0,3]
#     type_ = [ PRESSURE, VELOCITY ] 
#     return dof.DOF(ID_,dof_,type_)

def fluid_exc_dof(left_ID):
    
    ID_   = [left_ID]*2
    dof_  = [0,3]
    type_ = [ PRESSURE, VELOCITY] 
    return dof.DOF(ID_,dof_,type_)

def fluid_fluid_res_dof(left_ID):

    ID_   = [left_ID]*2
    dof_  = [0,3]
    type_ = [ PRESSURE, VELOCITY] 
    return dof.DOF(ID_,dof_,type_)

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
        self.type      = 'equivalent_fluid' 
     
    @staticmethod   
    def get_xdata(omega,kx):
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

    @property
    def isequivalentfluid(self):
        """
        Determine if layer is of type equivalent fluid

        Defauls parameter is True

        Returns
        -------
        bool
            True.

        """
        return True
    


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
    
    def __init__(self,thickness,fluid=mc.Fluid):
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
        # The ID is used only generically as this is set by layer position
        _left_dof = dof.DOF([0,0],[0,3],Tdof)
        _right_dof = dof.DOF([1,1],[0,3],Tdof)
                       
        super().__init__(thickness,_left_dof,_right_dof)
        self.fluid = fluid       # fluid in the layer
        

    def __str__(self):
        str_ = 'fluid layer of thickness {0}'.format(self.thickness)
        return str_
        
    def transfer_impedance(self,omega,kx=0,ID=[1,2]):
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
    
    def __init__(self,thickness,fluid=mc.Fluid):
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
        
    def transfer_impedance(self,omega,kx = 0, ID = [1,2]):
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
    
    def __init__(self,thickness,rho):
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
        _left_dof = dof.DOF([0,0],[0,3],Tdof)
        _right_dof = dof.DOF([1,1],[0,3],Tdof)
                       
        super().__init__(thickness,_left_dof,_right_dof)
        
        self.rho = rho
        
        # create lumped element component
        self._lumped = ac1D.MassStiffness(rho*thickness,0)
        
    def __str__(self):
        str_ = 'mass layer of thickness  {0}'.format(self.thickness)
        return str_

    def transfer_impedance(self,omega,kx = 0,ID=[1,2]):
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
    
    def __init__(self,plate_prop):
        """
        Class constructor for PlateLayer objects
        
        Parameters
        ----------
        plate_prop : PlateProp
            plate property of layer
        ID : list or tuple
            node IDs of left and right side
        """                
        # set DOF according to ID and natural DOF of the layer
        Tdof = ('pressure','velocity')
        _left_dof = dof.DOF([0,0],[0,3],Tdof)
        _right_dof = dof.DOF([1,1],[0,3],Tdof)
                       
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

    def transfer_impedance(self,omega,kx=0,ID = [1,2]):
        """
        Transferimpedance of plates
        
        Takes simple infinite plate model to calculate the transfer impedance
        
        Parameters
        ----------
        omega : float or ndarray
            scalar angular frequency
        kx : float or ndarray, optional
            In-plane wavenumber. The default is 0.
        ID : list of int,
            Left and right ID. The default is [1,2].
                      
        Returns
        -------
        DynamicMatrix
            [2 x 2] array of transferimpedance
        
        """
                
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
        
    def __init__(self,thickness,hole_radius,fluid=visc_air,pattern='square',alpha=2.,**kwargs):
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
        _left_dof = dof.DOF([0,0],[0,3],Tdof)
        _right_dof = dof.DOF([1,1],[0,3],Tdof)
         
        # assign attributes of mother classes
        super().__init__(thickness,_left_dof,_right_dof)
        self._lumped = ac1D.PerforatedLayer(thickness,hole_radius,1.,fluid=fluid, \
                                            pattern = pattern,alpha = alpha,**kwargs)
        
    def transfer_impedance(self,omega,kx = 0,ID=[1,2]):
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
        
                    
class SolidLayer(AcousticLayer):
    """
    Class for Modelling elastic solid material as InfiniteLayer.
    
    Attributes
    ----------
    thickness: float
        thickness of the perforted layer
    prop: PlateProp
        property of plate modelled as 3D layer
    """
    
    def __init__(self,plate_prop):
        """
        Class contructor for SolidLayer objects.

        Parameters
        ----------
        plate_prop : PlateProp
            plate property of layer
        ID : list or tuple
            node IDs of left and right side
        """
                
        # set DOF according to ID and natural DOF of the layer
        Tdof = ('velocity','velocity','stress','stress')
        _left_dof = dof.DOF([0,0,0,0],[1,3,3,1],Tdof)
        _right_dof = dof.DOF([1,1,1,1],[1,3,3,1],Tdof)
                       
        super().__init__(plate_prop.thickness,_left_dof,_right_dof)
        
        self.prop = plate_prop
        self.type = 'solid'
    
    def __repr__(self):
        """
        repr for SolidLayer

        Returns
        -------
        str_ : TYPE
            DESCRIPTION.

        """
        
        str_ = 'solid layer of thickness {0}'.format(self.thickness)
        return str_

    @property
    def isequivalentfluid(self):
        """
        Determine if layer is of type equivalent fluid

        Defauls parameter is False

        Returns
        -------
        bool
            False.

        """
        return False

    @property
    def mass_per_area(self):
        """
        Mass per area of SolidLayer.
        
        Returns
        -------
        float
            mass per area.

        """        
        return self.prop.mass_per_area
    
    def transfer_impedance(self,omega,kx=0, ID = [1,2], allard = False):
        """
        Transferimpedance of Solidlayer
        
        Implementation according to [1] except the not used denominator D1+D2*kx.
        [1] A. Uthayasuriyan, Advanced vibro-acoustic condensed models for multi-layer structures,
        PhD Thesis, University de Lyon, Lyon, France, 2021.

        Parameters
        ----------
        omega : float or ndarray
            scalar angular frequency
        kx : float or ndarray, optional
            In-plane wavenumber. The default is 0.
        ID : list of int
            Left and right ID. The default is [1,2].
                      
        Returns
        -------
        DynamicMatrix
            [2 x 2] array of transferimpedance         
        
        """
                    
        xdata = self.get_xdata(omega, kx)
                       
        # Bulk wavenmumber L(ongitudinal) and S(hear)
        k_L = self.prop.material.wavenumber_L(omega) # \delta_1^2 Allard
        k_S = self.prop.material.wavenumber_S(omega) # \delta_3^2 Allard
        k_L2 = k_L*k_L   
        k_S2 = k_S*k_S       
        
        # wavenumber in z for L and S
        kLz2 = k_L2-kx**2
        kSz2 = k_S2-kx**2
        kLz = np.sqrt(kLz2) #k_13
        kSz = np.sqrt(kSz2) #k_33

        #rho_om  = lambda omega_: self.prop.fluid.rho_freq(omega_)*omega_
        h  = self.thickness
        
        #lame1 = self.prop.material.lambda_lame
        mu    = self.prop.material.G_complex
        #mu    = self.prop.material.G
        D1  = mu*(kSz2-kx**2)
        #D1  = self.prop.material.lambda_lame()*(kx**2+kLz2)+2*mu*kLz2
        D2  = 2*mu*kx
        
        # if allard:

        # Gamma = np.zeros((len(xdata),4,4),dtype = np.complex128)
        # h = -h
        # # 1st row
        # Gamma[:,0,0] = omega*kx*np.cos(kLz*h)
        # Gamma[:,0,1] = -1j*omega*kx*np.sin(kLz*h)
        # Gamma[:,0,2] = 1j*omega*kSz*np.sin(kSz*h)
        # Gamma[:,0,3] = -omega*kSz*np.cos(kSz*h)
        # # 2nd row
        # Gamma[:,1,0] = -1j*omega*kLz*np.sin(kLz*h)
        # Gamma[:,1,1] = omega*kLz*np.cos(kLz*h)
        # Gamma[:,1,2] = omega*kx*np.cos(kSz*h)
        # Gamma[:,1,3] = -1j*omega*kx*np.sin(kSz*h)
        # # 3rd row
        # Gamma[:,2,0] = -D1*np.cos(kLz*h)
        # Gamma[:,2,1] = 1j*D1*np.sin(kLz*h)
        # Gamma[:,2,2] = 1j*D2*kSz*np.sin(kSz*h)
        # Gamma[:,2,3] = -D2*kSz*np.cos(kSz*h)
        # # 4th row
        # Gamma[:,3,0] = 1j*D2*kLz*np.sin(kLz*h)
        # Gamma[:,3,1] = -D2*kLz*np.cos(kLz*h)
        # Gamma[:,3,2] = D1*np.cos(kSz*h)
        # Gamma[:,3,3] = -1j*D1*np.sin(kSz*h)

       
        # Gamma0_inv = np.zeros((len(xdata),4,4),dtype = np.complex128) 



        # Gamma0_inv[:,0,0] = 2*kx/omega/k_S**2
        # Gamma0_inv[:,0,2] = -1/mu/k_S**2
        # Gamma0_inv[:,1,1] = (kSz2-kx**2)/omega/kLz/k_S**2
        # Gamma0_inv[:,1,3] = -kx/mu/kLz/k_S**2
        # Gamma0_inv[:,2,1] = kx/omega/k_S**2
        # Gamma0_inv[:,2,3] = 1/mu/k_S**2
        # Gamma0_inv[:,3,0] = -(kSz2-kx**2)/omega/kSz/k_S**2 # corrected in VAOne docs
        # Gamma0_inv[:,3,2] = -kx/mu/kSz/k_S**2
        
        # Gamma = np.matmul(Gamma,Gamma0_inv)
        # Gamma = np.moveaxis(Gamma,0,-1) # put xaxis last

        # else:
            
        Gamma = np.zeros((4,4,len(xdata)),dtype = np.complex128) 
        A = 1./(D1 + D2*kx) # Denominator not used!
        
        coshkLz = np.cos(kLz*h)
        sinhkLz = np.sin(kLz*h)
        coshkSz = np.cos(kSz*h)
        sinhkSz = np.sin(kSz*h)
      
        
        Gamma[0,0,:] = A*(D1*coshkSz+D2*kx*coshkLz)
        Gamma[3,3,:] = Gamma[0,0,:]
        Gamma[1,1,:] = A*(D1*coshkLz+D2*kx*coshkSz)
        Gamma[2,2,:] = Gamma[1,1,:]
        Gamma[0,1,:] = A*-1j*(D2*kLz*kSz*sinhkSz-D1*kx*sinhkLz)/kLz
        Gamma[2,3,:] = Gamma[0,1,:] 
        Gamma[1,0,:] = A*1j*(D2*kLz*kSz*sinhkLz-D1*kx*sinhkSz)/kSz
        #Gamma[1,0,:] = A*1j*kx*(D2*kSz*sinhkLz-D1*sinhkSz)/kSz # kLz replaced by kx 
        Gamma[3,2,:] = Gamma[1,0,:]
        # Buffer to avoid T31 and T42 singularity for kx = 0
        T13_ = (coshkSz - coshkLz)
        Gamma[0,2,:] = A*T13_*kx*omega
        Gamma[1,3,:] = Gamma[0,2,:]
        Gamma[2,0,:] = A*D1*D2/omega*T13_
        Gamma[3,1,:] = Gamma[2,0,:]
        Gamma[0,3,:] = A*-1j*omega*(kLz*kSz*sinhkSz+kx**2*sinhkLz)/kLz
        Gamma[1,2,:] = A*-1j*omega*(kLz*kSz*sinhkLz+kx**2*sinhkSz)/kSz
        #Gamma[1,2,:] = A*-1j*omega*kx*(kSz*sinhkLz+kx*sinhkSz)/kSz # kLz replaced by kx
        Gamma[2,1,:] = A*-1j*(D2**2*kLz*kSz*sinhkSz+D1**2*sinhkLz)/(omega*kLz)
        Gamma[3,0,:] = A*-1j*(D2**2*kLz*kSz*sinhkLz+D1**2*sinhkSz)/(omega*kSz)
            
            
            
        xdata = self.get_xdata(omega, kx )   
                       
        
        # Update ID in DOF attribute
        ldof = dof.DOF([ID[0],ID[0],ID[0],ID[0]], self.left_dof.dof,self.left_dof.type)
        rdof = dof.DOF([ID[1],ID[1],ID[1],ID[1]], self.right_dof.dof,self.right_dof.type)
        
        return mC.DynamicMatrix(Gamma,xdata,rdof,ldof)