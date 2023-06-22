# -*- coding: utf-8 -*-
"""
Module for deterministic, random and hybrid models
"""
import scipy.linalg as linalg
import scipy.integrate as integrate
import scipy.special as scl
import numpy as np

import pyva.data.matrixClasses as mC
import pyva.properties.materialClasses as matC
import pyva.loads.loadCase as lC
import pyva.data.dof as dof
import pyva.systems.infiniteLayers as iL
import pyva.systems.acousticRadiators as aR
import pyva.useful as uf

import copy


class FEM:
    
    """
    The FEM class is dedictated to deterministic discrete systems.
    This class manages the degrees of freedom of such systems and the matrix that
    describes the system.

    This class is the least mature class in this toolbox, because it was mainly 
    created fot training purpose and the examples from the book.
    """
    
    def __init__(self,ID,mesh,modes,damping_loss = 0.01,**kwargs):
        """
        Constructor of finite element models (FEM)

        Parameters
        ----------
        ID : int
            identifier of FE model.
        mesh : mesh
            mesh of FE model.
        modes : shape
            modes of the FE model
        damping_loss : Signal or float, optional
            global damping loss of FE model. The default is 0.01.
        **kwargs : Arbitrary keyword arguments
            DESCRIPTION.

        Returns
        -------
        None.

        """
         
        self.ID = ID
        self.mesh = mesh
        self.modes = modes
        self._damping_loss = damping_loss
        self._loads = {}
        self._results = {}
                    
        for kw in kwargs:
            if kw == 'loads':
                 self._loads.update = kwargs[kw]
            else:
                 print('Unkown argument')

    def __repr__(self):
        _str = 'FEM object\n with loads \n'
        for lkey in self._loads.keys():
            _str += "ID=" + str(lkey) + str(self._loads[lkey]) + "\n"
            
        for rkey in self._results.keys():
            _str += "ID=" + str(rkey) + str(self._results[rkey]) + "\n"
            
        return _str
    
    def add_load(self,load):
        """
        Adds load to FEM objects

        Parameters
        ----------
        load : dict or Load
            load of FE model

        Returns
        -------
        None.

        """
        if isinstance(load,(dict)):
            self._loads.update(load)
        elif isinstance(load,lC.Load):
            self._loads.update({load.name : load })   

    def add_result(self,result):
        """
        Add result to FE model

        Parameters
        ----------
        result : Signal
            result from solve.

        Returns
        -------
        None.

        """
   
        self._results.update(result)
        
    def plot(self,nfig=1,loadID=1,**kwargs):
        """
        plot results

        Parameters
        ----------
        nfig : int, optional
            figure identifier. The default is 1.
        loadID : int, optional
            identifier of load case. The default is 1.
        **kwargs : dict
             Arbitrary keyword arguments passed to plot method of Signal

        Returns
        -------
        None.

        """
        
        _res = self._results[loadID]
        _res.plot(nfig,**kwargs)
        
    def modal_dynamic_stiffness_diag(self,omega):
        """
        Diagonal dynamic stiffness matrix in modal coordinates

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        om_n2 : ndarray
            diagonal of modal dynamic stiffness.

        """
        
        # arrange modes on modal matrix
        om_n  = self.modes.xdata.angular_frequency
        om_n2 = om_n**2*(1+1j*self.damping_loss(om_n))-omega**2
        
        return om_n2
    
    def damping_loss(self,omega):
        """
        Damping loss of SEA systemns
        
        Arguments
        ---------
            omega : ndarray 
                angular frequency 
            
        """
        return np.array([self._damping_loss]*len(omega))
    
    
    @property
    def loads(self):
        """
        Property method for loads

        Returns
        -------
        laod
            load of FE model.

        """
        return self._loads

    @property
    def results(self):
        """
        Property method for results
        

        Returns
        -------
        Signal
            results of FEmoodel.

        """
        return self._results
    
    def modal_force(self,key):
        """
        Calculates forve vector in modal coordinates

        Parameters
        ----------
        key : int
            loadID.

        Returns
        -------
        Signal
            modal force vector.

        """
        if self.loads[key].dof.unique_type.typestr == 'force':
            F = self.loads[key]
        
        indexF = self.modes.dof.find(F.dof.ID,F.dof.dof)
        _F                = np.zeros(self.modes.Nsig, dtype = np.complex128)
        _R                = np.zeros((self.modes.Nxdata,F.Nxdata),dtype = np.complex128 )

        _M          = self.modes.ydata.conj().T
        
        #loop over frequency
        for ifreq in range(F.Nxdata):
            _F[indexF]  = F.ydata[:,ifreq]
            _R[:,ifreq] = _M.dot(_F) # works also for scalars
            
        resdof = self.modes.dof.unique_type*F.dof.unique_type
        return mC.Signal(F._xdata, _R, dof = dof.DOF(1+np.arange(self.modes.Nxdata),[0],resdof,repetition = True))
    
    def solve(self):
        """
        solver of FE model - not yet implemented

        Returns
        -------
        str
            DESCRIPTION.

        """
        return 'HUU'
    
    def modal_dof(self,doftype):
        """
        provides modal degress of freedom from mode set 
        
        Parameters
        ----------
        doftype : DOFtype
            type multiplier e.g. displacment.

        Returns
        -------
        DOF
            modal degrees of freedom.

        """
        
        modal_type = dof.DOFtype(typestr = 'mass', exponent = 0.5)*doftype
        return dof.DOF(1+np.arange(self.modes.Nxdata),[0],modal_type,repetition = True)
        
        
        
    def modal_vec_to_space(self,q_mod):
        """
        Transformation of model space to grid space

        Parameters
        ----------
        q_mod : mC.Signal
            Vector in modal coordinates.

        Returns
        -------
        mC.Signal
            Vector in mesh coordinates.

        """
 
        q = np.zeros((self.modes.Nsig,q_mod.Nxdata),dtype = np.complex128)
    
        for i_freq in range(q_mod.Nxdata):   
            q[:,i_freq] = self.modes.ydata @ q_mod[:,i_freq]
            
        q_type = self.modes.type*q_mod.type
        
        return mC.Signal(q_mod.xdata, q, dof.DOF(self.modes.ID,self.modes.dof,q_type))

    def modal_matrix_to_space(self,mat_mod):
        """
        Transforms modal matrices to grid space

        Parameters
        ----------
        q_mod : mC.Signal
            Vector in modal coordinates.

        Returns
        -------
        mC.Signal
            Vector in mesh coordinates.

        """
 
        q_ = np.zeros((self.modes.Nsig,self.modes.Nsig,mat_mod.shape[2]),dtype = np.complex128)
    
        for i_freq in range(mat_mod.shape[2]):   
            q_[:,:,i_freq] = self.modes.ydata @ mat_mod[:,:,i_freq] @ (self.modes.ydata.conj().T)
            
        return q_

    
    def rms_vec_from_modal(self,q_mod,ID = 0):
        """
        Calculates rms average over all nodes from modal result 

        Parameters
        ----------
        q_mod : Signal
            modal vevtor.
        ID : int, optional
            result ID of output. The default is 0.

        Returns
        -------
        Signal
            rms avrage of result.

        """
        
        q_ = np.zeros((1,q_mod.Nxdata),dtype = np.complex128)
    
        for i_freq in range(q_mod.Nxdata):   
            q_[0,i_freq] = np.sqrt(np.mean((self.modes.ydata @ q_mod.ydata[:,i_freq])**2))
            
        q_type = self.modes.dof.unique_type*q_mod.dof.unique_type
        
        if np.unique(self.modes.dof.dof).size == 1:
           dof_ =  np.unique(self.modes.dof.dof)
        else:
           dof_ = 0
        
        return mC.Signal(q_mod.xdata, q_, dof.DOF(ID,dof_,q_type))
    
    def rms_vec_from_modal_cpsq(self,q_mod,ID = 0,sqq_type = dof.DOFtype()):
        """
        Calculates rms average over all nodes from csd result 

        Parameters
        ----------
        q_mod : ndarray
            csd result.
        ID : int, optional
            ID for ouput Signal. The default is 0.
        sqq_type : DOFtype, optional
            type of result, eg. 'displacement'. The default is dof.DOFtype().

        Returns
        -------
        ndarray, DOF
            rms average of result, dofs

        """
        
        q_ = np.zeros((1,q_mod.shape[2]),dtype = np.complex128)
    
        for i_freq in range(q_mod.shape[2]):
            M_ = self.modes.ydata @ q_mod[:,:,i_freq] @ self.modes.ydata.conj().T
            q_[0,i_freq] = np.sqrt(np.mean(np.diag(M_)))
            
        q_type = self.modes.dof.unique_type*sqq_type
        
        if np.unique(self.modes.dof.dof).size == 1:
           dof_ =  np.unique(self.modes.dof.dof)
        else:
           dof_ = 0
        
#        return mC.Signal(q_mod.xdata, q_, dof.DOF(ID,dof_,q_type))
        return q_, dof.DOF(ID,dof_,q_type)
        
class TMmodel:
    
    """
    The TMmodel class deals with systems that are described best by transfer matrices
    
    The TM-method uses cascades of matrices to calculate the propagation through the sequentially 
    connected 1D subsystems with DOF pressure and velocity/volume velocity
    
    The RHS is considered as excitation port, the LHS as response port.  
    
    .. math:: 
        \\begin{Bmatrix} p_1 \\\\ v_1  \\end{Bmatrix} =
        \\begin{bmatrix} 
        T_{11} & T_{12} \\\\
        T_{21} & T_{22} 
        \\end{bmatrix} 
        \\begin{Bmatrix} p_2 \\\\ v_2 \\end{Bmatrix}
        
    Usally the coefficients :math:`T_{ij}(\\omega)` are functions of frequency 
    but for infinite layers they are also a function of wavenumber 
    :math:`T_{ij}(\\omega,k_x)`.
    
    """
    
    def __init__(self,layers,**kwargs):
        """
        Class constructor of TMmodel

        Parameters
        ----------
        layers : tuple or lost of AcousticLayer
            lay-up of different infinitly extende layers
        **kwargs : dict
            Arbitrary keyword agrguments

        Returns
        -------
        None.

        """
        
    
        #super().__init__(data,xdata,excdof,resdof,**kwargs)
        #self._loads = {}
        #self._results = {}
        #self._boundary_conditions = {}
        self._scalar_param = 0. # defauls kx = 0 
                    
        for kw in kwargs:
            if kw == 'loads':
                 self._loads.update = kwargs[kw]
            elif kw == 'results':
                 self._results.update = kwargs[kw]
            elif kw in {'BC','boundary_conditions'}:
                 self._boundary_conditions = kwargs[kw]
            else:
                 print('Unkown argument')
        
        if isinstance(layers,list):         
            self.layers = layers
        elif isinstance(layers,tuple):         
            self.layers = list(layers)
        else:
            raise(ValueError,'layers must be a list or tuple of layers')
    
    @property            
    def N(self):
        return len(self.layers)

    def __repr__(self):
        """
        Overloaded __repr__ method

        Returns
        -------
        _str : str
            string representation of TM_model.

        """
        
        #_str = super().__repr__()
        _str = 'TMmodel of layers:\n'
        for isys,lsys in enumerate(self.layers):
            _str += " {0:d}: {1} \n".format(isys,lsys)
            
        return _str
    
    def add_BC(self,BC):
        """
        Adds boundary condition

        Parameters
        ----------
        BC : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        self._boundary_conditions.update(BC)

            
    @property
    def loads(self):
        """
        Property methods for loads

        Returns
        -------
        load
            loads.

        """
        
        return self._loads

    def kx_DataAxis(self,kx):
        """
        Creates DataAxis from kx vector

        Parameters
        ----------
        kx : ndarray
            wavenumber in x-direction.

        Returns
        -------
        DataAxis
            wavenumber in x-direction.

        """
        
        return mC.DataAxis(kx,typestr = 'wavenumber')

    @property
    def boundary_conditions(self):
        """
        Property method for boundary conditions

        Returns
        -------
        TYPE
            boundary conditions.

        """
        return self._boundary_conditions

    @property
    def BC(self):
        """
        Property method for boundary conditions (BC) form

        Returns
        -------
        TYPE
            boundary conditions.

        """
        return self._boundary_conditions
    
    def V0(self,boundary_condition = 'equivalent_fluid'):
        """
        Multilayer state variable V_0 of Allards D0 matrix.
        
        The state variable is determined by the DOF left of the connection
        and the right side of the each layer.
        
        Layers of same nature are considered as one layer because they can be 
        calculated by simple transfermatrix multiplication.
        
        So, every change in layer nature require new DOFs.

        Returns
        -------
        v0 : mC.DOF
            State variable of multilayer.

        """
        ID_  = [0]*2 # we start with 0 for the entry fluid ID
        dof_ = [0,3]
        doftype_ = ('pressure','velocity')
 
        v0 = dof.DOF(ID_,dof_,doftype_)
        # first resdof is always 
        
        # Set current type from first layer  
        current_type = self.layers[0].type
        #ix_current = 0  
        # and set response dofs depending on the first layer nature
        if current_type == 'equivalent_fluid':
            v1 = iL.fluid_fluid_res_dof(0)
        elif current_type =='solid':
            v1 = iL.solid_fluid_res_dof(1) # for response, always
        elif current_type == 'poro_elastic':
            v1 = iL.porous_fluid_res_dof(1)
        
        for ix in range(1,self.N):
            # check if the next layer is of different type
            if current_type != self.layers[ix].type: #New DOFs are required
                if current_type == 'equivalent_fluid':
                    v0 = v0 + iL.fluid_exc_dof(ix*2)
                    if self.layers[ix].type == 'solid':
                        v1 = v1 + iL.solid_fluid_res_dof(ix*2+1) # fluid-solid takes right response ID
                    if self.layers[ix].type == 'poro_elastic':
                        v1 = v1 + iL.porous_fluid_res_dof(ix*2+1) # fluid-solid takes right response ID
                elif current_type == 'solid':
                    v0 = v0 + iL.solid_exc_dof(ix*2)
                    if self.layers[ix].type == 'equivalent_fluid':
                        v1 = v1 + iL.solid_fluid_res_dof(ix*2)
                    if self.layers[ix].type == 'poro_elastic':
                        v1 = v1 + iL.solid_porous_res_dof(ix*2+1) # solid-porous takes right response ID
                elif current_type == 'poro_elastic':
                    v0 = v0 + iL.porous_exc_dof(ix*2)
                    if self.layers[ix].type == 'equivalent_fluid':
                        v1 = v1 + iL.porous_fluid_res_dof(ix*2)
                    if self.layers[ix].type == 'solid':
                        v1 = v1 + iL.solid_porous_res_dof(ix*2) # porous-solid takes left response ID
                        
                #ix_current = ix      
                current_type = self.layers[ix].type
        
        # Deal with the last layer 
        if current_type == 'equivalent_fluid':
            v0 = v0 + iL.fluid_exc_dof(self.N*2)
            #v1 = v1 + iL.fluid_fluid_res_dof(self.N*2)
        elif current_type == 'solid':
            v0 = v0 + iL.solid_exc_dof(self.N*2)
            #v1 = v1 + iL.solid_fluid_res_dof(self.N*2)
        elif current_type == 'solid':
            v0 = v0 + iL.porous_exc_dof(self.N*2)

        # Add row depending on boundary condition
        if boundary_condition == 'fixed':
            if current_type == 'equivalent_fluid':
                v1 = v1 + dof.DOF([self.N*2+1],[3],['velocity'])
            elif current_type == 'solid':
                v1 = v1 + dof.DOF([self.N*2+1,self.N*2+1],[1,3],['velocity','velocity'])
            elif current_type == 'poro_elastic':
                v1 = v1 + dof.DOF([self.N*2+1,self.N*2+1],[1,3,2],['velocity','velocity','velocity'])
        elif boundary_condition == 'equivalent_fluid':
            # Deal with condition of eq. (11.83)
            v0 = v0 + iL.fluid_exc_dof(self.N*2+1)
            # set response dof regarding 
            if current_type == 'equivalent_fluid':
                v1 = v1 + iL.fluid_fluid_res_dof(self.N*2)
            elif current_type == 'solid':
                v1 = v1 + iL.solid_fluid_res_dof(self.N*2)
            elif current_type == 'poro_elastic':
                v1 = v1 + iL.porous_fluid_res_dof(self.N*2)
            # finally implement (11.84)
            v1 = v1 + dof.DOF([self.N*2+1],[0],['pressure'])

        return v0,v1
    
    def allard_matrix(self,omega,kx=0.,boundary_condition = 'equivalent_fluid', out_fluid = matC.Fluid(),reduced = False):
        """
        Calculate Allard matrix of all layers.

        Namely the matrix given by equation (11.79) and the terminiation given by (11.82) or (11.85)

        Parameters
        ----------
        omega : float
            angular frequency.
        kx : ndarray of float, optional
            wavenumber in x-direction. The default is 0..
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        TM_ : DynamicMatrix
            overall Allard Matrix D0 matrix.

        """
        # 
        #   #
        #   
        # Determine the xdata depending on the combination of kx and omega 
        xdata = iL.AcousticLayer.get_xdata(omega, kx)
        
        # Create empty dynamic matrix
        
        # The current DOF 1 is the DOF just left of the next interface
        current_left_ID = 0
        
        #SIF_in  = aR.HalfSpace(fluids[0])
        SIF_out = aR.HalfSpace(out_fluid)
        
        #Z_in  = SIF_in.radiation_impedance_wavenumber(omega,kx)
        Z_out = SIF_out.radiation_impedance_wavenumber(omega,kx)
        
        
        # Set current type from first layer  
        current_type = 'equivalent_fluid'
        layer_set_finished = False # no element creation as long as set is not finished
        layer_set_started  = False # no element creation as long as set is not finished
        
        V0,V1 = self.V0(boundary_condition=boundary_condition)
        D0 = mC.DynamicMatrix.zeros(xdata,V0,V1,0,dtype = np.complex)
        
        for ix in range(self.N):
            # check if the next layer is of different type
            if current_type != self.layers[ix].type or ix==0: #New DOFs are required
                if layer_set_started:
                    D0 += I12
                    D0 += J12
                else:
                    layer_set_started = True
                    
                    
                # Layer nature has changed
                if current_type == 'equivalent_fluid':
                    # fluid - fluid connection, only possible at first layer
                    if self.layers[ix].type == 'equivalent_fluid':
                        I12 = iL.I_fluid_fluid(xdata,2*ix)        # both M2
                        J12 = iL.J_fluid_fluid(xdata,2*ix,2*ix+1) # e.g. M2 and M3 
                                                                            
                    # fluid - solid connection
                    if self.layers[ix].type == 'solid':
                        I12 = iL.J_solid_fluid(xdata,2*ix+1,2*ix)
                        J12 = iL.I_solid_fluid(xdata,2*ix+1,2*ix+1) # e.g. M2 and M3 
                        
                elif current_type == 'solid':
                    # solid - fluid connection
                    if self.layers[ix].type == 'equivalent_fluid':
                        I12 = iL.I_solid_fluid(xdata,2*ix,2*ix)
                        J12 = iL.J_solid_fluid(xdata,2*ix,2*ix+1) # e.g. M2 and M3 
                elif current_type == 'poro_elastic':
                    pass
                        
                # set current type accordinmg to the new layer
                current_type = self.layers[ix].type
                # and multiply with tansfermatrix to reach the right side of the layer
                T0 = self.layers[ix].transfer_impedance(omega,kx,ID=[2*ix+1,2*ix+2])
                J12 = J12.dot(T0)                
                
            else: # go to next layer and use transfermatrix with out recreating I12 and J12
                # In a first step the transfermatrices are multiplied without I because
                # there is currently no porosity involved
                J12 = J12.dot(self.layers[ix].transfer_impedance(omega,kx,ID = [2*ix,2*ix+2])) # 1st ID is set to left 
          
            # consider last layer
            if ix == self.N-1:
                D0 += I12
                D0 += J12
            
        # Deal with boundary condition
        if boundary_condition == 'fixed':
            if current_type == 'equivalent_fluid':
                D0 += iL.Y_fluid_fixed(xdata,2*self.N)
            elif current_type == 'solid':
                D0 += iL.Y_solid_fixed(xdata,2*self.N)
        elif boundary_condition == 'equivalent_fluid': # (11.83)
            if current_type == 'equivalent_fluid': 
                D0 += iL.I_fluid_fluid(xdata,2*self.N)        # both M2
                D0 += iL.J_fluid_fluid(xdata,2*self.N,2*self.N+1) # e.g. M2 and M3 
            elif current_type == 'solid': 
                D0 += iL.I_solid_fluid(xdata,2*self.N,2*self.N)
                D0 += iL.J_solid_fluid(xdata,2*self.N,2*self.N+1) # e.g. M2 and M3 
            # Use SIF radiation impedanz to get the radiation condition
            data_ = np.zeros((1,2,len(xdata) ),dtype=np.complex128)
            data_[0,0,:] = -1
            data_[0,1,:] = Z_out
                            
            D0 += mC.DynamicMatrix(data_, xdata, iL.fluid_exc_dof(2*self.N+1), dof.DOF([self.N*2+1],[0],['pressure']))
        
            # Consider First Layer
            
        if reduced:
            # set the pressure in V to 1 and create an excitation vector from the first matrix column
            D1 = D0[:,1:,:]
            # Create force from first column
            F = -D0.signal(slice(0,None),0)
            return D1,F
        else:
            return D0

    #def local_allard_matrix(self,omega,kx=0.,**kwargs):
        
    def connect(self,omega,kx=0.,**kwargs):
        """
        Calculate TM of all layers by simple transfermatrix multiplication.
        
        Multiplication makes sense only for layers of same nature. If porous 
        layers are involved (and they are not modelled as equivalent fluid) the
        multiplication involves a connection matrix [I] that consideres the porosity.
    
        Parameters
        ----------
        omega : float
            angular frequency.
        kx : ndarray of float, optional
            wavenumber in x-direction. The default is 0..
        **kwargs : TYPE
            DESCRIPTION.
    
        Returns
        -------
        TM_ : DynamicMatrix
            overall transfer matrix.
    
        """
        # collects the TMM        
        
        if self.layers[0].isequivalentfluid:
            TM_ = self.layers[0].transfer_impedance(omega,kx,ID = [1,2])
        else:
            raise(ValueError,'Only equivalent fluid layer can be connected via transfermatrix muliplication')
            
        for ix,TM in enumerate(self.layers[1:]):
            # check is each layer is of type equivalent fluid
            if TM.isequivalentfluid:
                TM_ = TM_.dot(TM.transfer_impedance(omega,kx,ID = [ix+2,ix+3]))
            else:
                raise(ValueError,'Only equivalent fluid layer can be connected via transfermatrix muliplication')
    
                
        return TM_
        
            
    def impedance(self,omega,kx=0.,ID=1, boundary_condition = 'fixed',signal = True,out_fluid = matC.Fluid()):
        """
        impedance provides the surface impedance of a TMmodel with specified end condition
        
        Parameters
        ----------
        ID : int
            ID for surface condition (can be an inner ID)
        boundary_condition : str
            at end of layers
        signal : bool, optional
            Switch for Signal (True) or array (False) output. The default is True.

                
                
        Returns
        -------
        Signal 
            impedance at port one   
        """
        
        TM_ = self.connect(omega,kx)
        
        # get load
        if isinstance(boundary_condition, str ) and boundary_condition == 'fixed':
            _z = TM_.data[0,0,:]/TM_.data[1,0,:]
            _doftype = TM_.resdof.type[0]/TM_.resdof.type[1]
        elif isinstance(boundary_condition, str ) and boundary_condition == 'free':
            _z = TM_.data[1,0,:]/TM_._data[1,1,:]
            _doftype = TM_.resdof.type[0]/TM_.resdof.type[1]
        elif isinstance(boundary_condition, str ) and boundary_condition == 'equivalent_fluid':
            SIF_out = aR.HalfSpace(out_fluid)
            z_out = SIF_out.radiation_impedance_wavenumber(omega,kx)
            _z = (TM_.data[0,0,:]*z_out+TM_.data[0,1,:])/(TM_.data[1,0,:]*z_out+TM_.data[1,1,:])
            _doftype = TM_.resdof.type[0]/TM_.resdof.type[1]
            
            
            
        # @todo pick our any ID even in the middle, default 1
        # @todo 
        # @todo create impedance or admittance BC
        
        # xdata predifined by non-scalar element kx or omega
        if signal:
            xdata_ = TM_.xdata # turn single values into arrays
            return mC.Signal(xdata_,_z,dof.DOF(ID,3,_doftype))
        else:
            return _z

    def impedance_allard(self,omega,kx=0.,ID=1, boundary_condition = 'fixed', signal = True,out_fluid = matC.Fluid()):
        """
        impedance provides the surface impedance of a TMmodel with specified end condition
        
        Parameters
        ----------
        ID : int
            ID for surface condition (can be an inner ID)
        boundary_condition : str
            at end of layers
        signal : bool, optional
            Switch for Signal (True) or array (False) output. The default is True.

                
                
        Returns
        -------
        Signal 
            impedance at port one   
        """
        
        
        # Create reduces Allard matrix
        D1,F = self.allard_matrix(omega, kx = kx, reduced=True, boundary_condition = boundary_condition, out_fluid = out_fluid)
        V    = D1.solve(F) 
                    
        Z_s = 1./V.ydata[0,:]
        
        # xdata predifined by non-scalar element kx or omega
        if signal:
            _doftype = dof.DOFtype(typestr = 'pressure')/F.dof.type[0]
            xdata_ = D1.xdata # turn single values into arrays
            return mC.Signal(xdata_,Z_s,_doftype)
        else:
            # ix 0 is the index of the velcoity at the surface
            return Z_s
    
    def transmission_allard(self,omega,kx=0.,ID=1, fluids = (matC.Fluid(),matC.Fluid()), signal = True,plot = 0):
        """
        Calculate the transmission of TMM according to Allards method.
        
        Parameters
        ----------
        ID : int
            ID for surface condition (can be an inner ID)
        boundary_condition : str
            at end of layers
        signal : bool, optional
            Switch for Signal (True) or array (False) output. The default is True.
    
                
                
        Returns
        -------
        Signal 
            impedance at port one   
        """
        
        
        # Create reduces Allard matrix
        D1,F = self.allard_matrix(omega, kx = kx, reduced=True, 
                                  boundary_condition = 'equivalent_fluid',out_fluid=fluids[1])
        V    = D1.solve(F) 
                    
        Z_s = 1./V.ydata[0,:] # = 1/V_3(A)
        SIF_in = aR.HalfSpace(fluids[0])
        Zin = np.real(SIF_in.radiation_impedance_wavenumber(omega,kx))
              
        R = (Z_s-Zin)/(Z_s+Zin)
        #print(R)
        T = (1+R)*V.ydata[-2,:] # = (1+R)*P(B)
        # Validation plot if requestet
        if plot > 0:
            V[-2].plot(plot,marker='.',res = 'real')
            V[-2].plot(plot,marker='>',res = 'imag')
        
        tau = np.abs(T)**2
        
        # xdata predifined by non-scalar element kx or omega
        if signal:
            _doftype = dof.DOFtype(typestr = 'transmission')
            xdata_ = D1.xdata # turn single values into arrays
            return mC.Signal(xdata_,tau,_doftype)
        else:
            # ix 0 is the index of the velcoity at the surface
            return tau #,1-np.abs(R**2) 
    

    def stiffness_matrix(self,omega,distances,ks,dA,Nstep = 10):
        """
        Calculate stiffness matrix of TMmodel

        Parameters
        ----------
        omega : float
            angular frequency.
        distances : ndarray
            distances for discrete mesh.
        ks : fload
            maximum wavenumber integration limit
        dA : float
            element area
        Nstep : int, optional
            Number of intervals for numeric wavenumber integration. The default is 10.

        Returns
        -------
        D11 : ndarray
            1-1 coeffcient of stifness matrix.
        D12 : ndarray
            1-2 coeffcient of stifness matrix equal to D21.
        D22 : ndarray
            2-2 coeffcient of stifness matrix.

        """
        
        # \todo Calulate wavenumber sampling....
        J_1st_zero = 2.40482556
        dks = J_1st_zero/Nstep
        
        kx = np.linspace(0,ks,int(np.ceil(ks*np.max(distances)/dks))+1)
        fak = dA/2/np.pi
        
        _D = self.stiffness_matrix_wavenumber(omega,kx)
        
        
        # create besselmatrix with rows for distance and coloums for kx 
        J0 = scl.j0(distances.reshape(-1,1) @ kx.reshape(1,-1) )
        D11 = fak*integrate.trapz(_D[0,0,:].data.flatten()*J0*kx,kx)
        D12 = fak*integrate.trapz(_D[0,1,:].data.flatten()*J0*kx,kx)
        D22 = fak*integrate.trapz(_D[1,1,:].data.flatten()*J0*kx,kx) # ,axis = 1
        
        return (D11,D12,D22)
        
    
    def stiffness_matrix_mesh(self,omega,mesh):
        """
        Calculate stiffness matrix of TMmodel for mesh   

        Parameters
        ----------
        omega : float
            angular frequency.
        mesh : mesh
            grid points on surface.

        Returns
        -------
        LinearMatrix
            Stiffness matrix of TMmodel for mesh.

        """
        
        distances, index = mesh.distance()
        
        Ndist  = len(distances)
        Nx     = len(omega)
        Nmesh  = mesh.Nmesh
        dA     = mesh.dA
        ks     = mesh.ks
        
        D11 = np.zeros(((Nmesh*(Nmesh+1))//2,Nx),dtype=np.complex128)
        D12 = np.zeros(((Nmesh*(Nmesh+1))//2,Nx),dtype=np.complex128)
        D22 = np.zeros(((Nmesh*(Nmesh+1))//2,Nx),dtype=np.complex128)
        
        for ix in range(Nx):
            #for id in range(Ndist):
            _A,_B,_C = self.stiffness_matrix(omega[ix],distances,ks,dA)
            D11[:,ix], D12[:,ix],D22[:,ix] =  (_A[index],_B[index],_C[index])
                    
        return (mC.LinearMatrix(D11,sym = 1,shape=(Nmesh,Nmesh,Nx) ),
                mC.LinearMatrix(D12,sym = 1,shape=(Nmesh,Nmesh,Nx) ),
                mC.LinearMatrix(D22,sym = 1,shape=(Nmesh,Nmesh,Nx) ))
    
    def stiffness_matrix_wavenumber(self,omega,kx=0.):
        """
        provides the stiffness matrix of TMM
        
        Parameters
        ----------
        omega : float
            angular frequency.
        kx : ndarray or float, optional
            wavenumber in x-direction. The default is 0..

        Returns
        -------
        LinearMatrix
            stiffness matrix in wavenumber space.

        """
        
        
        TM_ = self.connect(omega,kx)
        
        # get coefficients
        T12 = TM_.data[0,1,:]
        T21 = TM_.data[1,0,:]
        T11 = TM_.data[0,0,:]
        T22 = TM_.data[1,1,:]
        
        det = T11*T22-T12*T21
        
        Dx = np.ones(np.shape(TM_.data),dtype = np.complex128)

        fak = 1j*omega/T21
        Dx[0,0,:] = fak*T11
        Dx[0,1,:] = -fak*det # shout be 1
        Dx[1,1,:] = fak*T22
        Dx[1,0,:] = -fak
        
        # Update ID in DOF attribute
        IDs = [TM_.resdof.ID[0],TM_.resdof.ID[1]]
        
        # set DOF according to ID and natural DOF of the layer
        resdof = dof.DOF(TM_.resdof.ID, 3 , dof.DOFtype(typestr='pressure'),repetition=True)
        excdof = dof.DOF(IDs, 3,  dof.DOFtype(typestr='displacement'),repetition=True)
        
        return mC.DynamicMatrix(Dx,TM_.xdata,excdof,resdof)
    
    def insertion_loss(self,omega,kx,fluid = matC.Fluid(),Signal = False):
        """
        calculates the insertion loss for plate radiation 

        Parameters
        ----------
        omega : float
            angular frequency.
        kx : ndarray of float
            wavenumber in x-direction.
        fluid : matC.Fluid
            fluid material that is connected to the outer trim surface
        Signal : bool, optional
            switch if ndarray or singal output is desired.

        Returns
        -------
        insertion loss

        """
        
        Dx = self.stiffness_matrix_wavenumber(omega,kx)
        
        SIF_out = aR.HalfSpace(fluid)
        D_dir   = SIF_out.radiation_stiffness_wavenumber(omega,kx)
        
        IC = np.abs(Dx.data[1,0,:]/(Dx.data[1,1,:]+D_dir))**2

        if Signal:
            return mC.Signal(Dx.xdata,IC,dof.DOF(0,0,'generic'))
        else:
            return IC
        
    def insertion_loss_diffuse(self,omega,plate,theta_max = 78/180*np.pi,theta_step = np.pi/180,fluid = matC.Fluid()):
        """
        calculates the diffurse field insertion loss for plate radiation 
        
        Parameters
        ----------
        omega : float
            angular frequency.
        plate : 2Dstructure
            reference plate with applied trim.
        theta_max : float, optional
            maximum integration angle. The default is 78/180*np.pi.
        theta_step : float, optional
            integral sampling step size. The default is np.pi/180.
        fluid : matC.Fluid
            fluid material that is connected to the outer trim surface

        Returns
        -------
        TYPE
            insertion loss.

        """

        omega_ = np.array(omega).flatten() # ensure np.array
        #IC_diffuse = np.zeros((omega_.size,)) 
        #denom = np.sin(theta_max)**2
        
        # call method with trim = False option for baseline
        tau_plate = plate.resonant_TMM(trim=False).transmission_diffuse(omega_,fluids=(fluid,fluid),signal = False)
        tau_trim  = plate.resonant_TMM().transmission_diffuse(omega_,fluids=(fluid,fluid),signal = False)
        
        return tau_trim/tau_plate
               
    
    def absorption(self,omega,kx,in_fluid = matC.Fluid(),ID=1, boundary_condition = 'fixed', allard = False, signal = True,out_fluid = matC.Fluid()):
        """
        Calculates the surface absorption with specified end condition and input fluid

        Parameters
        ----------
        omega : float
            angular frequency
        kx : ndarray of float
            surface wavenumber.
        in_fluid : fluid, optional
            irradiating fluid. The default is matC.Fluid().
        ID : int, optional
            ID for surface condition (can be an inner ID). The default is 1.
        boundary_condition : str or Signal, optional
            boundary_condition. The default is 'fixed'.
        signal : bool, optional
            Switch for Signal (True) or array (False) output. The default is True.

        Returns
        -------
        Signal
            absorption coefficient    
        """
        
        
        SIF = aR.HalfSpace(in_fluid)
        
        Zin = SIF.radiation_impedance_wavenumber(omega,kx)

        if allard:
            Zsurf  = self.impedance_allard(omega, kx, ID, boundary_condition, out_fluid=out_fluid)
        else:
            Zsurf  = self.impedance(omega, kx, ID, boundary_condition)
        
        xdata_ = Zsurf.xdata # turn single values into arrays
        
        refl_ = (Zsurf.ydata-Zin)/(Zsurf.ydata+Zin)
        abs_  = 1-np.abs(refl_)**2
        
        if signal:
            return mC.Signal(xdata_,abs_,dof.DOF(ID,0,dof.DOFtype(typestr='general')))
        else:
            return abs_.flatten()
   
    
    def absorption_diffuse(self,omega,theta_max = 78/180*np.pi,theta_step = np.pi/180,\
                           in_fluid = matC.Fluid(),ID=1,boundary_condition= 'fixed', allard = False, \
                           signal = True,out_fluid = matC.Fluid()):
        """
        diffuse surface absorption with specified end condition and input fluid

        Parameters
        ----------
        omega : ndarray or float
            angular frequency
        theta_max : float, optional
            Maximum angle for diffuse field integeration. The default is 78/180*np.pi.
        theta_step : float, optional
            Angle step for diffuse field integeration. The default is np.pi/180.
        in_fluid : fluid, optional
            Fluid of sound field. The default is matC.Fluid().
        ID : int, optional
            ID of considered layer. The default is 1.
        boundary_condition : str, optional
            End condtion. The default is 'fixed'.
        signal : bool
            Switch for Signal output,  The default is True.
        allard : bool
            Switch for calculation method, The default is False

        Returns
        -------
        Signal, ndarray
            Absorption coefficient   
    
        """
        
        # \todo convention for first layer is structure layer must be implemented
        
        omega_ = np.array(omega).flatten()
        abs_diffuse = np.zeros((omega_.size,))
        
        theta_ = np.linspace(0,theta_max,int(np.floor((theta_max+theta_step)/theta_step)))
        
        for ix,om in enumerate(omega_):

            k  = np.real(in_fluid.wavenumber(om))
            kx_    = np.sin(theta_)*k
            abs_kx = self.absorption(om,kx_,in_fluid,ID,boundary_condition,allard = allard,out_fluid = matC.Fluid()).ydata
            #remove nans
            abs_kx[np.isnan(abs_kx)] = 0.
            denom = np.sin(theta_max)**2
            abs_diffuse[ix] = 2*integrate.simps(abs_kx*np.sin(theta_)*np.cos(theta_), theta_)/denom

        if uf.isscalar(omega):
            return abs_diffuse[0]
        elif signal:
            xdata = mC.DataAxis(omega_, typestr = 'angular frequency')
            return mC.Signal(xdata,abs_diffuse,dof.DOF(ID, 0, dof.DOFtype(typestr = 'general')))
        else:
            return abs_diffuse
        
    def transmission_coefficient(self,omega,kx,fluids = (matC.Fluid(),matC.Fluid()),ID=0,signal = True):
        """
        provides the acoustic transmission coefficient of TMmodel objects
        

        Parameters
        ----------
        omega : TYPE
            angular frequency.
        kx : TYPE
            surface wavenumber.
        fluids : tuple of fluid, optional
            2x1 vector of fluids. The default is (matC.Fluid(),matC.Fluid()).
        ID : int, optional
            node ID for output. The default is 0.

        Returns
        -------
        Signal
            transmission coefficient.

        """
        
        if uf.isscalar(omega):
            xdata  = self.kx_DataAxis(kx)
        elif uf.isscalar(kx):
            xdata  = mC.DataAxis(np.array(omega).flatten(),type = dof.DOFtype(typestr='angular frequency'))
        elif ~(uf.isscalar(omega) or uf.isscalar(kx)):
            if len(omega) == len(kx):
                xdata = mC.DataAxis(np.array(omega).flatten(),type = dof.DOFtype(typestr='angular frequency'))
            else:
                raise(ValueError,'kx and omega must be one scalar or len(kx) == len(omega)')
            
        
        
        SIF_in  = aR.HalfSpace(fluids[0])
        SIF_out = aR.HalfSpace(fluids[1])
        
        Z_in  = SIF_in.radiation_impedance_wavenumber(omega,kx)
        Z_out = SIF_out.radiation_impedance_wavenumber(omega,kx)

                
        kx = np.array(kx).flatten() # turn single values into arrays
        TM_ = self.connect(omega,kx)
        T11 = TM_.data[0,0,:]
        T12 = TM_.data[0,1,:]
        T21 = TM_.data[1,0,:]
        T22 = TM_.data[1,1,:]
        
        tau_  = 4*np.real(1/Z_out)/np.real(1/Z_in)/np.abs(T11+T12/Z_out+T21*Z_in+T22*Z_in/Z_out)**2
        
        if signal:
            return mC.Signal(xdata,tau_,dof.DOF(ID,0,dof.DOFtype(typestr='transmission')))        
        else:
            return tau_
        
    def transmission_diffuse(self,omega,theta_max = 78/180*np.pi,theta_step = np.pi/180,\
                                  fluids = (matC.Fluid(),matC.Fluid()),ID = 0,signal = True,
                                  allard = False):
        """
        

        Parameters
        ----------
        omega : float
            angular frequency.
        plate : 2Dstructure
            reference plate with applied trim.
        theta_max : float, optional
            maximum integration angle. The default is 78/180*np.pi.
        theta_step : float, optional
            integral sampling step size. The default is np.pi/180.
        fluids : tuple of fluids, optional
            vector of coupled fluids. The default is (matC.Fluid(),matC.Fluid()).
        ID : int, optional
            node ID for output Signal. The default is 0.
        signal : bool, optional
            Switch for Signal (True) or array (False) output. The default is True.

        Returns
        -------
        ndarray of Signal
            diffuse field transmission coefficient.

        """
                
        omega_ = np.array(omega).flatten()
        tau_diffuse = np.zeros((omega_.size,))
        denom = np.sin(theta_max)**2

        theta_ = np.linspace(0,theta_max,int(np.floor((theta_max+theta_step)/theta_step)))
               
        for ix,om in enumerate(omega_):
            k  = np.real(fluids[0].wavenumber(om))
            kx_    = np.sin(theta_)*k
            if allard:
                tau_kx = self.transmission_allard(om,kx_,fluids).ydata
            else:
                tau_kx = self.transmission_coefficient(om,kx_,fluids).ydata
            #remove nans
            tau_kx[np.isnan(tau_kx)] = 0.
            tau_diffuse[ix] = 2*integrate.trapz(tau_kx*np.sin(theta_)*np.cos(theta_), theta_)/denom

        if uf.isscalar(omega):
            return tau_diffuse[0]
        else:
            if signal:
                xdata = mC.DataAxis(omega_, typestr = 'angular frequency')
                return mC.Signal(xdata,tau_diffuse,dof.DOF(ID, 0, dof.DOFtype(typestr = 'transmission')))
            else:
                return tau_diffuse
      
class VAmodel(mC.DynamicMatrix):
    
    """
    The VAmodel is an extension of the DynamicMatrix class with the 
    additional
    
    The VAmodel (vibroacoustic model) class is used for models that can be
    represented by one dynamic matrix, for example the stiffness matrix
    
    .. math:: 
       \\left[\\bm D \\right] \\lbrace {\\bm u} \\rbrace = \\lbrace {\\bm F} \\rbrace
    
    Different loads can by applied. The DOF handling is done by DynamicMatrix class.
    
    Attributes
    ----------
    results : dict of Signals
        results due to loads
    loads : dict of loads
        loads
    
    
    """
    
    def __init__(self,data,xdata,excdof,resdof,sym=1,**kwargs):
        """
        Constructor of VAmodel class 

        Parameters
        ----------
        data : ndarray,LinearMatrix or None
            matrix over frequency. In case of None - empty model with zeros is created
        xdata : DataAxis
            frequency values.
        excdof : dof
            excitation degrees of freedom.
        resdof : dof
            response degrees of freedom.
        sym : int
            Symmetry identifier.
       **kwargs : dict
            Arbitrary keyword arguments passed to DynamicMatrix constructor.

        Returns
        -------
        None.

        """
        
        # create empty model with empty matrix
        if data is None:
            Ndepth = len(xdata)
            Nrow   = resdof.size
            Ncol   = excdof.size
        
            shape  = (Nrow,Ncol,Ndepth)
            data = mC.LinearMatrix.zeros(sym, shape, **kwargs)
            
        super().__init__(data,xdata,excdof,resdof,sym=sym,**kwargs)
            
        self._loads = {}
        self._results = {}
                    
        for kw in kwargs:
            if kw == 'loads':
                 self._loads.update = kwargs[kw]
            elif kw == 'results':
                 self._results.update = kwargs[kw]
            else:
                 print('Unkown argument')

    def __repr__(self):
        _str = super().__repr__()
        for lkey in self._loads.keys():
            _str += "Load with ID=" + str(lkey) + " " +str(self._loads[lkey]) + "\n"
            
        for rkey in self._results.keys():
            _str += "Results with ID=" + str(rkey) + " " + str(self._results[rkey]) + "\n"
            
        return _str
    
    def add_load(self,load):
        """
        Adds load to VAmodel

        Parameters
        ----------
        load : dict of Load
            excitation load.

        Returns
        -------
        None.

        """
   
        self._loads.update(load)

    def add_result(self,result):
        """
        Adds result to VAmodel

        Method is intended fot internal use, e.g. storing results from solve

        Parameters
        ----------
        load : dict of load
            excitation load.

        Returns
        -------
        None.

        """   
        self._results.update(result)
        
    def plot(self,nfig=1,loadID=1,**kwargs):
        """
        Plot results

        Parameters
        ----------
        nfig : int, optional
            figure number. The default is 1.
        loadID : int, optional
            key for load/result to plot. The default is 1.
        **kwargs : 
            Arbitrary keyword arguments passed to Signal.plot() .

        Returns
        -------
        None.

        """
        
        
        _res = self._results[loadID]
        _res.plot(nfig,**kwargs)
    
    @property
    def loads(self):
        """
        property method for loads

        Returns
        -------
        dict of loads
            loads.

        """
        
        return self._loads

    @property
    def results(self):
        """
        property method for results

        Returns
        -------
        dict of results
            results.

        """
        return self._results
    
    def solve(self,**kwargs):
        """
        Solve linear system of equations

        Parameters
        ----------
        **kwargs : dict
            Arbitrary keyword arguments.
        loadresponse : True when the response shall overwrite the load     
        

        Returns
        -------
        None.

        """
        
        
        lresp = False
        
        for kw in kwargs:
            if kw=='loadresponse':
                lresp = True
                
        
        # loop over loads
        for iload in self._loads.keys():
            print("Solving Load Case ID = {:d}".format(iload))

            _load = self._loads[iload]

            # set up load vector according to dofs
            _F     = np.zeros(self.Nrow, dtype = complex)
            if lresp:
                _LR     = np.zeros((self.Nrow,self.Ndepth),dtype = complex )
                
            comdof,indexF,ix2 = self._resdof.intersect(_load.dof) #,ignore='doftype')
            _R     = np.zeros((self.Nrow,self.Ndepth),dtype = complex )
            
                       
            #loop over frequency
            for ifreq in range(len(self._xdata)):
                
                _M = self.Dindex(ifreq)
                _F[indexF] = _load.ydata[:,ifreq]
                
                _R[:,ifreq] = linalg.solve(_M,_F)
                
                if lresp:
                    # Overwrite load with recalculated response
                    _LR[:,ifreq] = _M.dot(_R[:,ifreq])
                
            self.add_result({iload:mC.Signal(self._xdata, _R, dof = self._excdof )})
            if lresp:
                self.add_load({iload:mC.Signal(self._xdata, _LR, dof = self._resdof )})
            
            
    def power(self,LID,NID,boundary = None):
        """
        Calculate powerflow into network

        Parameters
        ----------
        LID : int
            load ID.
        NID : int
            node ID.
        boundary : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        """
        
        
        # get load
        exc = self._loads[LID]
        res = self._results[LID]

        _P = res.pressure(NID)
        
        # check if excitation is available
        if boundary == None:
            _Q = exc.volume_flow(NID)
        else:
            _Q = boundary.dot(_P)
            
        
        #print("q=",q)
        #print("p=",p)
        
#        _Q.plot(100,res = 'real')
#        _Q.plot(101,res = 'imag')
#        _P.plot(102,res = 'real')
#        _P.plot(103,res = 'imag')
        
        
        return (0.5*_P.conj()*_Q).real()

    
class HybridModel:
    """
    The HybridModel class provides the framework for the random part of a vibroacoustic model.
    
    
    
    Attributes
    ----------
    systems : dict
        SEA systems of the SEA model
    FEsystems : dict
        FEM systems of the SEA model
    junctions : dict
        junction defining the coupling of systems
    hybrid_junctions : dict
        junction defining the hybrid (via FE) coupling of systems
    SEAmatrix = DynamicMatrix
        SEAmatrix of CLF and DLF of the model
    energy : Signal
        energy of subsystems
    result : dict
        physical unit of subsystens, e.g. velocity
    hybrid_result : dict
        physical unit of FE subsystens, e.g. velocity
    loads : dict
        loads of the model
    """
    
    def __init__(self,systems=tuple(),FEsystems=tuple(),sifs=tuple(),
                      xdata=mC.DataAxis.octave_band(), # set default spectrum
                      **kwargs):
        """
        Class constructor of HybridModel

        Parameters
        ----------
        systems : tuple, optional
            List of SEA systems. The default is tuple().
        FEsystems : tuple, optional
            List of FE systems. The default is tuple().
        xdata : DataAxis
            Global frequency of the model. The default is mC.DataAxis.octave_band().
        loads :  load
            loads of the model.

        Returns
        -------
        None.

        """
    
        self.systems = dict()
        self.FEsystems = dict()
        self.sifs    = dict()
        for i_sys,syst in enumerate(systems):
            ID = syst.wave_DOF.ID[0]
            self.systems.update( {ID : systems[i_sys]})
        
        for i_sys,syst in enumerate(FEsystems):
            ID = syst.ID
            self.FEsystems.update( {ID : systems[i_sys]})
            
        self.xdata   = xdata
        self._wave_DOF = self.get_SEA_model_DOFs()
        self.junctions = dict()
        self.hybrid_junctions = dict()
        self.loads     = dict()
        self.SEAmatrix = False
        self.energy  = []
        self.result  = []
        self.hybrid_result = []
        
        self._current = True # switch is all attributes are currect
                           
        for kw in kwargs:
            if kw == 'loads':
                 self._loads.update = kwargs[kw]
            else:
                 print('Unkown argument')

    def __repr__(self):
        """
        overloaded repr method

        Returns
        -------
        _str : str
            Text summary of the model.

        """
        
        _str = 'HybridModel with subsystems:\n'
        
        for lsys in self.systems:
            _str += str(self.systems[lsys]) + "\n"
            
        for load_key in self.loads.keys():
            _str += 'Load with name:'+str(load_key) + "\n"
            
                        
        return _str
    
    def add_system(self,system):
        """
        Adds systems to the SEA model

        Parameters
        ----------
        system : SEAsytem or FEsystem
            system to be added.

        Returns
        -------
        None.

        """
        
        # use dictionary to add system
        ID = system.wave_DOF.ID[0]
        self.systems.update({ID:system})
        self._current = False

    def add_SIF(self,sifdict):
        """
        Adds semi infinite fluid (SIF) to the SEA model
        

        Parameters
        ----------
        sifdict : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        # use dictionary to add system
        self.sifs.update(sifdict)
 
    def add_junction(self,junctions):
        """
        Adds junction to the SEA model
        
        Parameters
        ----------
        junctions : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        # use dictionary to add system
        self.junctions.update(junctions)
        self._current = False

    def add_hybrid_junction(self,junctions):
        """
        Adds junction to the SEA model
        
        Parameters
        ----------
        junctions : hybrid_junction
            hybrid junction of the SEA model.

        Returns
        -------
        None.

        """
        
        # use dictionary to add system
        self.hybrid_junctions.update(junctions)
        self._current = False

    def add_load(self,key,load):
        """
        Adds load to the HybridModel

        Parameters
        ----------
        key : int
            load key.
        load : load
            load of the HybridModel.

        Returns
        -------
        None.

        """
        
        # use dictionary to add system
        self.loads.update( { key : load })
        
    def update(self):
        """
        Updates all relevant attributes of HybridModel
        """
        
        self._wave_DOF = self.get_SEA_model_DOFs()
        
    @property
    def N_freq(self):
        """
        Determines number of frequency lines

        Returns
        -------
        uint16
            number of frequency lines

        """
        return len(self.xdata) 

    @property
    def wave_DOF(self):
        """
        property method that gives the wave_DOFs of the SEA model

        Returns
        -------
        dof
            wave DOF of SEAsystems.

        """
        if not(self._current):
            self.update
        return self._wave_DOF
       
    @property
    def N_sys(self):
        """
        Determines number of physical systems (not wavefields)

        Returns
        -------
        uint16
            number of physical systems

        """
        return self.systems.count()
    
    @property
    def N_sea(self):
        """
        Determines number of SEA systems (of wavefields)
        
        Structural SEA systems as plates or beams may have several wave fields that are considererd
        as different reverberant fields and thus as different subsystems. 
        The number of reverberant fields considered considered as independent is determined here.

        Returns
        -------
        uint16
            number of reverberant fields/SEA systems

        """
    
        N_rev = 0
        
        for syst in self.systems.values():
            N_rev += syst.N_wave_fields
            
        return N_rev
            
    def index(self,ID):
        """
        Determines index of ID in system list of HybridModel

        Parameters
        ----------
        ID : integer
            System ID

        Returns
        -------
        integer
            index of ID in system list (-1) if not found

        """
        
        ix = -1
        
        for ix in range(self.N):
            if ID == self.systems[ix].ID:
                return ix
        return -1
        
    
    def create_SEA_matrix(self,sym = 1,force = None):
        """
        Creates the SEA matrix from subystems and junction
        
        This method loops over subsystems and junctions to set up the SEA matrix, 
        including the hybrid junctions. In case of hybrid junction excitation 
        by external laods is included.
        (not from the reverberant SEA fields) 

        Parameters
        ----------
        sym : int, optional
            Symmmetry of SEA matrix. The default is 1 (symmetric).
        force : force, optional
            Force loads on FEsystems. The default is None.

        Returns
        -------
        int
            DESCRIPTION.

        """
        
        N_sea  = self.N_sea 
        N_freq = len(self.xdata)
        omega = self.xdata.angular_frequency
        
        # preset
        data = np.zeros((N_sea*(N_sea+1)//2,N_freq))
        # create symmetric LinearMatrix with zeros
        self.SEAmatrix = mC.DynamicMatrix(data, mC.DataAxis(omega), self.exc_DOF, self.res_DOF, sym = 1 , shape = (N_sea,N_sea,N_freq))
        
        if self._current == False:
            self.update()
        
        #modal density over all SEA DOFs 
        m_dens =  self.get_modal_density()
        
        
    
        # Deal with damping \eta_nn entries in diagonal 
        i_sea = 0 # set counter
       
        for syst in self.systems.values():
            i_sys = self._wave_DOF.index(syst.wave_DOF)
            
            for i_local,i_global in enumerate(i_sys):
                # pre-fill diagonal with n_i*eta_ii
                self.SEAmatrix[i_sea,i_sea,:] = m_dens[i_sea,:]*syst.damping_loss(omega,i_local)
                i_sea += 1
                
        # Deal with damping linked to sinks (SIF)        
        for sif in self.sifs.values():
            # loop over connected subsystem to sifs
            for syst in sif.systems:
                # pre-fill diagonal with n_i*eta_ii
                if syst.iscavity():
                    eta   = sif.non_resonant_dampingloss(omega)
                    i_sea = self._wave_DOF.index(syst.wave_DOF[0])     

                elif syst.isplate():
                    eta = sif.resonant_dampingloss(omega)
                    i_sea = self._wave_DOF.index(syst.wave_DOF[0])     

                self.SEAmatrix[i_sea,i_sea,:] = self.SEAmatrix[i_sea,i_sea,:].data +  m_dens[i_sea,:].data*eta
                
        # Loop over juncions using junction matrix JM
        for jun_key in self.junctions.keys():
            
            print('evaluating junction {0}'.format(jun_key))
            jun = self.junctions[jun_key]
            jm = jun.junction_matrix(omega)
            self.SEAmatrix += jm
            
        # Loop of hybrid junction (in most case one!)
        # Create input power load to SEA in case of external FE excitation 
        for jun_key in self.hybrid_junctions.keys():
            
            print('evaluating hybrid junction {0}'.format(jun_key))
            jun = self.hybrid_junctions[jun_key]
            # check if FE-system is excited by externl force
            if force in jun.fem.loads.keys():
                print('Apply FEM load:'+jun_key+' in hybrid Simulation response')
                femload_sw = True
                eta,eta_alpha,power_in,x_modal = jun.CLF(omega,force = force)
            else:
                femload_sw = False
                eta,eta_alpha          = jun.CLF(omega)
            
            # index of junction DOF in global SEA
            #id_global = self.wave_DOF.index(jun.wave_DOF)
            
            # do the FE damping
            #for ix,id_ in enumerate(id_global):
                # !-------------------- global index ------------------------------------!!-local index -!
                # self.SEAmatrix[id_,id_,:] = self.SEAmatrix[id_,id_,:].data + m_dens[id_,:]*eta_alpha[ix,:]
                
            
            # get all row and col dofs of upper triangular junction matrix
            jun_rows = jun.out_dofs
            jun_cols = jun.in_dofs
            # get indexes of jun_dof in overall wave_DOFs
            ir_global = self.wave_DOF.index(jun_rows)
            ic_global = self.wave_DOF.index(jun_cols)

            # and loop over dofs of the junction!
            for ix in range(len(jun_rows)):
                # determine ID of junction dofs
                ic = ic_global[ix] # index of currently treated column 
                ir = ir_global[ix]
                # add CLFs to the diagonal index ic,ic
                self.SEAmatrix[ic,ic,:] = self.SEAmatrix[ic,ic,:].data + m_dens[ir,:]*eta[ix,:]
                # and add to ofther diagonals ir,ir
                self.SEAmatrix[ir,ir,:] = self.SEAmatrix[ir,ir,:].data + m_dens[ir,:]*eta[ix,:]
                    
                # and fill the off diagonal use upper triangle, thus exchange ix_row and 
                self.SEAmatrix[ir,ic,:] = self.SEAmatrix[ir,ic,:].data - m_dens[ir,:]*eta[ix,:] 
        
 
            # create additional load case from FE excitation 
            if femload_sw:
                power_in_dof = dof.DOF(jun.wave_DOF.ID,jun.wave_DOF.dof,['power'])
                self.add_load('FE_'+force, lC.Load(self.xdata,power_in,power_in_dof,name='FE_'+force))
                
                # add hybrid result from power in.
                
                # create type of modal display resuls
                modal_dof = jun.fem.modal_dof(dof.DOFtype(typestr = 'displacement'))
                q_rms = jun.fem.rms_vec_from_modal(mC.Signal(self.xdata,x_modal,modal_dof))
                self.hybrid_result = q_rms*self.xdata
                
            
                        
            return 1
        
        
    def calculate_physical_units(self,**kwargs):
        """
        Create pysical units from energy results of system nature
        

        Parameters
        ----------
        **kwargs : dict
            Arbitrary keyword argmuments.

        

        Returns
        -------
        None.

        """
                
        i_sea = 0
        ydata =  np.zeros((self.N_sea,self.N_freq))       
        
        
        # loop over systems
        for i_sys,syst in enumerate(self.systems.values()):
            if syst.isSIF():
                pass
            else:
                # loop over wave_DOFs i.e. local SEA systems
                for ii,sea in enumerate(syst.wave_DOF):
                    #_DOF           = syst.wave_DOF[ii] # not needed
                    ydata[i_sea,:] = syst.physical_unit(self.xdata.angular_frequency, self.energy.ydata[i_sea,:], **kwargs)
                    i_sea += 1
               
        self.result = mC.Signal(self.xdata, ydata, self.wave_DOF)
                       
 
    def get_SEA_model_DOFs(self):
        """
        

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        """
        Generates the SEA DOFs of the full system setup
        
        These DOF are called wave_DOFs because they are not energy dofs but
        the according wave field, e.g. rms velocity and the ID describing the wave field 

        Returns
        -------
        dof
            wave_DOFs of SEA systems

        """
        
        # buffers of ID and DOF collection
        _ID  = np.zeros(self.N_sea,dtype=np.uint32)
        _dof = np.zeros(self.N_sea,dtype=np.uint32)
        _type = [0]*self.N_sea
        
        i_sea = 0
        
        
        # loop over systems
        for i_sys,syst in enumerate(self.systems.values()):
            
            # loop over wave_DOFs i.e. local SEA systems
            if not(syst.isSIF()):
                for ii,sea in enumerate(syst.wave_DOF):
                    _DOF           = syst.wave_DOF[ii]
                    _ID[i_sea]     = _DOF.ID[0]
                    _dof[i_sea]    = _DOF.dof[0]
                    _type[i_sea]   = _DOF.type[0]
                    i_sea += 1
                
        return dof.DOF(_ID,_dof,_type)
    
    @property
    def res_DOF(self):
        """
        property method for response dofs of SEA matrix (energy)

        Returns
        -------
        DOF
            reponse DOF of the SEA matrix
        """
        
        return dof.DOF(self.wave_DOF.ID,self.wave_DOF.dof,['energy'])
 
    @property        
    def exc_DOF(self):
        """
        property method for excitation dofs of SEA matrix (power)

        Returns
        -------
        DOF
            excitation DOF of the SEA matrix
        """
        
        return dof.DOF(self.wave_DOF.ID,self.wave_DOF.dof,['power'])

                
                
    def get_load_vector(self):
        """
        Creates load vector of SEA system

        Returns
        -------
        ndarray
            of power load

        """

        _Pi = np.zeros((self.N_sea,self.N_freq))
        
        for i_load,load in enumerate(self.loads.values()):
            # currently only power loads are consideres
            if load.dof.unique_type.typestr == 'power': # power
                i_sea = self.wave_DOF.find(load.dof.ID,load.dof.dof,repetition = False)
                _Pi[i_sea,:] += load.ydata
            elif load.dof.unique_type.typestr == 'force': # power
                # loop over force components
                for i_sig in range(load.Nsig):
                    if load[i_sig].dof.dof == 3:
                        _dof = load[i_sig].dof
                        # get vector index from DOF of excitation
                        i_sea = self.wave_DOF.find(_dof.ID,_dof.dof)
                        #i_sys = self.systems[load.dof.ID]
                              
                        ############## ID and system logics must be cleaned up!!!
                        F     = load[i_sig].ydata
                        omega = load[i_sig].xdata.angular_frequency 
    
                        _Pi[i_sea,:] = _Pi[i_sea,:] + self.systems[_dof.ID[0]].force_excitation_power(omega,F)
                
        return _Pi
    
    def get_modal_density(self):
        """
        modal densities in order of wave_DOF

        Returns
        -------
        m_dens : ndarray ( N_sea, N_freq)
            modal density of systems over frequency.

        """
        
        m_dens =  np.zeros((self.N_sea,self.N_freq))
        
        i_sea = 0
        
        # loop over sytems
        for syst in self.systems.values():

            if not(syst.isSIF()):            
                wave_dofs = syst.wave_DOF.dof
    
                for lw_dof in wave_dofs:
                    # get and prepare modal densities
                    m_dens[i_sea,:] = syst.modal_density(self.xdata.angular_frequency,lw_dof)
                    i_sea += 1
                
        return m_dens
    
    def get_CLF(self,m_dof,n_dof):
        """
        get sum of all CLF from off diagonal of existing SEA matrix 
        
        If two systems are connected via multiple paths, the SEA matrix contains
        the sum of all paths.
        For example two cavities separated by different panels contain all 
        sinlge CLFs of the panels in the off diagonal of the matrix

        Parameters
        ----------
        m_dof : dof.DOF 
            ... of from / source subsystem
        n_dof : dof.DOF 
            ... of to / receiver subsystem
 
        Returns
        -------
        n_eta : ndarry
            Sum of coupling loss factor

        """
        
        m = self.wave_DOF.index(m_dof)
        n = self.wave_DOF.index(n_dof)
        
        m_dens = self.get_modal_density()
        
        n_eta  = -1/m_dens[m]*self.SEAmatrix.data[n,m,:] #/self.xdata.angular_frequency
        
        return n_eta
        
    def solve(self,**kwargs):
        """
        Solve SEA model
        
        Solves the SEA model and provides the energy result for each subsystem
        In case of requested power the input of all power is

        Parameters
        ----------
        power_in: bool
            switch for input power calculation

        Returns
        -------
        None.

        """
        
        # for kw in kwargs:
        #     if kw=='loadresponse':
        #         lresp = True
                
        
        print("Solving SEA matrix")

        # set up load vector according to dofs
        Pi  = self.get_load_vector()
        E   = np.zeros((self.N_sea,self.N_freq),dtype = float )
        
        nn   = self.get_modal_density()    
        
        #loop over frequency
        for ifreq,omega in enumerate(self.xdata.data):
            self.SEAmatrix.Dindex(ifreq)
            E_per_n    = linalg.solve(omega*np.real(self.SEAmatrix.Dindex(ifreq)),Pi[:,ifreq])
            #print(Pi[:,ifreq])
            E[:,ifreq] = E_per_n*nn[:,ifreq]
        
        
        _dof = self.wave_DOF
        _edof = dof.DOF(_dof.ID,_dof.dof,dof.DOFtype(typestr = 'energy'))
        
        self.energy = mC.Signal(self.xdata,E,_edof)
        self.calculate_physical_units()
        
    def power_input(self,ID,wave_DOF = 0,result = 'detail'):
        """
        Calculates power input to specified subsystem
        

        Parameters
        ----------
        ID : int
            ID of systems.
        wave_DOF : int, optional
            wave_DOF of reverberant field. The default is 0 for all.
        result : str, optional
            Method of input calculation, 'detail' for single path contribution, ´total´ overall. The default is 'detail'.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        Signal
            power input of subsystems.

        """


        
        if ID in self.systems.keys():
            sys    = self.systems[ID]
            sif_sw = False
        elif ID in self.sifs.keys(): 
            sys    = self.sifs[ID]
            sif_sw = True
        else:
            raise ValueError('system or SIF of ID {0} not fount'.format(ID))

        
        omega = self.xdata.angular_frequency
 
        
        if result == 'total': # total power from system calculated not pathwise (if system is connectged via several junctions)
            source_dofs = self.wave_DOF
            in_dof   = sys.wave_DOF
            # find index of input system 
            ix_self = source_dofs.index(in_dof)
                 
            ix_sea  = np.arange(self.N_sea)
            # remoce index from possible input systems
            ix_sea  = ix_sea[ix_sea !=ix_self]
    
            # initialize loop indexes and buffers
            power_buf = []
            # Create index filter of DOFs
            ix_dof    = np.ones((self.N_sea,),dtype=bool)
            ix_dof[ix_self] = False
        
            # loop over sytems
            for ix in ix_sea:
                
                eta = self.get_CLF(source_dofs[ix],in_dof)
                
                if np.allclose(eta, 0.): 
                    # skip input power
                    ix_dof[ix] = False
                else:
                    power_buf.append(omega*eta*self.energy.ydata[ix,:])
     
            # cleanup and fill
            Nsig  = len(power_buf)
            power = np.zeros((Nsig,self.N_freq))
            ix = np.arange(self.N_sea)
            pow_in_dof = source_dofs[ix[ix_dof]]
            
            for i in range(Nsig):
                power[i,:] = power_buf[i]
                pow_in_dof[i].typestr = 'power'
        else: # pathwise calculation, single CLFs must be recalculated
            # specific treament of SIF. They have junction character
            if sif_sw:
                # loop over all subsystems connected to sif 
                power = np.zeros((sys.N,len(omega)))
                in_dof,_ = sys.get_wave_DOF()
                pow_in_dof = dof.DOF(in_dof.ID, in_dof.dof ,dof.DOFtype(typestr = 'power'))
                
                for ix,syst in enumerate(sys.systems):
                    if syst.iscavity():
                        i_sea = self.wave_DOF.index(syst.wave_DOF)
                        eta = sys.non_resonant_dampingloss(omega)
                    elif syst.isplate():
                        i_sea = self.wave_DOF.find(syst.ID,3)
                        eta = sys.resonant_dampingloss(omega)
                        
                    power[ix,:] = omega*eta*self.energy.ydata[i_sea,:]
            else:
                
                source_dofs = self.wave_DOF # all wave_dofs of the SEA model
                ix_self     = source_dofs.find(ID,wave_DOF) # index of system and wavedof in SEA model
                
                # initialize loop indexes and buffers
                power_buf = []
                IDs       = []
                dofs      = []
                              
                # loop over all junctions
                for junc_key in self.junctions.keys():
                    if sys in self.junctions[junc_key].systems:
                        jun = self.junctions[junc_key]
                        print ('Power input system found in junction {0}'.format(junc_key))
                        # get inxed of system in junction
                        ix_target_sys = jun.index(ID) # gives index into system list
                        # get all wave_dofs of junction
                        jun_dofs = jun.wave_DOF
                        # get indexes of jun_dof in overall wave_DOFs required for energy!
                        #ix_global = self.wave_DOF.index(jun_dofs)
                        # buffer for reusing values for
                        ix_local = jun_dofs.find(ID,dof=wave_DOF)
                        # index over all dofs of junction
                        ix_junc_in = np.arange(len(jun_dofs))
                        # remove index from possible input systems
                        ix_junc_in = ix_junc_in[ix_junc_in !=ix_local]
                        
                        for ix in ix_junc_in: # loop over all junction indexs except input
                            sID = jun_dofs[ix].ID
                            ix_source_sys = jun.index(sID)
                            s_wave_dof    = jun_dofs[ix].dof
                            eta = jun.CLF(omega,i_sys = (ix_source_sys,ix_target_sys),i_in_wave = (s_wave_dof,),i_out_wave = (wave_DOF,),Signal = False)
                            ix_global = self.wave_DOF.find(sID,dof=s_wave_dof)
                            E = self.energy.ydata[ix_global,:]                            
                            # source wave dofs and IDs to list
                            IDs.append(sID)
                            dofs.append(s_wave_dof)
                            power_buf.append(eta*omega*E)
                            
                # prepare details for the Signal
                pow_in_dof = dof.DOF(IDs,dofs,dof.DOFtype(typestr='power'))
                power = np.zeros((len(power_buf),len(omega)))
                for ix in range(len(power_buf)):
                    power[ix,:] = power_buf[ix]
                           
        return mC.Signal(self.xdata, power, pow_in_dof)
    
                          

        
        
        
        
        
        
            
            
                
                
                
                
            
            
            
            
            
            
        
                    

