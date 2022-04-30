"""
This module deals with the physics of junctions

In contrast to connections the junction classes deal with the physics
and the degrees of freedom/wavefields of the coupling.

The abstract junction class :class:`pyva.coupling.junctions.Junction` is 
extended by classes with specific geometry, namely:
    

- :class:`pyva.coupling.junctions.LineJunction` Plate line junction 
- :class:`pyva.coupling.junctions.AreaJunction` Cavity-plate junction
- :class:`pyva.coupling.junctions.HybridAreaJunction` Hybrid cavity-FEM(surface) junction
- :class:`pyva.coupling.junctions.SemiInfiniteFluid` Fluid half space junction 
    
"""



import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

import pyva.data.matrixClasses as mC
import pyva.data.dof as dof
import pyva.systems.acoustic3Dsystems as ac3Dsys
import pyva.systems.acousticRadiators as aR
import pyva.systems.structure2Dsystems as st2Dsys
import pyva.systems.infiniteLayers as iL
import pyva.models as mds

#import pyva.TMmodel as TMM
#import pyva.FEM
import pyva.useful as uf

# Debug switches used for detailed output
debug_sw = 0 # 2 global 3 amplitude 4 force_in 5 Pow_in 6 Pow_out
 

def edge_transform(theta):
    """
    helper function for edge coordinate transformation
    
    Calculates the coordinate transfer matrix for local edge coordinates 
    to global edge coordinates
    
    Parameters
    ----------
    theta : double
        plate angle

    Returns
    -------
    LinearMatrix 
        4x4 coordinate transformation matrix

    """

        
    cs = np.cos(theta)
    sn = np.sin(theta)
    
    return mC.LinearMatrix(np.array([[1., 0.,  0., 0.],
                                     [0., cs, -sn, 0.],
                                     [0., sn,  cs, 0.],
                                     [0., 0.,  0., 1.]],dtype = np.complex128))

def three_step(x,x1,x2):
    buf_ = np.ones(x.shape)
    buf_[x<x1] += 1
    buf_[x<x2] += 1
    return (buf_)
    
    

class Junction:
    """ 
    Junction defines the physics and systems of junctions
    
    Junction handle a specific set of dofs similar to the dynamix DOFS, but 
    the space direction and rotation DOFs are replaced by wave field DOFs for
    example logitudinal, shear and bending wave fields
    
    The DOF class used the DOF.ID for the identification of the system and the 
    wave_dof for the wave field that is considered as reverberant field
    
    Attributes
    ----------
    wave_DOF: DOF
        of the reverberant waves e.g. dof=3 and velocity for bending waves
    systems: list or tupel
        of systems
        
    """
    
    def __init__(self,systems):
        """
        Class contructor for junction
        
        Parameters
        ----------
        systems : tupel or list 
            systems connected by the junction

        Returns
        -------
        None.

        """
        
        self.systems = systems
        # get DOF and system pattern for junctin DOFs
        self.wave_DOF,self.wave_sys = self.get_wave_DOF()
        

    def __str__(self):
        """
        Implements __str__ 

        Returns
        -------
        text  : str
             Representing the junction.

        """
        
        _str = 'Junction with systems:\n'
        for sys in self.systems:
            _str = _str + '      {0}\n'.format(sys) 
        return _str
        
    def __repr__(self):
        """
        Implements __repr__ 

        Returns
        -------
        text  : str
             Representing the junction.

        """
        return self.__str__()
        
    @property    
    def N(self):
        """
        Number of connected physical SEA systems (not wave_fields)

        Returns
        -------
        int
            number of systems

        """
        return len(self.systems)
    
    @property
    def N_wave(self):
        """
        Number of wave fields 

        Returns
        -------
        int 
            Number of wave fields

        """
        
        return len(self.wave_DOF)
  
    def index(self,ID):
        """
        Determines index of ID in system list of junction

        Parameters
        ----------
        ID : int
            System ID

        Returns
        -------
        int
            index of ID in system list (-1) if not found

        """
        
        ix = -1
        
        for ix in range(self.N):
            if ID == self.systems[ix].ID:
                return ix
        return -1
        
    def get_wave_DOF(self,ix = slice(None,None)):
        """
        Provides the local wave DOFs of the junction
        
        Returns
        -------
        DOF,sys_list 
            dofs of junction, list of systems
        """
        
        N_sea = 0
        
        for syst in self.systems[ix]:
            N_sea += syst.N_wave_fields
            
        # buffers of ID and DOF collection
        _ID      = np.zeros(N_sea,dtype=np.uint32)
        _dof     = np.zeros(N_sea,dtype=np.uint32)
        _type    = [0]*N_sea
        sys_list = [0]*N_sea 
        
        # counter over wave fields
        i_sea = 0
        
        # loop over systems
        for i_sys,syst in enumerate(self.systems):
            
            # loop over wave_DOFs i.e. local SEA systems
            for ii,sea in enumerate(syst.wave_DOF):
                _DOF           = syst.wave_DOF[ii]
                _ID[i_sea]     = _DOF.ID[0]
                _dof[i_sea]    = _DOF.dof[0]
                _type[i_sea]   = _DOF.type[0]
                sys_list[i_sea]      = syst
                i_sea += 1
                
        return (dof.DOF(_ID,_dof,_type),sys_list)
    
    @property
    def res_DOF(self):
        """
        Provides 'response' DOFs of junction SEA matrix

        The 'response' and excitation DOFs are equal in terms of ID and 
        wave_dof but use 'energy' as physical output

        Returns
        -------
        DOF
            response DOFs of junction SEA matrix

        """
        
        return dof.DOF(self.wave_DOF.ID,self.wave_DOF.dof,['energy'])
 
    @property        
    def exc_DOF(self):
        """
        Provides 'excitation' DOFs of junction SEA matrix

        The 'excitation' and response DOFs are equal in terms of ID and 
        wave_dof but use 'power' as physical input

        Returns
        -------
        DOF
            excitation DOFs of junction SEA matrix

        """
        return dof.DOF(self.wave_DOF.ID,self.wave_DOF.dof,['power'])
 
    
    def modal_density(self,omega):
        """
        Provides modal density
        
        Calculates the modal density for all connected wavefields of the 
        junction.

        Parameters
        ----------
        omega : ndarray
            angular frequency.

        Returns
        -------
        res : ndarray of dimension N_wave x size(omega)
            modal density.

        """
        
        res = np.zeros((self.N_wave,omega.size))
        
        i_sea = 0
        for ix,syst in enumerate(self.systems):
            for iw,wdof in enumerate(self.wave_DOF[ix]):
                res[i_sea,:] = syst.modal_density(omega,wave_DOF = wdof.dof[iw])
                i_sea += 1
                
        return res
        
    
    @property
    def in_dofs(self):
        """
        Provides input wave dofs of junction for upper triangular input matrix
        
        These DOFs are requiered to determine the column index of 
        the upper triangular components in the global SEA matrix.  
        
        For example a plate - cavity junction has wave dofs [ ID=1,dof=3 ; 1,5 ; 2,0 ]
        
        The full of diag CLF in the SEA matrix are
        
        .. math:: 
            \\begin{bmatrix}
                \\Sigma & - n_{15}\\eta_{15,13} & -n_{20}\eta_{20,13} \\\\ \
                - n_{13}\\eta_{13,15} & \\Sigma & -n_{20}\\eta_{20,15} \\\\ \
                - n_{13}\\eta_{13,20} & - n_{15}\\eta_{15,20} & \\Sigma 
            \\end{bmatrix}
            
        The column dofs would be : [ 1,5 ; 2,0 ; 2,0]
        
        See also
        --------
        out_dofs
        
        Returns
        -------
        DOF
            wave_DOF of column.

        """
        
        i_rows,i_cols = np.triu_indices(self.N,1)
        return self.wave_DOF[i_cols]
    
    @property    
    def out_dofs(self):
        """
        Provides output wave dofs of junction for upper triangular input matrix
        
        These DOFs are requiered to determine the row index of 
        the upper triangular components in the global SEA matrix.  
        
        For example a plate - cavity junction has wave dofs [ 1,3 ; 1,5 ; 2,0 ]
        
        The full of diag CLF in the SEA matrix are
        
        .. math:: 
            \\begin{bmatrix}
                \\Sigma & - n_{15}\\eta_{15,13} & -n_{20}\eta_{20,13} \\\\ \
                - n_{13}\\eta_{13,15} & \\Sigma & -n_{20}\\eta_{20,15} \\\\ \
                - n_{13}\\eta_{13,20} & - n_{15}\\eta_{15,20} & \\Sigma 
            \\end{bmatrix}

        The column dofs would be : [ 1,3 ; 1,3 ; 1,5]
        
        See also
        --------
        in_dofs
        
        Returns
        -------
        DOF
            wave_DOF of rows.

        """

   



        i_rows,i_cols = np.triu_indices(self.N,1)
        return self.wave_DOF[i_rows]

    def junction_matrix(self,omega):
        """
        Create empty junction matrix for further calculation in daughter classes  
    
        Parameters
        ----------
        omega : float or ndarray
            angular frequency.
    
        Returns
        -------
        DynamicMatrix
            Empty junction matrix with correct in- and out-DOFs 
    
        """
        
        if uf.isscalar(omega):
            N_freq = 1
        else:
            N_freq = len(omega)
        
        data = np.zeros(((self.N_wave*(self.N_wave+1))//2,N_freq))
                        
        return mC.DynamicMatrix(data, mC.DataAxis(omega), self.exc_DOF, self.res_DOF, sym = 1, shape = (self.N_wave,self.N_wave,N_freq))

        
def all2array(*A):
    """
    Converts input to np.arrays

    Parameters
    ----------
    *A : list of tupel, np.arrays or lists

    Returns
    -------
    r : tuple of np.arrays of input A

    """
    
    r = ()
    
    for a in A:
        r += np.array(a),
    return r
        
class LineJunction(Junction) :
    """ 
    LineJunction defines the physics based on the connected systems of junctions 

    Attributes
    ----------
    systems : list or tuple
        of systems
    length : float
        length of LineJunction
    thetas : ndarray
        angles of connected plates
    
    """
 
       
    
    def __init__(self,systems,length,thetas):
        """
        Class contructor for LineJunction
        
        Parameters
        ----------
        systems : list or tuple
            systems connected to junction
        length : float
            length of LineJunction
        thetas : list or tuple
            angles of connected plates
                
                        
        """
        
        # Check if systems are allowed structur2DSystem sytems
        for sys_ in systems:
            if not(isinstance(sys_,(st2Dsys.Structure2DSystem,) )):
                raise ValueError('Line junction system must be of type structure2Dsytem or structure1Dsytem.')
               
        
        super().__init__(systems)
        self.length = length
        self.thetas = thetas
        # buffers fpr performance increas
        self._CLFs = dict()
        self._omega = dict()
        
    def __str__(self):
        """
        Implements __str__ 

        Returns
        -------
        text  : str
             Representing the line junction.

        """
        
        _str = 'LineJunction with systems:\n'
        for i_,sys in enumerate(self.systems):
            _str += '{0} angle: {1:.4f}\n'.format(sys,self.thetas[i_])
            
        _str += 'length       : {0}\n'.format(self.length) 
        
        return _str
        
    def __repr__(self):
        """
        Implements __repr__ 

        Returns
        -------
        text  : str
             Representing the junction.

        """
        return self.__str__()
                
    def total_radiation_stiffness_wavenumber(self,omega,wavenumber):
        """
        Calculates the total radiation stiffness of line junctions
        
        The result is given in wavenumber domain and global displacement 
        coordinate system 
        
        Parameters
        ----------
        omega : float
            angular frequency.
        wavenumber : ndarray
            wavenumber in x- (edge) direction.

        Returns
        -------
        LinearMatrix
            Total stiffness matrix of line junction 

        """
        
                
        # # manage datastorage
        # if omega in self._D_tot.keys():
        #     # check is frequency parameters have not changed
        #     if np.allclose(self._wavenumber[omega],wavenumber,atol = 1e-16):
        #         return self._D_tot[omega]
            
        # start with the first plate
        D_tot = self.systems[0].edge_radiation_stiffness_wavenumber(omega,wavenumber)
        # Rotation matrix from egde coordinate transformation
        T_rot = edge_transform(self.thetas[0])
        # Transform to global system
        D_tot = (T_rot.dot(D_tot)).dot(T_rot.transpose())
        
        for i_sys in range(1,self.N):
            _D_buf = self.systems[i_sys].prop.edge_radiation_stiffness_wavenumber(omega,wavenumber)
            # Rotation matrix from egde coordniate transformation
            T_rot = edge_transform(self.thetas[i_sys])
            
            _T_buf = (T_rot.dot(_D_buf)).dot(T_rot.transpose())
            
            D_tot += _T_buf
            
        #self._wavenumber.update({omega : wavenumber})
        #self._D_tot.update({omega : D_tot})
        
        return(D_tot)
    
    def transmission_wavenumber(self,omega,wavenumber,i_sys = (0,1),i_in_wave = (1,2,3),i_out_wave = (1,2,3),rad_sw = 'wave',Signal = True):
        """
        Calculates the transmission coefficient of line junctions
        
        This methods applies the hybrid CLF formulation from [1] but
        using the radiated power calcualted from the wave amplitude as 
        disscussed in the reference.
        
        The rad_sw argument is used for validation. For simulation always use 'wave'

        Parameters
        ----------
        omega : float or ndarray
            angular frequency.
        wavenumber : float or ndarray
            wavenumber in x- or edge-direction.
        i_sys : tuple or list, optional
            incident and radiating system index vector. The default is (0,1).
        i_in_wave : tuple or list, optional
            incident wave type. The default is (1,2,3).
        i_out_wave : tuple or list, optional
            radiating wave type. The default is (1,2,3).
        rad_sw : str, 'wave' or 'im_dir', optional
            identifier for radiated power method. The default is 'wave'.
        Signal : bool, optional
            switch for Signal output. The default is True.

        Raises
        ------
        ValueError
            When input values are not or correct type of not constent.

        Returns
        -------
        ndarray or Signal
            transmission coefficient for one frequency over wavenumber.

        """
        
        i_sys,i_in_wave,i_out_wave = all2array(i_sys,i_in_wave,i_out_wave)

        if i_out_wave.size != i_in_wave.size:
            raise ValueError('In and out wave must have same size')
   
        Nsig = len(i_in_wave)
        
        # convert to np.array to become independent from tupels
        i_sys,i_in_wave,i_out_wave = all2array(i_sys,i_in_wave,i_out_wave)

        D_tot = self.total_radiation_stiffness_wavenumber(omega,wavenumber)
   
        # Determine transformation matrices for in- and out-systems
        T_rot_1  = edge_transform(self.thetas[i_sys[0]])
        T_rot_1_T = T_rot_1.transpose()
        T_rot_2  = edge_transform(self.thetas[i_sys[1]])
        T_rot_2_T = T_rot_2.transpose()
        
        # Determine wave transformation for output wave amplitude determiniation
        T_wave_2 = self.systems[i_sys[1]].wave_transformation_matrix(omega,wavenumber)
        #T_wave_2_inv = self.systems[i_sys[1]].wave_transformation_matrix(omega,wavenumber,inv=True)

        D_mn_part = T_rot_1_T.dot(D_tot).dot(T_rot_2)
        
        Nsig = len(i_in_wave)
        _ydata= np.zeros((Nsig,np.size(wavenumber)))
        
        _tdof = dof.DOFtype(typestr='transmission')
        xdata = mC.DataAxis(wavenumber,typestr='wavenumber')
        
        
        kL = self.systems[i_sys[0]].prop.wavenumber_L(omega) #longitudinal wavenumber
        kS = self.systems[i_sys[0]].prop.wavenumber_T(omega) #shear wavenumber
        kB = self.systems[i_sys[0]].prop.wavenumber_B(omega) #shear wavenumber
        B  = self.systems[i_sys[0]].prop.B
        S  = self.systems[i_sys[0]].prop.S
    
    
        for i_in in range(Nsig):

            ii_wave = i_in_wave[i_in]
            # I suppose this is not used anymore
            if ii_wave == 0: # all S,L and B
                D_in  = self.systems[i_sys[0]].edge_imaginary_radiation_stiffness_wavenumber(omega,wavenumber,0)
                den_ = three_step(wavenumber,kL,kS)

            else:    
                D_in  = self.systems[i_sys[0]].edge_imaginary_radiation_stiffness_wavenumber(omega,wavenumber,ii_wave)
                den_ = np.ones(wavenumber.shape)

            Sqqe = D_in.HDH(D_mn_part)

            # Remove incoming motion if same system
            # Meaning in principle that we have to rely on the 
            if i_sys[0]==i_sys[1]: # and io_wave != ii_wave:
                # Here the single result is required
                
                # derive correction factor to reduce from Im Ddir normalisation
                if ii_wave == 1:
                    kyL = np.sqrt(kL**2 - wavenumber**2)
                    cor_fac = (2*np.sqrt(S*kyL)*kS) #np.real
                elif ii_wave == 2:
                    kyS = np.sqrt(kS**2 - wavenumber**2)
                    cor_fac = (2*np.sqrt(S*kyS)*kS)
                elif ii_wave in (3,4):
                    kyB = np.sqrt(kB**2 - wavenumber**2)
                    cor_fac = (2*np.sqrt(2*B*kyB)*kB)
                elif ii_wave == 5: # I suppose this does not work with a combined correction factor
                    kyS = np.sqrt(kS**2 - wavenumber**2)
                    kyL = np.sqrt(kL**2 - wavenumber**2)
                    cor_facS = (2*np.sqrt(S*kyS)*kS)
                    cor_facL = (2*np.sqrt(S*kyL)*kS)
                    cor_fac = cor_facS+cor_facL
                    
                    
                if ii_wave == 5:
                    # force from incoming waves
                    f_Sin  = self.systems[i_sys[0]].edge_wave_excitation_force(omega,wavenumber,2,True)/cor_facS
                    f_Lin  = self.systems[i_sys[0]].edge_wave_excitation_force(omega,wavenumber,1,True)/cor_facL
                    # global junction displacement due to f_in
                    qL0 = D_tot.solve(T_rot_1.dot(f_Lin))
                    qS0 = D_tot.solve(T_rot_1.dot(f_Sin))
                    # edge junction diplacement due to f_in
                    qLe = T_rot_2_T.dot(qL0)
                    qSe = T_rot_2_T.dot(qS0)
    
                    # edge displaqcement in local edge 1
                    qLe_in = self.systems[i_sys[0]].edge_wave_excitation_displacement(omega,wavenumber,1)/cor_facL
                    qSe_in = self.systems[i_sys[0]].edge_wave_excitation_displacement(omega,wavenumber,2)/cor_facS
                   
                    # Matrix version of incoming wave correction
                    # Correction follows from (qe-qe_in)*(qe-qe_in)^H
                    SqqeredL = -(qLe_in.dot(qLe.H()))-(qLe.dot(qLe_in.H()))+(qLe_in.dot(qLe_in.H()))
                    SqqeredS = -(qSe_in.dot(qSe.H()))-(qSe.dot(qSe_in.H()))+(qSe_in.dot(qSe_in.H()))
                    #Sqqered = -2*(qe_in.dot(qe.H())).real()+(qe_in.dot(qe_in.H()))
                    #MMred   = Sqq0red.HDH(T_rot_wave_2)
                    ix = wavenumber > kL
                    SqqeredL._data[:,ix] = 0. # set single S wave area to 2
                    Sqqe += SqqeredS
                    Sqqe += SqqeredL
                    
                else:      
                    # force from incoming wave
                    f_in  = self.systems[i_sys[0]].edge_wave_excitation_force(omega,wavenumber,ii_wave,True)/cor_fac
                    # global junction displacement due to f_in
                    q0 = D_tot.solve(T_rot_1.dot(f_in))
                    # edge junction diplacement due to f_in
                    qe = T_rot_2_T.dot(q0)
    
                    # edge displaqcement in local edge 1
                    qe_in = self.systems[i_sys[0]].edge_wave_excitation_displacement(omega,wavenumber,ii_wave)/cor_fac
                   
                    # Matrix version of incoming wave correction
                    # Correction follows from (qe-qe_in)*(qe-qe_in)^H
                    Sqqered = -(qe_in.dot(qe.H()))-(qe.dot(qe_in.H()))+(qe_in.dot(qe_in.H()))
                    #Sqqered = -2*(qe_in.dot(qe.H())).real()+(qe_in.dot(qe_in.H()))
                    #MMred   = Sqq0red.HDH(T_rot_wave_2)
                    Sqqe += Sqqered
                
            # Transform to Psi...
            Sqq_psi = Sqqe.HDH(T_wave_2)
            #MM = D_in.HDH.(D_mn_tot)

                
            #print('index {0}'.format(i_out))
            io_wave = i_out_wave[i_in]
            
            
            if io_wave == 3:
                io_wave = 4
                
            if rad_sw == 'im_dir': # just for presentation purpose works for Bending or SL i_in = 5 for k<= kL but not for kL < k < kS
                if io_wave == 0: # all not finally tested
                    D_out = self.systems[i_sys[1]].edge_imaginary_radiation_stiffness_wavenumber(omega,wavenumber,0)
                else:
                    # version delivering the imag dircetly
                    D_out = self.systems[i_sys[1]].edge_imaginary_radiation_stiffness_wavenumber(omega,wavenumber,io_wave)   
                    # version using the full stiff and perform imag later
                    D_out = self.systems[i_sys[1]].edge_radiation_stiffness_wavenumber(omega,wavenumber,io_wave)   
                    #D_out = D_out.imag()

                #_ydata[i_in*Nout+i_out,:] = (4*D_out.imag()*Sqqe).sum()# Imag over all
                _ydata[i_in,:] = (4*D_out*Sqqe).imag().sum() # Imag seperately assuming symmetry of D_out
                    
            else:

                if io_wave == 5: # all
                    Psi2 = np.abs(Sqq_psi.data[0,0,:])
                    Psi21 = np.abs(Sqq_psi.data[1,1,:])
                else:    
                    Psi2 = np.abs(Sqq_psi.data[io_wave-1,io_wave-1,:])



                if io_wave == 5: # all
                    # add in-plane waves
                    WQ0 = self.systems[i_sys[1]].edge_wave_amplitude_radiated_power(1.,omega,wavenumber,1)
                    WQ1 = self.systems[i_sys[1]].edge_wave_amplitude_radiated_power(1.,omega,wavenumber,2)
                    _ydata[i_in,:] = 8/omega*(np.abs(WQ0*Psi2)+np.abs(WQ1*Psi21))
                else:
                    WQ = self.systems[i_sys[1]].edge_wave_amplitude_radiated_power(1.,omega,wavenumber,io_wave)
                    _ydata[i_in,:] = 8/omega*WQ*Psi2 # Imag seperately assuming symmetry of D_out
                #tot_imag = np.sum(np.imag((4*D_out*MM).imag().sum()))
                #print('Total imaginary values {0:f}'.format(tot_imag))
        if Signal:
            return mC.Signal(xdata,_ydata,dof.DOF(i_out_wave,np.zeros((1,Nsig)),_tdof))
        else:
            return _ydata
            


    def transmission_wavenumber_wave(self,omega,wavenumber,i_sys = (0,1),i_in_wave = (1,)*3+(2,)*3+(3,)*3,i_out_wave = (1,2,3)*3,matrix = False):
        """
        transmission coefficient assuming wave transformations for radiation stiffness 
        calculations blocked forces from langley are used
        

        Parameters
        ----------
        omega : TYPE
            DESCRIPTION.
        wavenumber : TYPE
            DESCRIPTION.
        i_sys : TYPE, optional
            DESCRIPTION. The default is (0,1).
        i_in_wave : TYPE, optional
            DESCRIPTION. The default is (1,)*3+(2,)*3+(3,)*3.
        i_out_wave : TYPE, optional
            DESCRIPTION. The default is (1,2,3)*3.
        matrix : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        None.

        """
        
        print("DON'T USE! FOR DEMONSTRATION PURPOSE ONLY!")
        
        D_tot = self.total_radiation_stiffness_wavenumber(omega,wavenumber)
        #D_tot_inv = D_tot.inv()        
        
        #T_wave_1 = self.systems[i_sys[0]].wave_transform(omega,wavenumber)
        T_wave_2 = self.systems[i_sys[1]].wave_transform(omega,wavenumber)
        #T_wave_1_inv = T_wave_1.inv()
        T_wave_2_inv = self.systems[i_sys[1]].wave_transform(omega,wavenumber,inv=True)
        
        #D_dir_1_edge = self.systems[i_sys[0]].radiation_stiffness_wavenumber(omega,wavenumber)
        #D_dir_2_edge = self.systems[i_sys[1]].radiation_stiffness_wavenumber(omega,wavenumber)
        #D_dir_1_wave = T_wave_1_inv.dot(D_dir_1_edge).dot(T_wave_1)
        #D_dir_2_wave = T_wave_2_inv.dot(D_dir_2_edge).dot(T_wave_2)
            
        
        T_rot_1  = edge_transform(self.thetas[i_sys[0]])
        T_rot_1_T = T_rot_1.H() #transpose()
        T_rot_2  = edge_transform(self.thetas[i_sys[1]])
        T_rot_2_T = T_rot_2.H() #transpose()
        

        # Again radiation in wave coordinates
        #Dtot_interim = T_rot_1_T.dot(D_tot).dot(T_rot_2)
        #Dtot_interim0 = T_rot_1_T.dot(D_tot)
        #Dtot_        = Dtot_interim0.dot(T_rot_2).dot(T_wave_2)

        #D_mn_wave = D_mn


        Nin   = len(i_in_wave)
        Nout  = len(i_out_wave)
        Nsig  = Nin*Nout
        
        _ydata= np.zeros((Nsig,np.size(wavenumber) ))
        
        _tdof = dof.DOFtype(typestr='transmission')
        xdata = mC.DataAxis(wavenumber,typestr='wavenumber')

        plate1 = self.systems[i_sys[0]] 
        plate2 = self.systems[i_sys[1]] 
        rho_area = plate1.prop.mass_per_area
        
        
        for i_in in range(Nin):

            ii_wave = i_in_wave[i_in]
            if ii_wave == 2:
                    ii_wave = 3            




            Sff    = plate1.wave_excitation_force_cross_correlation(omega,wavenumber,ii_wave,matrix)
            
            # Plate forces in global coordinates
            Sff0   = T_rot_1.dot(Sff).dot(T_rot_1_T)
            # The use of D_tot_inv is not precise enough, so we switch to soilve
            #Sqq0   = D_tot.solve(Sff0)
            #Sqq0   = D_tot_inv.dot(Sff0).dot(D_tot_inv.H())
            Sqq0 = Sff0.HDH(D_tot)
            
            # Remove incoming motion if same system
            if i_sys[0]==i_sys[1]:
                # Here the single result is required
                f_in  = self.systems[i_sys[0]].wave_excitation_force(omega,wavenumber,ii_wave,matrix)
                q0    = D_tot.solve(T_rot_1.dot(f_in))
                q_in = self.systems[i_sys[0]].wave_excitation_displacement(omega,wavenumber,ii_wave)
                q_in = T_rot_1.dot(q_in)

                # Reflexion means coherent wave that must be removed
                Sqq0 += -(q_in.dot(q0.H()))-(q0.dot(q_in.H()))+(q_in.dot(q_in.H()))  # evtl wegen der -nomenklatur 
            #Sqq0   = D_tot_inv.dot(Sff0).dot(D_tot_inv.H())
            #Sqq0  = D_mn_interim0.dot(Sff).dot(D_mn_interim0.H())
            #Sqqe  = D_mn_interim.dot(Sff).dot(D_mn_interim.H())
            Sqqe  = T_rot_2_T.dot(Sqq0).dot(T_rot_2)
            #SqqPsi2 = T_wave_2_inv.dot(Sqqe).dot(T_wave_2_inv.H())
            SqqPsi = T_wave_2_inv.dot(Sqqe).dot(T_wave_2_inv.H())

            #i_test = 30
            
            if debug_sw in (3,5):
                plt.figure(100+ii_wave)
                
            for i_out in range(Nout):
                

                io_wave = i_out_wave[i_out]
                if io_wave == 2:
                    io_wave = 3
                    
                Psi     = SqqPsi.data[io_wave,io_wave,:]    
                #Psi      = np.abs(buf).flatten()
                Psi_comp = np.abs(Psi).flatten()
                
                pow_unit = plate2.wave_amplitude_radiated_power(np.sqrt(Psi_comp),omega,wavenumber,io_wave).flatten()
                
                if debug_sw == 6:
                    #plt.figure(300+ii_wave)
                    plt.plot(wavenumber,pow_unit,label = 'Pow in DOF:{0}'.format(io_wave))
                elif debug_sw == 3: 
                    plt.plot(wavenumber,Psi, label = 'Re {0}to{1}plate{2}'.format(ii_wave,io_wave,i_sys[1]))
#                    plt.plot(wavenumber,np.real(Psi), label = 'Re {0}to{1}plate{2}'.format(ii_wave,io_wave,i_sys[1]))
#                    plt.plot(wavenumber,np.imag(Psi), label = 'Im {0}to{1}plate{2}'.format(ii_wave,io_wave,i_sys[1]))


                if ii_wave == 0:
                    k_plate = np.real(plate1.prop.wavenumber_L(omega))
                    c_gr    = plate1.prop.c_L()
                elif ii_wave == 1:
                    k_plate = np.real(plate1.prop.wavenumber_T(omega))
                    c_gr = plate1.prop.c_T()
                elif ii_wave in (2,3):
                    k_plate = np.real(plate1.prop.wavenumber_B(omega))
                    c_gr = plate1.prop.c_B_group(omega)
                else:
                    raise ValueError('i_wave argument must be in [0,1,2]')

                phi = np.arccos(wavenumber/k_plate)
                sinphi = np.sin(phi)

                if ii_wave in (0,1):
                    _ydata[i_in*Nout+i_out,:] = 2*pow_unit/(k_plate**2*rho_area*c_gr*omega**2*sinphi)
                elif ii_wave in (2,3):
                    _ydata[i_in*Nout+i_out,:] = 2*pow_unit/(rho_area*c_gr*omega**2*sinphi)

            plt.legend()
            if debug_sw == 3:
                plt.title('Amplitude due to excitation from {0}'.format(ii_wave))
            if debug_sw == 6:
                plt.title('Out power due to wave {0}'.format(ii_wave))


        return mC.Signal(xdata,_ydata,dof.DOF(np.array(i_out_wave),np.zeros((1,Nsig)),_tdof))
     
                
   
                
    def transmission_wavenumber_langley(self,omega,wavenumber,i_sys = (0,1),i_in_wave = (1,)*3+(2,)*3+(3,)*3,i_out_wave = (1,2,3)*3,matrix=False,Signal = True):
        """
        transmission coefficient based on the wave transionmssion method
        
        The version is based on the wave trasmission theory from [2]
        The transmission is calculated for each specific wave type, this
        is not possible for the hybrid CLF as the diffuse field reciprocity
        is not valid for single in-plane waves.
        
        The matrix argument is for test pupose to check if the analytical expression 
        is correct.
        
        Parameters
        ----------
        omega : float or ndarray
            angular frequency.
        wavenumber : float or ndarray
            wavenumber in x- or edge-direction.
        i_sys : tuple or list, optional
            incident and radiating system index vector. The default is (0,1).
        i_in_wave : tuple or list, optional
            incident wave_DOF. The default is (1,1,1,2,2,2,3,3,3).
        i_out_wave : tuple or list, optional
            transmissted wave_DOF. The default is (1,2,3)*3.
        matrix : bool, optional
            Switch if the matrix D.f - q (True) or the expicity analytical version is used. The default is False.
        Signal : TYPE, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        ndarray
            transmission coefficient.

        """
      
        
        # convert to np.array to become independent from tupels
        # i_sys      = np.array(i_sys)
        # i_in_wave  = np.array(i_in_wave)
        # i_out_wave = np.array(i_out_wave)
        i_sys,i_in_wave,i_out_wave = all2array(i_sys,i_in_wave,i_out_wave)

        
        D_tot = self.total_radiation_stiffness_wavenumber(omega,wavenumber)
        D_tot_inv = D_tot.inv()        
        
        #T_wave_1 = self.systems[i_sys[0]].wave_transform(omega,wavenumber)
        #T_wave_2 = self.systems[i_sys[1]].wave_transform(omega,wavenumber)
        #T_wave_1_inv = T_wave_1.inv()
        T_wave_2_inv = self.systems[i_sys[1]].wave_transformation_matrix(omega,wavenumber,inv=True)
        
                    
        T_rot_1  = edge_transform(self.thetas[i_sys[0]])
        #T_rot_1_T = T_rot_1.transpose()
        T_rot_2  = edge_transform(self.thetas[i_sys[1]])
        T_rot_2_T = T_rot_2.transpose()
        
        
        
        D_mn  = T_wave_2_inv.dot(T_rot_2_T).dot(D_tot_inv).dot(T_rot_1) 
        #D_mn  = (T_rot_1.transpose().dot(D_tot).dot(T_rot_2).dot(T_wave_2)).inv()
        
        # Number of input and output Signals
        Nsig  = len(i_in_wave)
        
        _ydata= np.zeros((Nsig,np.size(wavenumber) ))
        
        _tdof = dof.DOFtype(typestr='transmission')
        xdata = mC.DataAxis(wavenumber,typestr='wavenumber')
        
        
        for i_in in range(Nsig):

            ii_wave = i_in_wave[i_in]
            
            #k_in = self.systems[i_sys[0]].plate_wavenumber(omega,ii_wave)

            #if ii_wave in (3,4):
            #    k_in = 1.
                        
            f_in  = self.systems[i_sys[0]].edge_wave_excitation_force(omega,wavenumber,ii_wave,matrix)
            #f_in0 = T_rot_1.dot(f_in)
            
            # Response in global coordinates, use solve for higher precision
            _buf = D_tot.solve(T_rot_1.dot(f_in))
            
            
            # Remove incoming motion if same system
            if i_sys[0]==i_sys[1]:
                _q = self.systems[i_sys[0]].edge_wave_excitation_displacement(omega,wavenumber,ii_wave)
                _buf -= T_rot_1.dot(_q) # evtl wegen der -nomenklatur 

            Pow_in   = self.systems[i_sys[0]].edge_wave_amplitude_radiated_power(1.,omega,wavenumber,ii_wave)
            
                
            # from global to edge 2
            qe       = T_rot_2_T.dot(_buf)
            # from edge 2 to psi
            wave_out = T_wave_2_inv.dot(qe)

            io_wave = i_out_wave[i_in]

            if io_wave == 3:
                io_psi = 3
            else:
                io_psi = io_wave-1 # because index 0,1,2,3 corresponds to 1,2,(not used),(3 and 4)
                
            Psi_out = wave_out[io_psi,0,:].data.flatten()
            Pow_out = self.systems[i_sys[1]].edge_wave_amplitude_radiated_power(Psi_out,omega,wavenumber,io_wave)
            
            #print('index {0}'.format(i_out))
            _ydata[i_in,:] = Pow_out/Pow_in

        #return mC.Signal(xdata,_ydata,dof.DOF(1+np.arange(Nsig),np.zeros((1,Nsig)),_tdof))
    
        if Signal:
            return mC.Signal(xdata,_ydata,dof.DOF(i_out_wave,np.zeros((1,Nsig)),_tdof))
        else:
            return _ydata
    
    
    def transmission_wavenumber_diffuse(self,omega,i_sys = (0,1),i_in_wave = (5,5),i_out_wave = (5,3), N_step = 100,method = 'diffuse',Signal = True):
        """
        diffuse transmission coefficient for line junctions
        

        Parameters
        ----------
        omega : float or ndarray
            angular frequency.
        i_sys : tuple or list, optional
            incident and radiating system index vector. The default is (0,1).
        i_in_wave : tuple or list, optional
            incident wave_DOF. The default is (5,5).
        i_out_wave : tuple or list, optional
            transmissted wave_DOF. The default is (5,3).
        N_step : int, optional
            Number of wavenumber integration steps. The default is 100.
        method : str, optional
            Mewthod seclector 'diffuse' for diffuse field reciprocity, 'langley' for wave transmission. The default is 'diffuse'.
        Signal : bool, optional
            Swith for return datatype. The default is True.

        Returns
        -------
        ndarray (Signal = False) or Signal
            diffuse field transmission coefficient.
        """

                
        # convert to np.array to become independent from tupels
        i_sys,i_in_wave,i_out_wave = all2array(i_sys,i_in_wave,i_out_wave)

        tau = np.zeros((len(i_out_wave),omega.size))


        
        for ifreq,om in enumerate(omega):
            
            kx   = self.kx(om,i_sys,i_in_wave[0],i_out_wave,N_step,method = 'in-plane')     
            k_in = self.systems[i_sys[0]].plate_wavenumber(om,i_in_wave[0])
            if method == 'diffuse':
                taus = self.transmission_wavenumber(om,kx,i_sys,i_in_wave,i_out_wave,Signal=False)#.ydata # check late with Signal option
            elif method == 'langley':
                taus = self.transmission_wavenumber_langley(om,kx,i_sys,i_in_wave,i_out_wave,Signal=False)#.ydata
            
            for i_wave,i_type in enumerate(i_out_wave):
                tau[i_wave,ifreq] = integrate.trapz(taus[i_wave,:],kx)/k_in

        if Signal:
                        
            xdata = mC.DataAxis(omega,typestr='angular frequency')
            tdof = dof.DOFtype(typestr='transmission')
            return mC.Signal(xdata,tau,dof.DOF((i_sys[1],),i_out_wave,tdof,repetition=True))
        else:
            return tau        


    
    def max_inplane_wavenumber(self,omega):
        """
        Determine the maximal in-plane wavenumber of all connected plates
        
        Bending wavenumber can be orders of magnitude higher. In this case the in-plane
        wave may be undersampled. Thus, for all wavenumber intergration the upper boundary
        of the in-plane wavenumber is required for better sampling of 0 > kx > kS. 
        
        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        k_max : float
            maximal wavenumber.
        
        
        """
        
        k_max = 0.
        
        for i_sys,syst in enumerate(self.systems):
            if k_max < syst.plate_wavenumber(omega,2): # 2 for shear!
                k_max = syst.plate_wavenumber(omega,2)
                
        return k_max
        
        

    def kx(self,omega,i_sys,i_in_wave,i_out_wave,Nstep,method = 'angle'):
        """
        wavenumber sampling for diffuse field intgegration
        
        The different wave number regimes require a well balanced sampling of the 
        wavnumber. High sampling in-plane waves and low sampling for bending waves
        
        Parameters
        ----------
        omega: double
            angular frequency
        i_sys : tupel or list of int
            system index
        i_in_wave : tupel or list of int
            wave_DOF of irradiating wave field
        i_out_wave : tupel or list of int
            wave_DOF of radiating wave field
        Nstep : int
            Integration steps
        method : str
            'angle' for angular integration with dense sampling at 90 degrees
            'in-plane' for two intervals, kx = [0 kS ) + [kS kB) 

        Returns
        -------
        kx : ndarray
            wavenumber in x

        """
        
        k_in = self.systems[i_sys[0]].plate_wavenumber(omega,i_in_wave)
        
        if method == 'angle':
            phi = np.linspace(0,0.99/2*np.pi,Nstep)
            kx  = k_in*(1-np.cos(phi))
        elif method == 'in-plane':
            k_max = self.max_inplane_wavenumber(omega)
            if k_max < k_in: 
                kx1 = np.linspace(0,k_max,Nstep)
                kx2 = np.linspace(k_max,k_in,Nstep)
                kx  = np.concatenate((kx1, kx2[1:]))
            else:
                kx  = np.linspace(0,k_in,2*Nstep)
        else:
            raise(ValueError,'unknown option for method')
        return kx[:-1]
        
        
        
        
        
        
        


    def CLF(self,omega,i_sys = (0,1),i_in_wave = (5,5),i_out_wave = (5,3), N_step = 100,method = 'diffuse',Signal = True):
        """
        coupling loss factor for line junctions
        
        This method calculates the CLF from system i_sys[0] to i_sys[1]. 
        The wavefields that are considered  are given by i_in_wave and i_out_wave.
        The CLF is calculated running through through both indexes
        

        Parameters
        ----------
        omega : ndarray or float
            angular frequency.
        i_sys : ndarray,list or tuple, optional
            indexes of incident and radiating system. The default is (0,1).
        i_in_wave : ndarray,list or tuple, optional
            DESCRIPTION. The default is (5,5).
        i_out_wave : ndarray,list or tuple, optional
            DESCRIPTION. The default is (5,3).
        N_step : int, optional
            Number of angle integration intervals. The default is 100.
        method : str, 'diffuse' or 'langley', optional
            identifier for calculation method. The default is 'diffuse'.
        Signal : bool, optional
            Switch for Signal output. The default is True.

        Raises
        ------
        ValueError
            For non consistent arguments.

        Returns
        -------
        ndarray or Signal
            coupling loss factor.

        """
                      
        # convert to np.array to become independent from tupels
        # i_sys      = np.array(i_sys)
        # i_in_wave  = np.array(i_in_wave)
        # i_out_wave = np.array(i_out_wave)
        i_sys,i_in_wave,i_out_wave = all2array(i_sys,i_in_wave,i_out_wave)
                
        if i_out_wave.size != i_in_wave.size:
            raise ValueError('In and out wave must have same size')
   
        Nsig = len(i_in_wave)
         
        _tdof = dof.DOFtype(typestr='general')
        etas  = np.zeros((Nsig,omega.size))

        # 1D eta-tau proportionality factor 
        fak1 = self.length/(2*np.pi**2 *omega) # /n(omega) missing 

        xdata = mC.DataAxis(omega,typestr='frequency')
        
        for ifreq,om in enumerate(omega):
            # specify wavenumber for integration
            modal_dens = self.systems[i_sys[0]].modal_density(om,i_in_wave[0])

            # The wavenumber samplingh is tricky because the radiation wavenumber is much 
            # smaller for in-plane waves. So we need special treatment of both regimes1

            #k_in = self.systems[i_sys[0]].plate_wavenumber(om,i_in_wave[0])
            
            kx   = self.kx(om,i_sys,i_in_wave[0],i_out_wave,N_step,method = 'in-plane')
            
            if method == 'diffuse':
                taus = self.transmission_wavenumber(om,kx,i_sys,i_in_wave,i_out_wave,Signal=False) # check late with Signal option
            elif method == 'langley':
                taus = self.transmission_wavenumber_langley(om,kx,i_sys,i_in_wave,i_out_wave,Signal=False)
                        
            for i_o_wave,i_type in enumerate(i_out_wave):
                etas[i_o_wave,ifreq] = fak1[ifreq]/modal_dens*integrate.trapz(taus[i_o_wave,:],kx)
                
                

        if Signal:
            xdata = mC.DataAxis(omega,typestr = 'angular frequency')
            return mC.Signal(xdata,etas,dof.DOF(i_out_wave,np.zeros((1,Nsig)),_tdof))
        else:
            return etas
    
    def CLF_angle(self,omega,i_sys = (0,1),i_in_wave = (5,3),i_out_wave = (5,3), N_step = 200,method = 'diffuse',Signal = True):
        """
        coupling loss factor for wave fields of plates
        
        This method uses the transmission_wave_number_diffuse function that is based in 
        angular integration. The normal CLF method uses the wavenumber integration.
        
        Options for i_in/out_wave are 
        0 : all  
        1 : longitudinal
        2 : shear
        3/4 : bending
        5 : in-plane (longitudinal + shear)
            

        Parameters
        ----------
        omega : np.array
            angular frequency
        L : double
            length of junction
        i_sys : tupel of integer
            pair of in- and out systems. The default is (0,1), if (0,0) only tau between waves is considered
        i_in_wave : TYPE, optional
            wave_fields index of exciting waves. The default is (5,3).
        i_out_wave : TYPE, optional
            wave_fields index of outgoing waves. The default is (5,3).
        N_step : TYPE, optional
            Number of steps for numerical angle integration. The default is 200.
        method : TYPE, optional
            Method of tau calculation. 'diffuse' for diffuse reciprocity, 
            'langley' for wave transmission. The default is 'diffuse'.

        Returns
        -------
        Signals with coupling loss factors

        """
      
        # convert to np.array to become independent from tupels
        i_sys,i_in_wave,i_out_wave = all2array(i_sys,i_in_wave,i_out_wave)
       
        # tau option only for waves with specific wacenumber
        taus = self.transmission_wavenumber_diffuse(omega,i_sys,i_in_wave,i_out_wave,N_step,method)
            
        
        # create data buffer
        if i_out_wave.size != i_in_wave.size:
            raise ValueError('In and out wave must have same size')
   
        Nsig = len(i_in_wave)
        
        etas = np.zeros(( Nsig,np.size(omega) ))
        _tdof = dof.DOFtype(typestr='general')

        
        i_sig = 0

        for i_in,in_wave in enumerate(i_in_wave):
            
            modal_dens = self.systems[i_sys[0]].modal_density(omega,in_wave)
            wavenumber = self.systems[i_sys[0]].plate_wavenumber(omega,in_wave)
            
            fak = wavenumber*self.length/(2*np.pi**2*modal_dens*omega)
            
            for i_out in i_out_wave:
                
                etas[i_sig,:] = fak*taus.ydata[i_sig,:]
                i_sig += 1
                
                
                
            
        if Signal:
            xdata = mC.DataAxis(omega,typestr = 'angular frequency')
            return mC.Signal(xdata,etas,dof.DOF(i_out_wave,np.zeros((1,Nsig)),_tdof))
        else:
            return etas
        
    def junction_matrix(self,omega,N_step = 90, method = 'diffuse'):
        """
        Creates junction matrix for LineJunction (all upper triangular)
        
        This method calculates the all CLF of junction matric [J] 
        The wavefields that are considered ar either the pressure wave or the bending wave
        The CLF is calculated running through through both indexes

        Parameters
        ----------
        omega : float
            angular frequency.
        N_step : int, optional
            Number of integration steps. The default is 90.
        method : str, optional
            CLF calculation method. The default is 'diffuse'.

        Returns
        -------
        JM : DynamiMatrix
            JunctionMatrix.

        """
                
        mod_dens = self.modal_density(omega)

        # prepare relevant wave_DOFs of junction only 3,4 couples to the fluid!          
        j_wave_DOF = self.wave_DOF
        Nw         = self.N_wave

    
        # .. indexes into upper triangular
        i_row,i_col = np.triu_indices(self.N_wave,1)    
 
        # initialise
        JM  = super().junction_matrix(omega)
        fak = self.length/2/np.pi**2/omega
        
        
        #buffer for h_buf
        h_buf = np.zeros((Nw-1,len(omega)))
        

        
        # loop over upper triangular
        # For plates it is better to use the lower triangulat coefficient, because in this
        # case the exciting wave field remains constant and the total stiffness matrix can be reused much better
        # ir_start = 0
           
        for ic in range(Nw-1): # loop over coloums

           # irows = i_row[ir_start:ir_start+self.N-ixc]  # take all row indexes of this coloums
           irows = np.arange(ic+1,Nw)         # take all row indexes of this coloums
           sys_ = self.wave_sys[ic]
           i_in_sys = self.systems.index(sys_)
           Nr = Nw-1-ic # number of rows
           
           # transmission call must be seperated into case of two same systems and
           # different systmes
           
           # get list of respose systems
           sys_IDs = self.wave_DOF[irows].ID
           # prepare input arguments of transmission_wavenumber for each system ID configuraiion
           out_IDs,ix_out,Ns_out= np.unique(sys_IDs, return_index = True,return_counts=True)
           
           #i0 = 0       
           # loop over system configurations
           for ii,ix_sys in enumerate(ix_out):
               Nr = Ns_out[ii]
               i0 = ix_out[ii] # the sort option requires this i0 may be jumping
               irows_ = np.arange(i0,(i0+Nr)) # row index with unique system combination
               i_out_sys=self.systems.index(self.wave_sys[ix_sys])
               h_buf[irows_,:] = -fak*self.transmission_wavenumber_diffuse(omega, (i_in_sys,i_out_sys), \
                                                                      np.tile(self.wave_DOF[ic].dof,(Nr,)), \
                                                                      self.wave_DOF[irows[irows_]].dof,Signal=False)
               #i0 += Nr
           
           # And fill the lower tri and diagonal entries
           for ix,ir in enumerate(irows):
               JM[ir,ic,:] = h_buf[ix,:]    
               JM[ir,ir,:] = JM[ir,ir,:].data - JM[ir,ic,:].data
               JM[ic,ic,:] = JM[ic,ic,:].data - JM[ir,ic,:].data
        
            
        # diagonal!
        return JM
        
class AreaJunction(Junction):
    """ 
    AreaJunction defines the dynamics and transmission of area junctions 

    Attributes
    ----------
    systems : list or tuple
        connected systems
    area : float
        area of AreaJunction
    double_cavity : bool
        switch for simple cavity-cavity connection
    cavity1 : Acoustic3DSystem
        first connected cavity
    cavity2 : Acoustic3DSystem
        second connected cavity
    plate : Structure2DSystem
        plate subsysten
    ix_plate: int
        index of plate subsystem in list
        
    """
 
       
    
    def __init__(self,systems,area=0):
        """
        Class contructor for AreraJunction
        
        Possible system configurations are:
            cavity - plate - cavity
            cavity - cavity
            plate - cavity
            cavity - plate
            
    
        Parameters
        ----------
        systems : list
            involved systems
        area : float
            area of area junction
                
        Examples
        --------
                        
        """
        
        super().__init__(systems)
        
        # Set switch for directly connected => tau = 1
        self.double_cavity = False
                
        if self.N == 3:
            if systems[1].isplate():
                if systems[0].iscavity() and systems[2].iscavity():
                    # create non resonent Transfer Matrix Model
                    self.non_res = systems[1].non_resonant_TMM()
                    self.cavity1 = systems[0]
                    self.cavity2 = systems[2]
                    self.plate   = systems[1]
                    self.ix_plate = 1
                else:
                    raise ValueError('The outer systems in AreaJunction of 3 systems must be a cavity system')
            else:
                raise ValueError('The center system in AreaJunctions of 3 systems must be a structure 2D system')
        elif self.N == 2: # all combinations are possible
            # in any case no non_res system
            self.non_res = mds.TMmodel((iL.MassLayer(1, 0.), ))
            if systems[0].iscavity:
                self.cavity1 = systems[0]
                if systems[1].iscavity():
                    self.cavity2 = systems[1]
                    self.ix_plate = 1E20 # rediculus value
                    self.double_cavity = True
                    
                else:
                    self.plate = systems[1]
                    self.ix_plate = 1
            else:
                self.plate = systems[0]
                self.ix_plate = 0
                if systems[1].iscavity():
                    self.cavity1 = systems[1]
 
                else:
                    raise ValueError('Area junctions must not have 2 plates')

        if area == 0:
            if self.ix_plate < 1E19: # with plate
                self.area = self.plate.area
            else:
                raise ValueError('Without plate in the junction, the area must be given')
        else:
            self.area = area                
                
    def __str__(self):
        """
        Implements __str__ 

        Returns
        -------
        text  : str
             Representing the line junction.

        """
        
        _str = 'AreaJunction with systems:\n'
        for i_,sys in enumerate(self.systems):
            _str += '{0} \n'.format(sys)
            
        _str += 'area       : {0}\n'.format(self.area) 
        
        return _str
        
    def __repr__(self):
        """
        Implements __repr__ 

        Returns
        -------
        text  : str
             Representing the junction.

        """
        return self.__str__()

    def get_wave_DOF(self,ix = slice(None,None)):
        """
        Provides the local wave DOFs of the area junction
        
        Parameters
        ----------
        ix : int, optional
            DESCRIPTION. The default is slice(None,None).

        Returns
        -------
        DOF
            ldofs of junction
            systems participating to junction
            DESCRIPTION.

        """
        
        
        if isinstance(ix,slice):
            N_sea = len(self.systems[ix]) 
            syslist = self.systems[ix]
        elif uf.isscalar(ix):
            N_sea = 1
            syslist = (self.systems[ix],)
            
        # buffers of ID and DOF collection
        _ID  = np.zeros(N_sea,dtype=np.uint32)
        _dof = np.zeros(N_sea,dtype=np.uint32)
        _type = [0]*N_sea
        
        # loop over systems
        for i_sys,syst in enumerate(syslist):
            if syst.iscavity():
                _DOF           = syst.wave_DOF[0]
                _ID[i_sys]     = _DOF.ID[0]
                _dof[i_sys]    = _DOF.dof[0]
                _type[i_sys]   = _DOF.type[0]
            else: 
                _DOF           = syst.wave_DOF[0]
                _ID[i_sys]     = _DOF.ID[0]
                _dof[i_sys]    = _DOF.dof[0] # bending or all must be first
                _type[i_sys]   = _DOF.type[0]
                
        
        return dof.DOF(_ID,_dof,_type),self.systems
                
    def CLF_fluid_fluid(self,omega,pos_dir=True,theta_max = np.pi/180*78):
        """
        Calculates the coupling loss factor of fluid fluid junctions
        
        The infinite layer method is used here for the mass law and non-resonant coupling
        
        Parameters
        ----------
        omega : TYPE
            angular frequency.
        pos_dir : bool, optional
            switch for direction of wave transfer. The default is True.
        theta_max : float, optional
            Maximum angle for diffuse field integration. The default is np.pi/180*78.

        Returns
        -------
        eta_non_res : ndarray
            CLF.

        """
        
        
        if pos_dir:
            cavity1 = self.cavity1
            cavity2 = self.cavity2
        else:
            cavity2 = self.cavity1
            cavity1 = self.cavity2
                
        
        # Calculate mass law using TMM theory
        if self.double_cavity:
            if cavity1.fluid == cavity2.fluid:
                tau_diff = 1 # @todo implement different fluids
            else:
                tau_diff = self.non_res.transmission_diffuse(omega,theta_max,\
                            fluids = (cavity1.fluid,cavity2.fluid),signal = False)
                
        else:
            tau_diff = self.non_res.transmission_diffuse(omega,theta_max,\
                            fluids = (cavity1.fluid,cavity2.fluid),signal = False)
        
        #Calculate eta2D
        k2 = np.real(cavity1.fluid.wavenumber(omega))**2
        #cc = cavity1.fluid.c_freq(omega)

        mod_dens = cavity1.modal_density(omega)
    
        eta_non_res = k2*self.area*tau_diff/(8*np.pi**2*mod_dens*omega)        
        #eta_non_res = cc*self.area*tau_diff/(4*omega*cavity1.volume)
        
        return eta_non_res

    def CLF_structure_fluid(self,omega,ix_cavity = 0,pos_dir=True):
        """
        Calculates the coupling loss factors of plate to fluid junctions

        ----------
        omega : float
            angular frequency.
        ix_cavity : int, optional
            index of cavity. The default is 0.
        pos_dir : bool, optional
            swith for direction. True means from plate to cavity. The default is True.

        Returns
        -------
        ndarray
            coupling loss factor.

        """
        
        
        if ix_cavity < 1:
            cavity = self.cavity1
        else: # 1 and 2 means secont cavity is taken 
            cavity = self.cavity2 # to b echanges later
        
        # check which side of plate must be considered for trim
        if ix_cavity < self.ix_plate: 
            ix_trim = 0
        else:
            ix_trim = 1
            
        z = np.real(cavity.fluid.impedance(omega))
        sigma = self.plate.radiation_efficiency(omega,cavity.fluid,Nstep = 90)
         
        eta_s_fluid = z/omega/self.plate.prop.mass_per_area*sigma
        
        if pos_dir:
            if self.plate.trim_sw[ix_trim]:
                kB = self.plate.prop.wavenumber_B(omega)
                #IC = self.plate.trim[1].insertion_loss(omega,kB,fluid = cavity.fluid)    
                IC = self.plate.trim[ix_trim].insertion_loss_diffuse(omega,self.plate,fluid = cavity.fluid)    
                return eta_s_fluid*IC
            else:
                return eta_s_fluid         
        else:
            modal_dens_cavity = cavity.modal_density(omega)
            modal_dens_plate  = self.plate.modal_density(omega)
            if self.plate.trim_sw[ix_trim]:
                IC = self.plate.trim[ix_trim].insertion_loss_diffuse(omega,self.plate,fluid = cavity.fluid)    
                return eta_s_fluid*IC*modal_dens_plate/modal_dens_cavity
            else:
                return eta_s_fluid*modal_dens_plate/modal_dens_cavity
            
        
            
                   
        
    def CLF(self,omega,i_sys = (0,1), i_in_wave = (0,0),i_out_wave = (0,0),N_step = 90, method = 'diffuse',Signal = True):
        """
        coupling loss factor for area junctions
        
        This method calculates the CLF from system i_sys[0] to i_sys[1]
        and calls the appropriate CLF_fluid_fluid or CLF_structure_fluid 
        method
        
        Parameters
        ----------
        omega : float
            angular frequency.
        i_sys : tuple, list of ndarray, optional
            2x1 vector of input output index. The default is (0,1).
        i_in_wave : int, optional
            index of in wave. The default is (0,0).
        i_out_wave : int, optional
            index of out wave. The default is (0,0).
        N_step : int, optional
            Integration steps for diffuse field integration. The default is 90.
        method : str, optional
            method str for selection of calculation method. The default is 'diffuse'.
        Signal : bool, optional
            switch for Signal output. The default is True.

        Returns
        -------
        Signal or ndarray
            CLF.

        """
                
        # convert to np.array to become independent from tupels
        i_sys, = all2array(i_sys)
           
         
        _tdof = dof.DOFtype(typestr='general')
        etas  = np.zeros((1,omega.size))

        # 2D eta-tau proportionality faltor 
        #fak1 = self.area/(2*np.pi**2 *omega) # /n(omega) missing 

        
        if i_sys[0] < i_sys[1]:
            pos_dir = True
        else:
            pos_dir = False
            
        # fluid to fluid
        if self.systems[i_sys[0]].iscavity():
            if self.systems[i_sys[1]].iscavity(): 
                etas = self.CLF_fluid_fluid(omega, pos_dir = pos_dir)
            else: # fluid to structure (in case of trim trim1)
                etas = self.CLF_structure_fluid(omega,pos_dir = False)
        else: # structure to fluid (in case of trim trim2)
            etas = self.CLF_structure_fluid(omega,pos_dir = True)
                    

        if Signal:
            xdata = mC.DataAxis(omega,typestr = 'angular frequency')
            return mC.Signal(xdata,etas,dof.DOF(0,np.zeros((1,1)),_tdof))
        else:
            return etas     

    def junction_matrix(self,omega,N_step = 90, method = 'diffuse',Signal = True):
        """
        Calculates the local SEA matrix for the junction DOFs
        
        This method calculates the all CLF of junction matric [J] 
        The wavefields that are considered ar either the pressure wave or the bending wave
        The CLF is calculated running through through both indexes
        

        Parameters
        ----------
        omega : float
            angular frequency.
        N_step : int, optional
            Integration steps for diffuse field integration. The default is 90.
        method : str, optional
            method str for selection of calculation method. The default is 'diffuse'.
        Signal : bool, optional
            switch for Signal output. The default is True.

        Returns
        -------
        JM : DynamicMatrix
            Junction SEA matrix.

        """
                
        mod_dens = self.modal_density(omega)

        # prepare relevant wave_DOFs of junction only 3,4 couples to the fluid!          
        j_wave_DOF = self.wave_DOF

    
        # .. indexes into upper triangular
        i_row,i_col = np.triu_indices(self.N_wave,1)    
 
        # initialise using mother class method
        JM = super().junction_matrix(omega)

        
        # loop over upper triangular
        for ix,ic in enumerate(i_col):
           ir = i_row[ix]
           if self.wave_sys[ic].iscavity():
               if self.wave_sys[ir].iscavity():
                   JM[ir,ic,:] = -mod_dens[ic]*self.CLF_fluid_fluid(omega, pos_dir = ic < ir)
               else: # fluid to structure (in case of trim trim1)
                   JM[ir,ic,:] = -mod_dens[ic]*self.CLF_structure_fluid(omega,ix_cavity = ic, pos_dir = False)
           else: # structure to fluid (in case of trim trim2)
                JM[ir,ic,:] = -mod_dens[ic]*self.CLF_structure_fluid(omega,ix_cavity = ir, pos_dir = True)
           # And fill the diagonal entries
           JM[ir,ir,:] = JM[ir,ir,:].data - JM[ir,ic,:].data
           JM[ic,ic,:] = JM[ic,ic,:].data - JM[ir,ic,:].data
            
            
        # diagonal!
        return JM










class HybridAreaJunction(Junction) :
    """ 
    HybridAreaJunction defines the physics and systems of hybrid area junctions
    for SEA-FEM or SEA-FEM-SEA configuration

    Attributes
    ----------
    systems: list of SEA systems
    FEMsys:  list of FEMsystems
    """
 
       
            
    def __init__(self,systems,fem,coupling='all',trim = None):
        """
        Class contructor for HbridAreaJunction
        

        Parameters
        ----------
        systems : tuple, list or ndarray
            vector of SEA systems.
        fem : fem
            FEM system.
        coupling : str,tuple, list or ndarray, optional
            indexes defining nodes coupled to SEA subsystem. The default is 'all'.
        trim : tuple, list of TMmodel, optional
            trim layer definition. The default is None.

        Returns
        -------
        None.

        """
        
        super().__init__(systems)
        # mesh and normal modes can be later combined to FEMresult object
        self.fem = fem
        self.coupling = []
        self.trim = []
        
        self.trim_sw = False # Global switch it trim is used
        self.trim_list_sw = [False]*self.N

        
        # Deal with trim
        if trim != None:
            if isinstance(trim,(tuple,list)) and len(trim)==self.N:
                
                for ix,_trim in enumerate(trim):
                    self.trim.append(_trim)
                    if isinstance(_trim,mds.TMmodel):
                        self.trim_list_sw[ix] = True
                        self.trim_sw = True
                    
        
        if isinstance(coupling,str) and coupling == 'all':
            for ix in enumerate(systems):
                self.coupling.append(fem.mesh)
                
    def CLF(self,omega,trim = None, force = None):
        """
        coupling loss factor for hybrid area junction
        
        This method calculates the CLF from system i_sys[0] to i_sys[1]. 
        The wavefields that are considered are either the pressure waves
        The CLF is calculated running through through both indexes
        
        This method provides the total stiffness, thus it is efficient to perform all solutions
        that require the total stiffness matrix  

        Parameters
        ----------
        omega : ndarray
            angular frequency.
        trim : TMmodel, optional
            trim as layers. The default is None.
        force : Load, optional
            force defined on nodes and DOFs. The default is None.
        Signal : bool, optional
            switch for Signal output. The default is True.

        Returns
        -------
        CLF : ndarray
            hybrid coupling loss factor         
        CLF_alpha : ndarray
            dymping of SEA systems due to loss in FE-system
        input_power : ndarray
            power input to the SEA subsystem due to the force load
        modal_displacement : ndarray
            modal displacement due to force load
        """
        
        
        
        
        # modal density of all connected systems
        damping_loss = self.fem.damping_loss(omega) 
               
        mod_dens = self.modal_density(omega)
        
        # initialize practical factors    
        alpha_fac  = 2/np.pi/omega
        eta_fac    = alpha_fac
        femload_sw = False
        
        # .. indexes into upper triangular
        i_row,i_col = np.triu_indices(self.N_wave,1)
        
           
        # prepare soluion arrays (type double)    0
        #_tdof = dof.DOFtype(typestr='general')
        etas       = np.zeros((self.N_wave*(self.N_wave-1)//2,omega.size))
        etas_alpha = np.zeros((self.N_wave,omega.size))
        power_in   = np.zeros((self.N_wave,omega.size))
        modal_disp = np.zeros((self.fem.modes.Nxdata,omega.size),dtype = np.complex128)
        
        # preallocate lists
        Ddir=self.N_wave*[None]
        DDD =self.N_wave*[None]
        
        if self.trim_sw:
            # create buffer for trim matrices
            D11 = self.N_wave*[None]
            D12 = self.N_wave*[None]
            D22 = self.N_wave*[None]
            Dred = self.N_wave*[None]
            D_F2_F1= self.N_wave*[None]
            
                    

        # collect radation stiffness matrices
        # Start with first mesh
        hs = aR.HalfSpace(self.systems[0].fluid)
        HS = [hs] # 1.radiation_stiffness_mesh(omega, self.mesh,method = 'wavelet')]

        if self.N_wave > 1:
            D_equal = np.zeros((self.N_wave-1,),dtype=bool)
            # @todo currently only two halfspaces are checked if equal
            for ix in range(1,self.N_wave):
                if self.coupling[0] == self.coupling[ix]:
                    HS.append(hs) # use the same
                    D_equal[ix-1] = True
                else:
                    HS.append(aR.HalfSpace(self.systems[ix].fluid))

        # modal base in short form
        MM          = self.fem.modes.ydata
        
        if force in self.fem.loads.keys():
            modal_force = self.fem.modal_force(force)
            femload_sw    = True
            
        print('Solving hybrid Junction...')
        
        
        
        # frequency loop
        for i_om,om in enumerate(omega):
            print("\r","Dealing with omega = {0:>8.1f} Hz\r".format(om), end="")

            # onvert radiation stiffness to modal space
            for idir,hs_ in enumerate(HS):
                if idir < 1:
                    D_ = HS[idir].radiation_stiffness_mesh([om], self.fem.mesh,method = 'wavelet').Dindex(0) # only one freq calc.
                    Ddir[idir] = np.conj(MM.T).dot(D_).dot(MM)
                    # Check for trim and calculate Reduced stiffness
                        
                else:
                    if D_equal[idir-1]:
                        Ddir[idir] = Ddir[0]
                    else:
                        D_ = HS[idir].radiation_stiffness_mesh([om], self.fem.mesh,method = 'wavelet').Dindex(0)
                        Ddir[idir] = np.conj(MM.T).dot(D_).dot(MM)

            # deal with trim
            if self.trim_sw:
                for idir,hs_ in enumerate(HS):
                    # Check for trim and calculate Reduced stiffness
                    if self.trim_list_sw[idir]:
                        D11_,D12_,D22_ = self.trim[idir].stiffness_matrix_mesh([om],self.fem.mesh)
                        D11[idir] = np.conj(MM.T).dot(D11_.Dindex(0)).dot(MM)
                        D12[idir] = np.conj(MM.T).dot(D12_.Dindex(0)).dot(MM)
                        D22[idir] = np.conj(MM.T).dot(D22_.Dindex(0)).dot(MM)
                        # 
                        D_F2_F1[idir] = D12[idir].dot(np.linalg.inv(D22[idir]+Ddir[idir]))
                        Dred[idir] = D_F2_F1[idir] @ Ddir[idir] @ D_F2_F1[idir].conj().T 


                        
            # from now D1_ and D2_ are the inverse
            Ds_mod    = np.diag(self.fem.modal_dynamic_stiffness_diag(om))
            D_tot_inv = Ds_mod
            for ix in range(self.N_wave):
                if self.trim_list_sw[ix]:
                    # Add reduced form or D_dir
                    #D_tot_inv -= D_F2_F1[idir]
                    D_tot_inv += D11[idir] - D_F2_F1[idir] @ D12[idir]
                else:
                    # use the normal one
                    D_tot_inv += Ddir[ix] # add the di
                    
                
            D_tot_inv = np.linalg.inv(D_tot_inv)
            
            if femload_sw:
                modal_disp [:,i_om] = D_tot_inv.dot(modal_force.ydata[:,i_om])
                for ix,_ in enumerate(self.wave_DOF):
                    if self.trim_list_sw[ix]:
                        power_in[ix,i_om] = om/2*np.imag(modal_disp [:,i_om].conj().T @ Dred[ix] @ modal_disp [:,i_om])
                    else:
                        power_in[ix,i_om] = om/2*np.imag(modal_disp [:,i_om].conj().T @ Ddir[ix] @ modal_disp [:,i_om])
            
            # get all D^-1 Im D_in D^-H for each D_in             
            for ix in range(self.N_wave):
                if self.trim_list_sw[ix]:
                    DDD[ix] = np.real(D_tot_inv.dot(np.imag(Dred[ix])).dot(np.conj(D_tot_inv.T)))
                else:
                    DDD[ix] = np.real(D_tot_inv.dot(np.imag(Ddir[ix])).dot(np.conj(D_tot_inv.T)))

            
            for ix,ir in enumerate(i_row):
                # coupling loss factors
                if self.trim_list_sw[ix]:
                    etas[ix,i_om] = (np.imag(Dred[i_col[ix]])*DDD[ir]).sum()
                else:
                    etas[ix,i_om] = (np.imag(Ddir[i_col[ix]])*DDD[ir]).sum()
            
            for ix,_ in enumerate(self.wave_DOF):
                etas_alpha[ix,i_om] = damping_loss[i_om]/mod_dens[ix,i_om]*\
                    (np.imag(Ds_mod)*DDD[ix]).sum()
                    
                    
             
        # correct for upper triangular
        for ix,ir in enumerate(i_row):
            etas[ix,:] *= (1/mod_dens[i_row,:]).flatten()

        print('\n...finished\n')
        if femload_sw:
            return etas*eta_fac,etas_alpha*alpha_fac,power_in,modal_disp
        else:
            return etas*eta_fac,etas_alpha*alpha_fac
            

    def FEM_response(self,omega,energy):
        """        
        Calculates the cross spectral resonse of the plate due to 
        reverberant fields in SEA systems

        Parameters
        ----------
        omega : ndarray
            angular frequency.
        energy : Signal
            energy in SEA systems.

        Returns
        -------
        ndarray : 
            CSD of FEM DOFs.

        """
      
        # modal density of all connected systems
        mod_dens = self.modal_density(omega)
        
        # initialize practical factors    
        fak  = 4/np.pi
        
        # .. indexes into upper triangular
        i_row,i_col = np.triu_indices(self.N_wave,1)
        
           
        # prepare soluion arrays (type double)    0
        #_tdof = dof.DOFtype(typestr='general')
        Sqq = np.zeros((self.fem.modes.Nxdata,self.fem.modes.Nxdata,omega.size),dtype = np.complex128)
        
        # preallocate lists
        Ddir=self.N_wave*[None]
                    

        # collect radation stiffness matrices
        # Start with first mesh
        hs = aR.HalfSpace(self.systems[0].fluid)
        HS = [hs] # 1.radiation_stiffness_mesh(omega, self.mesh,method = 'wavelet')]

        if self.N_wave > 1:
            D_equal = np.zeros((self.N_wave-1,),dtype=bool)
            # @todo currently only two hs are checked if equal
            for ix in range(1,self.N_wave):
                if self.coupling[0] == self.coupling[ix]:
                    HS.append(hs) # use the same
                    D_equal[ix-1] = True
                else:
                    HS.append(aR.HalfSpace(self.systems[ix].fluid))

        # modal base in short form
        MM          = self.fem.modes.ydata
        

        # frequency loop
        print('Calulating FEM response of hybrid junctions...')
        for i_om,om in enumerate(omega):
            print("\r","Dealing with omega = {0:>8.1f} Hz".format(om), end="")

            #convert to modal space
            for idir,hs_ in enumerate(HS):
                if idir < 1:
                    D_ = HS[idir].radiation_stiffness_mesh([om], self.fem.mesh,method = 'wavelet').Dindex(0) # only one freq calc.
                    Ddir[idir] = np.conj(MM.T).dot(D_).dot(MM)
                else:
                    if D_equal[idir-1]:
                        Ddir[idir] = Ddir[0]
                    else:
                        D_ = HS[idir].radiation_stiffness_mesh([om], self.fem.mesh,method = 'wavelet').Dindex(0)
                        Ddir[idir] = np.conj(MM.T).dot(D_).dot(MM)
                        
            # from now D1_ and D2_ are the inverse
            Ds_mod    = np.diag(self.fem.modal_dynamic_stiffness_diag(om))
            D_tot_inv = Ds_mod+Ddir[0]
            for ix in range(1,self.N_wave):
                D_tot_inv += Ddir[ix]
                
            D_tot_inv = np.linalg.inv(D_tot_inv)
            
            # get all D^-1 Im D_in D^-H for each D_in
            
            for ix,dof_ in enumerate(self.wave_DOF):
                ie = energy.dof.find(dof_.ID,dof_.dof)
                E  = energy.ydata[ie,i_om]
                Sqq[:,:,i_om] += fak/om/mod_dens[ie,i_om]*E*np.real(D_tot_inv.dot(np.imag(Ddir[int(ie)])).dot(np.conj(D_tot_inv.T)))


        print('\nfinished')

        return Sqq
            
        
        
class SemiInfiniteFluid(AreaJunction):
    """
    Class for non reverberant fluid sinks
    
    The aim of this class is to add absoption to the connected area junctions, to calculate
    the radiated power into the sink and to calculate the power at certain discance
    """
    

    def __init__(self,systems, fluid, area = 0.):
        """
        Constructor of simi infinite fluids

        Parameters
        ----------
        systems : tuple of list of SEA systems
            connected SEA systems.
        fluid : fluid
            DESCRIPTION.
        area : float, optional
            area of SIF. The default is 0..

        Returns
        -------
        None.

        """
        

        super().__init__(systems,area)
        
        self.fluid = fluid
        # create stange cavity to allow for non resonent methods of mother class
        self.cavity2 = ac3Dsys.Acoustic3DSystem(0,0,0,0,fluid)
        # in case of 2 connected systems the non_res attribute is not set correctly in the 
        # mother class constructor 
        if self.N == 2 and not(self.double_cavity):
            self.non_res = self.plate.non_resonant_TMM()

              
    def __repr__(self):
        _str = 'SEA semi infinite fluid {0}\n connected to :'.format(self.fluid)
        for sys in self.systems:
            _str += '{0}'.format(sys) 
        return _str
    
    def resonant_dampingloss(self,omega):
        """
        resonant DLF of connected plates

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        eta : ndarray
            damping loss factor.

        """
        
        eta = self.CLF_structure_fluid(omega,ix_cavity = 2,pos_dir=True)
        
        return eta
        
    def non_resonant_dampingloss(self,omega):
        """
        non-resonant DLF of connected cavities

        Parameters
        ----------
        omega : float
            angular frequency.

        Returns
        -------
        eta : ndarray
            damping loss factor.

        """
        eta = self.CLF_fluid_fluid(omega)
        return eta

    