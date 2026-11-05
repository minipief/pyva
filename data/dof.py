# -*- coding: utf-8 -*-
"""
Module for degrees of freedom management

The dof module deals with the handling of nodal and system (wave) degrees of freedom. 
This comprised the identification of node IDs and related local degrees of freedom.
An example would be a node with ID = 99 and displacment in x-, y- and z-direction.

In finite element methods this is called a grid or mesh when combined with the nodal
positions.

In addition every DOF is linked to a physical property, for example 
displacement or pressure. This is treated by the DOFtype class.
The aim is to automatically manage units and units systems, but also to
allow operations between physical quantities.

For example a pressure multiplied by an area leads into a force. 
"""
import numpy as np
import copy
#from pyva import ureg
from pint import UnitRegistry

ureg = UnitRegistry()
Q_   = ureg.Quantity


# useful methods
def shape(lst):
    dims = np.shape(lst)
    if len(dims)>1:
        return dims
    elif len(dims)==1:
        return(dims,1)
    else:
        raise ValueError('Input must not be scalar')
        
#    if isinstance(lst[0],str):
#        dim2 = 1
#    else:
#        dim2 = len(lst[0])
#        if dim1 > 1:
#            for idim1 in range(1,dim1):
#                if dim2 != len(lst[idim1]):
#    return np.shape(lst)




# Class definition and methods


class DOF:
    """
    Class for degrees of freedom
    
    The attributes are organise as vectors, so that a list of DOFs is represented
    by the node ID, the local orientation and physical unit.

    
    Attributes
    ----------
    ID : ndarray of int
        ID for nodes, elements, etc. 
    dof : ndarray of int
        local degree of freedom  e.g. 0 for sclaer 1,2,3 for x,y,z and 4,5,6 for RX,RY,RZ
    doftype : list
        DOFtype 
    
    

    """
    
    
    def __init__(self, ID, dof, doftype, repetition = False):
        """
        Costructor of DOF objects

        Parameters
        ----------
        ID : ndarray, list, tuple of integer or integer   
            ID for nodes, elements, etc. 
        dof : ndarray, list, tuple of integer or integer   (of same dimension as ID)
            degree of freedom e.g. 0 for sclaer 1,2,3 for x,y,z and 4,5,6 for RX,RY,RZ
        doftype : list of DOFtype, string or typeID
            type for DOF, e.g. pressure, displacent delivers unit information
        repetition : bool, optional
            switch for repetition, when True doftype is used as pattern, every dof is repeted for every ID

        Raises
        ------
        ValueError
            For incorrect input values and value combinations.

        Returns
        -------
        DOF object

        """

        
        ID = np.array(ID,dtype = np.dtype('u4')).flatten()
        self._ID = ID
        dof = np.array(dof,dtype = np.dtype('u4')).flatten()
        self._dof = dof 
        self._uniquetype_sw = False

        if isinstance(doftype,DOFtype):
            self._type = [doftype] # doftypes must be lists to include dimension information, how ever single doftypes are alloe
        elif isinstance(doftype,(list,tuple)):
            # Check if all members are of type DOFtype
            # and put values into list
            self._type = []
            for ilist in range(len(doftype)):
                if isinstance(doftype[ilist],DOFtype):
                    self._type.append(doftype[ilist])  
                elif isinstance(doftype[ilist],str):
                    self._type.append(DOFtype(typestr=doftype[ilist]))
                else:                    
                    raise ValueError('All doftype members must be of class DOFtype or str')
                    
                    
        else:
            raise ValueError('dof argument must by of class DOFtype')

        self._repetition = repetition
        
        if repetition:
            self._size = (1,len(ID)*len(dof))
            if len(doftype)==1: # single type for all
                #self._type = doftype # should be list already
                self._uniquetype_sw = True
            else:
                raise ValueError('If repetition is true, the doftype mulst be single')
            
        else:
            # check consistent length
            if self._ID.shape == self._dof.shape:
                #self._ID = np.array(ID,dtype = np.dtype('u4'))
                #self._dof = np.array(dof,dtype = np.dtype('u4'))
                if len(doftype)==1 and self._ID.size > 1: # single type for all and not single ID!!!
                    #self._type = [doftype]
                    self._uniquetype_sw = True
                elif shape(self._type) == self._ID.shape:
                    self._uniquetype_sw = False
                    # loop over all fields and allow str, ID and DOFtype
                    Nrow = self._ID.shape[0]
                    Ncol = self._ID.shape[1]
                    print(Nrow,Ncol)
                    # Initialies list for types
                    self._type = [[None]*Ncol for _ in range(Nrow)]
                    for irow in range(Nrow):
                        for icol in range(Ncol):
                            if isinstance(doftype[irow][icol],int):
                                self._type[irow][icol] = DOFtype(type=doftype[irow][icol])
                            elif isinstance(doftype[irow][icol],str):
                                self._type[irow][icol] = DOFtype(typestr=doftype[irow][icol])
                            elif isinstance(doftype[irow][icol],DOFtype):
                                self._type[irow][icol] = doftype[irow][icol]
                            else:
                                raise ValueError('type must be of same dimension as ID and dof')
                    
            else:
                raise ValueError('ID and dof must be of same dimension')
                    
            self._size = (1,len(self._ID))
            
           
            
    def __len__(self):
        """
        Oveloaded __len__ method

        Returns
        -------
        int
            length of DOF object.

        """
        
        return self._size[0]*self._size[1]
    
    @property    
    def size(self):
        """
        size of DOF object
        
        TODO size and shape are different things, should be corrected
        
        Returns
        -------
        int
            size.

        """
        return self._size[0]*self._size[1]
    
    @property
    def shape(self):
        """
        size of DOF object
        
        Returns
        -------
        int
            shape of DOF object.

        """
        return self._size
    
    
    def __str__(self):
        """
        overloaded __str__ method

        Returns
        -------
        str
            of DOF.

        """
        return "DOF object with ID {0}, DOF {1} of type {2}".format(self.ID,self.dof,self.type)

    def __repr__(self):
        """
        overloaded __repr__ method

        Should allow class construction, must be improved

        Returns
        -------
        str
            of DOF.

        """
        return str(self)
    
    def __eq__(self,other):
        """
        overloaded __eq__ method

        Parameters
        ----------
        other : DOF

        Returns
        -------
        bool

        """
        return all(self.ID==other.ID) and all(self.dof==other.dof) and (self.type==other.type)

    def __ne__(self,other):
        """
        overloaded __ne__ method

        Parameters
        ----------
        other : DOF

        Returns
        -------
        bool
        """
        return not(self == other)
    
    def __add__(self,other):
        
        if self.isuniform() and other.isuniform() and \
            self.type[0] == other.type[0]:
            return DOF(np.append(self.ID,other.ID),\
                   np.append(self.dof,other.dof), \
                   self._type)
        else:
            return DOF(np.append(self.ID,other.ID),\
                   np.append(self.dof,other.dof), \
                   self.type + other.type)
            
                
                
    
    def __getitem__(self, position):
        """
        Overlaoded __getitem__ method for indexing

        Parameters
        ----------
        position : index type
            index of DOF.

        Returns
        -------
        DOFtype
            of index.

        """
        if self._uniquetype_sw:
            return DOF(self.ID[position].flatten(), \
                       self.dof[position].flatten(),self._type[0])
        else:
            # extra indexing of lists
            if isinstance(position,(int,slice,np.integer)):
                type_buf = self.type[position]
            elif isinstance(position,(np.ndarray,list,tuple)):
                type_buf = []
                for i in range(len(position)):
                    type_buf.append(self.type[position[i]])
                    
            return DOF(self.ID[position].flatten(),
                       self.dof[position].flatten(),
                       type_buf) # self.type[position]) #.flatten()    
    def label(self):
        """
        Creates label
        
        Label can be used for plotting or in log files

        Returns
        -------
        str
            label for all dofs.

        """
        return "ID {0} DOF {1}".format(self.ID,self.dof)
        
    def index(self,other,**kwargs):
        """
        Get index of other into DOF of self
        
        Searches th positions of other in self and provides the index of
        these positions
        
        ix = index(A,B)
        
        A[ix] == B
        
        Parameters
        ----------
        other : DOF
            DOFs that are searched in self.
        **kwargs : dict
            Arbitrary keyword arguments.

        Returns
        -------
        ndarray
            Indexes.

        """
        
        
        iRes = np.zeros(len(other),dtype=np.int64)
                
        # Create tupelof hashes for fast ans easy search     
        sT = self.hashtupel(**kwargs)
        oT = other.hashtupel(**kwargs)
        
        for i in range(len(iRes)):
            try:
#                iRes[i] = list(sT).index(oT[i])
                iRes[i] = sT.index(oT[i])
            except: # no index found
                iRes[i] = -1

                #raise IndexError('No DOF found' )
    
        return iRes #list(iRes)
    
    def find(self,ID,dof =-1,repetition=True):
        """
        index method for finding index from ID and dof, ignoring doftype
        

        Parameters
        ----------
        ID : list, tuple or ndarray
            node identifiers.
        dof : degree of freedom, optional
            DESCRIPTION. The default is -1.
        repetition : bool, optional
            When true every ID and dof is repeated. The default is True.

        Raises
        ------
        ValueError
            For repetition == False ID and dof must have same dimentions.

        Returns
        -------
        ndarray of int
            index of found DOFs.

        """
        
        ID  = np.array(ID).flatten()
        dof = np.array(dof).flatten()
        
        
        if not(isinstance(dof,int)) and np.all(dof > -1):   
            dof = np.array(dof).flatten()
            if not(repetition):
                if len(ID)!=len(dof):
                    raise ValueError('In case of non-repetitive find arguments the length of ID and dof must be equal')
            return self.index(DOF(ID,dof,DOFtype(typestr='unknown'),repetition=repetition),ignore='doftype')
        else:

            return self.index(DOF(ID,[0],DOFtype(typestr='unknown'),repetition=repetition),ignore=('dof','doftype') )
    
    
    def hashtupel(self,**kwargs):
        """
        Returns a list of hash for fast sorting

        Parameters
        ----------
        **kwargs : dict
            Arbitrarty keyword arguments.

        Raises
        ------
        ValueError
            For wrong keyword argument.

        Returns
        -------
        list
            of hashes.

        """

        sID = self.ID;  sdof = copy.copy(self.dof);  sft = self.fulltype
           
        # set values to 0 if component should be ignored
        for kw in kwargs:
            if kw == 'ignore':
                if isinstance(kwargs[kw],str):
                    if kwargs[kw] == 'ID':
                        sID[:]  = 0
                    elif kwargs[kw] == 'dof':
                        sdof[:] = 0
                    elif kwargs[kw] == 'doftype':
                        sft[:]  = 0
                    else:
                        raise ValueError('Unknown type to ignore')
                elif isinstance(kwargs[kw],(list,tuple)):
                    ignore_list = kwargs[kw]
                    for itup in range(len(kwargs[kw])):
                        if ignore_list[itup] == 'ID':
                            sID[:]  = 0
                        elif ignore_list[itup] == 'dof':
                            sdof[:] = 0
                        elif ignore_list[itup] == 'doftype':
                            sft[:]  = 0
                        else:
                            raise ValueError('Unknown type to ignore')
                else:
                    raise ValueError('hashtupel kwargs must be str or tupel of str')
                            
        return [ hash((sID[i],sdof[i],sft[i])) for i in range(len(sID)) ]

    def tupel(self):
        """
        Is this really used

        Returns
        -------
        list
            DESCRIPTION.

        """
        sID = self.ID;  sdof = self.dof;  sft = self.fulltype
     
        return [ (sID[i],sdof[i],sft[i]) for i in range(len(sID)) ]


        
        
    @property
    def ID(self):
        """
        property method for ID

        Returns
        -------
        tuple of IDs
            DESCRIPTION.

        """
        if self._repetition:
            return np.repeat(self._ID,len(self._dof))
        else:
            return self._ID
            

    @property
    def dof(self):
        """
        property method for dof

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        if self._repetition:
            return np.tile(self._dof,len(self._ID))
        else:
            return self._dof

    @property
    def fulltype(self):
        """
        Provides the full list of type IDs for indexing purposes
        
        Extends the list if repetitive type descriptions are used

        Returns
        -------
        np.ndarray
            of typeIDs.

        """
        if self._uniquetype_sw:
            return np.tile(self._type[0].fulltype,len(self))
        else:
            if self._repetition:
                ftb = np.zeros(len(self._dof))
                for ir in range(len(self._dof)):
                    ftb[ir] = self._type[ir].fulltype
                return np.tile(ftb,len(self._ID))
            else:
                ftb = np.zeros(len(self))
                for ir in range(ftb.size):              
                    ftb[ir] = self._type[ir].fulltype
                return ftb
        
        
       
    @property
    def type(self):
        """
        property method for DOF(types)

        Returns
        -------
        list
            of DOFtypes.

        """
        # if self.isuniform():
        #     return self._type*self.size
        # else:
        #     return self._type
    
        return self._type
    
    @property
    def typeID(self):
        
        tID = np.zeros((self.size,),dtype = np.int32 )

        if self._repetition:
            tID[:] = self.type[0].type
        else:
            for ix,dt in enumerate(self.type):
                tID[ix] = dt.type

        return tID
    
    @property
    def unique_type(self):
        """
        Provides unique type if all DOfs are of same type

        Raises
        ------
        ValueError
            When DOFtype is not unique.

        Returns
        -------
        DOFtype
            unique DOFtype.

        """
        if self._uniquetype_sw:
            return self._type[0]
        else:
            _types = np.unique(self.fulltype)
            if len(_types) == 1:
                return self._type[0]
            else:
                raise ValueError('dof is not unique')

            
    def intersect(self,other,**kwargs):
        """
        Intersects DOF and provides indices into subparts
        
        Parameters
        ----------
        other : TYPE
            DESCRIPTION.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.
        idx1 : TYPE
            DESCRIPTION.
        idx2 : TYPE
            DESCRIPTION.

        """
        
        idx1 = self.index(other,**kwargs)
        idx2 = other.index(self,**kwargs)
        
        idx1 = [x for x in idx1 if x>= 0]
        idx2 = [x for x in idx2 if x>= 0]
        
               
        return (self[idx1],idx1,idx2)
    
    def union(self,other):
        """
        Unites DOF and provides indices into subparts


        Parameters
        ----------
        other : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.
        idx1 : TYPE
            DESCRIPTION.
        idx2 : TYPE
            DESCRIPTION.

        """
        
        idx1 = self.index(other)
        idx2 = other.index(self)
        
        idx1 = [x for x in idx1 if x>= 0]
        idx2 = [x for x in idx2 if x>= 0]
        
               
        return (self[idx1],idx1,idx2)
            
    def isuniform(self):
        """
        checks if all DOFs are of the same doftype


        Returns
        -------
        bool
            DESCRIPTION.

        """
        
        if self._uniquetype_sw:
            return True
        else:
            _types = np.unique(self.fulltype)
            if len(_types) == 1:
                return True
            else:
                return False
            
    @property
    def dBref(self):
        """
        Provides reference value for dB values
        
        For example the dBref of pressure is p0 = 20e-6 Pa

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        float
            reference value.

        """
        
        if self.isuniform():
        # @todo consider denominator and nominatr
            return self.type[0].dBref    
        else:
            raise ValueError('No dB reference available for multiDOF Signals')
        

def gettypeID(IDstr):
    """
    Get ID from typestr

    Parameters
    ----------
    IDstr : str
        type string.

    Raises
    ------
    ValueError
        When type string is not existing.

    Returns
    -------
    key : int
        typeID.

    """
    
    if IDstr in list(TYPEDICT.values()):
        for key in TYPEDICT:
            if TYPEDICT[key]==IDstr:
                return key
    else:
        raise ValueError('IDstr {0} not in DOF list \n {1}'.format(IDstr))
                

# Module constants used in DOF class 
        
TYPEDICT = {   0 : 'unknown',          1 : 'general',
               2 : 'stress',           3 : 'strain',
               5 : 'temperature',      6 : 'heat flux',
               8 : 'displacement',    10 : 'volume flow',
              11 : 'velocity',        12 : 'acceleration',
              13 : 'force',           15 : 'pressure',
              16 : 'mass',            17 : 'time',
              18 : 'frequency',       19 : 'rpm',
              20 : 'order',           21 : 'angular frequency', 
              22 : 'sound intensity',
              23 : 'power',           24 : 'area',
              25 : 'transmission' ,
              26 : 'wavenumber'   ,   27 : 'energy' ,
              28 : 'impedance' }

SIUNITDICT = { 0 : '',                1 : '',
               2 : 'pascal',          3 : 'pascal',
               5 : 'kelvin',          6 : '----',
               8 : 'meters',         10 : 'meter**3/second',
              11 : 'meter/second',   12 : 'meter/second**2',
              13 : 'newton',         15 : 'pascal',
              16 : 'kg',             17 : 'seconds',
              18 : 'hertz',          19 : '1/minute',
              20 : '',               21 : '1/second', 
              22 : 'watt/meter**2',
              23 : 'watt',           24 : 'meter**2' ,
              25 : '' ,
              26 : '1/meter'      ,  27 : 'joule' ,
              28 : 'pascal*second/meter'}

# Dictionary for L,M,T and related type
LMTDICT    = {( 0, 0, 0) : 1 ,
              (-1, 1,-2) : 15,
              ( 1, 0, 0) : 8 ,
              ( 3, 0,-1) : 10,
              ( 1, 0,-1) : 11,
              ( 1, 0,-2) : 12,
              ( 1, 1,-2) : 13,
              ( 0, 1, 0) : 16,
              ( 0, 0, 1) : 17,
              ( 0, 0,-1) : 18,
              ( 0, 1,-3) : 22,
              ( 2, 1,-3) : 23,
              ( 2, 0 ,0) : 24,
              (-1, 0, 0) : 26,
              ( 2, 1,-2) : 27,
              (-2, 1,-1) : 28 }

# Dictionary that provides the multiplication of types
MULDICT = { (15,24) : 13 , # pressure * area = force
            (24,15) : 13 , #
            (22,24) : 23 , # intensity * area = power
            (24,22) : 23 ,
            ( 8, 8) : 24 , # length **2
            (17,18) : 1  , # time * frequency
            (18,17) : 1  ,
            (11,24) : 10 ,
            (24,11) : 10}

DIVDICT = { (13,24) : 15 ,
            (23,24) : 22 ,
            (8 ,17) : 11 ,
            (10,24) : 11 }

DBREF   = {    8 : (1.E-12, 2,'pm') ,      
              11 : (1.E-9 , 2,'nm/s')  ,        12 : (1.E-6, 2, '\mu m') ,
              15 : (2.E-5 , 2, '20 \mu Pa')  ,  23 : (1.E-12,1, 'pW'),
              25 : (   1. ,-1 ,'1') } 



def typeIDfromLMT(LMT):
    return LMTDICT.get(tuple(list(LMT)),0)



class DOFtype:
    """
    The DOFtype class provides methods to deal with types of phyical quantities
    in combination with units. 
    
    Units are defined by the pint package, but this does
    not include a clear definition of the phyiscal quantity. For
    example a force can be a reaction force of a force itself.
    Stress and pressure have the same base units so there is a need 
    for managing the quantitiy.
    
    The ``exp`` and ``xexp`` attributes are aimed at spectral data 
    infomation. For example if a specta is an amplitude spectrum,
    a an amplitude per root frequency or the squared amplitud per 
    frequency. Thus, they are not inpluded in the mathematical 
    operations ``+, -, *, / `` 
    
    This class uses the pint module
    
    Attributes
    ----------
        type : int
            ID of DOF type, e.g. 8 for displacement
        xtype : int
            ID of nodal DOF for x axis, e.g. 1,2,3 for x,y,z direction 4 for 
            roation around x-axis
        exp :  float
            Exponent of the data, used for Signal power e.g. 2 for p^2
        xexp : floot
            Exponent of data-axis exponent, e.g. -1 spectral density p^2/f 
        unit : ureg
            Physical unit of the quantity
        
    
    .. _tab-DOFtype:IDs:
        
    .. table:: Idientifier of DOFtype

        ==== =================    
        type physical quantity
        ==== =================    
         0   'unknown',
         1   'general',
         2   'stress',
         3   'strain',
         5   'temperature',
         6   'heat flux',
         8   'displacement',
         9   'force',
        10   'volume flow',
        11   'velocity',
        12   'acceleration',
        13   'excitation force',
        15   'pressure',
        16   'mass',
        17   'time',
        18   'frequency',
        19   'rpm',
        20   'order',
        21   'angular frequency',
        22   'sound intensity',
        23   'power'
        24   'area'
        25   'transmission'
        26   'wavenumber'
        27   'energy'
        ==== =================    

    """    
 
    def __init__(self, **kwargs):
        """
        Constructor of DOFtype class

        Parameters
        ----------
        **kwargs : dict
            Arbitrary keyword arguemts list.
        typestr : int
            ID of DOF type, e.g. 8 for displacement
        typeID : int
            ID of nodal DOF for x axis, e.g. 1,2,3 for x,y,z direction 4 for 
            roation around x-axis
        xtypestr : str
            typestr, e.g. 'diplacement'
        xtypeID : int
            ID of xaxis, e.g. 18 for frequency
        exponent : float
            Exponent of the data, used for Signal power e.g. 2 for p^2
        xdata_exponent': floot
            Exponent of data-axis exponent, e.g. -1 spectral density p^2/f 


        Raises
        ------
        ValueError
            Unkown values.

        Returns
        -------
        None.

        """
        

        self._type = int(1)
        self._xtype = int(1)
        self._exp = 1
        self._xexp = 0
        self._ureg = ureg.dimensionless
        self._xureg = ureg.dimensionless
        xexp_sw = False # switch for making sure that xexp is initialized
    

        for kw in kwargs:
            if kw == 'typestr':
                self._type = gettypeID(kwargs[kw])
            elif kw == 'typeID':
                self._type = np.round(kwargs[kw]) # guarantee integers
            elif kw == 'xtypestr':
                self._xtype = gettypeID(kwargs[kw])
                if not(xexp_sw): # use 1 is xexp is not used
                    self._xexp = 1
            elif kw == 'xtypeID':
                self._xtype = int(kwargs[kw])
            elif kw == 'exponent':
                self._exp = kwargs[kw]
            elif kw == 'xdata_exponent':
                self._xexp = kwargs[kw]
                xexp_sw = True # avoid overriding by default value
            else:
                 raise ValueError('Unkown keyword argument {0}'.format(kw))
                 
        self._ureg = ureg.parse_expression(SIUNITDICT[self._type])
        self._xureg = ureg.parse_expression(SIUNITDICT[self._xtype])
                    
    # helpers
    def _equal_exp_args(self,other):
        """
        Ckecks if exp amd xexp attribute equal

        Parameters
        ----------
        other : DOFtype
            objekt to be compared.

        Returns
        -------
        Bool
            True if equal.

        """              
        
        return np.allclose(self._exp,other._exp) and np.allclose(self._xexp,other._xexp) 
        
    def _set_types(self):
        """
        Sets the types according to ureg values
        

        Returns
        -------
        None.

        """
        
        self._type = typeIDfromLMT(self.LMT)
        self._xtype = typeIDfromLMT(self.xLMT)
    
    
    def __str__(self):
        """
        overloaded __str__ method

        Returns
        -------
        str
            of DOF.

        """        
        
        _str = self.typestr
        if self._exp != 1:
            _str += '**{0}'.format(self._exp)
        if self._xexp != 0:
            _str += '/{0}**{1}'.format(self.xtypestr,self._xexp)
        if self.ureg != ureg.dimensionless:
            _str += ' in {0}'.format(str(self.ureg.units))
        if self.xureg != ureg.dimensionless:
            _str += '/{0}'.format(str(self.xureg.units))
            

        return _str

    def __repr__(self):
        """
        overloaded __repr__ method

        Should allow class construction, must be improved

        Returns
        -------
        str
            of DOF.

        """
        return "DOFtype(typestr='{0}')".format(self.typestr)

    
    def qlabel(self):
        """
        

        Returns
        -------
        _str : TYPE
            DESCRIPTION.

        """
        
        _str = self.typestr
        if self._exp != 1:
            _str += '**{0}'.format(self._exp)
        if self._xexp != 0:
            _str += '/{0}**{1}'.format(self.xtypestr,self._xexp)
        return _str

    def ulabel(self):
        """
        Provide unit label

        Returns
        -------
        _str : str
            unit str 'meter'.

        """
        
        _str = ""
        if self.ureg != ureg.dimensionless:
            _str += ' ({0}'.format(str(self.ureg.units))
        if self.xureg != ureg.dimensionless:
            _str += '/{0})'.format(str(self.xureg.units))
        else:
            _str += ')'
        return _str
    
    def label(self):
        """
        Provides label
        
        Returns
        -------
        _str : str
            Full label, e.g. 'displacement / meter' .

        """
        
        _str = self.qlabel()
        if self.ureg != ureg.dimensionless:
            _str += "/" + self.ulabel()
        return _str
    

        

    def __eq__(self,other):
        """
        Overloaded __eq__ 

        Parameters
        ----------
        other : DOFtype.

        Returns
        -------
        bool

        """
        
        return all(self.LMT==other.LMT) and all(self.xLMT==other.xLMT)
    
    def __neq__(self,other):
        """
        Overloaded __neq__ 

        Parameters
        ----------
        other : DOFtype.

        Returns
        -------
        bool

        """
        return not(self==other)
    
    def __len__(self):
        """
        Overloaded __len__

        Returns
        -------
        int
            always 1 because units are scalar.

        """
        
        return 1
                
    @property
    def typestr(self):
        """
        

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return TYPEDICT[self._type]
    @property
    def xtypestr(self):
        """
        property method

        Returns
        -------
        str
            xtypestr.

        """
        return TYPEDICT[self._xtype]
    
    @property
    def type(self):
        """
        property method
        
        Returns
        -------
        int
            type.

        """
        return self._type

    @property
    def ureg(self):
        """
        property method

        Returns
        -------
        ureg
            unit registry.

        """
        return self._ureg**self._exp
    
    @property
    def xureg(self):
        """
        property method
        
        Returns
        -------
        ureg
            unit registry for x-axis.

        """
        return self._xureg**self._xexp

    @property
    def fulltype(self):
        """
        Special ID that includes all types and exponents generated from LMT
        
        Returns
        -------
        int
            Integer with specific digits for L, M or T

        """
        LMT = self.fullLMT
        return int(10*(10000*(LMT[0]+50)+100*(LMT[1]+50)+LMT[2]+50))
 
    @property    
    def fullLMT(self):
        """
        LMT as combination of nominator and denominator

        Length - Mass - Time units

        Returns
        -------
        ndarray
            np.array([L,M,T]).

        """
        return self.LMT-self.xLMT
    
    @property
    def LMT(self):
        """
        Length, Mass, Time, exponents of instance

        Returns
        -------
        ndarray
            np.array([L,M,T]).

        """
        dim = self._ureg.dimensionality
        L = dim['[length]']
        M = dim['[mass]']
        T = dim['[time]']
        return np.array([L,M,T])*self._exp

    @property
    def xLMT(self):
        """
        Length, Mass, Time, exponents of x-axis of instance

        Returns
        -------
        ndarray
            np.array([L,M,T]).
        """
        
        dim = self._xureg.dimensionality
        L = dim['[length]']
        M = dim['[mass]']
        T = dim['[time]']
        return np.array([L,M,T])*self._xexp
    
    @property
    def dBref(self):
        """
        decibel reference value linked to unit

        Returns
        -------
        _ret : float
            reference value, e.g. 20e-6 Pa for pressure.

        """
        
        
        # @todo consider denominator and nominatr
        
        _ret =  DBREF[self._type]
        
        return _ret
    
    def __mul__(self,other):
        """
        Multiplication for DOFtype
        
        Takes care of unit multiplication e.g. 'force' * 'length' = 'work'

        Parameters
        ----------
        other : Doftype
            factor.

        Returns
        -------
        DOFtype
            self * other.

        """
        
        # check nominator denomnator cancellation
        if self==other:
            resDOF = copy.deepcopy(self)
            resDOF._exp *= 2
            resDOF._xexp *= 2
            resDOF._set_types()
            return resDOF
        # if self nominator is cancelled by other denominator  
        elif all(self.LMT==other.xLMT):
            return DOFtype(typeID=other._type, exponent=other._exp,xtypeID = self._xtype, xdata_exponent=self._xexp)
        # if self denominator is cancelled by other nominator  
        elif all(self.xLMT==other.LMT):
            return DOFtype(typeID=self._type, exponent=self._exp,xtypeID = other._xtype, xdata_exponent=other._xexp)
        else: # rely on the pint capapilities
            res_LMT  = self.LMT+other.LMT
            res_xLMT = self.xLMT+other.xLMT
            _type  = typeIDfromLMT(res_LMT)
            _xtype = typeIDfromLMT(res_xLMT)
            if _type == 0 or _xtype == 0: # one is unknown so try total
                res_LMT -= res_xLMT
                _type = typeIDfromLMT(res_LMT)
                if _type > 0:
                    return DOFtype(typeID=_type)
                else: # just rely on pint logics
                    resDOF = DOFtype()
                    resDOF._ureg = self._ureg**self._exp*other._ureg**other._exp/ \
                                  (self._xureg**self._xexp*other._xureg**other._xexp)
                                  
                    resDOF._set_types()
                    return resDOF
            else:
                return DOFtype(typeID=_type,xtypeID=_xtype)
            
    def __imul__(self,other):
        """
        Selfmultiplication of DOFtype

        Parameters
        ----------
        other : DOFtype
            factor.

        Returns
        -------
        DOFtype
            self *= other.

        """
        
        self = self*other
        return self
                                  
    def __truediv__(self,other):
        """
        Division for DOFtype
        
        Takes care of unit division e.g. 'force' / 'area' = 'pressure'

        Parameters
        ----------
        other : DOFtype
            divisor.

        Returns
        -------
        DOFtype
            self/other.
        """
        
        # if  not self._equal_exp_args(other):
        #    raise ValueError('Both arguments must have equal exponent attributes')

        if self==other:
            return DOFtype(typeID=0)
        elif all(self.LMT==other.LMT):
            return DOFtype(typeID=other._xtype, exponent=other._xexp, xtypeID = self._type, xdata_exponent=self._xexp)
        elif all(self.xLMT==other.LMT):
            return DOFtype(typeID=self._type, exponent=self._exp,xtypeID = other._type, xdata_exponent=other._exp)
        else: # rely on the pint capapilities
            res_LMT  = self.LMT-other.xLMT
            res_xLMT = other.LMT-self.xLMT
            _type  = typeIDfromLMT(res_LMT)
            _xtype = typeIDfromLMT(res_xLMT)
            if _type == 0 or _xtype == 0: # one is unknown so try total
                res_LMT -= res_xLMT
                _type = typeIDfromLMT(res_LMT)
                if _type > 0:
                    return DOFtype(typeID=_type)
                else: # just rely on pint logics
                    resDOF = DOFtype()
                    resDOF._ureg = self.ureg/other.ureg
                    # exponent logic cannot be kept with division
                    resDOF._exp = 1
                    if self._xexp != 0:
                        resDOF._xureg = self.xureg/other.xureg
                        resDOF._xexp = 1
                        
                    resDOF._set_types()
                    
                    return resDOF
            else:
                return DOFtype(typeID=_type,xtypeID=_xtype)

        #if all(self.xLMT==other.xLMT):
        #    return DOFtype(typeID=self._type, exponent=self._exp, xtypeID = other._type, xdata_exponent=other._exp)

                                  
    def __floordiv__(self,other):
        return self/other
    
    def __add__(self,other):
        if self != other:
            raise ValueError('Arguments must be of same DOFtype to be added')
        else:
            return self
        
    def __sub__(self,other):
        if self != other:
            raise ValueError('Arguments must be of same DOFtype to be subtracted')
        else:
            return self
        
    def __pow__(self,other):
        """
        power of DOFtype
        
        This method affects only the exponents

        Parameters
        ----------
        other : float
            exponent.

        Raises
        ------
        ValueError
            Wrong input type.

        Returns
        -------
        resDOF : DOFtype
            self**other.

        """
        
        if np.isscalar(other) and np.isreal(other):
            resDOF = DOFtype(typeID=self._type, exponent=self._exp*other, \
                             xtypeID=self._type, xdata_exponent=self._xexp*other)
            # make sure that ureg is done fine in all cases
            # if ID is unknown
            resDOF._ureg = self._ureg
            resDOF._xureg = self._xureg
            
            return resDOF
        else:
            raise ValueError('other must be a scalar real numeric')
            
            
                                  
        
        
        
        
        
    
    
    
    
     
        

    
    