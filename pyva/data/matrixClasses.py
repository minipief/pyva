# -*- coding: utf-8 -*-
"""
Module for matrix and vector handling

The matrix classes modules contains classes that deal with numerical tubular
data. This comprises state vectors and system matrices that can change over 
certain frequency or time. 

The data_axix class deals with the dimension over which the vectors and matrices
are changing. 

The Signal class deals with vectors where every coefficient represent the degree
of freedom, for example the complex pressure amplitude at different positions in 
space.

The LinearMatrix class deals with the matrices. The first two dimensions are the rows
and columns of system matrices, the third dimension is the DataAxis dimension.
One main task of the LinearMatrix is to pick out a system matrix at one specific
frequency.

In other words, the LinearMatrix is not real three dimensional dataset, a better
description would be the linear externsion of a set of matrices.
Most numpy methods are implemented and used in such a way, that the operations are
performed along the full depth index.

The DynamicMatrix is an extension of the LinearMatrix class adding excitation and 
response degrees of freedom.  

.. _fig-linmat-sig:
    
.. figure:: ./images/LinMat_Sig.*
   :align: center
   :width: 90%
   
   Sketch of numeric data classes in pyva.

"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import scipy.linalg as linalg
#import numpy.linalg as linalg
import copy
from cycler import cycler


import pyva.data.dof as dof
import pyva.useful as uf




RTOL = 1e-18
ATOL = 1e-24


# useful functions
def isdiagonal(data):
    """
    Check if all matrices in 3D data set are diagonal.

    Parameters
    ----------
    data : ndarry
        3D dataset array.

    Returns
    -------
    bOut : bool
        True for diagonal argument.

    """
    """
    Tests if first two dimensions of a matrix are diagonal
    """
    
    bOut = True
        
    iN = data.shape[2]
    for iz in range(iN):
        mBuf = np.array(data[:,:,iz])
        if np.any(mBuf != np.diag(np.diagonal(mBuf))):
            bOut = False
            return bOut
    return bOut


def linearIndex(N,irow,icol, btriu = True):
    """
    Provides the linear index of triangular matrices

    Parameters
    ----------
    N : int
        Number of columns.
    irow : int
        row index.
    icol : int
        column .
    btriu : bool, optional
        switch for upper triangular index. The default is True.

    Raises
    ------
    IndexError
        Wenn btriu = False the lower triangle is used, which is only possible when 
        the matrix is quadratic.

    Returns
    -------
    int
        linear index in upper triangula matrix.
            
    """
    if btriu:
        return N * irow - irow * (irow + 1) // 2 + icol
    else:
        if irow == icol:
            return irow
        else:
            raise IndexError('Diagonal matrices require irow == icol')

def hermitian(M):
    """
    Calculate hermitian of matrix along first dimension.

    Parameters
    ----------
    M : np.ndarray of size (N,M,M)
        complex 3D array.

    Returns
    -------
    Series of hermitian mastrices.

    """
    return np.conj(np.transpose(M,axes=(0,2,1)))

        
# Class definition and methods

class LinearMatrix:
    """
    Class for handling 3D data matrices in an efficient way.
    
    Attributes
    ----------
    data: ndarry
        3D matrix in 3D or 2D depending on symmetry
    sym: int
        Type of symmetry 0-no 1-sym 2-hermitian 3-diagonal
    shape: tuple
        shape of the matrix
    """
    
    def __init__(self, data, **kwargs):
        """
        Class contructor for LinearMatrix.
        
        When 3D data is provided the constructor checks if the matrix has any symmetry. When 
        2D data is provided the kwords must specify the detailed shape and symmetry
        
        When only 3D data is given the symmetry will be automatically detected.
        
        When the sym argemunet is used the data is given as 2D ndarray using a
        linear index for the upper triangle (sym in {1,2}) or the diagonal (sym = 3)
                
        Parameters
        ----------
        data : ndarray
            3D matrix or 2D when kwargs are used
        **kwargs :
            Arbitrary argument list.
        shape: tuple of int 
            (Nrow,Ncol,Ndepth) shape dimension of the 3D matrix
        sym: int 
            symmetry sym=0,1,2 or 3 for irregular, symmetric, hermitian and diagonal
        
        Raises
        ------
        ValueError
            When data, sym and shape are not consistent.

        Returns
        -------
        None.

        Examples
        --------
            import matrixClasses as mC
            mat3D    = mC.LinearMatrix(3Ddata)
            mat3Dsym = mc.LinearMatrix(2Dtriudata,sym = 1,shape=(3,3,4) )
        """
        if len(kwargs) == 0: # only D given, check dimension
            self._shape = ()
            self._sym = 0
            status = self._check_symmetry(data)
            
            if status == 1:
                raise ValueError('In case of 1dimensional input data, the shape must be given')
            
            #print('{0:13} LinarMatrix of shape {1} created'.format(self.symstr,self.shape))
            
        elif len(kwargs) == 2:
            for kw in kwargs:
                if kw == 'shape':
                     self._shape = kwargs[kw]
                elif kw == 'sym':
                     self._sym  = kwargs[kw]
                else:
                     print('Unkown argument')
                
            # chick input
            if self._sym == 0:
                if data.shape != self._shape:
                    raise ValueError('shape value and data.shape must be equal')
            elif self._sym in (1,2):
                Ntriu = self._shape[0]
                if data.shape != ( (Ntriu*(Ntriu+1))//2 , self._shape[2]):
                    raise ValueError('shape value and data.shape must fit the number of triu matrix coefficients')
            elif self._sym == 3:
                Ndia = self._shape[0]
                if data.shape != (Ndia, self._shape[2]):
                    raise ValueError('shape value and data.shape must fit the number of diagonals')
            
            
            if self._sym in (1,2,3):
                Nx = self._shape[0]
                Ny = self._shape[1]
                if Nx != Ny:
                    raise ValueError('Matrix must be square for sym 1,2 or 3')
                    
            self._data = data
            #print('{0:13} LinarMatrix of shape {1} created'.format(self.symstr,self.shape))
        else:
            raise ValueError('Zero or two kw-arguments required')

    def _check_symmetry(self,data,*shape):
        """
        Checks the symmetry of the input matrix.
        
        Parameter
        ---------
            data: ndarray
                3D numpy array of shape (rows,cols,depth)
            
        Returns
        -------
            status: hint for correct symmetry and size 
            
        """
        # check dimension of data
        if uf.isscalar(data):
            self._shape = (1,1,1)
            data = np.array([[[data]]])
            status = 0
        if len(data.shape) == 1:
            # in this case the dimension cannot be derived and shape must be provided
            status = 1
            raise ValueError('Symmetry of 1 dimensional data cannot by identified automatically')
        if len(data.shape) == 2:
            self._shape = (data.shape[0],data.shape[1], 1)
            data = np.expand_dims(data,2)
            status = 2
        elif len(data.shape) == 3:
            self._shape = data.shape
            status = 3           
        else:
            raise ValueError('Input array has more than 3 dimensions')
          
        # check principle possiblily of symmetry
        ix,iy = data.shape[0:2] 
        if ix != iy:
            self._sym = 0
            self._data = data
        else:
            
            if isdiagonal(data):
                sym = 3 # test this first because they are symetric too!
                iD  = np.diag_indices(ix)
            elif  np.allclose(data, data.transpose(1,0,2),rtol=RTOL, atol=ATOL):
                sym = 1
                iD  = np.triu_indices(ix)
            elif np.allclose(data, data.transpose(1,0,2).conj(),rtol=RTOL, atol=ATOL):
                sym = 2
                iD  = np.triu_indices(ix)
            else:
                sym = 0
                iD = tuple(slice(s) for s in data.shape)
                
        
            #self._data = np.array(data[iD])
            self._data = data[iD] # version where LM and np.array share same data!
            self._sym = sym

        self._shape = data.shape
        return status

            
    def __len__(self):
        return self._shape[0]*self._shape[1]*self._shape[2]
    
    def __str__(self):
        return "LinearMatrix of size {0}, sym: {1}".format(self.shape,self._sym)

    def __repr__(self):
        str = "LinearMatrix of size {0}, sym: {1}\n".format(self.shape,self._sym)
#        for irow in range(self.shape[0] if self.shape[0] < 5 else 5):
#            #str = str + "[" +"{0:6 .3f}, sym: {1}".format(self.data[irow,1:5,0])
#            str = str + "{:6.3f}".format(*self.data[irow,0:5,0]) + "\n"
#        return str + "First matrix up to index 5 at iz = 0 \n{0}".format(self.data[0:5,0:5,0])
        return str + "First matrix up to index 5 at iz = 0 \n{0}".format(np.round(self.data[0:5,0:5,0],2))


    @property
    def shape(self):
        """  Shape of LinearMatrix.

        Returns
        -------
        tuple
            (Nrow,Ncol,Ndepth) shape vector of Lineax matrix.
        """
        return self._shape
    
        
    @property
    def Nrow(self):
        """
        Number of rows

        Returns
        -------
        int
            Number of matrix rows.

        """
        
        return self._shape[0]
    
    @property
    def Ncol(self):
        """
        Number of columns

        Returns
        -------
        int
            Number of columns.

        """
        return self._shape[1]
        
    @property
    def Ndepth(self):
        """
        Number of matrices
        
        size of third dimension

        Returns
        -------
        int
            Number of matrices.

        """
        return self._shape[2]
    
    @staticmethod
    def zeros(sym, shape, **kwargs):
        """
        zeros matrix creator with defined dimension
        
        Parameters
        ----------
        sym : int
            Symmetry identifier.
        shape : tuple of int 
            (Nrow,Ncol,Ndepth) shape.
        **kwargs : dict
            Arbitrary argument list passed to np.zeros method.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        LinearMatrix
            zero matrix with defined dimension and shape.

        """
        
        if sym == 0:
            return LinearMatrix(np.zeros(shape,**kwargs),sym=sym,shape=shape)
        elif sym in (1,2):
            if shape[0] == shape[1]:
                nn = shape[0]
                N = nn*(nn+1)//2
                return LinearMatrix(np.zeros((N,shape[2]),**kwargs),sym=sym,shape=shape)
            else:
                raise ValueError('First two dimensions must be equal for sym = 1,2')
        elif sym == 3:
            if shape[0] == shape[1]:
                nn = shape[0]
                N = nn
                return LinearMatrix(np.zeros((N,shape[2]),**kwargs),sym=sym,shape=shape)
            else:
                raise ValueError('First two dimensions must be equal for sym = 3')
        else:
            raise ValueError('Unknown symmetrie argument')
            
        
    @property
    def sym(self):
        """
        Property method for sym

        Returns
        -------
        int
            Symmetry identifier.

        """
        return self._sym


    @property
    def symstr(self):
        """
        Property method for symstr

        Returns
        -------
        str
            Symmetry string.

        """
        if self._sym == 1:
            return 'symmetric'
        elif self._sym == 2:
            return 'hermitian'
        elif self._sym == 3:
            return 'diagonal'
        else:
            return 'non-symmetric'
            
    def copy(self):
        """
        copy method for LinearMatrix
        
        uses deepcopy

        Returns
        -------
        Linear Matrix
            Deepcopy instance of LinearMatrix instance.

        """
        
        return copy.deepcopy(self)
        
    
    @property
    def data(self):
        """
        data property method
        
        The method recreates the 3D data from the efficiently stored 2D data 
        in case of sym in (1,2,3)
        
        Returns
        -------
        ndarray
            3D data.

        """
        
        # @todo use the data format of data !!! dtype = self._data.dtype
        if self._sym == 0:
            return self._data
        else:
            mBuf = np.zeros(self.shape, dtype = np.complex128)
            lBuf = np.zeros(self.shape[0:2], dtype = np.complex128)
            id = np.diag_indices(self.shape[0])
            Nz = self.shape[2]
            if self._sym == 1:
                iu = np.triu_indices(self.shape[0])
                for iz in range(Nz): 
                    lBuf[iu] = self._data[:,iz]
                    mBuf[:,:,iz] = lBuf + lBuf.T - np.diag(np.diag(lBuf))
                
            elif self._sym == 2:
                iu = np.triu_indices(self.shape[0])
                for iz in range(Nz): 
                    lBuf[iu] = self._data[:,iz]
                    mBuf[:,:,iz] = lBuf + lBuf.T.conj() - np.diag(np.diag(lBuf)).conj()
            elif self._sym == 3:
                for iz in range(Nz): 
                    mBuf[id[0],id[1],iz] = self._data[:,iz]
            return mBuf

    def Dindex(self,iz):
        """
        Provides 2D matrix of in-depth index iz
        
        MM.Dindex(iz) == MM.data(:,:,iz)
    

        Parameters
        ----------
        iz : int
            index in depth.

        Returns
        -------
        ndarray
            2D matrix.

        """
        
        if self._sym == 0:
            return self._data[:,:,iz]
        else:
            lBuf = np.zeros(self.shape[0:2], dtype = np.complex128)
            if self._sym == 1:
                iu = np.triu_indices(self.shape[0])
                #il = np.tril_indices(self.shape[0],-1)
                lBuf[iu] = self._data[:,iz]
                lBuf += lBuf.T - np.diag(np.diag(lBuf))
                #lBuf[il] += lBuf.T[il]
                
            elif self._sym == 2:
                iu = np.triu_indices(self.shape[0])
                #il = np.tril_indices(self.shape[0],-1)
                lBuf[iu] = self._data[:,iz]
                lBuf += lBuf.T.conj() - np.diag(np.diag(lBuf)).real
                #lBuf[il] += lBuf.T[il].conj()
            elif self._sym == 3:
                ix = np.diag_indices(self.shape[0])
                lBuf[ix[0],ix[1]] = self._data[:,iz]
            return lBuf


    def __getitem__(self, position):
        """
        Overloaded method getitem for LinearMatrix
        
        Deals with specific indexing issues espacially when 1D data is accessed
        """
        
        if np.size(np.shape(self.data[position])) < 3:
            # Catch case when there are singleton dimensions and dimension must be extended
            _shape = [1,1,1]
            ix = 0
            _data = self.data[position]
            for ipos in position:
                _helpindex = [0,0,0] # np.array((0,0,0),dtype = np.int64)
                if isinstance(ipos,slice): # this is the dimension that carries the non-singleton data
                    _helpindex[ix] = ipos
                    _shape[ix] = np.size(self.data[tuple(_helpindex)])
                # old version leading to problems when 1 singleton dimansion was already kept and then doubled 
                #elif np.size(ipos) == 1: # is singleton fill dimension
                #    _data = np.expand_dims(_data,axis = ix)
                # else:
                #     raise IndexError('This should not happen')
                    
                ix += 1
                
            return LinearMatrix(_data.reshape(_shape),sym = 0,shape = tuple(_shape))
        else:
        
            return LinearMatrix(self.data[position])

    def __setitem__(self, position, value):
        try:
            index = self._get_index(position)
        except IndexError:
            # Create new matrix because of broken symmetry
            self._data  = self.data
            self._sym   = 0
            self._shape = self.shape
            index = self._get_index(position)
         
        
        # handlye specific behaviour of hermitian matrices
        if self._sym == 2:
            if position[0] > position[1]:
                print('Lower triangle used')
                self._data[index] = np.conj(value)
            elif position[0] == position[1]:
                if np.iscomplex(value):
                    print('Complex value on diagonal brakes symmetry')
                    # Create new matrix because of broken symmetry
                    self._data  = self.data
                    self._sym   = 0
                    self._shape = self.shape
                    index = self._get_index(position)
                    
                self._data[index] = value
        else:
            self._data[index] = value
            
        return self
            

    def _get_index(self, position):
        """
        Provides index into full,triangular or diagonal indexing
        """
        
        if len(position) < 3:
            row,column = position
            Npos = 2
        else:
            row, column, depth = position
            Npos = 3
            
        if self._sym == 0: # keep index as is
            return position
        elif self._sym in (1,2):
            if column < row: # make sure that only upper triangle is indexed
                row, column = column, row
                
            #index = self.Ncol * row - (0 + row) * (row + 1) // 2 + column
            index = linearIndex(self.Ncol,row,column)
            
            if Npos == 3:
                return index, depth
            else:
                return index
        elif self._sym == 3:
            if row == column:
                index = row 
                if Npos == 3:
                    return index, depth
                else:
                    return index
            else:
                raise IndexError('This index is not valid for diagonal matrices')
        else:
            pass
        
    def _check_dimension(self, other, cross_index=((0,0),(1,1),(2,2)) ):
        """
        This method determines type and dimension for several algebraic operation
        It returns a tuple for each dimension with         
           (sigrow,sigcol,sigdepth)
           
        Each sigX indicates...
            0 when both dimensions are equal
            1 when the first (self)  is singleton
            2 when the second (other) is singteton
           -1 when they don't agree at all

        This method is used for methods as matrix multiplication. Singleton dimensions are
        exended to the shape of the other argument.
        
        Some operations, e.g. the dot-product require that the number of columns of the first 
        matrix are equal to the rows of the second. (Inner dimensions) 
        This is considered  by the crossIndex argument that
        changes the index of the second argument to be compared.
        For example (1,0,2) would suit for the dot operation.
     
        
        Args:
            other:      2nd argument
            cross_index: definition of index of other that shall be compared
                               
        Returns:
            signature tupel (sigrow,sigcol,sigdepth)          
        
        """
        
        res = np.zeros(3, dtype = np.int32 )
        
        if uf.isscalar(other):
            other_shape = (1,1,1)
        else:
            other_shape = other.shape            
        
        for ishape,ci in enumerate(cross_index):
            if self.shape[ci[0]] != other_shape[ci[1]]:
                if self.shape[ci[0]] == 1:
                    res[ci[0]] = 1
                elif other_shape[ci[1]] == 1:
                    res[ci[0]] = 2
                else:
                    res[ci[0]] = -1
                    
        return res
        
    def __neg__(self):
        return -1*self
    
    def __pos__(self):
        return self
        
    def __add__(self,other):
        if uf.isscalar(other):
            return LinearMatrix(self.data + other)
        else:
            (sigrow,sigcol,sigdepth) = self._check_dimension(other)
            if sigcol == -1 or sigrow == -1:
                    raise ValueError('Row and column dimensions must fit')
            else:
                if sigdepth == 0:
                    return LinearMatrix(self.data + other.data)
                elif sigdepth == 1:
                    Ndepth = other.Ndepth
                    data   = np.zeros(other.shape,dtype = np.complex128) 
                    for i3 in range(Ndepth):
                        data[:,:,i3] = self.data[:,:,0] +  other.Dindex(i3)
                    return LinearMatrix(data)
                elif sigdepth == 2:
                    Ndepth = self.Ndepth
                    data   = np.zeros(self.shape,dtype = np.complex128) 
                    for i3 in range(Ndepth):
                        data[:,:,i3] = other.data[:,:,0] +  self.Dindex(i3)
                    return LinearMatrix(data)
                elif sigdepth == -1:
                    raise ValueError('Depth dimension must fit or 1 must be singular')
                else:
                    print(sigdepth)
                    raise ValueError('I have no idea what happened')
                

    def __radd__(self,other):
        if uf.isscalar(other):
            return LinearMatrix(self.data + other)
        else:
            return other + self  

    def __iadd__(self,other):
        if uf.isscalar(other):
            return LinearMatrix(self.data + other)
        else:
            return self + other

    def __sub__(self,other):
        return self + (-1)*other  

    def __mul__(self,other):
        if uf.isscalar(other):
            return LinearMatrix(self._data * other, sym=self._sym, shape=self._shape)
        else:
            return LinearMatrix(self.data * other.data)  

    def __truediv__(self,other):
        if uf.isscalar(other):
            return LinearMatrix(self._data / other, sym=self._sym, shape=self._shape)
        else:
            return LinearMatrix(self.data/other.data)  
        
    def __rmul__(self,other):
            return LinearMatrix(self._data * other, sym=self._sym, shape=self._shape)



    def dot(self,other):
        
        if uf.isscalar(other):
            print('I should not be in the dot produkt')
            return LinearMatrix(self._data * other, sym=self._sym, shape=self._shape)
        else:
            resT = self._check_dimension(other,cross_index=((1,0),(2,2)) )
        
            if np.any(resT < 0):
                raise TypeError("Inner matrix dimensions must agree")
            else:
                matbuf = np.zeros((self.shape[0],other.shape[1],max(self.shape[2],other.shape[2])),dtype = np.complex128)
                if resT[2]==0: # third dimension agrees
                    for iz in range(self.shape[2]):
                        #matbuf[:,:,iz] = self.Dindex(iz).dot(other.Dindex(iz))
                        matbuf[:,:,iz] = np.matmul(self.Dindex(iz),other.Dindex(iz))
                else: # singleton matrix is applied to all others
                    if resT[2]==1:
                        ix2 = range(other.shape[2])
                        ix1 = [0] * other.shape[2]
                        _N  = other.shape[2]
                    else:
                        ix1 = range(self.shape[2])
                        ix2 = [0] * self.shape[2]
                        _N  = self.shape[2]
                    for iz in range(_N):
                        #matbuf[:,:,iz] = self.Dindex(ix1[iz]).dot(other.Dindex(ix2[iz]))
                        matbuf[:,:,iz] = np.matmul(self.Dindex(ix1[iz]),other.Dindex(ix2[iz]))
                        
                return LinearMatrix(matbuf)
            
    def HDH(self,other):
        """
        HDH performt the left and right matrix multiplication by the inverse of other from
        left and right in a numerical stable way using solve 
        
        Args:
            other: NOT inverted matrix D 
                               
        Returns:
            D^-1 . self . D^-H

        """

        out = other.solve(self) # Solve D . X = self => out = D^-1 . self 
        return other.solve(out.H()).H()

            
    def sqrt(self):
        """
        square root.

        Returns
        -------
        LinearMatrix
            sqrt(self).

        """
        return LinearMatrix(np.sqrt(self.data))

    def solve(self,other):
        """
        Solve linear matrix equation.
        
        With self=A and other=B the method returns the solution of A.X = B

        Parameters
        ----------
        other : 
            dependent variable of matrix equation.

        Raises
        ------
        TypeError
            DESCRIPTION.

        Returns
        -------
        LinearMatrix 
            X, the solution of A.X = B

        """
              
        resT = self._check_dimension(other)
    
        #check if matrix is square
        if self.Nrow != self.Ncol:
            raise TypeError("Matrix must be square")
            
        if self.Nrow != other.Nrow:
            raise TypeError("other Ncols must equal self.Ncols")
        else:
            
            matbuf = np.zeros((other.shape[0],other.shape[1],max(self.shape[2],other.shape[2])),dtype = np.complex128)

            if resT[2]==0: # third dimension agrees
                for iz in range(self.shape[2]):
                    matbuf[:,:,iz] = np.linalg.solve(self.Dindex(iz),other.Dindex(iz))
            else: # singleton matrix is applied to all others
                if resT[2]==1:
                    ix2 = range(other.shape[2])
                    ix1 = [0] * other.shape[2]
                    _N  = other.shape[2]
                else:
                    ix1 = range(self.shape[2])
                    ix2 = [0] * self.shape[2]
                    _N  = self.shape[2]
                for iz in range(_N):
                    matbuf[:,:,iz] = np.linalg.solve(self.Dindex(ix1[iz]),other.Dindex(ix2[iz]))
                    
            return LinearMatrix(matbuf)

    def abs(self):
        """
        Abxolute value of Linear Matrix.

        Returns
        -------
        LinearMatrix
            with all absolute value of coefficients.

        """
        if self.sym in (1,2):
            # abs makes hermitian matrixes symmetric
            _sym = 1
        else:
            _sym = self.sym
            
        return LinearMatrix(np.abs(self._data),sym = _sym, shape = self.shape)
    
    def sum(self):
        
        """ 
        Calculates the full sum of all coefficients in the matrix
        """
        
        _sumdata = np.zeros((self.Ndepth,))
        
        for iz in range(self.Ndepth):
            _sumdata[iz] = np.sum(self.Dindex(iz))
            
        return _sumdata
    
    def mean(self):
        
        """ 
        Calculates the mean value of all coefficients in the matrix
        """
        
        _data = np.zeros((self.Ndepth,))
        
        for iz in range(self.Ndepth):
            _data[iz] = np.mean(self.Dindex(iz))
            
        return _data        
    

    def real(self):
        if self.sym in (2,):
            # abs makes hermitian matrixes symmetric
            _sym = 1
        else:
            _sym = self.sym
            
        return LinearMatrix(np.real(self._data),sym = _sym, shape = self.shape)

    def imag(self):
        if self.sym in (2,):
            # imag makes hermitian matrixes antisymmetric
            return LinearMatrix(np.imag(self.data))
           
        else:
            _sym = self.sym
            return LinearMatrix(np.imag(self._data),sym = _sym, shape = self.shape)



        
    def inv(self):
        """ 
        Inverse of linear matrix
        """
        
        _data = self._data.astype(np.complex128)
        
        if self.sym == 3:
            #return LinearMatrix(1/self._data,sym=3,shape=self.shape)
            _data = 1./self._data
        elif self.sym in [1,2]: # inverse if hermitian and symmetric is also symmetric
            iu = np.triu_indices(self.shape[0])
            for iz in range(self.Ndepth):
                lBuf   = np.linalg.inv(self.Dindex(iz))
                _data[:,iz] = lBuf[iu]
        else:
            for iz in range(self.Ndepth):
                _data[:,:,iz] = np.linalg.inv(self.Dindex(iz))
        
        return LinearMatrix(_data,shape = self.shape, sym = self.sym)

    def transpose(self):
        """ 
        Transpose of linear matrix
        
        Returns 
        """
        
        if self.sym in [1,3]:
            return self
        elif self.sym in [2]:
            _data = self._data.copy()
            _data = self._data.conj()
            return LinearMatrix(_data,shape = self.shape, sym = self.sym)
        else: # means 0
            _data = np.zeros((self.Ncol, self.Nrow, self.Ndepth),dtype = np.complex128)
            for iz in range(self.Ndepth):
                _data[:,:,iz] = self.Dindex(iz).transpose()
            
            return LinearMatrix(_data)

    def diag(self):
        """ 
        Transpose of linear matrix
        
        Returns 
        """
        _data = np.zeros((min((self.Ncol, self.Nrow)),self.Ndepth),dtype = np.complex128)
        for iz in range(self.Ndepth):
            _data[:,iz] = np.diag(self.Dindex(iz))
            
        return LinearMatrix(_data,sym=3,shape=self._shape)
        

    def H(self):
        """ 
        Hermitian of linear matrix
        """
        
        
        
        if self.sym in [2]:
            return self
        elif self.sym in [1,3]:
            # this works for symmetric too, because the diagonal must be real
            #_data = self._data.copy()
            _data = np.zeros(self._data.shape,dtype = np.complex128)            
            _data = self._data.conj()
            return LinearMatrix(_data,sym=self._sym,shape=self._shape)     

        else: # means 0
            _data = np.zeros((self.Ncol, self.Nrow, self.Ndepth),dtype = np.complex128)            
            for iz in range(self.Ndepth):
                _data[:,:,iz] = self.Dindex(iz).conj().transpose()
  
            return LinearMatrix(_data)     

    def cond(self,p=None):
        """ 
        cond of linear matrix
        """

        _ydata = np.zeros((1,self.Ndepth))
        
        for iz in range(self.Ndepth):
            _ydata[0,iz] = linalg.cond(self.Dindex(iz),p)

        return _ydata


            
    def isa(self,typestr):
        return typestr in ('LinearMatrix')
    
    def spy(self,res='mag',iz=0,**kwargs):
        """
        Show sparasity of matrix in a graphical way. Uses matplotlib.pyplot.spy
           
       
        Args:
            res: type of result 'mag','real','imag' or 'phase'
            iz:  index of z-dimension
            kwargs: for passing to matplotlib routine 
                               
        
        """        

        if res == 'mag':
            MM = np.abs(self.Dindex(iz))
        elif res == 'imag':
            MM = np.imag(self.Dindex(iz))
        elif res == 'real':
            MM = np.real(self.Dindex(iz))
        elif res == 'phase':
            MM = np.angle(self.Dindex(iz))*180/np.pi
        else:
            raise ValueError("res paramteter must be 'mag','real','imag' or 'phase'")
            
        plt.spy(MM,**kwargs)
            
            
        
     
        
class DynamicMatrix(LinearMatrix):
    """
        
           data: 3D matrix or 2D when kwargs are used
          xdata: 

            
                    
    """    
    def __init__(self,data,xdata,excdof,resdof,**kwargs):
        """
        The DynamicMatrix class extends LinearMatrix by difinitions of the 
        row and coloum degrees of freedom (excdof and resdof) as far as the 
        definition for the in-depth axis.
        
        Parameters
        ----------
        data : ndarray or LinearMatrix
            3D matrix or 2D when kwargs are used
        xdata : DataAxis
            in depth axis definition e.g. frequency, time.
        excdof : DOF
            DOF of excitatin (columns).
        resdof : DOF
             DOF of response  (rows).
        **kwargs : dict
            Arbitrary argument list passed to LinearMatrix constructor.
        shape: tuple of int 
            (Nrow,Ncol,Ndepth) shape dimension of the 3D matrix
        sym: int 
            symmetry sym=0,1,2 or 3 for irregular, symmetric, hermitian and diagonal

        Raises
        ------
        ValueError
            DESCRIPTION.

        Examples
        --------        
            ::
                
                import matrixClasses as mC
                dynmat3D = mC.DynamicMatrix(3Ddata,xdata,resdof)
                dynmat3Dsym = mc.DynamicMatrix(2Dtriudata,xdata,sym = 1,shape=(3,3,4) )


        """
        # data is already LinearMatrix
        if isinstance(data,LinearMatrix):
            super().__init__(data._data,sym=data.sym,shape=data.shape)
        else:
            super().__init__(data,**kwargs)
        
        if xdata.isa('DataAxis'):
            Nxdata = len(xdata)
            self._xdata = xdata
        else:
            raise ValueError('xdata must be an instance of DataAxis')
            
        (Nres,Nexc,Nx) = self.shape
        
        if Nx != Nxdata:
            raise ValueError('3rd dimension of LinearMatrix must equal length of xdata')
            
        if isinstance(excdof,dof.DOF):
            if len(excdof) == Nexc:
                self._excdof = excdof
            else:
                raise ValueError('excdof must be equal to columns of LinearMatrix')
        else:
            raise ValueError('excdof must be an instance of DOF')
            
        if isinstance(resdof,dof.DOF):
            if len(resdof) == Nres:
                self._resdof = resdof
            else:
                raise ValueError('resdof must be equal to rows of LinaerMatrix')
        else:
            raise ValueError('resdof must be an instance of DOF')
            
    @property
    def xdata(self):
        """
        Property method for xdata

        Returns
        -------
        DataAxis
            xdata.

        """
        
        return self._xdata

    @property
    def resdof(self):
        """
        Property method for resdof

        Returns
        -------
        DOF
            response DOFs.

        """
        return self._resdof

    @property
    def excdof(self):
        """
        Property method for excdof

        Returns
        -------
        DOF
            excitatin DOFs.

        """
        return self._excdof

    def __getitem__(self, position):
        """
        Overloaded method getitem for DynamicMatrix
        
        Deals with specific indexing issues espacially when 1D data is accessed
        """
        
        l_data= super().__getitem__(position)
        resdof_ = self.resdof[position[0]]
        excdof_ = self.excdof[position[1]]
        xdata_  = self.xdata[position[2]]

        return DynamicMatrix(l_data,xdata_,excdof_,resdof_)

            
    def __str__(self):
        """
        str method of DynamicMatrix

        Returns
        -------
        _str : str
            Text description of DynamicMatrix.

        """
        
        _str = super().__str__() + '\n'
        _str += str(self._xdata) + '\n'
        _str += 'resdof: '+ str(self._resdof) + '\n'
        _str += 'excdof: '+ str(self._excdof)  
        return _str
        
    def __repr__(self):
        return str(self) + '\n'
        
    def __add__(self,other):
        """
        Overloaded add method of Dynmatrix
        
        In contrast to standard matrix addition this method considers
        the degrees of freedom. Thus, only those coefficients are added to the
        coefficients of self that share the same DOFs.

        Parameters
        ----------
        other : DynamicMatrix
            summand.

        Raises
        ------
        ValueError
            When xaxis don't agree or not all DOF of other are found in self.

        Returns
        -------
        self + other.

        """
                
        # determine joined numer of DOFs
        outresdof,ires,_ = self._resdof.intersect(other._resdof)
        outexcdof,iexc,_ = self._excdof.intersect(other._excdof)
        
        # Throw error when some non common dofs are found
        # Debug commands
        #print(outresdof,ires)
        #print(outexcdof,iexc)
        
        # Check common xdata or interpolate
        if self.xdata != other.xdata:
            raise ValueError("Both dynamic matrix must have consistent xdata")

            
        # Check if symmetry is similar
        
        if self.sym == other.sym:
            
            # Create new data from self
            outObj = copy.copy(self)
            
            if self.sym == 0:
                                
                for iz in range(self.Ndepth):
                    for irow in range(len(ires)):
                        for icol in range(len(iexc)):
                            outObj._data[ires[irow],iexc[icol],iz] += other.data[irow,icol,iz]
            elif self.sym in (1,2):
                for iz in range(self.Ndepth):
                    for irow in range(len(ires)):
                        for icol in range(irow,len(iexc)):
                            outObj._data[ires[irow],iexc[icol],iz] += other.data[irow,icol,iz]
                
        else:
            if self.sym == 0:
                
                outObj = copy.copy(self)
                
                for iz in range(self.Ndepth):
                    for irow in range(len(ires)):
                        for icol in range(len(iexc)):
                            outObj._data[ires[irow],iexc[icol],iz] += other.data[irow,icol,iz]
                            
                return(outObj)                
    
    def __radd__(self,other):
        """
        Overloaded method radd for DynamicMatrix

        Should not be used!!!!!

        Parameters
        ----------
        other : DynamicMatrix
            summand.

        Returns
        -------
        DynamicMatrix
            self+other.

        """
        
        raise ValueError('This should not happen')
        
        return other+self

    def __iadd__(self, other):
        """
        Incremental add for DynamicMatrix
        

        Parameters
        ----------
        other : DynamicMatrix
            summand.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        .. todo::
            clarify the symmetry thing


        """
        

        # determine joined numer of DOFs
        outresdof,ires,_ = self._resdof.intersect(other._resdof)
        outexcdof,iexc,_ = self._excdof.intersect(other._excdof)
        
        # Throw error when some non common dofs are found
        #print(outresdof,ires)
        #print(outexcdof,iexc)
        
        # Check common xdata or interpolate
        if self.xdata != other.xdata:
            raise ValueError("Both dynamic matrix must have consistent xdata")

        # Check if symmetry is similar
        
        if self.sym == other.sym:
            if self.sym == 0:
                for iz in range(self.Ndepth):
                    for irow in range(len(ires)):
                        for icol in range(len(iexc)):
                            self._data[ires[irow],iexc[icol],iz] += other.data[irow,icol,iz]

            elif self.sym in (1,2):
                for iz in range(self.Ndepth):
                    for irow in range(len(ires)):
                        for icol in range(irow,len(iexc)):
                            self._data[self._get_index( (ires[irow],iexc[icol],iz) )] += other.data[irow,icol,iz]
                
        else:
            if self.sym == 0:
                for iz in range(self.Ndepth):
                    for irow in range(len(ires)):
                        for icol in range(len(iexc)):
                            self._data[ires[irow],iexc[icol],iz] += other.data[irow,icol,iz]
                            
        return self
            
        
    def isa(self,typestr):
        """
        Check if self is an insance of typestr 

        Parameters
        ----------
        typestr : str
            str of classtype.

        Returns
        -------
        bool
            True if typestr = 'DynamicMatrix'.

        """
        
        return typestr in ('DynamicMatrix')

            
    def dot(self,other):
        """
        Dot product method of DynamicMatrix
        
        This method allows flexible matrix multiplication. 
        If other is a DynamicMatrix the resdof of other must equal the excdof of self.
        If other is a Signal all DOFs of the Signal vector must occur in the excdof of self.

        Parameters
        ----------
        other : DynamicMatrix, Signal 
            factor.

        Raises
        ------
        ValueError
            When xdata or dofs are not consistent.

        Returns
        -------
        DynamicMatrix
            product.

        """
        

        if other.isa('DynamicMatrix'):
            # solve equation if xdata is a frequency type 
            if self.xdata == other.xdata and self._excdof == other._resdof:

                    self_ = super().dot(other)
                    return DynamicMatrix(self_.data,self.xdata,other.excdof,self.resdof)
            else:
                raise ValueError('The matrix multiplication of DynamicMatrix instances require self.excdof == self.resdof')
                                        
        elif other.isa('Signal'): # required for excitation 
            # todo check xdata!    
            comdof,indexF,ix2 = self._excdof.intersect(other.dof)
            _F                = np.zeros(self.Nrow, dtype = np.complex128)
            _R     = np.zeros((self.Nrow,other.Nxdata),dtype = np.complex128 )

            
             #loop over frequency
            for ifreq in range(len(other.xdata)):
                _M          = self.Dindex(ifreq)
                _F[indexF]  = other.ydata[:,ifreq]
                _R[:,ifreq] = _M.dot(_F) # works also for scalars
                
            return Signal(self._xdata, _R, dof = self._resdof)
        else:
            raise ValueError('DynamicMatrix require Signal or DynamicMatrix instances')
            
            
    def self_dot(self,other):
        """
        Incremental dot product
        
        Modifies self

        Parameters
        ----------
        other : DynamicMatrix
            factor matrix.

        Raises
        ------
        ValueError
            When xdata or dofs are not consistent.

        Returns
        -------
        DynamicMatrix
            matrix product.

        """
        
        if other.isa('DynamicMatrix'):
            # solve equation if xdata is a frequency type 
            if self.xdata == other.xdata and self._excdof == other._resdof:

                self = super().dot(other)
                self.excdof = other.excdof 
            else:
                raise ValueError("excdof of self must match resdof of other")
                
        return self
    
    @staticmethod
    def zeros(xdata,excdof,resdof,sym,**kwargs):
        """
        zeros matrix creator with defined dimension and DOFs
        
        Parameters
        ----------
        xdata : DataAxis
            in depth axis definition e.g. frequency, time.
        excdof : DOF
            DOF of excitatin (columns).
        resdof : DOF
             DOF of response  (rows).
        sym : int
            Symmetry identifier.
        **kwargs : dict
            Arbitrary argument list passed to np.zeros method.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        DynamicMatrix
            zero matrix with defined DOFs and symmetry.

        """
        
        Ndepth = len(xdata)
        Nrow   = resdof.size
        Ncol   = excdof.size
        
        shape  = (Nrow,Ncol,Ndepth)
        
        data = LinearMatrix.zeros(sym, shape, **kwargs)
        
        return DynamicMatrix(data, xdata, excdof, resdof)
    
    def signal(self,irow,icol,exc_val = 1.):
        """
        create a Signal object from 1-dim parts of the DynamicMatrix

        Parameters
        ----------
        icol : TYPE
            DESCRIPTION.
        irow : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        sig_data = self.data[irow,icol,:]*exc_val
        sig_DOFs = self.resdof[irow] 
        
        return Signal(self.xdata,sig_data,sig_DOFs)
            
    def inv(self):
        """
        Inverts dynamic matrix
        
        excitation and response dofs are exchanged

        Returns
        -------
        DynamicMatrix
            inverse(self).

        """
                
        super().inv()
        self._resdof,self._excdof = self._excdof,self._resdof
        return self

    def solve(self,other):
        """
        Solve dynamic matrix
        
        Returns the solutoin of A.x = B with A=self and B=other
        
        Excitation dofs are dofs of the solution. All DOFs of B must be a subset 
        of res_dof of the Dynamic Matrix


        Parameters
        ----------
        other : Signal
            vector B.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        Signal
            res.

        """
        
        
        # Code snipped from dot taking signal part only
        if self.xdata == other.xdata:
              
            comdof,indexF,ix2 = self._resdof.intersect(other.dof)
            _F  = np.zeros(self.Nrow, dtype = np.complex128)
            _R  = np.zeros((self.Nrow,other.Nxdata),dtype = np.complex128 )
        
            #loop over third dimension
            for ifreq in range(len(other.xdata)):
                _M          = self.Dindex(ifreq)
                _F[indexF]  = other.ydata[:,ifreq]
                #_R[:,ifreq] = np.linalg.solve(_M,_F)
                try:
                    _R[:,ifreq] = linalg.solve(_M,_F)
                except np.linalg.LinAlgError as le:
                    print('{0} at index {1} with data {2} = {3}'.format(le,ifreq,str(other.xdata.type.typestr),other.xdata.data[ifreq]))
                
            return Signal(self._xdata, _R, dof = self._excdof)
        else:
            raise ValueError('xdata of self and other must be the same')



# Class definition and methods

class DataAxis:
    """ 
    Class for handling 1D data axis e.g. the time of requency axis
    
    Attributes
    ----------
        data: ndarray
            data axis samples 
        type: of data e.g. time, frequency @type str
    """

    def __init__(self, data, **kwargs):
        """
        Class constructor of data:axis

        Parameters
        ----------
        data : ndarray
            data axis samples.
        **kwargs : dict
            Arbitrary argument list.
        typestr: str
            type string e.g. 'pressure
        typeID: int
            corresponding typeID
        type: DOFtype
            DOFtype
        
            
        Raises
        ------
        ValueError
            DESCRIPTION.


        """
        
        
        
        if isinstance(data,DataAxis):
            DataAxis(data._data,type=data._type)
        else:
            self._data = data
            self.__class__ = DataAxis
    
            if len(kwargs) == 0: # only D given, check dimension
                self._type = dof.DOFtype()
               
            else:      
                for kw in kwargs:
                    if kw == 'typestr':
                        self._type = dof.DOFtype(typestr=kwargs[kw])
                    elif kw == 'typeID':
                        self._type = dof.DOFtype(typeID=kwargs[kw])
                    elif kw == "type":
                        self._type = kwargs[kw]
                    else:
                         raise ValueError('Unkown keyword {0}'.format(kw))
                     
            
    @property
    def shape(self):
        """
        shape method of DataAxis object

        Returns
        -------
        tuple
            (1,Ndata) .

        """
        
        return (1,len(self._data))

    def __len__(self):
        """
        Overloaded len method of DataAxis

        Returns
        -------
        int
            length of DataAxis object.

        """
        
        return len(self._data)

    
    def __str__(self):
        return "DataAxis of {0} samples and type {1}".format(len(self),self.type)

    def __repr__(self):
        return str(self)

    def __setitem__(self, position, value):
        """
        __setitem__ of DataAxis

        Parameters
        ----------
        position : int or slice
            index.
        value : float
            values to set.

        Returns
        -------
        None.

        """
        
        self._data[position] = value
        
    def __eq__(self,other):
        sw = self._type == other._type
        return sw and np.allclose(self._data,other._data)
        
    def __ne__(self,other):
        return not(self==other)
                
    @property
    def data(self):
        """
        Property method for data

        Returns
        -------
        ndarray
            data axis samples.

        """
        
        return self._data

    @property
    def angular_frequency(self):
        """
        Provides the angular frequency independent from frequency type

        Returns
        -------
        ndarray
            angular frequency
            
        Raises
        ------
        ValueError
            When self is no frequency type.


        """
        if self.type.type == 18:
            return 2*np.pi*self._data
        elif self.type.type == 21:
            return self.data
        else:
            raise ValueError('Angular frequency requires a frequency data type')

    @property
    def frequency(self):
        """
        Provides the frequency independent from frequency type

        Returns
        -------
        ndarray
            frequency

        Raises
        ------
        ValueError
            When self is no frequency type.

        """
        if self.type.type == 21:
            return self._data/2/np.pi
        elif self.type.type == 18:
            return self.data
        else:
            raise ValueError('Angular frequency requires a frequency data type')


    @property
    def unit(self):
        """
        Provides units of DataAxis

        Returns
        -------
        units
            units of dataaxis.

        """
        
        return self._type._ureg.units

    @property
    def type(self):
        """
        DOFtype of DataAxis

        Returns
        -------
        DOFtype 
            DOFtype of DataAxis.

        """
        
        return self._type
    
    @property
    def typeID(self):
        """
        

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        
        return self._type.typeID
           
    def label(self):
        """
        Creates text label of DataAxis
            
        Returns
        -------
        str
            DataAxis label.

        """
        
        return self._type.label()
    
    def isa(self,typestr):
        """
        Check is self is of type typestr

        Parameters
        ----------
        typestr : str
            typestr.

        Returns
        -------
        bool 
            True if typestr = 'DataAxis'.

        """
        
        return typestr in ('DataAxis')
        
    def __getitem__(self, position):
        """
        getitem method for DataAxis

        Parameters
        ----------
        position : int or slice
            index.

        Returns
        -------
        DataAxis
            indexed result.

        """
        
        return DataAxis(self.data[position],type = self._type)
    
    @staticmethod
    def octave_band (f_min=2*np.pi*250, f_max=2*np.pi*4000, bands_per_octave = 3, typeID = 21):
        """
        Mathod for octave / log space creation
        
        Parameters
        ----------
        f_min : float
            lower frequency limit
        f_max : float
            upper frequency limit
        bands_per_octave : integer, optional
            frequency points per octave band. The default is 3.
        typeID : integer, optional
            typeID 18 frequency, 21 angular frequency. The default is 21.

        Returns
        -------
        DataAxis with center log spaces center frequencies

        """
        
        fa1 = 2**(1/bands_per_octave) 
        
        if typeID == 21:
            pi_fac = 2*np.pi
        elif typeID == 18:
            pi_fac = 1
        else:
            raise(ValueError,'typeID must be 18 or 21')
                
        f_min = f_min/pi_fac
        f_max = f_max/pi_fac
        
        # first index
        i0 = np.round(np.log10(f_min/1000)/np.log10(fa1)) # -3 from log10(1000)
        # number of lines
        N  = np.round((np.log10(f_max)-3)/np.log10(fa1))-i0+1
        ix = np.arange(0,N)
        freqc = pi_fac*1000*(fa1**(i0+ix))  
                      
        return DataAxis(freqc,typeID = typeID)
                      
class Signal:
    """
    Class for test and simulation results

    Attributes
    ----------
    xdata: DataAxis
        x-axis samples 
    ydata: ndarray or function 
        array of Signal samples with dimension Ndof x Nxdata with ndarray output of (Ndofs,Nxdata) output
    dof: DOF 
        object with ID, local DOF and DOFtype of ydata channels
    
    """
    def __init__(self, xdata, ydata, dof):
        """
        Constructor for Signal classes
        
        Parameters
        ----------
        xdata : DataAxis
            in depth data of Signals e.g. time, frequency
        ydata : ndarray of size (Ndof,Nxdata) or list of function handles
            DESCRIPTION.
        dof : TYPE, optional
            DESCRIPTION. The default is dof.DOF(0,0,dof.DOFtype(typeID = 0)).

        Raises
        ------
        ValueError
            input is checked for consistency.

        Returns
        -------
        None.

        """
        
        # Check consistency of x- and ydata
        if ydata.ndim < 2:
            Nsig   = 1
            Nxdata = len(ydata)
            ydata = np.array([ydata])
        else:
            Nsig   = ydata.shape[0]
            Nxdata = ydata.shape[1]
            
        if Nxdata == len(xdata) and Nsig == len(dof):
            self._ydata = ydata
            self._xdata = xdata
            self._Nsig   = Nsig
            self._Nxdata = Nxdata
        else:
            raise ValueError('ydata, xdata and dof must be of consistent dimenension')
            
        self._dof = dof
        
    def __str__(self):
        return "Signal of {0} samples and {1} DOFs".format(self._Nxdata,len(self._dof))

    def __repr__(self):
        _str = 'Signal object of {0} samples with {1} channels and properties ...\n'.format(self.Nxdata,self.Nsig)
        _str += str(self._xdata) + '\n'
        _str += str(self._dof) + '\n'
        
        return _str
    
    def __len__(self):
        """
        len of Signal

        Returns
        -------
        int
            length of Signal.

        """
        
        return len(self._xdata)
    
    @property
    def ydata(self):
        
        if isinstance(self._ydata, np.ndarray):
            return self._ydata
        else:
            _ydata = np.zeros((len(self._dof),len(self)),dtype = np.complex128)
            for idof in range(len(self._dof)):
                _ydata = self._ydata[idof](self._xdata.data)
                
            return _ydata
        
    @property
    def xdata(self):
        return self._xdata

    def getdof(self):
        return self._dof
    
    def setdof(self,value):
        if isinstance(value, dof.DOF):
            self._dof = value
 
    dof = property(getdof,setdof)
    
 
    @property
    def Nxdata(self):
        """
        Number of DataAxis samples

        Returns
        -------
        int
            Number of samples.

        """
        
        return len(self)
    
    @property
    def Nsig(self):
        """
        Number of Signal channels

        Returns
        -------
        int
            Number of Signal DOFs.

        """
        
        return self._Nsig

    def conj(self):
        _sig = copy.deepcopy(self)
        _sig._ydata = np.conj(self._ydata)
        return _sig

    def __abs__(self):
        _sig = copy.deepcopy(self)
        _sig._ydata = np.abs(self._ydata)
        return _sig

    def real(self):
        _sig = copy.deepcopy(self)
        _sig._ydata = np.real(self._ydata)
        return _sig

    def sum(self):
        """
        sum of Signals
        
        For reasonable results DOFtype must not change 

        Returns
        -------
        Signal
            sum of all DOFs per x-axis.

        """
        
        _ydata = np.sum(self._ydata,axis=0)
        return Signal(self.xdata, _ydata, dof = dof.DOF(0,0,self.dof.type[0]))
        

    def imag(self):
        _sig = copy.deepcopy(self)
        _sig._ydata = np.imag(self._ydata)
        return _sig
    
    def interp(self,xvalues,ix=0):
        """
        interpolation for Signals
        
        Provides interpolated values when sample is not given

        Parameters
        ----------
        xvalues : double
            xvalues for interpolation from Signal
        ix : integer, optional
            Signal index. The default is 0.

        Returns
        -------
        ndarray
            interpolated ydata values for Signal of index ix

        """
        
        return(np.interp(xvalues,self.xdata.data,self.ydata[ix,:]))
    
    def nan2zero(self):
        """
        sets all 'NaN's to zero
        

        Returns
        -------
        Signal
            with NaNs replaced by 0.

        """
        
        aydata = self._ydata
        index = np.isnan(aydata)
        
        self._ydata[index] = 0.
        
        return self
    
    def inf2zero(self):
        """
        sets all infs to zero
             

        Returns
        -------
        Signal
            with inf's replaced by 0.

        """
        
        aydata = self._ydata
        index = np.isinf(aydata)
        
        self._ydata[index] = 0.
        
        return self    
    
    def __pow__(self,exponent):
        _sig = copy.deepcopy(self)
        _sig._ydata = self._ydata**exponent
        for idof in range(len(self._dof._type)):
            _sig._dof._type[idof]   = self._dof._type[idof]**exponent
        return _sig
    
    def isa(self,typestr):
        return typestr in ('Signal')

    def __getitem__(self, position):
        """
        Overloaded method getitem for Signal
        
        Deals with specific indexing issues especially when 1D data is accessed
        """
        
        _ydata = self._ydata[position]
        
        if np.size(np.shape(position)) > 1:
            _xdata = self.xdata(position[1])
            _pos1  = position[0]
        else:
            _xdata = self.xdata
            _pos1  = position

        return Signal(_xdata,_ydata,dof = self._dof[_pos1])


    def __rmul__(self,other):

        # check xaxis dimension
        
        if uf.isscalar(other):
            _sig = copy.deepcopy(self)
            _sig._ydata = _sig._ydata * other

        return _sig
    
    def __neg__(self):
        return self*(-1)
    
    def __mul__(self,other):
        
        if uf.isscalar(other):
            return other*self # use already implemented rmul
        if isinstance(other, DataAxis):
            if self.xdata == other:
                _sig = copy.deepcopy(self)
                for isig in range(self.Nsig):
                    _sig._ydata[isig,:]*= other.data
                    _sig._dof.type[isig] *= other.type
            return _sig
        if self.xdata == other.xdata:
            # check number of Signals
            
            # 0 dof allows for multiplication because it is like scalar
            if all(self.dof.dof == 0) or all(other.dof.dof== 0):
                ignore_args = ('doftype','dof')
            else:
                ignore_args = 'doftype'
                    
            
            if self.Nsig == other.Nsig:
                _sig = copy.deepcopy(self)
                index = self.dof.index(other.dof,ignore=ignore_args)
                
                if any(index == -1):
                    raise ValueError('Not all DOFs in self consistent to those in other')
                
                for isig in range(self.Nsig):
                    _sig._ydata[isig,:]*= other._ydata[index[isig],:]
                    _sig._dof.type[isig] *= other._dof.type[index[isig]]
                    
        return _sig

    def __add__(self,other):

        if self.xdata != other.xdata:
            raise ValueError("data must (currently) be the same")
        if self._dof != other._dof:
            raise ValueError("DOF must (currently) be the same")
            
        _sig = copy.deepcopy(self)
                
        _sig._ydata += other._ydata
                    
        return _sig
    
    def __sub__(self,other):
        
        return self + -1*other

    
    def __truediv__(self,other):
        
        if other.isa('DynamicMatrix'):
            # solve equation if xdata is a frequency type 
            if self.xdata == other.xdata and self.xdata.type.type in (21,18):
                
                comdof,indexF,ix2 = other._excdof.intersect(self.dof)
                _F                = np.zeros(other.Nrow, dtype = np.complex128)
                _R     = np.zeros((other.Nrow,other.Ndepth),dtype = np.complex128 )

                
                 #loop over frequency
                for ifreq in range(len(other.xdata)):
                    _M = other.Dindex(ifreq)
                    _F[indexF] = self.ydata[:,ifreq]
                    _R[:,ifreq] = linalg.solve(_M,_F) # works also for scalars
                return Signal(self._xdata, _R, dof = other._resdof)
                    
        elif other.isa('Signal'):
            
            if self.xdata == other.xdata and self.xdata.type.type in (21,18):
                comdof,indexF,ix2 = self.dof.intersect(other.dof,ignore = 'doftype')
                _res = copy.deepcopy(self)
                
                if len(indexF) == len(ix2):
                    for isig in range(len(indexF)):
                        _res._ydata[indexF[isig],:] = self._ydata[indexF[isig],:]/other._ydata[isig,:]
#                        _res._dof[indexF[isig]]._type[0] = self._dof[indexF[isig]]._type[0]/other._dof[isig]._type[0]
                        _res._dof[indexF[isig]]._type[0] = self._dof._type[indexF[isig]]/other._dof._type[isig]
                elif other.Nsig == 1:
                    _res = self._ydata/other._ydata
                    _res._dof = self._dof/other.dof
                
                return _res                    

    
    def volume_flow(self,ID,dofarg=1,S=1.):
        """
        Exctracts the volume flow from Signals
        
        Parameters
        ----------
        ID : int
            node ID.
        dofarg : int, optional
            local degree of freedom. The default is 1.
        S : float, optional
            DESCRIPTION. The default is 1..

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        Signal
            of volume flow DOFs.

        
        """
        
        index = self.dof.index(dof.DOF(ID,dofarg,dof.DOFtype(typestr='volume flow')))
        if index < 0:
            raise ValueError('No volume flow for ID={0} and DOF={1} found'.format(ID,dofarg) )
                             
        return Signal(self._xdata, self.ydata[index,:] , dof=dof.DOF(ID,dofarg,dof.DOFtype(typestr='volume flow')))
        

        
    def pressure(self,ID):
        """
        Exctracts the pressure from Signals
        
        Parameters
        ----------
        ID : int
            node ID.

        Returns
        -------
        Signal
            extracted pressure Signal.

        """
        
        index = self.dof.index(dof.DOF(ID,0,dof.DOFtype(typestr='pressure')))
        return Signal(self._xdata, self.ydata[index,:] , dof=dof.DOF(ID,0,dof.DOFtype(typestr='pressure')))
    
    def transfer(self,other,IDs=None):
        """
        Calculates the transfer function of similar types between differnt DOFs
        
        Parameters
        ----------
        other : Signal
            input Signal.
        IDs : tuple or list of int, optional
            output and input ID. The default is [0,1].

        Raises
        ------
        ValueError
            IN case of not consitent data.

        Returns
        -------
        Signal
            transfer coefficient, self/other.

        """
        
        if self.xdata != other.xdata:
            raise ValueError('Both Signals must have consistent xdata')
        
        
        if IDs == None:
            ix1 = 0
            ix2 = 0
        else:
            ix1 = self.dof.find(ID = IDs[0])
            ix2 = other.dof.find(ID = IDs[1])
            
        print(ix1,ix2)
            
        # check equality of doftype
        if self.dof[ix1]._type == self.dof[ix2]._type:
            _ydata = self._ydata[ix1,:]/other._ydata[ix2,:]
        else:
            raise ValueError('Both Signals must have same DOFtype')
        
        return Signal(self.xdata,_ydata,dof.DOF(0,0,dof.DOFtype(typestr = 'transmission')))
        
        
        
            
        
        

        
    def plot(self,nfig,loc=1,res='mag',xscale='linear',yscale='linear', **kwargs):
        """
        Plot method for Signal obejcts.
        
        Parameters
        ----------
        nfig : TYPE
            number of figure.
        loc : int or str, optional
            location of legend. The default is 1.
        res : str, optional
            kind of result mag, phase, real, imag. The default is 'mag'.
        xscale : str, optional
            linear or log. The default is 'linear'.
        yscale : str, optional
            linear or log. The default is 'linear'.
                    ls : line style 
        **kwargs : 
            Arbitrary argument list.
        ID : IDs to plot
        dof : DOFs to plot e.g. [1,2]
        legstr : additional label added to automatic label 
        fulllegstr : list of labels overwriting the automatic label
        DX : int 
            value for skipping x samples DX=2 skips every second sample
        lw : list or tuple
            line width
        cs : list or tuple
            colorstyle
        ms : list or tuple
            markerstyle
        ls : list or tuple
            line_styles
            
        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        None.

        """
        #
        iDOF = range(self._Nsig)
        dofs = -1          # default value to ignore dofs in find
        repetition = True
        counter_ID_dof = 0
        leg_sw    = False # True when own labels are used
        style_sw  = False # True when line style is specified
        _lw = 2
        DX = 1

        ID_sw  = False       # When True, filter with dof and IDs
        IDs    = self.dof.ID # Preset IDs to all existing

        leg_auto_sw = True
        label_sw = True      # false if no label should be used
        # Default line styles applied in the cycler definition
        line_styles = ['-', '--', ':', '-.']
        marker_styles = ['', '.', 'o']
        color_styles = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', \
                        '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', \
                        '#bcbd22', '#17becf', '#17becf', '#17becf']
                
        plt.figure(nfig)
        ax = plt.gca()
                
        for kw in kwargs:
            if kw == 'ID':
                IDs = kwargs[kw]
                ID_sw = True
                counter_ID_dof += 1
            elif kw == 'dof':
                dofs = kwargs[kw]
                ID_sw = True
                counter_ID_dof += 1
            elif kw == 'legstr':
                leg_sw = True
                leg_pre_str = kwargs[kw]  
            elif kw == 'fulllegstr':
                leg_sw = True
                leg_auto_sw = False
                if kwargs[kw] == '':
                    leg_str = kwargs[kw] # not used
                    label_sw = False
                    leg_sw = False
                else:
                    leg_str = kwargs[kw]
            elif kw in ('ls','linestyle'): #second key for compatibility
                line_styles = kwargs[kw]
                style_sw = True
            elif kw in ('cs','color'): #second key for compatibility
                color_styles = kwargs[kw]
                style_sw = True
            elif kw in ('ms','marker'): #second key for compatibility
                marker_styles = kwargs[kw]
                style_sw = True
            elif kw =='lw':
                _lw = kwargs[kw]
            elif kw =='DX':
                DX = kwargs[kw]
            else:
                 raise ValueError('Unkown keyword {0}'.format(kw))
                 
        # create cycler from parameters
        if style_sw:
            color_c = cycler('color', color_styles)
            style_c = cycler('linestyle', line_styles)
            markr_c = cycler('marker', marker_styles)
            c_ms = markr_c * style_c
            min_len = min(len(c_ms),len(color_c))
            custom_cycler = color_c[:min_len] + c_ms[:min_len]
            
            # set ax style cycler
            ax.set_prop_cycle(custom_cycler)
        


        if counter_ID_dof > 1: # set repetition to True because ID and dof are selceted specifically
            repetition = False
 
        
        if ID_sw: #len(kwargs) > 0:
            iDOF = self.dof.find(IDs,dofs,repetition = repetition)
        
        # check consistent arguments
        if leg_sw:
            if leg_auto_sw:
                if len(iDOF) != len(leg_pre_str):
                    raise(ValueError('Number of legpre elements must be equal to plotted Signals'))
            else:
                if len(iDOF) != len(leg_str):
                    raise(ValueError('Number of legpre elements must be equal to plotted Signals'))        
        
        if not(leg_sw):
            leg_str = []
            leg_pre_str = []
            for ileg in range(len(iDOF)):
                leg_pre_str.append('')
                leg_str.append('')
                

        # Manage y-axis scaling and labeling
        ylabel = self.dof[iDOF[0]].type[0].qlabel()        
        if np.iscomplexobj(self._ydata):
            if res == 'mag':
                pfun = lambda x: np.abs(x)
                ylabel = '|{0}/{1}|'.format(ylabel,self._dof.type[0].ulabel())
            elif res == 'dB':
                dBref,dBpot,refstr = self.dof.dBref
                pfun = lambda x: dBpot*10*np.log10(np.abs(x)/dBref)
                ylabel = 'dB re {0}'.format(refstr)
            elif res == 'real':
                pfun = lambda x: np.real(x)
                ylabel = 'real({0}/{1})'.format(ylabel,self._dof.type[0].ulabel())
            elif res == 'imag':
                pfun = lambda x: np.imag(x)
                ylabel = 'imag({0}/{1})'.format(ylabel,self._dof.type[0].ulabel())
            elif res == 'phase':
                pfun = lambda x: np.angle(x)
                ylabel = 'phase({0})/deg'.format(ylabel)
            else:
                raise ValueError('Wrong res parameter {0}'.format(res))
        else:
            if res == 'dB':
                dBref,dBpot,refstr = self.dof.dBref
                pfun = lambda x: dBpot*10*np.log10(np.abs(x)/dBref)
                ylabel = 'dB re {0}'.format(refstr)
            else:
                pfun = lambda x: x
                #ylabel = self._dof.type[iDOF[0]].label()
                ylabel = self.dof[iDOF[0]].type[0].label()
            
        i_leg = 0
        

        
        for idof in iDOF:
            if leg_auto_sw:
                 _label = leg_pre_str[i_leg]+':'+self._dof[idof].label()
            else:
                 _label = leg_str[i_leg]
                
            line, = plt.plot(self._xdata.data[::DX],pfun(self._ydata[idof,::DX].flatten()), lw=_lw)

            if label_sw:
                line.set_label(_label)
                
            i_leg+= 1
        
        # Zero line
        plt.plot([self._xdata.data[0],self._xdata.data[-1]],[0,0],'k-',lw=0.5)

        
        plt.xlabel(self._xdata.label())
        plt.ylabel(ylabel)
        
        plt.gca().set_yscale(yscale)
        plt.gca().set_xscale(xscale)
        
        plt.legend(loc=loc)
        plt.tight_layout
        plt.show()
        
class ShapeSignal(Signal):
    """
    Class for Signals with link to shapes ydata can be reprented by ndarrays but also function definition

    Attributes:
        xdata:  x-axis samples @type DataAxis
        ydata:  array of Signal samples with dimension Ndof x Nxdata @type np.array or function with ndarray output of (Ndofs,Nxdata(]) output
          dof:  DOF object with ID, local DOF and type of ydata channels
    
    """
    
    def __init__(self,mesh,xdata, ydata, dof):
        """
        constructor for Signal classes with link to 2D meshes
        
        Parameters
        ----------
        mesh : meshtypes
            mesh that links the DOFs to positions 
        xdata : DataAxis
            in depth data of Signals e.g. time, frequency
        ydata : ndarray of size (Ndof,Nxdata) or function with input arguments (pos_x<,posy,posz,>,xdata) depending on mesh position arguments
            DESCRIPTION.
        dof : TYPE, optional
            DESCRIPTION. The default is dof.DOF(0,0,dof.DOFtype(typeID = 0)).

        Raises
        ------
        ValueError
            input is checked for consistency.

        Returns
        -------
        None.

        """     
        
        super().__init__(xdata, ydata, dof = dof )
        self._mesh = mesh
        
    @property
    def ydata(self):
        
        if isinstance(self._ydata, np.ndarray):
            return self._ydata
        else:
            _ydata = np.zeros((len(self._dof),len(self)),dtype = np.complex128)
            if self._mesh.ndim == 2:
                _ydata = self._ydata(self._mesh.X,self._mesh.Y,self._xdata.data)
            else:
                raise ValueError('Other dimensions for meshes as 2 not implemented')
                
            return _ydata
        
    def plot3d(self,fig=1,index = 1):
    
        # This import registers the 3D projection, but is otherwise unused.


        fig = plt.figure(fig)
        ax = fig.gca(projection='3d')

        X, Y = self._mesh.nodes()
        Z = self.ydata[:,index].reshape(X.shape)

        # Plot the surface.
        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0.02, antialiased=True)

        # Customize the z axis.
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

        ax.set_zlim(np.min(Z.flatten())*1.01, np.max(Z.flatten())*1.01)
        ax.set_xlim(self._mesh.X0-0.01, self._mesh.X1+.01)
        ax.set_ylim(self._mesh.Y0-0.01, self._mesh.Y1+.01)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        # Strange code taking care of scaling
        max_range = np.array([X.max()-X.min(), Y.max()-Y.min()]).max() / 2.0 # , Z.max()-Z.min()]

        mid_x = (X.max()+X.min()) * 0.5
        mid_y = (Y.max()+Y.min()) * 0.5
        mid_z = (Z.max()+Z.min()) * 0.5
        ax.set_xlim3d(mid_x - max_range, mid_x + max_range)
        ax.set_ylim3d(mid_y - max_range, mid_y + max_range)
        #ax.set_zlim3d(mid_z - max_range, mid_z + max_range)        

        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.5, aspect=5)

        # ax.auto_scale_xyz([self._mesh._X0, self._mesh._X1], [self._mesh._Y0, self._mesh._Y1], [Z.max(), Z.min()])
        # ax.pbaspect = [mid_x, mid_y, mid_z]
        #ax.set_aspect('equal')

        plt.show()

        
        
        
        

    
    
        


        
        
        
            
            
     
        

                        
            
    
        
        
    
    
        

    
    