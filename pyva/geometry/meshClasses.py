# -*- coding: utf-8 -*-
"""
This module provides methods for mesh handling. 
This comprised simple and regular meshes without an element definition 
but will later invlude complex meshes of irregular elements, different 
coordinate systems and node definitions.

The latter is not yet implemeneted. Feel free :)
"""

import numpy as np
import pandas as pd
import types

#from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy import integrate

#import pyva.properties.geometricalPropertyClasses as geoPC
import pyva.data.dof as dof

class RegMesh2D:

    """ 
    This class deals with regular and flat 2D meshes. 
    These meshes are very useful for radiation stiffness calculation of vibrating structures
    
    """
    
    def __init__(self,X0,Y0,X1,Y1,NX,NY,doftype = dof.DOFtype(typestr = 'general')):
        
        """ 
        Class constructor of RegMesh2D
        
        Parameters
        ----------
        X0 : float
            lowest X position
        Y0 : float
            lowest Y position
        X1 : float
            highest X position
        Y1 : float
            highest Y position
        NX : float
            number of nodes in X-direction
        NY : float
            number of nodes in Y-direction
            

        """
       
        self.X0 = X0
        self.Y0 = Y0
        self.X1 = X1
        self.Y1 = Y1
        self.NX = int(NX)
        self.NY = int(NY)
        self.doftype = doftype
        
        # update flags
        self._dist_update = True
                
        self.dX = (X1-X0)/(NX-1)
        self.dY = (Y1-Y0)/(NY-1)

    def __str__(self):
        """
        str for RegMesh2D

        Returns
        -------
        str

        """
        _str  = "RegMesh2D: \n"
        _str += "X0              : {0}\n".format(self.X0)
        _str += "Y0              : {0}\n".format(self.Y0)
        _str += "X1              : {0}\n".format(self.X1)
        _str += "Y1              : {0}\n".format(self.Y1)
        _str += "NX              : {0}\n".format(self.NX)
        _str += "NY              : {0}\n".format(self.NY)
        _str += "doftype         : {0}".format(self.doftype)
        
        return _str
            
    def __repr__(self):
        """
        repr of RegMesh2D

        Returns
        -------
        None.

        """
        
        return "RegMesh2D({0},{1},{2},{3},{4},{5},doftype ={6})"\
            .format(self.X0,self.Y0,self.X1,self.Y1,self.NX,self.NY,self.doftype.__repr__())
    
    @property    
    def area(self):
        """ 
        Returns area of full mesh. inlcuding edge elements

        Returns
        -------
        float
            mesh area.

        """
                
        return (self.X1-self.X0+self.dX)*(self.Y1-self.Y0+self.dY)
        
   
    @property    
    def shape(self):
        return (self.NY,self.NX)
       
    @property    
    def dA(self):
        return self.dX*self.dY
    
    @property    
    def X(self):
        """
        x-coordinates of mesh

        Returns
        -------
        ndarray
            x-node positions.

        """
        
        return np.linspace(self.X0,self.X1,self.NX)

    @property    
    def Y(self):
        """
        y-coordinates of mesh

        Returns
        -------
        ndarray
            y-node positions.

        """
        return np.linspace(self.Y0,self.Y1,self.NY)
            
    @property    
    def ks(self):
        """
        Wavenumber parameter from Langley [Lan2007]_

        Technically with is the wavenumber that korresponds to the 
        smallest allowed wavelength of the mesh represented by the
        diagonal length of the element 

        Returns
        -------
        float
            wavenumber.

        """
        
        return 2*np.pi/np.sqrt(self.dX**2+self.dY**2)
    

    @property    
    def Nmesh(self):
        """
        Number of nodes

        Returns
        -------
        int

        """
        
        return self.NX*self.NY

    
    def nodes(self,index='all'):
        """
        Node coordinates

        Parameters
        ----------
        index : int or str, optional
            index of node. The default is 'all'.

        Returns
        -------
        ndarray
            x,y coordinates of mesh or x of 1st index.
        float 
            y of 2nd index.

        """
        
        if isinstance(index, str ) and index =='all':
            return np.meshgrid(self.X,self.Y)
        else:
            return (self._X0+self._dX*index[0],self._Y0+self._dY*index[1])
                

    def distance(self):
        """
        Calculate the distances between nodes of mesh
        
        This method provides the occuring distances in a mesh in a reduced way.
        Regular meshes have many repeating distances, thus this method returns a unique
        array of distances and an index to reconstruct the upper triangular matrix
        
        dist[index] gives the full triangular matrix
        
        Returns
        -------
        ndarray 
            dist : unique distances 
        nd.adday of int
            index : index to reconstruct all distances from dist
        """
        
        if self._dist_update:
            
            N = self.Nmesh
            
            _X,_Y  = self.nodes()
            _X = _X.reshape(N,1)
            _Y = _Y.reshape(N,1)        
    
            
            # The distance matrix is symmetric, only upper triangle is calculated
            dist2 = np.zeros(((N+1)*N)//2,dtype=np.float32)
            
            i0 = 0
    
            for i in range(N):

                il = np.arange(i+1,N)
                dist2[i0+il-i] = (_X[il,0]-_X[i,0])**2+(_Y[il,0]-_Y[i,0])**2
                i0 += N-i
    
    
            # Filter first occurrenc
            dist2u = pd.unique(dist2)
            
            # Provide sorted index of distances
            idx_sorted = np.argsort(dist2u)
            # Provide index for occurence in full distance list
            ix_pos = np.searchsorted(dist2u[idx_sorted], dist2)
            self._dist  = np.sqrt(dist2u)
            self._index = np.uint32(idx_sorted[ix_pos]) 
    
        return self._dist,self._index
    
    def plot3d(self,fig=1):
        """
        3D plot

        Parameters
        ----------
        fig : int, optional
            figure identifier. The default is 1.

        Returns
        -------
        None.

        """
        
    
        # This import registers the 3D projection, but is otherwise unused.


        fig = plt.figure(fig)
        ax = fig.gca(projection='3d')

        X, Y = self.nodes()
        Z = np.zeros(np.shape(X))

        # Plot the surface.
        #surf = 
        ax.plot_wireframe(X, Y, Z,
                       linewidth=0.5, antialiased=True)

        # Customize the z axis.
        ax.set_zlim(-1.01, 1.01)
        ax.set_xlim(self.X0-0.01, self.X1+.01)
        ax.set_ylim(self.Y0-0.01, self.Y1+.01)

        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))


        plt.show()


class RegShape2D(RegMesh2D):

    """ 
    This class deals with shapes mapped on regular and flat 2D meshes. 
    These meshes are very useful for radiation stiffness calculation of vibrating structures
    
    Check if required, usefull for function based shapes. 
    Note, that there is no frequency depence!
    
    """
    
    def __init__(self,X0,Y0,X1,Y1,NX,NY, doftype = dof.DOFtype(typestr='displacement'), shape = 1.):
        """
        Class constructor of RegShape2D
        
        Parameters
        ----------                  
        X0 : float
            lowest X position
        Y0 : float
            lowest Y position
        X1 : float
            highest X position
        Y1 : float
            highest Y position
        NX : float
            number of nodes in X-direction
        NY : float
            number of nodes in Y-direction
        doftype : DOFtype, optional
            type of shape. The default is dof.DOFtype(typestr='displacement').
        shape : functin, optional
            Function for shape. The default is 1..

        Returns
        -------
        None.

            shape: shape of mesh, 2D fun or ndarray of shape NX,NY function is only defines in mesh area the rest is ignored
            
        Examples:
            >>> regMesh(0,0,0.8,0.9,10,10)

        """  
        
        super().__init__(X0,Y0,X1,Y1,NX,NY)
        
        if isinstance(shape,types.FunctionType):
            X,Y     = self.nodes()
            self.surfshape = shape(X,Y) 
            self._surfshapefun = shape
        else:
            self.surfshape = shape
            self._surfshapefun = 'None' # later interpoliation

    def shapefun(self,X,Y):
        """
        Shape function with consideratin of mesh limits 

        Parameters
        ----------
        X : float
            x-coordinate.
        Y : TYPE
            y-coordinate.

        Returns
        -------
        res : ndarray
            of shape dimension.
        """
        
        
        res = np.zeros(X.shape)
        # allow only values in mesh area
        ix = (X>self.X0 or X<self.X1 or Y>self.Y0 or Y< self.Y1)
        res[ix] = self._surfshapefun(X[ix],Y[ix])
        return res
            
    def normalised_shape(self,N=2):
        """ 
        Normalised reshape considering integration over element area
        
        Integrates the shape function around the node to get better
        results when nodal shape value if not representing the average
        value in the node vincinity
        
        Parameters
        ----------
        N : int
            oversamping for integration. The default is 2.

        Returns
        -------
        ndarray
            normalised shape            
        """
        
        _shape = np.zeros((self.NY,self.NX))
        X,Y = self.nodes()
        
        f = lambda x,y: self._surfshapefun(x,y) # x*y 
        dX2 = self.dX/2
        dY2 = self.dY/2
        
        for ix in range(self.NX):
            for iy in range(self.NY):
                X0 = X[iy,ix]-dX2
                X1 = X[iy,ix]+dX2
                xlim0 = lambda x: X0
                xlim1 = lambda x: X1
                _shape[iy,ix] = integrate.dblquad(f,Y[iy,ix]-dY2,Y[iy,ix]+dY2, xlim0, xlim1)[0]
            
        return _shape/self.dA
            
    
    def mean_square(self,normalise):
        """
        mean square value

        Parameters
        ----------
        normalise : bool
            switch if normalisation shall be used.

        Returns
        -------
        float
            mean square.
        
        """
        
        return 0.5*self.Zvec(normalise).conj() @ self.Zvec(normalise)
       
    def Zvec(self,normalise=False):
        """
        Vector of shape values

        Parameters
        ----------
        normalise : bool, optional
            switch if normalisation shall be used. The default is False.

        Returns
        -------
        ndarray
            shape vector.

        """
        if normalise:
            return self.normalised_shape().flatten()
        else:
            return self.shape.flatten()
            
    def plot3d(self,fig=1):
        """
        3D plot for shapes

        Parameters
        ----------
        fig : int, optional
            figure identifier. The default is 1.

        Returns
        -------
        None.

        """
        
    
        # This import registers the 3D projection, but is otherwise unused.


        fig = plt.figure(fig)
        ax = fig.gca(projection='3d')

        X, Y = self.nodes()
        Z = self.surfshape

        # Plot the surface.
        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0.02, antialiased=True)

        # Customize the z axis.
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

        ax.set_zlim(np.min(Z.flatten())*1.01, np.max(Z.flatten())*1.01)
        ax.set_xlim(self._X0-0.01, self._X1+.01)
        ax.set_ylim(self._Y0-0.01, self._Y1+.01)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        # Strage code taking care of scaling
        max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max() / 2.0

        mid_x = (X.max()+X.min()) * 0.5
        mid_y = (Y.max()+Y.min()) * 0.5
        mid_z = (Z.max()+Z.min()) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)        

        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.5, aspect=5)

        ax.auto_scale_xyz([self._X0, self._X1], [self._Y0, self._Y1], [Z.max(), Z.min()])
        ax.pbaspect = [mid_x, mid_y, mid_z]
        #ax.set_aspect('equal')

        plt.show()
        
