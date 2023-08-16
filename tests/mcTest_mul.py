# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 19:20:00 2018
"""

from pyva.data import matrixClasses as matC
import numpy as np
import random as rnd
import timeit 


Nx = 4
Ny = Nx
Nz = 500

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
    
    return np.array([[1., 0.,  0., 0.],
                     [0., cs, -sn, 0.],
                     [0., sn,  cs, 0.],
                     [0., 0.,  0., 1.]],dtype = np.complex128)


# Data for nonsymetric random matrix
Adata = np.round(40*((np.random.random_sample((Nz,Nx,Nx))-0.5)+1j*(np.random.random_sample((Nz,Nx,Nx))-0.5)))
#A = matC.LinearMatrix(Adata)
#Asafe = A.copy()

# Symmetric data
Bdata = (Adata + Adata.transpose(0,2,1))/2
#B     = matC.LinearMatrix(Bdata)
#Bsafe = B.copy()

# Hermitian data
Cdata = (Adata + Adata.transpose(0,2,1).conj())/2
#C = matC.LinearMatrix(Cdata)
#Csafe = C.copy()

Ddata = np.zeros((Nz,Nx,Nx), dtype = np.complex128)
iD = np.diag_indices(Adata.shape[1])

Ddata[:,iD,iD] = Adata[:,iD,iD]

def my_mul(A,B):
    C = np.zeros(A.shape,dtype=np.complex128)
    for i in range(A.shape[0]):
        C[i,:,:] = np.matmul(A[i,:,:],B[i,:,:])
    return C

t_0 = timeit.default_timer()

AmulB = np.matmul(Adata,Bdata)

t_1 = timeit.default_timer()

AmulBtest=my_mul(Adata,Bdata)

t_2 = timeit.default_timer()

# calculate elapsed time and print
elapsed_time = round((t_1 - t_0) * 10 ** 6, 3)
print(f"Elapsed time1: {elapsed_time} µs")
elapsed_time = round((t_2 - t_1) * 10 ** 6, 3)
print(f"Elapsed time2: {elapsed_time} µs")

# test if single mastrix along full Nz works

T = edge_transform(np.pi/2)

TA= np.matmul(T,Adata)


