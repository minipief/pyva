# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 19:20:00 2018
"""

from pyva.data import matrixClasses as matC
import numpy as np
import random as rnd
#import pyva.data.data1D as d1D

Adata = np.round(40*((np.random.random_sample((3,3,4))-0.5)+1j*(np.random.random_sample((3,3,4))-0.5)))

A = matC.LinearMatrix(Adata)
Asafe = A.copy()

Bdata = (Adata + Adata.transpose(1,0,2))/2
B     = matC.LinearMatrix(Bdata)
Bsafe = B.copy()

Cdata = (Adata + Adata.transpose(1,0,2).conj())/2
C = matC.LinearMatrix(Cdata)
Csafe = C.copy()

Ddata = np.zeros((3,3,4), dtype = np.complex64)
iD = np.diag_indices(Adata.shape[0])

Ddata[iD,iD,:] = Adata[iD,iD,:]

D = matC.LinearMatrix(Ddata)
Dsafe = D.copy()


AdotC = A.dot(B)
CdotD = C.dot(D)

Cdot3 = C.dot(3)

sqA = A.sqrt

Dteil = D[0:2,1:3,3]
DTeil = D[1,1,1]

Dinv = Dsafe.inv()
print(A.Dindex(0))
A.inv()
B.inv()
C.inv()
print(A.Dindex(0))
Aone = A.dot(Asafe)

Bone = B.dot(Bsafe)
Cone = C.dot(Csafe)
Done = D.dot(Dsafe)

