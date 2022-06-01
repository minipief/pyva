# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 19:20:00 2018
"""

from pyva.data import matrixClasses as matC
import numpy as np
import random as rnd

Nx = 4
Ny = Nx
Nz = 6


# Data for nonsymetric random matrix
Adata = np.round(40*((np.random.random_sample((Nx,Nx,Nz))-0.5)+1j*(np.random.random_sample((Nx,Nx,Nz))-0.5)))
A = matC.LinearMatrix(Adata)
Asafe = A.copy()

# Symmetric data
Bdata = (Adata + Adata.transpose(1,0,2))/2
B     = matC.LinearMatrix(Bdata)
Bsafe = B.copy()

# Hermitian data
Cdata = (Adata + Adata.transpose(1,0,2).conj())/2
C = matC.LinearMatrix(Cdata)
Csafe = C.copy()

Ddata = np.zeros((Nx,Nx,Nz), dtype = np.complex128)
iD = np.diag_indices(Adata.shape[0])

Ddata[iD,iD,:] = Adata[iD,iD,:]

D = matC.LinearMatrix(Ddata)
Dsafe = D.copy()

AdotB = A.dot(B)
CdotD = C.dot(D)
Cdot3 = C.dot(3)

# compare with np.results

AdotBdata = np.zeros((Nx,Nx,Nz),dtype = complex)
CdotDdata = np.zeros((Nx,Nx,Nz),dtype = complex)

for iz in range(4):
    AdotBdata[:,:,iz] = Adata[:,:,iz].dot(B.data[:,:,iz])
    CdotDdata[:,:,iz] = Cdata[:,:,iz].dot(D.data[:,:,iz])
    

sqA = A.sqrt

Ateil1 = A[1,1,:]

Dteil  = D[0:2,1:3,3]
DTeil  = D[1,1,1]
DTTeil = D[1,0,:]



print(Ateil1.data)
print(A.data[1,1,:])

if all(Ateil1.data.flatten()==A.data[1,1,:]):
    print('Indexing is working')




# Try setting single values
D[0,0,:] = -3.
D[0,1,0] = 100.

# Try setting values in symmetric matrices
B[0,2,0] = 100.
B[2,0,0] = 100.

C[2,0,0] = 333. + 444.j
C[2,2,0] = 200. + 111.j
# does not work! Kills symmetry and this is not coovered


Ei = matC.LinearMatrix(np.eye(3))

# This should by element wise multiplication
AxB   = A * B
Ax2   = A * 2
C2x   = 2 * C

A0 = A.Dindex(0)
B0 = B.Dindex(0)
A0xB0 = A0*B0

print(A0xB0)
print(AxB.Dindex(0))


A4sum = A.copy()
A4sum += A
A4sum += B
Atest = A4sum - B - A

A,B,C,D = Asafe,Bsafe,Csafe,Dsafe


print('Check A before')
print(A.Dindex(0))
AT = A.transpose()
BT = B.transpose()
CT = C.transpose()
DT = D.transpose()

Asum = A.sum()
ATsum = AT.sum()



ATH = A.H()
BTH = B.H()
CTH = C.H()
DTH = D.H()

# test inversion - results must be near to unity
unityA = A.dot(A.inv())
unityB = B.dot(B.inv())
unityC = C.dot(C.inv())
unityD = D.dot(D.inv())

# check single dimension dot
AsD  = A[:,:,2]
AsDI = A.inv()[:,:,2]

AsDtest = A.dot(AsD).dot(AsDI)

# test column selection
col1 = matC.LinearMatrix(np.array([[1,0,0],[1,0,0],[1,0,0]]))

#crazy things
links = A.transpose().dot(B).dot(C.H())
rechts1 = links.inv()
rechts2 = C.H().inv().dot(B.inv()).dot(A.transpose().inv())

one1 = links.dot(rechts1)
one2 = links.dot(rechts2)
 
print('one1 == one2 is {0}'.format(np.allclose(one1.data,one2.data)))

# check the vector issue
# Here Sff = f . f.H()
# 1: q = D.f Sqq = q.q.H()
# 2: Sqq = D.Sff.D.H()
# 1 == 2 ?

# create Vector
f = A[:,1,:]
Ainv = A# .inv()
q = Ainv.dot(f)
AA = Ainv.inv()

# 1:
Sqq1 = q.dot(q.H())

Sff = f.dot(f.H())
Sqq2 = Ainv.dot(Sff).dot(Ainv.H())

if np.allclose(Sqq1.data,Sqq2.data):
    print('juchhu beide gleich!')
else:
    print('schade!')
    
Sqq3 = Sff.HDH(AA)
    
if np.allclose(Sqq1.data,Sqq3.data):
    print('juchhu beide gleich!')
else:
    print('schade!')





time = matC.DataAxis(np.arange(0,1,0.01),typestr='time')
print(time)
part_time = time[1:40]



