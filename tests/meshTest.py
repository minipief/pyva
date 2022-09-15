"""
Test script for regular meshes
"""

import pyva.geometry.meshClasses as rC
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

Nx = 15 
Ny = 30
X0 = -0.2
Y0 = 0.1
Lx = 0.3
Ly = 0.8

NmodeX = 5
NmodeY = 1

myMesh = rC.RegMesh2D(X0,Y0,X0+Lx,Y0+Ly,Nx,Ny)

print('Area: {:.1f}'.format(myMesh.area))
print('dX= {:.3f} dY= {:.3f}'.format(myMesh.dX,myMesh.dY))
print('Element Area: {:.5f}'.format(myMesh.dA))

print(myMesh.X)
print(myMesh.Y)

X,Y     = myMesh.nodes()
X1,Y1   = myMesh.nodes((4,3))
X1s,Y1s = myMesh.nodes(np.array([[4,4,4,5,5,5],[2,3,4,2,3,4]]))

myMesh.plot3d(1)

dist,index = myMesh.distance()

shapefun = lambda x,y: np.sin(NmodeX*np.pi*(x-X0)/Lx)*np.sin(NmodeY*np.pi*(y-Y0)/Ly)

my_shape = rC.RegShape2D(X0,Y0,X0+Lx,Y0+Ly,Nx,Ny,shape=shapefun)

my_shape.plot3d(2)

fineshape = my_shape.normalised_shape()
my_shape = rC.RegShape2D(X0,Y0,X0+Lx,Y0+Ly,Nx,Ny,shape=fineshape)

my_shape.plot3d(3)

very_simple_shape = rC.RegShape2D(X0,Y0,X0+Lx,Y0+Ly,Nx,Ny)
very_simple_shape.plot3d(4)





