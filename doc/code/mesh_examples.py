# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 19:44:59 2022

"""

import pyva.geometry.meshClasses as meshC
import pyva.data.dof as dof
import matplotlib.pyplot as plt


# Dimensions
Lx = 1.2
Ly = 1.5
Nx = 6
Ny = 7

dispDOF = dof.DOFtype(typestr='displacement')

my_mesh = meshC.RegMesh2D(0., 0., Lx, Ly, Nx, Ny, dispDOF)

X,Y = my_mesh.nodes()

dist,index = my_mesh.distance()

dist.size

dist[index]

my_mesh.plot3d(1)

#plt.savefig('../source/images/mesh_plot.png')


