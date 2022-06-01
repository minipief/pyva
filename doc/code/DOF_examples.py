# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 20:51:13 2021

"""

import pyva.data.dof as dof
import pyva.data.matrixClasses as mC
import numpy as np
import copy

force = dof.DOFtype(typestr='force')
disp  = dof.DOFtype(typestr='displacement')


ID   = np.arange(1,3)
ldof = np.arange(1,4)

ID   = np.repeat(ID,3)
ldof = np.tile(ldof,2)

my_dof = dof.DOF(ID,ldof,disp)

my_dof = dof.DOF([1,2],[1,2,3],disp,repetition = True)

many_dof = dof.DOF([1,1,2,2],[1,2,3,1],[disp,disp,force,force])

my_part_dof = dof.DOF([1,2],[2],disp,repetition = True)

ix = my_dof.index(my_part_dof)

my_dof[ix]

ang_freq  = dof.DOFtype(typestr='angular frequency')

u_dof = dof.DOF([1,2],[1,2,3],disp,repetition=True)
f_dof = dof.DOF([1,2],[1,2,3],force,repetition=True)
xdata = mC.DataAxis([10.,20.,30], typestr = 'angular frequency')

data = np.arange(108).reshape((6,6,3))

freq_axis  = mC.DataAxis(np.arange(0.,2.,0.1),typestr = 'frequency')

freq_axis.data
freq_axis.type

omega = freq_axis.angular_frequency

p_dof = dof.DOF([1,2],[0,0],dof.DOFtype(typestr = 'pressure'))
ydata = np.array([np.sin(omega),np.cos(omega)])

sig1 = mC.Signal(freq_axis,ydata,p_dof)

sig1.plot(1)


u_dof = dof.DOF([1,2],[1,2,3],disp,repetition=True)
f_dof = dof.DOF([1,2],[1,2,3],force,repetition=True)
x_data = mC.DataAxis([10.,20.,30], typestr = 'angular frequency')

data = 40*np.random.random_sample((6,6,3))

DD   = mC.DynamicMatrix(data, x_data, u_dof, f_dof)

DDinv = DD.inv()

f_data = np.zeros((6,3))
f_data[1,:] = 1.
force_load = mC.Signal(x_data,f_data, f_dof)

u_res = DDinv.dot(force_load)

f_data = np.ones((1,3))
force_load = mC.Signal(x_data,f_data, f_dof[0])

u_res = DDinv.dot(force_load)












