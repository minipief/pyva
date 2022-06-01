# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 20:51:13 2021

"""

from pyva.data import matrixClasses as mC
import numpy as np
import random as rnd

data = np.arange(18).reshape(3,3,2)
data[:,:,0]

lin_data = mC.LinearMatrix(data) 

sym_data = (data + data.transpose(1,0,2))/2
sym_data[:,:,0]

lin_sym_data = mC.LinearMatrix(sym_data)

one = np.eye(3)
mC.LinearMatrix(one)

triu_data = np.arange(6).reshape(6,1)
mC.LinearMatrix(triu_data, sym = 1, shape = (3,3,1))

lin_sym_data[0:2,0:2,0]
lin_sym_data[0:2,1:3,0]

lin_data.cond()










