# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 19:20:00 2018

"""

from pyva.data import matrixClasses as matC
import numpy as np

N = 4

for irow in range(N):
    for icol in range(irow,N):
        print(irow,icol,matC.linearIndex(N,irow,icol))

# itr2 fits perfectly for seclecti uper right triangle        
itr1,itr2 = np.triu_indices(N,1)
itr       = np.triu_indices(N,1)

# create simpe 4x4 matrix
test = np.arange(16).reshape(4,4)

# make it flat and add singleton dimension
test1   = test.reshape(N*N,1)
# take the upper triangle and reshape it to linear and transpose
iutest2 = test[itr].reshape(1,((N-1)*N)//2)


