# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 20:08:09 2018

"""

import numpy as np

# Generate some example data...
x = np.array((4,3,2,3,2,1,1,2,3,4,1))
#np.random.shuffle(x)
y = np.array((4,2,1,3))

# Actually preform the operation...
xsorted = np.argsort(y)
xpos = np.searchsorted(y[xsorted], x)
indices = xsorted[xpos]