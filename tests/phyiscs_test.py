# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 21:34:09 2018

"""

from pyva.properties.physics import Q 


nix = 4.
print(nix)

dist = 200*ureg.meter
print(dist)

time = 8*ureg.second
print(time)

velocity = dist/time
print(velocity)
acceleration = dist/time**2


print(dist.dimensionality)
print(velocity.dimensionality)
print(acceleration.dimensionality)

#velocity.to_reduced_units()
print(velocity.to_base_units())

P = 20*ureg.pascal
area = 2*ureg.meter**2

F = P*area
print(F)
print(F.to_base_units())


