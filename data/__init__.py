# -*- coding: utf-8 -*-
"""
Subpackage for matrix, vector and degree of freedom handling

The data package include modules for handling degrees of freedom (DOF) and useful 
data and matrix classes that are required for vibroacoustic simuilation.

The DOF classes allow for handling the different physical quantities in vibroacoustics,
for example pressure, displacment, energy etc. as far as orientations, for example the
three room directions.
 
"""

__all__ = ['dof.py','matrixClasses.py']