# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 09:04:17 2022

"""

import numpy as np
import matplotlib.pyplot as plt

import pyva.systems.structure1Dsystems as st1Dsys
import pyva.properties.materialClasses as matC
import pyva.properties.geometricalPropertyClasses as geoPC
import pyva.properties.structuralPropertyClasses as stPC



# Define the properties
alu = matC.IsoMat() # Alu is default

# Beam constants
# First example for oscillator
h    = 0.0005
b    = 0.0005
L    = 2

# Cross section
beam_section = geoPC.RectBeam(h, b)
# Beam property
beam_prop    = stPC.BeamProp(beam_section,alu)

beam = st1Dsys.beam(L,beam_prop)

