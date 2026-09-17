# -*- coding: utf-8 -*-
"""
The pyva package is a framework for the simulation of vibroacoustic systems.
It comprises deterministic (FEM) and random (SEA) methods, as far as the combination 
of both the hybrid FEM/SEA method.

This includes classes ...

- for material properties of fluids, solids amd porous materials
- for physical properties of plates, shells and beams
- for systems and the coupling of them


"""

# from pint import UnitRegistry
# ureg = UnitRegistry()
# Q_   = ureg.Quantity

__all__ = ['models.py','useful.py']
