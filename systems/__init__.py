# -*- coding: utf-8 -*-
"""
Subpackage provides for dynamic systems

Globally, this are canonic systems which can be modelled using
analytical formulations. For example rectangular rooms, pipes, rectangular plates etc.

In addition the aproximattve and random description of SEA systems is also part of this 
package.  

"""

__all__ = ["acoustic1Dsystems.py", "acoustic3Dsystems.py",\
           "acousticRadiators.py", "lumpedSystems.py", \
           "structure1Dsystems.py","structure2Dsystems.py", \
           "SEAsystem.py"]