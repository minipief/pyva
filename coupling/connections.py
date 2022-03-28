"""
This module deals with the topology of connections

The connections define the geometry of the coupling, whereas the junction 
classes deals with the physics and the degrees of freedom/wavefields of the 
junction.

.. moduleauthor:: Alexander Peiffer <alexander.peiffer@pyva.eu>
"""

import numpy as np



# Debug switches used for detailed output
debug_sw = 0 # 2 global 3 amplitude 4 force_in 5 Pow_in 6 Pow_out
 
class Connection:
    
    """ 
    Connections define the contact between systems. 

    Attributes
    ----------
    geometries: list or tuple
        geometry describing the topology of the coupling

    """
    
    def __init__(self,geometries):
        """
        Class contructor for connection
        
        Parameters
        ----------
        geometries : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        self.geometries = geometries
