# -*- coding: utf-8 -*-
"""
The load module deals with dynamic loads, currently only there is only
an extension of the Signal class implemeted with an extra artribute for 
the name of the load case
"""

import pyva.data.matrixClasses as mC

class Load(mC.Signal):
    """ 
    Class Load management in dynamic model extends Signal
    
    Attributes
    ----------
    xdata: DataAxis
        x-axis samples 
    ydata: ndarray or function 
        array of Signal samples with dimension Ndof x Nxdata with ndarray output of (Ndofs,Nxdata) output
    dof: DOF 
        object with ID, local DOF and DOFtype of ydata channels
    name : str
        identifier for loadcase
    """
    
    def __init__(self,xdata, ydata, dof, name = 'Load', **kwargs):
        """
        Constructor of Load 

        Parameters
        ----------
        xdata: DataAxis
            x-axis samples 
        ydata: ndarray or function 
            array of Signal samples with dimension Ndof x Nxdata with ndarray output of (Ndofs,Nxdata) output
        dof: DOF 
            object with ID, local DOF and DOFtype of ydata channels
        name : str
            identifier for loadcase. The default is 'Load'.
        **kwargs : dict
            Arbitrary keyword list.

        Returns
        -------
        None.

        """
        
        super().__init__(xdata,ydata,dof)
        
        self.name = name
        
        
