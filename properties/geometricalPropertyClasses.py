# -*- coding: utf-8 -*-
"""
Module for handling of pure geometrical properties. 
"""

class CrossSection:
    """
    Class for beam cross section geometrical properties
    
    Attributes
    ----------
    Ix : float
        second moment of area of x-axis.
    Iy : float
        second moment of area of y-axis..
    Ixy : float
        product moment of area.
    area : float
        area of cross sectopm.    
    
    
    """    
    
    def __init__(self,Ix,Iy,Ixy,area):
        """
        Class constructor for beam crosSection

        Parameters
        ----------
        Ix : float
            second moment of area of x-axis.
        Iy : float
            second moment of area of y-axis..
        Ixy : float
            product moment of area.
        area : float
            area of cross sectopm.

        Returns
        -------
        None.

        """
        self.Ix = Ix
        self.Iy = Iy
        self.Ixy = Ixy
        self.area = area
        
    def __str__(self):
        """
        String output of CrossSection

        Returns
        -------
        str

        """
        _str  = "CrossSection: \n"
        _str += "Ix              : {0}\n".format(self.Ix)
        _str += "Iy              : {0}\n".format(self.Iy)
        _str += "Ixy             : {0}\n".format(self.Ixy)
        _str += "area            : {0}\n".format(self.area)
        
        return _str
            
    def __repr__(self):
        """
        Repr of CrossSection

        Returns
        -------
        None.

        """
        
        return "CrossSection({0},{1},{2},{3})".format(self.Ix,self.Iy,self.Ixy,self.area)
        
class RectBeam(CrossSection):
    """
    Daughter class of CrossSection for rectangular beams
    """
    
    def __init__(self,Lx,Ly):
        """
        Constructor for specific rectangular beam geometry

        Parameters
        ----------
        Lx : float
            length in x-direction.
        Ly : float
            lmegth in y-direction.

        Returns
        -------
        None.

        """
        
        area = Lx*Ly
        Ix   = Ly**3*Lx/12
        Iy   = Lx**3*Ly/12
        Ixy  = 0.
        
        self.Lx = Lx
        self.Ly = Ly
         
        
        super().__init__(Ix,Iy,Ixy,area)

    def __str__(self):
        """
        String output of CrossSection

        Returns
        -------
        str

        """
        _str  = "RectBeam: \n"
        _str += "Ix              : {0}\n".format(self.Ix)
        _str += "Iy              : {0}\n".format(self.Iy)
        _str += "Ixy             : {0}\n".format(self.Ixy)
        _str += "area            : {0}\n".format(self.area)
        _str += "Lx              : {0}\n".format(self.Lx)
        _str += "Ly              : {0}\n".format(self.Ly)
                
        return _str
            
    def __repr__(self):
        """
        Reps of acoutic tube

        Returns
        -------
        None.

        """
        
        return "RectBeam({0},{1})".format(self.Lx,self.Ly)
       
    
    
    
    
    
