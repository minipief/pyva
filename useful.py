# -*- coding: utf-8 -*-
"""
Collection of usefull functions, no classes. 
"""

import numpy as np
import pyva.data.matrixClasses as mC
import matplotlib.pyplot as plt
from cycler import cycler


from pathlib import Path
import sys
 
def set_vascript_graphics(bw = True,cycler_sw = True,form = 'lecture'):
    """
    

    Parameters
    ----------
    bw : boolean, optional
        Switch for black and white output. The default is True.
    cycler_sw : boolean, optional
        Swith for setting a new linestyle cyler. The default is True.
    form : str, optional
        'book' or 'lecture'. The default is 'book'.

    Returns
    -------
    figpath : TYPE
        DESCRIPTION.
    aprop : TYPE
        DESCRIPTION.
    aprop2 : TYPE
        DESCRIPTION.
    tlfz : TYPE
        DESCRIPTION.
    alfz : TYPE
        DESCRIPTION.
    lst : TYPE
        DESCRIPTION.
    pst : TYPE
        DESCRIPTION.
    lst2nd : TYPE
        DESCRIPTION.
    lst3rd : TYPE
        DESCRIPTION.

    """

    #plt.close('all')

    global figpath

    if sys.platform=='linux':
        figpath = '//home/alex//Documents//work//VA_script//'
    else:
        #figpath = 'C://Users//alex//Documents//VA_script//matlab//'
        #figpath = str(Path.home())+'//Documents//VA_script//'
        figpath = str(Path.home())+'//Documents//python_VA_lib//VAScript//scripts//figs//'

    # display constant
        
    #global tlfz,alfz  
    if form=='lecture':
        tlfz = 14
        alfz = 16
    
        #global aprop, aprop2  
        aprop  = dict(facecolor='black', shrink=0.01, width=1 )
        aprop2 = dict(arrowstyle="<->", facecolor='black' , lw=2 )
              
        font = {'family': 'sans-serif', 'sans-serif': ['DejaVu Sans'], 'size'   : 18}
                
    else:
        tlfz = 14
        alfz = 16
    
        #global aprop, aprop2  
        aprop  = dict(facecolor='black', shrink=0.01, width=1 )
        aprop2 = dict(arrowstyle="<->", facecolor='black' , lw=2 )
              
        font = {'family': 'sans-serif', 'sans-serif': ['DejaVu Sans'], 'size'   : 16}
                

    
    if bw:
        color_c = cycler('color', ['k'])
        style_c = cycler('linestyle', ['-', '--', ':', '-.'])
        markr_c = cycler('marker', ['', '.', 'o'])
        c_cms = color_c * markr_c * style_c
    else:
        # nice default colors
        mycolors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#17becf', '#17becf']
        color_c = cycler('color', mycolors)
        style_c = cycler('linestyle', ['-', '--', ':', '-.'])
        markr_c = cycler('marker', ['', '.', 'o'])
        
    if form == 'lecture':
        c_cms = style_c * color_c 
    else:
        c_ms = markr_c * style_c
        c_cms = color_c + c_ms
            
    plt.rc('font',**font)
    plt.rc('text', usetex=True)
    # Set parameter to use boldmath package bm
    plt.rcParams['text.latex.preamble'] = r'\usepackage{bm}'

    if cycler_sw:
        plt.rc('axes',prop_cycle = c_cms )
                
    #global lst,pst    
    lst = ('r-','b-','m-','g-','c-','k-')
    pst = ('d','o','x','<','>','+')
    lst2nd = ('r:','b:','m:','g:','c:','k:')
    lst3rd = ('r--','b--','m--','g--','c--','k--')

    return figpath,aprop,aprop2,tlfz,alfz,lst,pst,lst2nd,lst3rd


# helpful simple functions
def isscalar(data):
    if isinstance(data, np.ndarray) and data.size == 1:
        return True
    
    return np.isscalar(data) and not(isinstance(data,str))

def get_omega_values(omega):
    """
    provides angular frequency from DataAxis in any case 
    """
     
    if isinstance(omega, mC.DataAxis):
        if omega.type.type == 18:
            omega = omega.data*2*np.pi
        elif omega.type.type == 21:
            omega = omega.data
        else:
            raise ValueError('Omega must be an instance of DataAxis of type frequency')
    else:
        pass
    
    return omega

def get_omega_DataAxis(omega):
    """
    provides angular frequency xdata object in any case 
    """
    if isinstance(omega, mC.DataAxis):
        if omega.type.type == 18:
            omega = mC.DataAxis(omega.data*2*np.pi,typeID = 21)
        elif omega.type.type == 21:
            pass
        else:
            raise ValueError('Omega must be an instance of DataAxis of type frequency')
    else:
        omega = mC.DataAxis(omega,typeID = 21)
    
    return omega

def get_3rd_oct_axis_labels(range = 'SEA'):
    """
    sets typical frequency axis labels for 3rd oct
    """
    
    if range == 'SEA':
        fclabels = ['100','250','500','1k','2.5k','5k']
        fc       = np.array((100,250,500,1000,2500,5000))
    elif range == 'hybrid':
        fclabels = ['50','100','250','500','1k','2.5k']
        fc       = np.array((50,100,250,500,1000,2500))
        
        
    return (fc,fclabels)
    
def get_n_th_oct_frequencies(f_min,f_max,N=3):
    """
    Method for octave / log space creation.
 

    Parameters
    ----------
    f_min : float
        lower frequency limit
    f_max : float
        upper frequency limit
    N : integer, optional
        bands_per_octave. The default is 3.
 
    Returns
    -------
    fc,flim : ndarray
        center and limit frequencies - limit freqs have one more member
 
    """
 
    fa1 = 2**(1/N)
 
    pi_fac = 1
 
    f_min = f_min/pi_fac
    f_max = f_max/pi_fac
 
    # first index
    i0 = np.round(np.log10(f_min/1000)/np.log10(fa1)) # -3 from log10(1000)
    # number of lines
    NN  = np.round((np.log10(f_max)-3)/np.log10(fa1))-i0+1
    ix = np.arange(0,NN)
    freqc = pi_fac*1000*(fa1**(i0+ix)) 
    freql = freqc/np.sqrt(fa1)
 
    return freqc,np.append(freql,freqc[-1]*np.sqrt(fa1))



