# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 23:45:03 2018

"""

import pyva.data.dof as dof
import pyva.data.matrixClasses as mC
import numpy as np

gendof =  dof.DOFtype(typeID=1)
distdof = dof.DOFtype(typestr='displacement')

# errordof = dof.DOFtype(typestr='jhdfkh')

print(gendof)
print(distdof)

autopressure    = dof.DOFtype(typestr='pressure', exponent=2)
pressuredensity = dof.DOFtype(typestr='pressure', xtypestr='frequency', xdata_exponent=0.5)
psd             = dof.DOFtype(typestr='pressure', xtypeID=18 ,exponent=2, xdata_exponent=1)

listoftypestr   = ('pressure','force','pressure')
tistoftypeID    = (24,1,1)

#listoftypes     = dof.DOFtype(typestr=listoftypestr)
dofs3           = (3,3,3)
IDs3            = [1,2,3]


print(psd)
print(pressuredensity)

griddof  = dof.DOF(np.array([1,2,3,4]),np.array([1,2,3]),dof.DOFtype(typestr='displacement'),repetition=True)
vdof  = dof.DOF(np.array([1,2,3,4]),np.array([1,2,3]),dof.DOFtype(typestr='velocity'),repetition=True)
ogriddof = dof.DOF(np.array([2,3,5]),np.array([3,1]),dof.DOFtype(typestr='displacement'),repetition=True)

ogriddof.unique_type

print(griddof.typeID)


#sdof     = dof.DOF(1,3,dof.DOFtype(typestr='displacement'))

specdof = dof.DOF(IDs3,dofs3,listoftypestr)

partgriddof = griddof[3]

interdofs,_,_ = griddof.intersect(ogriddof)
print('Intersection-->')
print(interdofs)


# Show the details
print(griddof)
print(griddof.ID)
print(griddof.dof)
print(griddof.fulltype)


print(specdof)
print(specdof.ID)
print(specdof.dof)
print(specdof.fulltype)

print('Index of ogriddof into griddof')
print(griddof.index(ogriddof))


ht = griddof.hashtupel()
print(ht)
print(ht.sort())

# create full dof set



# DOFtype test

FF = dof.DOFtype(typestr='force')
AA = dof.DOFtype(typestr='area')
PP = dof.DOFtype(typestr='pressure')

FF2 = PP*AA
PP2 = PP*PP

FR = dof.DOFtype(typestr='frequency')
TT = dof.DOFtype(typestr='time')

PP2/FR
PP2/FR*FR

FplusF = FF+PP*AA

try:
    FF + PP
except ValueError:
    print("Add with wrong DOFtype doesn't work")
    
PP2 = PP**2
NIX = FR*TT

try:
    PP**TT
except ValueError:
    print("Pow with wrong argument doesn't work")
    
timeX = mC.DataAxis(np.arange(0,1,0.001),typestr='time')
print(timeX)

freqX = mC.DataAxis(np.arange(0,1,0.001),typestr='frequency')
print(freqX)
        
sinY   = 2*np.sin(2*np.pi*20*timeX.data)
sinY2  = 1.6*np.sin(2*np.pi*15*timeX.data)
sinfun = lambda timeX: 2.2*np.sin(2*np.pi*20*timeX)

s1 = mC.Signal(timeX,sinY,dof.DOF(1,0,['pressure']))
s1.plot(1)

#s2 = mC.Signal(timeX,np.stack((sinY,sinY2)),dof = dof.DOF([1,2],[0,0],dof.DOFtype(typestr='pressure')))
s2 = mC.Signal(timeX,np.stack((sinY,sinY2)),dof = dof.DOF([1,2],[0,0],['pressure']))
s2.plot(2)

s3 = s1[0]
s3.plot(1)


# test functoin as ydata argumenr - not yet implemented
#s3 = mC.Signal(timeX,sinfun)
#s1.plot(3)

# create full dof set



dofset = []






