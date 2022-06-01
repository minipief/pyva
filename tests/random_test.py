# -*- coding: utf-8 -*-
"""

My random experiment

Created on Wed Apr  5 09:44:50 2017
"""
import numpy as np
import matplotlib.pyplot as plt
import random as rnd

plt.close('all')

N  = 10000
f  = 200
Ni = 20
om = 2*np.pi*f
A  = 3

dt = 1/f/40

t  = np.arange(0.,Ni/f,dt)
Nt = len(t)

# random phase
def phasernd():
    return rnd.uniform(0,2*np.pi) 

def amprnd(Amax):
    return rnd.uniform(0,Amax)



    
Fmean = np.zeros(Nt,dtype=np.complex128)
FmsRe = np.zeros(Nt,dtype=np.complex128)
FmsCo = np.zeros(Nt,dtype=np.complex128)

for i in np.arange(1,N):
    sigRe = A*np.real(np.exp(1j*(om*t+phasernd())))
    sigCo = A*np.exp(1j*(om*t+phasernd()))
    #print(phasernd())
    #Fms   += sig**2/N # *random.uniform(1.-vc,1.+vc)  
    FmsRe   += sigRe*sigRe/N # *random.uniform(1.-vc,1.+vc)  
    FmsCo   += 0.5*np.real(sigCo*np.conj(sigCo))/N # *random.uniform(1.-vc,1.+vc)  
    Fmean += sigRe/N
    

plt.figure(1)
plt.plot(t,sigRe*np.conj(sig),'g',label='$A^2$')
plt.plot(t,FmsRe,'r-',label='mean square from Re')
plt.plot(t,FmsCo,'r:',label='mean square from Co')
plt.plot(t,Fmean,'b',label='mean')
plt.legend()
plt.xlabel('t')





