# Example for harmonic oscilator

import numpy as np

import pyva.systems.lumpedSystems as lSys
import matplotlib.pyplot as plt

# evenly sampled time at 200ms intervals
time = np.arange(0., 0.5, 0.001)

plt.close('all')

# First example for oscillator
mass   = 0.1
ks     = 98.696 # spring stiffness 

# Undamped HO
myHO = lSys.HarmonicOscillator(mass,ks)

# Derive all other constants from this
c_vc   = myHO.critical_viscous_damping

print('Critical damping from object is:'+str(myHO.critical_viscous_damping ))

# Initial conditions
x0 = 0.1
v0 = 1.4

# %% plot 1

plt.figure(1)
plt.plot(time, myHO.displacement(time,x0,v0),lw=2)
plt.xlabel('$t/$s')
plt.ylabel('$x$')
plt.ylim(-0.12,0.12)

#plt.savefig('../source/images/HO_oscillations.png')


# Damped HOs
cv1  = c_vc*3   # overdamped 
cv2  = c_vc/10  # underdamped

myHO_uD = lSys.HarmonicOscillator(mass,ks,cv2)
myHO_oD = lSys.HarmonicOscillator(mass,ks,cv1)
myHO_cD = lSys.HarmonicOscillator(mass,ks,c_vc)

# %% plot 2

plt.figure(2)
plt.plot(time, myHO_uD.displacement(time,x0,v0),lw=2,label = 'underdamped')
plt.plot(time, myHO_oD.displacement(time,x0,v0),lw=2,label = 'overdamped')
plt.plot(time, myHO_cD.displacement(time,x0,v0),lw=2,label = 'critically damped')
plt.xlabel('$t/$s')
plt.ylabel('$x$')
plt.ylim(-0.12,0.12)
plt.legend()

#plt.savefig('../source/images/HO_damped_oscillations.png')


# %%  plot 3

force = 10.0
omega = np.linspace(0,4*myHO.omega_mode,200)

plt.figure(3)
plt.plot(omega, np.abs(myHO_uD.u_force(omega,force)),lw=2,label = 'underdamped')
plt.plot(omega, np.abs(myHO_oD.u_force(omega,force)),lw=2,label = 'overdamped')
plt.plot(omega, np.abs(myHO_cD.u_force(omega,force)),lw=2,label = 'critically damped')
plt.ylabel('amplitude')
plt.xlabel('$\omega/$s$^{-1}$')
#plt.ylim(-0.12,0.12)
plt.legend()

#plt.savefig('../source/images/HO_forced_oscillations.png')

