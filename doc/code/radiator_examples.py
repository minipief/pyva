# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 09:04:17 2022

"""

import numpy as np
import matplotlib.pyplot as plt

import pyva.properties.materialClasses as matC
import pyva.systems.acousticRadiators as acR
import pyva.geometry.meshClasses as meshC

plt.close('all')

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = r'\usepackage{bm}'

#Define the fluid
air   = matC.Fluid(eta = 0.)

# Create HalfSpace without trim
HS = acR.HalfSpace(air)

omega = 2*np.pi*np.linspace(100,5000,100)
omega0 = 1000.*2*np.pi

# Generate wavemuber
kx_max = air.wavenumber(omega0)
kx = np.linspace(0,2*kx_max,200)

D_wavenumber = HS.radiation_stiffness_wavenumber(omega0, kx)


#%% plot radstiffness over kx

plt.figure(1)
plt.plot(kx,np.imag(D_wavenumber), label='$Im$')
plt.xlabel('$k_x/$m$^{-1}$')
plt.ylabel('${\\bm D}\' N/$m$^{-2}$')
plt.legend()

#plt.savefig('../source/images/half_space_stiffness_wavenumber.png')


# Radiator dimensions
Lx = 0.8
Ly = 0.5

# Number of half sine wave on rectangle
nx = 10
ny = 3
# Corresponding wavenumber
kx = np.pi*nx/Lx
ky = np.pi*ny/Ly

sigma_LEP = HS.radiation_efficiency_leppington(omega,kx,ky,Lx,Ly,simple_muGT1=True)


Nx = 40
Ny = 25
dA = Lx*Ly/(Nx-1)/(Ny-1)

mesh     = HS.get_mesh(omega[-1], Lx, Ly, N=4)

 
shapefun   = lambda x,y: np.sin(kx*x)*np.sin(ky*y)
my_shape    = meshC.RegShape2D(0,0,Lx,Ly,Nx,Ny,shape = shapefun)   
sigma_wavelet  = HS.shape_radiation_efficiency(omega,my_shape,'piston',signal = False)  






#%% plot sigma

plt.figure(2)
plt.plot(omega,sigma_LEP,label='Leppington')
plt.plot(omega,sigma_wavelet,label='piston')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('$\\sigma$')
#plt.ylim(-50,50 )
plt.legend()

#plt.savefig('../source/images/half_space_sigma.png')

dist = 0.4

dA = my_shape.dA

D12piston = HS.radiation_stiffness_piston(omega, dist, dA, dA )
D11piston = HS.radiation_stiffness_piston(omega, 0., dA, dA )

ks = my_shape.ks

D12wavelet = HS.radiation_stiffness_wavelet(omega, dist, ks)
D11wavelet = HS.radiation_stiffness_wavelet(omega, 0., ks)

#%% discrete stiffness D12

plt.figure(3)
plt.plot(omega,np.real(D12piston),label='Re piston')
plt.plot(omega,np.imag(D12piston),label='Im piston')
plt.plot(omega,np.real(D12wavelet),':',label='Re wavelet')
plt.plot(omega,np.imag(D12wavelet),':',label='Im wavelet')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('$D_{ij}/$N m$^{-1}$')
#plt.ylim(-50,50 )
plt.legend()

plt.savefig('../source/images/half_space_D12.png')

#%% discrete stiffness D11

plt.figure(4)
plt.plot(omega,np.real(D11piston),label='Re piston')
plt.plot(omega,np.imag(D11piston),label='Im piston')
plt.plot(omega,np.real(D11wavelet),':',label='Re wavelet')
plt.plot(omega,np.imag(D11wavelet),':',label='Im wavelet')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('$D_{ii}/$N m$^{-1}$')
#plt.ylim(-50,50 )
plt.legend()

plt.savefig('../source/images/half_space_D11.png')

D_1000 = HS.radiation_stiffness_mesh_single([omega0], mesh,method = 'wavelet')



#%% Piston

radius = 0.3
my_piston = acR.CircularPiston(radius,fluid=air)

rad_impedance = my_piston.acousticImpedance(omega)

#%% plot 5

plt.figure(5)
plt.plot(omega,np.real(rad_impedance),label='Re piston')
plt.plot(omega,np.imag(rad_impedance),label='Im piston')
plt.xlabel('$\omega/$s$^{-1}$')
plt.ylabel('${\\bm Z}_a$')
#plt.ylim(-50,50 )
plt.legend()

plt.savefig('../source/images/piston_acoustic_impedance.png')


    
    