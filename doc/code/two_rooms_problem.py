# -*- coding: utf-8 -*-
"""
Example script for lecture example
for Ulf Orenius
"""
import numpy as np
import matplotlib.pyplot as plt

# pyva packages
import pyva.models as mds
import pyva.coupling.junctions as con

import pyva.systems.structure2Dsystems as st2Dsys
import pyva.systems.acoustic3Dsystems as ac3Dsys
import pyva.loads.loadCase as lC


import pyva.properties.structuralPropertyClasses as stPC
import pyva.properties.materialClasses as matC

import pyva.data.dof as dof
import pyva.data.matrixClasses as mC

import pyva.useful as uf

mycolors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#17becf', '#17becf']
 
# x-axis tics
fc,fclabels = uf.get_3rd_oct_axis_labels(range='SEA')
# Print option
printsw = True
printsw = False
figpath = '../source/images/'
plt.close('all')

# wall dimensions
Ly = 4.
Lz = 2.5
S = Lz*Ly

# Additional room dimensions
Lx1 = 3.
Lx2 = 5.
# Absorption area
As1  = 8.
As2  = 10. # As2 = S

# wall thickness
h = 0.05

# junction properties
angle_S = 0
angle_R = 90*np.pi/180

# Create materials
concrete = matC.IsoMat(E=3.8e9,nu=0.33,rho0=1250.)
air = matC.Fluid()

# Create properties
concrete_5cm = stPC.PlateProp(h,concrete)

# Create plate subsystems
wall  = st2Dsys.RectangularPlate(2, Ly,Lz,prop=concrete_5cm, eta = 0.03)
room1 = ac3Dsys.RectangularRoom(1, Lx1, Ly, Lz, air, absorption_area = As1, damping_type= ['surface'] )
# use receiving room with absorption
room2 = ac3Dsys.RectangularRoom(3, Lx2, Ly, Lz, air, absorption_area = As2, damping_type= ['surface'])

J123 = con.AreaJunction((room1,wall,room2))#,area_wall)
#om = 2*np.pi*np.logspace(2,np.log10(10000),3)


# Frequency range
omega = mC.DataAxis.octave_band(f_max=2*np.pi*10000)
om    = omega.data
freq  = om/2/np.pi 

# define load
# room power
power1mWatt = lC.Load(omega, 0.001*np.ones(omega.shape), dof.DOF(1,0,dof.DOFtype(typestr = 'power')), name = '1Watt')

# wall power
#power1Watt = lC.Load(omega, np.ones(omega.shape), dof.DOF(2,3,dof.DOFtype(typestr = 'power')), name = '1Watt')

#create SEA model
two_rooms = mds.HybridModel((wall,room1,room2),xdata=omega)
#connect both
two_rooms.add_junction({'areaJ_12':J123})
two_rooms.add_load('1mWatt',power1mWatt) # add 1Watt per band 

#pdof = J123.get_wave_DOF()

# Calculate specific CLF factors
eta_13 = J123.CLF(om,(0,2))
eta_31 = J123.CLF(om,(2,0))
eta_12 = J123.CLF(om,(0,1),)
eta_21 = J123.CLF(om,(1,0))

eta_all = J123.junction_matrix(om)

# check tau
tau_diff = J123.non_res.transmission_diffuse(om,signal = False)



#%% plot the different modal densities
plt.figure(1)
plt.plot(freq,wall.modal_density(om,3),label = 'wall')
plt.plot(freq,room1.modal_density(om),label = 'room1')
plt.plot(freq,room2.modal_density(om),label = 'room2')
# plt.loglog(f,n_1,':',label = 'n1 VA1')
# plt.loglog(f,n_2,':',label = 'n2 VA1')
# plt.loglog(f,n_3,':',label = 'n3 VA1')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$f_c/$Hz')
plt.ylabel('$n(\omega)/$Hz$^{-1}$')
plt.legend()
plt.show

if printsw:
     plt.savefig(figpath+'two_rooms_modal_density.png')


#%% plot modal overlap 
plt.figure(2)
plt.plot(freq,wall.modal_overlap(om,3),label = 'wall B')
plt.plot(freq,wall.modal_overlap(om,5),label = 'wall LS')
plt.plot(freq,room1.modal_overlap(om),label = 'room1')
plt.plot(freq,room2.modal_overlap(om),label = 'room2')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$f_c/$Hz')
plt.ylabel('$M$')
plt.ylim(bottom=0.1)
plt.xticks(fc[2:],fclabels[2:])
plt.legend()
plt.show

if printsw:
    plt.savefig(figpath+'two_rooms_modal_overlap.png')

#%% plot modes in band 
plt.figure(3)
plt.plot(freq,wall.modes_in_band(om,3),label = 'wall B')
plt.plot(freq,wall.modes_in_band(om,5),label = 'wall LS')
plt.plot(freq,room1.modes_in_band(om),label = 'room1')
plt.plot(freq,room2.modes_in_band(om),label = 'room2')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$f_c/$Hz')
plt.ylabel('modes in band')
plt.ylim(bottom=0.1)
plt.xticks(fc[2:],fclabels[2:])
plt.legend()
plt.show

if printsw:
    plt.savefig(figpath+'two_rooms_modes_in_band.png')


 



#%% PLot radation efficiency
sigma = wall.radiation_efficiency(om,fluid = air)
sigma_simple = wall.radiation_efficiency_simple(om,fluid = air)

# Show coincidence frequency
f_c = concrete_5cm.coincidence_frequency(air.c0)/2/np.pi
print('coincidence frequency = {0}'.format(f_c))

plt.figure(4)
plt.loglog(freq,sigma,label = '[Lep1982]')
plt.loglog(freq,sigma_simple,label = 'ISO EN 12354-1')
#plt.loglog(om_VA1,sigma_VA1,label = 'VAOne')
plt.xlabel('$f_c/$Hz')
plt.ylabel('$\sigma$')
#plt.ylim(0,0.005)
#plt.xticks(2*np.pi*fc,fclabels)
plt.legend()
plt.show

if printsw:
    plt.savefig(figpath+'two_rooms_radiation_efficiency.png')


#%% plot4
def mass_law(omega,theta=0):
    z0 = air.z0
    m  = wall.prop.mass_per_area
    return 1/(1+(m*omega/2/z0)**2)

def eta_f(omega,V,area):
    c0 = air.c0
    return c0*area/4/V/omega*mass_law(omega)

plt.figure(5)
eta_13.plot(5,yscale = 'log',fulllegstr = ['$\eta_{13}$'],ls='-')
eta_31.plot(5,yscale = 'log',fulllegstr = ['$\eta_{31}$'],cs=[mycolors[1]],ls = '--')
#plt.loglog(om_VA1,eta_13_VA1,':',label = '13 VA1')
#plt.loglog(om_VA1,eta_31_VA1,':',label = '31 VA1')
# #%% plot5
#eta_12.plot(4,yscale = 'log',fulllegstr = ['$\eta_{1,2B}$'],cs=[mycolors[2]],ls = ':',marker = [' '])
#eta_21.plot(4,yscale = 'log',fulllegstr = ['$\eta_{2B,1}$'],cs=[mycolors[3]],ls = '-.',marker = [' '])
#plt.loglog(om_VA1,eta_12_VA1,':',label = '13 VA1')
#plt.loglog(om_VA1,eta_21_VA1,':',label = '31 VA1')
plt.xlabel('$f_c/$Hz')
plt.ylabel('$\eta$')
#plt.ylim(0,0.005)
ax = plt.gca()
ax.axes.xaxis.set_ticks([])
plt.xscale('log')

plt.xticks(2*np.pi*fc[2:],fclabels[2:])
plt.legend()
plt.show

if printsw:
    plt.savefig(figpath+'two_rooms_etas.png')

#%% Solve model
two_rooms.create_SEA_matrix(sym = 1)
two_rooms.solve()

# Derive paths of power input
pow_in_room2 = two_rooms.power_input(3)

#%% Enrgry result plots
two_rooms.energy.plot(20,xscale = 'log',yscale = 'log',ls = ['-','--',':','-.'],
                       fulllegstr = ('wall B','wall LS','room1','room2'))
plt.xlabel('$f_c/$Hz')
plt.ylabel('$E/$ J')
#plt.ylim(bottom=0.1)
plt.xticks(2*np.pi*fc[2:],fclabels[2:])
plt.legend()
plt.show()

if printsw:
    plt.savefig(figpath+'two_rooms_energy.png')



#%% Engineering units velocity
plt.figure(30,figsize=(5,5/4*3))
two_rooms.result.plot(30,ID=[2],dof=[3],xscale = 'log',yscale = 'log',
                       fulllegstr = ('wall B',))
#plt.plot(om_VA1,v_wall_B,':',label = 'VA1')
plt.xlabel('$f_c/$Hz')
plt.ylabel('$v_{\\rm rms}/$m s$^{-1}$')
#plt.ylim(bottom=0.1)
ax = plt.gca()
ax.axes.xaxis.set_ticks([])

plt.xticks(2*np.pi*fc[2:],fclabels[2:])
plt.legend()
plt.tight_layout()  	
plt.show()

if printsw:
    plt.savefig(figpath+'two_rooms_velocity.png')

#%% engineering units pressure
plt.figure(31,figsize=(5,5/4*3))

two_rooms.result.plot(31,ID=[1,3],xscale = 'log',yscale = 'log',ls = ['-','--',':','-.'],
                       fulllegstr = ('room 1','room 2'))
#plt.plot(om_VA1,p_room1,':',label = 'room1 VA1')
#plt.plot(om_VA1,p_room2,':',label = 'room2 VA1')
plt.xlabel('$f_c/$Hz')
plt.ylabel('$p_{\\rm rms}/$Pa')
#plt.ylim(bottom=0.1)
ax = plt.gca()
ax.axes.xaxis.set_ticks([])

plt.xticks(2*np.pi*fc[2:],fclabels[2:])
plt.legend()
plt.tight_layout()  	
plt.show()

if printsw:
    plt.savefig(figpath+'two_rooms_pressure.png')
    
# Calculate transmission loss from pressure difference
# Note: This is allowed because wall and absorption surface are equal S=A2
p1 = two_rooms.result[2].ydata.flatten()
p2 = two_rooms.result[3].ydata.flatten()
tau = (p2/p1)**2

#%% Plot transmission loss 

plt.figure(32)
plt.semilogx(freq,-10*np.log10(tau),label = 'wall TL')
#plt.semilogx(f,20*np.log10(tau_VA1),':',label = 'pure VA1')
plt.semilogx(freq,-10*np.log10(mass_law(om)),':',label = 'mass law')
#plt.semilogx(freq,20*np.log10(wall.prop.mass_per_area*freq)-48,':',label = 'mass law diffuse')
#plt.semilogx(freq,-10*np.log10(tau_diff),':',label = 'mass law TMM')
plt.xticks(fc[2:],fclabels[2:])
plt.xlabel('$f_c/$Hz')
plt.ylabel('$TL/$dB')
plt.legend()
plt.show()

if printsw:
    plt.savefig(figpath+'two_rooms_TL.png')

#%% Plot input power
    
pow_in_room2.plot(33,yscale = 'log',xscale = 'log')
#plt.loglog(om_VA1,pow_wall,':')
#plt.loglog(om_VA1,pow_room1,':')
plt.xticks(2*np.pi*fc,fclabels)
plt.xlabel('$f_c/$Hz')
plt.ylabel('$\Pi_{in}/$W')
plt.legend()
plt.show()

if printsw:
    plt.savefig(figpath+'two_rooms_power_in.png')






    





