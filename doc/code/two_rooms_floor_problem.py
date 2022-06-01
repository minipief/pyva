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
fc,fclabels = uf.get_3rd_oct_axis_labels(range='hybrid')
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
h   = 0.05
h_f = 0.17 

# junction properties
angle_S = 0
angle_R = 90*np.pi/180

# Create materials
concrete = matC.IsoMat(E=3.8e9,nu=0.33,rho0=1250.)
air = matC.Fluid()

# Create properties
concrete_5cm   = stPC.PlateProp(h,concrete)
concrete_17cm  = stPC.PlateProp(h_f,concrete)

# Create plate subsystems
wall   = st2Dsys.RectangularPlate(2, Ly,Lz,prop=concrete_5cm, eta = 0.03)
floor1 = st2Dsys.RectangularPlate(4, Lx1,Ly,prop=concrete_17cm, eta = 0.03)
floor2 = st2Dsys.RectangularPlate(5, Lx2,Ly,prop=concrete_17cm, eta = 0.03)
room1 = ac3Dsys.RectangularRoom(1, Lx1, Ly, Lz, air, absorption_area = As1, damping_type= ['surface'] )
# use receiving room with absorption
room2 = ac3Dsys.RectangularRoom(3, Lx2, Ly, Lz, air, absorption_area = As2, damping_type= ['surface'])

# Area Junctions
J123 = con.AreaJunction((room1,wall,room2))
J14  = con.AreaJunction((room1,floor1))
J35  = con.AreaJunction((room2,floor2))

# Line Junctions
J425 = con.LineJunction((floor1,wall,floor2),length = Ly, thetas = (0,90,180))


# Frequency range
omega = mC.DataAxis.octave_band(f_min=2*np.pi*50,f_max=2*np.pi*10000)
om    = omega.data
freq  = om/2/np.pi 

# define load
# room power
# power1mWatt = lC.Load(omega, 0.001*np.ones(omega.shape), dof.DOF(1,0,dof.DOFtype(typestr = 'power')), name = '1Watt')
force10Nrms = lC.Load(omega, 10*np.ones(omega.shape), dof.DOF(4,3,dof.DOFtype(typestr = 'force')), name = '10N')

# wall power
#power1Watt = lC.Load(omega, np.ones(omega.shape), dof.DOF(2,3,dof.DOFtype(typestr = 'power')), name = '1Watt')

#create SEA model
two_rooms = mds.HybridModel((wall,room1,room2,floor1,floor2),xdata=omega)
#connect all
two_rooms.add_junction({'areaJ_123':J123})
two_rooms.add_junction({'areaJ_14':J14})
two_rooms.add_junction({'areaJ_35':J35})
two_rooms.add_junction({'lineJ_425':J425})

#two_rooms.add_load('10Watt',power1mWatt) # add 1Watt per band 
two_rooms.add_load('10N',force10Nrms)# add force excitatio to wave_DOF 3 of system 4
#
# Calculate specific CLF factors
eta_13 = J123.CLF(om,(0,2))
eta_31 = J123.CLF(om,(2,0))
eta_12 = J123.CLF(om,(0,1),)
eta_21 = J123.CLF(om,(1,0))

eta_all = J123.junction_matrix(om)

# check tau
tau_diff = J123.non_res.transmission_diffuse(om,signal = False)

# read reference results (keep if ulf will check with VA1)
#f,n_1,n_2,n_3 = np.loadtxt ('2rooms_m_dens.csv',
#                    unpack = True,
#                    usecols = (0,2,1,5), skiprows = 1,
#                    delimiter = ',')

# f,sigma_VA1 = np.loadtxt ('2rooms_sigma.csv',
#                     unpack = True,
#                     usecols = (0,1), skiprows = 1,
#                     delimiter = ',')

# f,p_room1,p_room2 = np.loadtxt ('2rooms_p_rms.csv',
#                     unpack = True,
#                     usecols = (0,1,2), skiprows = 1,
#                     delimiter = ',')

# f,v_wall_B = np.loadtxt ('2rooms_v_rms.csv',
#                     unpack = True,
#                     usecols = (0,1), skiprows = 1,
#                     delimiter = ',')

# f,pow_wall,pow_room1 = np.loadtxt ('2rooms_pow_in.csv',
#                     unpack = True,
#                     usecols = (0,3,2), skiprows = 1,
#                     delimiter = ',')

# f,tau_VA1,tau_eff_VA1 = np.loadtxt ('2rooms_tau.csv',
#                     unpack = True,
#                     usecols = (0,1,2), skiprows = 1,
#                     delimiter = ',')

# f,eta_13_VA1,eta_31_VA1,eta_12_VA1,eta_21_VA1,eta_23_VA1,eta_32_VA1 = np.loadtxt ('2rooms_eta.csv',
#                     unpack = True,
#                     usecols = (0,1,2,6,5,3,4), skiprows = 1,
#                     delimiter = ',')


# om_VA1 = f*2*np.pi


#%% plot the different modal densities
plt.figure(1)
plt.plot(freq,wall.modal_density(om,3),label = 'wall')
plt.plot(freq,floor1.modal_density(om,3),label = 'floor1')
plt.plot(freq,floor2.modal_density(om,3),label = 'floor2')
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
     plt.savefig(figpath+'two_rooms_floor_modal_density.png')


#%% plot modal overlap 
plt.figure(2)
plt.plot(freq,wall.modal_overlap(om,3),label = 'wall B')
#plt.plot(freq,wall.modal_overlap(om,5),label = 'wall LS')
plt.plot(freq,floor1.modal_overlap(om,3),label = 'floor1')
plt.plot(freq,floor2.modal_overlap(om,3),label = 'floor2')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$f_c/$Hz')
plt.ylabel('$M$')
plt.ylim(bottom=0.1)
plt.xticks(fc[2:],fclabels[2:])
plt.legend()
plt.show

if printsw:
    plt.savefig(figpath+'two_rooms_floor_modal_overlap.png')

#%% plot modes in band 
plt.figure(3)
plt.plot(freq,wall.modes_in_band(om,3),label = 'wall B')
#plt.plot(freq,wall.modes_in_band(om,5),label = 'wall LS')
plt.plot(freq,floor1.modes_in_band(om,3),label = 'floor1')
plt.plot(freq,floor2.modes_in_band(om,3),label = 'floor2')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$f_c/$Hz')
plt.ylabel('modes in band')
plt.ylim(bottom=0.1)
plt.xticks(fc[2:],fclabels[2:])
plt.legend()
plt.show

if printsw:
    plt.savefig(figpath+'two_rooms_floor_modes_in_band.png')


 



#%% PLot radation efficiency
sigma1 = floor1.radiation_efficiency(om,fluid = air)
sigma2 = floor2.radiation_efficiency(om,fluid = air)

# Show coincidence frequency
f_cw = concrete_5cm.coincidence_frequency(air.c0)/2/np.pi
f_cf = concrete_17cm.coincidence_frequency(air.c0)/2/np.pi
print('coincidence frequency wall = {0}'.format(f_cw))
print('coincidence frequency floor = {0}'.format(f_cf))

plt.figure(4)
plt.semilogx(freq,sigma1,label = 'floor1')
plt.semilogx(freq,sigma2,label = 'floor2')
#plt.loglog(om_VA1,sigma_VA1,label = 'VAOne')
plt.xlabel('$f_c/$Hz')
plt.ylabel('$\sigma$')
#plt.ylim(0,0.005)
#plt.xticks(2*np.pi*fc,fclabels)
plt.legend()
plt.show

if printsw:
    plt.savefig(figpath+'two_rooms_floor_radiation_efficiency.png')


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
    plt.savefig(figpath+'two_rooms_floor_etas.png')

#%% Solve model
two_rooms.create_SEA_matrix(sym = 1)
two_rooms.solve()

# Derive paths of power input
pow_in_room1 = two_rooms.power_input(1)
pow_in_room2 = two_rooms.power_input(3)


#%% Enrgry result plots
two_rooms.energy.plot(20,xscale = 'log',yscale = 'log',ls = ['-','--',':','-.'],
                        fulllegstr = ('wall B','wall LS','room1','room2','floor1 B','floor1 LS','floor2 B','floor2 LS'))
plt.xlabel('$f_c/$Hz')
plt.ylabel('$E/$ J')
#plt.ylim(bottom=0.1)
plt.xticks(2*np.pi*fc[2:],fclabels[2:])
plt.legend()
plt.show()

if printsw:
    plt.savefig(figpath+'two_rooms_floor_energy.png')



#%% Engineering units velocity
plt.figure(30,figsize=(5,5/4*3))
two_rooms.result.plot(30,ID=[2,4,5],dof=[3,3,3],xscale = 'log',yscale = 'log',
                       fulllegstr = ('wall B','floor1 B','floor2 B'))
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
    plt.savefig(figpath+'two_rooms_floor_velocity.png')

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
    plt.savefig(figpath+'two_rooms_floor_pressure.png')
    
# Calculate transmission loss from pressure difference
# Note: This is allowed because wall and absorption surface are equal S=A2
p1 = two_rooms.result[2].ydata.flatten()
p2 = two_rooms.result[3].ydata.flatten()
tau = (p2/p1)**2


#%% Plot input power
    
pow_in_room1.plot(32,yscale = 'log',xscale = 'log')
#plt.loglog(om_VA1,pow_wall,':')
#plt.loglog(om_VA1,pow_room1,':')
plt.xticks(2*np.pi*fc,fclabels)
plt.xlabel('$f_c/$Hz')
plt.ylabel('$\Pi_{in}/$W')
plt.legend()
plt.show()

if printsw:
    plt.savefig(figpath+'two_rooms_floor_power_in1.png')

pow_in_room2.plot(33,yscale = 'log',xscale = 'log')
#plt.loglog(om_VA1,pow_wall,':')
#plt.loglog(om_VA1,pow_room1,':')
plt.xticks(2*np.pi*fc,fclabels)
plt.xlabel('$f_c/$Hz')
plt.ylabel('$\Pi_{in}/$W')
plt.legend()
plt.show()

if printsw:
    plt.savefig(figpath+'two_rooms_floor_power_in2.png')





    





