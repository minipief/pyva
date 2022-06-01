# -*- coding: utf-8 -*-
"""
Created on Sat Oct 10 16:55:57 2020
two plate script for book presentation 
"""
import numpy as np
import matplotlib.pyplot as plt


import pyva.coupling.junctions as jun
import pyva.properties.structuralPropertyClasses as stPC
import pyva.systems.structure2Dsystems as st2Dsys
import pyva.systems.acoustic3Dsystems as ac3Dsys
import pyva.loads.loadCase as lC
import pyva.models as mds
import pyva.useful as uf
import pyva.data.dof as dof
import pyva.data.matrixClasses as mC
import pyva.properties.materialClasses as matC

# %% Initialisation

plt.close('all')

# x-axis tics
fc,fclabels = uf.get_3rd_oct_axis_labels(range='hybrid')
#fc,fclabels = fc[:],fclabels[:]

# Plate dimensions
Lx = 0.8
Ly = 0.5
# Mesh constants
Nx = 80
Ny = 50

# Plate thickness
t = 4.0E-3;

# Create materials
alu = matC.IsoMat(nu=0.3,eta = 0.0)
air = matC.Fluid()
# Create props
alu4mm = stPC.PlateProp(t,alu)

#Frequencies
omega_max = 2*np.pi*4000
omega = mC.DataAxis(np.geomspace(2*np.pi*25,omega_max,150), typestr = 'angular frequency')
freq  = omega.frequency 

f_c = alu4mm.coincidence_frequency(air.c0)/2/np.pi
print('coincidence frequency = {0}'.format(f_c))

# Create plate subsystem
plate = st2Dsys.RectangularPlate(2,Lx,Ly,prop=alu4mm,wave_DOF = [3],eta = 0.02)
# Create plate as FE-Model
modes,mesh = plate.normal_modes(omega_max*1.2,mapping = 'mesh')
plateFE    = mds.FEM(2,mesh,modes,damping_loss=0.02)


# Create rooms
room1 = ac3Dsys.Acoustic3DSystem(1, 64 , 96, 48, air)
room2 = ac3Dsys.Acoustic3DSystem(3, 80 ,112, 52, air, absorption_area = Lx*Ly, damping_type= ['surface'])

# %% Model set-up
HJ123 = jun.HybridAreaJunction((room1,room2),plateFE)#,area_plate)

# define loads
# plate force
forceID = 199
force10N = lC.Load(omega, 10*np.sqrt(2)*np.ones(omega.shape), \
                          dof.DOF(forceID,3,dof.DOFtype(typestr = 'force')), name = '10N@Node'+str(forceID))
# check position
X,Y = mesh.nodes()
print('Excitation at X={0:.2f}, Y={1:.2f}'.format(X.flatten()[forceID],Y.flatten()[forceID]))

# Add load to FEmodel
plateFE.add_load(force10N) 
# Get modal forec for debugging
#FM      = plateFE.modal_force('10N@Node200')

# Calculate CLF for debugging
#eta, eta_alpha,power_in,x_modal = HJ123.CLF(omega.angular_frequency,Signal = False,force='10N@Node200')
eta, eta_alpha = HJ123.CLF(omega.angular_frequency) #,Signal = False,force='10N@Node200')
tau_from_eta = 8*np.pi**2*room1.modal_density(omega.data,0)*omega.data/room1.fluid.wavenumber(omega.data)**2/plate.area*eta
tau_from_eta = tau_from_eta.flatten() 


#print_sw = False
#mat_sw = False

#create hybrid SEA model
RPR_FE_force = mds.HybridModel((room1,room2),FEsystems = (plateFE,),xdata=omega)

#connect both
RPR_FE_force.add_hybrid_junction({'HareaJ_12':HJ123})

pdof = HJ123.get_wave_DOF()

# %% solve hybrid SEA models
RPR_FE_force.create_SEA_matrix(sym = 1,force = '10N@Node'+str(forceID))
RPR_FE_force.solve()


#%% exctract FEM resonse using junction methods
sqq_type = dof.DOFtype(typestr='displacement',exponent = 2)
# Determine CSPD of FEM system
Sqq_F = HJ123.FEM_response(omega.angular_frequency , RPR_FE_force.energy)
# Determine nodal average from modal response
x2rms_F,v_type = plateFE.rms_vec_from_modal_cpsq(Sqq_F,sqq_type = sqq_type)
# Convert into velocity
v2rms_F = (omega.angular_frequency*x2rms_F).flatten()

#%% energy plots
RPR_FE_force.energy.plot(2,xscale = 'log',yscale = 'log',
                       fulllegstr = ('room1','room2'))
plt.xlabel('$f_c/$Hz')
plt.ylabel('$E/$ J')
#plt.ylim(bottom=0.1)
ax = plt.gca()
ax.axes.xaxis.set_ticks([])
plt.xticks(2*np.pi*fc,fclabels)
plt.legend()
plt.show()


#%% eta_alpha
plt.figure(3)
plt.loglog(freq,eta_alpha[0,:],label = 'room 1')
plt.loglog(freq,eta_alpha[1,:],label = 'room 2')
plt.xlabel('$f_c/$Hz')
plt.ylabel('$\eta_{\\alpha}$')
#plt.ylim(bottom=0.1)
plt.xticks(fc,fclabels)
plt.legend()
plt.show()


#%% engineering units force source
plt.figure(5)
RPR_FE_force.result.plot(5,ID=[1,3],xscale = 'log',yscale = 'log',
                       fulllegstr = ('room 1','room 2'))
plt.xlabel('$f_c/$Hz')
plt.ylabel('$p_{\\rm rms}/$Pa')
#plt.ylim(bottom=0.1)
ax = plt.gca()
ax.axes.xaxis.set_ticks([])
plt.xticks(2*np.pi*fc,fclabels)
plt.legend(loc = 3)
plt.tight_layout()
plt.show()

plt.savefig('../source/images/hybrid_RPR_force_pressure.png')

    

# %% power input from force at FE 
plt.figure(8)
RPR_FE_force.loads['FE_10N@Node'+str(forceID)].plot(8,ID=1,xscale = 'log',yscale = 'log',fulllegstr = ['room 1\&2'])
#plt.semilogx(freq,-10*np.log10(tau_diff),':',label = 'mass law TMM')
ax = plt.gca()
ax.axes.xaxis.set_ticks([])
plt.xticks(2*np.pi*fc,fclabels)
plt.xlabel('$f_c/$Hz')
plt.ylabel('$\Pi_{\\rm in}/$W')
plt.legend()
plt.tight_layout()
plt.show()

    
# %% engineering units input power

plt.figure(10)
RPR_FE_force.hybrid_result.plot(10,ID=2,xscale = 'log',yscale = 'log',fulllegstr = ['$S_{qq}$'])
plt.loglog(2*np.pi*freq,v2rms_F,'r',label = '$S_{qq,rev}$')
ax = plt.gca()
ax.axes.xaxis.set_ticks([])
plt.xticks(2*np.pi*fc,fclabels)
plt.xlabel('$f_c/$Hz')
plt.ylabel('$v_{\\rm rms}/$m s$^{-1}$')
plt.tight_layout()
plt.legend(loc=3)
plt.show()

#plt.savefig('../source/images/hybrid_RPR_force_velocity.png')    





