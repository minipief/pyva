# -*- coding: utf-8 -*-
"""
Two rooms and FE-plate under sound source excitation 
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
HJ123 = jun.HybridAreaJunction((room1,room2),plateFE)

# define loads
# room power
power1Watt = lC.Load(omega, np.ones(omega.shape), dof.DOF(1,0,dof.DOFtype(typestr = 'power')), name = '1Watt')

# Calculate CLF for debugging
#eta, eta_alpha,power_in,x_modal = HJ123.CLF(omega.angular_frequency,Signal = False,force='10N@Node200')
eta, eta_alpha = HJ123.CLF(omega.angular_frequency) #,Signal = False,force='10N@Node200')
tau_from_eta = 8*np.pi**2*room1.modal_density(omega.data,0)*omega.data/room1.fluid.wavenumber(omega.data)**2/plate.area*eta
tau_from_eta = tau_from_eta.flatten() 

# Create hybrid SEA model
RPR_SEA_exc = mds.HybridModel((room1,room2),FEsystems = (plateFE,),xdata=omega)
# connect and add load
RPR_SEA_exc.add_hybrid_junction({'HareaJ_12':HJ123})
RPR_SEA_exc.add_load('1Watt',power1Watt) 

# %% Solve hybrid SEA models
RPR_SEA_exc.create_SEA_matrix(sym = 1)
RPR_SEA_exc.solve()

#%% exctract FEM resonse using junction methods
sqq_type = dof.DOFtype(typestr='displacement',exponent = 2)
# Determine CSPD of FEM system
Sqq_P = HJ123.FEM_response(omega.angular_frequency , RPR_SEA_exc.energy)
# Determine nodal average from modal response
x2rms_P,v_type = plateFE.rms_vec_from_modal_cpsq(Sqq_P,sqq_type = sqq_type)
# Convert into velocity
v2rms_P = (omega.angular_frequency*x2rms_P).flatten()


#%% energy plots
RPR_SEA_exc.energy.plot(1,xscale = 'log',yscale = 'log', 
                       fulllegstr = ('room1','room2'))
plt.xlabel('$f_c/$Hz')
plt.ylabel('$E/$ J')
#plt.ylim(bottom=0.1)
ax = plt.gca()
ax.axes.xaxis.set_ticks([])
plt.xticks(2*np.pi*fc,fclabels)
plt.legend()
plt.show()

plt.savefig('../source/images/hybrid_RPR_SEA_energy.png')

#%% eta_alpha
plt.figure(3)
plt.loglog(freq,eta_alpha[0,:],label = 'cavity 1')
plt.loglog(freq,eta_alpha[1,:],label = 'cavity 2')
plt.xlabel('$f_c/$Hz')
plt.ylabel('$\eta_{\\alpha}$')
#plt.ylim(bottom=0.1)
plt.xticks(fc,fclabels)
plt.legend()
plt.show()


#%% engineering units power source
#RPR_SEA_exc.result.plot(4,ID=[1,3],xscale = 'log',yscale = 'log',
#                       fulllegstr = ('room1','room2'))
RPR_SEA_exc.result.plot(4,ID=[1,3],xscale = 'log',yscale = 'log',
                       fulllegstr = ('room 1','room 2',))
#plt.plot(SEA_omega,PR_ydata[1],color = 'grey')
plt.xlabel('$f_c/$Hz')
plt.ylabel('$p_{\\rm rms}/$Pa')
#plt.ylim(bottom=0.1)
ax = plt.gca()
ax.axes.xaxis.set_ticks([])
plt.xticks(2*np.pi*fc,fclabels)
plt.tight_layout()
plt.legend(loc = 3)
plt.show()

plt.savefig('../source/images/hybrid_RPR_SEA_pressure.png')

  
# %% engineering units velocity
plt.figure(9)
plt.loglog(2*np.pi*freq,v2rms_P,label = '$S_{qq,rev}$')
ax = plt.gca()
ax.axes.xaxis.set_ticks([])
plt.xticks(2*np.pi*fc,fclabels)
plt.xlabel('$f_c/$Hz')
plt.ylabel('$v_{\\rm rms}/$m s$^{-1}$')
plt.tight_layout()
plt.legend(loc=3)
plt.show()

plt.savefig('../source/images/hybrid_RPR_SEA_velocity.png')
  
# %% TL from delta Lp
    
p1 = RPR_SEA_exc.result[0].ydata.flatten()
p2 = RPR_SEA_exc.result[1].ydata.flatten()

tau = (p2/p1)**2
plt.close(7)
plt.figure(7)
plt.semilogx(freq,-10*np.log10(tau),label = 'plate TL')
#plt.semilogx(freq,-10*np.log10(tau_from_eta),label = 'plate TL fron tau')
ax = plt.gca()
ax.axes.xaxis.set_ticks([])
plt.xticks(fc,fclabels)
plt.xlabel('$f_c/$Hz')
plt.ylabel('TL/dB')
plt.legend()
plt.tight_layout()
plt.show()
plt.savefig('../source/images/hybrid_RPR_SEA_TL.png')
     