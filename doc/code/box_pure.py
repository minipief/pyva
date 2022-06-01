# -*- coding: utf-8 -*-
"""
Created on Sat Oct 10 16:55:57 2020

two plate script for book presentation 

"""
import numpy as np
import matplotlib.pyplot as plt

import pyva.models as mds
import pyva.coupling.junctions as jun
import pyva.properties.structuralPropertyClasses as stPC
import pyva.systems.structure2Dsystems as st2Dsys
import pyva.systems.acoustic3Dsystems as ac3Dsys
import pyva.loads.loadCase as lC

import pyva.useful as uf

# my packages

import pyva.data.dof as dof
import pyva.data.matrixClasses as mC
import pyva.properties.materialClasses as matC

plt.close('all')

# x-axis tics for better readability
fc,fclabels = uf.get_3rd_oct_axis_labels()
fc,fclabels = fc[1:],fclabels[1:]

# Frequency range
omega = mC.DataAxis.octave_band(f_max=2*np.pi*10000)

# Plate dimensions
Lx = 1.2
Ly = 1
Lz = 1.1

# Box dimensions
V = Lx*Ly*Lz
A = 2*(Lx*Ly+Ly*Lz+Lx*Lz)
P = 4*(Lx+Ly+Lz)

# Plates thickness
t = 2.0E-3;

# junction angle
angle_R = (0.,90*np.pi/180)

# Create materials
steel = matC.IsoMat(E=210e9,nu=0.3,rho0=7800, eta = 0.0)
air   = matC.Fluid() 

# Create props
steel2mm = stPC.PlateProp(t,steel)


area_dof = dof.DOF(0,0,dof.DOFtype(typestr='area'))


# %% create models
f_c = steel2mm.coincidence_frequency(air.c0)/2/np.pi
print('coincidence frequency = {0}'.format(f_c))

# Create plate subsystems
plate1 = st2Dsys.RectangularPlate(1,Lx,Lz,prop = steel2mm)
plate2 = st2Dsys.RectangularPlate(2,Ly,Lz,prop = steel2mm)
plate3 = st2Dsys.RectangularPlate(3,Lx,Lz,prop = steel2mm)
plate4 = st2Dsys.RectangularPlate(4,Ly,Lz,prop = steel2mm)
plate5 = st2Dsys.RectangularPlate(5,Lx,Ly,prop = steel2mm)

room     = ac3Dsys.Acoustic3DSystem(6, V , A, P, air)

# create semi infinite fluids
sif1 = jun.SemiInfiniteFluid((room,plate1), air)
sif2 = jun.SemiInfiniteFluid((room,plate2), air)
sif3 = jun.SemiInfiniteFluid((room,plate3), air)
sif4 = jun.SemiInfiniteFluid((room,plate4), air)
sif5 = jun.SemiInfiniteFluid((room,plate5), air)

juncs =  { 'j61' : jun.AreaJunction((room,plate1)), 
           'j62' : jun.AreaJunction((room,plate2)),
           'j63' : jun.AreaJunction((room,plate3)),
           'j64' : jun.AreaJunction((room,plate4)),
           'j65' : jun.AreaJunction((room,plate5)),
           'j12' : jun.LineJunction((plate1,plate2), Lz, angle_R),
           'j23' : jun.LineJunction((plate2,plate3), Lz, angle_R),
           'j34' : jun.LineJunction((plate3,plate4), Lz, angle_R),
           'j41' : jun.LineJunction((plate4,plate2), Lz, angle_R),
           'j15' : jun.LineJunction((plate1,plate5), Lx, angle_R),
           'j25' : jun.LineJunction((plate2,plate5), Ly, angle_R),
           'j35' : jun.LineJunction((plate3,plate5), Lx, angle_R),
           'j45' : jun.LineJunction((plate4,plate5), Ly, angle_R)
            }



#om = 2*np.pi*np.logspace(2,np.log10(10000),3)




# define load
# room power
power1Watt = lC.Load(omega, np.ones(omega.shape), dof.DOF(6,0,dof.DOFtype(typestr = 'power')), name = '1Watt')
# plate power
#power1Watt = lC.Load(omega, np.ones(omega.shape), dof.DOF(2,3,dof.DOFtype(typestr = 'power')), name = '1Watt')

#create SEA model
box = mds.HybridModel((plate1,plate2,plate3,plate4,plate5,room),xdata=omega)

#connect both
box.add_junction(juncs)
box.add_SIF({'sif1' : sif1, 
             'sif2' : sif2, 
             'sif3' : sif3, 
             'sif4' : sif4, 
             'sif5' : sif5})

box.add_load('1Watt',power1Watt) # add 1mWatt per band 


#%% solving
box.create_SEA_matrix(sym = 1)
box.solve()

#%% plotting 1

box.result.plot(1,ID = [1,5],xscale = 'log',yscale = 'log',fulllegstr=('a)','a)'))
plt.figure(1)
plt.yscale('log')
plt.ylim(1e-5,0.01)
plt.xticks(2*np.pi*fc,fclabels)
plt.xlabel('$f_c/$Hz')
#plt.ylabel('$TL/$dB')
plt.legend()
plt.show()

#%% plotting 2
 
box.result.plot(2,ID = [6],xscale = 'log',yscale = 'linear',fulllegstr=('a)',))   
plt.xticks(2*np.pi*fc,fclabels)
plt.xlabel('$f_c/$Hz')
plt.ylabel('$p_{\\rm rms}/$Pa')
plt.tight_layout()
plt.legend()
plt.show()

# plt.savefig('../source/images/box_pure_pressure.png')

#%% more info
sif1_in = box.power_input('sif1')
sif2_in = box.power_input('sif2')
sif3_in = box.power_input('sif3')
sif4_in = box.power_input('sif4')
sif5_in = box.power_input('sif5')

sif_all = sif1_in.sum()+sif2_in.sum()+sif3_in.sum()+sif4_in.sum()+sif5_in.sum()


#%% plot 3

sif1_in.plot(3,xscale='log',yscale='linear')
plt.figure(3)
#plt.plot(om_VA1,power_in_1_res,':',label = 'VA1 res')
#plt.plot(om_VA1,power_in_1_nonres,':',label = 'VA1 non-res')
plt.yscale('log')
plt.xlabel('$f_c/$Hz')
plt.xticks(2*np.pi*fc,fclabels)
plt.ylabel('$\Pi_{in}/$W')
plt.legend()
plt.tight_layout()
plt.show()

#plt.savefig('../source/images/box_pure_SIF_power_in.png')


#%% plot 4

sif_all.plot(4,xscale='log',res = 'dB',fulllegstr=('a)',))
plt.figure(4)
plt.xlabel('$f_c/$Hz')
plt.xticks(2*np.pi*fc,fclabels)
#plt.yticks([40,60,80,100,],[40,60,80,100,])
plt.ylim((60,110))
plt.legend()
plt.tight_layout()
plt.show()

#plt.savefig('../source/images/box_pure_SIF_total.png')

#%% plot 5    
NR = 1/sif_all.ydata[0,:]

plt.figure(5)
plt.semilogx(omega.data,10*np.log10(NR),label = 'a)')
plt.xlabel('$f_c/$Hz')
plt.xticks(2*np.pi*fc,fclabels)
plt.ylabel('$IL/$dB')
plt.legend()
plt.tight_layout()
plt.show()




    





