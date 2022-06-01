# -*- coding: utf-8 -*-
"""
Created on Sat Oct 10 16:55:57 2020

two plate script for book presentation 

"""
import numpy as np
import matplotlib.pyplot as plt

import pyva.models as mds
import pyva.coupling.junctions as con
import pyva.systems.infiniteLayers as iL
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

# x-axis tics
fc,fclabels = uf.get_3rd_oct_axis_labels()
fc,fclabels = fc[1:],fclabels[1:]

# Frequency range
omega = mC.DataAxis.octave_band(f_max=2*np.pi*10000)
om    = omega.data

# Plate dimensions
Lx = 1.2
Ly = 1
Lz = 1.1

# Room dimensions
V = Lx*Ly*Lz
A = 2*(Lx*Ly+Ly*Lz+Lx*Lz)
P = 4*(Lx+Ly+Lz)


# Plate thickness
t = 2.0E-3;

# junction properties
angle_R = (0.,90*np.pi/180)

# Create materials
steel = matC.IsoMat(E=210e9,nu=0.3,rho0=7800, eta = 0.0)
air   = matC.Fluid()

# Fibre material
rho_bulk = 0.98*1.20 + 30.
fibre1 = matC.EquivalentFluid(porosity = 0.98, \
                               flow_res = 25000.,\
                               tortuosity = 1.02, \
                               length_visc = 90.e-6, \
                               length_therm = 180.e-6,\
                               rho_bulk = rho_bulk , \
                               rho0 = 1.208, \
                               dynamic_viscosity = 1.81e-5 )

# Thickness of fibre layer
h1 = 0.05

# Create props
steel2mm = stPC.PlateProp(t,steel)
# Create Layer

fibre_5cm = iL.FluidLayer(h1,fibre1)
heavy_1kg = iL.MassLayer(0.001, 1000)
# and noise control treatement
nct_mass = mds.TMmodel((fibre_5cm,heavy_1kg)) 
# change order for absorption
nct_mass4abs = mds.TMmodel((heavy_1kg,fibre_5cm))
nct          = mds.TMmodel((fibre_5cm,)) 



# absorption due to trim
alpha_nct = nct.absorption_diffuse(om,in_fluid=air)
alpha_nct_mass = nct_mass4abs.absorption_diffuse(om,in_fluid=air)

# %% plot abs fig 10 
alpha_nct.plot(10,xscale='log',yscale='linear',fulllegstr = ('5cm fibre',))
alpha_nct_mass.plot(10,xscale='log',yscale='linear',fulllegstr = ('5cm fibre + 1kg',), ls = '--')
plt.xlabel('$f_c/$Hz')
plt.ylabel('$\\alpha_s$')
plt.xticks(2*np.pi*fc,fclabels)
plt.tight_layout()
plt.legend()
plt.show()

# %% create absorption
area_dof = dof.DOF(0,0,dof.DOFtype(typestr='area'))
abs_area = alpha_nct*(Lx*Ly)+alpha_nct_mass*(A-Lx*Ly)
# This is a signal those dof must be set to the
abs_area.dof = area_dof

# %% plot 1

abs_area.plot(1,xscale = 'log')

plt.savefig('../source/images/box_isolated_abs_area.png')

# %% plates

f_c = steel2mm.coincidence_frequency(air.c0)/2/np.pi
print('coincidence frequency = {0}'.format(f_c))

# Create plate subsystems
plate1 = st2Dsys.RectangularPlate(1,Lx,Lz,prop = steel2mm,trim=(nct_mass,'none'))
plate2 = st2Dsys.RectangularPlate(2,Ly,Lz,prop = steel2mm,trim=(nct_mass,'none'))
plate3 = st2Dsys.RectangularPlate(3,Lx,Lz,prop = steel2mm,trim=(nct_mass,'none'))
plate4 = st2Dsys.RectangularPlate(4,Ly,Lz,prop = steel2mm,trim=(nct_mass,'none'))
plate5 = st2Dsys.RectangularPlate(5,Lx,Ly,prop = steel2mm,trim=(nct_mass,'none'))

room   = ac3Dsys.Acoustic3DSystem(6, V , A, P, air,\
                        absorption_area = abs_area ,\
                        damping_type= ['eta','surface']    )



# use receiving room with abso
sif1 = con.SemiInfiniteFluid((room,plate1), air)
sif2 = con.SemiInfiniteFluid((room,plate2), air)
sif3 = con.SemiInfiniteFluid((room,plate3), air)
sif4 = con.SemiInfiniteFluid((room,plate4), air)
sif5 = con.SemiInfiniteFluid((room,plate5), air)


juncs =  { 'j61' : con.AreaJunction((room,plate1)), 
           'j62' : con.AreaJunction((room,plate2)),
           'j63' : con.AreaJunction((room,plate3)),
           'j64' : con.AreaJunction((room,plate4)),
           'j65' : con.AreaJunction((room,plate5)),
           'j12' : con.LineJunction((plate1,plate2), Lz, angle_R),
           'j23' : con.LineJunction((plate2,plate3), Lz, angle_R),
           'j34' : con.LineJunction((plate3,plate4), Lz, angle_R),
           'j41' : con.LineJunction((plate4,plate2), Lz, angle_R),
           'j15' : con.LineJunction((plate1,plate5), Lx, angle_R),
           'j25' : con.LineJunction((plate2,plate5), Ly, angle_R),
           'j35' : con.LineJunction((plate3,plate5), Lx, angle_R),
           'j45' : con.LineJunction((plate4,plate5), Ly, angle_R)
           }

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

box.add_load('1Watt',power1Watt) # add 1Watt per band 



#%% Solve
box.create_SEA_matrix(sym = 1)
box.solve()


#%% plotting 3

box.result.plot(2,ID = [1,5],xscale = 'log',yscale = 'log',fulllegstr=('a)','a)'))
plt.figure(2)
plt.yscale('log')
plt.ylim(1e-5,0.01)
plt.xticks(2*np.pi*fc,fclabels)
plt.xlabel('$f_c/$Hz')
#plt.ylabel('$TL/$dB')
plt.legend()
plt.show()

#%% plotting 3
 
box.result.plot(3,ID = [6],xscale = 'log',yscale = 'linear',fulllegstr=('a)',))   

plt.xticks(2*np.pi*fc,fclabels)
plt.xlabel('$f_c/$Hz')
plt.ylabel('$p_{\\rm rms}/$Pa')
plt.tight_layout()
plt.legend()
plt.show()

#%% more info
sif1_in = box.power_input('sif1')
sif2_in = box.power_input('sif2')
sif3_in = box.power_input('sif3')
sif4_in = box.power_input('sif4')
sif5_in = box.power_input('sif5')

sif_all = sif1_in.sum()+sif2_in.sum()+sif3_in.sum()+sif4_in.sum()+sif5_in.sum()


#%% plot 4

sif1_in.plot(4,xscale='log',yscale='linear')
plt.figure(4)
#plt.plot(om_VA1,power_in_1_res,':',label = 'VA1 res')
#plt.plot(om_VA1,power_in_1_nonres,':',label = 'VA1 non-res')
plt.yscale('log')
plt.xlabel('$f_c/$Hz')
plt.xticks(2*np.pi*fc,fclabels)
plt.ylabel('$\Pi_{in}/$W')
plt.legend()
plt.tight_layout()
plt.show()

#%% plot 5

sif_all.plot(5,xscale='log',res = 'dB',fulllegstr=('b)',))
plt.figure(5)
plt.xlabel('$f_c/$Hz')
plt.xticks(2*np.pi*fc,fclabels)
#plt.yticks([40,60,80,100,],[40,60,80,100,])
#plt.ylabel('$\Pi_{in}/$W')
plt.legend()
plt.tight_layout()
plt.show()

#plt.savefig('../source/images/box_isolated_SIF_total.png')



#%% plot 6    
NR = 1/sif_all.ydata[0,:]

plt.figure(6)
plt.semilogx(om,10*np.log10(NR),label = 'b)')
plt.xlabel('$f_c/$Hz')
plt.xticks(2*np.pi*fc,fclabels)
plt.ylabel('$IL/$dB')
plt.legend()
plt.tight_layout()
plt.show()




    





