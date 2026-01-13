# -*- coding: utf-8 -*-
"""
Script related to the paper from Price and Crocker [1] on the sound transmission of double panels.

[1] A. J. Price and M. J. Crocker, “Sound Transmission through Double Panels Using 
    Statistical Energy Analysis, ”The Journal of the Acoustical Society of America, 
    vol. 47, no. 3A, pp. 683–693, Mar. 1970, doi: 10.1121/1.1911951.
"""

import numpy as np
import matplotlib.pyplot as plt

# pyva main model module 
import pyva.models as mds

# Coupling module
import pyva.coupling.junctions as con    # junctions module

# Property modules
import pyva.properties.structuralPropertyClasses as stPC
import pyva.properties.materialClasses as matC

# Layer module for transfermatrix methods
import pyva.systems.infiniteLayers as iL

# SEA system modules
import pyva.systems.structure2Dsystems as st2Dsys # plate systems module
import pyva.systems.acoustic3Dsystems as ac3Dsys  # cavity systems moduke 
import pyva.loads.loadCase as lC                  # loadcase module

# Modules for system matrices and degree of freedom
import pyva.data.dof as dof
import pyva.data.matrixClasses as mC

# Pyva utilities
import pyva.useful as uf

aprop  = dict(facecolor='black', shrink=0.01, width=1 )

#%% Model settings

# Create frequency axis labels
fc,fclabels = uf.get_3rd_oct_axis_labels(range='SEA')
plt.close('all')

#%% Frequency range
omega = mC.DataAxis.octave_band(f_min=2*np.pi*100,f_max=2*np.pi*10000,bands_per_octave = 3)
# Create according ndarrays for plotting
om    = omega.data
freq  = om/2/np.pi 

# Panel dimensions
L1 = 1.55   # panel/DW cavity width
L2 = 0.071  # DW cavity thickness
L3 = 1.97   # panel/DW cavity height

L1_fig3 = 2.2 # height of cavity from fig 3

S = L1*L3   # panel surface
S_dwc_edge = L2*2*(L1+L3)# DW cavity edge surface

# Cavity volumes
V_dwc  = L1*L2*L3 # Volume DW cavity
V1 = 100 # room 1
V2 = 100 # room 2

# Absorption areas in large cavities
As1  = 2. # arbitrary room absorption area
As2  = S  # secong room absorption area equal to junction area for log10(S/As2)=0

# Materials
alu = matC.IsoMat(eta=0) # Alu is the default isomaterial
air = matC.Fluid(eta=0)  # Air is the default fluid

# Plate properties
alu_3_18mm = stPC.PlateProp(0.00318,alu) # 0.85 g/cm^2
alu_6_36mm = stPC.PlateProp(0.00636,alu) # 1.71 g/cm^2
alu_0_88mm = stPC.PlateProp(0.00088,alu) # 0.24 g/cm^2
alu_9_26mm = stPC.PlateProp(0.00925,alu) # 2.5  g/cm^2

#%% Tables for DW cavity damping
omega_min = 2*np.pi*100     # in all arguments pyva uses angular frequency
omega_max = 2*np.pi*10000
freq_tab_BI = mC.DataAxis.octave_band(omega_min,omega_max,typeID = 21)
 
# Create table until 10000Hz
eta_DWC_data = np.zeros(np.shape(freq_tab_BI))
# Fill data until 2000 Hz
eta_DWC_data[0,:14] = 0.01*np.array([[3.18, 2.29, 2.38, 2.89, 2.68, 
                         2.72, 2.79, 3.27, 3.09, 2.80,
                         2.38, 2.00, 1.60, 1.27 ]])

# Fill data from 3200 - 10000 Hz with equation  B11
def B11(omega):
    """
    Equation B11 of [1].

    Parameters
    ----------
    omega : float
        angular frequency.

    Returns
    -------
    float
        damping loss of cavity.

    """
    # alpha = 1
    return air.c0*S_dwc_edge/6/omega/V_dwc

eta_DWC_data[0,14:] = B11(freq_tab_BI.data[14:])
# Create damping loss as Signal 
eta_DWC = mC.Signal(freq_tab_BI, eta_DWC_data, dof.DOF(0,0,dof.DOFtype(typeID = 1)))



# Crossmode Frequency f_d
f_d = air.c0/2/L2 
print('Crossmode frequency: {0:.1f}'.format(f_d))

#%% Double wall resonance
m1 = alu_3_18mm.mass_per_area
m2 = alu_6_36mm.mass_per_area
m3 = alu_0_88mm.mass_per_area
m4 = alu_9_26mm.mass_per_area

def DW_frequency (m1,m2,d):
    return 0.5/np.pi*np.sqrt(air.rho0*air.c0**2/d/m1/m2*(m1+m2))

print('Double wall resonance Fig.5: {0:.1f}Hz'.format(DW_frequency(m1,m1,L2)))
print('Double wall resonance Fig.6: {0:.1f}Hz'.format(DW_frequency(m1,m2,L2)))
print('Double wall resonance Fig.7: {0:.1f}Hz'.format(DW_frequency(m1,m3,L2)))



#%% Subsystem definition
# Plate subsystems
plate1_fig_5 = st2Dsys.RectangularPlate(2, L1,L3,prop=alu_3_18mm, eta = 0.005) 
plate2_fig_5 = st2Dsys.RectangularPlate(4, L1,L3,prop=alu_3_18mm, eta = 0.005)
 
plate1_fig_6 = st2Dsys.RectangularPlate(2, L1,L3,prop=alu_6_36mm, eta = 0.005) 
plate2_fig_6 = st2Dsys.RectangularPlate(4, L1,L3,prop=alu_3_18mm, eta = 0.005) 

plate1_fig_7 = st2Dsys.RectangularPlate(2, L1,L3,prop=alu_0_88mm, eta = 0.005) 
plate2_fig_7 = st2Dsys.RectangularPlate(4, L1,L3,prop=alu_0_88mm, eta = 0.005) 

plate1_fig_10 = st2Dsys.RectangularPlate(2, L1,L3,prop=alu_9_26mm, eta = 0.005) 
plate2_fig_10 = st2Dsys.RectangularPlate(4, L1,L3,prop=alu_9_26mm, eta = 0.005) 


# Cavity subsystems
room1 = ac3Dsys.Acoustic3DSystem(1, V1, 0, 0, air, absorption_area = As1, damping_type= ['surface'] )
room2 = ac3Dsys.Acoustic3DSystem(5, V2, 0, 0, air, absorption_area = As2, damping_type= ['surface'])

# Double wall cavity
# The flat_cavity_sw triggers the doubled radiation efficiency from plate to DWC
dwc     = ac3Dsys.RectangularRoom(3, L1, L2, L3, air, eta = eta_DWC,flat_cavity_sw = True)
# Extra cavity of Fig. 3 that uses other dimensions
dwc_fig3 = ac3Dsys.RectangularRoom(3, L1_fig3, L2, L3, air, eta = eta_DWC,flat_cavity_sw = True)

#%% Modal densities (of other modified cavity) Fig. 3
plt.figure(1)
plt.plot(freq,dwc_fig3.modal_density(om),label = 'Maa')
# Calculate modal density from mode counting
n_num,om_lim = dwc_fig3.modal_density_precise(om)

plt.plot(freq[1:],n_num, label = 'Mode counting')
#plt.plot(freq,dwc_fig3.modal_density(om,method = 'Price'),label = 'Price mod dens')
plt.plot(freq,(L1_fig3*L3*om/2/np.pi/air.c0**2),label = 'Eq (28)')
plt.plot(freq,(L1_fig3*L3*L2*om**2/2/np.pi**2/air.c0**3),label = 'Eq (29)')
#plt.plot(freq,dwc_fig3.modal_density(om,method = 'Price'),label = 'Price mod dens')

plt.annotate(r'$c_0/2L_2$', xy=(f_d,L1_fig3*L3*f_d/air.c0**2), xytext=(5000,0.05),
              arrowprops=aprop)
plt.title("Fig.3")
plt.xlim((100,1e4))
plt.ylim((1e-2,1))
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$f/$Hz')
plt.ylabel('$n(\omega)/$s')
plt.xticks(fc,fclabels)
plt.legend()
plt.grid(which = 'both')
plt.show


#%% Junctions Fig 5 Model
J123_fig_5 = con.AreaJunction((room1,plate1_fig_5,dwc))
J345_fig_5 = con.AreaJunction((dwc,plate2_fig_5,room2))

#%% Junctions Fig 6 Model
J123_fig_6 = con.AreaJunction((room1,plate1_fig_6,dwc))
J345_fig_6 = con.AreaJunction((dwc,plate2_fig_6,room2))

#%% Junctions Fig 7 Model
J123_fig_7 = con.AreaJunction((room1,plate1_fig_7,dwc))
J345_fig_7 = con.AreaJunction((dwc,plate2_fig_7,room2))

#%% Junctions Fig 10 Model
J123_fig_10 = con.AreaJunction((room1,plate1_fig_10,dwc))
J345_fig_10 = con.AreaJunction((dwc,plate2_fig_10,room2))

#%% Load definition
# Source in cavity 1
power1mWatt = lC.Load(omega, 0.001*np.ones(omega.shape), 
                      dof.DOF(1,0,dof.DOFtype(typestr = 'power')), name = '1mWatt')


#%% Select to Plot Fig 5 (and 8),6 7 or 10
fig_nr = 10

if fig_nr == 5:
    # SEA Model used in Fig 5
    SEA_model = mds.HybridModel((room1,plate1_fig_5,dwc,plate2_fig_5,room2),xdata=omega)
    # Add junctions
    SEA_model.add_junction({'areaJ_123':J123_fig_5})
    SEA_model.add_junction({'areaJ_345':J345_fig_5})
elif fig_nr == 6:
    # SEA Model used in Fig 6
    SEA_model = mds.HybridModel((room1,plate1_fig_6,dwc,plate2_fig_5,room2),xdata=omega)
    # Add junctions
    SEA_model.add_junction({'areaJ_123':J123_fig_6})
    SEA_model.add_junction({'areaJ_345':J345_fig_6})
elif fig_nr == 7:
    # SEA Model used in Fig 6
    SEA_model = mds.HybridModel((room1,plate1_fig_7,dwc,plate2_fig_7,room2),xdata=omega)
    # Add junctions
    SEA_model.add_junction({'areaJ_123':J123_fig_7})
    SEA_model.add_junction({'areaJ_345':J345_fig_7})
elif fig_nr == 10:
    # SEA Model used in Fig 6
    SEA_model = mds.HybridModel((room1,plate1_fig_10,dwc,plate2_fig_10,room2),xdata=omega)
    # Add junctions
    SEA_model.add_junction({'areaJ_123':J123_fig_10})
    SEA_model.add_junction({'areaJ_345':J345_fig_10})
else:
    raise ValueError

# Add loads (same for all cases)
SEA_model.add_load('1mWatt',power1mWatt) # add 1Watt per band 
#%% Solve models Fig 5-7
if fig_nr in {5,6,7}:
    SEA_model.create_SEA_matrix(sym = 1)
    SEA_model.solve()

    # Derive paths of power input
    pow_in_room2 = SEA_model.power_input(5)

    #%% Engineering units velocity - plate results
    plt.figure(2)
    plt.title('Plate velocity of plates from fig {0}'.format(fig_nr))
    
    SEA_model.result.plot(2,ID=[2,4],dof=[3,3],xscale = 'log',yscale = 'log',
                           fulllegstr = ('plate B1','plate B2'),
                           xlabel = '$f/$Hz', ylabel = '$v_{\\rm rms}/$m s$^{-1}$',
                           xticks = 2*np.pi*fc[2:], xticklabels = fclabels[2:])
    
    
    #%% Engineering units pressure - cavity results
    plt.figure(3)
    plt.title('Cavity pressure from fig {0}'.format(fig_nr))
    SEA_model.result.plot(3,ID=[1,3,5],xscale = 'log',yscale = 'log',ls = ['-','--',':','-.'],
                        fulllegstr = ('room 1','DWC','room 2'),
                        xlabel = '$f/$Hz', ylabel = '$p_{\\rm rms}/$Pa',
                        xticks = 2*np.pi*fc[2:], xticklabels = fclabels[2:])
    
    
    
    #%% plot comparative modes in band 
    plt.figure(4)
    plt.plot(freq,plate1_fig_5.modes_in_band(om,3),label = 'plate 1 Fig 5')
    plt.plot(freq,plate1_fig_6.modes_in_band(om,3),label = 'plate 1 Fig 6')
    plt.plot(freq,plate1_fig_7.modes_in_band(om,3),label = 'plate 1 Fig 7')
    #plt.plot(freq,plate2.modes_in_band(om,3),label = 'plate 2')
    plt.plot(freq,dwc.modes_in_band(om),label = 'DWC')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$f/$Hz')
    plt.ylabel('N')
    plt.title('Modes in 3rd oct band')
    plt.ylim(bottom=0.1)
    plt.xticks(fc,fclabels)
    plt.legend()
    plt.grid(which = 'both')
    plt.show
 
#%% Plot Fig 9 radation efficiency
sigma1 = plate1_fig_5.radiation_efficiency(om,fluid = air)
sigma2 = plate1_fig_5.radiation_efficiency_simple(om,fluid = air)

# Show coincidence frequencies
f_c1 = alu_3_18mm.coincidence_frequency(air.c0)/2/np.pi
print('3.18 mm alu plate coincidence frequency = {0:.1f}'.format(f_c1))
f_c2 = alu_6_36mm.coincidence_frequency(air.c0)/2/np.pi
print('6.36 mm alu plate coincidence frequency = {0:.1f}'.format(f_c2))
f_c3 = alu_0_88mm.coincidence_frequency(air.c0)/2/np.pi
print('0.88 mm alu plate coincidence frequency = {0:.1f}'.format(f_c3))
f_c4 = alu_9_26mm.coincidence_frequency(air.c0)/2/np.pi
print('0.88 mm alu plate coincidence frequency = {0:.1f}'.format(f_c4))

plt.figure(5)
plt.plot(freq,10*np.log10(sigma1),'-.',label = 'Leppington')
plt.plot(freq,10*np.log10(sigma2),'-.',label = 'ISO EN 12354-1')
plt.title('Fig. 9')
plt.xlabel('$f/$Hz')
plt.ylabel('$\sigma/dB$')
plt.xscale('log')
plt.ylim(-50,25)
plt.xticks(fc,fclabels)
plt.legend()
plt.grid(which = 'both')


#%% Calculate transmission loss fig 5-7
# Calculate transmission loss from pressure difference
# Note: This is allowed because wall and absorption surface are equal S=As2
#%% Solve models Fig 5-7
if fig_nr in {5,6,7}:
    p1 = SEA_model.result[0].ydata.flatten()
    p2 = SEA_model.result[6].ydata.flatten()
    tau = (p2/p1)**2

#%% Plot transmission loss fig 5-7
    plt.figure(6)
    plt.semilogx(freq,-10*np.log10(tau),'-.',label = 'Pyva',color = 'C1')
    plt.xticks(fc,fclabels)
    plt.ylim((0,75))
    plt.xlabel('$f/$Hz')
    plt.ylabel('$TL/$dB')
    plt.title('Fig. {0}'.format(fig_nr))
    plt.grid(which = 'both')
    plt.legend()
    plt.show()

#%% Plot input power
    plt.figure(7)
    plt.title('Powert input into room 2 from fig {0}'.format(fig_nr))
    
    pow_in_room2.plot(7,yscale = 'log',fulllegstr = ['$\Pi_{DWC}$','$\Pi_{plate2}$'],
                      xscale = 'log',
                      xlabel = '$f/$Hz', ylabel = '$\Pi_{in}/$W',
                      xticks = 2*np.pi*fc[2:], xticklabels = fclabels[2:])

#%% Redo calculatiwith 5% damping if Fig 5 is selected to generate Fig 8
if fig_nr == 5: # Initiates Fig 8 calulation 
    plate1_fig_5.eta = 0.05
    plate2_fig_5.eta = 0.05

    SEA_model.create_SEA_matrix(sym = 1)
    SEA_model.solve()  

    p1 = SEA_model.result[0].ydata.flatten()
    p2 = SEA_model.result[6].ydata.flatten()
    tau_eta_0_05 = (p2/p1)**2
 
    plt.figure(8)
    plt.semilogx(freq,-10*np.log10(tau),'-.',label = r'$\eta = 0.5$%', color = 'C1')
    plt.semilogx(freq,-10*np.log10(tau_eta_0_05),'-.',label = r'$\eta = 5$%', color = 'C2')
    plt.xticks(fc[2:],fclabels[2:])
    plt.ylim((0,75))
    plt.xlabel('$f/$Hz')
    plt.ylabel('$TL/$dB')
    plt.title('Fig 8')
    plt.grid(which = 'both')
    plt.legend()
    plt.show()


#%% calculate DWC thickness variation of Fig 10 

thicknesses = [0.01,0.2,0.4]
tau = [None]*3
if fig_nr == 10:
    plt.figure(6)

    for i_cav,t in enumerate(thicknesses):
        print('Evaluating thickness {0:.1f}mm'.format(t*1000))
        dwc.Ly = t
        dwc.update()
        print('Double wall frequency for L2={0:.0f}cm : {1:.1f}Hz'.format(t*100,DW_frequency(m4,m4,t)))
        
        SEA_model.create_SEA_matrix(sym = 1)
        SEA_model.solve()  

        p1 = SEA_model.result[0].ydata.flatten()
        p2 = SEA_model.result[6].ydata.flatten()
        tau[i_cav] = (p2/p1)**2
        
        plt.semilogx(freq,-10*np.log10(tau[i_cav])+0*i_cav,label = '$L_2=${0:.0f}cm'.format(t*100))
    
    plt.xticks(fc,fclabels)
    plt.ylim((0,75))
    plt.xlabel('$f/$Hz')
    plt.ylabel('$TL/$dB')
    plt.title('Fig. {0}'.format(fig_nr))
    plt.grid(which = 'both')
    plt.legend()
    plt.show()
    

#%% Plot damping
#%% Engineering units pressure - cavity results
plt.figure(7)
plt.title('Double wall cavity damping')
plt.grid(which='both')

eta_DWC.plot(7,xscale = 'log', fulllegstr = ('$\eta_3$',),
                    xlabel = '$f/$Hz', ylabel = r'$\eta$',
                    xticks = 2*np.pi*fc, xticklabels = fclabels)




#%% Create DW cavity and plate2 as trim
air_damped = matC.Fluid(eta = eta_DWC)
il_air_7_11cm = iL.FluidLayer(L2,air_damped)
il_plate2     = iL.PlateLayer(alu_3_18mm)


