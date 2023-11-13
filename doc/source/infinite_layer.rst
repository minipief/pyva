.. _sec-TMM:

TMM - infinite Layer
--------------------

In contrast to pipe networks the chance for side paths is much lower in infinite layer systems.
A typical example application is the design of sound absorption systems. 

Infinite layers are a special version of one-dimensional system. Such system are one-dimensional in that 
sense, that the properties change in one dimension and remain constant for the other two space dimensions.
The systems or layers are supposed to be infinite in those remaining two dimensions.

In contrast to acoustic networks from section :ref:`sec-acoustic-network` there is an additional 
variable required to describe the propagation parallel to the layers: the wavenumber perpendicular to the main dimension.

.. _fig-infinite-layer:
    
.. figure:: ./images/infinite_layer.*
   :align: center
   :width: 80%
   
   Sketch of connected infinite layers.

In figure :ref:`fig-infinite-layer` such a set-up is shown. The use and application of the
:class:`pyva.models.TMmodel` is shown is this section.

In the first version of pyva TMM was restricted to fluid layer, i.e. layers that can be described by fluid state variables
pressure and velocity. Note, that a this plate can also be described by fluid variable when no friction is considered.
The current version of pyva has implemented :class:`~pyva.systems.infiniteLayers.SolidLayer` and 
:class:`~pyva.systems.infiniteLayers.PoroElasticLayer`. In addition, the Allard approach is implemented that 
allows the combination of multiple layers of different nature. In the last examples of this section 
the new possibilities will be shown. 
 
Absorber design
+++++++++++++++


.. figure:: ./images/absorber.*
   :align: center
   :width: 40%
   
   Absorber configurations
   
In our example we start with the definition of air and the usual fibre material ::

    import pyva.TMmodel as TMM
    import pyva.systems.infiniteLayers as iL

    import pyva.properties.materialClasses as matC
    import pyva.useful as uf

    air    = matC.Fluid(eta = 0.0)
    fibre1 = matC.EquivalentFluid(porosity = 0.98, \
                                   flow_res = 25000.,\
                                   tortuosity = 1.02, \
                                   length_visc = 90.e-6, \
                                   length_therm = 180.e-6,\
                                   rho_bulk = 0.98*1.20 + 30.,\
                                   rho0 = 1.208, \
                                   dynamic_viscosity = 1.81e-5 ) 
                                   
We will design three absorbers:

#. A pure fibre absorber with large thickness (and lot of required space) 
#. The pure thinner fibre absorber
#. A thinner fibre absorber with a perforate front layer

The components or layers of our absorber are created as infinite-layer::

    fibre_10cm  = iL.FluidLayer(h1,fibre1)
    fibre_20cm  = iL.FluidLayer(h2,fibre1)
    perforate   = iL.PerforatedLayer(0.005, 0.001, distance = 0.02) 

And the set-up are created as transfer matrix ::

    TMM_fibre_10      = mods.TMmodel((fibre_10cm,))
    TMM_fibre_20      = mods.TMmodel((fibre_20cm,))
    TMM_perf_fibre_20 = mods.TMmodel((perforate, fibre_10cm,))

The diffuse field absorption is determined with the appropriate method ::

    alpha_fibre_10      = TMM_fibre_10.absorption_diffuse(omega,theta_max=np.pi/2,signal = False)
    alpha_fibre_20      = TMM_fibre_20.absorption_diffuse(omega,theta_max=np.pi/2,signal = False)
    alpha_perf_fibre_10 = TMM_perf_fibre_20.absorption_diffuse(omega,theta_max=np.pi/2,signal = False)

The result can be plotted and shows the typical effects. The thick fibre layer provides excellent absorption 
even at low frequencies. When we take the thinner version we loose low frequency performance. 
In order to improve this we take the perforate cover, but here we loose the high frequency performance due to impedance mismatch at higher
frequencies. 
However, the design of absorbers at given space and weight restrictions is an art. Pyva can support here and provide the tools for absorber design.

.. _fig-TMM_absorber_abs_diffuse:
    
.. figure:: ./images/TMM_absorber_abs_diffuse.*
   :align: center
   :width: 70%
   
   Diffuse absorption of different absorber configurations.

Double walls
++++++++++++

Double walls are the master tool set for efficient acoustic isolation.
Infinite layers are an excellent option to calculate the performance of double walls even though the
theory fails for example for single plates at the coincidence frequency. In such cases a full SEA model 
in twin chamber configuration as in :ref:`sec-two-rooms` section is recommended.

However, for fast and comparative studies on sound isolation the infinite layer in combination with the transfer matrices 
is very useful.

.. figure:: ./images/DW.*
   :align: center
   :width: 40%
   
   Double wall set-up.
   
In this example the materials from above are used with an additional aluminium plate layer ::

    # Create props
    alu = matC.IsoMat()
    alu1mm = stPC.PlateProp(0.001,alu)
    
that is used in combination with several :class:`~pyva.systems.infiniteLayers.FluidLayer` and :class:`~pyva.systems.infiniteLayers.MassLayer` ::

    # Fluid layer
    air_5cm    = iL.FluidLayer(0.05,air)
    fibre_5cm  = iL.FluidLayer(0.05,fibre1)
    # Limp mass layer
    heavy_1kg  = iL.MassLayer(0.001, 1000)
    heavy_2kg7 = iL.MassLayer(0.001, 2700) 
    # Plate layer
    iL_alu1mm  = iL.PlateLayer(alu1mm,)
    
All these layers are compiled in four different versions::

    T_alu            = mds.TMmodel((iL_alu1mm,))
    T_alu_air_mass   = mds.TMmodel((iL_alu1mm,air_5cm,heavy_1kg))
    T_alu_fibre_mass10 = mds.TMmodel((iL_alu1mm,fibre_5cm,heavy_1kg))
    T_alu_fibre_mass27 = mds.TMmodel((iL_alu1mm,fibre_5cm,heavy_2kg7))

Once the TMmodels are created the transmission coefficient is calculated by the :meth:`~pyva.models.TMmodel.transmission_diffuse` method::

    tau_alu = T_alu.transmission_diffuse(omega,signal = False)
    tau_alu_air_mass = T_alu_air_mass.transmission_diffuse(omega,theta_step=np.pi/1000,signal = False)
    tau_alu_fibre_mass10 = T_alu_fibre_mass10.transmission_diffuse(omega,signal = False)
    tau_alu_fibre_mass27 = T_alu_fibre_mass27.transmission_diffuse(omega,signal = False)

Plotting the results in the following figure shows the tremendous increase of isolation above the
double wall resonance. However, in addition is becomes obvious that the inner cavity requires damping and that 
the lower resonance frequencies provide a much better performance in the mid-frequency.

.. figure:: ./images/TMM_DW_transmission.*
   :align: center
   :width: 70%
   
   Transmission loss of different single and double wall configurations.
  
Absorber Design with Poroelastic Materials
++++++++++++++++++++++++++++++++++++++++++

The implementation of Brouards and Allards theory that allows layups of different nature.
Thus examples for such set-ups are required. 
We used the carpet-impervious screen-fibre layup of [All2009]_ in section 11.7.2.
The layup is shown in the following figure.

  .. figure:: ./images/carpet_screen_fibre_layup.*
   :align: center
   :width: 70%   
   
   Layup from [All2009]_ figure 11.16
   
The material parameters and thicknesses are taken from Table 11.7 but with corrected thickness of the fibre layer.

+-----------+-----------+--------++---+----------------+-----------------------+-----------------+-----------------+----------------+------+----------------+
| Material  | :math:`h` | :math:`\Phi`| :math:`\sigma` | :math:`\alpha_\infty` | :math:`\Lambda` | :math:`\Lambda'`| :math:`\rho_1` | E    | :math:`\eta_s` |
|  Unit     |  mm       |  .          |  Ns/m^4        |  .                    |  µm             |  µm             |  kg/m³         |  kPa |                |
+===========+===========+=============+================+=======================+=================+=================+================+======+================+
| Carpet 1  |    3.5    |   0.99      |   5000         |   1                   |        23.      |      28.        |    60.         | 20.  |    0.5         |
+-----------+-----------+-------------+----------------+-----------------------+-----------------+-----------------+----------------+------+----------------+
| Carpet 2  |    3.5    |   0.99      |   5000         |   1                   |        23.      |      28.        |    60.         | 20.  |    0.5         |
+-----------+-----------+-------------+----------------+-----------------------+-----------------+-----------------+----------------+------+----------------+
| Screen    |    3.0    |             |                |                       |                 |                 |  2000.         | 30.  |                |
+-----------+-----------+-------------+----------------+-----------------------+-----------------+-----------------+----------------+------+----------------+
| Carpet 1  |   12.5    |   0.98      |  33000         |   1.1                 |        50.      |     110.        |    60.         | 100. |    0.88        |
+-----------+-----------+-------------+----------------+-----------------------+-----------------+-----------------+----------------+------+----------------+

Two carpet layers are mentioned but with similar properties. Thus, the original paper should be checked.

The code is given in :ref:`sec-abs-poro-examples` and explained in detail here. The carpets material is created by ::

    # Carpet as Poroelastic
    carpet_solid = matC.IsoMat(E=20000.,rho0=60,nu=0.,eta=0.5)
    carpet_mat   = matC.PoroElasticMat(carpet_solid, \
                                flow_res = 5000., \
                                porosity = 0.99, \
                                tortuosity = 1., \
                                length_visc = 23.E-6, length_therm = 28.E-6)

The impervious screen uses the plate property ::

    # Impervious scree as PlateProp
    screen_mat  = matC.IsoMat(E=30000.,rho0 = 2000, nu=0.49)
    screen_prop = stPC.PlateProp(0.003, screen_mat)

And the fibre material is also defined as poroelastic material. ::

    # Fibre as Poroelastic
    fibre_solid = matC.IsoMat(E=100000.,rho0=60,nu=0.,eta=0.88)
    fibre_mat   = matC.PoroElasticMat(fibre_solid, \
                                flow_res = 33000., \
                                porosity = 0.98, \
                                tortuosity = 1.1, \
                                length_visc = 50.E-6, length_therm = 110.E-6)
                                
The single layers are created using the :mod:`~pyva.systems.infiniteLayers` 
module with the related classes ::

    # Define infiniteLayers
    carpet1  = iL.PoroElasticLayer(carpet_mat, 0.0035)
    carpet2  = iL.PoroElasticLayer(carpet_mat, 0.0035)
    screen   = iL.ImperviousScreenLayer(screen_prop)
    fibre    = iL.PoroElasticLayer(fibre_mat, 0.0125)
    
and are connected as :class:`~pyva.models.TMmodel` object::

    # Create lay-up as TMmodel
    TMM_layup = mds.TMmodel((carpet1,carpet2,screen,fibre))
    
The impedance is calculated using the Allard version of the surface impedance calculation::

    # Calculate normal surface impedance
    Z          = TMM_layup.impedance_allard(omega,kx=0,signal = False)
    
which gives the following result in accordance with [All2009]_ except the frequency unit.

  .. figure:: ./images/carpet_fibre_impedance.*
   :align: center
   :width: 70%   
   
   Normal impedance of layup from [All2009]_ as in figure 11.17
   
With the following command the normal and diffuse absorption can be caluclated ::

    # Calculate normal absorption
    alpha0     = TMM_layup.absorption(omega,kx=0.0,signal = False,allard=True)
    # Calculate diffure absorption
    alpha_diff = TMM_layup.absorption_diffuse(omega,theta_max=np.pi*78/180,theta_step = np.pi/180, signal = False, allard=True)
    
Note the ``allard = True`` keyword argument required to use the Allard method of the absorption calculation.
The figure reveals that there are better absorbers in the world.

  .. figure:: ./images/carpet_fibre_absorption.*
   :align: center
   :width: 70%   
   
   Normal and diffuse absorption of layup
   
 
Acoustic transmission design with poroelastic foam and rubber
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

A practical version double wall systems is a so called mass-spring system. 
Such layup consist of a soft foam covered by a heavy layer, e.g. rubber.
With the Allard theory poroelaxtic foams can be investigated in detail.
An example layup is shown in the following figure.

  .. figure:: ./images/alu_melamin_rubber_layup.*
   :align: center
   :width: 70%   
   
   Mass-spring layup

In this example the impact of different modelling approaches and assumptions is shown.
The code is given in :ref:`sec-TL-poro-examples`, please refer to the code for the material details.

The rubber and the Aluminium of the base plate are given as :class:`~pyva.properties.materialClasses.IsoMat`. ::

    # Isotropic materials
    alu     = matC.IsoMat(eta = 0.1)
    rubber  = matC.IsoMat(E=2.6e6,rho0=1200,nu=0.49,eta=0.00)

The melamin foam is given by ::

    # Melamin foam paramters
    melamin_vac = matC.IsoMat(E=3.0e5,rho0=12.0,nu=0.4,eta=0.1) # Frame in vaccuum    E6
    melamin = matC.PoroElasticMat(melamin_vac, \
                                flow_res = 30000., \
                                porosity = 0.99, \
                                tortuosity = 1.01, \
                                length_visc = 250.E-6, length_therm = 550.E-6)
                                
Every layer is given as plate, so that we can try different options::

    # plate properties
    alu_1mm        = stPC.PlateProp(0.001, alu)
    rubber_2mm     = stPC.PlateProp(0.002, rubber)
    # test the foam as solid
    foam_3cm_solid = stPC.PlateProp(0.03, melamin_vac) 

We define the melamin foam either as :class:`~pyva.infiniteLayers.PoroElasticLayer` or :class:`~pyva.infiniteLayers.SolidLayer` ::

    # Foam Layers
    iL_foam_3cm = iL.PoroElasticLayer(melamin, 0.03)
    iL_foam_3cm_solid = iL.SolidLayer(foam_3cm_solid)

whereas the aluminium plate and the rubber layer are given as :class:`~pyva.infiniteLayers.SolidLayer` or :class:`~pyva.infiniteLayers.ImperviousScreenLayer` ::

    # rubber and alu as solid- and screen layer
    iL_rubber_solid_2mm = iL.SolidLayer(rubber_2mm)
    iL_rubber_imper_2mm = iL.ImperviousScreenLayer(rubber_2mm)
    iL_alu_solid_1mm    = iL.SolidLayer(alu_1mm)
    iL_alu_imper_1mm    = iL.ImperviousScreenLayer(alu_1mm)
    
For the decoupling we create a super thin layer of low mass ::

    # Mass of Fluid as gap
    iL_nothing    = iL.MassLayer(1e-6,1.)

First we are interested in the impact of the solid or screen layer modelling ::

    # TMmpodel using the solid or impervious screen formulation 
    alu_melamin_rubber_solid = mds.TMmodel((iL_alu_solid_1mm,iL_foam_3cm,iL_rubber_solid_2mm))
    alu_melamin_rubber_imper = mds.TMmodel((iL_alu_imper_1mm,iL_foam_3cm,iL_rubber_imper_2mm))

and two other variations that decouple the alu plate of model the melamin foam as solid - neglecting the fluid waves. ::

    alu_melamin_rubber_decoup = mds.TMmodel((iL_alu_imper_1mm,iL_nothing,iL_foam_3cm,iL_rubber_imper_2mm))
    alu_melamin_rubber_foam_as_solid = mds.TMmodel((iL_alu_imper_1mm,iL_foam_3cm_solid,iL_rubber_imper_2mm))
    
Calculating the transmission loss of all options is done by::

    tau_imper = alu_melamin_rubber_imper.transmission_diffuse(omega,theta_max=78/180*np.pi,allard=True,signal=False)
    tau_solid = alu_melamin_rubber_solid.transmission_diffuse(omega,theta_max=78/180*np.pi,allard=True,signal=False)
    tau_decoup = alu_melamin_rubber_decoup.transmission_diffuse(omega,theta_max=78/180*np.pi,allard=True,signal=False)
    tau_foam_as_solid = alu_melamin_rubber_foam_as_solid.transmission_diffuse(omega,theta_max=78/180*np.pi,allard=True,signal=False)

The results are shown in the following figure.

  .. figure:: ./images/allard_DW_TL.*
   :align: center
   :width: 70%   
   
   Transmission loss of alu + mass-spring system using different modelling options
   
First, the difference between the versions modelling both skins as solid or screen is low.
Second, decoupling the alu from the foam leads to very different results. As a consequence this means that the 
connection iof the noise control treatment must be exactly known to get correct results. This is also true for 
experiments where the treatment must be carefully glued to the plate to get correct results.
Third, when the foam is modelled as solid the result is different, but the difference is lower than compared to the decoupling
effect.




   


  