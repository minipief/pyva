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


  