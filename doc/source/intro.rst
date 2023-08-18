About pyva
==========

Vibroacoustic simulation is a challenging and demanding task. 
Students in this field of acoustics must get used to several methods and tools. 
A full frequency vibroacoustic simulation requires very different tools
ranging from deterministic approaches, for example finite element methods (FEM), to random methods as the statistical energy analysis (SEA)
In the mid frequency range deterministic and random approaches must be used simultaneously.
This is called the hybrid FEM/SEA-method.

Students (and small acoustic consultancies) usually cannot afford the licence costs of typical 
vibroacoustic  simulation software.  Thus, the **pyva** toolbox aims at providing a software framework
to apply all the above mentioned methods in a python environment for free.
The idea is to allow students to exercise vibroacoustic  simulation using a free software without the need to implement
everything from scratch.

Objectives
----------

So what are the main goals of the pyva package?
The main goals are ...

- to establish a data structure that deals with typical data of vibroacoustic simulation as transfer matrices, 
  system matrices, signals, shape function etc.
- to implement an API for SEA-methods that can be easily used and extended
- to implement an API for hybrid FEM/SEA-methods
- to implement an API for acoustic and elastodynamic material models
- to implement system descriptions for deterministic and random systems
- to provide an open framework to allow further extensions

Modules
-------

.. figure:: ./images/pyva_packages.*
   :align: center
   :width: 100%
   
   Overview of pyva sub-packages and modules 

The pyva package has several modules in a hierarchical way. 

The main module
+++++++++++++++

The top-level module :mod:`pyva.models` consists of classes that describe
a deterministic, random or hybrid model.

Those are classes are:

- Deterministic Finite Element Models :class:`pyva.models.FEM`
- Hybrid Models (that may only contain SEA Systems) :class:`pyva.models.HybridModel`
- Vibro Acoustic Models :class:`pyva.models.VAmodel`

The transfermatrix class is a deterministic model class but with the specific task
of multiple and infinite layer simulation.

- Transfermatrix Models :class:`pyva.models.TMmodel`

Subpackages
+++++++++++

The submodules have the following purposes:

1. To describe and model the dynamics of systems
2. To populate the database for system properties (materials, properties, contact dynamics)
3. To deal efficiently with the required system matrices and state vectors

The first task is performed by specific system classes in the modules of :mod:`pyva.systems` subpackage,
the second by the modules of the :mod:`pyva.properties` subpackage and the latter by the 
modules of the :mod:`pyva.data` subpackage.






