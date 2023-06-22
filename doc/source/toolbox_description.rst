Toolbox description
===================

This section introduces the modules and sub-packages of pyva from the [All2009] perspective.
It presents each module in the context of its [All2009] purpose. 

.. figure:: ./images/pyva_packages.*
   :align: center
   :width: 100%
   
   Overview of pyva sub-packages and modules 

The toolbox is organised in accordance with the setup of vibroacoustic systems. 
The top level module of the pyva package is the :mod:`pyva.models` module providing methods and classes for the joined [All2009] 
of deterministic and random systems. 

All sub-packages and their included modules populate the database to describe the systems or to provide classes and
methods to simulate the dynamics and coupling of systems. This includes the physical models but also data structures
to efficiently handle system matrices, vectors and degrees of freedom. 

This section gives an overview about the functionality and features of the implemented modules and classes. 
For details of classes, methods and attributes please refer to the API section :ref:`sec-pyva_API`, for the application 
of the classes in a [All2009] context you should refer to the :ref:`sec-creating_models` section of the :ref:`sec-user_guide` 


.. toctree::
   :maxdepth: 4

   model_classes
   property_classes
   system_classes
   coupling_classes
   geometry_classes
   data_structures


    
    

 
 

    
    








