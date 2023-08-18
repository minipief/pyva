Models module
=============

The model classes are part of the models :mod:`pyva.models` module. 

.. _fig-models:
    
.. figure:: ./images/models.*
   :align: center
   :width: 90%
   
   Global setup of a vibroacoustic model.
   
Finite Element Model (FEM)
--------------------------

In general, finite element models are the result of specific numerical method that models a structure or a fluid 
by subdividing the system into several finite elements. These methods are subject of many books, projects and 
research. In [Bat1982]_ a good overview about the method is for example given.

The current implementation of finite element models is very basic. Currently, only analytical mode shapes are mapped to 
two dimensional meshes in order to provide the modal database for the hybrid models.

A full scale FE-model comes with a mesh and coordinate system definition and the related mass and stiffness matrices that 
are required to determine the dynamic stiffness matrix for every frequency or to calculate the modal solutions.
However, there is the intention to include this into later versions of pyva preferably from other packages that are already available.

Vibro Acoustic Model VAmodel
----------------------------

The VAmodel class :class:`pyva.models.VAmodel` describes a deterministic system by a stiffness matrix for a given set of
degrees of freedom.

One typical example is given in equation :eq:`dynamicStiffnessEOM_VA`.

.. math:: 
    :label: dynamicStiffnessEOM_VA
	
	\begin{bmatrix} 
	D_{11}(\omega) & D_{12}(\omega) & \cdots & D_{1M}(\omega) \\
	D_{21}(\omega) & D_{22}(\omega) &        & D_{2M}(\omega) \\
	\vdots & \vdots & \ddots & \vdots  \\
	D_{N1}(\omega) & D_{N2}(\omega) & \cdots & D_{NM}(\omega) \\
	\end{bmatrix}  
	\begin{Bmatrix} q_1(\omega) \\ q_2(\omega)  \\ \vdots \\ q_N(\omega) \end{Bmatrix} =
	\begin{Bmatrix} F_1(\omega) \\ F_2(\omega)  \\ \vdots \\ F_N(\omega) \end{Bmatrix}
	
:math:`D_{ij}(\omega)` are the coefficients of the dynamic stiffness matrix, 
:math:`q_{i}(\omega)` the qeneralised displacement degrees of freedom and  
:math:`F_{i}(\omega)` the generalized forces.

In FEmodels the mass and stiffness matrices are not frequency dependent.
In contrast to that, the dynamic stiffness matrix is frequency dependent matrix.
Naturally, the VAmodel class is an extension of the DynamicMatrix class with additional results and loads.
Many system classes from the systems module :mod:`pyva.systems` provide methods to generate entries into the matrix in 
order to model such systems as part of a VAmodel.

For example the mass and stiffness matrices of an FE-model matrix equation can be converted into a dynamic stiffness representation by:

.. math:: 
    :label: dynamicStiffnessEOM_VA
	
	\left( \begin{bmatrix} K \end{bmatrix} - \omega^2  \begin{bmatrix} M \end{bmatrix} \right) 
    \begin{Bmatrix} \bm{q} \end{Bmatrix} =  
    \begin{bmatrix} \bm{D} \end{bmatrix}
    \begin{Bmatrix} \bm{q} \end{Bmatrix} =  
 	\begin{Bmatrix} \bm{F} \end{Bmatrix}



HybridModel
-----------

The HybridModel class :class:`pyva.models.HybridModel` serves as a model container for SEA- and FEM-subsystems. 
The subsystems can be connected via junctions and excited by specific loads.
The class attributes include SEA- and FEM-matrices that are required to calculate the 
SEA and FEM response of the system.

The SEA part is given by the SEA matrix:

.. math::
    :label: SEA-matrix

    \omega
	\begin{bmatrix} 
	n_1 \sum_{n=1}^{N} \eta_{1n} & -n_2(\omega)\eta_{21}        & \cdots & -n_N(\omega)\eta_{N1} \\
	-n_1(\omega)\eta_{12}        & n_2 \sum_{n=1}^{N} \eta_{2n} &         & \vdots \\
    \vdots                       &                              & \ddots  &        \\
	-n_1(\omega)\eta_{1N}        & \cdots                       &         & n_N \sum_{n=1}^{N} \eta_{Nn} 
    \end{bmatrix}  
	\begin{Bmatrix} \frac{E_1}{n_1(\omega)} \\ \frac{E_2}{n_2(\omega)}  \\ \vdots \\ \frac{E_N}{n_N(\omega)} \end{Bmatrix} =
	\begin{Bmatrix} \Pi_1(\omega) \\ \Pi_2(\omega)  \\ \vdots \\ \Pi_N(\omega) \end{Bmatrix}
	
:math:`\eta_{ij}(\omega)` are the coupling loss factors, 
:math:`n_{i}(\omega)` the modal density of the *wave* system, 
:math:`E_{i}(\omega)` the energy of the *wave* system and  
:math:`\Pi_{i}(\omega)` the input power.   

Please note that only SEA Systems are obligatory. FEM-subsystems are not required, but can be included. 

Transfer matrix model
---------------------

Especially the  simulation of flat noise control treatment can be performed efficiently using the 
transfer matrix method. The :class:`pyva.models.TMmodel` class provides methods to deal with system configurations that can be 
treated by transfer matrices. 
Obviously, this is restricted to sequential arrangements, for example a series of tubes or infinite layers. 
Furthermore, systems with complex degrees of freedom (for example the fluid and solid phase of porous layer) cannot yet be modelled by 
the transfer matrix method.  

  
