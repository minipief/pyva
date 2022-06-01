Geometry package
================

Currently the geometry classes is implemented very simplistic. The only module
is the mod:`pyva.geometry.meshClasses` module. 

The required import is ::

    import pyva.geometry.meshClasses as meshC

Mesh classes
------------

The only yet implemented mesh is a regular two-dimensional meshClass that provides methods for 
regular meshes, used for example for the discrete radiation stiffness. 

The mesh is created with ::

    # Dimensions
    Lx = 1.2
    Ly = 1.5
    Nx = 3
    Ny = 7

    my_mesh = meshC.RegMesh2D(0., 0., Lx, Ly, Nx, Ny)
    
There are many properties that can be found in the API description of class:`pyva.geometry.meshClasses.RegMesh2D`.
An important method that gives the X and Y positions is the nodes method ::

    X,Y = my_mesh.nodes()
    >>> X
    array([[0.  , 0.24, 0.48, 0.72, 0.96, 1.2 ],
           [0.  , 0.24, 0.48, 0.72, 0.96, 1.2 ],
           [0.  , 0.24, 0.48, 0.72, 0.96, 1.2 ],
           [0.  , 0.24, 0.48, 0.72, 0.96, 1.2 ],
           [0.  , 0.24, 0.48, 0.72, 0.96, 1.2 ],
           [0.  , 0.24, 0.48, 0.72, 0.96, 1.2 ],
           [0.  , 0.24, 0.48, 0.72, 0.96, 1.2 ]])
       
    >>> Y
    array([[0.  , 0.  , 0.  , 0.  , 0.  , 0.  ],
           [0.25, 0.25, 0.25, 0.25, 0.25, 0.25],
           [0.5 , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 ],
           [0.75, 0.75, 0.75, 0.75, 0.75, 0.75],
           [1.  , 1.  , 1.  , 1.  , 1.  , 1.  ],
           [1.25, 1.25, 1.25, 1.25, 1.25, 1.25],
           [1.5 , 1.5 , 1.5 , 1.5 , 1.5 , 1.5 ]])

The radiation efficiency required the distance between all nodes. This means for a mesh of N nodes N*(N+1)/2 distances
of the upper triangular matrix because the distances are symmetric. 
The good thing with regular meshes is that that most distances have similar values. The distance method has two 
outputs ::

    dist,index = my_mesh.distance()
    
The number of mesh nodes is ::

    >>> my_mesh.Nmesh
    42
    
And the size of the distance is also 42 ::

    >>> dist.size
    42
    >>> dist
    array([0.        , 0.24      , 0.48      , 0.72      , 0.96      ,
           1.2       , 0.25      , 0.34655446, 0.54120237, 0.762168  ,
           0.99201816, 1.2257651 , 0.5       , 0.554617  , 0.6931089 ,
           0.8765843 , 1.0824047 , 1.3000001 , 0.75      , 0.7874643 ,
           0.89044935, 1.0396634 , 1.2182364 , 1.4150972 , 1.        ,
           1.0283968 , 1.109234  , 1.2322338 , 1.3862178 , 1.56205   ,
           1.25      , 1.2728314 , 1.3389921 , 1.4425324 , 1.5761029 ,
           1.7327724 , 1.5       , 1.5190787 , 1.5749286 , 1.6638509 ,
           1.7808987 , 1.9209373 ], dtype=float32)

With the index we can reconstruct all distances ::

    >>> index.size
    903
       
    >>> index
    array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
           17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
           34, 35, 36, 37, 38, 39, 40, 41,  0,  1,  2,  3,  4,  7,  6,  7,  8,
            9, 10, 13, 12, 13, 14, 15, 16, 19, 18, 19, 20, 21, 22, 25, 24, 25,
           26, 27, 28, 31, 30, 31, 32, 33, 34, 37, 36, 37, 38, 39, 40,  0,  1 ... ])

The upper trangle can be reconstructed by ::

    dist_triu = dist[index]
    
Thus, with the indexing trick the number of data is tremendously reduces. As the distance 
is used in many stiffness matrix formula this increases the calculation speed also.
Finally, the mesh can be plotted with ::

    my_mesh.plot3D(1)
    
.. figure:: ./images/mesh_plot.*
   :align: center
   :width: 70%
   
   3D plot of RegMesh2D 


    




