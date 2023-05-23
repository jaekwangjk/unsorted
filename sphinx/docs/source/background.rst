Background
==========


Data structures
---------------

* For efficient implementation of the threshold-dynamics algorithms, following *MATLAB* data structures are used. They are widely used as inputs and outputs of modules.   

.. _grains:

**grains** 
 *grains* is  a *matlab cell array* with indexed grain data containers 
 It stores the following information: 
 
 grains{k,1} : List of grid points that belongs to grain k. Each point is indicated by 1d array
 
 grains{k,2} : Level set values. If the point is on the value is -1, otherwise +1 
 
 grains{k,3} : Convolution vals. of alpha kernel
 
 grains{k,4} : Convolution vals. of beta kernel
 
 .. note::
 	Each grain is defined by the list of grid indices that belong to (a) the grain itself 
 	and (b) grids that are near to outside the grain boundary. 
 	The buffering region (b) is where the results of the convolution kernel are stored.  
 	This advantageous for saving computational memeory, since we do not need to save the result
 	of kernel in the entire domain.  
 
 
.. _family:

**family** 

 Each family is a collection of grains that are located far apart. The concept of familiy
 allows to use a single kernel to all grains in the family, since their operations will not interfere each other. 

.. _presence:

**presence**
  
 *presence* is an dim1*dim2-by-3 cell array:
 Each element is related to information on grid data point. 
 Specifically, its k-th element stores the following information.
 
 presence{k,1} : Grains in the neighborhood of pixel at loc. index. 
 
 presence{k,2} : convolu vals. of those grains, in same order. 
 
 presence{k,3} : necessary information used to locate this pixel from each grains's data structure.
 


Global Variables
----------------

In the code, following global variables also are frequently used. 

.. _ori:

**ori** 

  *ori* is a 1d array that stores orientation of 2D grain respect to the reference state. 
  Its k-th element ori(k) is that of k-grain. This must be stated in radians.

.. _dims:

**dims** 
 *dims* is  a *cell array* that saves size of the domain, i.e. [nx,ny]




Level set and volume of fluid representation of grains 
------------------------------------------------------


Family 
------


