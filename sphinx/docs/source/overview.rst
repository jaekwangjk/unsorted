Overview
========

* We divided the original source files into two types. :doc:`generalModule` are functions that implement the fundamental algorithm of the thersholding-dynamics. :doc:`extensionModule` are parts of implementation functions. User might have to edit these functions for their own problems.

* In both modules, some functions are orignally written in **C** for the speed of computation. However, these functions still need to be called in MATLAB environment. To this end, these functions must be compiled in advance using the *mex* command, which turn them as *MEX functions*. For example, 

.. code-block:: console

   mex CGModule_pgrow.c

* To edit this c-functions, a user may need to refer the source `MEX functions <https://www.mathworks.com/help/matlab/call-mex-file-functions.html>`_ to be familarize with its own API.

* :doc:`background` provides some additional knowledge and concepts that are useful how the code works. This includes data structure as well as mathematical theories. :doc:`examples` are driver functions that stands alone. 


