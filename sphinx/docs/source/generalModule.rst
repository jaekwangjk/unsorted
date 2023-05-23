General Modules
===============


.. note::
	General Modules are functions that implement the fundamental algorithms of the thersholding-dynamics. It is unlikley that a user edits these modules, unless one develops a new thresholding rule. 

.. mat:automodule:: m_functions

M file Functions
----------------

.. mat:autofunction:: GModule_grains2families2d
.. mat:autofunction:: GModule_convolvefamily2d
.. mat:autofunction:: GModule_Thresholding

C file Functions
----------------

.. mat:automodule:: c_functions

.. mat:function:: CGModule_pgrow( x, y, Rgrains, WORKSPACE, work_x, work_y, work_bx, work_by, work_candidate_x, work_candidate_y)

	* Takes the grid points occuppied by a grain and dilates grains by *Rgrains* pixels outward with the periodic boundary conditions. These are the points where the level set value of the grain is positive value.
	
	In MATLAB interface, the function can be used like::
	
		[x2,y2] = CGModule_pgrow3(int32(x),int32(y),Rgrains,WORKSPACE,work_x,work_y,work_bx,work_by,work_candidate_x,work_candidate_y);
		
	:param x: List of indices in x grid occuppied by the grain. 
	:param y: List of indices in y grid occuppied by the grain.
	:param Rgrain: Width of periodic growth (window size) 
	:param WORKSPACE: pre-declared variables for memoery allocation 
	:param work_x: pre-declared variables for memoery allocation. Same for other similarly named variables 
	:return x2 y2: new x (and y) coord to list of x (and y) coords occupied by the grain 

.. mat:function:: CGModule_get_nhd_grains2d(grains,dims(1)*dims(2))

	* Identify :ref:`presence <presence>` from :ref:`grains <grains>`. 
	
	In MATLAB interface, the function can be used like::
	
		presence = CGModule_get_nhd_grains2d(grains,dims(1)*dims(2));
	
	:param grains: :ref:`grains <grains>`
	:param dims: 1d array dims
	:return presence: data structure presence 
	
.. mat:function:: CGModule_loc_levset_to_volfluid(int32(x),int32(y),val,Z)

	* Convert the level represenation of grain to characteristic function (volume of fluid) on the union of grains from its level-set representation 
	
	This function is used the matlab function *Module_convolvefamily2d*
	
	In MATLAB interface, the function can be used like::
	
		vf = GModule_loc_levset_to_volfluid(int32(x),int32(y),val,Z);


