Extension Modules
==================

.. note::
	Extension Modules are parts of implementation codes which need to be edited for different uses of the algorithm. 
	A user are suggested to develop edit these modules for their own purpose. 
	

M file Functions
----------------

.. mat:automodule:: m_functions
.. mat:autofunction:: EModule_savegrainsfigure_withLabelRGB_ori
.. mat:autofunction:: EModule_initialdoubleVoronoidata2d
.. mat:autofunction:: EModule_kernel_widths2d
.. mat:autofunction:: EModule_kernel_widths2d_ABC

C file Functions
----------------

.. mat:automodule:: c_functions

.. _CModule_updatelevelset_ABC:
.. mat:function:: CEModule_updatelevelsetdata2d_ABC(presence,grains,ID,ori_in_deg,alpha,beta,angBrandon,option)

	* Update (threshold) the grain labels for ABC-group example 
	
	:param presence: :ref:`presence <presence>`
	:param grains: :ref:`grains <grains>`
	:param ID: grain ID
	:param ori_in_deg: orientation angels for grains. Here, it must be stated in `deg`.
	:param alpha: longer Gaussian kerenl width 
	:param beta: shorter Gaussian kerenl width 
	:param angBrandon: Brandon angle in degree.  For example, '30'. (Actually this parameter is unused in this function)
	:param option: choice of mobility option. If option=1, then a constant mobility is used. 
		If option=2, then a higher mobility for high angle graion boundary (misorientation angle larger than 10 degree).
		Mobility functions can be further generalized by changing the internal mobility functions defined in this c-file.  
	
	surfacetension::

		double surfacetension(double ori1, double ori2)
		{
		  int id1Type, id2Type; 
  
		  if(ori1< 5)
		    id1Type=0; // type A 
		  else if(ori1<20)
		    id1Type=1; 
		  else
			id1Type=2;
  
		  if(ori2< 5)
		    id2Type=0; // type A 
		  else if(ori2<20)
		    id2Type=1; 
		  else
			id2Type=2;
    
		  double st; 
		  double gamma1=1.0; 
		  double gamma2=1.5;
		  double gamma3=2.0; 
  
		  if(id1Type==id2Type)
		      st = gamma1; 
		  else if (id1Type== 2 || id2Type==2)
		      st = gamma2; // interaction A-C and B-C
		  else 
		      st = gamma3;  // interaction A-B 
 
		  return st;
		}
		
.. _CModule_updatelevelset_CVE:
.. mat:function:: CEModule_updatelevelsetdata2d_CVE(presence,grains,ID,ori_in_deg,alpha,beta,angBrandon,option,CVE)

	* Update (threshold) the grain labels when grain boundary energy is given by an external data 
	
	:param presence: :ref:`presence <presence>`
	:param grains: :ref:`grains <grains>`
	:param ID: grain ID
	:param ori_in_deg: orientation angels for grains. Here, it must be stated in `deg`.
	:param alpha: longer Gaussian kerenl width 
	:param beta: shorter Gaussian kerenl width 
	:param angBrandon: Brandon angle in degree.  For example, '30'. (Actually this parameter is unused in this function)
	:param option: choice of mobility option. If option=1, then a constant mobility is used. 
		If option=2, then a higher mobility for high angle graion boundary (misorientation angle larger than 10 degree).
	:param CVE: a string of the external data file 
	
	The external data file should look like::
	
		#*misorientation values* *grain boundary energy data* 
	
		0.0 0.1
		
		0.5 0.10458
		
		1.0 0.203463
	
		1.5 0.296963 

		2 0.385376
		... ...