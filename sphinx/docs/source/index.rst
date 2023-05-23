.. Lumache documentation master file, created by
   sphinx-quickstart on Mon Mar 27 16:48:06 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation of statisticsGBM
======================================================

**statisticsGBM** is a *MATLAB*-based *threshold-dynamics* code for simulating the motion of grain boundaries during
grain growth.   
The original code was first published in `statsticsGBM <https://github.com/tiagosalvador/statisticsGBM>`_ 
together with `Ref[1] <https://arxiv.org/abs/1907.11574>`_, written by *Tiago Salvador and Selim Esedoglu* in 2019.  
The algorithm is a recent modification of threshold-dynamics of `Ref[2] <https://doi.org/10.1002/cpa.21527>`_
that allows anisotropic grain boundary mobilities without significant complications in algorithm. 
  
In this document, we discuss its modified version that is used for the subsequent study on influence on grain
boundary character anisotropy and grain statistics. 
The document is structured as follow.  

.. toctree::
   overview
   background
   generalModule
   extensionModule
   examples
   :maxdepth: 2
   :caption: Contents:

Start with :doc:`overview` for general introduction on of the source files and documentation. 

