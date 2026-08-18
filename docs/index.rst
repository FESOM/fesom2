.. fesom2 documentation master file, created by
   sphinx-quickstart on Sat Sep 28 22:37:42 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

FESOM2 documentation
====================

The Finite volumE Sea Ice-Ocean Model (FESOM2).

Multi-resolution ocean general circulation model that solves
the equations of motion describing the ocean and sea ice using
finite-volume methods on unstructured computational grids. The
model is developed and supported by researchers at the Alfred
Wegener Institute, Helmholtz Centre for Polar and Marine
Research (AWI), in Bremerhaven, Germany.

Authors
   -------

   Sergey Danilov, Dmitry Sidorenko, Nikolay Koldunov, Patrick Scholz, Qiang Wang, Thomas Rackow, Helge Goessling and Lorenzo Zampieri


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting_started/getting_started
   general_configuration/general_configuration
   ocean_configuration/ocean_configuration
   forcing_configuration
   output_configuration
   seaice_configuration
   icepack_in_fesom
   data_processing/data_processing
   geometry
   meshes/meshes
   main_equations
   vertical_discretization
   temporal_discretization
   spatial_discretization
   time_stepping_transport
   subcycling_instead_solver
   isoneutral_diffusion_triangular_prisms
   zreferences
   code_documentation/io_mesh_info
   developer_documentation/github_actions_test_suite
   developer_documentation/running_the_test_suite_offline
   

Proposed structure:

::

    Introduction
    Getting started
        TL;DR version for supported HPC systems
        Detailed steps of compiling and runing the code
        Ubuntu based Docker container (to get first impression of the model)
        Troubleshooting
    Looking at the results
        Notebooks that comes with the model
        pyfesom2 (short description with link to documentation)
    Tutorials
        Build model on Ubuntu (Video)
        Add new output variable (Video)
        First look at model output (Video)
    General configuration (namelist.config)
        Time stepping
        Restarts
        ALE options
        Mesh geometry and partitioning
    Ocean configuration (namelist.oce, namelist.dyn, namelist.tra, namelist.cvmix, namelist.transit)
        Ocean dynamics and GM/Redi choices
        Momentum advection and viscosity options
        Tracer advection/diffusion and initialisation
        Vertical mixing (CVMix) and transient tracers
    Sea ice configuration (namelist.ice)
        Ice dynamics
        Ice thermodynamics
    Atmospheric forcing (namelist.forcing)
    Output (namelist.io)
        Adding new output variable
    Meshes
        Mesh format
        Mesh generation?
    Partitioning
        MESSY
        Hierarchical partitioning
    Data pre/post processing
        Initial conditions
        Convert grid to netCDF that CDO understands
    Discretizations and Algorithms
    Coupling interfaces
        To atmosphere
        To ocean biogeochemistry?
    Example experiments
    FAQ
    History

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
