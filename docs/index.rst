.. fesom2 documentation master file, created by
   sphinx-quickstart on Sat Sep 28 22:37:42 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

FESOM2 documentation
====================

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
    General configuration (namelist.configure)
        Time stepping
        Restarts
        ALE options
        Mesh geometry and partitioning
    Ocean configuration (namelist.oce)
        Ocean dynamics
        Ocean tracers
        Adding a new tracer
        Initial conditions
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
