.. _icepack_in_fesom:

Icepack sea ice configuration
*****************************

This section describes the implementation of the Icepack sie ice column physics package in the FESOM2 model. The scope of the following paragraphs is to provide a practical guide to users interested in detailed simulations of the sea ice system with FESOM2, and not to describe the scientific features of Icepack. A detailed description of the sea ice parameterizations here implemented can be found on the website of the `CICE Consortium <https://github.com/CICE-Consortium/Icepack/wiki/Icepack-Release-Table>`_, which maintains and continuously develops this package. 

.. attention::
   The Icepack implementation in FESOM2 is still in a testing phase and we cannot guarantee a bugfree code nor good scientific results.

.. note::
   To get more information regardng the implementation of Icepack in FESOM2, to report bugs, or to get advice regarding the model setup do not hesitate to open an issue on the FESOM2 GitHub repository or to contact Lorenzo Zampieri at lorenzo(dot)zampieri(at)awi(dot)de. 
   
   You are invited to update and develop further this documentation by pushing your changes to the `FESOM2 Documentation repository <https://github.com/FESOM/fesom2_doc>`_ on GitHub.

General information
===================

Icepack–the column physics package of the sea-ice model CICE–is a collection of physical parameterizations that account for thermodynamic and mechanic sub-grid processes not explicitly resolved by the hosting sea-ice model, in our case FESOM2. The modular implementation of Icepack allows the users to vary substantially the complexity of the sea-ice model, with the possibility of choosing between several schemes and a broad set of active and passive tracers that describe the sea-ice state. Icepack v1.2.1 has been implemented in FESOM2 and can be used as an alternative to the standard FESIM thermodynamic module. As the standard FESIM implementation, the Icepack column-physics subroutines run every ocean time step. All the Icepack variables are defined directly on the nodes of the FESOM2 mesh, ensuring an optimal consistency between the ocean and the sea-ice components of the model. The inclusion of Icepack in FESOM2 required a revision of the calling sequence within the sea-ice model, which now follows that of the CICE model as illustrated in :numref:`call_seq`.

.. _call_seq:
.. figure:: img/call_seq.png

   Schematic describing the calling sequences of the Standard FESOM2 and FESOM2-Icepack implementations.

Icepack is licensed for use through the CICE Consortium. Therefore, we encourage the FESOM2 userbase interested in the Icepack features to be aware of the `License <https://github.com/CICE-Consortium/Icepack/blob/master/LICENSE.pdf>`_ when working with this model configuration. We report here a disclaimer from the `Icepack website <https://github.com/CICE-Consortium/Icepack/wiki>`_.

.. important::  
   Icepack releases are “functional releases” in the sense that the code runs, does not crash, passes various tests, and requires further work to establish its scientific validity. In general, users are not encouraged to use any of the CICE Consortium’s model configurations to obtain “scientific” results. The test configurations are useful for model development, but sea ice models must be evaluated from a physical standpoint in a coupled system because simplified configurations do not necessarily represent what is actually happening in the fully coupled system that includes interactive ocean and atmosphere components.

How to cite
"""""""""""

The current Icepack version implemented in FESOM2 is Icepack 1.2.1. To acknowledge the development work behind the implementation of Icepack in FESOM2 please cite `Zampieri et al. (2021) <https://search.proquest.com/docview/2469422827?fromopenview=true&pq-origsite=gscholar>`_, part of which used to compile this documentation, and `Hunke et al. (2020) <https://zenodo.org/record/3712299#.Xvn3DPJS9TZ>`_, in addition to the usual FESOM2 papers.

Implementation
==============

The implementation of Icepack in FESOM2 is fully modular, meaning that the users are free to vary the configuration via namelist parameters. When Icepack is used, ``namelist.icepack`` controls all settings related to the sea ice subgrid parameterizations, thus overriding the content of ``namelist.ice``. The dynamics (EVP) and advection schemes are still controlled by the standard ``namelist.ice``. Below we describe some of the most important namelist parameters, while we recommend consulting the `official Icepack documentation <https://buildmedia.readthedocs.org/media/pdf/cice-consortium-icepack/icepack1.2.1/cice-consortium-icepack.pdf>`_ for a more comprehensive description.

Namelist section &env_nml
"""""""""""""""""""""""""

- **nicecat** Defines the number of sea ice thickness categories.
- **nfsdcat** Defines the number of categories of the floe size distribution. This parameter should be set to 1 as the floe size distribution has never been tested in FESOM2.
- **nicelyr** and **nsnwlyr** Defines the number of vertical layers in sea ice and snow.

.. attention::
   Increasing substantially the number of thickness classes and vertical layers can lead to numerical instabilities (very thin vertical layers), memory issues, very large output files, and finally to a substantial slow down of the model because of the high number of tracers that need to be advected.

Namelist section &grid_nml
""""""""""""""""""""""""""

- **kcatbound** Specifies which criteria is followed to discretize the Ice Thickness Distribution (ITD). Setting **kcatbound** equal to 0, 1, or 3 gives lower thickness boundaries for any number of thickness categories. Setting **kcatbound=2** corresponds to the World Meteorological Organization ITD classification, and it is compatible only with **nicecat=5,6,7**.

Namelist section &tracer_nml
""""""""""""""""""""""""""""

Logical parameters to specify parameterizations and passive tracers. Only **tr_pond_cesm** has been tested extensively.  

Namelist section &nml_list_icepack
""""""""""""""""""""""""""""""""""

It regulates the type, frequency, and precision of the output for Icepack variables. Most of the Icepack variables can be defined as average over the grid cell (e.g. **aice**: average sea ice area fraction – 2D variable), or separately for each thickness class (e.g. **aicen**: sea ice area fraction in each thickness class – 3D variable), with the ITD information saved as a vertical dimension in the netCDF file. At the moment, variables defined over multiple vertical layers are output in separated files. For example, in a model configuration with **n** sea ice vertical layers, activating the **qice** output stream will lead to **n** files where ``qice_i.fesom.yyyy.nc`` contains the sea ice enthalpy of the **i**-*th* vertical layer averaged over the ITD. Similarly, activating the **qicen** output stream will lead to **n** files where ``qicen_i.fesom.yyyy.nc`` contains the sea ice enthalpy of the **i**-*th* sea ice vertical layer for each thickness class.

Compatibility with FESOM2 configurations
""""""""""""""""""""""""""""""""""""""""

In `Zampieri et al. (2020) <https://search.proquest.com/docview/2469422827?fromopenview=true&pq-origsite=gscholar>`_ the model was run with linear free surfaces (**which_ALE=’linfs’**), and other ALE coordinates have not been tested. In principle, Icepack should be independent of the scheme used to solve the sea ice dynamics. However, at the moment only the standard EVP is supported, while the mEVP and aEVP still show some strange behaviors. We are working on solving this issue as well as on testing further setups, and we will update this document as soon as progress is made.

Compilation
===========

Compiling FESOM2 with Icepack is very easy if you are already used to the FESOM2 workflow. After cloning fesom2 from the GitHub repository, download the Icepack single column package:
::

   cd src/icepack_drivers/
   bash -l download_icepack.sh
The next step is to activate the Icepack flag in ``CMakeLists.txt`` by setting **USE_ICEPACK** from **OFF** to **ON**. At this point, you can proceed with the usual compilation via
::

   bash -l configure.sh   
The compilation of this FESOM2 version with the ESM Tools is not yet supported.

Running the model
=================

Running FESOM2 with Icepack is not different from the standard case. Make sure to add the ``namelist.icepack`` file to your ``work`` directory. Two diagnostic files are generated in addition to the standard ``fesom2.0.out``. ``icepack.diagnostics`` contains information about the Icepack configuration such as the value of some parameters, the tracers employed, and the boundaries of the ITD. ``icepack.errors`` possibly contains diagnostic information about errors in Icepack that can occur during the model run. Information about the running time are given in ``fesom2.0.out`` with the usual division in **dynamics**, **thermodynamics**, and **advection**.

The model output is saved in the result folder together with the standard ocean output. Note that outputting sea ice information using the standard FESIM variables (**a_ice**, **m_ice**, **m_snow**, etc.) is still possible also when using Icepack. These variables are consistent with the Icepack sea ice description (**a_ice** = **aice**, **m_ice** = **vice**, **m_snow** = **vsno**). An additional restart file is generated for Icepack, ``fesom.yyyy.icepack.restart.nc``, and it is written with the same frequency as ``fesom.yyyy.oce.restart.nc`` and ``fesom.yyyy.ice.restart.nc``.

.. attention::
   Restarting the model after changing the number of ice thickness classes, the vertical discretization of ice and/or snow, and the number of passive tracers is currently not possible. Also, changing the thermodynamic and melt pond schemes during the run is not recommended. In these cases consider a cold start and repeat your spinup run.

Code structure
==============

Icepack is a single column model and therefore its subroutines act on one grid cell. The Icepack code is downloaded from a separate repository (see instructions on how to compile the model) and is located in ``src/icepack_drivers/Icepack/columnphysics/``. To integrate this code in a host General Circulation Model (GCM), in our case FESOM2, additional instructions are needed to define an interface between the two systems and to drive the Icepack subroutines. This interface is contained in the ``src/icepack_drivers/icedrv_*.F90`` files, which are part of the FESOM2 repository, and will be briefly described in the following section.

Icepack drivers
"""""""""""""""

- ``icedrv_main.F90`` This file contains the main module of the Icepack drivers. All the variables are declared here, together with the interface of the subroutines contained in various submodules. If new variables or subroutines need to be added to the code, this is a good place to start. Try to maintain all the variables private to increase the modularity of the code, and use the transfer interface to exchange variables with FESOM2. 

- ``icedrv_set.F90`` This file contains few subroutines that initialize the model parameters by reading the Icepack namelists or alternatively by extracting default values from the Icepack package. Furthermore, ``icepack.diagnostics`` is written here, and the sea ice state is initialized in case of a cold start of the model.  

- ``icedrv_allocate.F90`` This file contains subroutines that allocate the Icepack variables declared in ``icedrv_main.F90``. 

- ``icedrv_init.F90`` This file contains subroutines that initialize the Icepack variables declared in ``icedrv_main.F90`` and allocated in ``icedrv_allocate.F90``.

- ``icedrv_step.F90`` This file contains few subroutines that describe the calling sequence of the sea ice model when Icepack is used in FESOM2.  

- ``icedrv_advection.F90`` This file contains few subroutines that advect the Icepack tracers. If new parameterization or options are explored, you should check if the relative tracers are advected properly.  

- ``icedrv_transfer.F90`` This file contains subroutines that describe the procedure to pass information between FESOM2 and Icepack.

- ``icedrv_io.F90`` This file contains subroutines that describe the I/O streams for the Icepack variables, including restart procedures. If new parameterization or options are explored, you should check if the relative tracers are restarted properly. 

- ``icedrv_kinds.F90`` This file declares some standard types for variable declarations. 

- ``icedrv_system.F90`` This file contains subroutines that handle model errors inside Icepack, possibly stopping the model run, and that output warning messages when appropriate.

- ``icedrv_constants.F90`` This file defines some constants that are used in the Icepack drivers.

Communication between Icepack and FESOM2
""""""""""""""""""""""""""""""""""""""""

The Icepack environment is separated from the rest of FESOM2 and consists of a single big module with multiple submodules. Almost all the variables are private and are not visible by the FESOM2 code. The variables exchange between Icepack and FESOM2 takes place through the passing subroutines ``fesom_to_icepack`` and ``icepack_to_fesom``.

Frequently asked questions
==========================

Should I use Icepack for my simulations?
""""""""""""""""""""""""""""""""""""""""

It depends on your scientific questions. Icepack might be a good option if you are interested in sea ice processes in polar regions. In principle, the employment of Icepack should not negatively affect the ocean state but could make FESOM2 slower.

Is FESOM2 slower when run with Icepack?
"""""""""""""""""""""""""""""""""""""""

Yes, the model integration is slower for two reasons: 1. The sea ice subgrid parameterizations are more complex compared to the standard FESIM. 2. Much more sea-ice tracers need to be advected. Overall, the sea ice component of FESOM2 becomes approximately four times slower with Icepack. Including additional output related to a more complex sea ice description can also contribute to deteriorating the model performances.    

Which EVP scheme should I use with Icepack?
""""""""""""""""""""""""""""""""""""""""""

In principle, Icepack should be independent of the scheme used to solve the sea ice dynamics. However, at the moment only the standard EVP is supported, while the mEVP and aEVP still exhibit some strange behaviors. We are working on solving this issue and we will update this document as soon as progress is made.

Can Icepack be configured as the standard FESIM?
""""""""""""""""""""""""""""""""""""""""""""""""

Yes, in principle it is possible to run Icepack with a single thickness class and with the 0-layer thermodynamics. However, the results obtained during the testing phase with this configuration were not very convincing and they seemed not compatible with the standard FESOM2 results. More investigations are needed to understand the cause of this behavior, which is likely related to a different implementation of the thermodynamic processes in the model.   

Can Icepack be used in coupled configurations?
""""""""""""""""""""""""""""""""""""""""""""""

No, at the moment FESOM2 with Icepack has not been coupled with atmospheric models. A coupling with OpenIFS is planned and might be released in the upcoming months.

Can Icepack be used with data assimilation?
"""""""""""""""""""""""""""""""""""""""""""

No, at the moment FESOM2 with Icepack has not been equipped with data assimilation capabilities. 
