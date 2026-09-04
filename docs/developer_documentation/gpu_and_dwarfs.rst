.. _gpu_and_dwarfs:

GPU support and compute dwarfs
******************************

Parts of FESOM2 carry OpenACC directives so that the most expensive kernels
can run on a GPU. The port is not complete, so this page describes what
exists rather than a supported workflow.

Building with OpenACC
=====================

OpenACC is switched on at configure time. Only the Cray and NVHPC compiler
paths set the corresponding compiler flags:

::

    cmake -DENABLE_OPENACC=ON -DGPU_COMPUTE_CAPABILITY=cc80 ..

The relevant options are:

- **ENABLE_OPENACC=OFF** compile the OpenACC code paths. The directives are inactive without it, and the model builds and runs exactly as a CPU code.
- **DISABLE_OPENACC_ATOMICS=ON** compile the kernels that would otherwise use ``atomic`` statements in a form that does not depend on the order in which threads accumulate. Atomics make the result depend on that order, so leaving this option at its default keeps GPU runs reproducible. Switch it off only when measuring the performance the atomics would buy.
- **GPU_COMPUTE_CAPABILITY=cc80** intended target architecture, combined with the CUDA version into a ``GPU_FLAGS`` cache variable. Both are currently unused: the NVHPC flags in ``src/CMakeLists.txt`` name ``cc80`` and ``cuda12.2`` directly, so changing these options has no effect and a different architecture has to be set in the build file itself.

With the Cray compiler the build adds ``-hacc``. With NVHPC it adds
``-acc=verystrict`` together with the ``-gpu`` flags derived from the
options above. Note that the Cray CPU build explicitly passes ``-hnoacc``,
because that compiler otherwise activates OpenACC by itself.

What is ported
==============

OpenACC directives cover the sea-ice dynamics and transport
(``ice_EVP.F90``, ``ice_fct.F90``, ``ice_setup_step.F90``), the tracer
advection chain (``oce_adv_tra_driver.F90``, ``oce_adv_tra_hor.F90``,
``oce_adv_tra_ver.F90``, ``oce_adv_tra_fct.F90``, ``oce_ale_tracer.F90``,
``oce_tracer_mod.F90``), the halo exchange (``gen_halo_exchange.F90``) and
the data movement around the time loop in ``fesom_module.F90``. The ocean
dynamics and the I/O layer have no directives. Sections that had to be
restructured for the accelerator are guarded with ``#ifndef ENABLE_OPENACC``,
so the CPU formulation stays in the file next to the ported one. To find the
current extent of the port, search the sources for ``acc`` directives rather
than relying on this list.

Compute dwarfs
==============

The ``dwarf`` directory holds two standalone mini-applications that build and
run a single kernel outside the full model. ``dwarf_ice`` exercises the
sea-ice dynamics and ``dwarf_tracer`` the tracer advection, which are the two
areas the GPU port targets. A dwarf isolates the kernel from the rest of the
time loop, so a change in performance or in results can be attributed to the
kernel alone.

Each dwarf directory contains a link script and a ``work`` directory to run
in:

::

    cd dwarf/dwarf_ice
    ./dwarf_linkfiles.sh
    ./configure.sh
    cd work

``dwarf_linkfiles.sh`` symlinks the environment files and ``configure.sh``
from the repository root, then builds a local ``src`` directory whose entries
are symlinks into the model's own ``src``. A dwarf therefore always compiles
the current state of the model code, and editing a kernel in ``src`` is
picked up by the dwarf without any copying. The driver program that replaces
the full time loop lives in the dwarf's ``dwarf_ini`` directory, and ``work``
holds the job script to run from.
