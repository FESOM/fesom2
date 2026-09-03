.. _gpu_and_dwarfs:

GPU support and compute dwarfs
******************************

Parts of FESOM2 carry OpenACC directives so that the most expensive kernels
can run on a GPU. The port is not complete and is not enabled in production
setups, so this page describes what exists rather than a supported workflow.

Building with OpenACC
=====================

OpenACC is switched on at configure time and is only wired up for the two
compilers that support the directives used here, Cray and NVHPC:

::

    cmake -DENABLE_OPENACC=ON -DGPU_COMPUTE_CAPABILITY=cc80 ..

The relevant options are:

- **ENABLE_OPENACC=OFF** compile the OpenACC code paths. The directives are inactive without it, and the model builds and runs exactly as a CPU code.
- **DISABLE_OPENACC_ATOMICS=ON** compile the kernels that would otherwise use ``atomic`` statements in a form that does not depend on the order in which threads accumulate. Atomics make the result depend on that order, so leaving this option at its default keeps GPU runs reproducible. Switch it off only when measuring the performance the atomics would buy.
- **GPU_COMPUTE_CAPABILITY=cc80** target architecture for the NVHPC compiler. It is combined with the CUDA version into ``GPU_FLAGS``, which defaults to ``cuda12.2,cc80``.

With the Cray compiler the build adds ``-hacc``. With NVHPC it adds
``-acc=verystrict`` together with the ``-gpu`` flags derived from the
options above. Note that the Cray CPU build explicitly passes ``-hnoacc``,
because that compiler otherwise activates OpenACC by itself.

What is ported
==============

OpenACC directives are present in eleven source files. They cover the sea-ice
dynamics and transport (``ice_EVP.F90``, ``ice_fct.F90``,
``ice_setup_step.F90``), the tracer advection chain
(``oce_adv_tra_driver.F90``, ``oce_adv_tra_hor.F90``, ``oce_adv_tra_ver.F90``,
``oce_adv_tra_fct.F90``, ``oce_ale_tracer.F90``, ``oce_tracer_mod.F90``) and
the halo exchange (``gen_halo_exchange.F90``). Everything else, including the
ocean dynamics and the whole I/O layer, runs on the host. Sections that had
to be restructured for the accelerator are guarded with
``#ifndef ENABLE_OPENACC``, so the CPU code path is preserved verbatim
alongside the ported one.

Compute dwarfs
==============

The ``dwarf`` directory holds two standalone mini-applications that build and
run a single kernel outside the full model. ``dwarf_ice`` exercises the
sea-ice dynamics and ``dwarf_tracer`` the tracer advection, which are the two
areas the GPU port targets. A dwarf is the practical way to work on a kernel:
it starts in seconds, needs no forcing data, and isolates the kernel from the
rest of the time loop, so a change in performance or in results can be
attributed to the kernel alone.

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
