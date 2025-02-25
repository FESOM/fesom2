==============================================
Running FESOM's CI Jobs without GitHub Actions
==============================================

To run and fine-tune the FESOM CI jobs, you can use docker/singularity to launch the same
container that the jobs run in.

Example on Albedo
-----------------

.. code-block::

  $ module load apptainer
  Module for apptainer version 1.1.7 loaded
  $ git clone https://github.com/FESOM/fesom2.git
  Cloning into 'fesom2'...
  remote: Enumerating objects: 26648, done.
  remote: Counting objects: 100% (400/400), done.
  remote: Compressing objects: 100% (179/179), done.
  remote: Total 26648 (delta 260), reused 247 (delta 204), pack-reused 26248 (from 3)
  Receiving objects: 100% (26648/26648), 135.92 MiB | 41.51 MiB/s, done.
  Resolving deltas: 100% (19037/19037), done.
  Updating files: 100% (688/688), done.
  # In the next line, replace computing.computing with your SLURM account
  $ srun -A computing.computing singularity pull fesom2-ci.sif docker://ghcr.io/fesom/fesom2_docker:fesom2_test_refactoring-master
  INFO:    Converting OCI blobs to SIF format
  INFO:    Starting build...
  2025/02/14 13:59:23  info unpack layer: sha256:31e907dcc94a592a57796786399eb004dcbba714389fa615f5efa05a91316356
  2025/02/14 13:59:24  info unpack layer: sha256:d63063514200d35fd93aabad3af880cadda8716a3c35f86c66d6e93e0fcda303
  2025/02/14 13:59:34  warn rootless{usr/lib/x86_64-linux-gnu/gstreamer1.0/gstreamer-1.0/gst-ptp-helper} ignoring (usually) harmless EPERM on setxattr "security.capability"
  2025/02/14 13:59:55  info unpack layer: sha256:6fe2fecbd7fe377d431f3ffec826b32068530fdb7f127e4fc1410233a0a0059c
  2025/02/14 13:59:55  info unpack layer: sha256:b575551b9253c2a15679059d003335d83a63bb42a182b10c2415d4f7e4dfc516
  2025/02/14 14:00:05  info unpack layer: sha256:45d449d7d266f873dcfe8ded709f9205ba2b7bfe723406064f746da70149126d
  2025/02/14 14:00:05  info unpack layer: sha256:4f4fb700ef54461cfa02571ae0db9a0dc1e0cdb5577484a6d75e68dc38e8acc1
  2025/02/14 14:00:05  info unpack layer: sha256:f8cf99cd14bc14cbf0e4fed16961c0d5ecc48a47dccf89be8e8079a4129a59af
  2025/02/14 14:00:08  info unpack layer: sha256:4f4fb700ef54461cfa02571ae0db9a0dc1e0cdb5577484a6d75e68dc38e8acc1
  2025/02/14 14:00:08  info unpack layer: sha256:91697b90f8ddf44ede09eb97fd8cf9c8fa6f9b4ed00c5671a481683a95f3bed8
  INFO:    Creating SIF file...
  $ singularity run -B $(pwd)/fesom2:/fesom/fesom2 fesom2-ci.sif /bin/bash

This gets you into the container! Now, you can follow the rest of the steps here: https://github.com/FESOM/FESOM2_Docker/blob/master/fesom2_test/Readme.md

From inside the container
-------------------------

.. code-block::

  $ cd /fesom/fesom2
  $ bash -l configure.sh ubuntu  # Compile the model (ubuntu settings used since the container is ubuntu based)

Prepare the run with whatever test you want to run:

.. code-block::

  $ mkrun pi test_pi -m docker

Now, run the test simulation:

.. code-block::

  $ cd work_pi
  $ ./job_docker_new

