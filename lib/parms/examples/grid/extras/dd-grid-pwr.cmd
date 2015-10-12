#!/usr/bin/tcsh
#@ initialdir = /homes/r31/lizhongz/parms/examples/grid
#@resources = ConsumableMemory(1024)
#@blocking  = unlimited
#@requirements = (Feature == "Fastp655")
#@ output       = test
#@ error        = test.err
#@ job_type     = parallel
#@ wall_clock_limit = 1:00:00
#@ network.MPI = css0,shared,US
#@ node = 4
#@ total_tasks = 4
#@ node_usage = shared
#@ queue
poe ./dd-grid.ex
