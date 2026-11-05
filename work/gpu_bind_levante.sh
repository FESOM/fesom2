#!/bin/bash
# gpu_bind_levante.sh <cmd...>: with all GPUs of the node visible (job_gpu_levante_p2p), select the GPU on this
# rank's NUMA node for OpenACC (ACC_DEVICE_NUM). Falls back to the local rank id.
cpu=$(taskset -cp $$ | cut -d: -f2 | tr ',' '\n' | head -1 | cut -d- -f1 | tr -d ' ')
mynuma=$(basename /sys/devices/system/cpu/cpu$cpu/node* | sed 's/node//')
dev=$SLURM_LOCALID; i=0
for bus in $(nvidia-smi --query-gpu=pci.bus_id --format=csv,noheader); do
  b=$(echo $bus | tr 'A-F' 'a-f' | sed 's/^0000//'); f=$(ls -d /sys/bus/pci/devices/*$b 2>/dev/null | head -1)
  [ -n "$f" ] && [ "$(cat $f/numa_node)" = "$mynuma" ] && dev=$i
  i=$((i+1))
done
export ACC_DEVICE_NUM=$dev
echo "rank $SLURM_PROCID cpu $cpu numa $mynuma -> ACC_DEVICE_NUM=$dev" >&2
exec "$@"
