#!/bin/bash
set -euo pipefail
root=${XML_POC_ROOT}
cd "$root"
export LD_LIBRARY_PATH="$root/dependencies/install/lib:/scratch/kohlbach/contrib/build/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OPENMS_TMPDIR="$root/tmp"
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1
(while test ! -f results/MATRICES_COMPLETE; do date -Is; cat /proc/loadavg; grep -E 'MemAvailable|SwapFree|SwapTotal' /proc/meminfo; grep -E '^pswp(in|out)' /proc/vmstat; sleep 15; done) >> results/resource-monitor.log &
monitor=$!
trap 'kill "$monitor" 2>/dev/null || true' EXIT
PYTHONPATH="$root" python3 staging/run_openmp_scaling.py --matrix results/benchmarks.json --threads 1,2,4,8,16,32,64,96 --output results/openmp-scaling.json > logs/openmp-runs.log 2>&1
python3 staging/environment.py after
touch results/MATRICES_COMPLETE
