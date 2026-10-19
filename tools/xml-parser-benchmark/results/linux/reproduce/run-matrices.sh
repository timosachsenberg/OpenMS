#!/bin/bash
set -euo pipefail
root=${XML_POC_ROOT}
cd "$root"
test -f results/ALL_BUILDS_COMPLETE
export LD_LIBRARY_PATH="$root/dependencies/install/lib:/scratch/kohlbach/contrib/build/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OPENMS_TMPDIR="$root/tmp"
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1
mkdir -p "$root/results"
python3 staging/verify-inputs.py
python3 staging/environment.py before
available=$(awk '/MemAvailable:/ {print $2}' /proc/meminfo)
if (( available < 600000000 )); then echo "Less than600GB memory available before96worker tests: ${available}KiB"; exit 1; fi
(while test ! -f "$root/results/MATRICES_COMPLETE"; do date -Is; cat /proc/loadavg; grep -E 'MemAvailable|SwapFree|SwapTotal' /proc/meminfo; grep -E '^pswp(in|out)' /proc/vmstat; sleep 15; done) > "$root/results/resource-monitor.log" &
monitor=$!
trap 'kill "$monitor" 2>/dev/null || true' EXIT
backends=(--backend "xerces=$root/snapshots/xerces/lib" --backend "libxml2=$root/snapshots/libxml2/lib" --backend "expat=$root/snapshots/expat/lib" --backend "pugixml=$root/snapshots/pugixml/lib")
inputs=(--input "$root/data/small.mzML" --input "$root/data/medium.mzML" --input "$root/data/large.mzML" --input "$root/data/medium-zlib.mzML" --input "$root/data/medium-zlib.mzML.gz" --input "$root/source/share/OpenMS/examples/FRACTIONS/BSA1_F1.mzML")
workers=1,2,4,8,16,32,64,96
python3 run_benchmarks.py --binary "$root/benchmark-build/xml_parser_benchmark" "${backends[@]}" "${inputs[@]}" --workers "$workers" --modes stream,indexed,index-open --repetitions 5 --queries 2048 --timeout 600 --output results/benchmarks.json > logs/benchmark-runs.log 2>&1
python3 run_benchmarks.py --binary "$root/benchmark-build/xml_parser_benchmark" "${backends[@]}" --input "$root/data/large.mzML" --input "$root/data/medium-zlib.mzML" --workers "$workers" --modes load,contexts,idle --repetitions 5 --queries 2048 --timeout 600 --output results/memory-scaling.json > logs/memory-runs.log 2>&1
PYTHONPATH="$root" python3 staging/run_openmp_scaling.py --matrix results/benchmarks.json --threads "$workers" --output results/openmp-scaling.json > logs/openmp-runs.log 2>&1
python3 staging/environment.py after
touch results/MATRICES_COMPLETE
