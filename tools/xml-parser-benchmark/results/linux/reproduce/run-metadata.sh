#!/bin/bash
set -euo pipefail
root=${XML_POC_ROOT}
cd "$root"
test -f results/MATRICES_COMPLETE
export LD_LIBRARY_PATH="$root/dependencies/install/lib:/scratch/kohlbach/contrib/build/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OPENMS_TMPDIR="$root/tmp"
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1
cmake -S staging/metadata-benchmark-source -B benchmark-metadata-build -DOpenMS_DIR="$root/build" -DCMAKE_BUILD_TYPE=Release > logs/metadata-benchmark-configure.log 2>&1
cmake --build benchmark-metadata-build -j8 > logs/metadata-benchmark-build.log 2>&1
for backend in xerces libxml2 expat pugixml; do
 LD_LIBRARY_PATH="$root/snapshots/$backend/lib:$LD_LIBRARY_PATH" ctest --test-dir benchmark-metadata-build --output-on-failure > "results/$backend-metadata-benchmark-tests.log" 2>&1
done
python3 staging/environment.py metadata-before
(while test ! -f "$root/results/METADATA_COMPLETE"; do date -Is; cat /proc/loadavg; grep -E 'MemAvailable|SwapFree|SwapTotal' /proc/meminfo; grep -E '^pswp(in|out)' /proc/vmstat; sleep 15; done) > results/metadata-resource-monitor.log &
monitor=$!
trap 'kill "$monitor" 2>/dev/null || true' EXIT
python3 staging/metadata-benchmark-source/run_benchmarks.py --binary "$root/benchmark-metadata-build/xml_parser_benchmark" --backend "xerces=$root/snapshots/xerces/lib" --backend "libxml2=$root/snapshots/libxml2/lib" --backend "expat=$root/snapshots/expat/lib" --backend "pugixml=$root/snapshots/pugixml/lib" --input "$root/data/large.mzML" --input "$root/data/medium-zlib.mzML" --workers 1,2,4,8,16,32,64,96 --modes indexed-metadata,contexts-metadata --repetitions 5 --queries 2048 --timeout 600 --output results/metadata-scaling.json > logs/metadata-runs.log 2>&1
python3 staging/environment.py metadata-after
touch results/METADATA_COMPLETE
