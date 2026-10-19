#!/bin/bash
set -euo pipefail
root=${XML_POC_ROOT}
cd "$root"
while ! test -f results/METADATA_COMPLETE; do
 for marker in MATRICES_EXIT_STATUS METADATA_EXIT_STATUS; do
  if test -f "results/$marker"; then
   status=$(cat "results/$marker")
   if [[ "$status" != 0 ]]; then echo "$marker failed: $status"; exit "$status"; fi
  fi
 done
 sleep 15
done
export LD_LIBRARY_PATH="$root/dependencies/install/lib:/scratch/kohlbach/contrib/build/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cmake -S staging/validator-probe-source -B validator-probe-build -DOpenMS_DIR="$root/build" -DCMAKE_BUILD_TYPE=Release > logs/validator-probe-configure.log 2>&1
cmake --build validator-probe-build -j4 > logs/validator-probe-build.log 2>&1
python3 staging/validate-huge.py > logs/huge-schema-api.log 2>&1
touch results/EXTRA_VALIDATION_COMPLETE
