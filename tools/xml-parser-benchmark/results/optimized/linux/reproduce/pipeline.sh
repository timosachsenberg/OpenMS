#!/bin/bash
set -euo pipefail
original=${XML_POC_ROOT}
root="$original/optimized"
cd "$root"
while ! test -f "$original/results/EXTRA_VALIDATION_COMPLETE"; do
 for marker in MATRICES_EXIT_STATUS METADATA_EXIT_STATUS; do
  if test -f "$original/results/$marker"; then
   code=$(cat "$original/results/$marker")
   if [[ "$code" != 0 ]]; then echo "Original $marker failed: $code"; exit "$code"; fi
  fi
 done
 sleep 15
done
bash staging/prepare-source.sh > logs/prepare-source.log 2>&1
bash staging/run-builds.sh > logs/build-pipeline.log 2>&1
taskset -c 0-95 bash staging/run-matrices.sh > logs/matrix-pipeline.log 2>&1
