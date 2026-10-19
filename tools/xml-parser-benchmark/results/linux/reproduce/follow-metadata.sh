#!/bin/bash
set -euo pipefail
root=${XML_POC_ROOT}
cd "$root"
while ! test -f results/MATRICES_COMPLETE; do
 if test -f results/MATRICES_EXIT_STATUS; then
  status=$(cat results/MATRICES_EXIT_STATUS)
  if [[ "$status" != 0 ]]; then echo "Primary matrix failed: $status"; exit "$status"; fi
 fi
 sleep 15
done
set +e
taskset -c 0-95 bash staging/run-metadata.sh
status=$?
printf '%s\n' "$status" > results/METADATA_EXIT_STATUS
exit "$status"
