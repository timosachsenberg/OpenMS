#!/bin/bash
root=${XML_POC_ROOT}
cd "$root" || exit 1
taskset -c 0-95 bash staging/run-matrices.sh
status=$?
printf '%s\n' "$status" > results/MATRICES_EXIT_STATUS
exit "$status"
