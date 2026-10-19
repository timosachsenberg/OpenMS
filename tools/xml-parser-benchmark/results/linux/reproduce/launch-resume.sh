#!/bin/bash
root=${XML_POC_ROOT}
cd "$root"
taskset -c 0-95 bash staging/resume-matrices.sh
status=$?
echo "$status" > results/MATRICES_EXIT_STATUS
exit "$status"
