#!/bin/bash
root=${XML_POC_ROOT}/optimized
cd "$root"
bash staging/pipeline.sh
code=$?
echo "$code" > results/PIPELINE_EXIT_STATUS
exit "$code"
