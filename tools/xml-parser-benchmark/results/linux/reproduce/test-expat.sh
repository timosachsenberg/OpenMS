#!/bin/bash
set -euo pipefail
root=${XML_POC_ROOT}
cmake -S "$root/dependencies/expat-2.9.0" -B "$root/dependencies/expat-build" -DEXPAT_BUILD_TESTS=ON
cmake --build "$root/dependencies/expat-build" -j12
ctest --test-dir "$root/dependencies/expat-build" --output-on-failure
