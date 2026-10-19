#!/bin/bash
set -euo pipefail
root=${XML_POC_ROOT}
cd "$root/dependencies"
curl --fail --location --retry 3 --connect-timeout 20 https://github.com/facebook/zstd/releases/download/v1.5.7/zstd-1.5.7.tar.gz -o zstd-1.5.7.tar.gz
tar xf zstd-1.5.7.tar.gz
cmake -S zstd-1.5.7/build/cmake -B zstd-build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$root/dependencies/install" -DZSTD_BUILD_PROGRAMS=OFF -DZSTD_BUILD_TESTS=OFF -DZSTD_BUILD_STATIC=OFF -DZSTD_BUILD_SHARED=ON > "$root/logs/zstd-configure.log" 2>&1
cmake --build zstd-build -j12 > "$root/logs/zstd-build.log" 2>&1
cmake --install zstd-build > "$root/logs/zstd-install.log" 2>&1
sha256sum zstd-1.5.7.tar.gz > "$root/results/zstd-archive.sha256"
touch "$root/dependencies/EXTRAS_COMPLETE"
