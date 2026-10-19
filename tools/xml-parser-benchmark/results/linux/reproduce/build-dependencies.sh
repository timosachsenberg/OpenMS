#!/bin/bash
set -euo pipefail
root=${XML_POC_ROOT}
cd "$root/dependencies"
prefix="$root/dependencies/install"
mkdir -p "$prefix"
fetch() { test -f "$1" || curl --fail --location --retry 3 --connect-timeout 20 "$2" -o "$1"; }
fetch xerces-c-3.3.0.tar.gz https://archive.apache.org/dist/xerces/c/3/sources/xerces-c-3.3.0.tar.gz
fetch libxml2-2.15.4.tar.xz https://download.gnome.org/sources/libxml2/2.15/libxml2-2.15.4.tar.xz
fetch expat-2.9.0.tar.xz https://github.com/libexpat/libexpat/releases/download/R_2_9_0/expat-2.9.0.tar.xz
fetch pugixml-1.16.tar.gz https://github.com/zeux/pugixml/releases/download/v1.16/pugixml-1.16.tar.gz
for archive in xerces-c-3.3.0.tar.gz libxml2-2.15.4.tar.xz expat-2.9.0.tar.xz pugixml-1.16.tar.gz; do tar xf "$archive"; done
cmake -S xerces-c-3.3.0 -B xerces-build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$prefix" -DBUILD_SHARED_LIBS=ON -Dnetwork-accessor=socket -Dtranscoder=gnuiconv -DBUILD_TESTING=OFF > "$root/logs/xerces-dependency-configure.log" 2>&1
cmake --build xerces-build -j24 > "$root/logs/xerces-dependency-build.log" 2>&1
cmake --install xerces-build > "$root/logs/xerces-dependency-install.log" 2>&1
cmake -S libxml2-2.15.4 -B libxml2-build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$prefix" -DBUILD_SHARED_LIBS=ON -DLIBXML2_WITH_PYTHON=OFF -DLIBXML2_WITH_PROGRAMS=OFF -DLIBXML2_WITH_TESTS=OFF > "$root/logs/libxml2-dependency-configure.log" 2>&1
cmake --build libxml2-build -j24 > "$root/logs/libxml2-dependency-build.log" 2>&1
cmake --install libxml2-build > "$root/logs/libxml2-dependency-install.log" 2>&1
cmake -S expat-2.9.0 -B expat-build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$prefix" -DEXPAT_SHARED_LIBS=ON -DEXPAT_BUILD_DOCS=OFF -DEXPAT_BUILD_EXAMPLES=OFF -DEXPAT_BUILD_TESTS=OFF -DEXPAT_BUILD_TOOLS=OFF > "$root/logs/expat-dependency-configure.log" 2>&1
cmake --build expat-build -j24 > "$root/logs/expat-dependency-build.log" 2>&1
cmake --install expat-build > "$root/logs/expat-dependency-install.log" 2>&1
cmake -S pugixml-1.16 -B pugixml-build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$prefix" -DBUILD_SHARED_LIBS=ON > "$root/logs/pugixml-dependency-configure.log" 2>&1
cmake --build pugixml-build -j24 > "$root/logs/pugixml-dependency-build.log" 2>&1
cmake --install pugixml-build > "$root/logs/pugixml-dependency-install.log" 2>&1
sha256sum *.tar.* > "$root/results/dependency-archives.sha256"
touch "$root/dependencies/BUILD_COMPLETE"
