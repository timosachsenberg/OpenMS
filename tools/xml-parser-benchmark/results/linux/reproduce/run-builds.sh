#!/bin/bash
set -euo pipefail
root=${XML_POC_ROOT}
cd "$root"
export LD_LIBRARY_PATH="$root/dependencies/install/lib:/scratch/kohlbach/contrib/build/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OPENMS_TMPDIR="$root/tmp"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
prefix="$root/dependencies/install"
contrib=/scratch/kohlbach/contrib/build
for marker in BUILD_COMPLETE EXTRAS_COMPLETE; do
  if ! test -f "$root/dependencies/$marker"; then echo "Dependency marker missing: $marker"; exit 1; fi
done
mkdir -p source
if ! test -f source/CMakeLists.txt; then tar xzf staging/openms-base.tar.gz -C source; fi
cd source
if ! test -d .git; then
 git init -q
 git -c user.name='OpenMS XML POC' -c user.email='poc@localhost' add -A
 git -c user.name='OpenMS XML POC' -c user.email='poc@localhost' commit -qm 'Import exact OpenMS77cccda7 archive'
 git tag imported-base
fi
existing=(OnDiscMSExperiment ConsensusXMLFile FeatureXMLFile IndexedMzMLDecoder IndexedMzMLFile IndexedMzMLFileLoader MzMLFile MzMLSpectrumDecoder TransformationXMLFile XMLFile XMLHandler XMLValidator CVMappingFile IdXMLFile MzDataFile MzIdentMLFile MzDataValidator MzMLValidator MzXMLFile ParamXMLFile PepXMLFile ProtXMLFile SemanticValidator TraMLFile UnimodXMLFile)
targets=()
regex='^('
for name in "${existing[@]}"; do targets+=("${name}_test"); regex+="${name}|"; done
regex="${regex%|})_test$"
for backend in xerces libxml2 expat pugixml; do
 echo "BEGIN $backend $(date -Is)"
 git checkout -q -B "$backend-poc" imported-base
 git apply "$root/staging/$backend.patch"
 git -c user.name='OpenMS XML POC' -c user.email='poc@localhost' add -A
 git -c user.name='OpenMS XML POC' -c user.email='poc@localhost' commit -qm "Import original $backend feature snapshot"
 cmake -S "$root/source" -B "$root/build" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-I$root/dependencies/header-priority" -DCMAKE_CXX_FLAGS_RELEASE='-O3 -DNDEBUG' -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DCMAKE_PREFIX_PATH="$prefix;$contrib" -DOPENMS_USE_VCPKG=OFF \
  -DXercesC_INCLUDE_DIR="$prefix/include" -DXercesC_LIBRARY_RELEASE="$prefix/lib/libxerces-c.so" \
  -DXercesC_LIBRARY_DEBUG="$prefix/lib/libxerces-c.so" \
  -DLIBXML2_INCLUDE_DIR="$prefix/include/libxml2" -DLIBXML2_LIBRARY="$prefix/lib/libxml2.so" \
  -DEXPAT_INCLUDE_DIR="$prefix/include" -DEXPAT_LIBRARY="$prefix/lib/libexpat.so" \
  -UCURL_INCLUDE_DIR -UCURL_LIBRARY_RELEASE -DCURL_NO_CURL_CMAKE=OFF -DCURL_DIR="$contrib/lib/cmake/CURL" \
  -DBUILD_TOPP_TOOLS=OFF -DWITH_GUI=OFF -DPYOPENMS=OFF -DWITH_OPENTIMS=OFF -DWITH_THERMO_RAW=OFF \
  -DENABLE_DOCS=OFF -DENABLE_CLASS_TESTING=ON -DENABLE_TOPP_TESTING=OFF -DENABLE_UPDATE_CHECK=OFF \
  > "$root/logs/$backend-configure.log" 2>&1
 focused=()
 case "$backend" in libxml2) focused=(XMLParserBackend_test);; expat) focused=(ExpatXMLFile_test);; pugixml) focused=(XMLBackendPrototype_test);; esac
 cmake --build "$root/build" --target "${targets[@]}" "${focused[@]}" -j24 > "$root/logs/$backend-build.log" 2>&1
 ctest --test-dir "$root/build" -R "$regex" --output-on-failure -j1 > "$root/results/$backend-tests.log" 2>&1
 if ((${#focused[@]})); then
  ctest --test-dir "$root/build" -R "^${focused[0]}$" --output-on-failure -j1 >> "$root/results/$backend-tests.log" 2>&1
 fi
 mkdir -p "$root/snapshots/$backend/lib"
 cp -a "$root/build/lib/"libOpenMS.so* "$root/build/lib/"libOpenSwathAlgo.so* "$root/snapshots/$backend/lib/"
 cp "$root/build/CMakeCache.txt" "$root/results/$backend-build-cache.txt"
 cp "$root/build/compile_commands.json" "$root/results/$backend-compile_commands.json"
 sha256sum "$root/snapshots/$backend/lib/"* > "$root/results/$backend-libraries.sha256"
 cmake -S "$root/source/tools/xml-parser-benchmark" -B "$root/benchmark-build" -DOpenMS_DIR="$root/build" -DCMAKE_BUILD_TYPE=Release > "$root/logs/$backend-benchmark-configure.log" 2>&1
 cmake --build "$root/benchmark-build" -j8 > "$root/logs/$backend-benchmark-build.log" 2>&1
 ctest --test-dir "$root/benchmark-build" --output-on-failure > "$root/results/$backend-benchmark-tests.log" 2>&1
 LD_LIBRARY_PATH="$root/snapshots/$backend/lib:$LD_LIBRARY_PATH" "$root/benchmark-build/xml_benchmark_launcher" "$root/benchmark-build/xml_parser_benchmark" verify "$root/data/huge-array.mzML" 1 2048 > "$root/results/$backend-huge-array-equivalence.json" 2> "$root/logs/$backend-huge-array.log"
 if [[ "$backend" == libxml2 ]]; then
  cp -a "$root/source/tools/xml-poc" "$root/schema-tools"
  g++ -std=c++17 -I"$prefix/include" "$root/source/tools/xml-poc/xerces_schema_probe.cpp" -L"$prefix/lib" -Wl,-rpath,"$prefix/lib" -lxerces-c -o "$root/results/xerces_schema_probe"
  python3 "$root/source/tools/xml-poc/check_schemas.py" --library "$prefix/lib/libxml2.so" --xerces-probe "$root/results/xerces_schema_probe" --allow-known-differences > "$root/results/schema-results.json" 2> "$root/results/schema-diagnostics.txt"
 fi
 echo "DONE $backend $(date -Is)"
 touch "$root/results/$backend-BUILD_TEST_COMPLETE"
done
cp "$root/source/tools/xml-parser-benchmark/run_benchmarks.py" "$root/run_benchmarks.py"
cp "$root/staging/source-manifest.json" "$root/results/source-manifest.json"
touch "$root/results/ALL_BUILDS_COMPLETE"
