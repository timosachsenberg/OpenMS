# OpenMS XML parser benchmarks

Compare the actual OpenMS mzML paths on `feature/xml-libxml2-poc`, `feature/xml-expat-poc`, and `feature/xml-pugixml-poc` against upstream develop `77cccda7caf40b46a4919d33748450270eccd565` (3.7). See [HANDOFF.md](HANDOFF.md) for the measured results and remaining migration work.

## Build and correctness

Use the same compiler, optimization, dependencies, and OpenMP configuration for every branch. The macOS run used Apple Clang 21, Release `-O3 -DNDEBUG`, C++23, Xerces 3.3.0, libxml2 2.15.4, Expat 2.9.0, and pugixml 1.16. Expat was built from upstream tag `R_2_9_0` into a private prefix; its upstream test suite passed. Linux used GCC 13.3 on Ubuntu 24.04 with the same optimization, C++ standard and parser versions. No global dependencies were modified.

The following configuration disables unrelated optional dependencies. Homebrew-specific paths should be replaced on other systems. The explicit CURL paths avoid an obsolete system framework on the measured host.

```sh
cmake -S . -B ../OpenMS-build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DOPENMS_USE_VCPKG=OFF \
  -DWITH_GUI=OFF -DBUILD_TOPP_TOOLS=OFF -DENABLE_DOCS=OFF \
  -DWITH_OPENTIMS=OFF -DWITH_THERMO_RAW=OFF -DWITH_ONNX=OFF \
  -DWITH_WNETALIGN=OFF -DGIT_TRACKING=OFF -DENABLE_UPDATE_CHECK=OFF \
  -DMT_REQUIRE_OPENMP=ON \
  '-DCMAKE_PREFIX_PATH=/opt/homebrew;/opt/homebrew/opt/libxml2;/opt/homebrew/opt/libomp' \
  -DCURL_INCLUDE_DIR=/opt/homebrew/opt/curl/include \
  -DCURL_LIBRARY=/opt/homebrew/opt/curl/lib/libcurl.dylib \
  -DCURL_LIBRARY_RELEASE=/opt/homebrew/opt/curl/lib/libcurl.dylib
```

For libxml2, set `LIBXML2_INCLUDE_DIR` and `LIBXML2_LIBRARY` explicitly if the system provides an older version; the prototype requires >=2.13. For Expat, set `EXPAT_INCLUDE_DIR` and `EXPAT_LIBRARY` to the private 2.9.0 prefix. pugixml uses its CMake config package. The vcpkg manifests declare the additional Expat/pugixml dependencies on their respective branches.

| Branch | Focused class test |
|---|---|
| libxml2 | `XMLParserBackend_test` |
| Expat | `ExpatXMLFile_test` |
| pugixml | `XMLBackendPrototype_test` |

Build the selected test target and the existing affected class tests, then run CTest. The exact 25 existing test names and the per-branch focused test outcomes are recorded in [results](results/). The libxml2 [schema comparator](https://github.com/okohlbacher/OpenMS/tree/feature/xml-libxml2-poc/tools/xml-poc) separately compiles all bundled schemas and compares 195 instance/schema pairs; its strict mode deliberately fails on the recorded URI difference.

Build the external benchmark against the same OpenMS build tree:

```sh
cmake -S tools/xml-parser-benchmark -B ../benchmark-build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DOpenMS_DIR="$PWD/../OpenMS-build" \
  '-DOpenMP_CXX_FLAGS=-Xclang -fopenmp' -DOpenMP_CXX_LIB_NAMES=omp \
  -DOpenMP_omp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib \
  -DOpenMP_CXX_INCLUDE_DIR=/opt/homebrew/opt/libomp/include
cmake --build ../benchmark-build
ctest --test-dir ../benchmark-build --output-on-failure
```

The benchmark smoke test generates plain/zlib mzML and compares full loading, event delivery, and indexed access, including coordinates, intensities, native IDs, and selected metadata. Auxiliary arrays are compared when present, but the synthetic generator does not populate them; existing OpenMS tests supply broader format coverage. It checks independent workers, failure propagation, and the stale OpenMP registration regression.

## Generate the shared inputs

The commands in this section and the run command below assume the workspace directory containing `OpenMS/`, `benchmark-build/`, and the snapshot directories. Use `OpenMS/tools/xml-parser-benchmark/run_benchmarks.py` and prefix the tracked example path with `OpenMS/` from that directory.

Generate once with the Xerces baseline and give every backend the identical files. Arrays use a fixed random seed; peak coordinates and intensities vary across spectra. The writer streams output instead of retaining the full experiment. Save the generated file hashes with every run.

```sh
mkdir -p benchmark-data
benchmark-build/xml_benchmark_launcher benchmark-build/xml_parser_benchmark generate benchmark-data/small.mzML 256 1024 plain
benchmark-build/xml_benchmark_launcher benchmark-build/xml_parser_benchmark generate benchmark-data/medium.mzML 4096 2048 plain
benchmark-build/xml_benchmark_launcher benchmark-build/xml_parser_benchmark generate benchmark-data/large.mzML 32768 2048 plain
benchmark-build/xml_benchmark_launcher benchmark-build/xml_parser_benchmark generate benchmark-data/medium-zlib.mzML 4096 2048 zlib
benchmark-build/xml_benchmark_launcher benchmark-build/xml_parser_benchmark generate benchmark-data/huge-array.mzML 1 1048576 plain
```

The run also uses gzip around `medium-zlib.mzML` and the tracked real-data file `share/OpenMS/examples/FRACTIONS/BSA1_F1.mzML`. The huge-array file is a correctness regression, not part of the timed random-query matrix.

## Run the comparison

Freeze each built OpenMS library and its OpenSwath companion in a separate directory. Preserve all runtime dependencies; the measured Expat snapshot also contains its exact `libexpat.1.dylib`. On macOS the runner selects snapshots with `DYLD_LIBRARY_PATH`; on Linux it uses `LD_LIBRARY_PATH`. Every child reports its actual loaded OpenMS path, which must match the selected snapshot. Binary and library SHA256 hashes are recorded.

```sh
python3 OpenMS/tools/xml-parser-benchmark/run_benchmarks.py \
  --binary benchmark-build/xml_parser_benchmark \
  --backend xerces=/absolute/snapshots/xerces/lib \
  --backend libxml2=/absolute/snapshots/libxml2/lib \
  --backend expat=/absolute/snapshots/expat/lib \
  --backend pugixml=/absolute/snapshots/pugixml/lib \
  --input benchmark-data/small.mzML --input benchmark-data/medium.mzML \
  --input benchmark-data/large.mzML --input benchmark-data/medium-zlib.mzML \
  --input benchmark-data/medium-zlib.mzML.gz \
  --input OpenMS/share/OpenMS/examples/FRACTIONS/BSA1_F1.mzML \
  --workers 1,2,4,8 --repetitions 5 --queries 2048 --output results/benchmarks.json
```

## Measurement meaning

- `stream`: `MzMLFile::transform` with a single-pass, non-retaining consumer. The pugixml implementation builds a whole document before delivering events; it is not a streaming parser.
- `indexed`: `OnDiscMSExperiment::openFile(..., true)` followed by 2048 deterministic random spectrum requests per worker. Report opening and query phases separately. All three prototypes replace the footer and fragment decoder; Expat does not retain Xerces for the measured indexed path.
- `index-open`: index parsing/counting without retrieving peaks. Index counts must match across workers and backends.
- `load`: full `MzMLFile::load`, retaining each complete `MSExperiment` until all workers finish. This measures whole-dataset loading; libxml2/Expat still parse XML incrementally, while pugixml additionally builds the XML DOM.
- `indexed-metadata`: default `OnDiscMSExperiment::openFile` including the full metadata scan, followed by the same random spectrum requests. `contexts-metadata` retains each default reader and its loaded metadata until all workers finish. These paths are measured in `metadata-scaling.json`, separately from the fast index-only path.
- `contexts`: open one `OnDiscMSExperiment` per worker with global metadata skipped, retaining every reader and its index until all workers finish. `idle` measures common initialization and worker creation without reading the input.
- `max_rss_bytes` is peak whole-process RSS from `wait4`. Newer binaries also report `resident_rss_bytes`, sampled after joining while retained reader contexts/experiments remain alive (`/proc/self/statm` on Linux; Mach task info on macOS). The baseline-adjusted resident estimate is `(RSS_contexts(N) - RSS_idle(N)) / N`. It includes index storage and allocator retention; it is not the byte size of a native parser context. Linux `/proc/self/statm` uses approximate RSS accounting, so sampled resident values can occasionally exceed `wait4`’s reported peak; keep the counters separate ([Linux man-pages](https://man7.org/linux/man-pages/man5/proc_pid_statm.5.html)).
- Each worker reads the same cached file and, for indexed access, the same query sequence. This is independent-parser throughput/weak scaling, not fixed-work speedup or independent-disk I/O scaling. `OMP_NUM_THREADS=1` avoids nested OpenMP work.
- Timed work includes actual binary decoding plus a consumer hashing all peak coordinates, intensities, and native IDs. This is an OpenMS workload, not a tokenizer-only benchmark. Cross-backend digests must agree.
- A common tiny in-memory mzML initializes shared CV data before timing. It never reads the timed input. RSS is whole-process peak RSS from `wait4`, including initialization. Tables report the median of five per-process peaks.
- Workload order is randomized with a fixed seed. Inputs are hashed before timing, warming the file cache. Raw records retain repetitions, min/max/median, per-worker timings, CPU time, and RSS.
- The POSIX launcher handles a measured-host libomp startup defect: only its own old, empty, regular `/tmp/__KMP_REGISTERED_LIB_<pid>` file owned by the current UID may be removed before `exec`. No nonempty files or other PIDs are touched. `exec` preserves the PID and avoids Python-fork RSS contamination. Startup failures are retained separately from parser failures; see [the diagnosis](results/macos/runtime-startup-failure.json).

These branches retain Xerces for mzIdentML DOM and shared input/transcoding infrastructure; Expat/pugixml also retain Xerces XSD validation. `XMLAttributes` gains storage and changes layout, so SDK consumers must be rebuilt. The recorded builds and tests cover macOS ARM64 and Linux x86-64. The benchmark runner/launcher are POSIX-only and need a Windows port before measuring that platform. Windows runtime behavior remains untested; concurrent-reader results do not establish shared-object safety or race freedom.

## Memory and single-file parallelism

Run memory probes separately from the primary matrix with `--modes load,contexts,idle`, the same backends, and inputs `large.mzML` and `medium-zlib.mzML`. On Linux use `--workers 1,2,4,8,16,32,64,96`; macOS measurements stop at eight independent readers. The macOS residency-only supplement uses `--modes contexts,idle` and a newer small benchmark binary so the original timing binary remains frozen.

For fixed-work OpenMP scaling, keep one reader and vary binary-decoding threads:

```sh
python3 OpenMS/tools/xml-parser-benchmark/run_openmp_scaling.py \
  --matrix results/benchmarks.json --output results/openmp-scaling.json
```

On Linux add `--threads 1,2,4,8,16,32,64,96`. This reuses the completed primary matrix's exact binary, inputs, snapshots and expected decoded digests. It changes `OMP_NUM_THREADS`, leaves `OPENBLAS_NUM_THREADS=1`, and explicitly disables dynamic OpenMP thread counts. Each combination has five randomized repetitions. The x-axis gives the requested per-reader OpenMP team size, not the total native process-thread count. Tiny-input warmup runs on the main thread and the timed reader on a separate C++ thread. Linux observations found 128 process threads at OpenMP=64 and 192 at OpenMP=96, consistent with a retained initialization team plus the reader team; no inherited OpenMP thread limit was found. RSS includes this runtime overhead; see [runtime observations](results/linux/openmp-runtime-observation.json). Independent-reader matrices use OpenMP=1. Syntax parsing and the digest consumer remain serial; OpenMS's array decoding supplies OpenMP work.

The Linux run uses CPUs 0–95 on the selected host: 64 physical cores plus 32 SMT siblings, within one NUMA node. Every tested thread count uses that same allowed CPU set. The report marks the physical-core boundary; these data cannot establish scaling on 96 physical cores. Test/build work and timings run separately. Inputs are node-local in `/scratch`; temporary files use a private directory inside that job’s `/scratch` folder. The tested gzip path decompresses incrementally and does not unpack a ZIP archive.

`plot_results.py --mac MAC_RESULT_DIRECTORY --linux LINUX_RESULT_DIRECTORY --output FIGURE_DIRECTORY` regenerates PNG/SVG figures (concurrent throughput, peak memory, default/fast indexed paths, resident context overhead, and fixed-work OpenMP speed/memory) using matplotlib. Both directories contain `benchmarks.json`, `memory-scaling.json`, `openmp-scaling.json`, and `metadata-scaling.json`; a separate `context-residency.json` overrides the memory report for resident-context estimates when present. Graphs use within-host comparisons and median-of-five results; different CPUs/toolchains preclude attributing differences solely to the operating system.


## Recorded evidence and validation probes

The [handoff](HANDOFF.md) integrates 13,760 final measurements: 7,040 original and 6,720 reviewed/optimized. The original phase comprises the following evidence. [Mac results](results/macos/) and [Linux results](results/linux/) include raw records, input/library/executable SHA256, summaries and test logs. [comparison.csv](results/comparison.csv) contains all modes and thread counts across both hosts. Regenerate the six PNG/SVG figures with:

```sh
python3 plot_results.py --mac results/macos --linux results/linux --output figures
```

Run that command from this tool directory in a Python environment with matplotlib installed. Plotting versions used for the recorded figures are in `plotting-requirements.txt`.

Each implementation stays fixed within its original or reviewed measurement phase. In the original phase, the macOS primary/OpenMP executable uses the harness at libxml2-branch ancestor `715d804`; retained whole loads use `17a8a3ef`; resident/default-metadata probes use `28a160db`. Linux uses the harness at `cf45d453` for primary/memory/OpenMP and a separate metadata executable from `28a160db`. These harness updates add measurement modes or retention/resident sampling; they do not change the OpenMS parser implementation. Raw reports pin the actual executable, launcher and library hashes. The Linux [source manifest](results/linux/reproduce/source-manifest.json) and patches retain the exact initial build inputs; supplemental scripts retain later harness setup. The public branch tips add documentation, plotting and evidence to those tested ancestors.

The `xml_schema_validation_probe` target exercises the actual OpenMS `XMLValidator::isValid` API:

```sh
../benchmark-build/xml_schema_validation_probe INPUT_XML BUNDLED_XSD
```

Select the matching library snapshot through the same runtime library path used for timing. The probe accepts exactly those two positional arguments and returns zero for validity. The recorded runs validate the huge-array input and the two minimal URI cases on all four implementations. libxml2's standalone schema comparator additionally tests all bundled schemas and the larger fixture corpus; Expat/pugixml still route validation through Xerces.

Private Linux builds/runs use the `ibmi-hpc` workflow and a task-specific `/scratch` tree; exact dependency setup, topology, CPU affinity, resource monitor and cleanup evidence accompany the [Linux reproduction notes](results/linux/README.md). The source archive, build trees and dataset binaries are deliberately excluded from the public evidence package. Inputs can be regenerated with this harness; the exact SHA256 manifest identifies the recorded bytes.


## Reviewed implementations and before/after evidence

The reviewed phase adds 2,240 macOS and 4,480 Linux measurements using the same inputs, timed work, digests, worker counts and five-repetition policy. macOS no longer duplicates its resident-context matrix: the current memory runner records both peak and retained RSS. An additional 288-run preliminary screen and 48-run unchanged-Xerces calibration are separate from final medians. [Selection and calibration](results/optimized/macos/SELECTION.md) explain why small apparent changes are not assigned to a particular patch.

[The integrated before/after CSV](results/comparison-before-after.csv) includes all original/reviewed cells, retained reader overhead and idle-adjusted excess peak RSS per reader where matching idle cases exist. These estimates include reader indexes, metadata, allocator retention, decoding buffers and retained data; they do not isolate native parser structs.

Source/library fingerprints and test logs are in [reviewed macOS results](results/optimized/macos/) and [reviewed Linux results](results/optimized/linux/). The optimized Xerces control applies shared OpenMS fixes and retains native Xerces parsing; it is not a fourth migration branch. Its complete patch against the same develop base is under `results/optimized/macos/reproduce/xerces-control.patch`. The manifest records its source tree ID so a fresh base checkout plus that patch can verify identical source with `git write-tree`, independently of reconstructing local commit metadata.

The additive `XMLHandler` virtual hook changes ABI; rebuild the benchmark and SDK consumers against the matching headers. Keep original executables with original snapshots and reviewed executables with reviewed snapshots. Mixing old/new vtables can invalidate the experiment.

The source-only review records contain actual Vibe/Kimi/Claude final responses, their fixed input, CLI settings and [evidence-based triage](reviews/REVIEW-DECISIONS.md). The Vibe run uses a task-local thinking-off profile after high-thinking timeouts; no private reasoning output is published. Reviewer suggestions are not benchmark evidence.

`xml_conformance_probe` exercises generic XML event delivery on small malformed/encoding cases. The corpus/runner is supplied in the reviewed reproduction files. Acceptance is recorded separately from decoded UTF-16 values: accepting malformed XML is not a test pass. The existing `xml_schema_validation_probe` separately repeats the huge-array and URI compatibility cases through the actual validator API.

Regenerate reviewed and combined figures from this tool directory:

```sh
python3 plot_results.py --mac results/optimized/macos --linux results/optimized/linux --output figures/optimized
python3 plot_optimization.py --results results --output figures
```

The archived sources/scripts may require adapting workspace/dependency paths. `${XML_POC_ROOT}` inside a Python string is a publication placeholder, not automatic environment expansion: replace it with the chosen absolute scratch path in a copied reproduction directory before executing. Git branch sources and parser-dependency archive hashes remain pinned.

Transfer checksum manifests describe the temporary remote workspace before publication path sanitization; they may name build artifacts excluded from this source-only package. `EVIDENCE-SHA256.txt` instead hashes the final published files. Historical validation records may say cleanup was pending at transfer time; the later [cleanup record](results/linux/scratch-cleanup.json) records the completed outcome.

The [cross-clock audit](results/optimized/macos/clock-consistency.json) preserves eight discrepant macOS elapsed-time observations (seven original OpenMP, one optimized metadata) and records separate confirmation runs; no Linux case has the discrepancy. C++ uses `steady_clock` and Python uses `perf_counter`; the reason for their divergence is not established. Original rows and final five-run medians remain unchanged. These clock anomalies must not be interpreted as parser CPU stalls.
