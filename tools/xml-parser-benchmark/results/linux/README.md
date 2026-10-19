# Linux reproduction and provenance

Host: `ibminode05`, Ubuntu 24.04, AMD EPYC 7763, 64 physical/128 logical CPUs, one NUMA node, approximately 995 GiB RAM. Effective locale is `C.UTF-8`. Build/compiler: GCC 13.3.0, CMake 3.28.3, C++23, Release `-O3 -DNDEBUG`.

The baseline archive is the exact OpenMS 3.7 develop commit `77cccda7caf40b46a4919d33748450270eccd565`. Its 9022 non-export-subst tracked blobs were checked against the original Git object IDs; `.git_archival.txt` intentionally expands archive metadata. Feature patches and original revision IDs are in `reproduce/source-manifest.json`. Imported local Git branches preserve archived file bytes using `.git/info/attributes` with `* -text`.

Private scratch parser builds use Xerces 3.3.0, libxml2 2.15.4, Expat 2.9.0, pugixml 1.16 and zstd 1.5.7. The original tarball checksums are retained. Existing read-only dependencies from `/scratch/kohlbach/contrib/build` supply Arrow 23.0.0, Boost 1.87, Eigen 5.0.1, COIN-OR, libSVM, curl 8.12.1 and compression dependencies. No installed/global dependency was modified.

Two configuration issues were resolved before measurements:

* Arrow's existing static package requires its matching `CURL::libcurl_static` imported target. The build uses the original curl 8.12.1 CMake package instead of substituting the system curl target.
* That prefix's broad include directory contains old Xerces 3.2 headers. An isolated `dependencies/header-priority` include directory symlinks only the private parser headers and precedes the broad dependency include. A compile-time probe verifies Xerces 3.3 headers; all OpenMS objects were rebuilt afterward. This was a dependency-prefix collision, rather than a parser-source portability failure.

`reproduce/build-dependencies.sh`, `build-extras.sh`, `prepare-source-and-headers.sh`, `run-builds.sh`, `test-expat.sh`, `verify-source.py`, and `verify-inputs.py` contain the exact setup and checks. The build script runs 25 affected class tests for each backend, one additional focused test for each prototype, the shared benchmark smoke checks, and huge-array equivalence checks. libxml2 additionally runs 22 schema compilations and 195 instance comparisons, allowing only the two recorded anyURI differences. Expat's upstream `runtests` suite also runs.

All six timing inputs plus the huge-array regression input match the macOS bytes and SHA256. Input verification warms the file cache. Each backend has a separate copied library snapshot; the measurement runner checks the actual loaded library path and compares decoded hashes/counts across backends, repetitions, and worker counts.

To reproduce on the same node/dependency prefix, choose a fresh private scratch directory and export `XML_POC_ROOT` to that absolute path. The published assembler substitutes that variable for the original job root in these exact executed helpers. Create `dependencies`, `staging`, `results`, `logs`, `tmp`, and `data`, copy `reproduce/` into `staging/`, then run `build-dependencies.sh`, `build-extras.sh`, `prepare-source-and-headers.sh`, `run-builds.sh`, and `test-expat.sh` in that order. The helpers retain the original read-only contrib prefix; using another installation requires corresponding CMake dependency paths. Populate `data/` with the five generated timing files and `huge-array.mzML` from the main reproduction recipe, keeping the recorded hashes; BSA comes from the exact source archive. The build script creates frozen backend snapshots and the primary benchmark executable. The metadata source directory compiles separately after the main matrices.

The intended measurement launch is:

```sh
taskset -c 0-95 bash staging/run-matrices.sh
```

CPUs 0–95 include 64 physical cores and 32 SMT siblings; 96 workers must not be described as 96 physical cores. Environment snapshots include actual affinity and topology. There are no concurrent builds during timing. The script runs 2560 primary measurements, 960 whole-load/context/idle measurements, and 320 fixed-work OpenMP measurements, using five randomized repetitions and worker/thread counts 1, 2, 4, 8, 16, 32, 64, 96. Each concurrent reader repeats the same workload; OpenMP measurements instead keep one reader and one fixed input workload. A separately compiled metadata supplement then adds 640 measurements of default OnDisc opens, using `indexed-metadata` and `contexts-metadata`, for 4480 Linux measurements overall. It preserves the original running executable.

A read-only observation of live fixed-work runs found 192 native process threads at `OMP_NUM_THREADS=96` (and 128 at 64). No inherited `OMP_THREAD_LIMIT` or GOMP/KMP cap was present; the linked runtime is system libgomp. Warm-up on the main thread can retain one OpenMP pool while the reader thread creates another, so native process-thread count and RSS include parked warm-up threads as well as the requested per-reader team. The raw observation and its helper are retained.

Linux `wait4.ru_maxrss` is converted from KiB to bytes. `resident_rss_bytes` samples `/proc/self/statm` after all workers join while whole experiments or OnDisc contexts remain retained. Linux RSS counters can lag or account differently with many threads, so the sampled resident value can slightly exceed the separately reported peak. Neither value is treated as byte-exact native parser allocation. Peak and resident values are distinct; baseline-adjusted context estimates include indexes, allocator and thread costs rather than only a native parser structure. Raw JSON, build caches, compiler commands, test logs, dependency checksums and resource monitoring remain the evidence.

The first OpenMP launch stopped before any measurement because the staged Python helper predated its `--threads` option. The completed 2560 primary and 960 memory measurements were preserved. `resume-matrices.sh` uses the current shared helper to run all 320 requested OpenMP measurements, followed by the unchanged metadata supplement. The corrected helper SHA256 and original argument-error log are retained; the benchmark binary and OpenMS snapshots were unchanged.

Scratch was retained for the user-requested external reviews, optimized rebuilds and before/after reruns. Cleanup runs only after every required measurement succeeds and the local evidence passes integrity checks. The actual cleanup outcome and preserved shared contrib prefix are recorded in [scratch-cleanup.json](scratch-cleanup.json).
