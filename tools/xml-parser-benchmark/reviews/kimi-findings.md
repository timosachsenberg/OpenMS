• # Adversarial Review: libxml2 / Expat / pugixml XML-parser prototypes

  Scope note: findings are based solely on the embedded source and the reported measurements. Items marked **Hypothesis** are plausible but not provable from the supplied material; each carries a falsification test.

  ---

  ## A. High-confidence findings

  ### A1. Per-lookup heap-allocating transcoder in `XMLAttributes::value(const char*)` — all three branches
  **Files/functions:** `libxml2/src/.../XMLAttributes.cpp`, `expat/src/.../XMLAttributes.cpp`, `pugixml/src/.../XMLAttributes.cpp`, `value(const char* qname)`.
  **Mechanism:** Every attribute lookup by narrow name constructs a `xercesc::TranscodeFromStr` (pugi branch: unconditionally; libxml2/expat: on the native path). This is the hottest accessor in every handler (`attributeAsString_`/`attributeAsInt_`/etc. in `XMLHandler.h` call it once per attribute per element; a cvParam pays 4+ lookups). Each call allocates a transcoder object plus a transcoded `XMLCh` buffer, per lookup, per element, for a 1.05 GiB file with millions of elements.
  **Minimal change:** Pre-transcode the (small, fixed) set of attribute names once — e.g. `constexpr char16_t` literals as already done in `MzMLHandler.cpp` lines 810–850 (`s_count`, `s_default_array_length`) — or cache per-`XMLAttributes` transcoded names. Do not transcode on every query.
  **Risk:** None functional; identical comparison semantics.
  **Falsify:** Allocation profiler (e.g. `malloc` flamegraph / `heaptrack`) on the single-reader stream workload; if `TranscodeFromStr` construction is <1% of samples, this finding is wrong. Compare stream seconds before/after on the Mac 1-reader case (lib2 5.536 / E 6.776).

  ### A2. Per-SAX-event allocation storm in the libxml2 and Expat adapters
  **Files/functions:** `libxml2/.../XMLFile.cpp` `LibxmlState::start/characters`; `expat/.../XMLFile.cpp` `ExpatCallbacks::start/end/characters`, `utf16()`.
  **Mechanism:** libxml2: one `std::u16string` per element name, two per attribute, plus a `std::vector<std::pair<...>>` per start-tag; no `reserve` in `utf16()`. Expat: a `std::vector<std::u16string>` (one heap string per attribute), a second pointer vector, and `utf16(name)` computed twice (start and end). `characters` chunks (base64 text dominates byte volume) each allocate a fresh `u16string` before `appendASCII` copies it again in `MzMLHandler::onCharacters`.
  **Minimal change:** Add a reusable scratch arena to the parser-state struct (already per-parse, so thread-safe): one grow-only `std::u16string` buffer for names/attrs and one `std::vector` for the attribute pairs, cleared per event. Reserve `utf16()` result to `length` (UTF-16 length ≤ UTF-8 length).
  **Risk:** Scratch lifetime must not outlive the callback (handlers receive `c_str()` pointers); reuse is safe only because callbacks are synchronous — a test asserting digests across all 25 XML classes covers this.
  **Falsify:** Allocation count per stream run should drop by ~O(elements+attributes); if Mac 1-reader stream time does not improve measurably, the hypothesis that conversion/allocations explain the ~2× gap vs Xerces (2.926 vs 5.536/6.776) is weakened.

  ### A3. pugixml branch: function-local `static` enforced encoding is a one-shot global
  **File/function:** `pugixml/src/.../XMLFile.cpp`, `parse_()` and `parseBuffer_()`: `static const XMLCh* s_enc = xercesc::XMLString::transcode(enforced_encoding_.c_str());`
  **Mechanism:** The `static` is initialized exactly once; the first `XMLFile` instance's `enforced_encoding_` wins for the lifetime of the process, and every later instance with a different encoding is silently ignored. Also a deliberate leak. Both sibling branches handle this per-instance.
  **Minimal change:** Transcode per call (as libxml2/expat do) or store the transcoded value as an instance member when `enforceEncoding_` is called.
  **Risk:** None; restores per-instance semantics.
  **Falsify:** Two `XMLFile` objects with different enforced encodings (`UTF-8`, `ISO-8859-1`) used in sequence on files requiring the respective override; the second file must parse identically to the libxml2/expat branches.

  ### A4. pugixml whole-file path: unbounded recursive `visit` + acceptance of documents sibling parsers reject
  **File/function:** `pugixml/src/.../XMLFile.cpp` `parse()`, `std::function<void(pugi::xml_node)> visit`.
  **Mechanism:** (a) Traversal recursion depth equals XML nesting depth, with a large `std::function` frame per level; a well-formed document nested ~10⁵–10⁶ deep overflows the thread stack — a crash on hostile input, not a clean `ParseError`. libxml2 caps depth (256 default) and Expat reports an error via accounting; pugi accepts arbitrary depth and the user code then recurses. (b) `std::function` type erasure adds an indirect call per node.
  **Minimal change:** Convert `visit` to an explicit stack/iterative traversal (siblings pointers make this trivial), which removes both the stack-depth hazard and the `std::function` overhead.
  **Risk:** Event order must remain identical (pre-order, text nodes interleaved); digest equality across branches falsifies regressions.
  **Falsify:** Generated file with 500k nested elements: pugi branch currently crashes or at minimum diverges from libxml2/expat rejection; after the fix it must deliver all events without crash.

  ### A5. pugixml does not validate character encoding — malformed-input rejection is not preserved
  **File/function:** `pugixml/src/.../XMLFile.cpp` `parse()` (`load_buffer_inplace`) and `pugixml/.../MzMLSpectrumDecoder.cpp` `domParseString_`.
  **Mechanism:** pugixml performs structural well-formedness checks but does not validate that text is legal UTF-8 (invalid continuation bytes, overlong encodings, stray 0xFF pass through). Xerces, libxml2 and Expat all reject such input. The stated constraint "preserve malformed input rejection" is violated for encoding-level malformation, and `utf16()` then decodes garbage into the handler with no error.
  **Minimal change:** Add a UTF-8 validation pass over `input` before `load_buffer_inplace` when encoding is UTF-8/auto (a linear SIMD-able scan; must be measured, not assumed cheap), or document the divergence explicitly.
  **Risk:** Rejecting inputs previously accepted by this prototype only; must match Xerces accept/reject behavior exactly on the existing test corpus.
  **Falsify:** Corpus of UTF-8 malformations (truncated sequences, overlongs, isolated continuations) injected into element text and attribute values: assert identical accept/reject across the three prototypes and the Xerces baseline. If pugixml 1.16 actually rejects these (verifiable with a 20-line standalone), this finding is wrong.

  ### A6. `findIndexListOffset`: UB on regex mismatch and on files smaller than the tail buffer — all three branches
  **File/function:** `*/src/openms/source/FORMAT/HANDLERS/IndexedMzMLDecoder.cpp`, `findIndexListOffset`.
  **Mechanism:** (a) `boost::regex_search(buffer.get(), matches, listoffset_rx)` result is discarded; `std::string thismatch(matches[1].first, matches[1].second)` is executed unconditionally. On no match, `matches[1]`'s iterators are null → constructing a `std::string` from null pointers is undefined behavior (crashes observed in practice with libstdc++ debug, silent garbage otherwise). (b) `f.seekg(-buffersize, f.end)` on a file smaller than `buffersize` sets failbit; the subsequent `read` leaves the heap buffer uninitialized and the regex runs on garbage.
  **Minimal change:** Check the `regex_search` return value before touching `matches[1]`; clamp the seek/read length to file size and check stream state after `read`.
  **Risk:** None; only error paths change. Note the error-message branch (`Could not find element indexListOffset...`) is currently unreachable dead code precisely because of (a).
  **Falsify:** Under ASan/UBSan: (1) a 100-byte file, (2) a valid mzML whose last `buffersize` bytes lack `<indexListOffset>`. Expected: clean `-1`/exception; observed (pre-fix): UB/crash.

  ### A7. Data race in `MzMLHandler::populateSpectraWithData_` exception path
  **File/function:** common `MzMLHandler.cpp` lines ~200–340, `#pragma omp parallel for` block.
  **Mechanism:** The `catch (...)` handler executes `++errCount; error_message = "...";` with **no** `omp critical`/`atomic`, while the `BaseException` handler uses `omp critical(MZMLErrorHandling)`. Concurrent non-OpenMS exceptions from two threads race on both `errCount` (lost updates → `errCount != 0` check unreliable) and `error_message` (torn `std::string`, potential heap corruption). The chromatogram variant uses `omp atomic` for the counter in the catch-all — the spectrum variant does not.
  **Minimal change:** Wrap the catch-all body in the same `omp critical(MZMLErrorHandling)` region (or drop the message write there and only count, atomically).
  **Risk:** None; error-path only.
  **Falsify:** TSan build, input with ≥2 spectra whose binary data triggers non-OpenMS exceptions (e.g. `std::bad_alloc` from corrupt declared lengths), OMP threads ≥ 2: TSan must report the race pre-fix and none post-fix.

  ### A8. Double copy of the index tail in `parseOffsets` — all three branches
  **File/function:** `*/IndexedMzMLDecoder.cpp::parseOffsets`.
  **Mechanism:** The file tail is read into a raw `char[]`, then copied again into `tmp_fixed_xml = "<indexedmzML>" + std::string(buffer) + "\n"`, then (libxml2/pugi) a DOM is built on top. For a large index (the code itself warns about index >10% of file) this is 2× the index size in transient memory plus a full extra copy pass. libxml2 branch additionally uses raw `new(std::nothrow)[]`/`delete[]` instead of RAII (leak on exception from `domParseIndexedEnd_` — e.g. `stringToStreampos` throws `ConversionError`, `buffer` leaks; Expat branch already uses `unique_ptr`).
  **Minimal change:** Allocate one buffer sized `prefix + readl + suffix`, `read` into it at offset `prefix.size()`, and parse it in place; use RAII in the libxml2 branch.
  **Risk:** None; byte content identical.
  **Falsify:** Synthetic mzML with a 100 MB index: peak RSS during `parseOffsets` and wall time vs before; also leak-check (ASan LSan) the conversion-error path.

  ### A9. libxml2 fragment decode builds a full DOM per query; Expat streams and pugixml parses leaner
  **File/function:** `libxml2/src/.../HANDLERS/MzMLSpectrumDecoder.cpp::domParseString_`.
  **Mechanism:** Each indexed `getSpectrum()` constructs a libxml2 DOM (node structs, property strings, `xmlGetProp` allocations per cvParam attribute) only to walk it once and free it. Measured: Mac 2048-query fragment time lib2 0.205 s vs E 0.187 / P 0.144 — libxml2 is the slowest prototype here, consistent with DOM overhead.
  **Minimal change:** Parse the fragment with the libxml2 push-parser + SAX callbacks (the same machinery already in this branch's `XMLFile.cpp`), populating `BinaryData` directly — mirroring the Expat `Fragment` struct's logic, including stop-at-root-end via `xmlStopParser`. This stays within libxml2 (no parser replacement).
  **Risk:** Must preserve the single-text-node-child rule for `<binary>` (track CDATA/comment callbacks as errors), root-name/defaultArrayLength checks, and entity semantics; covered by the existing focused indexed-read tests plus digest equality.
  **Falsify:** Mac/Linux fragment-query benchmark: lib2 should move from 0.205 s toward ≤0.19 s with identical digests; if it does not, DOM cost is not the differentiator and the gap lies elsewhere (e.g. `xmlProperty` allocations — profile will show).

  ### A10. `XMLValidator` recompiles the schema grammar on every call
  **File/function:** `common/reference .../VALIDATORS/XMLValidator.cpp::isValid`.
  **Mechanism:** `xmlSchemaNewParserCtxt` + `xmlSchemaParse` run per `isValid()` invocation. For the 22-XSD suite and any repeated validation of many files against one schema, schema compilation (file I/O + grammar build) dominates per-call cost. The streaming `xmlSchemaValidateFile` itself is appropriately bounded-memory.
  **Minimal change:** Cache the compiled `xmlSchema*` keyed by canonical schema path, behind a mutex (the cache is shared state; this is a bounded cache lock, not a serialization of parsing — validation itself still runs in per-call `xmlSchemaValidCtxt`, which is already the case). Make cache opt-out to preserve current per-call freshness semantics for tests that mutate schemas.
  **Risk:** Stale-schema reuse if a schema file changes in-process (documented; opt-out covers); thread-safety of shared compiled grammars is guaranteed read-only in libxml2 ≥2.9 but must be asserted with a concurrent-validation test.
  **Falsify:** Validate N=100 files against one schema, workers=8: pre/post wall time; if schema parse is <5% of per-call time, drop the proposal.

  ### A11. Harness: identical query sequence in every worker biases indexed-read measurements
  **File/function:** `common/reference tools/xml-parser-benchmark/benchmark.cpp`, `run()` "indexed" mode: `std::mt19937_64 random(20261019);` per worker.
  **Mechanism:** All W workers issue the *same* 2048 offsets in the *same* order. After the first worker touches a page, the remaining W−1 workers are guaranteed cache hits (context already notes cached input, but this is a stronger, per-offset aliasing effect). Linux 96-worker indexed numbers therefore measure near-pure CPU decode with maximal page-cache sharing and understate per-worker I/O divergence; they also cannot detect offset-dependent concurrency bugs.
  **Minimal change:** Seed per worker (`seed ^ (i+1)`) or use disjoint query partitions; keep one fixed-seed configuration for cross-branch digest comparability (digests must be compared against a canonical digest, not worker-vs-worker, if seeds differ).
  **Risk:** Changes comparability with the frozen snapshots — keep the fixed-seed mode as a flag and re-baseline.
  **Falsify:** Compare per-worker-seed vs fixed-seed at 96 workers on Linux; if files/s moves <2%, the bias is negligible and this is downgraded to a nit.

  ---

  ## B. Hypotheses (need measurement to confirm)

  **H1. UTF-16 conversion + per-event allocation explains most of the single-reader stream gap.** Xerces 2.926 s vs lib2 5.536 / E 6.776 / P 5.902 on Mac 1-reader. All three prototypes transcode every event to UTF-16 and heap-allocate per event (A1/A2); Xerces delivers UTF-16 natively with parser-owned buffers. If after A1+A2 the gap does not close materially, the residual lies in `MzMLHandler`'s own consumption (common code) or in `BinInputStream` chunking. Falsify: `perf record`/Instruments profile of one stream run before and after A1+A2.

  **H2. Digest hashing is a non-trivial fraction of indexed-query time, compressing branch differences.** 2048 queries × 2048 points × 2 words × 8 FNV iterations ≈ 67M multiply-xor steps per worker on top of decode. Mac fragment times (.293/.205/.187/.144) may understate true parser differences. Falsify: a `--no-digest` variant run; report parser-only deltas.

  **H3. pugixml internal-DTD entity divergence.** pugixml does not process internal DTD subsets the way Xerces/libxml2-with-`NOENT`/Expat do; documents relying on internal general entities will be rejected or deliver different text in the pugi branch (whichever, it diverges). The benchmark generator never emits DTDs, so digests agree on tested inputs. Falsify: one test file with an internal `<!ENTITY>` used in an attribute and in text; compare across branches.

  **H4. libxml2 `domParseString_` SAX-table surgery is fragile across libxml2 versions.** Overriding only `parser->sax->endElementNs` relies on the default SAX2 table's initialization and on `ctxt->node == NULL` exactly at root end; correct on 2.15.4 (tests pass) but a silent behavior change in another 2.x would turn "stop at root end" into "parse trailing garbage" (which *would* then throw, so failure is loud, not silent). Falsify: run the indexed test suite against the oldest supported libxml2; also add a unit test asserting the trailing parent-closing-tags case explicitly.

  **H5. Expat fragment parser: `XML_Parse` copies each 64 KB chunk.** Using `XML_GetBuffer`/`XML_ParseBuffer` parses the caller's memory in place. For the fragment path (0.187 s / 2048 queries) the copy is small; unlikely to matter. Falsify: fragment benchmark before/after; keep only if ≥2%.

  **H6. `load_buffer_inplace` in the pugi whole-file path makes `input` non-shareable and forces one full extra residency of the file on top of the DOM** — unavoidable in pugi, but `input.reserve(file_size)` (available via `File::size` before streaming) would remove ~24 realloc/memcpy doublings for a 1.05 GiB file. Falsify: time the read loop alone; likely <1% — listed for completeness, not as a recommendation.

  ---

  ## C. Top-3 optimizations ranked by expected impact

  1. **Eliminate per-event UTF-16 conversion and allocation in the SAX adapters (A1 + A2).** Targets the dominant workload (streaming full-file reads) where every prototype is ~1.9–2.3× slower than Xerces (Mac 1-reader: 2.926 vs 5.536/6.776/5.902 s). Mechanism is certain; attribution to the gap is H1 pending a profile. Applies to libxml2 and Expat directly and to pugi's traversal (`as_wide`→`wstring`→`u16string` double conversion there too).
  2. **Resolve the pugixml whole-file memory blowup or formally scope the pugi branch to indexed fragments (A4/A5 + context measurements).** Mac 8-reader stream RSS 26 427 MiB vs ~134–154 MiB for all others; Linux 96-reader 151 GiB vs ~0.9 GiB; whole-load 8-reader peak 33.86 GiB vs ~8.35 and 96-worker ~209.8 GiB vs ~99.8. This is two orders of magnitude and is disqualifying for the stream/load paths; it is architectural (pugi is DOM-only, acknowledged by the "ponytail" comment), so the honest minimal actions are: fix the correctness defects (A3, A4, A5), shave constant factors (`reserve`, single-conversion UTF-16), and gate `XMLFile::parse_` in the pugi branch to non-streaming uses. No within-constraints change makes pugi streaming bounded-memory.
  3. **Replace per-query DOM construction in the libxml2 fragment decoder with its own SAX machinery (A9).** Indexed fragment reads are the only scenario where all prototypes already beat Xerces (0.144–0.205 vs 0.293 s); libxml2 is the slowest of the three there. Reusing the branch's existing push-parser/SAX code (parser unchanged) should close most of the 0.205→~0.19 s gap and, more importantly, removes per-query DOM allocation — the same defect class that makes pugi untenable at scale.

  ## D. Explicitly not recommended
  - Any global lock around parser calls (independent-reader thread safety is a hard requirement; the measured 96-reader scaling of libxml2/Expat confirms none is needed).
  - `XML_PARSE_HUGE`-style safety-limit relaxation beyond the existing bounded use, or disabling amplification limits.
  - Replacing pugixml or hand-rolling a streaming tokenizer (the pugi fragment-end scanner in `domParseString_` — A12-adjacent — is already at the edge of this; its DOCTYPE/`<!…>` handling gap should be closed by explicit rejection, not by extending the scanner).

