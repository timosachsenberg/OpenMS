# libxml2 proof of concept

Base: OpenMS develop `77cccda7` (OpenMS 3.7). libxml2 >=2.13 is required for `XML_PARSE_NO_XXE`; tested library is Homebrew 2.15.4.

Implemented paths: common XMLFile SAX parsing, XMLValidator streaming XSD validation, IndexedMzMLDecoder footer, MzMLSpectrumDecoder spectrum/chromatogram fragments. Existing decompression InputSources, StringManager UTF16 helpers, initialization, and mzIdentML DOM still use Xerces. This branch does **not** eliminate the Xerces dependency.

`XMLParserBackend_test` exercises Unicode, entity expansion, CDATA, crossing the 64 KiB input boundary, malformed XML, handler exceptions, early cancellation, and eight independent parser contexts, per-reader encoding overrides, external entity blocking, huge binary text, entity-amplification rejection, and propagation of errors from a throwing validator output stream. Existing format and indexed-reader class tests exercise actual format decoding and compressed input. libxml2 contexts are local to each operation and no global cleanup is called.

The SAX1 adapter preserves the current UTF16 callback ABI. It converts each event and attribute; this is a compatibility POC, not a direct UTF8 handler redesign. It expands internal entities while preventing external entity loading (`NO_XXE`) and network access (`NONET`). Unlike the old parser, documents depending on external entities are intentionally unsupported. DOM fragment size is limited to `INT_MAX` by the memory API. Spectrum/chromatogram fragments use `HUGE` to accommodate base64 text larger than 10 MB (covered by a 13 MiB decoder regression). This also raises element depth from 256 to 2048, entity depth from 20 to 40, names from 50 KB to 10 MB, and text from 10 MB to 1 GB. The independent entity amplification factor is explicitly fixed at 5 and tested with an exponential entity input; external entities remain disabled. Footer parsing keeps the ordinary limits. Indexed fragment DOM stops at the completed XML root using a SAX2 end callback, because offset spans can contain following parent-close/list-open tags; malformed content inside the root is rejected. These values follow the [upstream parser documentation](https://gnome.pages.gitlab.gnome.org/libxml2/html/parser_8h.html).

## Reproduce schema compatibility check

```sh
c++ -std=c++17 -I/opt/homebrew/include tools/xml-poc/xerces_schema_probe.cpp -L/opt/homebrew/lib -lxerces-c -o /tmp/xerces_schema_probe
python3 tools/xml-poc/check_schemas.py --xerces-probe /tmp/xerces_schema_probe > schema-results.json 2> schema-diagnostics.txt
```

The checked result is `schema-results.json`: all 22 bundled XSDs compile with both libraries. The original 13 fixtures with known class-test validation expectations all agree. The expanded corpus covers 193 existing XML instances across 16 bundled schemas (mzML/indexed mzML, idXML, mzData, mzXML, mzIdentML 1.1/1.2/1.3, pepXML, protXML, TraML, parameters, featureXML, consensusXML, trafoXML, and validator fixtures). Of those 193, both accept 127 and both reject 65; one disagrees. Historical fixtures can target older schema versions, so shared rejection does not indicate a reader failure. Additional entries have `expected_valid: null` and compare the two validators rather than assuming that every readable fixture is schema-valid.

**Validation compatibility difference:** `TargetedSpectraExtractor_13C1_spectra0to100.mzML:16` has a Windows-style `sourceFile@location` beginning `file://F:\Users\...` with spaces/backslashes. Xerces accepts it; libxml2 reports error 1824 for `xs:anyURI`. This affects explicit XSD validation, not ordinary streaming/indexed reading. The two small instances in `fixtures/` isolate the same difference and show that a normalized `file:///F:/Users/...` URI passes both. The complete matrix therefore has 195 cases and two disagreements representing the same issue.

The [XSD 1.0 lexical rule](https://www.w3.org/TR/xmlschema-2/#anyURI) applies URI escaping before testing RFC 2396/2732 syntax; it explicitly allows spaces in principle. The drive/authority/backslash form needs a standards/conformance decision before choosing whether to retain Xerces tolerance, normalize output, or accept stricter validation. We do not classify the entire input file as invalid XML or relax validation to hide the difference.

The strict comparison exits 1 for the recorded mismatch. A regression gate allowing only the explicitly listed differences (and only with libxml2 status 1824 versus Xerces status 0) is:

```sh
python3 tools/xml-poc/check_schemas.py --xerces-probe /tmp/xerces_schema_probe --allow-known-differences
```

That gate passes; a new difference, reversed result, any failed schema compilation, or failed known expectation still fails. This establishes compatibility for this corpus, not exhaustive XSD equivalence. Diagnostic wording/error codes differ. The standalone comparator uses namespace/schema/full-checking settings matching the original XMLValidator, but disables loading additional instance-referenced remote schemas after caching the supplied grammar. This ensures both validators compare the explicitly supplied bundled schema without network fallback.

Configure OpenMS with:

```sh
cmake -S . -B ../OpenMS-build -DLIBXML2_INCLUDE_DIR=/opt/homebrew/opt/libxml2/include/libxml2 -DLIBXML2_LIBRARY=/opt/homebrew/opt/libxml2/lib/libxml2.dylib
cmake --build ../OpenMS-build --target XMLParserBackend_test XMLValidator_test MzMLFile_test MzMLSpectrumDecoder_test IndexedMzMLDecoder_test
ctest --test-dir ../OpenMS-build -R '^(XMLParserBackend|XMLValidator|MzMLFile|MzMLSpectrumDecoder|IndexedMzMLDecoder)_test$' --output-on-failure
```

Additional generated-file and scaling benchmarks are coordinated by the shared `tools/xml-parser-benchmark` harness and the parent handoff document. Raw benchmark results are reported there rather than fabricated in this branch.
