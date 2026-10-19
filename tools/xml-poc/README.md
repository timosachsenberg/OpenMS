# libxml2 proof of concept

Base: OpenMS develop `77cccda7` (OpenMS 3.7). libxml2 >=2.13 is required for `XML_PARSE_NO_XXE`; tested library is Homebrew 2.15.4.

Implemented paths: common XMLFile SAX parsing, XMLValidator streaming XSD validation, IndexedMzMLDecoder footer, MzMLSpectrumDecoder spectrum/chromatogram fragments. Existing decompression InputSources, StringManager UTF16 helpers, initialization, and mzIdentML DOM still use Xerces. This branch does **not** eliminate the Xerces dependency.

`XMLParserBackend_test` exercises Unicode, entity expansion, CDATA, crossing the 64 KiB input boundary, malformed XML, handler exceptions, early cancellation, and eight independent parser contexts, per-reader encoding overrides, external entity blocking, huge binary text, entity-amplification rejection, and propagation of errors from a throwing validator output stream. Existing format and indexed-reader class tests exercise actual format decoding and compressed input. libxml2 contexts are local to each operation and no global cleanup is called.

The SAX1 adapter preserves the current UTF16 callback ABI. It converts each event and attribute; this is a compatibility POC, not a direct UTF8 handler redesign. It expands internal entities while preventing external entity loading (`NO_XXE`) and network access (`NONET`). Unlike the old parser, documents depending on external entities are intentionally unsupported. DOM fragment size is limited to `INT_MAX` by the memory API. Spectrum/chromatogram fragments use `HUGE` to accommodate base64 text larger than 10 MB (covered by a 13 MiB decoder regression). This also raises element depth from 256 to 2048, entity depth from 20 to 40, names from 50 KB to 10 MB, and text from 10 MB to 1 GB. The independent entity amplification factor is explicitly fixed at 5 and tested with an exponential entity input; external entities remain disabled. Footer parsing keeps the ordinary limits. These values follow the [upstream parser documentation](https://gnome.pages.gitlab.gnome.org/libxml2/html/parser_8h.html).

## Reproduce schema compatibility check

```sh
c++ -std=c++17 -I/opt/homebrew/include tools/xml-poc/xerces_schema_probe.cpp -L/opt/homebrew/lib -lxerces-c -o /tmp/xerces_schema_probe
python3 tools/xml-poc/check_schemas.py --xerces-probe /tmp/xerces_schema_probe > schema-results.json 2> schema-diagnostics.txt
```

The checked result is `schema-results.json`: all 22 bundled XSDs compile with both libraries; both agree with the existing class-test expectations for all 13 selected fixtures (valid and invalid XMLValidator, featureXML, consensusXML, trafoXML, indexed mzML). This establishes compatibility for this corpus, not exhaustive XSD equivalence. Diagnostic wording/error codes differ. The standalone comparator uses the same namespace/schema/full-checking settings as the original XMLValidator.

Configure OpenMS with:

```sh
cmake -S . -B ../OpenMS-build -DLIBXML2_INCLUDE_DIR=/opt/homebrew/opt/libxml2/include/libxml2 -DLIBXML2_LIBRARY=/opt/homebrew/opt/libxml2/lib/libxml2.dylib
cmake --build ../OpenMS-build --target XMLParserBackend_test XMLValidator_test MzMLFile_test MzMLSpectrumDecoder_test IndexedMzMLDecoder_test
ctest --test-dir ../OpenMS-build -R '^(XMLParserBackend|XMLValidator|MzMLFile|MzMLSpectrumDecoder|IndexedMzMLDecoder)_test$' --output-on-failure
```

Additional generated-file and scaling benchmarks are coordinated by the shared `tools/xml-parser-benchmark` harness and the parent handoff document. Raw benchmark results are reported there rather than fabricated in this branch.
