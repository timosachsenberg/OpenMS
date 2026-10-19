# libxml2 proof of concept

Base: OpenMS develop `77cccda7` (OpenMS 3.7). libxml2 >=2.13 is required for `XML_PARSE_NO_XXE`; tested library is Homebrew 2.15.4.

Implemented paths: common XMLFile SAX parsing, XMLValidator streaming XSD validation, IndexedMzMLDecoder footer, MzMLSpectrumDecoder spectrum/chromatogram fragments. Existing decompression InputSources, StringManager UTF16 helpers, initialization, and mzIdentML DOM still use Xerces. This branch does **not** eliminate the Xerces dependency.

`XMLParserBackend_test` exercises Unicode, entity expansion, CDATA, crossing the 64 KiB input boundary, malformed XML, handler exceptions, early cancellation, and eight independent parser contexts. Existing format and indexed-reader class tests exercise actual format decoding and compressed input. libxml2 contexts are local to each operation and no global cleanup is called.

The SAX1 adapter preserves the current UTF16 callback ABI. It converts each event and attribute; this is a compatibility POC, not a direct UTF8 handler redesign. It expands internal entities while preventing external entity loading (`NO_XXE`) and network access (`NONET`). Unlike the old parser, documents depending on external entities are intentionally unsupported. DOM fragment size is limited to `INT_MAX` by `xmlReadMemory`; the libxml2 default text-node resource limits also apply. No relaxed `HUGE` parsing mode is used. These limits need a policy decision before adoption for very large single-spectrum binary arrays.

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
