#!/usr/bin/env python3
"""Compile every bundled XSD with the selected libxml2; emit JSON for comparison."""
import argparse
import ctypes
import json
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--library', default='/opt/homebrew/opt/libxml2/lib/libxml2.dylib')
parser.add_argument('--xerces-probe', type=Path, help='Compile xerces_schema_probe.cpp and pass its executable')
parser.add_argument('--allow-known-differences', action='store_true', help='Accept only the two documented anyURI disagreements with their recorded outcomes')
parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[2])
args = parser.parse_args()
if args.allow_known_differences and not args.xerces_probe:
    parser.error("--allow-known-differences requires --xerces-probe")
lib = ctypes.CDLL(args.library)
lib.xmlSchemaNewParserCtxt.argtypes = [ctypes.c_char_p]
lib.xmlSchemaNewParserCtxt.restype = ctypes.c_void_p
lib.xmlSchemaParse.argtypes = [ctypes.c_void_p]
lib.xmlSchemaParse.restype = ctypes.c_void_p
lib.xmlSchemaFree.argtypes = [ctypes.c_void_p]
lib.xmlSchemaFreeParserCtxt.argtypes = [ctypes.c_void_p]
lib.xmlSchemaNewValidCtxt.argtypes = [ctypes.c_void_p]
lib.xmlSchemaNewValidCtxt.restype = ctypes.c_void_p
lib.xmlSchemaValidateFile.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int]
lib.xmlSchemaFreeValidCtxt.argtypes = [ctypes.c_void_p]
results = []
for path in sorted((args.root / 'share/OpenMS/SCHEMAS').rglob('*.xsd')):
    ctx = lib.xmlSchemaNewParserCtxt(bytes(path))
    schema = lib.xmlSchemaParse(ctx)
    results.append({'schema': str(path.relative_to(args.root)), 'compiled': bool(schema), 'xerces_compiled': subprocess.run([str(args.xerces_probe), str(path)]).returncode == 0 if args.xerces_probe else None})
    if schema:
        lib.xmlSchemaFree(schema)
    lib.xmlSchemaFreeParserCtxt(ctx)
# Exercise the repository's validator fixtures, including malformed XML.
data = args.root / 'src/tests/class_tests/openms/data'
corpus = [('XMLValidator.xsd', path.name, path.stem == 'XMLValidator_valid') for path in sorted(data.glob('XMLValidator_*.xml'))]
corpus += [
    ('FeatureXML_1_9.xsd', 'FeatureXMLFile_1.featureXML', True),
    ('FeatureXML_1_9.xsd', 'FeatureXMLFile_2_options.featureXML', True),
    ('ConsensusXML_1_7.xsd', 'ConsensusXMLFile_1.consensusXML', True),
    ('ConsensusXML_1_7.xsd', 'ConsensusXMLFile_2_options.consensusXML', True),
    ('TrafoXML_1_1.xsd', 'TransformationXMLFile_1.trafoXML', True),
    ('TrafoXML_1_1.xsd', 'TransformationXMLFile_2.trafoXML', True),
    ('TrafoXML_1_1.xsd', 'TransformationXMLFile_3.trafoXML', False),
    ('TrafoXML_1_1.xsd', 'TransformationXMLFile_4.trafoXML', True),
    ('mzML_idx_1_10.xsd', 'MzMLFile_4_indexed.mzML', True),
]
# Additional instances are observational: their historical schema versions can
# differ from the bundled current grammar. Xerces is the baseline oracle; do
# not label an old fixture schema-valid solely because OpenMS can read it.
patterns = {
    '*.idXML': 'IdXML_1_5.xsd', '*.mzData': 'mzData_1_05.xsd',
    '*.mzXML': 'mzXML_idx_3.1.xsd', '*.mzid': 'mzIdentML1.3.0.xsd',
    '*.pepXML': 'pepXML_v114.xsd', '*.pepxml': 'pepXML_v114.xsd',
    '*.protXML': 'protXML_v6.xsd', '*.traML': 'TraML1.0.0.xsd',
    '*.mzML': 'mzML_1_10.xsd', '*.ini': 'Param_1_8_0.xsd',
    'ParamXMLFile*.xml': 'Param_1_8_0.xsd',
}
seen = {filename for _, filename, _ in corpus}
for pattern, default_schema in patterns.items():
    for path in sorted(data.glob(pattern)):
        if path.name in seen:
            continue
        schema_name = default_schema
        try:
            root = next(ET.iterparse(path, events=('start',)))[1]
            root_name = root.tag.rsplit('}', 1)[-1]
            if root_name == 'indexedmzML':
                schema_name = 'mzML_idx_1_10.xsd'
            elif root_name == 'MzIdentML':
                version = root.get('version', '1.3.0')
                candidate = f'mzIdentML{version}.xsd'
                if (args.root / 'share/OpenMS/SCHEMAS' / candidate).exists():
                    schema_name = candidate
        except (ET.ParseError, StopIteration):
            pass
        corpus.append((schema_name, path.name, None))
        seen.add(path.name)
corpus.extend([('tools/xml-poc/fixtures/anyuri.xsd', 'tools/xml-poc/fixtures/anyuri-windows.xml', None),
               ('tools/xml-poc/fixtures/anyuri.xsd', 'tools/xml-poc/fixtures/anyuri-normalized.xml', True)])
validations = []
for schema_name, filename, expected in corpus:
    schema_path = args.root / schema_name if '/' in schema_name else data / schema_name if schema_name == 'XMLValidator.xsd' else args.root / 'share/OpenMS/SCHEMAS' / schema_name
    path = args.root / filename if '/' in filename else data / filename
    ctx = lib.xmlSchemaNewParserCtxt(bytes(schema_path))
    schema = lib.xmlSchemaParse(ctx)
    if not schema:
        raise RuntimeError(f'Cannot compile {schema_path}')
    valid_ctx = lib.xmlSchemaNewValidCtxt(schema)
    status = lib.xmlSchemaValidateFile(valid_ctx, bytes(path), 0)
    xerces = subprocess.run([str(args.xerces_probe), str(schema_path), str(path)]).returncode if args.xerces_probe else None
    validations.append({'file': filename, 'schema': schema_name, 'status': status, 'expected_valid': expected, 'xerces_status': xerces})
    lib.xmlSchemaFreeValidCtxt(valid_ctx)
    lib.xmlSchemaFree(schema)
    lib.xmlSchemaFreeParserCtxt(ctx)
mismatches = [v for v in validations if v['xerces_status'] is not None and (v['status'] == 0) != (v['xerces_status'] == 0)]
known = {
    'TargetedSpectraExtractor_13C1_spectra0to100.mzML': 'Windows file:// drive/backslash xs:anyURI lexical handling',
    'tools/xml-poc/fixtures/anyuri-windows.xml': 'Minimal reproducer of the same xs:anyURI difference',
}
for mismatch in mismatches:
    mismatch['known_reason'] = known.get(mismatch['file']) if mismatch['status'] == 1824 and mismatch['xerces_status'] == 0 else None
summary = {'schema_count': len(results), 'instance_count': len(validations), 'libxml_valid': sum(v['status'] == 0 for v in validations),
           'xerces_valid': sum(v['xerces_status'] == 0 for v in validations), 'mismatches': mismatches}
print(json.dumps({'library': args.library, 'summary': summary, 'schemas': results, 'validation': validations}, indent=2))
raise SystemExit(any(not r['compiled'] or r['xerces_compiled'] is False for r in results) or any((v['status'] == 0) != v['expected_valid'] for v in validations if v['expected_valid'] is not None)
                 or any(not (args.allow_known_differences and v['known_reason']) for v in mismatches))
