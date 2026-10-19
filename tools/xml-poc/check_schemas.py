#!/usr/bin/env python3
"""Compile every bundled XSD with the selected libxml2; emit JSON for comparison."""
import argparse
import ctypes
import json
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--library', default='/opt/homebrew/opt/libxml2/lib/libxml2.dylib')
parser.add_argument('--xerces-probe', type=Path, help='Compile xerces_schema_probe.cpp and pass its executable')
parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[2])
args = parser.parse_args()
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
validations = []
for schema_name, filename, expected in corpus:
    schema_path = data / schema_name if schema_name == 'XMLValidator.xsd' else args.root / 'share/OpenMS/SCHEMAS' / schema_name
    path = data / filename
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
print(json.dumps({'library': args.library, 'schemas': results, 'validation': validations}, indent=2))
raise SystemExit(any(not r['compiled'] for r in results) or any((v['status'] == 0) != v['expected_valid'] for v in validations)
                 or any((v['status'] == 0) != (v['xerces_status'] == 0) for v in validations if v['xerces_status'] is not None))
