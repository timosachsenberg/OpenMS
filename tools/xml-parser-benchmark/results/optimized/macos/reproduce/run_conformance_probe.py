"""Compare generic OpenMS XML event adapters on a small adversarial corpus."""
import argparse,json,os,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--snapshots',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
cases={'entity-valid':b'<r>A&amp;B</r>','undefined-entity':b'<r>&undefined;</r>',
       'malformed-reference':b'<r>A&-B</r>','duplicate-attribute':b'<r a="1" a="2"/>',
       'invalid-utf8':b'<r>\xff</r>','overlong-utf8':b'<r>\xc0\xaf</r>',
       'truncated-utf8':b'<r>\xe2\x82</r>',
       'latin1':b'<?xml version="1.0" encoding="ISO-8859-1"?><r>\xe9</r>',
       'windows1252':b'<?xml version="1.0" encoding="windows-1252"?><r>\x80</r>',
       'unicode': '<r>β🧬</r>'.encode(),'mismatched-tag':b'<r>\n<x></r>'}
fixtures=a.output.parent/'conformance-fixtures';fixtures.mkdir(parents=True,exist_ok=True)
rows=[]
for name,xml in cases.items():
 fixture=fixtures/(name+'.xml');fixture.write_bytes(xml)
 for backend in ['xerces','libxml2','expat','pugixml']:
  env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1');env['DYLD_LIBRARY_PATH' if os.uname().sysname=='Darwin' else 'LD_LIBRARY_PATH']=str((a.snapshots/backend/'lib').resolve())
  result=subprocess.run([str(a.binary.resolve()),str(fixture.resolve())],env=env,capture_output=True,text=True,timeout=30)
  assert result.returncode==0,(backend,name,result.returncode,result.stderr)
  value=json.loads([s for s in result.stdout.splitlines() if s.startswith('{')][-1])
  rows.append(dict(backend=backend,case=name,stderr=result.stderr,**value))
a.output.write_text(json.dumps({'method':'OpenMS XMLFile generic SAX callbacks, isolated process per case','results':rows},indent=2)+'\n')
