"""Repeat actual OpenMS XMLValidator API compatibility probes on frozen snapshots."""
import argparse,hashlib,json,os,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--binary',required=True,type=Path);p.add_argument('--snapshots',required=True,type=Path);p.add_argument('--cases',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
rows=json.loads(a.cases.read_text());rows=rows['results'] if isinstance(rows,dict) else rows
results=[]
for row in rows:
 env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1');env['DYLD_LIBRARY_PATH' if os.uname().sysname=='Darwin' else 'LD_LIBRARY_PATH']=str((a.snapshots/row['backend']/'lib').resolve())
 outcome=subprocess.run([str(a.binary.resolve()),row['file'],row['schema']],env=env,capture_output=True,text=True,timeout=120)
 value=json.loads(next(s for s in outcome.stdout.splitlines() if s.startswith('{')))
 expected=not(row['backend']=='libxml2' and row['case']=='anyuri-windows')
 assert value['valid']==expected and outcome.returncode==int(not expected),(row,outcome)
 results.append(dict(backend=row['backend'],case=row['case'],returncode=outcome.returncode,stderr=outcome.stderr,**value))
a.output.write_text(json.dumps({'binary_sha256':hashlib.sha256(a.binary.read_bytes()).hexdigest(),'results':results},indent=2)+'\n')
