"""Read-only observation of native threads and effective OMP environment."""
import json,time
from pathlib import Path
root=Path('${XML_POC_ROOT}')
out=root/'results/openmp-runtime-observation.json'
samples=json.loads(out.read_text())['samples'] if out.exists() else []
for attempt in range(55):
 for p in Path('/proc').iterdir():
  if not p.name.isdigit(): continue
  try:
   args=(p/'cmdline').read_bytes().split(b'\0')
   if not args or args[0].decode(errors='replace')!=str(root/'benchmark-build/xml_parser_benchmark'): continue
   env=dict(x.split(b'=',1) for x in (p/'environ').read_bytes().split(b'\0') if b'=' in x)
   omp={k.decode():v.decode() for k,v in env.items() if k.startswith((b'OMP_',b'GOMP_',b'KMP_'))}
   threads=int(next(x.split()[1] for x in (p/'status').read_text().splitlines() if x.startswith('Threads:')))
   row=dict(pid=int(p.name),command=[a.decode(errors='replace') for a in args if a],openmp_environment=omp,threads_observed=threads)
   samples.append(row)
   if omp.get('OMP_NUM_THREADS')=='96' and threads>=96:
    out.write_text(json.dumps(dict(samples=samples),indent=2)+'\n'); print(json.dumps(row)); raise SystemExit
  except (FileNotFoundError,ProcessLookupError,PermissionError): pass
 time.sleep(1)
out.write_text(json.dumps(dict(samples=samples),indent=2)+'\n')
print('No OMP=96 row observed in sampling window')
