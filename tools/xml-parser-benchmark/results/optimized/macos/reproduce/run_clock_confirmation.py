"""Keep anomalous clock observations, then confirm affected cells separately."""
import hashlib,json,os,random,statistics,sys
from pathlib import Path
sys.path.insert(0,str(Path('OpenMS/tools/xml-parser-benchmark').resolve()))
from run_benchmarks import measured,sha256
root=Path.cwd(); output=root/'results/optimized/macos/clock-consistency.json'
report={'method':'C++ steady_clock seconds > Python launcher perf_counter seconds + 0.5s. Cause unestablished. Frozen matrices unchanged; independent three-repetition confirmations are not mixed with final medians.','anomalies':[],'confirmations':[]}
cells={}
for directory in ['results','results/linux','results/optimized/macos','results/optimized/linux']:
 for name in ['benchmarks','memory-scaling','metadata-scaling','openmp-scaling','context-residency']:
  p=root/directory/(name+'.json')
  if not p.is_file():continue
  data=json.loads(p.read_text())
  for row in data['records']:
   if row.get('seconds',0)>row.get('process_seconds',float('inf'))+.5:
    key=(str(p.relative_to(root)),row['backend'],row['input'],row.get('mode','stream'),row.get('workers',1),row.get('openmp_threads',1))
    report['anomalies'].append(dict(matrix=key[0],**row));cells[key]=data
jobs=[(key,rep) for key in cells for rep in range(3)];random.Random(20261019).shuffle(jobs)
for n,(key,rep) in enumerate(jobs,1):
 path,backend,filename,mode,workers,threads=key; data=cells[key]; binary=Path(data['binary']); launcher=binary.with_name('xml_benchmark_launcher')
 assert sha256(binary)==data['binary_sha256'] and sha256(launcher)==data['launcher_sha256']
 env=dict(os.environ,DYLD_LIBRARY_PATH=data['backends'][backend],OMP_NUM_THREADS=str(threads),OPENBLAS_NUM_THREADS='1')
 if 'openmp' in path:env['OMP_DYNAMIC']='FALSE'
 row=measured([str(launcher),str(binary),mode,filename,str(workers),'0' if 'openmp' in path else '2048'],env,600)
 expected=next(r for r in data['records'] if r['backend']==backend and r['input']==filename and r.get('mode','stream')==mode)
 assert row['ok'] and row['hash_per_worker']==expected['hash_per_worker']
 assert Path(row['openms_library']).resolve().parent==Path(data['backends'][backend]).resolve()
 assert row['seconds']<=row['process_seconds']+.5,row
 report['confirmations'].append(dict(matrix=path,backend=backend,input=filename,openmp_threads=threads,repetition=rep,**row))
 output.write_text(json.dumps(report,indent=2)+'\n');print(f'{n}/{len(jobs)} {backend} {mode} OMP{threads}: {row["seconds"]:.3f}s',flush=True)
report['summary']=[]
for key,data in cells.items():
 path,b,f,m,w,t=key
 rows=[r for r in report['confirmations'] if (r['matrix'],r['backend'],r['input'],r['mode'],r['workers'],r['openmp_threads'])==key]
 old=[r for r in data['records'] if r['backend']==b and r['input']==f and r.get('mode','stream')==m and r.get('workers',1)==w and r.get('openmp_threads',1)==t]
 report['summary'].append(dict(matrix=path,backend=b,input=f,mode=m,workers=w,openmp_threads=t,original_median=statistics.median(r['seconds'] for r in old),confirmation_median=statistics.median(r['seconds'] for r in rows),confirmation_min=min(r['seconds'] for r in rows),confirmation_max=max(r['seconds'] for r in rows),repetitions=3))
report['complete']=True;output.write_text(json.dumps(report,indent=2)+'\n')
