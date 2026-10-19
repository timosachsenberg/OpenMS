"""Run frozen optimized macOS snapshots; never mutate original measurements."""
import argparse,json,subprocess,sys,datetime,hashlib,os
from pathlib import Path
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/optimized/macos'
TOOLS=ROOT/'OpenMS/tools/xml-parser-benchmark'
BACKENDS=['xerces','libxml2','expat','pugixml']

def run(args,log):
    with (OUT/log).open('w') as handle:subprocess.run(list(map(str,args)),cwd=ROOT,stdout=handle,stderr=subprocess.STDOUT,check=True)

def environment(phase):
    record=dict(timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),phase=phase,scope='Boundary snapshot; not a continuous thermal trace.')
    for key,args in [('thermal_status',['pmset','-g','therm']),('swap_status',['sysctl','vm.swapusage'])]:
        value=subprocess.run(args,capture_output=True,text=True);record[key]=dict(returncode=value.returncode,stdout=value.stdout,stderr=value.stderr)
    (OUT/('environment-'+phase+'.json')).write_text(json.dumps(record,indent=2)+'\n')

def main():
    p=argparse.ArgumentParser();p.add_argument('--screen',action='store_true');a=p.parse_args()
    manifest=json.loads((OUT/'source-manifest.json').read_text());assert set(manifest)==set(BACKENDS)
    for backend in BACKENDS:
        for name,digest in manifest[backend]['libraries'].items():
            assert hashlib.sha256((ROOT/'snapshots-optimized'/backend/'lib'/name).read_bytes()).hexdigest()==digest
    # Older build-driver logs combined the startup marker with JSON; keep marker as a separate diagnostic.
    for backend in BACKENDS:
        path=OUT/(backend+'-huge-array-equivalence.json')
        text=path.read_text()
        try:json.loads(text)
        except json.JSONDecodeError:
            (OUT/(backend+'-huge-array-stderr.log')).write_text('\n'.join(s for s in text.splitlines() if not s.startswith('{'))+'\n')
            path.write_text(json.dumps(json.loads(next(s for s in text.splitlines() if s.startswith('{'))),indent=2)+'\n')
    base=[sys.executable,TOOLS/'run_benchmarks.py','--binary',ROOT/'benchmark-optimized-build/xml_parser_benchmark']
    for backend in BACKENDS:base+=['--backend',backend+'='+str(ROOT/'snapshots-optimized'/backend/'lib')]
    inputs=[ROOT/'benchmark-data'/name for name in ['small.mzML','medium.mzML','large.mzML','medium-zlib.mzML','medium-zlib.mzML.gz']]
    inputs.append(ROOT/'OpenMS/share/OpenMS/examples/FRACTIONS/BSA1_F1.mzML')
    selected=[ROOT/'benchmark-data/large.mzML',ROOT/'benchmark-data/medium-zlib.mzML']
    def matrix(files,modes,output,workers='1,2,4,8',reps=5):
        command=base+['--workers',workers,'--modes',modes,'--repetitions',str(reps),'--output',OUT/(output+'.json')]
        for filename in files:command+=['--input',filename]
        run(command,output+'-runs.log')
    if a.screen:
        matrix(selected,'stream,load,indexed,indexed-metadata,contexts,idle','screen',workers='1,8',reps=3)
        return
    environment('before')
    matrix(inputs,'stream,indexed,index-open','benchmarks')
    matrix(selected,'load,contexts,idle','memory-scaling')
    matrix(selected,'indexed-metadata,contexts-metadata','metadata-scaling')
    run([sys.executable,TOOLS/'run_openmp_scaling.py','--matrix',OUT/'benchmarks.json','--output',OUT/'openmp-scaling.json'],'openmp-scaling-runs.log')
    environment('after')
    (OUT/'resource-monitor.log').write_text('macOS thermal and swap boundary snapshots: environment-before.json and environment-after.json. No continuous trace.\n')
    print('Completed optimized macOS matrices.',flush=True)

if __name__=='__main__':main()
