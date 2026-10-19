import datetime,hashlib,json,os,platform,subprocess,sys
from pathlib import Path
root=Path('${XML_POC_ROOT}')
phase=sys.argv[1] if len(sys.argv)>1 else 'before'
def command(*args):
    r=subprocess.run(args,capture_output=True,text=True)
    return {'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr}
report={'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'phase':phase,'host':platform.node(),'platform':platform.platform(),'machine':platform.machine(),'python':sys.version,'locale_environment':{key:os.environ.get(key) for key in ('LANG','LC_ALL','LC_CTYPE')},'locale':command('locale'),'thread_environment':{key:os.environ.get(key) for key in ('OMP_NUM_THREADS','OMP_DYNAMIC','OPENBLAS_NUM_THREADS','OMP_PROC_BIND','OMP_PLACES','GOMP_CPU_AFFINITY','KMP_AFFINITY')},'cpu_affinity':sorted(os.sched_getaffinity(0)),'cpu_count':os.cpu_count(),'lscpu':command('lscpu','-J'),'topology':command('lscpu','-e=CPU,CORE,SOCKET,NODE,ONLINE'),'compiler':command('g++','--version'),'cmake':command('cmake','--version'),'uptime':command('uptime'),'loadavg':Path('/proc/loadavg').read_text(),'memory':Path('/proc/meminfo').read_text(),'vmstat':Path('/proc/vmstat').read_text(),'mounts':command('df','-h','/scratch','/dev/shm'),'source_manifest':json.loads((root/'staging/source-manifest.json').read_text()),'versions':{'xerces':'3.3.0','libxml2':'2.15.4','expat':'2.9.0','pugixml':'1.16','zstd':'1.5.7','arrow':'23.0.0','boost':'1.87.0','eigen':'5.0.1','curl':'8.12.1'},'optimization':'Release -O3 -DNDEBUG, C++23; GCC13.3.0','benchmark_affinity_plan':'CPUs0-95:64physicalcores+32SMTsiblings,oneNUMAnode;96workersincludesSMT','temporary_root':str(root)}
(root/'results'/f'environment-{phase}.json').write_text(json.dumps(report,indent=2)+'\n')
print(report['timestamp_utc'],phase,report['loadavg'].strip(),report['cpu_affinity'])
