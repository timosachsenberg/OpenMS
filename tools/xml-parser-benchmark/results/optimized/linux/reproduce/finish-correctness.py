import hashlib,json,os,subprocess
from pathlib import Path
root=Path('${XML_POC_ROOT}/optimized')
manifest=json.loads((root/'staging/source-manifest.json').read_text())
for backend in manifest:
 manifest[backend]['libraries']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (root/'snapshots'/backend/'lib').iterdir() if p.is_file()}
for backend in manifest:
 manifest[backend]['source_tree']=subprocess.check_output(['git','-C',str(root/'source'),'rev-parse',manifest[backend]['source_commit']+'^{tree}'],text=True).strip()
 manifest[backend]['tracked_source_changes']=(root/'results'/f'{backend}-source-status.txt').read_text()
 manifest[backend]['timed_source_sha256']={name:hashlib.sha256((root/'staging/benchmark-source'/name).read_bytes()).hexdigest() for name in ('benchmark.cpp','launcher.cpp','run_benchmarks.py','run_openmp_scaling.py')}
 manifest[backend]['build_configuration']={'flags':'Release -O3 -DNDEBUG, C++23; -I private dependencies/header-priority', 'source_patches_outside_commit':False, 'dependency_workarounds':['Private parser header priority overlay ahead of broad contrib include','Original contrib CURL::libcurl_static config for Arrow'], 'shared_dependency_prefix':'Original private dependencies installation, unchanged'}
(root/'results/source-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
subprocess.run(['python3',str(root/'staging/run_conformance_probe.py'),'--binary',str(root/'benchmark-build/xml_conformance_probe'),'--snapshots',str(root/'snapshots'),'--output',str(root/'results/conformance-results.json')],check=True)
subprocess.run(['python3',str(root/'staging/validate-schema-api.py')],check=True)
