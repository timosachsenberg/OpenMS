import hashlib,json,os
from pathlib import Path
root=Path('${XML_POC_ROOT}')
entries=json.loads((root/'staging/base-tree-objects.json').read_text())
failures=[]
for mode,expected,name in entries:
    if name=='.git_archival.txt':
        continue  # git archive's intentional export-subst version expansion
    path=root/'source'/name
    if not path.exists() and not path.is_symlink():
        failures.append([name,'missing']);continue
    content=os.fsencode(os.readlink(path)) if mode=='120000' else path.read_bytes()
    observed=hashlib.sha1(f'blob {len(content)}\0'.encode()+content).hexdigest()
    if expected!=observed:failures.append([name,expected,observed])
report=dict(base='77cccda7caf40b46a4919d33748450270eccd565',checked=len(entries)-1,skipped_export_subst=['.git_archival.txt'],failures=failures)
(root/'results/baseline-source-verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
raise SystemExit(bool(failures))
