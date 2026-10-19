import hashlib,json
from pathlib import Path
root=Path('${XML_POC_ROOT}')
expected=json.loads((root/'staging/expected-inputs.json').read_text())
observed={}
for name,data in expected.items():
    path=root/'data'/name
    if name=='BSA1_F1.mzML':path=root/'source/share/OpenMS/examples/FRACTIONS'/name
    digest=hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda:handle.read(1024*1024),b''):digest.update(chunk)
    value={'bytes':path.stat().st_size,'sha256':digest.hexdigest()}
    observed[name]=value
    assert value==data,(name,data,value)
(root/'results/input-verification.json').write_text(json.dumps(observed,indent=2)+'\n')
print('Verified all',len(observed),'inputs match macOS bytes/SHA256')
