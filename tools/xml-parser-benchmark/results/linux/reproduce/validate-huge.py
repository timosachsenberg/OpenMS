import json,os,subprocess,xml.etree.ElementTree as ET
from pathlib import Path
root=Path('${XML_POC_ROOT}')
filename=root/'data/huge-array.mzML'
name=next(ET.iterparse(filename,events=('start',)))[1].tag.rsplit('}',1)[-1]
schema=root/'source/share/OpenMS/SCHEMAS'/('mzML_idx_1_10.xsd' if name=='indexedmzML' else 'mzML_1_10.xsd')
probe=root/'validator-probe-build/xml_schema_validation_probe'
rows=[]
cases=[('huge-array',filename,schema),
       ('anyuri-windows',root/'staging/schema-fixtures/anyuri-windows.xml',root/'staging/schema-fixtures/anyuri.xsd'),
       ('anyuri-normalized',root/'staging/schema-fixtures/anyuri-normalized.xml',root/'staging/schema-fixtures/anyuri.xsd')]
for backend in ('xerces','libxml2','expat','pugixml'):
    snapshot=root/'snapshots'/backend/'lib'
    env=dict(os.environ,LD_LIBRARY_PATH=str(snapshot)+os.pathsep+os.environ.get('LD_LIBRARY_PATH',''))
    linkage=subprocess.run(['ldd',str(probe)],env=env,capture_output=True,text=True,check=True).stdout
    openms_line=next(line for line in linkage.splitlines() if 'libOpenMS.so' in line)
    loaded=Path(openms_line.split('=>',1)[1].strip().split()[0]).resolve()
    assert loaded.parent==snapshot.resolve(),openms_line
    for case,xml,xsd in cases:
        command=[str(probe),str(xml),str(xsd)]
        result=subprocess.run(command,env=env,capture_output=True,text=True)
        row=dict(backend=backend,case=case,command=command,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr,openms_library=str(loaded))
        if result.stdout.strip(): row.update(json.loads(result.stdout.strip().splitlines()[-1]))
        expected=not (backend=='libxml2' and case=='anyuri-windows')
        row['expected_valid']=expected
        rows.append(row)
report=dict(kind='Actual OpenMS XMLValidator API on generated huge mzML and anyURI cases',input=str(filename),schema=str(schema),probe=str(probe),results=rows)
(root/'results/huge-schema-api-results.json').write_text(json.dumps(report,indent=2)+'\n')
assert all(r.get('valid')==r['expected_valid'] and r['returncode']==(0 if r['expected_valid'] else 1) for r in rows),report
print('Actual XMLValidator: all four accept huge-array and normalized URI; only libxml2 rejects Windows URI as expected')
