import datetime,hashlib,json,os,subprocess,time
from pathlib import Path
root=Path(__file__).resolve().parent
prompt=(root/'vibe-focused-input.txt').read_text();start=time.time()
command=['vibe','--agent','ask','--legacy-harness','--disabled-tools','*','--max-turns','1','--time-limit','600','--output','streaming','--prompt',prompt]
with (root/'vibe-review-output.txt').open('w') as out,(root/'vibe-review-stderr.txt').open('w') as err:
 result=subprocess.run(command,cwd=root,stdout=out,stderr=err,env=dict(os.environ,VIBE_CLI='python'),timeout=900)
record=dict(reviewer='vibe',version='2.26.1',runtime='supported Python CLI, legacy harness',mode='ask agent, all tools disabled',started_utc=datetime.datetime.fromtimestamp(start,datetime.timezone.utc).isoformat(),elapsed_seconds=time.time()-start,returncode=result.returncode,input_sha256=hashlib.sha256(prompt.encode()).hexdigest(),input_bytes=len(prompt.encode()),scope='Focused current optimization helpers across all three adapters; earlier full212224byte review exceeded time limit without final findings. Private reasoning output excluded from publication.')
(root/'vibe-review-run.json').write_text(json.dumps(record,indent=2)+'\n');print('Vibe focused exit',result.returncode,'elapsed',record['elapsed_seconds'],flush=True)
