import datetime,hashlib,json,os,subprocess,time
from pathlib import Path
root=Path(__file__).resolve().parent
prompt=(root/'review-input.txt').read_text();start=time.time()
command=['vibe','--agent','optimized-review','--legacy-harness','--disabled-tools','*','--max-turns','1','--time-limit','600','--output','streaming','--prompt',prompt]
with (root/'vibe-review-output.txt').open('w') as out,(root/'vibe-review-stderr.txt').open('w') as err:
 result=subprocess.run(command,cwd=root,stdout=out,stderr=err,env=dict(os.environ,VIBE_CLI='python',VIBE_AGENT_PATHS=json.dumps([str(root/'.vibe/agents')])),timeout=900)
record=dict(reviewer='vibe',version='2.26.1',runtime='supported Python CLI, legacy harness',mode='safe custom agent, all tools disabled; same installed default model alias with thinking off',started_utc=datetime.datetime.fromtimestamp(start,datetime.timezone.utc).isoformat(),elapsed_seconds=time.time()-start,returncode=result.returncode,input_sha256=hashlib.sha256(prompt.encode()).hexdigest(),input_bytes=len(prompt.encode()),scope='Same full212224byte frozen input as Claude/Kimi; earlier high-thinking attempts timed out without a final answer. Private reasoning output excluded from publication.',requested_model_name='mistral-vibe-cli-latest',requested_model_alias='mistral-medium-3.5',requested_thinking='off')
(root/'vibe-review-run.json').write_text(json.dumps(record,indent=2)+'\n');print('Vibe focused exit',result.returncode,'elapsed',record['elapsed_seconds'],flush=True)
