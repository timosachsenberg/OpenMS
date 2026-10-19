"""Run the user-requested reviewers with a fixed self-contained, review-only input."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

root=Path(__file__).resolve().parent
name=sys.argv[1]
prompt=(root/'review-input.txt').read_text()
commands={
    'vibe':['vibe','--agent','ask','--legacy-harness','--disabled-tools','*','--max-turns','2','--time-limit','900','--output','streaming','--prompt',prompt],
    'kimi':['kimi','--agent-file',str(root/'kimi-reviewer.md'),'--output-format','text','--prompt',prompt],
    'claude':['claude','--print','--safe-mode','--permission-mode','plan','--permission-prompts','none',
              '--tools','','--strict-mcp-config','--mcp-config','{"mcpServers":{}}',
              '--no-session-persistence','--output-format','json'],
}
version=subprocess.run([name,'--version'],capture_output=True,text=True)
started=time.time()
record={'reviewer':name,'version':version.stdout.strip(),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'input_sha256':hashlib.sha256(prompt.encode()).hexdigest(),'input_bytes':len(prompt.encode()),
        'mode':'self-contained review; Claude/Kimi tools disabled, Vibe ask mode; no automatic permission bypass'}
with (root/(name+'-review-output.txt')).open('w') as out, (root/(name+'-review-stderr.txt')).open('w') as err:
    result=subprocess.run(commands[name],input=prompt if name=='claude' else None,text=True,cwd=root,
                          stdout=out,stderr=err,timeout=1800,env=dict(os.environ,VIBE_CLI="python") if name=="vibe" else None)
record.update(returncode=result.returncode,elapsed_seconds=time.time()-started)
(root/(name+'-review-run.json')).write_text(json.dumps(record,indent=2)+'\n')
print(name,'exit',result.returncode,'elapsed',round(record['elapsed_seconds'],1),'output_bytes',(root/(name+'-review-output.txt')).stat().st_size,flush=True)
