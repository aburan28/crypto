#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,os,subprocess,time
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
rows=[]
for batch in (24,32):
 for tool in ('memcheck','initcheck','synccheck'):
  binary=ROOT/'build'/f'g7-split-batch{batch}';log=OUT/f'batch{batch}-{tool}.log';assert not log.exists()
  command=['compute-sanitizer','--tool',tool,str(binary),'--packed','--threads','129','--steps','2','--launches','1','--verify','0','--run-id',str(58341+batch),'--bench']
  start=time.monotonic()
  with log.open('x') as f:r=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,timeout=300,env=dict(os.environ,CUDA_DISABLE_PTX_JIT='1',OMP_NUM_THREADS='8'))
  text=log.read_text();row={'batch':batch,'tool':tool,'command':command,'returncode':r.returncode,'elapsed_seconds':time.monotonic()-start,'passed':r.returncode==0 and 'ERROR SUMMARY: 0 errors' in text and '0 dropped' in text,'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest()};rows.append(row);print(json.dumps(row),flush=True);assert row['passed'],text[-5000:]
(OUT/'sanitizers.json').write_text(json.dumps(rows,indent=2)+'\n')
