#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,subprocess,time
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
rows=[]
for batch in (24,32):
 for unit in ('packedstate','blockinverse','sharedx'):
  binary=ROOT/'build'/f'g7-batch{batch}-{unit}-test';log=OUT/f'batch{batch}-{unit}.log';assert not log.exists();start=time.monotonic()
  with log.open('x') as f:r=subprocess.run([str(binary)],stdout=f,stderr=subprocess.STDOUT,timeout=300)
  text=log.read_text();row={'batch':batch,'unit':unit,'returncode':r.returncode,'elapsed_seconds':time.monotonic()-start,'passed':r.returncode==0 and 'PASS' in text and 'FAIL' not in text,'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest()};rows.append(row);print(json.dumps(row),flush=True);assert row['passed'],text[-4000:]
(OUT/'unit-runs.json').write_text(json.dumps(rows,indent=2)+'\n')
