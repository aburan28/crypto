#!/usr/bin/env python3
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib,json,subprocess,time
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
SRC=ROOT/'build/g7-batch32-split-src'
UNITS={'packedstate':'testpackedstatecuda.cu','blockinverse':'testblockinversecuda.cu','sharedx':'testsharedxcuda.cu'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def task(item):
 batch,unit=item;build=json.loads((OUT/f'batch{batch}-build.json').read_text());cmd=list(build['command'])
 cmd[cmd.index(str(SRC/'src/main.cu'))]=str(SRC/'src'/UNITS[unit]);binary=ROOT/'build'/f'g7-batch{batch}-{unit}-test';cmd[cmd.index(str(ROOT/'build'/f'g7-split-batch{batch}'))]=str(binary)
 log=OUT/f'batch{batch}-{unit}-build.log';receipt=OUT/f'batch{batch}-{unit}-build.json';assert not binary.exists() and not log.exists() and not receipt.exists()
 start=time.monotonic()
 with log.open('x') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=600)
 row={'batch':batch,'unit':unit,'command':cmd,'returncode':r.returncode,'elapsed_seconds':time.monotonic()-start,'log_sha256':sha(log),'binary_sha256':sha(binary) if binary.exists() else None};receipt.write_text(json.dumps(row,indent=2)+'\n');return row
items=[(b,u) for b in (24,32) for u in UNITS]
with ThreadPoolExecutor(max_workers=2) as pool: rows=list(pool.map(task,items))
print(json.dumps(rows,indent=2));assert all(r['returncode']==0 for r in rows)
