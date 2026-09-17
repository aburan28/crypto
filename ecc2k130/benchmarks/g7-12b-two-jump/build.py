from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib,json,subprocess,time
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
BASE=ROOT/'benchmarks/g7-12b-seven-node/builds.json'
SRC=ROOT/'build/g7-two-jump-src/src/main.cu'

def sha(p):
 h=hashlib.sha256();h.update(Path(p).read_bytes());return h.hexdigest()
base=json.loads(BASE.read_text())['rows'][0]['command']

def task(mode):
 label=f'two-jump-{mode}'
 cmd=list(base)
 cmd[cmd.index('/home/ubuntu/crypto/ecc2k130/build/g7-seven-node-src/src/main.cu')]=str(SRC)
 cmd[-2]=str(ROOT/'build'/f'g7-{label}')
 cmd.insert(cmd.index('-Xptxas'),f'-DECC_PACKED_TWO_JUMP_MODE={1 if mode=="indexed" else 2}')
 if mode=='immediate':
  cmd[cmd.index('-DECC_PACKED_SHARED_SIGMA=1')]='-DECC_PACKED_SHARED_SIGMA=0'
 t=time.monotonic();p=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=600)
 log=OUT/f'{label}-build.log';log.write_text(p.stdout)
 binary=ROOT/'build'/f'g7-{label}'
 return {'label':label,'command':cmd,'returncode':p.returncode,'elapsed_seconds':time.monotonic()-t,
         'log':log.name,'log_sha256':sha(log),'binary':str(binary.relative_to(ROOT)),
         'binary_sha256':sha(binary) if p.returncode==0 else None}
with ThreadPoolExecutor(max_workers=2) as pool: rows=list(pool.map(task,['indexed','immediate']))
result={'passed':all(r['returncode']==0 for r in rows),'rows':rows}
(OUT/'builds.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
raise SystemExit(0 if result['passed'] else 1)
