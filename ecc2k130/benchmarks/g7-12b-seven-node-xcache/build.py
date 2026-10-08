#!/usr/bin/env python3
import hashlib,json,pathlib,subprocess,time
ROOT=pathlib.Path(__file__).resolve().parents[2]
EXP=pathlib.Path(__file__).resolve().parent
SRC=ROOT/'build/g7-cache-six-src'
NVCC=ROOT/'build/assembler134-toolchain/bin/nvcc'
common=[str(NVCC),'-O3','-std=c++17','-gencode','arch=compute_120,code=sm_120',
 '-DECC_BATCH=16','-DECC_THREADS=256','-DECC_MINBLOCKS=3','-DECC_PACKED_SINGLE_PRODUCT=0',
 '-DECC_PACKED_CACHE_DENOM=1','-DECC_PACKED_BY_VALUE=1','-DECC_PACKED_PERM_SIGMA=3',
 '-DECC_PACKED_POLY_CHAIN=1','-DECC_PACKED_UNROLL_INV=1','-DECC_PACKED_PAIR_PRODUCTS=1',
 '-DECC_PACKED_POLY_STATE=1','-DECC_PACKED_DIRECT_REDUCE=1','-DECC_PACKED_GENERATED_PRODUCT=1',
 '-DECC_PACKED_CLMAD=1','-DECC_PACKED_WEIGHTED_PREFIX=2','-DECC_PACKED_COMPACT_STATE=1',
 '-DECC_PACKED_SHARED_SIGMA=1','-DECC_PACKED_STATE_TILE=256','-DECC_PACKED_INLINE=3',
 '-DECC_PACKED_BLOCK_INVERSE=1','-DECC_PACKED_BATCH_SPLIT=2','-DECC_PACKED_PARTIAL_SIGMA=1',
 '-DECC_PACKED_LAST_SLOT_CACHE=1','-DECC_PACKED_TAIL_LAYOUT=1','-DECC_PACKED_SIGMA_ORDER=1',
 '-DECC_PACKED_BYTE_SIGMA=1','-DECC_PACKED_FAST_CONVERT=1','-DECC_PACKED_LOGICAL_PAIR_INVERSE=2',
 '-Xptxas','-v','-lineinfo','-Xcompiler','-O3','-Xcompiler','-fopenmp']
rows=[]
for slots in (4,5,6):
 label=f'xcache-{slots}'
 out=ROOT/'build'/f'g7-{label}'
 cmd=common+[f'-DECC_PACKED_SHARED_X_SLOTS={slots}',str(SRC/'src/main.cu'),'-o',str(out),'-lgomp']
 t=time.monotonic(); p=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=600)
 log=EXP/f'{label}-build.log';log.write_text(p.stdout)
 row={'label':label,'slots':slots,'command':cmd,'returncode':p.returncode,'elapsed_seconds':time.monotonic()-t,
      'log':log.name,'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest()}
 if p.returncode==0:row['binary']=str(out.relative_to(ROOT));row['binary_sha256']=hashlib.sha256(out.read_bytes()).hexdigest()
 rows.append(row);print(json.dumps(row),flush=True)
 if p.returncode:break
(EXP/'builds.json').write_text(json.dumps({'passed':len(rows)==3 and all(r['returncode']==0 for r in rows),'rows':rows},indent=2)+'\n')
raise SystemExit(0 if len(rows)==3 and all(r['returncode']==0 for r in rows) else 1)
