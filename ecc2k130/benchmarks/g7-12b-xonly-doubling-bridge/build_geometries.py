#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,subprocess,time
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent;SRC=ROOT/'build/g7-xonly-doubling-src';NVCC=ROOT/'build/assembler134-toolchain/bin/nvcc'
FLAGS={'SINGLE_PRODUCT':0,'CACHE_DENOM':1,'BY_VALUE':1,'PERM_SIGMA':3,'POLY_CHAIN':1,'UNROLL_INV':1,'PAIR_PRODUCTS':1,'POLY_STATE':1,'DIRECT_REDUCE':1,'GENERATED_PRODUCT':1,'CLMAD':1,'WEIGHTED_PREFIX':2,'COMPACT_STATE':1,'SHARED_SIGMA':1,'STATE_TILE':256,'INLINE':3,'BLOCK_INVERSE':1,'BATCH_SPLIT':2,'PARTIAL_SIGMA':1,'LAST_SLOT_CACHE':1,'TAIL_LAYOUT':1,'SIGMA_ORDER':1,'BYTE_SIGMA':1,'FAST_CONVERT':1,'LOGICAL_PAIR_INVERSE':0,'SHARED_X_SLOTS':4,'XONLY_23':1,'XONLY_DOUBLE_ONLY':1}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
for batch in (24,32):
 binary=ROOT/f'build/g7-xonly-doubling-b{batch}';log=OUT/f'b{batch}-build.log'
 cmd=[str(NVCC),'-O3','-std=c++17','-gencode','arch=compute_120,code=sm_120',f'-DECC_BATCH={batch}','-DECC_THREADS=256','-DECC_MINBLOCKS=3',*[f'-DECC_PACKED_{k}={v}' for k,v in FLAGS.items()],'-Xptxas','-v','-lineinfo','-Xcompiler','-O3','-Xcompiler','-fopenmp',str(SRC/'src/main.cu'),'-o',str(binary),'-lgomp']
 t=time.monotonic()
 with log.open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=900)
 row={'batch':batch,'command':cmd,'returncode':r.returncode,'elapsed_seconds':time.monotonic()-t,'binary_sha256':sha(binary) if binary.exists() else None,'log_sha256':sha(log)}
 (OUT/f'b{batch}-build.json').write_text(json.dumps(row,indent=2)+'\n');print(json.dumps(row),flush=True)
 if r.returncode:raise SystemExit(r.returncode)
