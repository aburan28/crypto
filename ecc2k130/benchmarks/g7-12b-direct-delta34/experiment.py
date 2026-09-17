from pathlib import Path
import hashlib,json,os,signal,subprocess,time
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
SRC=ROOT/'build/g7-seven-node-src';PRIOR=ROOT/'benchmarks/g7-2xlarge-assembler-confirmation'
OLD=ROOT/'build/cuda133/nvidia/cu13';MIXED=ROOT/'build/assembler134-toolchain';NEW=MIXED
REFERENCE_SHA='8c4ed76a152cc135280a68fd3559c3bf71bef726ab73eae0e04014c2e27c6d02'
FLAGS=dict(SINGLE_PRODUCT=0,CACHE_DENOM=1,BY_VALUE=1,PERM_SIGMA=3,POLY_CHAIN=1,UNROLL_INV=1,PAIR_PRODUCTS=1,POLY_STATE=1,DIRECT_REDUCE=1,GENERATED_PRODUCT=1,CLMAD=1,WEIGHTED_PREFIX=2,COMPACT_STATE=1,SHARED_SIGMA=1,STATE_TILE=256,INLINE=3,BLOCK_INVERSE=1,BATCH_SPLIT=2,PARTIAL_SIGMA=1,LAST_SLOT_CACHE=1,TAIL_LAYOUT=1,SIGMA_ORDER=1,BYTE_SIGMA=1,SHARED_X_SLOTS=4,FAST_CONVERT=1)
CONFIGS={'seven-node-control':dict(MINBLOCKS=3,LOGICAL_PAIR_INVERSE=0),'seven-node-three':dict(MINBLOCKS=3,LOGICAL_PAIR_INVERSE=2),'seven-node-four':dict(MINBLOCKS=4,LOGICAL_PAIR_INVERSE=2)}
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def put(n,x):
 with (OUT/n).open('x') as f:json.dump(x,f,indent=2);f.write('\n')
def command(label,source,binary):
 return [str(NEW/'bin/nvcc'),'-O3','-std=c++17','-gencode','arch=compute_120,code=sm_120','-DECC_BATCH=16','-DECC_THREADS=256',f'-DECC_MINBLOCKS={CONFIGS[label]["MINBLOCKS"]}',*[f'-DECC_PACKED_{k}={v}' for k,v in FLAGS.items()],'-DECC_PACKED_LOGICAL_PAIR_INVERSE='+str(CONFIGS[label]['LOGICAL_PAIR_INVERSE']),'-Xptxas','-v','-lineinfo','-Xcompiler','-O3','-Xcompiler','-fopenmp',str(SRC/source),'-o',str(binary),'-lgomp']
def finite(command,log,timeout=600,cwd=Path('/tmp'),env=None):
 start=time.monotonic();timed_out=False;assert not log.exists()
 with log.open('x') as f:
  p=subprocess.Popen(command,cwd=cwd,env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
  try:p.wait(timeout=timeout)
  except subprocess.TimeoutExpired:
   timed_out=True;os.killpg(p.pid,signal.SIGTERM)
   try:p.wait(timeout=5)
   except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
 return dict(command=command,returncode=p.returncode,timed_out=timed_out,timeout_seconds=timeout,elapsed_seconds=time.monotonic()-start,log=str(log.relative_to(OUT)),log_sha256=sha(log))
def build(task):
 label,suffix,source=task;name=label+suffix;binary=ROOT/'build'/('g7-'+name);assert not binary.exists()
 p=finite(command(label,source,binary),OUT/(name+'-build.log'));row=dict(p,label=name,configuration=label,binary=str(binary.relative_to(ROOT)),binary_sha256=sha(binary) if binary.exists() else None)
 put(name+'-build.json',row);print(name,'exit',p['returncode'],flush=True);return row
