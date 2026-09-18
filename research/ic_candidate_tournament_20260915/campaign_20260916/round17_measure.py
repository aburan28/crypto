#!/usr/bin/env python3
"""Byte-verify the round-0017 challenger against the oracle on every confirmation
fixture of round 0016, then measure Ir (callgrind) and native wall per cell for
challenger, scan_io and rho.  Usage: measure17.py CHALLENGER_BIN [reps]"""
import json,os,sys,subprocess,statistics,time,re,collections,tempfile
_OUT=os.path.join(tempfile.gettempdir(),'round17_measure.json')  # raw rows; the table below is the record,random
sys.path.insert(0,'.')
import tournament as T
from oracle import verify
CH=sys.argv[1]; REPS=int(sys.argv[2]) if len(sys.argv)>2 else 5
S='runs/round-0016/source_candidates/scan_io/worker'; W='runs/round-0016/worker'
fx=json.load(open('runs/round-0016/fixtures.json'))['confirmation']
cfg=json.load(open('runs/round-0016/candidates.json'))[0]['config']; env=T.child_env()
def run(b,job): return subprocess.run([b],input=json.dumps(job),text=True,capture_output=True,env=env,timeout=60)
def ir(b,job):
    r=subprocess.run(['valgrind','--tool=callgrind','--callgrind-out-file=/dev/null',b],input=json.dumps(job),text=True,capture_output=True,env=env,timeout=180)
    return int(re.search(r'Collected\s*:\s*(\d+)',r.stderr).group(1))
# 1. byte-verify on ALL 96 fixtures: oracle accepts, same base hash and same solutions as scan_io
bad=0; bytes_ch=[]; bytes_s=[]
for c in fx:
    j=dict(c['job'],mode='ic',config=cfg)
    a=run(CH,j); b=run(S,j)
    assert a.returncode==0 and b.returncode==0,(c['id'],a.returncode,b.returncode,a.stderr[-200:])
    ra=json.loads(a.stdout); rb=json.loads(b.stdout)
    pa=verify(ra,c['fixture'],expected_mode='ic',summands=3); pb=verify(rb,c['fixture'],expected_mode='ic',summands=3)
    ok = pa['factor_base_sha256']==pb['factor_base_sha256'] and pa['solutions']==pb['solutions'] and ra['relations']==rb['relations'] and ra['column_logs']==rb['column_logs'] and 'factor_base' not in ra and 'factor_base_orbits' in ra
    bad+= not ok; bytes_ch.append(len(a.stdout)); bytes_s.append(len(b.stdout))
print(f"verify: {len(fx)-bad}/{len(fx)} fixtures identical certificate (base hash, solutions, relations, column logs); output bytes median {statistics.median(bytes_ch):.0f} vs scan_io {statistics.median(bytes_s):.0f}")
assert bad==0
# 2. measurement, 4 fixtures per cell, interleaved arms
cells=sorted({c['cell'] for c in fx},key=lambda s:(int(s[1:s.index('a')]),s))
rows=collections.defaultdict(list)
for cell in cells:
    for c in [x for x in fx if x['cell']==cell][:4]:
        jobs={'ch':(CH,dict(c['job'],mode='ic',config=cfg)),'scan_io':(S,dict(c['job'],mode='ic',config=cfg)),'rho':(W,dict(c['job'],mode='rho',config=cfg))}
        irs={k:ir(b,j) for k,(b,j) in jobs.items()}
        walls={k:[] for k in jobs}
        for _ in range(REPS):
            order=list(jobs); random.shuffle(order)
            for k in order:
                b,j=jobs[k]; payload=json.dumps(j).encode()
                t0=time.perf_counter(); p=subprocess.run([b],input=payload,capture_output=True,env=env); walls[k].append(time.perf_counter()-t0)
        rows[cell].append((irs,{k:statistics.median(v) for k,v in walls.items()}))
print(f"\n{'cell':7} {'ch/rho Ir':>10} {'sio/rho Ir':>11} {'ch/sio Ir':>10} | {'ch/rho wall':>12} {'sio/rho wall':>13} {'ch/sio wall':>12}")
gm=lambda xs: statistics.geometric_mean(xs)
for cell in cells:
    R=rows[cell]
    f=lambda a,b,key: gm([r[key][a]/r[key][b] for r in R])
    print(f"{cell:7} {f('ch','rho',0):10.4f} {f('scan_io','rho',0):11.4f} {f('ch','scan_io',0):10.4f} | {f('ch','rho',1):12.4f} {f('scan_io','rho',1):13.4f} {f('ch','scan_io',1):12.4f}")
json.dump({c:[(i,w) for i,w in rows[c]] for c in cells},open(_OUT,'w'))
print('raw rows written to',_OUT)
