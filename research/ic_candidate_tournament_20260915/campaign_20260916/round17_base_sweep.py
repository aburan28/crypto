#!/usr/bin/env python3
"""Sweep the factor-base `points` parameter over frozen round-0016 fixtures of the
losing cells with the frozen round-0016 worker, recording Ir (callgrind), native
wall, base size, columns, trials and which pipeline ran. The table in
ROUND17-single-target.md is this script's output. Run from the campaign root."""
import json,os,sys,subprocess,statistics,time,re,collections
sys.path.insert(0,'.')
import tournament as T
W='runs/round-0016/worker'
fx=json.load(open('runs/round-0016/fixtures.json'))['confirmation']
cfg=json.load(open('runs/round-0016/candidates.json'))[0]['config']
env=T.child_env()
CELLS=['n19a0','n23a0','n23a1','n29a1','n31a0']
MULTS=[4,6,8,10,12,16,20,24,32]
NFIX=4
def ir_of(job):
    r=subprocess.run(['valgrind','--tool=callgrind','--callgrind-out-file=/dev/null',W],
        input=json.dumps(job),text=True,capture_output=True,env=env,timeout=120)
    m=re.search(r'Collected\s*:\s*(\d+)',r.stderr)
    out=json.loads(r.stdout) if r.stdout.strip().startswith('{') else {}
    return (int(m.group(1)) if m else None), r.returncode, out
def wall_of(job,reps=5):
    ts=[]
    for _ in range(reps):
        t0=time.perf_counter(); r=subprocess.run([W],input=json.dumps(job),text=True,capture_output=True,env=env,timeout=60); ts.append(time.perf_counter()-t0)
    return statistics.median(ts)*1e6
rows=[]
for cell in CELLS:
    cases=[c for c in fx if c['cell']==cell][:NFIX]
    n=cases[0]['job']['degree']
    for c in cases:
        rj=dict(c['job']); rj['mode']='rho'; rj['config']=cfg
        rir,rrc,_=ir_of(rj); rw=wall_of(rj)
        rows.append(dict(cell=cell,case=c['id'],arm='rho',mult=None,points=None,ir=rir,rc=rrc,wall=rw,path='rho',fb=None,cols=None,trials=None))
        for m in MULTS:
            j=dict(c['job']); j['mode']='ic'; j['config']=cfg
            j['factor_base']=dict(c['job']['factor_base'],points=m*n)
            ir,rc,out=ir_of(j); w=wall_of(j) if rc==0 else float('nan')
            path='tiny' if out.get('pair_table') else ('general' if out.get('sparse_report') is not None else '?')
            rows.append(dict(cell=cell,case=c['id'],arm='ic',mult=m,points=m*n,ir=ir,rc=rc,wall=w,path=path,
                             fb=len(out.get('factor_base',[])) if out else None,cols=out.get('columns'),trials=out.get('trials')))
        print(cell,c['id'],'done',flush=True)
json.dump(rows,open('/tmp/claude-0/-home-user/d851d573-3be8-5245-8e4e-603b7615435c/scratchpad/sweep17.json','w'))
# summary: per cell, median over fixtures of ic/rho ratio per multiplier
by=collections.defaultdict(dict)
for r in rows: by[(r['cell'],r['case'])][r['arm'],r['mult']]=r
print(f"\n{'cell':7} {'mult':>4} {'Ir ic/rho':>10} {'wall ic/rho':>12} {'fb pts':>7} {'cols':>5} {'trials':>7} {'path':>8} {'ok':>3}")
for cell in CELLS:
    for m in MULTS:
        irs=[];ws=[];fbs=[];cols=[];tr=[];paths=set();ok=0
        for (cl,case),d in by.items():
            if cl!=cell: continue
            rho=d[('rho',None)]; ic=d.get(('ic',m))
            if not ic or ic['rc']!=0 or not ic['ir']: continue
            ok+=1; irs.append(ic['ir']/rho['ir']); ws.append(ic['wall']/rho['wall']); fbs.append(ic['fb']); cols.append(ic['cols'] or 0); tr.append(ic['trials'] or 0); paths.add(ic['path'])
        if irs: print(f"{cell:7} {m:>4} {statistics.median(irs):10.4f} {statistics.median(ws):12.4f} {statistics.median(fbs):7.0f} {statistics.median(cols):5.0f} {statistics.median(tr):7.0f} {','.join(sorted(paths)):>8} {ok:>3}")
        else: print(f"{cell:7} {m:>4} {'FAIL':>10}")
