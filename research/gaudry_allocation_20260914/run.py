#!/usr/bin/env python3
"""Matched prime-cubic Gaudry replay; preserve failures and all input/output evidence."""
import argparse, hashlib, json, math, os, random, statistics, subprocess, time
from pathlib import Path
HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]

def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def scrub(v):
    if isinstance(v,dict): return {k:scrub(x) for k,x in v.items() if k not in {'wall_s','wall_ms','linear_algebra_ms'}}
    if isinstance(v,list): return [scrub(x) for x in v]
    return v

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--baseline',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--candidate',type=Path,default=REPO)
    a=parser.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    contract=json.loads((HERE/'contract.json').read_text())
    bins={'baseline':a.baseline.resolve()/'target/release/examples','candidate':a.candidate.resolve()/'target/release/examples'}
    cpu=min(os.sched_getaffinity(0))
    provenance={'cpu_affinity':[cpu],'contract_sha256':digest(HERE/'contract.json'),
        'baseline_source_sha256':digest(a.baseline/'src/cryptanalysis/gaudry_cubic.rs'),
        'candidate_source_sha256':digest(a.candidate/'src/cryptanalysis/gaudry_cubic.rs'),
        'runner_sha256':digest(Path(__file__)),
        'stage_harness_sha256':digest(REPO/'examples/gaudry_allocation_stage.rs'),
        'full_harness_sha256':digest(REPO/'examples/gaudry_cubic_bench.rs'),
        'binaries':{v:{name:digest(d/name) for name in ['gaudry_allocation_stage','gaudry_cubic_bench']} for v,d in bins.items()},
        'rustc':subprocess.check_output(['rustc','--version'],text=True).strip(),
        'host':os.uname().machine,
        'cargo_lock_sha256':digest(a.baseline/'Cargo.lock'),
        'cargo_manifest_sha256':digest(a.baseline/'Cargo.toml')}
    assert (a.baseline/'Cargo.lock').read_bytes()==(a.candidate/'Cargo.lock').read_bytes()
    assert (a.baseline/'Cargo.toml').read_bytes()==(a.candidate/'Cargo.toml').read_bytes()
    (a.output/'cargo-lock.txt').write_bytes((a.baseline/'Cargo.lock').read_bytes())
    (a.output/'cargo-manifest.toml').write_bytes((a.baseline/'Cargo.toml').read_bytes())
    (a.output/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    (a.output/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    schedule=[(kind,p,seed,rep) for kind in ['stage','full_dlp'] for p,seed in contract[kind+'_cases'] for rep in range(contract['repetitions'])]
    rng=random.Random(contract['timing_order_seed']);rng.shuffle(schedule)
    rows=[]
    with (a.output/'raw.jsonl').open('w') as out:
        for kind,p,seed,rep in schedule:
            variants=list(bins);rng.shuffle(variants)
            for variant in variants:
                name=f'{kind}-p{p}-seed{seed}-r{rep}-{variant}'
                result=a.output/(name+'.json')
                command=['taskset','-c',str(cpu)]
                if kind=='stage': command += [str(bins[variant]/'gaudry_allocation_stage'),str(p),str(seed),str(contract['stage_residuals']),str(result.resolve())]
                else: command += [str(bins[variant]/'gaudry_cubic_bench'),'--p',str(p),'--seed',str(seed),'--groebner','--cross-check','--max-residuals',str(contract['max_residuals']),'--json',str(result.resolve())]
                row=dict(kind=kind,p=p,seed=seed,repetition=rep,variant=variant,command=command,result=result.name)
                start=time.perf_counter()
                try:
                    proc=subprocess.run(command,capture_output=True,text=True,timeout=contract['timeout_seconds'])
                    row.update(returncode=proc.returncode,stdout=proc.stdout,stderr=proc.stderr)
                except subprocess.TimeoutExpired as e: row.update(error='timeout',stdout=str(e.stdout),stderr=str(e.stderr))
                row['process_wall_s']=time.perf_counter()-start
                if result.exists(): row['result_sha256']=digest(result)
                rows.append(row);out.write(json.dumps(row)+'\n');out.flush()
    failures=[];pairs=[]
    for row in rows:
        if row.get('returncode')!=0 or 'result_sha256' not in row:
            failures.append(row['result']);continue
        data=json.loads((a.output/row['result']).read_text())
        if row['kind']=='full_dlp':
            g,r=data[0]['gaudry'],data[0]['rho']
            if not (g['correct'] and r['correct'] and g['cross_check_mismatches']==0 and g['cross_checked']==g['residuals']): failures.append(row['result'])
        if row['variant']!='candidate':continue
        base=next(b for b in rows if all(b[k]==row[k] for k in ['kind','p','seed','repetition']) and b['variant']=='baseline')
        if base.get('returncode')!=0 or 'result_sha256' not in base:continue
        before=json.loads((a.output/base['result']).read_text())
        equal=scrub(before)==scrub(data)
        if not equal:failures.append('output/counter mismatch: '+row['result'])
        pairs.append(dict(kind=row['kind'],p=row['p'],seed=row['seed'],repetition=row['repetition'],equal=equal,
                          baseline_s=base['process_wall_s'],candidate_s=row['process_wall_s'],ratio=row['process_wall_s']/base['process_wall_s']))
    summary={'complete':not failures,'runs':len(rows),'failures':failures,'pairs':pairs,'aggregates':[],
             'primary_total_operations':None,'full_cost_floor_ratio':None,'full_cost_rho_ratio':None}
    for kind in ['stage','full_dlp']:
        ps=[p for p in pairs if p['kind']==kind];groups=sorted(set((p['p'],p['seed']) for p in ps))
        # Cluster by input (not by repetition); timings across seed/size cases
        # are paired. Small corpus: intervals describe only this frozen suite.
        logs=[statistics.mean(math.log(p['ratio']) for p in ps if (p['p'],p['seed'])==g) for g in groups]
        boot=sorted(math.exp(statistics.mean(rng.choices(logs,k=len(logs)))) for _ in range(10000)) if logs else []
        summary['aggregates'].append(dict(kind=kind,pairs=len(ps),geometric_time_ratio=math.exp(statistics.mean(logs)) if logs else None,
            paired_cluster_bootstrap_95=[boot[249],boot[9749]] if boot else None))
    (a.output/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary['aggregates'],indent=2));print('failures:',failures)
    return bool(failures)
if __name__=='__main__':raise SystemExit(main())
