#!/usr/bin/env python3
"""Full-DLP confirmation, retaining CPU time to diagnose descheduling outliers."""
import argparse, json, math, os, random, resource, subprocess, time
from pathlib import Path
from run import digest, scrub, interval
HERE=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--build-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    c=json.loads((HERE/'confirmation-contract.json').read_text());cpu=min(os.sched_getaffinity(0))
    provenance={'contract_sha256':digest(HERE/'confirmation-contract.json'),'runner_sha256':digest(Path(__file__)),
        'helper_sha256':digest(HERE/'run.py'),'cpu_affinity':[cpu],
        'binaries':{v:digest(args.build_root/'bin'/v/'gaudry_cubic_bench') for v in c['variants']}}
    (args.output/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    (args.output/'contract.json').write_text(json.dumps(c,indent=2)+'\n')
    rng=random.Random(c['order_seed']);schedule=[(p,s,r) for p,s in c['cases'] for r in range(c['repetitions'])];rng.shuffle(schedule)
    rows=[];data={};failures=[];pairs=[]
    with (args.output/'raw.jsonl').open('w') as out:
        for p,seed,rep in schedule:
            order=list(c['variants']);rng.shuffle(order)
            for v in order:
                name=f'p{p}-s{seed}-r{rep}-{v}.json';result=args.output/name
                cmd=['taskset','-c',str(cpu),str(args.build_root/'bin'/v/'gaudry_cubic_bench'),'--p',str(p),'--seed',str(seed),'--groebner','--cross-check','--max-residuals','20000','--json',str(result.resolve())]
                before=resource.getrusage(resource.RUSAGE_CHILDREN);t=time.perf_counter()
                row=dict(p=p,seed=seed,repetition=rep,variant=v,result=name,command=cmd)
                try:
                    proc=subprocess.run(cmd,capture_output=True,text=True,timeout=c['timeout_seconds'])
                    row.update(returncode=proc.returncode,stdout=proc.stdout,stderr=proc.stderr)
                except subprocess.TimeoutExpired as e:row.update(error='timeout',stdout=str(e.stdout),stderr=str(e.stderr))
                row['wall_s']=time.perf_counter()-t;after=resource.getrusage(resource.RUSAGE_CHILDREN)
                row['cpu_s']=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime
                if result.exists():
                    row['result_sha256']=digest(result);data[p,seed,rep,v]=json.loads(result.read_text())
                rows.append(row);out.write(json.dumps(row)+'\n');out.flush()
    # Get the actual per-call charge for the two fresh curves as well as old
    # ones with an untimed micro run; all polynomial digests must agree.
    charges={}
    for p,seed in c['cases']:
        oracle=None
        for v in c['variants']:
            result=args.output/f'charge-p{p}-s{seed}-{v}.json'
            proc=subprocess.run([str(args.build_root/'bin'/v/'summation_polynomial_bench'),str(p),str(seed),'256','1','1',str(result.resolve())],capture_output=True,text=True,timeout=180)
            if proc.returncode:
                failures.append('charge failure '+str(result));continue
            d=json.loads(result.read_text());actual=d['specialization_fp_muls']//256
            charges[p,seed,v]=15*(actual//11-4) if v=='legacy' else actual
            if oracle is None:oracle=d['polynomial_digest']
            if oracle!=d['polynomial_digest']:failures.append('polynomial mismatch '+str(result))
    for row in rows:
        p,seed,rep,v=[row[k] for k in ['p','seed','repetition','variant']]
        if row.get('returncode')!=0 or (p,seed,rep,v) not in data:
            failures.append(row['result']);continue
        try:
            d=data[p,seed,rep,v];b=data[p,seed,rep,'reference'];g=d[0]['gaudry'];bg=b[0]['gaudry']
            assert scrub(d,{'fp_muls','oracle_fp_muls','total_ops','s','fp_muls_per_pair_test'})==scrub(b,{'fp_muls','oracle_fp_muls','total_ops','s','fp_muls_per_pair_test'})
            assert g['correct'] and d[0]['rho']['correct'] and g['cross_check_mismatches']==0 and g['cross_checked']==g['residuals']
            delta=g['solve_stats']['solves']*(charges[p,seed,v]-charges[p,seed,'reference'])
            assert g['solve_stats']['fp_muls']-bg['solve_stats']['fp_muls']==delta
            assert g['oracle_fp_muls']-bg['oracle_fp_muls']==delta
            assert math.isclose(g['total_ops']-bg['total_ops'],delta/g['fp_muls_per_add'],abs_tol=1e-8)
            base=next(r for r in rows if r['p']==p and r['seed']==seed and r['repetition']==rep and r['variant']=='reference')
            for metric in ['wall_s','cpu_s']:
                pairs.append(dict(p=p,seed=seed,repetition=rep,variant=v,metric=metric,ratio=row[metric]/base[metric],baseline_s=base[metric],candidate_s=row[metric]))
        except (AssertionError,KeyError):failures.append('comparison '+row['result'])
    summary={'complete':not failures,'failures':failures,'runs':len(rows),'aggregates':[],'pairs':pairs}
    if not failures:
        for v in c['variants']:
            for metric in ['wall_s','cpu_s']:
                ps=[p for p in pairs if p['variant']==v and p['metric']==metric];ratio,ci=interval(ps,rng)
                summary['aggregates'].append(dict(variant=v,metric=metric,ratio=ratio,paired_95=ci,sum_s=sum(p['candidate_s'] for p in ps)))
    (args.output/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary['aggregates'],indent=2));print('failures:',failures)
    return bool(failures)
if __name__=='__main__':raise SystemExit(main())
