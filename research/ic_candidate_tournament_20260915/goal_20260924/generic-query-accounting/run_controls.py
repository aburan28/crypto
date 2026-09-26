import argparse, copy, hashlib, json, os, subprocess, sys
from pathlib import Path
root=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(root/'research/ic_candidate_tournament_20260915'))
from generic_queries import verify_queries
from oracle import verify
from identity import factor_base_inventory
parser=argparse.ArgumentParser(description='Toy query accounting controls, no performance comparison')
parser.add_argument('--worker',type=Path,required=True)
parser.add_argument('--out',type=Path,required=True)
args=parser.parse_args()
worker=args.worker.resolve()
out=args.out.resolve()
out.mkdir(exist_ok=False)
solvers=('pair_table','enumerate','f4','f5','inherited_f4','sat_xor','sat_cnf')
jobs=[]
for a in (0,1):
    for solver in solvers:
        for la in ('dense','sparse'):
            jobs.append((f'a{a}-{solver}-{la}',dict(mode='ic',degree=9,curve_a=a,target_seeds=[2026092554],algorithm_seed=2026092554,factor_base=dict(kind='factor',index=0),config=dict(solver=solver,linear_algebra=la,batch_trials=8,max_trials=256,summands=2)),True))
for solver in solvers:
    job=copy.deepcopy(jobs[0][1]);job['config'].update(solver=solver,max_trials=1,batch_trials=1)
    jobs.append((f'incomplete-{solver}',job,False))
(out/'inputs.json').write_text(json.dumps([dict(name=name,job=job,complete=complete) for name,job,complete in jobs],indent=2)+'\n')
results=[]
for name,job,complete in jobs:
    run=out/name;run.mkdir()
    cp=subprocess.run([str(worker)],input=json.dumps(job),text=True,capture_output=True,timeout=60,env={**os.environ,'RAYON_NUM_THREADS':'1'})
    (run/'stdout.json').write_text(cp.stdout);(run/'stderr.log').write_text(cp.stderr)
    report=json.loads(cp.stdout)
    assert cp.returncode==(0 if complete else 2),(name,cp.returncode,report)
    assert report['status']==('complete' if complete else 'incomplete'),(name,report)
    queries=verify_queries(report,report['fixture'],2)
    assert report['solve_attempts']==report['log_table_report']['solve_attempts']
    if complete:
        proof=verify(report,report['fixture'],expected_mode='ic',summands=2)
        inventory=factor_base_inventory(report,report['fixture'])
    else:
        proof=None
        inventory=factor_base_inventory(report,report['fixture'])
    receipt=dict(test=name,status='PASS',expected_complete=complete,queries=queries,certificate=proof,inventory=inventory)
    (run/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    results.append(receipt)
    print(name,'PASS',flush=True)
summary=dict(schema_version=1,tests=len(results),passed=len(results),scope='toy query/certificate accounting; no performance comparison',promotion_eligible=False,online_ms=None,cold_ms=None,speedup=None,worker_sha256=hashlib.sha256(worker.read_bytes()).hexdigest(),results=results)
(out/'summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
