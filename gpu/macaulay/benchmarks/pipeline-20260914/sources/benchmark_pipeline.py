#!/usr/bin/env python3
"""Replay a frozen matched transport benchmark and retain failures, never overwrite."""
import argparse
import hashlib
import json
import platform
import random
import statistics
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--contract', type=Path, default=ROOT/'benchmarks/pipeline-20260914/contract.json')
    args = parser.parse_args()
    contract = json.loads(args.contract.read_text())
    args.output.mkdir(parents=True, exist_ok=False)
    executable = ROOT/'pipeline'
    files = ['pipeline.cu', 'pipeline.cuh', 'kernels.cuh', 'macaulay.cuh', 'fpmod.cuh', 'params.h', 'pipeline', 'benchmark_pipeline.py']
    provenance = {
        'contract_sha256': sha(args.contract),
        'sha256': {name: sha(ROOT/name) for name in files},
        'git_head': subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'platform': platform.platform(),
        'gpu': subprocess.check_output(['nvidia-smi','--query-gpu=name,driver_version','--format=csv,noheader'],text=True).strip(),
    }
    (args.output/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    (args.output/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    pairs = [(split,seed,rep) for split in ['training','holdout']
             for seed in contract[split+'_seeds'] for rep in range(contract['repetitions'])]
    rng = random.Random(contract['run_order_seed']); rng.shuffle(pairs)
    records = []
    with (args.output/'raw.jsonl').open('w') as out:
        for split,seed,rep in pairs:
            variants=list(contract['variants'].items()); rng.shuffle(variants)
            for variant,slots in variants:
                cmd=[str(executable),str(contract['jobs']),str(contract['batch']),str(slots),str(seed),str(contract['min_gpu_batch'])]
                r={'split':split,'seed':seed,'repetition':rep,'variant':variant,'command':cmd}
                start=time.monotonic()
                try:
                    p=subprocess.run(cmd,text=True,capture_output=True,timeout=contract['timeout_seconds'])
                    r.update(returncode=p.returncode,stdout=p.stdout,stderr=p.stderr)
                    if p.returncode == 0:
                        m=json.loads(p.stdout)
                        assert m['verified']==m['jobs']==contract['jobs']
                        assert m['p']==contract['prime'] and [m['rows'],m['cols']]==contract['shape']
                        assert m['seed']==seed and m['slots']==slots and m['batch']==contract['batch']
                        assert m['cutoff']==contract['min_gpu_batch']
                        assert m['cpu_jobs']+m['gpu_jobs']==m['jobs']
                        r['metrics']=m
                except (subprocess.TimeoutExpired,ValueError,AssertionError) as e:
                    r['error']=str(e) or type(e).__name__
                r['process_wall_s']=time.monotonic()-start
                out.write(json.dumps(r)+'\n'); out.flush(); records.append(r)
    failures=[r for r in records if r.get('returncode')!=0 or 'metrics' not in r or 'error' in r]
    summary={'classification':contract['classification'],'complete':not failures,'runs':len(records),
             'failed_runs':len(failures),'primary_total_common_operations':None,'rho_ratio':None,
             'floor_ratio':None,'full_dlp_speedup':None,'rows':[]}
    if not failures:
        for split in ['training','holdout']:
            selected=[r for r in records if r['split']==split]
            for variant in contract['variants']:
                times=[r['metrics']['total_s'] for r in selected if r['variant']==variant]
                ratios=[]
                for r in selected:
                    if r['variant']!=variant: continue
                    base=next(b for b in selected if b['variant']=='serialized' and b['seed']==r['seed'] and b['repetition']==r['repetition'])
                    assert all(r['metrics'][k]==base['metrics'][k] for k in ['jobs','verified','cpu_jobs','gpu_jobs','transfer_bytes'])
                    ratios.append(r['metrics']['total_s']/base['metrics']['total_s'])
                summary['rows'].append({'split':split,'variant':variant,'runs':len(times),
                    'median_total_s':statistics.median(times),'median_paired_time_ratio':statistics.median(ratios),
                    'verified_matrices':sum(r['metrics']['verified'] for r in selected if r['variant']==variant)})
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    return bool(failures)

if __name__=='__main__':
    raise SystemExit(main())
