#!/usr/bin/env python3
"""Summarize a fully replayed high-K screen against the frozen K3072 row."""
import argparse
import json
import math
import statistics
from pathlib import Path

import run

HERE = Path(__file__).resolve().parent


def summarize(root):
    plan = run.frozen(Path(json.loads((HERE/'plan.json').read_text())['binary_path']))
    receipt = json.loads((root/'plan_receipt.json').read_text())
    if (receipt['plan_sha256'] != run.digest(HERE/'plan.json') or
            receipt['plan'] != plan or receipt['raw_root'] != str(root)):
        raise ValueError('raw plan receipt differs')
    baseline = json.loads(run.BASELINE.read_text())
    reference = baseline['per_target_online_ms']
    if (baseline['status'] != 'PASS_INDEPENDENT_BASE_RELATION_SCALAR_REPLAY' or
            len(reference) != 128 or
            run.digest(run.BASELINE) != plan['sha256']['baseline_k3072']):
        raise ValueError('frozen K3072 reference differs')
    rows = []
    for k in run.ARMS:
        folder = root/f'k{k}'
        if not folder.exists():
            rows.append({'k':k,'status':'NOT_RUN_AFTER_STOP_RULE'})
            continue
        producer = json.loads((folder/'run_receipt.json').read_text())
        if producer['producer_status'] != 'EXIT_0':
            rows.append({'k':k,'status':'CENSORED_PRODUCER',
                         'producer_status':producer['producer_status'],
                         'outer_wall_seconds':producer['outer_wall_seconds']})
            continue
        analysis = json.loads((folder/'analysis.json').read_text())
        if analysis['status'] != 'PASS_INDEPENDENT_BASE_RELATION_SCALAR_REPLAY':
            rows.append({'k':k,'status':'CENSORED_VERIFICATION',
                         'analysis':analysis})
            continue
        if analysis['raw_receipt_sha256'] != run.digest(folder/'run_receipt.json'):
            raise ValueError(f'K{k} raw analysis receipt differs')
        times = analysis['per_target_online_ms']
        if len(times) != 128 or any(not math.isfinite(x) or x <= 0 for x in times):
            raise ValueError(f'K{k} target intervals differ')
        differences = [math.log(b/a) for a,b in zip(times,reference)]
        mean = statistics.mean(differences)
        width = 1.979*statistics.stdev(differences)/math.sqrt(128)
        rows.append({'k':k,'status':analysis['status'],
                     'analysis_sha256':run.digest(folder/'analysis.json'),
                     'actual_usable_points':analysis['actual_usable_points'],
                     'effective_columns':analysis['effective_columns'],
                     'geometric_target_ms':analysis['online_geometric_mean_ms'],
                     'median_target_ms':analysis['online_median_ms'],
                     'mean_probes':analysis['mean_probes'],
                     'cold_process_seconds':analysis['cold_process_ms']/1000,
                     'index_allocated_gib':(
                         analysis['root_state_allocated_bytes']+
                         analysis['root_table_allocated_bytes'])/2**30,
                     'paired_geometric_speedup_vs_k3072':math.exp(mean),
                     'descriptive_95_interval_vs_k3072':[
                         math.exp(mean-width),math.exp(mean+width)],
                     'faster_point_count_vs_k3072':sum(a<b for a,b in zip(times,reference))})
    completed=[row for row in rows if row['status']=='PASS_INDEPENDENT_BASE_RELATION_SCALAR_REPLAY']
    return {'schema':'n83.root11-highk-screen-summary/v1',
            'role':'descriptive_candidate_discovery',
            'hardware_tier':'ordinary_local_exploratory',
            'plan_sha256':run.digest(HERE/'plan.json'),
            'baseline_k3072_sha256':plan['sha256']['baseline_k3072'],
            'baseline_k3072_geometric_target_ms':baseline['online_geometric_mean_ms'],
            'rows':rows,
            'lowest_observed_target_ms_k':min(completed,key=lambda x:x['geometric_target_ms'])['k'] if completed else None,
            'status':'COMPLETE_FOUR_ARMS' if len(completed)==4 else 'PARTIAL_OR_CENSORED'}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('root',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=summarize(args.root.resolve())
    blob=json.dumps(result,sort_keys=True,indent=2)+'\n'
    if args.output.exists() and args.output.read_text()!=blob:
        raise ValueError('existing summary differs')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(blob)
    print(json.dumps({'status':result['status'],
                      'lowest_observed_target_ms_k':result['lowest_observed_target_ms_k']},
                     sort_keys=True))


if __name__=='__main__':
    main()
