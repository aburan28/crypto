#!/usr/bin/env python3
"""Summarize the exact K4096 five-byte/four-byte ID30 comparison."""
import argparse
import json
import math
import statistics
from pathlib import Path

import run

HERE = Path(__file__).resolve().parent
PASS = 'PASS_INDEPENDENT_BASE_RELATION_SCALAR_REPLAY'


def summarize(root):
    plan = run.frozen(Path(json.loads((HERE/'plan.json').read_text())['binary_path']))
    receipt = json.loads((root/'plan_receipt.json').read_text())
    if (receipt['plan_sha256'] != run.digest(HERE/'plan.json') or
            receipt['plan'] != plan or receipt['raw_root'] != str(root)):
        raise ValueError('raw plan receipt differs')
    rows=[]
    analyses={}
    for arm in run.ARMS:
        folder=root/arm
        if not folder.exists():
            rows.append({'arm':arm,'status':'NOT_RUN_AFTER_STOP_RULE'})
            continue
        producer=json.loads((folder/'run_receipt.json').read_text())
        if producer['producer_status'] != 'EXIT_0':
            rows.append({'arm':arm,'status':'CENSORED_PRODUCER',
                         'producer_status':producer['producer_status'],
                         'outer_wall_seconds':producer['outer_wall_seconds']})
            continue
        analysis=json.loads((folder/'analysis.json').read_text())
        if analysis['status'] != PASS:
            rows.append({'arm':arm,'status':'CENSORED_VERIFICATION',
                         'analysis':analysis})
            continue
        if (analysis['raw_receipt_sha256'] != run.digest(folder/'run_receipt.json') or
                analysis['arm'] != arm or analysis['k'] != run.K or
                len(analysis['per_target_online_ms']) != 128):
            raise ValueError(f'{arm}: analysis receipt or coverage differs')
        analyses[arm]=analysis
        rows.append({'arm':arm,'status':PASS,
                     'analysis_sha256':run.digest(folder/'analysis.json'),
                     'actual_usable_points':analysis['actual_usable_points'],
                     'effective_columns':analysis['effective_columns'],
                     'geometric_target_ms':analysis['online_geometric_mean_ms'],
                     'median_target_ms':analysis['online_median_ms'],
                     'mean_probes':analysis['mean_probes'],
                     'cold_process_seconds':analysis['cold_process_ms']/1000,
                     'index_allocated_gib':(
                         analysis['root_state_allocated_bytes']+
                         analysis['root_table_allocated_bytes'])/2**30})
    comparison=None
    if len(analyses)==2:
        if analyses['indexed5_id30']['base_hash'] != analyses['indexed4_id30']['base_hash']:
            raise ValueError('indexed layouts produced different factor bases')
        five=analyses['indexed5_id30']['per_target_online_ms']
        four=analyses['indexed4_id30']['per_target_online_ms']
        if any(not math.isfinite(x) or x <= 0 for x in five+four):
            raise ValueError('invalid paired interval')
        delta=[math.log(a/b) for a,b in zip(five,four)]
        mean=statistics.mean(delta)
        width=1.979*statistics.stdev(delta)/math.sqrt(128)
        comparison={'ratio_definition':'five_byte_online_ms_over_four_byte_online_ms',
                    'paired_geometric_ratio':math.exp(mean),
                    'descriptive_95_interval':[math.exp(mean-width),
                                               math.exp(mean+width)],
                    'four_byte_faster_points':sum(a>b for a,b in zip(five,four)),
                    'paired_points':128,
                    'cold_seconds_five_minus_four':(
                        analyses['indexed5_id30']['cold_process_ms']-
                        analyses['indexed4_id30']['cold_process_ms'])/1000}
    return {'schema':'n83.indexed4-id30-screen-summary/v1',
            'role':'ordinary_local_candidate_diagnostic',
            'hardware_tier':'ordinary_local_exploratory',
            'plan_sha256':run.digest(HERE/'plan.json'),
            'rows':rows,'comparison':comparison,
            'status':'COMPLETE_BOTH_ARMS' if comparison else 'PARTIAL_OR_CENSORED'}


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
                      'paired_geometric_ratio':None if result['comparison'] is None
                      else result['comparison']['paired_geometric_ratio']},sort_keys=True))


if __name__=='__main__':main()
