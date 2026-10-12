#!/usr/bin/env python3
"""Bind completed evidence and produce a compact summary; never edit raw runs."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import statistics

HERE=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    assert not (HERE/'run_02').exists(), 'This report covers only the primary; bind any successor separately.'
    artifacts=[]
    for name,kind in [('resource_probe_01','discovery only'),('resource_probe_02','discovery only'),('resource_probe_03','discovery only'),('run_01','complete fixed-arm comparison')]:
        root=HERE/name
        for file,digest in read(root/'manifest.json')['files'].items():assert sha(root/file)==digest,(name,file)
        assert read(root/'metadata.json')['complete']
        artifacts.append({'path':name,'kind':kind,'manifest_sha256':sha(root/'manifest.json')})
    census=HERE/'support_census_01'
    for file,digest in read(census/'manifest.json')['files'].items():assert sha(census/file)==digest,file
    assert read(census/'results.json')['input_manifest_sha256']==sha(HERE/'resource_probe_02/manifest.json')
    artifacts.append({'path':census.name,'kind':'discovery structure only','manifest_sha256':sha(census/'manifest.json')})
    root=HERE/'run_01';result=read(root/'results.json');metadata=read(root/'metadata.json');protocol=read(root/'protocol.json')
    for name,digest in metadata['source_hashes'].items():assert sha(HERE/name)==digest,name
    gate_counts={}
    for candidate in protocol['new_candidates']:
        values=[d for d in result['comparisons'] if d['candidate']==candidate and d['reference']=='retained_frontier']
        ratios=[d['paired_ratio_median'] for d in values if d['complete']]
        gate_counts[candidate]={
            'dramatic_passes':sum(d['pass'] for d in values),'required':len(values),
            'incremental_passes':sum(d['incremental_pass'] for d in values),
            'paired_median_range':[min(ratios),max(ratios)] if ratios else None,
            'gates':result['gates'][candidate]}
    case_ratios=[];fingerprints=[]
    for n in protocol['variables']:
        cells={}
        for line in (root/f'raw-n{n}.jsonl').read_text().splitlines():
            row=json.loads(line)
            if row['type']=='fixture':
                cells[row['cell']]={'fixture':row,'samples':{}}
                fingerprints.append(sha_bytes(json.dumps([n,row['polys']],separators=(',',':')).encode()))
            else:cells[row['cell']]['samples'][row['rep'],row['variant']]=row
        for cell,values in cells.items():
            fixture=values['fixture'];pairs=values['samples']
            for candidate in protocol['new_candidates']:
                ratios=[]
                for rep in range(protocol['repetitions']):
                    samples=[pairs[rep,a] for a in protocol['reference_arms']+[candidate]]
                    if all(s['verified'] and s['outcome']!='UNKNOWN' for s in samples):
                        ratios.append(min(pairs[rep,a]['total_ns'] for a in protocol['reference_arms'])/pairs[rep,candidate]['total_ns'])
                case_ratios.append({'cell':cell,'n':n,'seed':fixture['seed'],'family':fixture['family'],'split':fixture['split'],
                    'candidate':candidate,'paired_median_ratio':statistics.median(ratios) if len(ratios)==protocol['repetitions'] else None})
    summary={'schema_version':1,'classification':'engineering','scope':protocol['scope'],
        'cells':result['cells'],'distinct_systems':len(set(fingerprints)),'observations':result['observations'],
        'all_complete_verified':result['all_complete_verified'],'gate_counts':gate_counts,
        'outcomes_per_cell':dict(Counter(c['arms']['gray_simd']['outcome'] for c in result['cells_detail'])),
        'fiber_dimension_per_cell':dict(Counter(c['arms']['fiber_simd']['fiber_dimension'] for c in result['cells_detail'])),
        'fiber_work_one_solve_per_fixture':{k:sum(c['arms']['fiber_simd'][k] for c in result['cells_detail']) for k in ['fiber_prefixes','fiber_filter_rounds','fiber_queries','fiber_zero_rejected','fiber_linear_rejected','fiber_extensions_rejected']},
        'peak_worker_rss_bytes':max(c['peak_rss_bytes'] for c in result['cells_detail']),
        'campaign_seconds':metadata['campaign_seconds'],
        'confirmation_eligible':result['all_complete_verified'] and any(v['retained_frontier']=='PASS' for v in result['gates'].values()),
        'confirmation_executed':False,'case_ratios':case_ratios,
        'calibrated_operation_ratio':None,'production_solver_cost':None,'full_ic_cost':None,'rho_ratio':None}
    (HERE/'SUMMARY.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    ledger={'schema_version':1,'classification':'engineering','scope':protocol['scope'],
        'artifacts':artifacts,'source_lineage_sha256':sha(HERE/'SOURCE_LINEAGE.json'),
        'summary_sha256':sha(HERE/'SUMMARY.json'),'reporter_sha256':sha(Path(__file__)),
        'primary_results_sha256':sha(root/'results.json'),'confirmation_eligible':summary['confirmation_eligible'],
        'confirmation_executed':False,'full_ic_cost':None,'production_solver_cost':None,'rho_ratio':None}
    (HERE/'RUN_LEDGER.json').write_text(json.dumps(ledger,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:summary[k] for k in ['cells','observations','all_complete_verified','gate_counts','confirmation_eligible']}))

def sha_bytes(data):return hashlib.sha256(data).hexdigest()
if __name__=='__main__':main()
