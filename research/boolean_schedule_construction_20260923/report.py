#!/usr/bin/env python3
"""Bind the preserved execution, corrected analysis and structural census."""
from collections import Counter
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    artifacts=[]
    for name,kind in [('run_01','original execution with analysis failure'),('analysis_01','additive corrected analysis'),('quotient_census_01','discovery algebra only')]:
        root=HERE/name
        for file,digest in read(root/'manifest.json')['files'].items():assert sha(root/file)==digest,(name,file)
        artifacts.append({'path':name,'kind':kind,'manifest_sha256':sha(root/'manifest.json')})
    meta=read(HERE/'run_01/metadata.json');analysis_meta=read(HERE/'analysis_01/metadata.json');result=read(HERE/'analysis_01/results.json')
    assert meta['complete'] and analysis_meta['complete']
    assert analysis_meta['input_manifest_sha256']==sha(HERE/'run_01/manifest.json')
    assert read(HERE/'quotient_census_01/results.json')['input_manifest_sha256']==sha(HERE/'run_01/manifest.json')
    for name,digest in meta['source_hashes'].items():
        if name.endswith('.rs'):assert sha(HERE/name)==digest,name
    summary={'schema_version':1,'classification':'accounting','protocol_classification':read(HERE/'protocol.json')['classification'],
        'scope':result['scope'],'cells':result['cells'],'observations':result['observations'],'all_complete_verified':result['all_complete_verified'],
        'decisions':{kind:dict(Counter(b['decision'] for b in result['bounds'] if b['n']>=16 and b['kind']==kind)) for kind in ['16','64','dispatch']},
        'campaign_seconds':meta['campaign_seconds'],'peak_worker_rss_bytes':max(c['peak_rss_bytes'] for c in result['cells_detail']),
        'constructor_optimization_implemented':False,'performance_promotion':None,'full_ic_cost':None,'production_solver_cost':None,'calibrated_operation_ratio':None,'rho_ratio':None}
    table=[]
    for arm in read(HERE/'protocol.json')['variants']:
        row={'variant':arm,'families':{}}
        for family in ['planted','cross_planted','unplanted']:
            c=next(c for c in result['costs'] if c['n']==24 and c['family']==family and c['variant']==arm)
            row['families'][family]=c['completion_ns']
        reference=next(c for c in result['costs'] if c['n']==24 and c['family']=='planted' and c['variant']=='word_dispatch')['completion_ns']
        row['dispatcher_over_arm_planted']=reference/row['families']['planted'] if reference is not None and row['families']['planted'] is not None else None
        table.append(row)
    summary['n24_cost_table']=table
    (HERE/'SUMMARY.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    ledger={'schema_version':1,'classification':'accounting','artifacts':artifacts,'summary_sha256':sha(HERE/'SUMMARY.json'),
        'reporter_sha256':sha(Path(__file__)),'source_lineage_sha256':sha(HERE/'SOURCE_LINEAGE.json'),
        'corrected_analyzer_sha256':sha(HERE/'corrected_analysis.py'),'replay_controller_sha256':sha(HERE/'run.py'),
        'timed_rust_sources_unchanged':True,'constructor_optimization_implemented':False,'performance_promotion':None,
        'full_ic_cost':None,'production_solver_cost':None,'rho_ratio':None}
    (HERE/'RUN_LEDGER.json').write_text(json.dumps(ledger,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:summary[k] for k in ['cells','observations','all_complete_verified','decisions','campaign_seconds','peak_worker_rss_bytes']}))
if __name__=='__main__':main()
