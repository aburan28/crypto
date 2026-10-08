#!/usr/bin/env python3
"""Verify complete half-word solves against retained controls and frozen inputs."""
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import statistics
import sys
sys.dont_write_bytecode=True
import reference as ref

read=ref.read
sha=ref.sha
dump=ref.dump
HALF=['half16_scalar','half64_scalar','half16_native','half64_native','half_dispatch']


def width(arm,n):
    return 64 if '64' in arm or (arm=='half_dispatch' and n>20) else 16


def check_half(sample,reference,n,arm):
    work=sample['half_work'];status=sample['outcome']
    assert set(work)=={'lane_bits','block_width','points','batches','projected_hits','original_checks','rejected_hits','verified_hits'}
    assert all(type(x) is int and x>=0 for x in work.values())
    assert work['lane_bits']==16 and work['block_width']==width(arm,n)
    assert work['points']==work['block_width']*work['batches']<=1<<n
    assert work['projected_hits']==work['original_checks']==work['rejected_hits']+work['verified_hits']<=work['points']
    assert work['verified_hits']==(1 if status=='SAT' else 0)
    assert work['points']==sample['logical']['enumeration_points']
    assert work['batches']==sample['logical']['enumeration_batches']
    if status=='UNSAT':assert work['points']==1<<n
    if status!='UNKNOWN':
        assert sample['outcome']==reference['outcome'] and sample['model']==reference['model']
        assert sample['logical']==reference['logical']
    assert sample['trace'] is None


def validate(root):
    protocol,meta=read(root/'protocol.json'),read(root/'metadata.json')
    assert meta['complete']
    for name,digest in meta['source_hashes'].items():assert sha(root/name)==digest,name
    prior=read(root/'REFERENCE_PROTOCOL.json')
    assert protocol['retained_arms']==protocol['reference_arms']==prior['variants']
    assert protocol['new_candidates']==HALF and protocol['variants']==prior['variants']+HALF
    assert protocol['variables']==[12,16,20,24] and protocol['repetitions']==8
    assert protocol['families']==['planted','cross_planted','unplanted'] and protocol['discovery_seeds']==[17,937]
    phase=protocol['phase'];assert phase in ['discovery','full']
    if phase=='discovery':assert protocol['splits']==['discovery']
    else:
        assert protocol['splits']==['discovery','regression','holdout']
        assert protocol['regression_seeds']==prior['regression_seeds']+prior['holdout_seeds']
        assert protocol['holdout_seeds']==[20261026,3407873]
        seeds=sum((protocol[s+'_seeds'] for s in protocol['splits']),[])
        assert len(seeds)==len(set(seeds))==20
    lineage=read(root/'SOURCE_LINEAGE.json')
    assert sha(root/'reference.py')==lineage['reference_python']['sha256']
    for name,digest in lineage['files'].items():
        data=(root/name).read_bytes()
        if name=='worker.rs':
            suffix=b'\ninclude!("halfword.rs");\ninclude!("halfword_tests.rs");\ninclude!("halfword_main.rs");\n'
            assert data.endswith(suffix);data=data[:-len(suffix)]
        assert hashlib.sha256(data).hexdigest()==digest,name
    test=read(root/'test_receipt.json');assert test['exit_code']==0 and not test['timed_out']
    for stream in ['stdout','stderr']:assert hashlib.sha256(test[stream].encode()).hexdigest()==test[stream+'_sha256']
    prior_inputs=read(root/'REFERENCE_FIXTURES.json')['fixtures']
    prior_inputs={(r['n'],r['seed'],r['family']):r for r in prior_inputs};assert len(prior_inputs)==216
    expected={(n,split,seed,f):f'n{n}-{split}-{seed}-{f}' for n in protocol['variables'] for split in protocol['splits'] for seed in protocol[split+'_seeds'] for f in protocol['families']}
    receipts=read(root/'receipts.json');assert len(receipts)==len(expected)==(24 if phase=='discovery' else 240)
    receipts={r['cell']:r for r in receipts};assert set(receipts)==set(expected.values())
    summaries=[];groups=defaultdict(list);all_complete=True;count=0
    for (n,split,seed,family),cell in expected.items():
        receipt=receipts[cell];assert receipt['exit_code']==0 and not receipt['timed_out']
        assert (receipt['n'],receipt['split'],receipt['seed'],receipt['family'])==(n,split,seed,family)
        for suffix,stream in [('.jsonl','stdout'),('.stderr','stderr')]:assert sha(root/(cell+suffix))==receipt[stream+'_sha256']
        assert (root/(cell+'.stderr')).read_bytes()==b''
        order_seed=(protocol['order_seed']^(n<<48)^(seed<<8)^protocol['families'].index(family))&((1<<64)-1)
        assert receipt['order_seed']==order_seed
        assert receipt['command'][1:]==[str(n),str(seed),family,'8','200000',str(order_seed)]
        rows=[json.loads(line) for line in (root/(cell+'.jsonl')).read_text().splitlines()]
        source,samples=rows[0],rows[1:]
        assert source['type']=='fixture' and (source['n'],source['seed'],source['family'])==(n,seed,family)
        polys,witness=ref.fixture(n,seed,family)
        assert source['polys']==polys and source['planted_witness']==witness
        if split!='holdout':assert prior_inputs[n,seed,family]=={k:source[k] for k in ['n','seed','family','polys','planted_witness']}
        if witness is not None:assert ref.satisfies(polys,witness)
        names=protocol['variants'];assert len(samples)==8*len(names)
        pairs={(s['rep'],s['variant']):s for s in samples};assert len(pairs)==len(samples)
        stable={}
        for rep in range(8):
            chunk=samples[rep*len(names):(rep+1)*len(names)]
            assert [s['variant'] for s in chunk]==[names[i] for i in ref.paired_order(len(names),order_seed,rep)]
            assert [s['order'] for s in chunk]==list(range(len(names)))
            assert all(s['rep']==rep and s['type']=='sample' for s in chunk)
            for sample in chunk:
                arm,status=sample['variant'],sample['outcome']
                assert all(type(sample[k]) is int and sample[k]>=0 for k in ['solve_ns','validation_ns','total_ns'])
                assert 0<sample['solve_ns']+sample['validation_ns']<=sample['total_ns']
                assert status in ['SAT','UNSAT','UNKNOWN']
                if status=='SAT':
                    assert type(sample['model']) is int and 0<=sample['model']<1<<n
                    assert sample['verified'] and sample['reason'] is None and ref.satisfies(polys,sample['model'])
                    assert source['search_reference']!='UNSAT'
                elif status=='UNSAT':
                    assert sample['model'] is None and sample['verified'] and sample['reason'] is None
                    assert source['search_reference']=='UNSAT' and witness is None
                else:
                    assert sample['model'] is None and not sample['verified'] and sample['reason'] is not None;all_complete=False
                if arm in HALF:
                    control=protocol['matched_controls'][arm]
                    check_half(sample,pairs[rep,control],n,arm)
                    assert sample['projected_work'] is None
                else:assert sample['half_work'] is None
                if arm.startswith('projected'):
                    ref.check_work(sample['projected_work'],polys,n,int(arm[-1]),status)
                else:assert sample['projected_work'] is None
                signature={k:v for k,v in sample.items() if k not in ['rep','order','solve_ns','validation_ns','total_ns']}
                if arm in stable:assert stable[arm]==signature
                stable[arm]=signature
                groups[split,n,family,arm].append(sample['total_ns'])
            if all(pairs[rep,a]['outcome']!='UNKNOWN' for a in HALF):
                counters=['projected_hits','original_checks','rejected_hits','verified_hits']
                for key in counters:assert len({pairs[rep,a]['half_work'][key] for a in HALF})==1
            assert pairs[rep,'search']['outcome']==source['search_reference']
        count+=len(samples)
        summaries.append({'cell':cell,'n':n,'split':split,'seed':seed,'family':family,'status':source['search_reference'],
            'complete':all(s['outcome']!='UNKNOWN' for s in samples),
            'medians_ns':{a:statistics.median(pairs[rep,a]['total_ns'] for rep in range(8)) for a in names},
            'half_work':{a:pairs[0,a]['half_work'] for a in HALF}})
    gates=[];native_attribution=[]
    splits=['discovery'] if phase=='discovery' else ['regression','holdout']
    for split in splits:
        for n in [16,20,24]:
            for family in protocol['families']:
                for arm in HALF:
                    values=groups[split,n,family,arm]
                    fastest=[min(groups[split,n,family,r][j] for r in protocol['reference_arms']) for j in range(len(values))]
                    ratios=[a/b for a,b in zip(fastest,values)];ci=ref.interval(ratios)
                    matched=[a/b for a,b in zip(groups[split,n,family,protocol['matched_controls'][arm]],values)]
                    gates.append({'candidate':arm,'split':split,'n':n,'family':family,'observations':len(values),
                        'median':statistics.median(ratios),'ci95':ci,'dramatic_group':all_complete and ci[0]>2,'incremental_group':all_complete and ci[0]>1,
                        'matched_control':protocol['matched_controls'][arm],'matched_median':statistics.median(matched),'matched_ci95':ref.interval(matched)})
                for size in [16,64]:
                    arm=f'half{size}_native';control=f'half{size}_scalar'
                    ratios=[a/b for a,b in zip(groups[split,n,family,control],groups[split,n,family,arm])]
                    native_attribution.append({'candidate':arm,'split':split,'n':n,'family':family,
                        'scalar_over_native_median':statistics.median(ratios),'ci95':ref.interval(ratios)})
    decisions={arm:{'dramatic_groups':sum(g['dramatic_group'] for g in gates if g['candidate']==arm),
        'incremental_groups':sum(g['incremental_group'] for g in gates if g['candidate']==arm),
        'groups':len(splits)*9,'dramatic_pass':phase=='full' and all(g['dramatic_group'] for g in gates if g['candidate']==arm),
        'incremental_pass':phase=='full' and all(g['incremental_group'] for g in gates if g['candidate']==arm)} for arm in HALF}
    table_split='discovery' if phase=='discovery' else 'holdout'
    return {'schema_version':1,'phase':phase,'cells':len(summaries),'observations':count,'all_complete':all_complete,'summaries':summaries,
        'gates':gates,'native_attribution':native_attribution,'decisions':decisions,
        'n24_ms':[{'variant':arm,**{f:statistics.median(groups[table_split,24,f,arm])/1e6 for f in protocol['families']}} for arm in names],
        'campaign_seconds':meta['campaign_seconds'],'peak_worker_rss_bytes':max(r['peak_rss_bytes'] for r in receipts.values()),
        'production_solver_cost':None,'full_ic_cost':None,'rho_ratio':None,'calibrated_operation_ratio':None}


def main(root):
    r=validate(root);dump(root/'results.json',r)
    print(json.dumps({k:r[k] for k in ['phase','cells','observations','all_complete','decisions','campaign_seconds','peak_worker_rss_bytes']}))


if __name__=='__main__':main(Path(sys.argv[1]).resolve())
