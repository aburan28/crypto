#!/usr/bin/env python3
"""Verify complete half-word solves against retained controls and frozen inputs."""
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys
sys.dont_write_bytecode=True
import reference as ref

read=ref.read
sha=ref.sha
dump=ref.dump
LEGACY_HALF=['half16_scalar','half64_scalar','half16_native','half64_native','half_dispatch']
HALF=LEGACY_HALF+['half16_eor3','half64_eor3']
WORD_EOR3=['eor3_word16','eor3_word64']


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


def verify_isolation(root,stage,cpus,cell,mode):
    resource_mode=stage.get('resource_mode',mode)
    assert resource_mode in [mode,'paired']
    name=f'{cell}-{resource_mode}-conditions.jsonl'
    assert stage['conditions_file']==name and sha(root/name)==stage['conditions_sha256']
    assert stage['qualified'] and stage['exit_code']==0 and not stage['timed_out']
    records=[json.loads(line) for line in (root/name).read_text().splitlines()]
    assert len(records)==1
    record=records[0];run=record['run'];pre=record['preflight']
    assert record['schema']=='isolated-bench/1' and record['mode']=='run'
    assert record['reserved_cpus']==run['cpus']==cpus and run['command']==stage['worker_command']
    assert run['exit_status']==0 and run['contended'] is False
    assert run['wall_seconds']>0 and 0<=run['other_cpu_seconds']<=.10*run['wall_seconds']
    assert pre['settle']['seconds']==2 and 0<=pre['settle']['other_cpu_seconds']<=.20
    assert len(pre['conditions']['loadavg'])==3
    for name in ['psi_cpu','psi_memory']:
        assert isinstance(pre['conditions'][name],dict)
        assert 0<=pre['conditions'][name]['some']['avg10']<=5
    for phase in ['before','after']:
        assert len(run[phase]['loadavg'])==3
        assert isinstance(run[phase]['psi_cpu'],dict) and isinstance(run[phase]['psi_memory'],dict)
    for key in ['voluntary_switches','involuntary_switches','minor_faults','major_faults','max_rss_kib']:
        assert type(run[key]) is int and run[key]>=0
    assert record['host']['logical_cpus']>len(cpus) and record['host']['machine'] and record['host']['kernel']
    assert stage['worker_peak_rss_bytes']==run['max_rss_kib']*1024
    command=stage['command']
    for flag,value in [('--settle','2'),('--max-other-cpu','0.10'),('--max-psi','5.0')]:
        assert command[command.index(flag)+1]==value


def verify_aa(root,stage,source,pairs,order_seed):
    cell=f"n{source['n']}-"  # Full identity is recovered from the recorded AA file name.
    name=stage['conditions_file'].removesuffix('-aa-conditions.jsonl')+'-aa.jsonl'
    assert name.startswith(cell)
    path=root/name
    assert sha(path)==stage['stdout_sha256']
    assert sha(path.with_suffix('.stderr'))==stage['stderr_sha256']
    assert path.with_suffix('.stderr').read_bytes()==b''
    rows=[json.loads(line) for line in path.read_text().splitlines()]
    fixture,samples=rows[0],rows[1:]
    assert fixture=={**source,'mode':'aa'}
    reps=int(stage['worker_command'][4]);assert reps in [8,16]
    mode='paired' if stage.get('resource_mode')=='paired' else 'aa'
    assert len(samples)==2*reps
    assert stage['worker_command'][1:]==[str(source['n']),str(source['seed']),source['family'],str(reps),'200000',str(order_seed),mode]
    aa={(s['rep'],s['variant']):s for s in samples};assert len(aa)==len(samples)
    ratios=[]
    for rep in range(reps):
        chunk=samples[2*rep:2*rep+2]
        assert [s['variant'] for s in chunk]==[['aa_a','aa_b'][j] for j in ref.paired_order(2,order_seed,rep)]
        assert [s['order'] for s in chunk]==[0,1]
        for sample in chunk:
            assert sample['rep']==rep and sample['type']=='sample' and sample['verified']
            assert sample['half_work'] is None and sample['projected_work'] is None
            assert 0<sample['solve_ns']+sample['validation_ns']<=sample['total_ns']
            baseline=pairs[rep,'word_dispatch']
            assert (sample['outcome'],sample['model'],sample['reason'],sample['logical'])==(baseline['outcome'],baseline['model'],baseline['reason'],baseline['logical'])
        ratios.append(aa[rep,'aa_a']['total_ns']/aa[rep,'aa_b']['total_ns'])
    return {'paired_ratios':ratios,'symmetric_ratios':[max(x,1/x) for x in ratios],
        'median':statistics.median(ratios),'minimum':min(ratios),'maximum':max(ratios)}


def verify_paired_stream(root,entry,cell):
    receipt=entry['paired'];path=root/(cell+'-paired.jsonl')
    assert sha(path)==receipt['stdout_sha256']
    assert sha(path.with_suffix('.stderr'))==receipt['stderr_sha256']
    assert path.with_suffix('.stderr').read_bytes()==b''
    data=path.read_bytes();parts=[];offset=0
    for mode in ['aa','ab']:
        stage=entry['stages'][mode];origin=stage['derived_from']
        assert stage['resource_mode']=='paired' and origin['file']==path.name
        assert origin['sha256']==receipt['stdout_sha256'] and origin['byte_start']==offset
        part=data[offset:offset+origin['byte_count']]
        expected=root/(cell+('-aa' if mode=='aa' else '')+'.jsonl')
        assert expected.read_bytes()==part and sha(expected)==stage['stdout_sha256']
        assert stage['worker_command']==receipt['worker_command']
        assert stage['conditions_file']==receipt['conditions_file']
        assert stage['conditions_sha256']==receipt['conditions_sha256']
        parts.append(part);offset+=origin['byte_count']
    assert offset==len(data) and b''.join(parts)==data


def validate(root):
    protocol,meta=read(root/'protocol.json'),read(root/'metadata.json')
    assert meta['complete']
    for name,digest in meta['source_hashes'].items():assert sha(root/name)==digest,name
    prior=read(root/'REFERENCE_PROTOCOL.json')
    assert protocol['retained_arms']==protocol['reference_arms']==prior['variants']
    isolated=protocol.get('schema_version',1)>=2
    paired=protocol.get('schema_version',1)>=3
    half=HALF if isolated else LEGACY_HALF
    candidates=half+WORD_EOR3 if isolated else half
    assert protocol['new_candidates']==candidates and protocol['variants']==prior['variants']+candidates
    if isolated:
        assert meta['qualified'] and meta['threads']==1 and meta['cpus']
        assert meta['host']['system']=='Linux' and meta['host']['cpu_info'] and meta['host']['memory_total_kib']>0
        for binary in ['baseline','worker','tests']:assert sha(root/binary)==meta[binary+'_sha256']
    assert protocol['variables']==[12,16,20,24] and protocol['repetitions']==8
    if paired:assert protocol['repetitions_by_n']=={'12':16,'16':8,'20':8,'24':8}
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
    if isolated:assert sha(root/'baseline_worker.rs')==lineage['files']['worker.rs']
    assert sha(root/'reference.py')==lineage['reference_python']['sha256']
    for name,digest in lineage['files'].items():
        data=(root/name).read_bytes()
        if name=='worker.rs':
            suffix=lineage.get('worker_suffix','\ninclude!("halfword.rs");\ninclude!("halfword_tests.rs");\ninclude!("halfword_main.rs");\n').encode()
            assert data.endswith(suffix);data=data[:-len(suffix)]
        assert hashlib.sha256(data).hexdigest()==digest,name
    test=read(root/'test_receipt.json');assert test['exit_code']==0 and not test['timed_out']
    for stream in ['stdout','stderr']:assert hashlib.sha256(test[stream].encode()).hexdigest()==test[stream+'_sha256']
    prior_inputs=read(root/'REFERENCE_FIXTURES.json')['fixtures']
    prior_inputs={(r['n'],r['seed'],r['family']):r for r in prior_inputs};assert len(prior_inputs)==216
    expected={(n,split,seed,f):f'n{n}-{split}-{seed}-{f}' for n in protocol['variables'] for split in protocol['splits'] for seed in protocol[split+'_seeds'] for f in protocol['families']}
    receipts=read(root/'receipts.json');assert len(receipts)==len(expected)==(24 if phase=='discovery' else 240)
    receipts={r['cell']:r for r in receipts};assert set(receipts)==set(expected.values())
    summaries=[];groups=defaultdict(list);noise=defaultdict(list);aa_summaries=[];all_complete=True;count=0
    for (n,split,seed,family),cell in expected.items():
        entry=receipts[cell]
        assert (entry['n'],entry['split'],entry['seed'],entry['family'])==(n,split,seed,family)
        reps=protocol.get('repetitions_by_n',{}).get(str(n),protocol['repetitions'])
        if paired:verify_paired_stream(root,entry,cell)
        if isolated:
            assert set(entry['stages'])=={'aa','ab'}
            for mode in ['aa','ab']:verify_isolation(root,entry['stages'][mode],meta['cpus'],cell,mode)
            receipt={**entry,**entry['stages']['ab'],'command':entry['stages']['ab']['worker_command']}
        else:receipt=entry
        assert receipt['exit_code']==0 and not receipt['timed_out']
        for suffix,stream in [('.jsonl','stdout'),('.stderr','stderr')]:assert sha(root/(cell+suffix))==receipt[stream+'_sha256']
        assert (root/(cell+'.stderr')).read_bytes()==b''
        order_seed=(protocol['order_seed']^(n<<48)^(seed<<8)^protocol['families'].index(family))&((1<<64)-1)
        assert receipt['order_seed']==order_seed
        assert receipt['command'][1:]==[str(n),str(seed),family,str(reps),'200000',str(order_seed)]+(['paired' if paired else 'ab'] if isolated else [])
        rows=[json.loads(line) for line in (root/(cell+'.jsonl')).read_text().splitlines()]
        source,samples=rows[0],rows[1:]
        assert source['type']=='fixture' and (source['n'],source['seed'],source['family'])==(n,seed,family)
        polys,witness=ref.fixture(n,seed,family)
        assert source['polys']==polys and source['planted_witness']==witness
        if split!='holdout':assert prior_inputs[n,seed,family]=={k:source[k] for k in ['n','seed','family','polys','planted_witness']}
        if witness is not None:assert ref.satisfies(polys,witness)
        if isolated:
            assert source['mode']=='ab' and type(source['eor3_available']) is bool
            assert source['architecture']==meta['capabilities']['architecture']
            assert source['eor3_available']==meta['capabilities']['eor3_available']
        names=protocol['variants'];assert len(samples)==reps*len(names)
        pairs={(s['rep'],s['variant']):s for s in samples};assert len(pairs)==len(samples)
        stable={}
        for rep in range(reps):
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
                if arm in half:
                    control=protocol['matched_controls'][arm]
                    check_half(sample,pairs[rep,control],n,arm)
                    assert sample['projected_work'] is None
                else:assert sample['half_work'] is None
                if arm in WORD_EOR3 and isolated:
                    old=pairs[rep,protocol['matched_controls'][arm]]
                    assert (sample['outcome'],sample['model'],sample['logical'])==(old['outcome'],old['model'],old['logical'])
                    assert sample['trace'] is None
                if arm.startswith('projected'):
                    ref.check_work(sample['projected_work'],polys,n,int(arm[-1]),status)
                else:assert sample['projected_work'] is None
                signature={k:v for k,v in sample.items() if k not in ['rep','order','solve_ns','validation_ns','total_ns']}
                if arm in stable:assert stable[arm]==signature
                stable[arm]=signature
                groups[split,n,family,arm].append(sample['total_ns'])
            if all(pairs[rep,a]['outcome']!='UNKNOWN' for a in half):
                counters=['projected_hits','original_checks','rejected_hits','verified_hits']
                for key in counters:assert len({pairs[rep,a]['half_work'][key] for a in half})==1
            assert pairs[rep,'search']['outcome']==source['search_reference']
        if isolated:
            aa=verify_aa(root,entry['stages']['aa'],source,pairs,order_seed)
            noise[split,n,family].extend(aa['symmetric_ratios'])
            aa_summaries.append({'cell':cell,**aa})
        count+=len(samples)
        summaries.append({'cell':cell,'n':n,'split':split,'seed':seed,'family':family,'status':source['search_reference'],
            'complete':all(s['outcome']!='UNKNOWN' for s in samples),
            'medians_ns':{a:statistics.median(pairs[rep,a]['total_ns'] for rep in range(reps)) for a in names},
            'half_work':{a:pairs[0,a]['half_work'] for a in half}})
    gates=[];native_attribution=[]
    splits=['discovery'] if phase=='discovery' else ['regression','holdout']
    for split in splits:
        for n in [16,20,24]:
            for family in protocol['families']:
                for arm in candidates:
                    values=groups[split,n,family,arm]
                    spread=sorted(noise[split,n,family])
                    floor=spread[math.ceil(.975*(len(spread)-1))] if isolated else 1.0
                    hardware_eligible=not isolated or ('eor3' not in arm or meta['capabilities']['eor3_available'])
                    fastest=[min(groups[split,n,family,r][j] for r in protocol['reference_arms']) for j in range(len(values))]
                    ratios=[a/b for a,b in zip(fastest,values)];ci=ref.interval(ratios)
                    matched=[a/b for a,b in zip(groups[split,n,family,protocol['matched_controls'][arm]],values)]
                    gates.append({'candidate':arm,'split':split,'n':n,'family':family,'observations':len(values),
                        'median':statistics.median(ratios),'ci95':ci,'noise_floor':floor,'hardware_eligible':hardware_eligible,'dramatic_group':all_complete and hardware_eligible and ci[0]>max(2,floor),'incremental_group':all_complete and hardware_eligible and ci[0]>max(1,floor),
                        'matched_control':protocol['matched_controls'][arm],'matched_median':statistics.median(matched),'matched_ci95':ref.interval(matched)})
                for size in [16,64]:
                    arm=f'half{size}_native';control=f'half{size}_scalar'
                    ratios=[a/b for a,b in zip(groups[split,n,family,control],groups[split,n,family,arm])]
                    native_attribution.append({'candidate':arm,'split':split,'n':n,'family':family,
                        'scalar_over_native_median':statistics.median(ratios),'ci95':ref.interval(ratios)})
    decisions={arm:{'dramatic_groups':sum(g['dramatic_group'] for g in gates if g['candidate']==arm),
        'incremental_groups':sum(g['incremental_group'] for g in gates if g['candidate']==arm),
        'groups':len(splits)*9,'dramatic_pass':phase=='full' and all(g['dramatic_group'] for g in gates if g['candidate']==arm),
        'incremental_pass':phase=='full' and all(g['incremental_group'] for g in gates if g['candidate']==arm)} for arm in candidates}
    table_split='discovery' if phase=='discovery' else 'holdout'
    result={'schema_version':protocol.get('schema_version',1),'phase':phase,'qualified':isolated,'aa_summaries':aa_summaries,'cells':len(summaries),'observations':count,'all_complete':all_complete,'summaries':summaries,
        'gates':gates,'native_attribution':native_attribution,'decisions':decisions,
        'n24_ms':[{'variant':arm,**{f:statistics.median(groups[table_split,24,f,arm])/1e6 for f in protocol['families']}} for arm in names],
        'campaign_seconds':meta['campaign_seconds'],'peak_worker_rss_bytes':max(stage['worker_peak_rss_bytes'] for r in receipts.values() for stage in r['stages'].values()) if isolated else max(r['peak_rss_bytes'] for r in receipts.values()),
        'production_solver_cost':None,'full_ic_cost':None,'rho_ratio':None,'calibrated_operation_ratio':None}
    if not isolated:
        # Preserve exact replay of the historical schema rather than rewriting
        # archived results when resource qualification was added later.
        result.pop('qualified');result.pop('aa_summaries')
        for gate in result['gates']:
            gate.pop('noise_floor');gate.pop('hardware_eligible')
    return result


def main(root):
    r=validate(root);dump(root/'results.json',r)
    print(json.dumps({k:r[k] for k in ['phase','cells','observations','all_complete','decisions','campaign_seconds','peak_worker_rss_bytes']}))


if __name__=='__main__':main(Path(sys.argv[1]).resolve())
