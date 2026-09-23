#!/usr/bin/env python3
"""Verify complete-solve evidence without promoting censored outcomes."""
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import random
import statistics
import sys

def read(path):return json.loads(path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def bootstrap(values):
    rng=random.Random(20260923)
    medians=sorted(statistics.median(rng.choices(values,k=len(values))) for _ in range(4000))
    return [medians[100],medians[3899]]
def satisfies(polys,model):return all(sum((model&m)==m for m in p)%2==0 for p in polys)

def main(root):
    protocol,meta=read(root/'protocol.json'),read(root/'metadata.json');assert meta['complete']
    for name,digest in meta['source_hashes'].items():assert sha(root/name)==digest,name
    test=read(root/'test_receipt.json');assert test['exit_code']==0 and not test['timed_out']
    for stream in ['stdout','stderr']:assert hashlib.sha256(test[stream].encode()).hexdigest()==test[stream+'_sha256']
    expected={f'n{n}-{split}-{seed}-{family}' for n in protocol['variables'] for split in ['discovery','holdout'] for seed in protocol[split+'_seeds'] for family in protocol['families']}
    receipts=read(root/'receipts.json');assert len(receipts)==len(expected)
    receipts={r['cell']:r for r in receipts};assert set(receipts)==expected
    raw=defaultdict(list)
    for n in protocol['variables']:
        for line in (root/f'raw-n{n}.jsonl').read_text().splitlines():
            row=json.loads(line);assert json.loads(row['raw_line'])=={k:v for k,v in row.items() if k not in ['cell','split','raw_line']}
            raw[row['cell']].append(row)
    assert set(raw)==expected
    summaries=[];groups=defaultdict(list);all_complete=True
    logical=['outcome','model','reason','nodes','kernel_calls','decisions','forced','specialized_terms','source_rows','source_columns','max_depth','trace']
    for cell,records in sorted(raw.items()):
        receipt=receipts[cell];assert receipt['exit_code']==0 and not receipt['timed_out']
        assert hashlib.sha256(''.join(r['raw_line'] for r in records).encode()).hexdigest()==receipt['stdout_sha256']
        assert receipt['stderr_sha256']==hashlib.sha256(b'').hexdigest()
        fixture=records[0];assert fixture['type']=='fixture'
        for key in ['n','seed','split','family']:assert fixture[key]==receipt[key]
        n=fixture['n'];polys=fixture['polys'];assert len(polys)==n+2
        for p in polys:
            assert p==sorted(set(p),key=lambda m:(-m.bit_count(),m))
            assert len(p)<=32 and all(0<=m<1<<n and m.bit_count()<=2 for m in p)
        if fixture['planted_witness'] is not None:assert satisfies(polys,fixture['planted_witness'])
        samples=records[1:];names=protocol['variants'];pairs={(s['rep'],s['variant']):s for s in samples}
        assert len(samples)==len(names)*protocol['repetitions'] and len(pairs)==len(samples)
        assert set(pairs)=={(rep,arm) for rep in range(protocol['repetitions']) for arm in names}
        for rep in range(protocol['repetitions']):
            ordered=[s for s in samples if s['rep']==rep]
            assert [s['variant'] for s in ordered]==[names[(rep+i)%len(names)] for i in range(len(names))]
            assert [s['order'] for s in ordered]==list(range(len(names)))
        kernel_signature=None
        for sample in samples:
            arm=sample['variant'];status=sample['outcome'];model=sample['model']
            assert 0<=sample['kernel_ns']<=sample['solve_ns']<=sample['total_ns']
            assert sample['solve_ns']+sample['validation_ns']<=sample['total_ns']
            assert sample['nodes']<=protocol['limits']['nodes'] and sample['max_depth']<=n
            assert sum(sample['kernel_calls_by_active'])==sample['kernel_calls']
            assert len(sample['kernel_calls_by_active'])==n+1
            assert sample['flat_calls']<=sample['kernel_calls']
            if arm=='search':assert sample['kernel_calls']==sample['kernel_ns']==0
            if arm=='flat':assert sample['flat_calls']==sample['kernel_calls']
            if arm=='bucket':assert sample['flat_calls']==0
            if status=='SAT':
                assert isinstance(model,int) and 0<=model<1<<n and satisfies(polys,model)
                assert sample['reason'] is None and sample['verified']
            elif status=='UNSAT':
                assert model is None and sample['reason'] is None and fixture['planted_witness'] is None
                assert sample['verified']==(fixture['search_reference']=='UNSAT')
            else:
                assert status=='UNKNOWN' and model is None and not sample['verified']
                assert sample['reason'] in ['NODE_CAP','ROW_CAP','COLUMN_CAP']
            if status=='UNKNOWN' or not sample['verified']:all_complete=False
            if arm!='search':
                signature=tuple(sample[k] for k in logical)
                if kernel_signature is None:kernel_signature=signature
                else:assert signature==kernel_signature,(cell,arm)
            groups[fixture['split'],n,fixture['family'],arm].append(sample)
        arms={}
        for arm in names:
            values=[s for s in samples if s['variant']==arm]
            assert len({tuple(s[k] for k in logical) for s in values})==1
            completed=all(s['outcome']!='UNKNOWN' and s['verified'] for s in values)
            arms[arm]={'outcome':values[0]['outcome'],'reason':values[0]['reason'],'verified':values[0]['verified'],
                       'completion_ns':statistics.median(s['total_ns'] for s in values) if completed else None,
                       'observed_total_ns':statistics.median(s['total_ns'] for s in values),
                       'median_solve_ns':statistics.median(s['solve_ns'] for s in values),'median_kernel_ns':statistics.median(s['kernel_ns'] for s in values),
                       'kernel_fraction':statistics.median(s['kernel_ns']/s['solve_ns'] for s in values),
                       **{k:values[0][k] for k in logical if k not in ['outcome','reason']},
                       'kernel_calls_by_active':values[0]['kernel_calls_by_active'],'flat_calls':values[0]['flat_calls']}
        summaries.append({'cell':cell,'n':n,'family':fixture['family'],'seed':fixture['seed'],'split':fixture['split'],'arms':arms,'peak_rss_bytes':receipt['peak_rss_bytes']})
    gates=[]
    for n in protocol['variables']:
        if n==12:continue
        for family in protocol['families']:
            cells=[c for c in summaries if c['n']==n and c['family']==family and c['split']=='holdout']
            for candidate in ['search','bucket','hybrid']:
                ratios=[]
                for cell in cells:
                    pairs={(s['rep'],s['variant']):s for s in raw[cell['cell']][1:]}
                    for rep in range(protocol['repetitions']):
                        a,b=pairs[rep,'flat'],pairs[rep,candidate]
                        if all(s['outcome']!='UNKNOWN' and s['verified'] for s in [a,b]):ratios.append(a['total_ns']/b['total_ns'])
                complete=len(ratios)==len(cells)*protocol['repetitions'];ci=bootstrap(ratios) if complete else None
                gates.append({'n':n,'family':family,'candidate':candidate,'control':'flat','paired_ratio_median':statistics.median(ratios) if complete else None,
                              'ci95_paired_median':ci,'threshold':2.0,'pass':complete and ci[0]>2.0,'complete':complete})
    strongest=[]
    for n in protocol['variables']:
        if n==12:continue
        for family in protocol['families']:
            cells=[c for c in summaries if c['n']==n and c['family']==family and c['split']=='holdout'];ratios=[]
            for cell in cells:
                pairs={(s['rep'],s['variant']):s for s in raw[cell['cell']][1:]}
                for rep in range(protocol['repetitions']):
                    values=[pairs[rep,a] for a in protocol['variants']]
                    if all(v['outcome']!='UNKNOWN' and v['verified'] for v in values):ratios.append(min(pairs[rep,a]['total_ns'] for a in ['search','flat','bucket'])/pairs[rep,'hybrid']['total_ns'])
            complete=len(ratios)==len(cells)*protocol['repetitions'];ci=bootstrap(ratios) if complete else None
            strongest.append({'n':n,'family':family,'paired_ratio_median':statistics.median(ratios) if complete else None,'ci95_paired_median':ci,'pass':complete and ci[0]>2.0,'complete':complete})
    result={'schema_version':1,'scope':protocol['scope'],'correctness':'PASS','cells':len(summaries),'samples':sum(len(v)-1 for v in raw.values()),
            'all_results_completed_and_verified':all_complete,
            'full_solve_gates':{arm:'PASS' if all_complete and all(g['pass'] for g in gates if g['candidate']==arm) else 'REJECTED' for arm in ['search','bucket','hybrid']},
            'hybrid_vs_fastest_control_gate':'PASS' if all_complete and all(g['pass'] for g in strongest) else 'REJECTED',
            'gate_details':gates,'strongest_control_details':strongest,'cells_detail':summaries,
            'classification':protocol['classification'],'production_solver_cost':None,'full_ic_cost':None,'rho_ratio':None}
    (root/'results.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    lines=['# Complete generic Boolean solve comparison','',f"Evidence integrity: **PASS** across {result['cells']} cells and {result['samples']} full-solve observations.",'',
           f"All results completed and independently verified: **{all_complete}**. Full-solve gates against flat: **{result['full_solve_gates']}**. Hybrid versus fastest control: **{result['hybrid_vs_fastest_control_gate']}**.",'',
           'The table reports cold solve plus result-validation milliseconds over two holdout seeds and four rotated repetitions. A completion cost is null if any sample in the group is censored or lacks the required verification. Observed capped work is retained separately. Ratios use pooled medians; acceptance uses paired intervals.','',
           '| Variables | Family | Arm | Complete cost (ms) | Observed cost (ms) | Flat / arm | Median nodes | Median kernel calls | Median kernel fraction |',
           '|---:|---|---|---:|---:|---:|---:|---:|---:|']
    for n in protocol['variables']:
        for family in protocol['families']:
            base=groups['holdout',n,family,'flat'];base_ok=all(s['outcome']!='UNKNOWN' and s['verified'] for s in base)
            baseline=statistics.median(s['total_ns'] for s in base)
            for arm in protocol['variants']:
                values=groups['holdout',n,family,arm];okay=all(s['outcome']!='UNKNOWN' and s['verified'] for s in values);time=statistics.median(s['total_ns'] for s in values)
                cost=f'{time/1e6:.6f}' if okay else 'null';ratio=f'{baseline/time:.3f}' if okay and base_ok else 'null'
                lines.append(f"| {n} | {family} | {arm} | {cost} | {time/1e6:.6f} | {ratio} | {statistics.median(s['nodes'] for s in values):.0f} | {statistics.median(s['kernel_calls'] for s in values):.0f} | {statistics.median(s['kernel_ns']/s['solve_ns'] for s in values):.3f} |")
    lines+=['','Kernel-backed solvers have equal outcomes/models, logical counters and trace digests on every matched cell. The search-only control may explore a different tree. SAT models are checked against original equations; UNSAT verification requires a completed independent search reference. UNKNOWN is censored, not UNSAT.','',
            'This driver solves bounded generated Boolean systems. It does not benchmark the repository inherited-F4 implementation, accept curve targets, recover scalars or establish index-calculus performance. Production and cryptanalytic costs stay null.']
    (root/'RESULT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:result[k] for k in ['correctness','cells','samples','all_results_completed_and_verified','full_solve_gates','hybrid_vs_fastest_control_gate']},sort_keys=True))

if __name__=='__main__':main(Path(sys.argv[1]).resolve())
