#!/usr/bin/env python3
"""Verify the discovery-only cost split and conditional unchanged-scan ceilings."""
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import random
import statistics
import sys
sys.dont_write_bytecode=True

PHASES=['encoding_ns','construction_ns','scan_ns','release_ns']
KINDS=['16','64','dispatch']
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def satisfies(polys,x):return all(sum((x&m)==m for m in p)%2==0 for p in polys)
def interval(values):
    if not values or any(v is None for v in values):return None
    rng=random.Random(20260924)
    medians=sorted(statistics.median(rng.choices(values,k=len(values))) for _ in range(4000))
    return [medians[100],medians[3899]]
def median(values):return None if not values or any(v is None for v in values) else statistics.median(values)
def order(count,seed,rep):
    mask=(1<<64)-1;state=seed^(((rep//2)*0xd1b54a32d192ed03)&mask);items=list(range(count))
    for i in range(count-1,0,-1):
        state=(state+0x9e3779b97f4a7c15)&mask;z=state
        z=((z^(z>>30))*0xbf58476d1ce4e5b9)&mask;z=((z^(z>>27))*0x94d049bb133111eb)&mask;z^=z>>31
        j=z%(i+1);items[i],items[j]=items[j],items[i]
    return items[::-1] if rep%2 else items
def original(kind,n):return 'word_dispatch' if kind=='dispatch' else 'word16_unrolled' if kind=='16' else 'wide64_unrolled'
def adjusted_scan(scan_ns,profile_total,rebound_total,clock_max):
    return max(0,scan_ns-max(0,profile_total-rebound_total)-clock_max)

def construction_arm(kind,suffix):return ('construction_dispatch_' if kind=='dispatch' else 'construction'+kind+'_')+suffix

def main(root, output):
    protocol=read(root/'protocol.json');meta=read(root/'metadata.json');assert meta['complete']
    for name,digest in meta['source_hashes'].items():assert sha(root/name)==digest,name
    predecessor=read(root/'REFERENCE_PROTOCOL.json')
    extra=[construction_arm(k,suffix) for k in KINDS for suffix in ['rebound','profile']]
    assert protocol['retained_arms']==predecessor['variants']
    assert protocol['variants']==predecessor['variants']+extra
    assert protocol['reference_arms']==predecessor['variants']+extra[::2]
    assert protocol['variables']==[12,16,20,24] and protocol['seeds']==[17,937] and protocol['repetitions']==8
    assert protocol['clock_samples']==4096
    check=read(root/'test_receipt.json');assert check['exit_code']==0 and not check['timed_out']
    for stream in ['stdout','stderr']:assert hashlib.sha256(check[stream].encode()).hexdigest()==check[stream+'_sha256']
    expected={f'n{n}-{seed}-{family}' for n in protocol['variables'] for seed in protocol['seeds'] for family in protocol['families']}
    receipts=read(root/'receipts.json');assert len(receipts)==len(expected)
    receipts={r['cell']:r for r in receipts};assert set(receipts)==expected
    cells=[];groups=defaultdict(list);costs=defaultdict(list);all_complete=True
    for cell,receipt in sorted(receipts.items()):
        path=root/(cell+'.jsonl');assert sha(path)==receipt['stdout_sha256']
        assert sha(root/(cell+'.stderr'))==receipt['stderr_sha256']==hashlib.sha256(b'').hexdigest()
        assert receipt['exit_code']==0 and not receipt['timed_out']
        rows=[json.loads(line) for line in path.read_text().splitlines()];fixture=rows[0];samples=rows[1:]
        assert fixture['type']=='fixture'
        n,seed,family=[fixture[k] for k in ['n','seed','family']]
        for k in ['n','seed','family']:assert fixture[k]==receipt[k]
        seed_order=(protocol['order_seed']^(n<<48)^(seed<<8)^protocol['families'].index(family))&((1<<64)-1)
        assert receipt['order_seed']==seed_order and receipt['command'][-1]==str(seed_order)
        assert fixture['clock_samples']==4096 and 0<=fixture['clock_median_ns']<=fixture['clock_p99_ns']<=fixture['clock_max_ns']
        polys=fixture['polys'];assert len(polys)==n+2 and all(all(0<=m<1<<n and m.bit_count()<=2 for m in p) for p in polys)
        if fixture['planted_witness'] is not None:assert satisfies(polys,fixture['planted_witness'])
        names=protocol['variants'];pairs={(r['rep'],r['variant']):r for r in samples}
        assert len(samples)==len(pairs)==len(names)*8 and set(pairs)=={(rep,name) for rep in range(8) for name in names}
        for rep in range(8):
            observed=[r for r in samples if r['rep']==rep]
            assert [r['variant'] for r in observed]==[names[i] for i in order(len(names),seed_order,rep)]
            assert [r['order'] for r in observed]==list(range(len(names)))
        for row in samples:
            name=row['variant'];status=row['outcome'];model=row['model']
            assert row['type']=='sample' and 0<=row['solve_ns']<=row['total_ns']
            assert 0<=row['validation_ns'] and row['solve_ns']+row['validation_ns']<=row['total_ns']
            if status=='SAT':assert row['verified'] and isinstance(model,int) and 0<=model<1<<n and satisfies(polys,model) and row['reason'] is None
            elif status=='UNSAT':assert row['verified'] and model is None and row['reason'] is None and fixture['search_reference']=='UNSAT'
            else:assert status=='UNKNOWN' and model is None and not row['verified'] and row['reason'] is not None
            if not row['verified'] or status=='UNKNOWN':all_complete=False
            phases=row['phases']
            if name.endswith('_profile'):
                assert phases is not None and set(phases)==set(PHASES)
                assert all(isinstance(v,int) and v>=0 for v in phases.values()) and sum(phases.values())<=row['solve_ns']
            else:assert phases is None
            if name.startswith('construction'):assert row['trace'] is None
            signature={k:row[k] for k in ['outcome','model','reason','verified','logical','trace']}
            baseline=pairs[0,name]
            assert signature=={k:baseline[k] for k in signature}
            costs[n,family,name].append(row)
        details=[]
        for kind in KINDS:
            for rep in range(8):
                measured=pairs[rep,construction_arm(kind,'profile')];rebound=pairs[rep,construction_arm(kind,'rebound')];old=pairs[rep,original(kind,n)]
                for key in ['outcome','model','reason','logical','trace']:
                    assert measured[key]==rebound[key]==old[key],(cell,kind,key)
                complete=all(row['verified'] and row['outcome']!='UNKNOWN' for row in [measured,rebound,old]+[pairs[rep,a] for a in protocol['reference_arms']])
                floor=adjusted_scan(measured['phases']['scan_ns'],measured['total_ns'],rebound['total_ns'],fixture['clock_max_ns'])
                reference=min(pairs[rep,a]['total_ns'] for a in protocol['reference_arms']) if complete else None
                detail={'kind':kind,'rep':rep,'complete':complete,'reference_ns':reference,'phases':measured['phases'],
                    'profile_over_rebound':measured['total_ns']/rebound['total_ns'] if complete else None,
                    'rebound_over_original':rebound['total_ns']/old['total_ns'] if complete else None,
                    'raw_ceiling':reference/measured['phases']['scan_ns'] if complete and measured['phases']['scan_ns']>0 else None,
                    'adjusted_scan_ns':floor,'optimistic_ceiling':reference/floor if complete and floor>0 else None}
                details.append(detail);groups[n,family,kind].append(detail)
        cells.append({'cell':cell,'n':n,'seed':seed,'family':family,'fixture':fixture,'details':details,'peak_rss_bytes':receipt['peak_rss_bytes']})
    bounds=[]
    for (n,family,kind),rows in sorted(groups.items()):
        p=interval([v['profile_over_rebound'] for v in rows]);r=interval([v['rebound_over_original'] for v in rows]);c=interval([v['optimistic_ceiling'] for v in rows])
        comparable=all(v['complete'] for v in rows) and p is not None and r is not None and p[0]>=0.8 and p[1]<=1.25 and r[0]>=0.8 and r[1]<=1.25
        decision='RULED_OUT_UNDER_MEASURED_SCAN' if comparable and c is not None and c[1]<2 else 'NOT_RULED_OUT' if comparable and c is not None else 'INCONCLUSIVE'
        bounds.append({'n':n,'family':family,'kind':kind,'profile_over_rebound_ci95':p,'rebound_over_original_ci95':r,'optimistic_ceiling_median':median([v['optimistic_ceiling'] for v in rows]),
            'optimistic_ceiling_ci95':c,'raw_ceiling_median':median([v['raw_ceiling'] for v in rows]),'comparable':comparable,'decision':decision,
            'phase_medians_ns':{k:statistics.median(v['phases'][k] for v in rows) for k in PHASES}})
    table=[]
    for (n,family,name),rows in sorted(costs.items()):
        complete=all(r['verified'] and r['outcome']!='UNKNOWN' for r in rows)
        table.append({'n':n,'family':family,'variant':name,'completion_ns':statistics.median(r['total_ns'] for r in rows) if complete else None,
            'observed_ns':statistics.median(r['total_ns'] for r in rows),'correct':complete})
    result={'schema_version':1,'scope':protocol['scope'],'classification':'engineering diagnostic','all_complete_verified':all_complete,
        'cells':len(cells),'observations':len(cells)*len(protocol['variants'])*8,'bounds':bounds,'cells_detail':cells,'costs':table,
        'constructor_optimization_implemented':False,'performance_promotion':None,'full_ic_cost':None,'production_solver_cost':None,'calibrated_operation_ratio':None,'rho_ratio':None}
    (output/'results.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    lines=['# Discovery cost partition and unchanged-scan ceilings','',f"{result['cells']} cells, {result['observations']} observations; complete and verified: {all_complete}.",'',
        'Ceilings are conditional timing diagnostics. No constructor optimization or speedup is measured. Null ceilings and failed comparability guards remain inconclusive.','',
        '| Variables | Family | Policy | Median optimistic ceiling | 95% interval | Comparability | Decision |','|---:|---|---|---:|---|---|---|']
    for v in bounds:lines.append(f"| {v['n']} | {v['family']} | {v['kind']} | {v['optimistic_ceiling_median']} | {v['optimistic_ceiling_ci95']} | {v['comparable']} | {v['decision']} |")
    (output/'RESULT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'cells':len(cells),'observations':result['observations'],'all_complete_verified':all_complete,
        'relevant_decisions':{k:{d:sum(v['n']>=16 and v['kind']==k and v['decision']==d for v in bounds) for d in ['RULED_OUT_UNDER_MEASURED_SCAN','NOT_RULED_OUT','INCONCLUSIVE']} for k in KINDS}}))

if __name__=='__main__':
    output=Path(sys.argv[2]);output.mkdir(parents=True,exist_ok=False);main(Path(sys.argv[1]),output)
