#!/usr/bin/env python3
"""Replay bounded necessary filters and their conditional preprocessing costs."""
from collections import defaultdict
from functools import lru_cache
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


def mobius(value, variables):
    full=(1<<(1<<variables))-1
    for j in range(variables):
        shift=1<<j
        upper=(full//((1<<shift)+1))<<shift
        value ^= (value<<shift)&upper
    return value&full


def pack(words):
    assert all(type(w) is int and 0<=w<1<<64 for w in words)
    return sum(word<<(64*j) for j,word in enumerate(words))


def unpack(value,h):
    return [(value>>(64*j))&((1<<64)-1) for j in range(((1<<h)+63)//64)]


def coefficient_projection(polys):
    coefficients={}
    for e,poly in enumerate(polys):
        for m in poly:coefficients[m]=coefficients.get(m,0)^(1<<e)
    basis=[]
    def project(value):
        for marker,row in basis:
            if value&marker:value^=row
        return value
    for j in range(4):
        for i in range(j):
            value=project(coefficients.get((1<<i)|(1<<j),0))
            if value:basis.append((value&-value,value))
    assert all(project(coefficients.get((1<<i)|(1<<j),0))==0 for j in range(4) for i in range(j))
    return basis,project


def determinant(rows):
    # Independent column-basis rank rather than the worker's row elimination.
    basis={}
    for col in range(5):
        value=sum(((row>>col)&1)<<i for i,row in enumerate(rows))
        while value:
            pivot=value.bit_length()-1
            if pivot not in basis:
                basis[pivot]=value
                break
            value^=basis[pivot]
    return int(len(basis)==5)


def consistent(columns,rhs):
    basis={}
    for value in columns:
        while value:
            pivot=value.bit_length()-1
            if pivot not in basis:
                basis[pivot]=value
                break
            value^=basis[pivot]
    while rhs:
        pivot=rhs.bit_length()-1
        if pivot not in basis:return False
        rhs^=basis[pivot]
    return True


@lru_cache(maxsize=16)
def oracle(polys,n):
    # Cache only this pure validation oracle under the complete ordered input.
    # The measured compiler and scan have no cache across coefficient assignments.
    basis,project=coefficient_projection(polys)
    pivots=sum(marker for marker,_ in basis)
    available=[e for e in range(len(polys)) if not pivots&(1<<e)]
    assert len(available)>=5
    rows=[[available[(offset+j)%len(available)] for j in range(5)] for offset in range(4)]
    h=n-4
    original_truth=[]
    original_unsatisfied=0
    for poly in polys:
        coefficients=0
        for m in poly:coefficients^=1<<m
        truth=mobius(coefficients,n)
        original_truth.append(truth)
        original_unsatisfied|=truth
    original_solutions=((1<<(1<<n))-1)^original_unsatisfied
    minors=[0]*4
    affine=0
    def direct(point):
        return project(sum(((truth>>point)&1)<<e for e,truth in enumerate(original_truth)))
    for outside in range(1<<h):
        point=outside<<4
        rhs=direct(point)
        columns=[direct(point|(1<<i))^rhs for i in range(4)]
        for r,selected in enumerate(rows):
            matrix=[sum(((c>>e)&1)<<i for i,c in enumerate(columns))|(((rhs>>e)&1)<<4) for e in selected]
            minors[r]|=determinant(matrix)<<outside
        if consistent(columns,rhs):affine|=1<<outside
    full=(1<<(1<<h))-1
    rejected=0
    for m in minors:rejected|=m
    accepted=full^rejected
    assert affine&~accepted==0
    solutions=original_solutions
    while solutions:
        bit=solutions&-solutions;point=bit.bit_length()-1;solutions^=bit
        assert accepted&(1<<(point>>4)) and affine&(1<<(point>>4))
    return {'rows':rows,'annihilated_rank':len(basis),'quotient_dimension':len(polys)-len(basis),
        'original_solutions':original_solutions.bit_count(),
        'oracle_minor':{'accepted':unpack(accepted,h),'minors':[unpack(m,h) for m in minors]},
        'oracle_affine':{'accepted':unpack(affine,h),'minors':[]}}


def filter_statistics(value,h):
    truth=[pack(m) for m in value['minors']]
    full=(1<<(1<<h))-1
    union=0
    for m in truth:
        assert m<=full
        union|=m
    accepted=pack(value['accepted'])
    assert accepted<=full
    if truth:assert accepted==full^union
    coeff=[mobius(m,h) for m in truth]
    degrees=[]
    for c in coeff:
        degree=None
        while c:
            bit=c&-c;c^=bit
            d=(bit.bit_length()-1).bit_count()
            degree=d if degree is None else max(degree,d)
        degrees.append(degree)
    return {'accepted_prefixes':accepted.bit_count(),'individual_rejected':[m.bit_count() for m in truth],
        'joint_rejected':union.bit_count() if truth else (1<<h)-accepted.bit_count(),
        'overlap':[[ (x&y).bit_count() for y in truth] for x in truth],
        'zero_minors':sum(m==0 for m in truth),'duplicate_minors':len(truth)-len(set(truth)),
        'degrees':degrees,'supports':[c.bit_count() for c in coeff]}


def check_work(work,arm,expected,h,stats,complete):
    for key in ['annihilated_rank','quotient_dimension']:assert work[key]==expected[key]
    assert all(type(v) is int and v>=0 for k,v in work.items() if k not in ['degrees','supports'])
    assert work['product_toggles']<=50000000
    assert work['product_toggles']-2*work['cancelled_toggles']==work['support_sum']
    if not complete:return
    if arm=='symbolic_minor4':
        assert work['dp_states']==124
        assert work['dp_payload_peak_bytes']==32*((1<<h)//64)*8
        assert work['truth_word_xors']==4*(6*((1<<h)//64)+(h-6)*((1<<h)//128))
        assert work['degrees']==stats['degrees'] and work['supports']==stats['supports']
        assert all(d is None or 0<=d<=6 for d in work['degrees'])
        assert max(work['supports'])<=work['max_support']<=1<<h
        assert work['prefixes']==work['numeric_determinants']==work['affine_queries']==0
    else:
        assert not work['degrees'] and not work['supports']
        assert work['product_toggles']==work['dp_states']==work['truth_word_xors']==0
        assert work['prefixes']==1<<h
        if arm=='numeric_minor4':
            assert work['numeric_determinants']==4*(1<<h)
            assert work['affine_queries']==work['affine_rejected']==work['screen_rejected']==0
        else:
            assert work['numeric_determinants']==0
            assert work['screen_rejected']+work['affine_queries']==1<<h
            assert work['affine_queries']-work['affine_rejected']==stats['accepted_prefixes']
            assert work['rank_sum']<=4*work['affine_queries']


def validate(root):
    protocol,meta=read(root/'protocol.json'),read(root/'metadata.json')
    assert meta['complete']
    for name,digest in meta['source_hashes'].items():assert sha(root/name)==digest,name
    previous=read(root/'REFERENCE_PROTOCOL.json')
    assert protocol['retained_arms']==protocol['reference_arms']==previous['variants']
    assert protocol['filter_arms']==['numeric_minor4','symbolic_minor4','affine_filter']
    assert protocol['variants']==previous['variants']+protocol['filter_arms']
    assert protocol['variables']==[12,16] and protocol['seeds']==[17,937] and protocol['repetitions']==8
    assert protocol['families']==['planted','cross_planted','unplanted']
    assert protocol['low_variables']==4 and protocol['minors']==4
    lineage=read(root/'SOURCE_LINEAGE.json')
    assert sha(root/'reference.py')==lineage['reference_python']['sha256']
    for name,digest in lineage['files'].items():
        data=(root/name).read_bytes()
        if name=='worker.rs':
            suffix=b'\ninclude!("minors.rs");\ninclude!("minors_main.rs");\n'
            assert data.endswith(suffix);data=data[:-len(suffix)]
        assert hashlib.sha256(data).hexdigest()==digest,name
    test=read(root/'test_receipt.json')
    assert test['exit_code']==0 and not test['timed_out']
    for stream in ['stdout','stderr']:assert hashlib.sha256(test[stream].encode()).hexdigest()==test[stream+'_sha256']
    receipts=read(root/'receipts.json')
    expected={(n,seed,f):f'n{n}-{seed}-{f}' for n in [12,16] for seed in [17,937] for f in protocol['families']}
    assert len(receipts)==len(expected)==12
    receipts={r['cell']:r for r in receipts};assert set(receipts)==set(expected.values())
    summaries=[];groups=defaultdict(list);bounds=defaultdict(list);all_complete=True;count=0
    for (n,seed,family),cell in expected.items():
        receipt=receipts[cell]
        assert receipt['exit_code']==0 and not receipt['timed_out']
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
        exact=oracle(tuple(tuple(p) for p in polys),n)
        for key in ['rows','original_solutions','oracle_minor','oracle_affine']:assert source[key]==exact[key],(cell,key)
        stats={arm:filter_statistics(exact['oracle_affine' if arm=='affine_filter' else 'oracle_minor'],n-4) for arm in protocol['filter_arms']}
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
                assert all(type(sample[k]) is int and sample[k]>=0 for k in ['compute_ns','validation_ns','total_ns'])
                assert 0<sample['compute_ns']+sample['validation_ns']<=sample['total_ns']
                if arm in protocol['filter_arms']:
                    assert sample['kind']=='filter' and status in ['COMPLETE','CAPPED']
                    assert sample['model'] is None and sample['logical'] is None and sample['projected_work'] is None and sample['trace'] is None
                    phases=sample['phases'];assert all(type(v) is int and v>=0 for v in phases.values())
                    assert set(phases)=={'setup_ns','construction_ns','evaluation_ns'}
                    assert sum(phases.values())<=sample['compute_ns']
                    if status=='COMPLETE':
                        assert sample['verified'] and sample['reason'] is None
                        assert sample['filter']==exact['oracle_affine' if arm=='affine_filter' else 'oracle_minor']
                    else:
                        assert not sample['verified'] and sample['filter'] is None and sample['reason'] is not None
                        all_complete=False
                    check_work(sample['work'],arm,exact,n-4,stats[arm],status=='COMPLETE')
                else:
                    assert sample['kind']=='solver' and status in ['SAT','UNSAT','UNKNOWN']
                    assert sample['phases'] is None and sample['work'] is None and sample['filter'] is None
                    if status=='SAT':
                        assert type(sample['model']) is int and 0<=sample['model']<1<<n
                        assert sample['verified'] and sample['reason'] is None and ref.satisfies(polys,sample['model']) and exact['original_solutions']>0
                    elif status=='UNSAT':
                        assert sample['model'] is None and sample['reason'] is None and sample['verified'] and exact['original_solutions']==0
                    else:
                        assert sample['model'] is None and sample['reason'] is not None and not sample['verified'];all_complete=False
                    if arm.startswith('projected'):
                        ref.check_work(sample['projected_work'],polys,n,int(arm[-1]),status)
                    else:assert sample['projected_work'] is None
                signature={k:v for k,v in sample.items() if k not in ['rep','order','compute_ns','validation_ns','total_ns','phases']}
                if arm in stable:assert stable[arm]==signature
                stable[arm]=signature
                groups[n,family,arm].append(sample['total_ns'])
            best=min(pairs[rep,arm]['total_ns'] for arm in protocol['reference_arms'])
            bounds[n,family].append(best/pairs[rep,'symbolic_minor4']['compute_ns'])
        count+=len(samples)
        summaries.append({'cell':cell,'n':n,'seed':seed,'family':family,'original_solutions':exact['original_solutions'],
            'annihilated_rank':exact['annihilated_rank'],'quotient_dimension':exact['quotient_dimension'],
            'statistics':stats,'symbolic_work':pairs[0,'symbolic_minor4']['work'],
            'complete':all(s['outcome'] not in ['UNKNOWN','CAPPED'] for s in samples),
            'medians_ns':{arm:statistics.median(pairs[rep,arm]['total_ns'] for rep in range(8)) for arm in names}})
    gates=[]
    for (n,family),values in sorted(bounds.items()):
        ci=ref.interval(values)
        decision='INCONCLUSIVE' if not all_complete else ('BLOCKS_MEASURED_FIXED_PREPROCESSING' if ci[1]<2 else ('OPPORTUNITY_ONLY' if ci[0]>2 else 'INCONCLUSIVE'))
        gates.append({'n':n,'family':family,'median':statistics.median(values),'ci95':ci,'decision':decision})
    return {'schema_version':1,'cells':len(summaries),'observations':count,'all_complete':all_complete,'summaries':summaries,'gates':gates,
        'complete_candidate_justified':all_complete and all(g['decision']=='OPPORTUNITY_ONLY' for g in gates),
        'complete_candidate_implemented':False,'n16_ms':[{'variant':a,'kind':'filter' if a in protocol['filter_arms'] else 'solver',**{f:statistics.median(groups[16,f,a])/1e6 for f in protocol['families']}} for a in names],
        'campaign_seconds':meta['campaign_seconds'],'peak_worker_rss_bytes':max(r['peak_rss_bytes'] for r in receipts.values()),
        'full_ic_cost':None,'production_solver_cost':None,'rho_ratio':None,'calibrated_operation_ratio':None}


def main(root):
    result=validate(root);dump(root/'results.json',result)
    print(json.dumps({k:result[k] for k in ['cells','observations','all_complete','gates','complete_candidate_justified','campaign_seconds','peak_worker_rss_bytes']}))


if __name__=='__main__':main(Path(sys.argv[1]).resolve())
