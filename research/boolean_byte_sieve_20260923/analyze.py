#!/usr/bin/env python3
"""Verify complete-solve evidence without promoting censored outcomes."""
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import random
import statistics
import sys
import importlib.util
sys.dont_write_bytecode=True

def read(path):return json.loads(path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def bootstrap(values):
    rng=random.Random(20260923)
    medians=sorted(statistics.median(rng.choices(values,k=len(values))) for _ in range(4000))
    return [medians[100],medians[3899]]
def satisfies(polys,model):return all(sum((model&m)==m for m in p)%2==0 for p in polys)

QUIET=['packed_untraced','gray_quiet','wide64_quiet','wide64_unrolled','word16_unrolled','word_dispatch','leaf16_quiet','leaf16_word_unrolled','byte_scalar','byte_simd','byte_planes','byte_unrolled','byte_single_quiet','leaf16_byte_scalar','leaf16_byte_simd','leaf16_byte_planes','leaf16_byte_unrolled','leaf16_single_quiet']
def is_partial(arm):return arm.startswith('byte_') or arm.startswith('leaf16_byte_') or arm=='leaf16_single_quiet'
def is_enumeration(arm):return arm.startswith(('gray_','packed_gray','initial_')) or arm in ['wide64_quiet','wide64_unrolled','word16_unrolled','word_dispatch','leaf16_quiet','leaf16_word_unrolled']
def paired_order(count,seed,rep):
    mask=(1<<64)-1;state=seed^(((rep//2)*0xd1b54a32d192ed03)&mask);order=list(range(count))
    for i in range(count-1,0,-1):
        state=(state+0x9e3779b97f4a7c15)&mask;z=state
        z=((z^(z>>30))*0xbf58476d1ce4e5b9)&mask;z=((z^(z>>27))*0x94d049bb133111eb)&mask;z^=z>>31
        j=z%(i+1);order[i],order[j]=order[j],order[i]
    return order[::-1] if rep%2 else order

def main(root):
    protocol,meta=read(root/'protocol.json'),read(root/'metadata.json');assert meta['complete']
    for name,digest in meta['source_hashes'].items():assert sha(root/name)==digest,name
    test=read(root/'test_receipt.json');assert test['exit_code']==0 and not test['timed_out']
    for stream in ['stdout','stderr']:assert hashlib.sha256(test[stream].encode()).hexdigest()==test[stream+'_sha256']
    expected={f'n{n}-{split}-{seed}-{family}' for n in protocol['variables'] for split in protocol.get('splits',['discovery','holdout']) for seed in protocol[split+'_seeds'] for family in protocol['families']}
    receipts=read(root/'receipts.json');assert len(receipts)==len(expected)
    receipts={r['cell']:r for r in receipts};assert set(receipts)==expected
    raw=defaultdict(list)
    for n in protocol['variables']:
        for line in (root/f'raw-n{n}.jsonl').read_text().splitlines():
            row=json.loads(line);assert json.loads(row['raw_line'])=={k:v for k,v in row.items() if k not in ['cell','split','raw_line']}
            raw[row['cell']].append(row)
    assert set(raw)==expected
    predecessor=read(root/'REFERENCE_PROTOCOL.json')
    assert protocol['retained_arms']==predecessor['variants']
    assert protocol['variants']==protocol['retained_arms']+['gray_quiet', 'wide64_quiet', 'byte_scalar', 'byte_simd', 'leaf16_quiet', 'leaf16_byte_scalar', 'leaf16_byte_simd', 'byte_single_quiet', 'leaf16_single_quiet', 'byte_planes', 'leaf16_byte_planes', 'byte_unrolled', 'wide64_unrolled', 'leaf16_byte_unrolled', 'word16_unrolled', 'word_dispatch', 'leaf16_word_unrolled']
    assert protocol['reference_arms']==protocol['retained_arms']+['gray_quiet', 'wide64_quiet', 'byte_scalar', 'leaf16_quiet', 'leaf16_byte_scalar', 'byte_single_quiet', 'leaf16_single_quiet']
    assert protocol['new_candidates']==['wide64_quiet', 'wide64_unrolled', 'word16_unrolled', 'word_dispatch', 'leaf16_word_unrolled', 'byte_simd', 'byte_planes', 'byte_unrolled', 'leaf16_byte_simd', 'leaf16_byte_planes', 'leaf16_byte_unrolled']
    assert protocol['matched_controls']=={'wide64_quiet': 'gray_quiet', 'wide64_unrolled': 'gray_quiet', 'word16_unrolled': 'gray_quiet', 'word_dispatch': 'gray_quiet', 'leaf16_word_unrolled': 'leaf16_quiet', 'byte_simd': 'byte_scalar', 'byte_planes': 'byte_scalar', 'byte_unrolled': 'byte_scalar', 'leaf16_byte_simd': 'leaf16_byte_scalar', 'leaf16_byte_planes': 'leaf16_byte_scalar', 'leaf16_byte_unrolled': 'leaf16_byte_scalar'}
    assert protocol['repetitions']==8 and protocol['variables']==[12,16,20,24]
    assert protocol['splits']==['discovery','regression','holdout']
    assert protocol['discovery_seeds']==predecessor['discovery_seeds']
    assert protocol['regression_seeds']==predecessor['regression_seeds']+predecessor['holdout_seeds']
    assert len(protocol['holdout_seeds'])==len(set(protocol['holdout_seeds']))==2
    assert not set(protocol['holdout_seeds']).intersection(protocol['discovery_seeds']+protocol['regression_seeds'])
    index=0
    for n in protocol['variables']:
        for split in protocol['splits']:
            for seed in protocol[split+'_seeds']:
                for family in protocol['families']:
                    cell=f'n{n}-{split}-{seed}-{family}'
                    assert receipts[cell]['order_offset']==0
                    assert receipts[cell]['order_seed']==(protocol['order_seed']^(n<<48)^(seed<<8)^protocol['families'].index(family))&((1<<64)-1)
                    index+=1
    spec=importlib.util.spec_from_file_location('selection_oracle',root/'selection_reference.py')
    selection=importlib.util.module_from_spec(spec);spec.loader.exec_module(selection)
    summaries=[];groups=defaultdict(list);all_complete=True
    logical=['screen_calls', 'screen_points', 'screen_batches', 'screen_passes', 'screen_second_checks', 'screen_second_rejected', 'screen_full_checks', 'screen_full_rejected', 'screen_fallback_points', 'screen_fallback_batches','fiber_filter_rounds', 'fiber_selected', 'fiber_dimension', 'fiber_selection_states', 'fiber_prefixes', 'fiber_batches', 'fiber_queries', 'fiber_zero_rejected', 'fiber_linear_rejected', 'fiber_extensions_rejected', 'fiber_rank_sum','projection_rows','projection_xors','projection_certified','outcome','model','reason','nodes','kernel_calls','decisions','forced','specialized_terms','source_rows','source_columns','max_depth','affine_eliminated','derived_rows','enumeration_points','enumeration_batches','enumeration_leaves','trace']
    for cell,records in sorted(raw.items()):
        receipt=receipts[cell];assert receipt['exit_code']==0 and not receipt['timed_out']
        assert hashlib.sha256(''.join(r['raw_line'] for r in records).encode()).hexdigest()==receipt['stdout_sha256']
        assert receipt['stderr_sha256']==hashlib.sha256(b'').hexdigest()
        assert receipt['command'][-3:] == ['with-final','0',str(receipt['order_seed'])]
        fixture=records[0];assert fixture['type']=='fixture'
        for key in ['n','seed','split','family']:assert fixture[key]==receipt[key]
        n=fixture['n'];polys=fixture['polys'];assert len(polys)==n+2
        for p in polys:
            assert p==sorted(set(p),key=lambda m:(-m.bit_count(),m))
            assert len(p)<=32 and all(0<=m<1<<n and m.bit_count()<=2 for m in p)
        if fixture['planted_witness'] is not None:assert satisfies(polys,fixture['planted_witness'])
        graph=selection.interaction_graph(polys,n);selected,selection_work=selection.maximum_independent(graph)
        selection_states=1 if not any(graph) else selection_work['right_subset_states']+selection_work['left_subsets']
        samples=records[1:];names=protocol['variants'];pairs={(s['rep'],s['variant']):s for s in samples}
        assert len(samples)==len(names)*protocol['repetitions'] and len(pairs)==len(samples)
        assert set(pairs)=={(rep,arm) for rep in range(protocol['repetitions']) for arm in names}
        for rep in range(protocol['repetitions']):
            ordered=[s for s in samples if s['rep']==rep]
            assert [s['variant'] for s in ordered]==[names[i] for i in paired_order(len(names),receipt['order_seed'],rep)]
            assert [s['order'] for s in ordered]==list(range(len(names)))
        kernel_signature={}
        for sample in samples:
            arm=sample['variant'];status=sample['outcome'];model=sample['model']
            assert 0<=sample['kernel_ns']<=sample['solve_ns']<=sample['total_ns']
            assert sample['solve_ns']+sample['validation_ns']<=sample['total_ns']
            for field,enabled in [('substitution_ns',arm.startswith(('affine_','initial_'))),('enumeration_ns',is_enumeration(arm)),('selection_ns',arm.startswith('fiber_')),('fiber_setup_ns',arm.startswith('fiber_')),('fiber_ns',arm.startswith('fiber_')),('partial_ns',is_partial(arm))]:
                assert isinstance(sample[field],int) and sample[field]>=0 if enabled else sample[field] is None
            assert sample['kernel_ns']+sum(sample[k] or 0 for k in ['substitution_ns','enumeration_ns','selection_ns','fiber_setup_ns','fiber_ns','partial_ns'])<=sample['solve_ns']
            assert sample['enumeration_points']<=1<<n
            assert sample['enumeration_batches']<=sample['enumeration_points']<=(64 if arm in ['wide64_quiet','wide64_unrolled'] or (arm=='word_dispatch' and n>20) else 16)*sample['enumeration_batches']
            if not is_enumeration(arm):
                assert sample['enumeration_points']==sample['enumeration_batches']==0
                if not is_partial(arm):assert sample['enumeration_leaves']==0
            if not arm.startswith(('affine_','initial_')):assert sample['affine_eliminated']==sample['derived_rows']==0
            assert sample['affine_eliminated'] <= (n if arm.startswith('initial_') else n*sample['nodes'])
            if arm.startswith('initial_'):
                assert sample['nodes']==sample['decisions']==sample['forced']==sample['max_depth']==0
                assert sample['projection_certified'] in [0,1]
                assert sample['kernel_calls']==2-sample['projection_certified']
                assert sample['kernel_calls_by_active'][n]==sample['kernel_calls']
                assert 0<=sample['projection_xors']<=64*sample['projection_rows']
                assert 0<sample['projection_rows']<=len(polys)
                if sample['projection_certified']:assert sample['affine_eliminated']==sample['derived_rows']==0
                assert sample['affine_eliminated']<=sample['derived_rows']<=n+1
                assert sample['enumeration_points']<=1<<(n-sample['affine_eliminated'])
            if not arm.startswith('initial_'):assert sample['projection_rows']==sample['projection_xors']==sample['projection_certified']==0
            if arm in QUIET:assert sample['trace'] is None
            else:assert isinstance(sample['trace'],int)
            screen_fields=['screen_calls','screen_points','screen_batches','screen_passes','screen_second_checks','screen_second_rejected','screen_full_checks','screen_full_rejected','screen_fallback_points','screen_fallback_batches']
            if is_partial(arm):
                assert isinstance(sample['full_update_words'],int) and sample['full_update_words']>=0
                assert isinstance(sample['secondary_update_words'],int) and sample['secondary_update_words']>=0
                assert sample['screen_points']==64*sample['screen_batches']
                assert sample['screen_points']+sample['screen_fallback_points']<=1<<n
                assert sample['screen_fallback_batches']<=sample['screen_fallback_points']<=16*sample['screen_fallback_batches']
                assert sample['screen_full_rejected']<=sample['screen_full_checks']<=sample['screen_passes']<=sample['screen_points']
                assert sample['enumeration_points']==sample['enumeration_batches']==0
                if arm.startswith('byte_'):
                    assert sample['screen_calls']==1 and sample['nodes']==0
                    assert sample['screen_fallback_points']==sample['screen_fallback_batches']==0
                else:assert sample['screen_calls']==sample['enumeration_leaves']
                if arm in ['byte_single_quiet','leaf16_single_quiet']:
                    assert sample['screen_second_checks']==sample['screen_second_rejected']==sample['secondary_update_words']==0
                    assert sample['full_update_words']<=6*sample['screen_batches']
                else:
                    assert sample['screen_second_checks']==sample['screen_second_rejected']+sample['screen_full_checks']
                    assert sample['screen_second_checks']<=sample['screen_passes']
                    assert sample['full_update_words']==0 and sample['secondary_update_words']<=sample['screen_batches']
                if status=='UNSAT':
                    assert sample['screen_full_checks']==sample['screen_full_rejected']
                    if arm in ['byte_single_quiet','leaf16_single_quiet']:assert sample['screen_full_checks']==sample['screen_passes']
                    else:assert sample['screen_second_checks']==sample['screen_passes']
                    if arm.startswith('byte_'):assert sample['screen_points']==1<<n
                elif status=='SAT':
                    assert 0<=sample['screen_full_checks']-sample['screen_full_rejected']<=1
                    if arm.startswith('byte_'):assert sample['screen_full_checks']==sample['screen_full_rejected']+1
            else:
                assert all(sample[k]==0 for k in screen_fields)
                assert sample['full_update_words'] is None and sample['secondary_update_words'] is None
            fiber_fields=['fiber_filter_rounds','fiber_selected','fiber_dimension','fiber_selection_states','fiber_prefixes','fiber_batches','fiber_queries','fiber_zero_rejected','fiber_linear_rejected','fiber_extensions_rejected','fiber_rank_sum']
            if arm.startswith('fiber_'):
                assert sample['fiber_selected']==selected and sample['fiber_dimension']==selected.bit_count()
                assert sample['fiber_selection_states']==selection_states
                h=n-selected.bit_count();block=1<<min(h,4)
                assert sample['fiber_prefixes']==sample['fiber_batches']*block and sample['fiber_prefixes']<=1<<h
                assert sample['fiber_prefixes']<=sample['fiber_filter_rounds']<=8*sample['fiber_prefixes']
                if arm=='fiber_zero_simd':assert sample['fiber_filter_rounds']==sample['fiber_prefixes']
                assert sample['fiber_zero_rejected']+sample['fiber_linear_rejected']<=sample['fiber_prefixes']
                assert sample['fiber_linear_rejected']<=sample['fiber_queries']<=sample['fiber_prefixes']-sample['fiber_zero_rejected']
                assert sample['fiber_rank_sum']<=sample['fiber_queries']*sample['fiber_dimension']
                assert sample['fiber_extensions_rejected']==(sample['fiber_zero_rejected']+sample['fiber_linear_rejected'])<<sample['fiber_dimension']
                assert sample['kernel_calls']==sample['nodes']==sample['kernel_ns']==sample['decisions']==sample['forced']==0
                if status=='UNSAT':
                    assert sample['fiber_prefixes']==1<<h and sample['fiber_extensions_rejected']==1<<n
                    assert sample['fiber_queries']==sample['fiber_linear_rejected']
                elif status=='SAT':assert sample['fiber_queries']==sample['fiber_linear_rejected']+1
            else:assert all(sample[k]==0 for k in fiber_fields)
            assert sample['nodes']<=protocol['limits']['nodes'] and sample['max_depth']<=n
            assert sum(sample['kernel_calls_by_active'])==sample['kernel_calls']
            assert len(sample['kernel_calls_by_active'])==n+1
            assert sample['flat_calls']<=sample['kernel_calls']
            if arm in ('search','merge_search','quadratic_state','packed_state'):assert sample['kernel_calls']==sample['kernel_ns']==0
            if arm in ('flat','small_flat'):assert sample['flat_calls']==sample['kernel_calls']
            if arm in ('bucket','word_tail'):assert sample['flat_calls']==0
            if arm in ('small_flat','word_tail'):assert all(count==0 for width,count in enumerate(sample['kernel_calls_by_active']) if width>10)
            if status=='SAT':
                assert isinstance(model,int) and 0<=model<1<<n and satisfies(polys,model)
                assert sample['reason'] is None and sample['verified']
            elif status=='UNSAT':
                assert model is None and sample['reason'] is None and fixture['planted_witness'] is None
                assert sample['verified']==(fixture['search_reference']=='UNSAT')
            else:
                assert status=='UNKNOWN' and model is None and not sample['verified']
                assert sample['reason'] in ['NODE_CAP','ROW_CAP','COLUMN_CAP','ENUM_CAP','FIBER_CAP','SCREEN_CAP']
            if sample['verified'] and fixture['search_reference']!='UNKNOWN':assert status==fixture['search_reference']
            if status=='UNKNOWN' or not sample['verified']:all_complete=False
            if arm!='search' or protocol.get('enable_merge'):
                signature=tuple(sample[k] for k in logical)
                group=('screen_leaf_single' if arm=='leaf16_single_quiet' else 'screen_single' if arm=='byte_single_quiet' else 'screen_leaf' if arm.startswith('leaf16_byte_') else 'screen' if arm.startswith('byte_') else arm if arm in QUIET else None) or ('fiber_zero' if arm=='fiber_zero_simd' else 'fiber' if arm.startswith('fiber_') else 'untraced' if arm=='packed_untraced' else 'initial' if arm.startswith('initial_') else 'gray' if arm.startswith('gray_') else 'leaf12' if arm.startswith('packed_gray12') else 'leaf16' if arm.startswith('packed_gray16') else 'affine_sl_basis' if arm in ('affine_sl_basis_list','affine_sl_basis_fast') else 'basis' if arm in ('basis_list','basis_wide') else 'tail' if arm in ('tail_list','tail_wide') else 'search' if arm in ('search','merge_search','quadratic_state','packed_state') else 'degree2' if arm in ('small_flat','word_tail') else 'degree3')
                if group not in kernel_signature:kernel_signature[group]=signature
                else:assert signature==kernel_signature[group],(cell,arm)
            groups[fixture['split'],n,fixture['family'],arm].append(sample)
        baseline=next(s for s in samples if s['variant']=='packed_state')
        for sample in samples:
            if sample['variant']=='packed_untraced':assert all(sample[k]==baseline[k] for k in logical if k!='trace')
        arms={}
        for arm in names:
            values=[s for s in samples if s['variant']==arm]
            assert len({tuple(s[k] for k in logical) for s in values})==1
            completed=all(s['outcome']!='UNKNOWN' and s['verified'] for s in values)
            arms[arm]={'outcome':values[0]['outcome'],'reason':values[0]['reason'],'verified':values[0]['verified'],
                       'completion_ns':statistics.median(s['total_ns'] for s in values) if completed else None,
                       'observed_total_ns':statistics.median(s['total_ns'] for s in values),
                       'median_solve_ns':statistics.median(s['solve_ns'] for s in values),'median_kernel_ns':statistics.median(s['kernel_ns'] for s in values),
                       'median_substitution_ns':None if values[0]['substitution_ns'] is None else statistics.median(s['substitution_ns'] for s in values),
                       'median_enumeration_ns':None if values[0]['enumeration_ns'] is None else statistics.median(s['enumeration_ns'] for s in values),
                       **{'median_'+field:(None if values[0][field] is None else statistics.median(s[field] for s in values)) for field in ['selection_ns','fiber_setup_ns','fiber_ns','partial_ns']},
                       'full_update_words':values[0]['full_update_words'],'secondary_update_words':values[0]['secondary_update_words'],
                       'kernel_fraction':statistics.median(s['kernel_ns']/s['solve_ns'] for s in values),
                       **{k:values[0][k] for k in logical if k not in ['outcome','reason']},
                       'kernel_calls_by_active':values[0]['kernel_calls_by_active'],'flat_calls':values[0]['flat_calls']}
        for rep in range(protocol['repetitions']):
            for arm in ['gray_quiet','wide64_quiet','wide64_unrolled','word16_unrolled','word_dispatch','byte_scalar','byte_simd','byte_planes','byte_unrolled','byte_single_quiet']:
                a=pairs[rep,arm];b=pairs[rep,'gray_simd']
                if a['verified'] and b['verified'] and a['outcome']!='UNKNOWN':assert (a['outcome'],a['model'])==(b['outcome'],b['model'])
            for arm in ['leaf16_quiet','leaf16_word_unrolled','leaf16_byte_scalar','leaf16_byte_simd','leaf16_byte_planes','leaf16_byte_unrolled','leaf16_single_quiet']:
                a=pairs[rep,arm];b=pairs[rep,'packed_gray16_simd']
                if a['verified'] and b['verified'] and a['outcome']!='UNKNOWN':
                    for key in ['outcome','model','nodes','decisions','forced','specialized_terms','source_rows','source_columns','max_depth','enumeration_leaves']:assert a[key]==b[key],(cell,arm,key)
            a=pairs[rep,'word_dispatch'];b=pairs[rep,'word16_unrolled' if n<=20 else 'wide64_unrolled']
            assert all(a[k]==b[k] for k in logical),(cell,'dispatch work')
        summaries.append({'cell':cell,'n':n,'family':fixture['family'],'seed':fixture['seed'],'split':fixture['split'],'arms':arms,'peak_rss_bytes':receipt['peak_rss_bytes']})
    details=[]
    for candidate in protocol['new_candidates']:
        for split in ['regression','holdout']:
            for n in [16,20,24]:
                for family in protocol['families']:
                    cells=[c for c in summaries if (c['split'],c['n'],c['family'])==(split,n,family)]
                    assert cells
                    for reference,controls,threshold in [('retained_frontier',protocol['reference_arms'],2.0),('same_policy',[protocol['matched_controls'][candidate]],1.05)]:
                        controls=[a for a in controls if a!=candidate]
                        ratios=[]
                        for cell in cells:
                            pairs={(s['rep'],s['variant']):s for s in raw[cell['cell']][1:]}
                            for rep in range(protocol['repetitions']):
                                values=[pairs[rep,a] for a in controls+[candidate]]
                                if all(s['verified'] and s['outcome']!='UNKNOWN' for s in values):
                                    ratios.append(min(pairs[rep,a]['total_ns'] for a in controls)/pairs[rep,candidate]['total_ns'])
                        complete=len(ratios)==len(cells)*protocol['repetitions']
                        ci=bootstrap(ratios) if complete else None
                        details.append({'candidate':candidate,'split':split,'n':n,'family':family,'reference':reference,'controls':controls,'complete':complete,
                            'paired_ratio_median':statistics.median(ratios) if complete else None,'ci95_paired_median':ci,'threshold':threshold,
                            'pass':complete and ci[0]>threshold,'incremental_pass':complete and ci[0]>1.0})
    gates={c:{ref:'PASS' if all_complete and all(d['pass'] for d in details if d['candidate']==c and d['reference']==ref) else 'REJECTED'
              for ref in ['retained_frontier','same_policy']} for c in protocol['new_candidates']}
    for c in gates:
        gates[c]['incremental']='PASS' if all_complete and all(d['incremental_pass'] for d in details if d['candidate']==c and d['reference']=='retained_frontier') else 'REJECTED'
    costs=[]
    for n in protocol['variables']:
        for arm in protocol['variants']:
            row={'n':n,'variant':arm,'families':{}}
            for family in protocol['families']:
                values=groups['holdout',n,family,arm]
                complete=all(s['verified'] and s['outcome']!='UNKNOWN' for s in values)
                row['families'][family]={
                    'completion_ns':statistics.median(s['total_ns'] for s in values) if complete else None,
                    'observed_ns':statistics.median(s['total_ns'] for s in values),
                    'median_kernel_ns':statistics.median(s['kernel_ns'] for s in values),
                    'median_substitution_ns':None if values[0]['substitution_ns'] is None else statistics.median(s['substitution_ns'] for s in values),
                    'median_enumeration_ns':None if values[0]['enumeration_ns'] is None else statistics.median(s['enumeration_ns'] for s in values),
                    **{field:statistics.median(s[field] for s in values) for field in ['fiber_filter_rounds', 'fiber_selected', 'fiber_dimension', 'fiber_selection_states', 'fiber_prefixes', 'fiber_batches', 'fiber_queries', 'fiber_zero_rejected', 'fiber_linear_rejected', 'fiber_extensions_rejected', 'fiber_rank_sum']},
                    **{'median_'+field:(None if values[0][field] is None else statistics.median(s[field] for s in values)) for field in ['selection_ns','fiber_setup_ns','fiber_ns','partial_ns']},
                    **{k:statistics.median(s[k] for s in values) for k in ['screen_calls', 'screen_points', 'screen_batches', 'screen_passes', 'screen_second_checks', 'screen_second_rejected', 'screen_full_checks', 'screen_full_rejected', 'screen_fallback_points', 'screen_fallback_batches']},
                    'affine_rank':statistics.median(s['affine_eliminated'] for s in values),
                    'enumeration_points':statistics.median(s['enumeration_points'] for s in values)}
            costs.append(row)
    for row in costs:
        before=next(r for r in costs if r['n']==row['n'] and r['variant']=='gray_simd')['families']['planted']['completion_ns']
        after=row['families']['planted']['completion_ns']
        row['display_direct_simd_ratio']=before/after if before is not None and after is not None else None
    result={'schema_version':1,'scope':protocol['scope'],'classification':'engineering','all_complete_verified':all_complete,
        'cells':len(summaries),'observations':sum(len(v)-1 for v in raw.values()),'cells_detail':summaries,'comparisons':details,'gates':gates,'costs':costs,
        'full_ic_cost':None,'production_solver_cost':None,'rho_ratio':None,'calibrated_operation_ratio':None}
    (root/'results.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    lines=['# Complete standalone projected-syndrome and compiled-schedule comparison','',f"{result['cells']} fixed cells, {result['observations']} observations; all completed and verified: {all_complete}.",'',
        'The current reference includes all 32 predecessor methods plus the declared quiet, full-word and single-stage controls. Every cold solve and validation cost is charged. Ratios below use matched repetitions; censored costs remain null.','',
        '| Candidate | Split | Variables | Family | Reference | Median ratio | 95% interval | Pass |','|---|---|---:|---|---|---:|---|---|']
    for d in details:
        lines.append(f"| {d['candidate']} | {d['split']} | {d['n']} | {d['family']} | {d['reference']} | {d['paired_ratio_median']} | {d['ci95_paired_median']} | {d['pass']} |")
    lines+=['','All gates: `'+json.dumps(gates,sort_keys=True)+'`','',
        'These are bounded generated Boolean solves. No calibrated-operation, production index-calculus, asymptotic or rho claim is established.']
    (root/'RESULT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'cells':result['cells'],'observations':result['observations'],'gates':gates}))

if __name__=='__main__':main(Path(sys.argv[1]))
