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
    expected={f'n{n}-{split}-{seed}-{family}' for n in protocol['variables'] for split in protocol.get('splits',['discovery','holdout']) for seed in protocol[split+'_seeds'] for family in protocol['families']}
    receipts=read(root/'receipts.json');assert len(receipts)==len(expected)
    receipts={r['cell']:r for r in receipts};assert set(receipts)==expected
    raw=defaultdict(list)
    for n in protocol['variables']:
        for line in (root/f'raw-n{n}.jsonl').read_text().splitlines():
            row=json.loads(line);assert json.loads(row['raw_line'])=={k:v for k,v in row.items() if k not in ['cell','split','raw_line']}
            raw[row['cell']].append(row)
    assert set(raw)==expected
    summaries=[];groups=defaultdict(list);all_complete=True
    logical=['outcome','model','reason','nodes','kernel_calls','decisions','forced','specialized_terms','source_rows','source_columns','max_depth','affine_eliminated','derived_rows','enumeration_points','enumeration_batches','enumeration_leaves','trace']
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
        kernel_signature={}
        for sample in samples:
            arm=sample['variant'];status=sample['outcome'];model=sample['model']
            assert 0<=sample['kernel_ns']<=sample['solve_ns']<=sample['total_ns']
            assert sample['solve_ns']+sample['validation_ns']<=sample['total_ns']
            for field,enabled in [('substitution_ns',arm.startswith('affine_')),('enumeration_ns',arm.startswith('gray_') or arm.startswith('packed_gray'))]:
                assert isinstance(sample[field],int) and sample[field]>=0 if enabled else sample[field] is None
            assert sample['kernel_ns']+sum(sample[k] or 0 for k in ['substitution_ns','enumeration_ns'])<=sample['solve_ns']
            assert sample['enumeration_points']<=1<<n
            assert sample['enumeration_batches']<=sample['enumeration_points']<=16*sample['enumeration_batches']
            if not arm.startswith(('gray_','packed_gray')):assert sample['enumeration_points']==sample['enumeration_batches']==sample['enumeration_leaves']==0
            if not arm.startswith('affine_'):assert sample['affine_eliminated']==sample['derived_rows']==0
            assert sample['affine_eliminated']<=n*sample['nodes']
            if arm=='packed_untraced':assert sample['trace'] is None
            else:assert isinstance(sample['trace'],int)
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
                assert sample['reason'] in ['NODE_CAP','ROW_CAP','COLUMN_CAP','ENUM_CAP']
            if sample['verified'] and fixture['search_reference']!='UNKNOWN':assert status==fixture['search_reference']
            if status=='UNKNOWN' or not sample['verified']:all_complete=False
            if arm!='search' or protocol.get('enable_merge'):
                signature=tuple(sample[k] for k in logical)
                group='untraced' if arm=='packed_untraced' else 'gray' if arm in ('gray_scalar','gray_simd') else 'leaf12' if arm.startswith('packed_gray12') else 'leaf16' if arm.startswith('packed_gray16') else 'affine_sl_basis' if arm in ('affine_sl_basis_list','affine_sl_basis_fast') else 'basis' if arm in ('basis_list','basis_wide') else 'tail' if arm in ('tail_list','tail_wide') else 'search' if arm in ('search','merge_search','quadratic_state','packed_state') else 'degree2' if arm in ('small_flat','word_tail') else 'degree3'
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
                       'kernel_fraction':statistics.median(s['kernel_ns']/s['solve_ns'] for s in values),
                       **{k:values[0][k] for k in logical if k not in ['outcome','reason']},
                       'kernel_calls_by_active':values[0]['kernel_calls_by_active'],'flat_calls':values[0]['flat_calls']}
        summaries.append({'cell':cell,'n':n,'family':fixture['family'],'seed':fixture['seed'],'split':fixture['split'],'arms':arms,'peak_rss_bytes':receipt['peak_rss_bytes']})
    candidates=protocol.get('candidate_arms',['search','bucket','hybrid'])
    gates=[]
    for n in protocol['variables']:
        if n==12:continue
        for family in protocol['families']:
            cells=[c for c in summaries if c['n']==n and c['family']==family and c['split']=='holdout']
            for candidate in candidates:
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
            'full_solve_gates':{arm:'PASS' if all_complete and all(g['pass'] for g in gates if g['candidate']==arm) else 'REJECTED' for arm in candidates},
            'hybrid_vs_fastest_control_gate':'PASS' if all_complete and all(g['pass'] for g in strongest) else 'REJECTED',
            'gate_details':gates,'strongest_control_details':strongest,'cells_detail':summaries,
            'classification':protocol['classification'],'production_solver_cost':None,'full_ic_cost':None,'rho_ratio':None}
    if protocol.get('enable_word'):
        word_gates=[]
        for n in protocol['variables']:
            if n==12:continue
            for family in protocol['families']:
                cells=[c for c in summaries if c['n']==n and c['family']==family and c['split']=='holdout']
                for label,controls,threshold in [('best_prior',['search','flat','bucket','hybrid'],2.0),('same_policy',['small_flat'],1.05)]:
                    ratios=[]
                    for cell in cells:
                        pairs={(s['rep'],s['variant']):s for s in raw[cell['cell']][1:]}
                        for rep in range(protocol['repetitions']):
                            values=[pairs[rep,a] for a in controls+['word_tail']]
                            if all(v['outcome']!='UNKNOWN' and v['verified'] for v in values):
                                ratios.append(min(pairs[rep,a]['total_ns'] for a in controls)/pairs[rep,'word_tail']['total_ns'])
                    complete=len(ratios)==len(cells)*protocol['repetitions'];ci=bootstrap(ratios) if complete else None
                    word_gates.append({'n':n,'family':family,'control':label,'paired_ratio_median':statistics.median(ratios) if complete else None,
                                       'ci95_paired_median':ci,'threshold':threshold,'pass':complete and ci[0]>threshold,'complete':complete})
        result['word_solve_details']=word_gates
        result['word_solve_gate']='PASS' if all_complete and all(g['pass'] for g in word_gates) else 'REJECTED'
    if protocol.get('enable_merge'):
        merge_gates=[]
        for n in protocol['variables']:
            if n==12:continue
            for family in protocol['families']:
                cells=[c for c in summaries if c['n']==n and c['family']==family and c['split']=='holdout'];ratios=[]
                for cell in cells:
                    pairs={(s['rep'],s['variant']):s for s in raw[cell['cell']][1:]}
                    for rep in range(protocol['repetitions']):
                        values=[pairs[rep,a] for a in protocol['variants']]
                        if all(v['outcome']!='UNKNOWN' and v['verified'] for v in values):
                            ratios.append(min(pairs[rep,a]['total_ns'] for a in ['search','flat','bucket','hybrid','small_flat','word_tail'])/pairs[rep,'merge_search']['total_ns'])
                complete=len(ratios)==len(cells)*protocol['repetitions'];ci=bootstrap(ratios) if complete else None
                merge_gates.append({'n':n,'family':family,'paired_ratio_median':statistics.median(ratios) if complete else None,
                                    'ci95_paired_median':ci,'threshold':2.0,'pass':complete and ci[0]>2.0,'complete':complete})
        result['merge_search_details']=merge_gates
        result['merge_search_gate']='PASS' if all_complete and all(g['pass'] for g in merge_gates) else 'REJECTED'
    if protocol.get('enable_quadratic'):
        quadratic_gates=[]
        controls=['search','flat','bucket','hybrid','small_flat','word_tail','merge_search']
        assert len(controls)==7 and 'merge_search' in controls
        for n in protocol['variables']:
            if n==12:continue
            for family in protocol['families']:
                cells=[c for c in summaries if c['n']==n and c['family']==family and c['split']=='holdout'];ratios=[]
                for cell in cells:
                    pairs={(s['rep'],s['variant']):s for s in raw[cell['cell']][1:]}
                    for rep in range(protocol['repetitions']):
                        values=[pairs[rep,a] for a in protocol['variants']]
                        if all(v['outcome']!='UNKNOWN' and v['verified'] for v in values):
                            ratios.append(min(pairs[rep,a]['total_ns'] for a in controls)/pairs[rep,'quadratic_state']['total_ns'])
                complete=len(ratios)==len(cells)*protocol['repetitions'];ci=bootstrap(ratios) if complete else None
                quadratic_gates.append({'n':n,'family':family,'controls':controls,'paired_ratio_median':statistics.median(ratios) if complete else None,
                                        'ci95_paired_median':ci,'threshold':2.0,'pass':complete and ci[0]>2.0,'complete':complete})
        result['quadratic_state_details']=quadratic_gates
        result['quadratic_state_gate']='PASS' if all_complete and all(g['pass'] for g in quadratic_gates) else 'REJECTED'
    if protocol.get('enable_packed'):
        packed_gates=[]
        for split in ['regression','holdout']:
            for n in protocol['variables']:
                if n==12:continue
                for family in protocol['families']:
                    cells=[c for c in summaries if c['n']==n and c['family']==family and c['split']==split]
                    for label,controls in [('current_frontier',['search','flat','bucket','hybrid','small_flat','word_tail','merge_search','quadratic_state']),('historical_seven',['search','flat','bucket','hybrid','small_flat','word_tail','merge_search'])]:
                        ratios=[]
                        for cell in cells:
                            pairs={(s['rep'],s['variant']):s for s in raw[cell['cell']][1:]}
                            for rep in range(protocol['repetitions']):
                                values=[pairs[rep,a] for a in controls+['packed_state']]
                                if all(v['outcome']!='UNKNOWN' and v['verified'] for v in values):
                                    ratios.append(min(pairs[rep,a]['total_ns'] for a in controls)/pairs[rep,'packed_state']['total_ns'])
                        complete=len(ratios)==len(cells)*protocol['repetitions'];ci=bootstrap(ratios) if complete else None
                        packed_gates.append({'split':split,'n':n,'family':family,'reference':label,'controls':controls,
                                             'paired_ratio_median':statistics.median(ratios) if complete else None,'ci95_paired_median':ci,
                                             'threshold':2.0,'pass':complete and ci[0]>2.0,'improvement_pass':complete and ci[0]>1.0,'complete':complete})
        result['packed_details']=packed_gates
        for reference in ['current_frontier','historical_seven']:
            selected=[g for g in packed_gates if g['reference']==reference]
            result['packed_'+reference+'_gate']='PASS' if all_complete and all(g['pass'] for g in selected) else 'REJECTED'
        result['packed_frontier_improvement_gate']='PASS' if all_complete and all(g['improvement_pass'] for g in packed_gates if g['reference']=='current_frontier') else 'REJECTED'
    if protocol.get('enable_basis'):
        basis_gates=[]
        for candidate,matched in [('basis_wide','basis_list'),('tail_wide','tail_list')]:
            for split in ['regression','holdout']:
                for n in protocol['variables']:
                    if n==12:continue
                    for family in protocol['families']:
                        cells=[c for c in summaries if c['n']==n and c['family']==family and c['split']==split]
                        for label,controls,threshold in [('strongest',[a for a in protocol['variants'] if a!=candidate],2.0),('same_policy',[matched],1.05)]:
                            ratios=[]
                            for cell in cells:
                                pairs={(s['rep'],s['variant']):s for s in raw[cell['cell']][1:]}
                                for rep in range(protocol['repetitions']):
                                    values=[pairs[rep,a] for a in controls+[candidate]]
                                    if all(v['outcome']!='UNKNOWN' and v['verified'] for v in values):
                                        ratios.append(min(pairs[rep,a]['total_ns'] for a in controls)/pairs[rep,candidate]['total_ns'])
                            complete=len(ratios)==len(cells)*protocol['repetitions'];ci=bootstrap(ratios) if complete else None
                            basis_gates.append({'candidate':candidate,'split':split,'n':n,'family':family,'reference':label,'controls':controls,
                                                'paired_ratio_median':statistics.median(ratios) if complete else None,'ci95_paired_median':ci,
                                                'threshold':threshold,'pass':complete and ci[0]>threshold,'improvement_pass':complete and ci[0]>1.0,'complete':complete})
        result['basis_details']=basis_gates
        for candidate in ['basis_wide','tail_wide']:
            strongest=[g for g in basis_gates if g['candidate']==candidate and g['reference']=='strongest']
            matched=[g for g in basis_gates if g['candidate']==candidate and g['reference']=='same_policy']
            result[candidate+'_dramatic_gate']='PASS' if all_complete and all(g['pass'] for g in strongest) else 'REJECTED'
            result[candidate+'_improvement_gate']='PASS' if all_complete and all(g['improvement_pass'] for g in strongest) else 'REJECTED'
            result[candidate+'_backend_gate']='PASS' if all_complete and all(g['pass'] for g in matched) else 'REJECTED'
    if protocol.get('final_selection'):
        new_gates=[]
        for candidate in protocol['new_candidates']:
            for split in ['regression','holdout']:
                for n in protocol['variables']:
                    if n==12:continue
                    for family in protocol['families']:
                        cells=[c for c in summaries if c['n']==n and c['family']==family and c['split']==split]
                        references=[('retained_frontier',protocol['reference_arms'],2.0),('same_policy',[protocol['matched_controls'][candidate]],1.05)]
                        for label,controls,threshold in references:
                            controls=[a for a in controls if a!=candidate];ratios=[]
                            for cell in cells:
                                pairs={(s['rep'],s['variant']):s for s in raw[cell['cell']][1:]}
                                for rep in range(protocol['repetitions']):
                                    values=[pairs[rep,a] for a in controls+[candidate]]
                                    if all(v['outcome']!='UNKNOWN' and v['verified'] for v in values):
                                        ratios.append(min(pairs[rep,a]['total_ns'] for a in controls)/pairs[rep,candidate]['total_ns'])
                            complete=len(ratios)==len(cells)*protocol['repetitions'];ci=bootstrap(ratios) if complete else None
                            new_gates.append({'candidate':candidate,'split':split,'n':n,'family':family,'reference':label,'controls':controls,
                                              'paired_ratio_median':statistics.median(ratios) if complete else None,'ci95_paired_median':ci,
                                              'threshold':threshold,'pass':complete and ci[0]>threshold,'complete':complete})
        result['new_details']=new_gates
        result['new_gates']={candidate:{reference:'PASS' if all_complete and all(g['pass'] for g in new_gates if g['candidate']==candidate and g['reference']==reference) else 'REJECTED' for reference in ['retained_frontier','same_policy']} for candidate in protocol['new_candidates']}
    (root/'results.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    lines=['# Complete generic Boolean solve comparison','',f"Evidence integrity: **PASS** across {result['cells']} cells and {result['samples']} full-solve observations.",'',
           f"All results completed and independently verified: **{all_complete}**. Full-solve gates against flat: **{result['full_solve_gates']}**. Hybrid versus fastest control: **{result['hybrid_vs_fastest_control_gate']}**.",'',
           'The table reports cold solve plus result-validation milliseconds over two holdout seeds and the declared rotated repetitions. A completion cost is null if any sample in the group is censored or lacks the required verification. Observed capped work is retained separately. Ratios use pooled medians; acceptance uses paired intervals.','',
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
    lines+=['','Kernel-backed solvers have equal outcomes/models, logical counters and trace digests within each matched policy group. The search-only control may explore a different tree. SAT models are checked against original equations; UNSAT verification requires a completed independent search reference. UNKNOWN is censored, not UNSAT.','',
            'This driver solves bounded generated Boolean systems. It does not benchmark the repository inherited-F4 implementation, accept curve targets, recover scalars or establish index-calculus performance. Production and cryptanalytic costs stay null.']
    if not protocol.get('enable_word'):
        lines=[line.replace('Kernel-backed solvers have equal outcomes/models, logical counters and trace digests within each matched policy group.','Kernel-backed solvers have equal outcomes/models, logical counters and trace digests on every matched cell.') for line in lines]
        lines=[line.replace("over two holdout seeds and the declared rotated repetitions.",'over two holdout seeds and four rotated repetitions.') for line in lines]
    else:
        lines += ['',f"Selective one-word degree-2 full-solve gate: **{result['word_solve_gate']}**. The small_flat and word_tail methods share a policy and must match exact traces/counters. They are compared separately with the always-degree-3 methods and search-only."]
    if protocol.get('enable_merge'):
        lines+=['',f"Ordered-specialization complete-solve gate: **{result['merge_search_gate']}**. The search and merge_search arms have identical outcomes/models, logical work and trace digests. Only representation work in specialization changes. The gate compares merge_search with the fastest of all six prior full-solve methods."]
    if protocol.get('enable_quadratic'):
        lines+=['',f"Fixed-quadratic complete-solve gate: **{result['quadratic_state_gate']}**. Compilation is charged to each cold solve. Search, merge_search and quadratic_state must match outcomes/models, logical counters and trace digests. The reference is the pointwise fastest of all seven prior methods.",'',
                '| Variables | Family | Fastest prior / quadratic | 95% paired interval | Gate |','|---:|---|---:|---|---|']
        for gate in result['quadratic_state_details']:
            ratio='null' if gate['paired_ratio_median'] is None else f"{gate['paired_ratio_median']:.3f}"
            interval='null' if gate['ci95_paired_median'] is None else '–'.join(f'{v:.3f}' for v in gate['ci95_paired_median'])
            lines.append(f"| {gate['n']} | {gate['family']} | {ratio} | {interval} | {'PASS' if gate['pass'] else 'REJECTED'} |")
        lines+=['','These paired intervals describe repeated timings on two fixed holdout fixtures per size/family, not a confidence interval over a population of Boolean systems. Fresh-worker RSS includes all arms and is not candidate-specific memory.']
    if protocol.get('enable_packed'):
        lines+=['',f"Packed-state current-frontier 2x gate: **{result['packed_current_frontier_gate']}**. Incremental improvement gate (>1x lower bound): **{result['packed_frontier_improvement_gate']}**. Cumulative historical-seven 2x gate: **{result['packed_historical_seven_gate']}**.",'',
                'Both fresh holdouts and the entire prior confirmation regression grid must pass. The current frontier includes the already faster quadratic_state method; the historical-seven comparison cannot substitute for it.', '',
                '| Split | Variables | Family | Reference | Paired ratio | 95% interval | >2x gate |','|---|---:|---|---|---:|---|---|']
        for gate in result['packed_details']:
            ratio='null' if gate['paired_ratio_median'] is None else f"{gate['paired_ratio_median']:.3f}"
            interval='null' if gate['ci95_paired_median'] is None else '–'.join(f'{v:.3f}' for v in gate['ci95_paired_median'])
            lines.append(f"| {gate['split']} | {gate['n']} | {gate['family']} | {gate['reference']} | {ratio} | {interval} | {'PASS' if gate['pass'] else 'REJECTED'} |")
    if protocol.get('enable_basis'):
        lines+=['','## Transported coefficient-row-space comparison','',
                'basis_list/basis_wide share a full-RREF policy and branch on that basis. tail_list/tail_wide keep quadratic rows in echelon form, reduce only the affine tail, and branch using carried original-equation residuals. The backends in each pair must match every logical counter, model and trace. New inference policies need not match the older search tree.', '',
                'Source-row and distinct-column counts price each reduction input, including zero rows. Specialization counts cover every carried representation. Basis histograms use nominal unassigned-variable width. All costs, including the extra original-equation state in the tail policy, are inside cold solve timing.', '',
                '| Candidate | Dramatic >2x strongest-reference gate | Any improvement strongest-reference gate | Matched backend >1.05x gate |','|---|---|---|---|']
        for candidate in ['basis_wide','tail_wide']:
            lines.append(f"| {candidate} | {result[candidate+'_dramatic_gate']} | {result[candidate+'_improvement_gate']} | {result[candidate+'_backend_gate']} |")
        lines+=['','| Candidate | Split | Variables | Family | Reference | Paired ratio | 95% interval | Gate |','|---|---|---:|---|---|---:|---|---|']
        for g in result['basis_details']:
            ratio='null' if g['paired_ratio_median'] is None else f"{g['paired_ratio_median']:.3f}"
            interval='null' if g['ci95_paired_median'] is None else '–'.join(f'{v:.3f}' for v in g['ci95_paired_median'])
            lines.append(f"| {g['candidate']} | {g['split']} | {g['n']} | {g['family']} | {g['reference']} | {ratio} | {interval} | {'PASS' if g['pass'] else 'REJECTED'} |")
    if protocol.get('final_selection'):
        lines+=['','## Frozen finalist comparison','',
                'Reference arms include all thirteen previously retained methods, the packed method without diagnostic hashing, the affine finalist pair, and all scalar enumeration controls. The three SIMD treatments are compared with this fixed reference roster. Comparisons among the SIMD treatments do not redefine that preregistered baseline. The affine finalist is also shown; its matched list implementation is a backend control.', '',
                'Enumeration uses assignment and block counters, not search-node counts. Its checksum folds each block index and the wrapping sum of all sixteen equation-syndrome words; it is a diagnostic, not a certificate or a canonical search-tree trace. Scalar and SIMD modes must match models, all logical counters and checksums. Tree hybrids also preserve their prefix traces. The untraced packed control must preserve models and logical counters; its absent trace is null.', '',
                '| Candidate | >2x retained-frontier gate | >1.05x matched-backend gate |','|---|---|---|']
        for candidate,verdict in result['new_gates'].items():lines.append(f"| {candidate} | {verdict['retained_frontier']} | {verdict['same_policy']} |")
        lines+=['','| Candidate | Split | Variables | Family | Reference | Paired ratio | 95% interval | Gate |','|---|---|---:|---|---|---:|---|---|']
        for g in result['new_details']:
            ratio='null' if g['paired_ratio_median'] is None else f"{g['paired_ratio_median']:.3f}"
            interval='null' if g['ci95_paired_median'] is None else '–'.join(f'{v:.3f}' for v in g['ci95_paired_median'])
            lines.append(f"| {g['candidate']} | {g['split']} | {g['n']} | {g['family']} | {g['reference']} | {ratio} | {interval} | {'PASS' if g['pass'] else 'REJECTED'} |")
    (root/'RESULT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:result[k] for k in ['correctness','cells','samples','all_results_completed_and_verified','full_solve_gates','hybrid_vs_fastest_control_gate','quadratic_state_gate','packed_current_frontier_gate','packed_frontier_improvement_gate','packed_historical_seven_gate','basis_wide_dramatic_gate','basis_wide_improvement_gate','basis_wide_backend_gate','tail_wide_dramatic_gate','tail_wide_improvement_gate','tail_wide_backend_gate','new_gates'] if k in result},sort_keys=True))

if __name__=='__main__':main(Path(sys.argv[1]).resolve())
