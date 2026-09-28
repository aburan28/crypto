#!/usr/bin/env python3
"""Compare saved processes without discarding incomplete runs."""
import json,pathlib,sys,statistics,itertools,math
folder=pathlib.Path(sys.argv[1]);contract=json.loads((folder/'contract.json').read_text())
processes={r['name']:r for r in json.loads((folder/'processes.json').read_text())}
def read(name):
    p=folder/(name+'.stdout')
    if not p.exists():return []
    text=p.read_text().strip()
    if not text:return []
    return [json.loads(text)] if name.startswith('dlp-') else [json.loads(x) for x in text.splitlines() if x.strip()]
def complete(name):
    r=processes.get(name,{})
    return r.get('status')=='finished' and r.get('returncode')==0
rows=[];stage_correct=True
for n,ell,m in contract['stages']['cases']:
    for variant in contract['stages']['variants']:
        samples=[];ok=True;checked=0;statuses=[]
        for rep in range(3):
            name=f'stage-n{n}-l{ell}-m{m}-r{rep}-{variant}'
            baseline=read(f'stage-n{n}-l{ell}-m{m}-r{rep}-baseline')
            raw=read(name);statuses.append(processes.get(name,{}).get('status','missing'))
            base={x['target']:x for x in baseline if x['phase']=='target'}
            targets=[x for x in raw if x['phase']=='target'];checked+=len(targets)
            good=all(x['correct'] and x['input_hash']==base.get(x['target'],{}).get('input_hash') and x['solutions']==base[x['target']]['solutions'] for x in targets)
            expected=(1<<n) if n<=5 else 32
            ok &= good and complete(name) and len(targets)==expected and len({x['target'] for x in targets})==expected
            if not good:stage_correct=False
            setup=next((x for x in raw if x['phase']=='setup'),None)
            if complete(name) and setup and len(targets)==expected:
                setup_ns=raw[0]['field_setup_ns']+setup['template_setup_ns']+setup['parameter_basis_setup_ns']
                samples.append({'encoding_ns':setup['template_setup_ns']+sum(x['instantiate_ns'] for x in targets),
                    'cold_f4_ns':setup_ns+sum(x['instantiate_ns']+x['solve_ns'] for x in targets),
                    'cold_sat_ns':setup_ns+sum(x['instantiate_ns']+x['sat_ns'] for x in targets),
                    'online_f4_ns':sum(x['instantiate_ns']+x['solve_ns'] for x in targets),
                    'setup_ns':setup_ns,'template_bytes':setup['template_bytes'],'parameter_basis_bytes':setup['parameter_basis_bytes'],
                    'root_xors':sum(x['f4_root_elimination_word_xors'] or 0 for x in targets),'sat_conflicts':sum(x['sat_conflicts'] for x in targets),
                    'targets':expected,'holdouts':sum(x['split']=='holdout' for x in targets)})
        row={'n':n,'ell':ell,'m':m,'variant':variant,'complete_and_correct':bool(ok),'checked_target_repetitions':checked,'process_statuses':statuses}
        row.update({k:statistics.median(x[k] for x in samples) for k in samples[0]} if len(samples)==3 else {})
        rows.append(row)
for row in rows:
    baseline=next(x for x in rows if (x['n'],x['ell'],x['m'],x['variant'])==(row['n'],row['ell'],row['m'],'baseline'))
    row['cold_f4_candidate_over_baseline']=row.get('cold_f4_ns',0)/baseline['cold_f4_ns'] if row['complete_and_correct'] and baseline.get('cold_f4_ns') else None
    row['cold_sat_candidate_over_baseline']=row.get('cold_sat_ns',0)/baseline['cold_sat_ns'] if row['complete_and_correct'] and baseline.get('cold_sat_ns') else None
    row['total_calibrated_operations']=None;row['S']=None;row['rho_ratio']=None;row['floor_ratio']=None
    if row['complete_and_correct'] and baseline.get('online_f4_ns'):
        saving=(baseline['online_f4_ns']-row['online_f4_ns'])/row['targets']
        row['estimated_f4_break_even_targets']=math.ceil(max(0,row['setup_ns']-baseline['setup_ns'])/saving) if saving>0 else None

dlp=[];matches=True;redis_errors=0;warm_hits=0
for n,solver,seed in itertools.product(contract['dlp']['degrees'],contract['dlp']['solvers'],contract['dlp']['seeds']):
    prefix=f'dlp-n{n}-{solver}-s{seed}';reference=read(prefix+'-reference');reference=reference[0] if reference else {}
    for variant in ['reference',*contract['dlp']['variants']]:
        name=prefix+'-'+variant;raw=read(name);r=raw[0] if raw else {};counts=r.get('counts',{})
        proc=processes.get(name,{});refproc=processes.get(prefix+'-reference',{})
        match=bool(reference) and proc.get('status')=='finished' and proc.get('returncode')==refproc.get('returncode') and r.get('status')==reference.get('status') and r.get('result')==reference.get('result') and all(counts.get(k)==v for k,v in reference.get('counts',{}).items())
        matches &= match
        stats=counts.get('algebra_cache_current_thread',[])
        errs=sum(x['errors']+x['invalid'] for x in stats);redis_errors+=errs
        hits=sum(x['redis_hits'] for x in stats)
        if variant=='both-redis-warm':warm_hits+=hits
        dlp.append({'n':n,'solver':solver,'seed':seed,'variant':variant,'matches_reference':match,'status':r.get('status','missing'),
            'verified':r.get('result',{}).get('verified',False),'elapsed_seconds':r.get('elapsed_seconds'),'cache':stats,'redis_hits':hits,'cache_errors':errs})
summary={'stage_rows':rows,'stage_observed_answers_correct':stage_correct,'dlp_rows':dlp,'all_dlp_outcomes_match_reference':bool(matches),
    'redis_errors':redis_errors,'new_process_warm_redis_hits':warm_hits,'expected_processes':159,'recorded_processes':len(processes),
    'limitations':['Timing medians are diagnostics on a shared host, not calibrated total operations or an end-to-end speedup claim.','Incomplete parameter runs remain in the table.','SAT exact answers are not cached; the template is shared by SAT and Groebner encoders.']}
(folder/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ['stage_rows','dlp_rows']},indent=2))

admitted = len(processes)==159 and stage_correct and matches and redis_errors==0 and warm_hits>0 and all(r['complete_and_correct'] for r in rows if r['variant']!='parameterized')
if not admitted:raise SystemExit('Comparison failed: missing baseline/candidate evidence, mismatched answers, or unhealthy Redis')
