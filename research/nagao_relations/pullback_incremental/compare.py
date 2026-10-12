"""Compare every iteration with frozen evidence and matched fresh holdouts."""
import argparse
import hashlib
import json
from pathlib import Path

import run as campaign


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def sumCounts(rows, key):
    return {name:sum(row[key][name] for row in rows) for name in rows[0][key]}


def summarize(folder):
    rows = campaign.readRows(folder/'raw.jsonl')
    require(not any(r['kind'] in ('error','validation_error') for r in rows), 'campaign contains errors')
    provenance = json.loads((folder/'provenance.json').read_text())
    contract = provenance['contract']
    stage = [r for r in rows if r['kind']=='stage']
    dlp = [r for r in rows if r['kind']=='e2e']
    require(len(stage)==len(provenance['targets']['targets'])*len(contract['modes'])*len(contract['variants']), 'stage run count mismatch')
    require(len(dlp)==len(contract['full_dlp']['panels'])*len(contract['full_dlp']['seeds'])*len(contract['full_dlp']['variants']), 'DLP run count mismatch')
    require(all(r['status']=='verified' for r in dlp), 'incomplete scalar recovery')
    key = lambda r:(r['n'],r['d'],tuple(r['target']),r['mode'],r['variant'])
    require(len({key(r) for r in stage})==len(stage), 'duplicate stage cell')
    byKey = {key(r):r for r in stage}
    previousStage = [r for r in campaign.readRows(campaign.PRIOR/'results/raw.jsonl')
                     if r['kind']=='trial' and r['variant']=='coefficient-pullback' and r['status']!='timeout']
    replays = 0
    censoredReplays = []
    for original in previousStage:
        replay = byKey[key(original)]
        if replay['status']=='timeout':
            censoredReplays.append(key(original))
            continue
        require(replay['solutions']==original['solutions'], 'frozen baseline solutions drifted')
        require(replay['field_api_counts']==original['field_api_counts'], 'frozen baseline field counters drifted')
        replays += 1
    stageGroups = []
    for corpus in contract['corpora']:
        for panel in contract['panels']:
            for variant in contract['variants']:
                for mode in contract['modes']:
                    rr=[r for r in stage if r['corpus']==corpus and r['n']==panel['n'] and r['variant']==variant and r['mode']==mode]
                    field = None
                    if 'field_api_counts' in rr[0]:
                        field = {k:sum(r['field_api_counts']['totals'][k] for r in rr)
                                 for k in ('additions','multiplications','squarings','fieldOperations')}
                    stageGroups.append({'corpus':corpus,**panel,'variant':variant,'mode':mode,
                        'runs':len(rr),'resolved':sum(r['status']!='timeout' for r in rr),
                        'resolved_within_budget':sum(r['status']!='timeout' and r['all_phase_seconds']<=contract['seconds_per_instance'] for r in rr),
                        'verified_relations':sum(r['verified_unique_relations'] for r in rr),
                        'field_api_counts_including_censored':field,
                        'all_phase_seconds_including_censored':sum(r['all_phase_seconds'] for r in rr)})
    gates = []
    for corpus in contract['corpora']:
        subset = [r for r in stage if r['corpus']==corpus and r['n']==11 and r['mode']=='enumerate']
        matched = []
        for b in [r for r in subset if r['variant']=='coefficient-pullback']:
            n = byKey[(*key(b)[:-1],'normalized-pullback')]
            c = byKey[(*key(b)[:-1],'cached-pullback')]
            if any(r['status']!='complete' for r in (b,n,c)):
                continue
            require(b['solutions']==n['solutions']==c['solutions'], 'complete outputs changed')
            require(n['field_api_counts']==c['field_api_counts'], 'cache changed arithmetic accounting')
            matched.append((b,n,c))
        vectors = {}
        for i,variant in enumerate(('coefficient-pullback','normalized-pullback','cached-pullback')):
            vectors[variant] = {k:sum(p[i]['field_api_counts']['totals'][k] for p in matched)
                               for k in ('additions','multiplications','squarings','fieldOperations')}
        before = vectors['coefficient-pullback']['multiplications']
        after = vectors['normalized-pullback']['multiplications']
        reduction = 1-after/before if before else None
        gates.append({'corpus':corpus,'complete_equal_output_targets':len(matched),'field_api_counts':vectors,
                      'multiplication_reduction_fraction':reduction,
                      'diagnostic_gate_passed':len(matched)==4 and reduction>=.10})
    oracles = {};calls = 0;distinct = set();dlpReplays = 0;paths = []
    oldDlp = campaign.readRows(campaign.PRIOR/'e2e_results/raw.jsonl')
    dlpGroups = []
    holdoutTargets = []
    for row in dlp:
        n,d = row['panel']['n'],row['panel']['d']
        if (n,d) not in oracles:
            oracles[n,d]=campaign.previous.previous.PairOracle(n,d)
        oracle = oracles[n,d];f,curve = oracle.f,oracle.curve
        G = tuple(f.fromCoords(x) for x in row['generator'])
        Q = tuple(f.fromCoords(x) for x in row['target'])
        require(curve.mul(G,row['recovered_scalar'])==Q, 'independent scalar verification failed')
        for attempt in row['attempts']:
            r=attempt['result'];point=tuple(f.fromCoords(x) for x in r['target'])
            expected=oracle.expected(point);got={tuple(xs) for xs in r['solutions']}
            require(got<=expected, 'DLP oracle emitted invalid relation')
            if r['status']=='complete':
                require(got==expected,'DLP oracle lost a relation')
            calls+=1;distinct.add((n,d,tuple(r['target'])))
        if row['corpus']=='frozen' and row['variant']=='coefficient-pullback':
            original=next(r for r in oldDlp if r['panel']==row['panel'] and r['seed']==row['seed'] and r['variant']=='coefficient-pullback')
            require(original['target']==row['target'] and original['all_phase_field_api_counts']==row['all_phase_field_api_counts']
                    and original['scalar_modular_counts']==row['scalar_modular_counts'], 'cold baseline counters drifted')
            dlpReplays+=1
        if row['variant'] in ('normalized-pullback','cached-pullback'):
            before=next(r for r in dlp if r['panel']==row['panel'] and r['seed']==row['seed'] and r['variant']=='coefficient-pullback')
            require(row['target']==before['target'] and row['generator']==before['generator'],'paired DLP instance changed')
            path=lambda r:[(a['phase'],a['coefficient'],a['result']['target'],a['result']['solutions']) for a in r['attempts']]
            paths.append({'n':n,'seed':row['seed'],'variant':row['variant'],'identical_attempt_path':path(row)==path(before)})
    for corpus in contract['corpora']:
        for panel in contract['full_dlp']['panels']:
            for variant in contract['full_dlp']['variants']:
                rr=[r for r in dlp if r['corpus']==corpus and r['panel']==panel and r['variant']==variant]
                dlpGroups.append({'corpus':corpus,**panel,'variant':variant,'verified':len(rr),
                    'field_api_counts':sumCounts(rr,'all_phase_field_api_counts'),
                    'scalar_modular_counts':sumCounts(rr,'scalar_modular_counts'),
                    'phase_seconds':sumCounts(rr,'phase_seconds'),
                    'cold_seconds':sum(r['cold_wall_seconds'] for r in rr),
                    'attempt_count':sum(len(r['attempts']) for r in rr),
                    'full_dlp_S':None,'rho_ratio':None,'cost_over_floor':None})
    for panel in contract['full_dlp']['panels']:
        frozen = {tuple(r['target']) for r in dlp if r['panel']==panel and r['corpus']=='frozen'}
        fresh = {tuple(r['target']) for r in dlp if r['panel']==panel and r['corpus']=='fresh'}
        holdoutTargets.append({**panel,'distinct_frozen_targets':len(frozen),'distinct_fresh_seed_targets':len(fresh),
                               'target_overlap':len(frozen&fresh),'new_distinct_targets':len(fresh-frozen)})
    return {'stage_run_count':len(stage),'dlp_run_count':len(dlp),'stage_baseline_exact_counter_replays':replays,
            'stage_baseline_censored_replays':censoredReplays,
            'dlp_baseline_exact_counter_replays':dlpReplays,'stage_groups':stageGroups,'stage_diagnostic_gates':gates,
            'dlp_oracle_checks':calls,'dlp_distinct_oracle_targets':len(distinct),'dlp_attempt_paths':paths,'dlp_groups':dlpGroups,
            'dlp_holdout_target_distinctness':holdoutTargets,
            'classification':'Engineering candidate; normalized-cost classification remains pending.',
            'end_to_end_speedup_established':False}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=summarize(args.input)
    result['input_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                           for p in [args.input/'raw.jsonl',args.input/'provenance.json',Path(__file__)]}
    with args.output.open('x') as out:
        json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('stage_groups','dlp_groups','dlp_attempt_paths','input_sha256')},indent=2))
