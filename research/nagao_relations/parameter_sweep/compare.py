"""Audit a frozen sweep and emit component-wise comparisons, never calibrated S."""
import argparse
from collections import Counter, defaultdict
import hashlib
import itertools
import json
from pathlib import Path
import statistics

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
FIELD=('additions','multiplications','squarings')
SCALAR=('additions','multiplications','inversions')


def read(path):return json.loads(path.read_text())
def rows(path):return [json.loads(line) for line in path.read_text().splitlines()]
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def total(rs,kind,keys):
    return {key:sum(r[kind]['totals'][key] if kind=='field_api_counts' else r[kind][key] for r in rs) for key in keys}
def smaller(a,b,keys):return all(a[k]<=b[k] for k in keys) and any(a[k]<b[k] for k in keys)


def compare(directory,targetsPath,regressionPath):
    raw=rows(directory/'raw.jsonl');contract=read(HERE/'contract.json');targets=read(targetsPath)
    provenance=read(directory/'provenance.json');summary=read(directory/'summary.json')
    assert not summary['errors'] and not [r for r in raw if r['kind']=='error']
    assert targets['contract_sha256']==digest(HERE/'contract.json')
    for path,sha in provenance['source_sha256'].items():
        assert digest(ROOT/path)==sha,('changed measured source',path)
    census=[r for r in raw if r['kind']=='census'];stages=[r for r in raw if r['kind'] in ('stage','algebraic')]
    dlps=[r for r in raw if r['kind']=='dlp'];rhos=[r for r in raw if r['kind']=='rho']
    stageKeys=set();expectedStage=set();expectedDlp=set();expectedCensus=set()
    for p in contract['profiles']:
        selected=next(t['targets'] for t in targets['profiles'] if t['profile']==p['name'])
        assert len({tuple(t['target']) for t in selected})==24
        for ell,k in itertools.product(p['dimensions'],contract['summands']):
            expectedCensus.add((p['name'],ell,k))
            variants=contract['stage_variants']+(['cached-pullback','s4-symmetric','chained-s3'] if p['q_bits']==1 and k==3 else [])
            expectedStage.update((p['name'],ell,k,v,t['stratum'],tuple(t['target'])) for v in variants for t in selected)
            expectedDlp.update((p['name'],ell,k,v,s) for v in contract['e2e']['variants'] for s in contract['e2e']['seeds'])
    for r in stages:
        key=(r['profile'],r['ell'],r['k'],r['variant'],r['stratum'],tuple(r['target']))
        assert key not in stageKeys;stageKeys.add(key)
    assert stageKeys==expectedStage
    assert len(dlps)==len(expectedDlp) and {(r['profile'],r['ell'],r['k'],r['variant'],r['seed']) for r in dlps}==expectedDlp
    assert len(census)==len(expectedCensus) and {(r['profile'],r['ell'],r['k']) for r in census}==expectedCensus
    assert len(rhos)==8 and all(r['status']=='verified' for r in rhos)
    assert {(r['profile'],r['seed']) for r in rhos}=={(p['name'],s) for p in contract['profiles'] for s in contract['e2e']['seeds']}
    for r in stages+dlps+rhos:
        assert r['full_dlp_S'] is None and r['rho_ratio'] is None
        report=r['field_api_counts']
        if report is not None:
            for key in FIELD:
                assert report['totals'][key]==sum(phase[key] for phase in report['phases'].values())
    for r in dlps+rhos:
        assert abs(sum(r['phase_seconds'].values())-r['cold_seconds'])<1e-8
        if r['status']=='verified':assert r['recovered_scalar']==r['planted_scalar']
        ref=next(x for x in rhos if x['profile']==r['profile'] and x['seed']==r['seed'])
        assert (r['target'],r['generator'],r['subgroup_order'])==(ref['target'],ref['generator'],ref['subgroup_order'])

    groups=defaultdict(list)
    for r in stages:groups[(r['profile'],r['ell'],r['k'],r['variant'],r['stratum'])].append(r)
    stageSummary=[]
    for (p,e,k,v,s),rs in sorted(groups.items()):
        stageSummary.append({'profile':p,'ell':e,'k':k,'variant':v,'stratum':s,'attempted':len(rs),
            'statuses':dict(Counter(r['status'] for r in rs)),'within_budget':sum(r['within_budget'] for r in rs),
            'cold_seconds_sum':sum(r['cold_seconds'] for r in rs),'cold_seconds_median':statistics.median(r['cold_seconds'] for r in rs),
            'field_api_totals':total(rs,'field_api_counts',FIELD) if all(r['field_api_counts'] for r in rs) else None})
    # Pair exact same target, not just equal counts of successes.
    paired=[]
    for c in census:
        for stratum in ('development','holdout'):
            rs={v:{tuple(r['target']):r for r in groups[(c['profile'],c['ell'],c['k'],v,stratum)]}
                for v in contract['stage_variants']}
            baseline=rs['enumerate-last'];candidate=rs['pair-mitm'];assert baseline.keys()==candidate.keys()
            complete=[t for t in baseline if baseline[t]['status']!='timeout' and candidate[t]['status']!='timeout']
            assert all(baseline[t]['status']==candidate[t]['status'] for t in complete)
            sums={v:total([rs[v][t] for t in complete],'field_api_counts',FIELD) for v in rs}
            paired.append({'profile':c['profile'],'ell':c['ell'],'k':c['k'],'stratum':stratum,
                'equal_completed_targets':len(complete),'censored_targets':len(baseline)-len(complete),
                'field_api_totals':sums,'pair_mitm_field_vector_improves':smaller(sums['pair-mitm'],sums['enumerate-last'],FIELD),
                'calibrated_speedup':None})

    dlpGroups=defaultdict(list)
    for r in dlps+rhos:dlpGroups[(r['profile'],r['ell'],r['k'],r['variant'])].append(r)
    dlpSummary=[]
    # None-valued rho parameters are sorted with a string key.
    for (p,e,k,v),rs in sorted(dlpGroups.items(),key=lambda item:str(item[0])):
        complete=[r for r in rs if r['status']=='verified']
        dlpSummary.append({'profile':p,'ell':e,'k':k,'variant':v,'attempted':len(rs),'verified':len(complete),
            'errors':dict(Counter(r.get('error') for r in rs if r['status']!='verified')),
            'field_api_all_attempts':total(rs,'field_api_counts',FIELD),'scalar_all_attempts':total(rs,'scalar_modular_counts',SCALAR),
            'cold_seconds_all_attempts':sum(r['cold_seconds'] for r in rs),
            'field_api_verified_only':total(complete,'field_api_counts',FIELD),'scalar_verified_only':total(complete,'scalar_modular_counts',SCALAR),
            'phase_seconds_all_attempts':{phase:sum(r['phase_seconds'][phase] for r in rs) for phase in rs[0]['phase_seconds']},
            'attempt_statuses':dict(Counter(a['status'] for r in rs for a in r['attempts'])),
            'columns':[r.get('columns') for r in rs],'ranks':[r.get('rank') for r in rs],
            'full_dlp_S':None,'rho_ratio':None,'cost_over_floor':None})
    # Freeze-gate test: BOTH individual seeds, not just an aggregate or a selected seed.
    changes=[];variantPairs=[]
    def gate(a,b):
        a={r['seed']:r for r in a};b={r['seed']:r for r in b}
        complete=all(a[s]['status']==b[s]['status']=='verified' for s in contract['e2e']['seeds'])
        return {'both_seeds_verified':complete,'field_vector_gate':complete and all(smaller(b[s]['field_api_counts']['totals'],a[s]['field_api_counts']['totals'],FIELD) for s in a),
            'scalar_nonincrease':complete and all(all(b[s]['scalar_modular_counts'][x]<=a[s]['scalar_modular_counts'][x] for x in SCALAR) for s in a),
            'calibrated_speedup':None}
    for p in contract['profiles']:
        for v in contract['e2e']['variants']:
            configs=list(itertools.product(p['dimensions'],contract['summands']))
            for a,b in itertools.permutations(configs,2):
                if abs(a[0]-b[0])+abs(a[1]-b[1])!=1:continue
                g=gate(dlpGroups[(p['name'],*a,v)],dlpGroups[(p['name'],*b,v)])
                changes.append({'profile':p['name'],'variant':v,'baseline':list(a),'candidate':list(b),**g})
        for ell,k in itertools.product(p['dimensions'],contract['summands']):
            variantPairs.append({'profile':p['name'],'ell':ell,'k':k,
                **gate(dlpGroups[(p['name'],ell,k,'enumerate-last')],dlpGroups[(p['name'],ell,k,'pair-mitm')])})
    regression=read(regressionPath)
    assert regression['stage_run_count']==240 and regression['dlp_run_count']==48
    assert regression['dlp_oracle_checks']==384
    return {'source':str(directory.resolve().relative_to(ROOT)),'raw_sha256':digest(directory/'raw.jsonl'),
        'comparison_source_sha256':digest(Path(__file__)),'expected_cell_coverage_verified':True,
        'counts':{**summary['counts'],'rho':len(rhos),'census':len(census),
            'independently_audited_dlp_oracle_calls':sum(len(r['attempts']) for r in dlps)},
        'validation':summary['validation'],'statuses':{kind:dict(Counter(r['status'] for r in raw if r['kind']==kind)) for kind in ('stage','algebraic','dlp','rho')},
        'failure_reasons':dict(Counter(r.get('error') for r in dlps if r['status']!='verified')),
        'census':census,'stage_summary':stageSummary,'matched_stage_comparisons':paired,'dlp_summary':dlpSummary,
        'neighbor_comparisons':changes,'passing_neighbor_gates':[r for r in changes if r['field_vector_gate']],
        'matched_e2e_solver_comparisons':variantPairs,'regression_replay':regression,
        'classification':'Accounting/coverage pilot; component-wise engineering candidates only, pending calibration.',
        'end_to_end_speedup_established':False,'full_dlp_S':None,'rho_ratio':None,'cost_over_floor':None}


def markdown(report):
    c=report['counts'];lines=['# Joint dimension / decomposition-length results','',
        'Frozen run: `results_v1/`; source and inputs pinned in its provenance. This is an accounting/coverage pilot, not a calibrated ECDLP speedup.','',
        f"{c['stage']} generic stage cells, {c['algebraic']} algebraic controls, {c['dlp']} cold IC runs, {c['rho']} matched rho runs, {c['census']} exact coverage cells. "
        f"{c['verified_dlp']} IC scalars verified; {c['incomplete_dlp']} incomplete runs retained. "
        f"All {c['independently_audited_dlp_oracle_calls']} full-run decomposition calls were independently audited, including empty results.",'',
        'The full previous benchmark replay retains 240 stage cells and 48 verified DLP runs; `regression_comparison.json` checks exact completed-baseline operation counters and 384 oracle calls.','',
        '## Exact coverage, including failures','',
        'Each probability covers **all affine targets** in the actual distinct-abscissa domain, not just successful sampled targets. M is the actual rational abscissa count; B is the number of distinct nonzero columns after cofactor projection. The paper expression is an intensity, not a probability above one.','',
        '| Profile | ell | k | M | B | Covered / affine | Exact coverage | Paper intensity | Poisson(actual mean) |',
        '|:--|--:|--:|--:|--:|:--|--:|--:|--:|']
    for r in report['census']:
        lines.append(f"| {r['profile']} | {r['ell']} | {r['k']} | {r['raw_abscissae']} | {r['projected_columns']} | {r['covered_targets']}/{r['affine_targets']} | {r['exact_uniform_coverage']:.4f} | {r['paper_intensity']:.4f} | {r['poisson_from_actual_mean']:.4f} |")
    lines+=['','## Full cold pipeline: separate units, all attempts','',
        'Each row sums **both fixed seeds**, including incomplete work. Compare cost only where both workloads complete. A/M/Q are binary-field API additions/multiplications/squarings; scalar A/M/I are separate modulo-subgroup operations. Inversions in the binary field are already expanded into field operations. No conversion prices Python control, lookups, field operations and scalar arithmetic together. Consequently all S/rho/floor ratios are N/A; do not turn the vector into a single operation score.','',
        '| Profile | ell/k | Variant | Verified | Field A | Field M | Field Q | Scalar A/M/I | Cold seconds | S | /rho | /floor | Class |',
        '|:--|:--|:--|:--|--:|--:|--:|:--|--:|:--|:--|:--|:--|']
    for r in report['dlp_summary']:
        f=r['field_api_all_attempts'];s=r['scalar_all_attempts'];param='—' if r['ell'] is None else f"{r['ell']}/{r['k']}"
        label='Reference' if r['variant'] in ('rho','enumerate-last') else 'Engineering candidate'
        if r['verified']<r['attempted']:label+='; incomplete'
        lines.append(f"| {r['profile']} | {param} | {r['variant']} | {r['verified']}/{r['attempted']} | {f['additions']:,} | {f['multiplications']:,} | {f['squarings']:,} | {s['additions']}/{s['multiplications']}/{s['inversions']} | {r['cold_seconds_all_attempts']:.6f} | N/A | N/A | N/A | {label} |")
    lines+=['','## Predeclared neighboring-parameter gate','',
        'For each directed neighboring change, both individual seeds must verify and each field component must not rise, with at least one decreasing. Scalar nonincrease is reported separately; neither flag establishes a common-unit end-to-end gain. All passing cases are listed; no selected winner is promoted as universal.','',
        '| Profile | Variant | Baseline ell/k | Candidate ell/k | Field-vector gate | Scalar nonincrease |',
        '|:--|:--|:--|:--|:--|:--|']
    for r in report['passing_neighbor_gates']:
        lines.append(f"| {r['profile']} | {r['variant']} | {'/'.join(map(str,r['baseline']))} | {'/'.join(map(str,r['candidate']))} | Pass on both seeds | {r['scalar_nonincrease']} |")
    lines+=['','## Retained limits','',json.dumps(report['failure_reasons'],indent=2),'',
        'Single stage repetition; two tiny DLP seeds. Wall times are descriptive, not a statistically established runtime gain. Generic controls cover k=2,3,4 and GF(8) supports; cached pullback/S4/S3 are measured only at q=2,k=3. No general-k Nagao, GF(8) polynomial adapter, F4 expansion, large-prime variation, prime-field result or asymptotic exponent fit is claimed.','',
        'Detailed stage timings, operation vectors, paired cells, phase totals, ranks, certificates and all failed attempts are in `comparison.json` and `results_v1/raw.jsonl`. See `README.md` for the interpretation and rerun commands.','']
    return '\n'.join(lines)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True)
    p.add_argument('--targets',type=Path,required=True);p.add_argument('--regression',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--markdown',type=Path);args=p.parse_args()
    report=compare(args.input,args.targets,args.regression)
    with args.output.open('x') as out:json.dump(report,out,indent=2);out.write('\n')
    if args.markdown:
        with args.markdown.open('x') as out:out.write(markdown(report))
    print(json.dumps({k:report[k] for k in ('counts','statuses','failure_reasons','passing_neighbor_gates','end_to_end_speedup_established')},indent=2))


if __name__=='__main__':main()
