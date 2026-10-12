#!/usr/bin/env python3
"""Render static report/scoreboard tables from the archived comparison."""
import argparse
from collections import Counter, defaultdict
import html
import json
from pathlib import Path
import statistics

HERE = Path(__file__).resolve().parent


def table(headers, rows, web=False):
    if web:
        return '<table><thead><tr>' + ''.join('<th>' + html.escape(str(x)) + '</th>' for x in headers) + '</tr></thead><tbody>\n' + '\n'.join('<tr>' + ''.join('<td>' + html.escape(str(x)) + '</td>' for x in row) + '</tr>' for row in rows) + '\n</tbody></table>\n'
    return '| ' + ' | '.join(headers) + ' |\n|' + '|'.join('---' for _ in headers) + '|\n' + '\n'.join('| ' + ' | '.join(map(str, row)) + ' |' for row in rows) + '\n'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('directory', type=Path)
    args = ap.parse_args()
    comp = json.loads((args.directory / 'comparison.json').read_text())
    structure = json.loads((args.directory / 'structure.json').read_text())
    records = [json.loads(line) for line in (args.directory / 'processes.jsonl').read_text().splitlines()]
    ledger = []
    for variant in ['ambient_f4', 'charts_f4', 'charts_linear', 'enumerate', 'pair_table', 'rho']:
        rows = [s for s in comp['summaries'] if s['variant'] == variant and s['kind'] == ('stage' if variant == 'charts_f4' else 'rho' if variant == 'rho' else 'dlp')]
        successes = sum(r['verified_targets'] for r in rows)
        attempts = sum(r['processes'] for r in rows)
        label = f'{attempts * 8} verified stage queries; {successes} witnesses' if variant == 'charts_f4' else f'{successes} verified DLPs / {attempts} attempts'
        ledger.append([variant, 'engineering diagnostic', 'null', 'null', 'null', 'null', label])
    headers = ['Variant', 'Class', 'Total operations', 'S', 'Rho ratio', 'Floor ratio', 'Verified workload']
    census_rows = []
    seen = set()
    for key, c in sorted(comp['census'].items(), key=lambda kv:tuple(map(int,kv[0].split(':')))):
        setup, census = c['setup'], c['census']
        identity = (setup['case'], setup['domain_hash'])
        if identity in seen:
            continue
        seen.add(identity)
        config = setup['configuration']
        ranks = census['relative_product_ranks']
        census_rows.append([setup['case'], f'{config["n"]}/{config["k"]}', config['family'], config['ell'],
            '17,937' if config['family'] != 'random' else setup['seed'], census['component_count'],
            f'{min(ranks)}–{max(ranks)}', f'{census["pair_weighted_product_rank"]:.3f}', setup['points'],
            census['projected_columns'], f'{census["exact_pair_coverage"]}/{census["coverage_denominator"]}',
            census['same_component_coverage'], census['tensor_bytes']])
    stage_rows = []
    stage_detail = []
    for ci in range(10):
        for variant in ['ambient_f4', 'charts_f4', 'charts_linear', 'enumerate']:
            runs = [r for r in records if r['kind']=='stage' and r['case']==ci and r['variant']==variant and r['seed']==937]
            vals, wins, reductions, complete = [], 0, 0, 0
            for r in runs:
                rows = [json.loads(x) for x in r['stdout'].splitlines() if x.startswith('{')]
                p = {x['phase']:x for x in rows}
                targets = [x for x in rows if x['phase']=='target']
                wins += sum(x['witness'] is not None for x in targets)
                reductions += sum(x.get('reductions',0) for x in targets)
                if r['status']=='finished' and len(targets)==8 and not any(t['exhausted'] for t in targets):
                    ns=p['setup']['setup_ns']+sum(p['target_setup'][k] for k in ['field_setup_ns','target_generation_ns'])+sum(t['solve_ns']+t['verification_ns'] for t in targets)
                    vals.append(ns/1e6)
                    complete += 1
            median = statistics.median(vals) if vals else None
            stage_rows.append([ci,variant,f'{complete}/3',wins,reductions,'—' if median is None else f'{median:.3f}'])
            stage_detail.append(dict(case=ci,variant=variant,complete=complete,witnesses=wins,reductions=reductions,median_cold_ms=median))
    dlp_rows=[]
    for s in comp['summaries']:
        if s['kind'] not in ('dlp','rho'):
            continue
        ns=s['median_cold_ns_completed']
        dlp_rows.append([s['case'],s['variant'],f'{s["verified_targets"]}/{s["processes"]}',
            s['statuses'].get('timeout',0),s['completed_processes']-s['verified_targets'],
            '—' if ns is None else f'{ns/1e6:.3f}'])
    paired_rows=[]
    for p in comp['comparisons']:
        if p['split']!='fresh' or p['candidate']!='charts_linear':
            continue
        t=p['timing']
        paired_rows.append([p['kind'],p['case'],p['reference'],'0' if t is None else t['pairs'],
            '—' if t is None else f'{t["candidate_over_reference_geomean"]:.4f}',
            '—' if t is None else f'[{t["paired_bootstrap_95"][0]:.4f}, {t["paired_bootstrap_95"][1]:.4f}]'])
    prime_rows=[]
    for r in structure['rows']:
        if r['n']!=131 or r['seed']!=17 or 'unavailable' in r:
            continue
        ranks=r['relative_product_ranks']
        prime_rows.append([r['ell'],r['family'],ranks[0],min(ranks[1:]),max(ranks[1:]),f'{r["pair_weighted_product_rank"]:.3f}'])
    raw_status=Counter((r['kind'],r['status']) for r in records)
    status_rows=[[kind,raw_status[kind,'finished'],raw_status[kind,'timeout']] for kind in ['native','stage','dlp','rho']]
    introduction=(f'Implemented an exact component-coordinate S3 solver with optional linear projection. '
        f'The saved run contains {len(records)} processes, {comp["native_root_set_comparisons"]:,} matched native root-set comparisons '
        f'({2*comp["native_root_set_comparisons"]:,} reference/candidate target executions), '
        f'{comp["complete_root_set_checks"]:,} complete component root-set checks, and '
        f'{comp["full_dlp_oracle_checks"]:,} full-DLP relation-attempt cross-checks. '
        'All completed checks pass. Four new unit tests and the no-default-features library/ic build also pass. '
        'Timeouts and incomplete attempts remain in the tables. No calibrated total-cost or exponent advance is established.')
    md='# Weil-friendly factor-base composition: measured results\n\n'+introduction+'\n\n'
    md+='The common operation-unit ledger stays explicit. Partial field/XOR/reduction counters cannot price the complete attack, so total operations, `S` and boundary ratios are null for every variant, including rho. Runtime tables below are secondary shared-host diagnostics.\n\n'+table(headers,ledger)
    md+='\n## Process outcomes\n\n'+table(['Kind','Finished','Timed out'],status_rows)
    md+='\nA finished DLP process may be incomplete; only an independently verified recovered scalar counts as a solve. Controls with zero projected columns are disqualified factor bases, not speedup wins. Every timed-out process keeps its partial stdout and stderr.\n'
    md+='\n## Actual factor-base census\n\nCoverage counts distinct nonzero subgroup points reachable by two factor-base summands. Same-component coverage is the incomplete ablation. Dimensions and point counts are not interchangeable.\n\n'+table(['Case','n/k','Seed family','ell','Seed','Components','Product ranks','Pair-weighted rank','Points','Projected columns','Exact coverage','Same-component coverage','Tensor bytes'],census_rows)
    md+='\nThe degree-nine ordinary seed and unscaled-subfield control both have zero projected columns. Scaling the subfield produces two useful columns and full coverage of the 126 nonzero subgroup points. At degree fifteen over GF(8), the scaled GF(8) components keep every product rank at three; their coverage is 70/660, versus 60/660 for the ordinary seed. The complementary GF(32) control has three columns and coverage 30/660, so it avoids the zero-column trap but has lower yield. These families change the base and its counting boundary; they are not an equal-base performance ratio.\n'
    md+='\n## Fresh solver workloads\n\nEach completed process covers the same eight targets. Three repetitions use seed 937. Cold milliseconds include setup and first-witness queries; complete root enumeration is separate validation. Witnesses/reductions in timeout rows are partial observations, not comparable totals. F4 reduction counts are exact calls to differently sized systems, not common-cost operation units.\n\n'+table(['Case','Variant','Complete processes','Witnesses observed','F4 reductions observed','Median cold ms / 8 targets'],stage_rows)
    md+='\nBoth chart variants finish all 60 stage processes across development and holdout seeds. The ambient implementation times out on 33/60. Enumeration remains stronger in several binary-curve cases; the projected adapter reduces its overhead on the GF(8)-defined cases. No finite speedup ratio is assigned to a timed-out baseline. Zero-yield holdouts cannot meet a cost-per-verified-relation gate.\n'
    md+='\n## Full-DLP practicality\n\nAll cases and all four seed/scalar combinations are included. Every row has 12 attempts. Cold time is the median among verified completions only; use the completion and timeout columns to avoid survivor bias. Rho is run once per distinct curve (case 1 also covers cases 2–3; case 4 covers case 5; case 7 covers cases 8–9).\n\n'+table(['Case','Variant','Verified / attempts','Timeouts','Finished incomplete','Median cold ms, verified only'],dlp_rows)
    md+='\n## Fresh paired timing ratios\n\nCandidate is projected charts throughout. Ratios below one favor the candidate. Only matched verified DLP completions or complete eight-target stage workloads enter a pair; missing pairs block an unqualified end-to-end claim. These are descriptive paired-bootstrap intervals on a shared host, not calibrated attack-cost ratios. Repeated seeds are not independent algorithm instances.\n\n'+table(['Workload','Case','Reference','Pairs','Candidate/reference time','95% interval'],paired_rows)
    md+='\n## Degree 131: the mixed-product obstruction\n\nThe field polynomial `x^131+x^13+x^2+x+1` passes the full irreducibility check. The following exact field-algebra diagnostics use seed 17; all offset lists and seed-937 controls are in `structure.json`. No proper subfield of dimension 8, 16, 32 or 44 exists at prime degree 131.\n\n'+table(['ell','Family','Self-product rank','Mixed minimum','Mixed maximum','Pair-weighted rank'],prime_rows)
    md+='\nAt ell=44 the power seed has rank 87 only on its self-product. Relative offsets 1 and 130 have rank 130; every other nonzero offset has rank 131. The complete cover has 8,646 unordered component pairs and weighted rank 130.318. Nearby offsets already explain the loss algebraically: products with the squared seed span 130 consecutive powers, and offset two contains a full field basis. This rejects the naive assumption that the low self-product dimension survives Frobenius composition. It does not rule out every other seed construction or fast algorithm.\n'
    md+='\n## Decision and remaining boundary\n\nThe frozen 20% stage gate is met on the positive-yield fresh cases with complete ambient baselines (cases 0 and 2): paired cold cost per witness is 0.01080 and 0.001213 times ambient F4. This is an encoding-stage improvement. Zero-column cases 1 and 3 are rejected as candidate factor bases, and ambient timeouts block finite matched ratios for cases 4–9.\n\nThe stronger reference gate is not met. No fresh full-DLP chart/pair-table interval establishes a chart win, and the fresh chart/rho time ratios range from 3.19 to 17.33 across the eight useful-base cases. Rho verifies all 60 attempts. These runtime ratios remain secondary; calibrated attack-cost ratios are still null.\n\nKeep the complete-cover adapter opt-in. Prefer common-subfield scalar components as the concrete composite-degree experiment; admission still depends on projected columns, relation yield and total verified work. Mixing GF(8) and GF(32) abscissa spaces would fill all 15 dimensions in their mixed product and lose the shared-subfield benefit.\n\nAt prime degree 131, a distinct future experiment could deliberately collect only selected component combinations and price the lost relation yield; that would be an incomplete relation oracle and must report unknown for the unsearched cover. The present implementation does not make that substitution. Higher arity, degree-131 solving, chart-plan Redis serialization, calibrated operation conversions and attack-exponent claims remain outside this measured implementation.\n\nDerivations and primary references are in [README.md](README.md). The frozen thresholds and negative controls are in [contract.json](contract.json). Raw data, source/binary hashes and every paired comparison are in [results/run-001](results/run-001).\n'
    (HERE/'RESULTS.md').write_text(md)
    summary=dict(ledger=ledger,census_rows=census_rows,fresh_stage=stage_detail,dlp_rows=dlp_rows,prime_rows=prime_rows,status_rows=status_rows,
                 source_comparison='comparison.json',total_calibrated_operations=None,S=None,rho_ratio=None,floor_ratio=None)
    (args.directory/'display.json').write_text(json.dumps(summary,indent=2)+'\n')
    section='<section id="weil-factor-composition-20260914" style="margin:2rem 0;overflow-x:auto">\n<h2>Weil-friendly factor-base composition <span class="chip">engineering diagnostic</span></h2>\n<p>'+html.escape(introduction)+'</p>\n'
    section+='<p>The positive-yield cases with completed ambient controls meet the frozen stage gate. The stronger baseline gate fails: no fresh full-DLP interval establishes a win over pair tables, and projected-chart time is 3.19–17.33× matched rho across the eight useful-base cases. These are secondary runtime ratios, not calibrated attack-cost ratios.</p>\n'
    section+=table(headers,ledger,True)
    section+='<h3>Fresh solver workloads: cold milliseconds</h3><p>Eight targets per process, three repetitions; complete workloads only. Partial timeout counters remain in the report. Zero-yield controls do not establish a gain.</p>\n'+table(['Case','Variant','Complete processes','Witnesses observed','F4 reductions observed','Median cold ms / 8 targets'],stage_rows,True)
    section+='<h3>Full DLP: verified completions and cold milliseconds</h3><p>12 attempts per row. Medians exclude incomplete attempts; their counts and all timeouts are retained. Rho shares the identical curve and target, without a factor base. Calibrated S and boundary ratios remain null.</p>\n'+table(['Case','Variant','Verified / attempts','Timeouts','Finished incomplete','Median cold ms, verified only'],dlp_rows,True)
    section+='<h3>Prime degree 131: exact product ranks</h3><p>Field algebra only. At dimension 44, a power seed has self-product rank 87 but pair-weighted rank 130.318 across the complete 8,646 component pairs. No nontrivial proper subfield exists. No degree-131 DLP or exponent inference is made.</p>\n'+table(['ell','Family','Self-product rank','Mixed minimum','Mixed maximum','Pair-weighted rank'],prime_rows,True)
    section+='<p>Common-subfield scalar components preserve small mixed products at composite degree. Every completed root set and point witness is checked; a small product rank alone does not establish a better attack. Earlier calibrated rows, exponents and attack verdicts remain unchanged.</p>\n<p>Sources: <a href="../research/weil_factor_composition_20260914/RESULTS.md">results and paired intervals</a>, <a href="../research/weil_factor_composition_20260914/results/run-001/comparison.json">comparison data</a>, <a href="../research/weil_factor_composition_20260914/results/run-001/structure.json">structural diagnostics</a>, <a href="../research/weil_factor_composition_20260914/contract.json">frozen contract</a>.</p>\n</section>\n'
    page=HERE.parents[1]/'docs/index-calculus-scoreboard.html'
    text=page.read_text()
    marker='<section id="weil-factor-composition-20260914"'
    if marker in text:
        begin=text.index(marker);end=text.index('</section>',begin)+len('</section>\n')
        text=text[:begin]+text[end:]
    assert '<script>\n  (function ()' in text
    text=text.replace('<script>\n  (function ()',section+'<script>\n  (function ()',1)
    page.write_text(text)
    print('Rendered RESULTS.md, display.json and canonical scoreboard from archived measurements.')


if __name__=='__main__':
    main()
