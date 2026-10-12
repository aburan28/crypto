#!/usr/bin/env python3
"""Render the note and canonical scoreboard only from saved evidence."""
from collections import Counter, defaultdict
import gzip
import html
import json
import math
from pathlib import Path
import random
import statistics

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CASES = json.loads((HERE.parent/'weil_factor_composition_20260914/contract.json').read_text())['stage_cases']
ORDERS = [71,127,127,127,991,991,2003,661,661,661]
# Immutable commit containing the complete accepted evidence, before link/CI
# integration fixes. Absolute sources survive the Pages /scoreboard/ layout.
EVIDENCE_URL = 'https://github.com/aburan28/crypto/blob/7259f1f6593078bb5d44299b0169b2079b0c0758/'


def table(head, rows, web=False):
    if web:
        return '<table><thead><tr>'+''.join('<th>'+html.escape(str(c))+'</th>' for c in head)+'</tr></thead><tbody>'+''.join(
            '<tr>'+''.join('<td>'+html.escape(str(c))+'</td>' for c in row)+'</tr>' for row in rows)+'</tbody></table>\n'
    return '| '+' | '.join(head)+' |\n|'+'|'.join(['---']*len(head))+'|\n'+''.join('| '+' | '.join(map(str,r))+' |\n' for r in rows)


def paired_mean_ratio(pairs):
    if not pairs:
        return None
    def ratio(ps):
        return sum(x for x,y in ps)/sum(y for x,y in ps)
    rng = random.Random(20260915)
    boot = sorted(ratio(rng.choices(pairs,k=len(pairs))) for _ in range(4000))
    return dict(ratio=ratio(pairs), ci95=[boot[100],boot[3899]], pairs=len(pairs))


def fmt(x, digits=3):
    return 'null' if x is None else f'{x:.{digits}f}'


def main():
    comparisons = {p.parent.name:json.loads(p.read_text()) for p in sorted((HERE/'results').glob('*/comparison.json'))}
    rounds = [k for k in comparisons if k.startswith('iteration-')]
    counts = {p.parent.name:json.loads(p.read_text()) for p in (HERE/'counts').glob('*/counts.json')
              if not (p.parent/'REJECTED.json').exists()}
    latest = max(rounds, key=lambda k:int(k.split('-')[-1]))
    final = 'confirmation' if 'confirmation' in comparisons else latest
    instruction_candidate = 'confirmation-candidate' if 'confirmation-candidate' in counts else latest
    instruction_reference = 'confirmation-reference' if 'confirmation-reference' in counts else 'reference'
    required = {'reference', *rounds}
    if final == 'confirmation':
        required.update(['confirmation-reference', 'confirmation-candidate'])
    assert required <= counts.keys(), 'incomplete instruction corpus'
    audit = json.loads((HERE/'instruction-audit.json').read_text())
    assert required <= {r['corpus'] for r in audit['accepted']}, 'native/instruction audit required'
    totals = Counter()
    for c in comparisons.values():
        totals.update(c['checks'])
    runtime_rows = []
    for name in rounds+[k for k in ['confirmation'] if k in comparisons]:
        for ci in range(10):
            rs = [r for r in comparisons[name]['comparisons'] if r['case']==ci and r['variant']=='charts_linear' and r['split']=='fresh']
            row = {r['reference']:r['runtime'] for r in rs}
            ratio = row['rho']
            runtime_rows.append([name,ci,fmt(row['reference']['ratio']) if row['reference'] else '—',
                                 fmt(ratio['ratio']) if ratio else '—',
                                 '['+', '.join(fmt(x) for x in ratio['ci95'])+']' if ratio else '—'])
    instruction_rows, instruction_data, instruction_fresh, phase_rows = [], [], [], []
    cohorts = [('validation', latest, 'reference', rounds)]
    if final == 'confirmation':
        cohorts.append(('confirmation', instruction_candidate, instruction_reference, []))
    for cohort, candidate_name, reference_name, revisions in cohorts:
        candidate = counts[candidate_name]
        reference = counts[reference_name]
        for ci,c in enumerate(CASES):
            rc = next(j for j,cc in enumerate(CASES) if all(c[k]==cc[k] for k in ['n','k','a','b']))
            cal = {r['mode']:r for r in candidate if r['case']==rc and r['mode'].startswith('calibrate')}
            per_add = (cal['calibrate']['scopes']['calibration_loop']-cal['calibrate_empty']['scopes']['calibration_loop'])/4096
            assert per_add > 0
            rho = {r['seed']:r for r in candidate if r['case']==rc and r['mode']=='rho'}
            variants = [('original charts',reference,'charts_linear')]
            variants += [(name+' charts',counts[name],'charts_linear') for name in revisions] if revisions else [('candidate charts',candidate,'charts_linear')]
            variants += [
                        ('original pair table',reference,'pair_table'),('candidate pair table',candidate,'pair_table'),
                        ('rho',candidate,'rho')]
            for label, corpus, variant in variants:
                selected = [r for r in corpus if r['case']==(rc if variant=='rho' else ci) and r['variant']==variant]
                verified = [r for r in selected if r['verified'] and r['status']=='finished']
                assert all(r['known_log']==rho[r['seed']]['known_log'] for r in selected)
                pairs = [(r['total_instructions'],rho[r['seed']]['total_instructions']) for r in verified if rho[r['seed']]['verified']]
                ratio = paired_mean_ratio(pairs)
                mean_i = statistics.mean(r['total_instructions'] for r in verified) if verified else None
                s = mean_i/per_add/math.sqrt(ORDERS[ci]) if mean_i else None
                original = {r['seed']:r for r in reference if r['case']==ci and r['variant']==variant}
                improvement = paired_mean_ratio([(original[r['seed']]['total_instructions'],r['total_instructions'])
                    for r in verified if r['seed'] in original and original[r['seed']]['verified']])
                row = dict(cohort=cohort,case=ci,variant=label,attempts=len(selected),verified=len(verified),mean_instructions=mean_i,
                           instructions_per_addition=per_add,S_I=s,rho_instruction_ratio=ratio,floor_ratio=None,
                           baseline_cost_speedup=improvement,
                           classification='engineering; machine-specific instruction model')
                instruction_data.append(row)
                if corpus is candidate and variant == 'charts_linear' and ci not in [1,3] and candidate_name == instruction_candidate:
                    fresh = [r for r in verified if r['seed']>=1000]
                    instruction_fresh.append(dict(case=ci, ratio=paired_mean_ratio([
                        (r['total_instructions'],rho[r['seed']]['total_instructions']) for r in fresh]),
                        baseline_cost_speedup=paired_mean_ratio([(original[r['seed']]['total_instructions'],r['total_instructions']) for r in fresh])))
                    phase_rows.append([ci]+[fmt(statistics.mean(r['scopes'][scope] for r in fresh)/1e6) for scope in
                        ['curve_setup','target_generation','factor_base_plan_index','driver_and_final_verification']])
                instruction_rows.append([cohort,ci,label,f'{len(verified)}/{len(selected)}',fmt(mean_i/1e6) if mean_i else 'null',
                                         fmt(s,2),fmt(ratio['ratio']) if ratio else 'null',
                                         fmt(improvement['ratio']) if improvement else 'null','null','engineering'])
    correct_rows = [[name,c['processes'],c['checks']['native_targets'],c['checks']['complete_root_sets'],c['checks']['oracle_attempts']]
                    for name,c in comparisons.items()]
    target_rows, instruction_steps = [], []
    for i, name in enumerate(rounds):
        before = counts['reference' if i == 0 else rounds[i-1]]
        after = counts[name]
        instruction_met, runtime_met = [], []
        for ci in [0,2,4,5,6,7,8,9]:
            x = {r['seed']:r for r in before if r['mode']=='dlp' and r['variant']=='charts_linear' and r['case']==ci and r['seed']>=1000}
            y = {r['seed']:r for r in after if r['mode']=='dlp' and r['variant']=='charts_linear' and r['case']==ci and r['seed']>=1000}
            assert x.keys() == y.keys() and all(r['verified'] for r in [*x.values(), *y.values()])
            ratio = paired_mean_ratio([(y[s]['total_instructions'], x[s]['total_instructions']) for s in x])
            instruction_steps.append(dict(iteration=name, case=ci, candidate_predecessor_ratio=ratio))
            if ratio['ratio'] <= .8 and ratio['ci95'][1] < 1:
                instruction_met.append(ci)
            timing = next(r['runtime'] for r in comparisons[name]['comparisons'] if r['case']==ci and r['variant']=='charts_linear' and r['split']=='fresh' and r['reference']=='reference')
            if timing and timing['ratio'] <= .8 and timing['ci95'][1] < 1:
                runtime_met.append(ci)
        target_rows.append([name, ', '.join(map(str,instruction_met)) or 'none', ', '.join(map(str,runtime_met)) or 'none'])
    final_comparisons = [r for r in comparisons[final]['comparisons'] if r['variant']=='charts_linear' and r['split']=='fresh' and r['reference']=='rho' and r['case'] not in [1,3]]
    parity = all(r['missing']==0 and r['runtime'] and r['runtime']['ci95'][1]<=1 for r in final_comparisons)
    instruction_parity = len(instruction_fresh)==8 and all(r['ratio'] and r['ratio']['pairs']==4 and r['ratio']['ci95'][1]<=1 for r in instruction_fresh)
    matrix_times = defaultdict(list)
    with gzip.open(HERE/'results'/final/'processes.jsonl.gz', 'rt') as f:
        for line in f:
            r = json.loads(line)
            if r['kind'] == 'dlp' and r['revision'] == 'candidate' and r['variant'] == 'charts_linear' and r['seed'] >= 1000:
                terminal = next(json.loads(s) for s in r['stdout'].splitlines() if s.startswith('{') and json.loads(s)['phase'] == 'dlp')
                if terminal['verified']:
                    matrix_times[r['case']].append((terminal['linear_algebra_ns'], terminal['cold_ns']))
    matrix_bounds = [dict(case=ci, timed_matrix_fraction=sum(x for x,y in times)/sum(y for x,y in times),
                         free_timed_matrix_speedup=sum(y for x,y in times)/sum(y-x for x,y in times))
                     for ci,times in sorted(matrix_times.items())]
    verdict = 'The frozen toy-suite parity gates are met in both measured models.' if parity and instruction_parity else 'Rho parity is not established: the frozen all-cases gates are not met.'
    intro = f'{verdict} {len(rounds)} engineering hypotheses were implemented and compared against their predecessors. '
    intro += 'The strongest structural change compresses the rational support of the degree-eleven ordinary seed from three coordinates to one while retaining every factor-base point. '
    intro += 'All claims below concern the measured small fields, with degree-131 solving and asymptotic parity unestablished.'
    if instruction_fresh:
        ratios = [r['ratio']['ratio'] for r in instruction_fresh]
        gains = [r['baseline_cost_speedup']['ratio'] for r in instruction_fresh]
        intro += f' On fresh confirmation targets, the candidate uses {min(ratios):.2f}–{max(ratios):.2f} times rho\'s instructions. Its full instruction cost is {100*(1-1/min(gains)):.1f}–{100*(1-1/max(gains)):.1f}% below the original chart implementation.'
    md = '# Full-cost index-calculus iterations against rho\n\n'+intro+'\n\n'
    md += 'The public algebraic visitor still returns every original S3 root. Group-witness search may use the smaller rational support. Point lifts, Frobenius coefficients and fixed-base point powers are reused; small projected systems use packed checks and fixed-capacity row elimination. All setup is charged, direct-relation shortcuts and caches remain disabled.\n\n'
    case_rows = [[ci,c['n'],c['k'],c['ell'],c['family'],ORDERS[ci], 'zero-column control' if ci in [1,3] else 'useful'] for ci,c in enumerate(CASES)]
    md += table(['Case','n','k','Seed dimension','Family','Subgroup order N','Role'],case_rows)+'\n'
    md += '## Correctness and archived workload\n\nEvery completed full-DLP oracle input is checked against independent pair truth, every reported logarithm verifies `[k]G=Q`, and matched coefficient/target prefixes agree. Incomplete zero-column controls remain in every suite; they never count as speedups.\n\n'
    md += table(['Comparison','Processes','Matched native targets','Complete algebraic root sets','Oracle checks'],correct_rows)
    md += '\n## Twenty-percent engineering targets\n\nThese lists identify useful cases meeting a candidate/predecessor cost ratio at most 0.8, with the paired 95% upper limit below 1. All other useful cases miss that numeric target or its confidence gate. The original four new inputs are reused validation inputs after iteration 1. These engineering gates are separate from the all-cases rho parity gate, which remains unmet.\n\n'
    md += table(['Iteration','Cases meeting instruction target','Cases meeting native runtime target'],target_rows)
    md += '\n## Complete instruction accounting\n\nThe table uses one machine-specific instruction model. `S_I = mean cold Ir / (measured Ir per addition × sqrt(N))`. The conversion uses 4,096 public random point additions minus a loop control on the same curve; all variants share the candidate calibration. `Ir/rho` is the ratio of arithmetic mean cold instruction counts on identical targets. This is distinct from the historical algebraic-operation S; it does not establish an attack exponent or a calibrated mathematical floor. Kernel execution and physical memory latency are outside this model. Empty controls have no verified-cost result.\n\n'
    cost_header = ['Cohort','Case','Variant','Verified','Mean cold million Ir','S_I','Ir/rho','Original/variant Ir','Floor ratio','Class']
    md += table(cost_header,instruction_rows)
    md += '\nThe validation cohort includes the four frozen targets and four reused validation targets. Every retained candidate is shown on those identical inputs; earlier costs are preserved. The confirmation cohort uses the four frozen targets and four new uniformly sampled nonzero scalars. The fresh-only parity gate uses the latter four, and does not pool them with tuning inputs.\n\n'
    md += '### Exclusive candidate phases on fresh confirmation targets\n\n'
    md += table(['Case','Curve setup million Ir','Target generation million Ir','Base/plan/index million Ir','Driver/verification million Ir'],phase_rows)
    md += '\n### Remaining gap\n\nThe original factor base and counting bound are unchanged. Rational-support compression helps the degree-eleven ordinary seed but does not reduce the degree-thirteen seed dimension. The latter still pays substantial work for failed chart queries and relation generation. The final relation-matrix timer is a small part of cold native runtime; setting just that measured timer to zero gives the following Amdahl limits. These are runtime bounds for that scope, not operation-count floors and not bounds on all polynomial preprocessing or chart projection.\n\n'
    md += table(['Case','Timed relation-matrix fraction','Maximum speedup from removing that timer'],
                [[r['case'],f"{100*r['timed_matrix_fraction']:.2f}%",fmt(r['free_timed_matrix_speedup'])] for r in matrix_bounds])
    md += '\nScopes include curve setup, target generation, factor-base/plan/index setup, and the entire driver plus final verification. Diagnostic JSON and independent exhaustive oracle replay are outside those scopes. Raw profiles, exclusive phase totals, compiler/binary hashes, and calibration controls are archived in `counts/`. No Valgrind runtime is used as a native timing result.\n\n'
    md += '## Secondary native runtime comparisons\n\nRatios below one favor the candidate. Intervals are paired 95% seed-cluster bootstrap intervals; three repetitions share a seed and are not treated as three independent instances. Iteration rows compare against their immediate predecessor. The confirmation compares directly against the original revision on newly sampled uniform targets. The original four new seed/log pairs became reused validation inputs after iteration 1. No cross-round timing ratios are multiplied.\n\n'
    md += table(['Comparison','Case','Candidate/predecessor','Candidate/rho','Rho ratio 95% interval'],runtime_rows)
    md += '\n## Limits and reproduction\n\nThe two zero-column controls cannot recover a logarithm. Cold single-target costs are reported; there is no hidden warm amortization, precomputed target logarithm, or rho fallback. The small-coordinate adapter remains opt-in, supports at most seven seed coordinates, and uses n<=63 field arithmetic (packed point arithmetic through n<=62). No Redis/network latency or degree-131 solve is inferred from this cache-off local experiment.\n\n'
    md += 'The instruction audit checks every non-timing field emitted by the immutable cost harness against the matching native run, including F4 reductions and every oracle witness. A stale shared-target reference build and two truncated raw-profile archives were rejected and retained with reasons. Accepted replacements use isolated package builds and validated, completed profile archives. The first iteration-5 native run overlapped profiling and is retained as a correctness diagnostic; its timing values are excluded from the tables. All native timings remain shared-host diagnostics.\n\n'
    md += 'See [README.md](README.md) for hypotheses frozen before each change, [COST_PROTOCOL.md](COST_PROTOCOL.md) for the instruction model, and [support-census.json](support-census.json) for the independent rational-support census. Use `build.py` with a new target directory per revision. Run `run.py`, then `compare.py` (the witness-order flag is required for the support/visitation experiments); `count.py` requires the identically generated `rho_parity_cost` executable and Valgrind. Run `audit.py` before `report.py`, which renders this note and the canonical scoreboard from archived comparisons.\n'
    (HERE/'RESULTS.md').write_text(md)
    display = dict(verdict=verdict,runtime_parity=parity,instruction_parity=instruction_parity,
                   instruction_fresh=instruction_fresh,confirmation=final=='confirmation',correctness=dict(totals),
                   instruction_rows=instruction_data,instruction_steps=instruction_steps,engineering_targets=target_rows,
                   phase_rows=phase_rows,matrix_bounds=matrix_bounds,runtime_rows=runtime_rows,sources=list(comparisons))
    (HERE/'display.json').write_text(json.dumps(display,indent=2)+'\n')
    section = '<section id="rho-parity-20260915" style="margin:2rem 0;overflow-x:auto">\n<h2>Full-cost rho parity iterations <span class="chip">engineering</span></h2>\n<p>'+html.escape(intro)+'</p>\n'
    section += '<p>S_I is a machine-specific instruction model, normalized by a measured addition cost. It is distinct from the historical algebraic-operation S. The mathematical floor ratio remains uncalibrated. All precomputation and verification are included; zero-column controls are incomplete.</p>\n'
    section += '<details><summary>Case definitions</summary>\n'+table(['Case','n','k','Seed dimension','Family','Subgroup order N','Role'],case_rows,True)+'</details>\n'
    section += '<details><summary>Complete instruction table: every retained iteration and independent confirmation</summary>\n'+table(cost_header,instruction_rows,True)+'</details>\n'
    section += '<h3>Exclusive candidate phases on fresh confirmation targets</h3>\n'+table(['Case','Curve setup million Ir','Target generation million Ir','Base/plan/index million Ir','Driver/verification million Ir'],phase_rows,True)
    section += '<p>Removing only the timed final relation-matrix work cannot close the measured gap. Chart projection, failed queries and relation generation remain separate costs; see the saved Amdahl diagnostics in the results note.</p>\n'
    section += '<h3>Paired native runtime diagnostics</h3><p>Three repetitions per seed; 95% seed-cluster intervals. Confirmation uses previously unseen uniform targets. Earlier iterations retain their prior values.</p>\n'
    section += table(['Comparison','Case','Candidate/predecessor','Candidate/rho','95% interval'],runtime_rows,True)
    section += '<p>Sources: <a href="../research/rho_parity_20260915/RESULTS.md">results</a>, <a href="../research/rho_parity_20260915/display.json">saved display values</a>, <a href="../research/rho_parity_20260915/COST_PROTOCOL.md">accounting protocol</a>. Earlier calibrated rows and exponent claims are unchanged; this round makes no exponent claim.</p>\n</section>\n'
    page = ROOT/'docs/index-calculus-scoreboard.html'
    content = page.read_text()
    old_verdict = '      do not establish an end-to-end speedup or belong in that normalized claim.</p>'
    if old_verdict in content:
        content = content.replace(old_verdict, '      do not establish an end-to-end speedup or belong in that normalized claim.\n      The new full-cost Weil-chart iterations use a separately labeled instruction model;\n      they improve the implementation but still do not establish rho parity.</p>', 1)
    marker = '<section id="rho-parity-20260915"'
    if marker in content:
        begin=content.index(marker);end=content.index('</section>',begin)+len('</section>\n')
        content=content[:begin]+content[end:]
    assert '<script>\n  (function ()' in content
    content = content.replace('<script>\n  (function ()',section+'<script>\n  (function ()',1)
    content = content.replace('href="../research/', 'href="'+EVIDENCE_URL+'research/')
    page.write_text(content)
    print(verdict)


if __name__ == '__main__':
    main()
