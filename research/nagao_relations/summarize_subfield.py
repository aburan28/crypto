"""Render the frozen subfield stage diagnostic; no extrapolated costs."""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FOLDER = HERE / 'subfield_01'


def render():
    raw = [json.loads(s) for s in (FOLDER / 'raw.jsonl').read_text().splitlines()]
    audit = json.loads((FOLDER / 'audit.json').read_text())
    trials = [r for r in raw if r['kind'] == 'trial']
    supplemental = [r for r in raw if r['kind'] == 'supplemental']
    panels = sorted({(r['n'], r['d']) for r in trials})
    variants = ['hybrid-image', 'chained-s3', 's4-symmetric']
    tables = {}
    htmlTables = []
    for mode in ('first', 'enumerate'):
        headings = ['Variant'] + ['n%d d%d' % p for p in panels] + ['Class', 'S / rho / cost speedup']
        rows = []
        for variant in variants:
            row = [variant]
            for n, d in panels:
                rr = [r for r in trials if (r['n'], r['d'], r['variant'], r['mode']) == (n, d, variant, mode)]
                row.append('%d/%d' % (sum(r['status'] != 'timeout' and r['within_budget'] for r in rr), len(rr)))
            row += ['engineering diagnostic', 'unmeasured']
            rows.append(row)
        tables[mode] = '\n'.join(['| ' + ' | '.join(headings) + ' |',
                                  '|' + '|'.join(['---'] * len(headings)) + '|'] +
                                 ['| ' + ' | '.join(r) + ' |' for r in rows])
        htmlTables.append('<div class="table-scroll"><table><caption>' + mode +
                          ': resolved within five seconds / attempts; includes uniform negative targets.</caption>' +
                          '<thead><tr>' + ''.join('<th>' + h + '</th>' for h in headings) + '</tr></thead><tbody>' +
                          ''.join('<tr>' + ''.join('<td>' + x + '</td>' for x in r) + '</tr>' for r in rows) +
                          '</tbody></table></div>')
    complete = sum(r['status'] == 'complete' and r['within_budget'] for r in supplemental)
    maximum = max(r['all_phase_seconds'] for r in supplemental)
    matched = {v: {'attempts': len([r for r in trials if r['variant'] == v]),
                   'resolved_within_budget': sum(r['status'] != 'timeout' and r['within_budget']
                                               for r in trials if r['variant'] == v),
                   'verified_relations': sum(r['verified_unique_relations'] for r in trials if r['variant'] == v)}
               for v in variants}
    lines = ['# Subfield curves with dimensions 6, 7, and 8', '',
             'The solver now handles binary curves with non-F2 coefficients in F4, evaluated over '
             'GF(2^18) and GF(2^30). The dimension-eight supplemental runs completed %d/%d exact '
             'enumerations within 60 seconds; maximum %.6f seconds. This extends the implementation '
             'and measures a larger-base obstacle. It does not establish a calibrated Semaev or ECDLP speedup.'
             % (complete, len(supplemental), maximum), '',
             'Frozen source and proofs: [subfield_01/README.md](subfield_01/README.md). '
             'Evidence: [raw.jsonl](subfield_01/raw.jsonl), [summary.json](subfield_01/summary.json), '
             '[audit.json](subfield_01/audit.json), and [contract.json](subfield_01/contract.json).', '',
             '## Matched five-second comparison', '',
             'Every table entry counts a verified first relation or a proved complete enumeration '
             '(including zero relations), as appropriate for its mode. A partial enumeration is a '
             'timeout. Budgets cover cold setup, search, extraction and verification. CryptoMiniSat '
             'uses a soft time limit; overruns remain recorded and do not count as within-budget completions.', '',
             'First relation / negative-target resolution:', '', tables['first'], '',
             'Complete enumeration:', '', tables['enumerate'], '',
             'These are 144 matched trials on 24 curve/base/target instances. For each field the '
             'same four targets are used at all dimensions: one uniform and one known-decomposable '
             'target for each of the predeclared development and holdout seeds. The supported stratum '
             'was sampled from d6 triples and is not a natural-yield estimate. Each stratum has only '
             'two targets per field, so no timing confidence interval or population success claim is made.', '',
             '## Larger bases and complete relation sets', '',
             'The following expected counts come from independent group-law pair enumeration, '
             'not from the tested solvers. Columns keep the same targets as the dimension increases.', '',
             '| Field bits | d | Signed base points | Development uniform | Holdout uniform | Development supported | Holdout supported |',
             '|---|---|---|---|---|---|---|']
    for n, d in panels:
        oracle = next(r for r in raw if r['kind'] == 'oracle' and (r['n'], r['d']) == (n, d))
        counts = []
        for st in ('uniform', 'known_decomposable'):
            for cohort in ('development', 'holdout'):
                r = next(r for r in raw if r['kind'] == 'instance' and
                         (r['n'], r['d'], r['stratum'], r['cohort']) == (n, d, st, cohort))
                counts.append(str(len(r['expected'])))
        lines.append('| %d | %d | %d | %s |' % (n, d, oracle['signed_base_size'], ' | '.join(counts)))
    lines += ['', 'Supplemental cold hybrid enumeration at d8, with a separate 60-second budget. '
              'No Semaev comparison uses these unmatched runs:', '',
              '| Field bits | Stratum | Cohort | Status | Verified relations | Seconds | Field API operations |',
              '|---|---|---|---|---|---|---|']
    for r in supplemental:
        lines.append('| %d | %s | %s | %s | %d | %.6f | %d |' %
                     (r['n'], r['stratum'], r['cohort'], r['status'], r['verified_unique_relations'],
                      r['all_phase_seconds'], r['field_api_counts']['totals']['fieldOperations']))
    lines += ['', '## What is proved and what remains open', '',
              '- The general-coefficient norm formula, quadratic reduction, and image support test '
              'are exact in the declared distinct-abscissa domain. The general-B symmetric S4 '
              'identity is checked symbolically against the resultant of two S3 polynomials.',
              '- All 53 affine GF64 targets agree across the hybrid, both Semaev controls, '
              'signed-triple enumeration, and the pair oracle. All 7,582 candidate support checks '
              'agree with modular H|L_V. All 576 AS inputs over GF64/GF512 pass; the GF64 point '
              'lifting set matches all 4,096 possible coordinate pairs.',
              '- The final audit replays %d signed-point certificates and independently recomputes '
              'all 24 larger oracle sets. Source and instance hashes, matched coverage, fixed targets '
              'across dimensions, exclusive timers, and summary aggregates pass with zero validation failures.'
              % audit['signed_certificates_replayed'],
              '- Curves over F4 are supported here; odd-characteristic subfield curves and F4-linear '
              'factor bases are not tested. The inherited coordinate Frobenius is not used for '
              'automorphism compression: only powers fixing the curve coefficients act on the same curve.', '',
              'The branch search remains O(2^(2d)) before field costs. Increasing d moves the '
              'counting bound 8 binom(M/2,3)/(#E-1); extra relations are not evidence of beating a '
              'fixed counting boundary. This round is classified as a functionality extension and '
              'engineering diagnostic. All common-operation speedups, full-DLP S, rho ratios and '
              'floor ratios remain null. The SAT field counters omit SAT and encoding work, so '
              'comparing them to the hybrid counters would be invalid. Only field additions, '
              'multiplications, squarings and inversionCalls are instrumented. Inherited curve-call '
              'fields in the raw CountedField reports are zero placeholders, not measurements; '
              'the audit marks those call counts null and all aggregates exclude them.', '',
              'The parent accounting/promotion contract remains unmet: this small stage panel '
              'does not replace the required solver regression, calibrated full-pipeline comparison, '
              'independent-curve holdouts or paired timing repetitions. No relation matrix or final '
              'scalar recovery is run. The user’s 20% all-cost goal remains open.', '',
              'Next useful experiment: reduce the quadratic dependence on base size or safely '
              'amortize image setup over a declared target batch, then run the broader matched '
              'regression with calibrated Boolean/field accounting. The present data alone do '
              'not justify an exponent fit.', '']
    (HERE / 'subfield_dimension_results.md').write_text('\n'.join(lines))
    result = {'matched': matched, 'supplemental_completions': complete,
              'supplemental_attempts': len(supplemental), 'supplemental_max_seconds': maximum,
              'audit': audit, 'classification': 'engineering diagnostic',
              'common_operation_speedup': None, 'full_dlp_S': None, 'rho_ratio': None, 'floor_ratio': None}
    (HERE / 'subfield_dimension_results.json').write_text(json.dumps(result, indent=2) + '\n')
    panel = ('  <div class="panel" id="nagao-subfield-dimension"><div class="panel-head">'
             '<h2>Non-F2 subfield curves, with dimensions 6, 7 and 8</h2>'
             '<p>Frozen source: <code>research/nagao_relations/subfield_dimension_results.json</code> '
             'and <code>subfield_01/raw.jsonl</code>. 144 matched five-second trials on F4-defined '
             'curves over GF(2^18) and GF(2^30), plus eight explicitly unmatched 60-second hybrid '
             'enumerations at dimension eight. All %d/8 supplemental enumerations completed; '
             'maximum %.6f seconds. Zero validation failures. This is a functionality extension '
             'and engineering diagnostic; calibrated cost speedup and full-DLP ratios remain unmeasured.</p>'
             '</div>%s<div class="panel-head"><p>First-relation and complete-enumeration results '
             'are separate. Bases use 64/128/256 candidate abscissas before curve membership and '
             'exclusions. The same targets are retained across dimensions. Larger support moves '
             'the counting boundary; it does not prove an attack advance. The required broader '
             'regression and full-pipeline performance gate remain unmet.</p></div></div>\n\n'
             % (complete, maximum, ''.join(htmlTables)))
    scoreboard = ROOT / 'docs/index-calculus-scoreboard.html'
    text = scoreboard.read_text()
    if 'id="nagao-subfield-dimension"' in text:
        raise ValueError('scoreboard panel already exists; do not duplicate it')
    marker = '  <!-- ============ NOTES ============ -->'
    if text.count(marker) != 1:
        raise ValueError('scoreboard insertion marker mismatch')
    scoreboard.write_text(text.replace(marker, panel + marker))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    render()
