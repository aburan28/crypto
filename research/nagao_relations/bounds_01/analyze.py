"""Audit saved bound calculations and render literal research/scoreboard tables."""
from collections import Counter
from fractions import Fraction
import hashlib
import json
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def decode(row):
    value = Fraction(row['numerator'], row['denominator'])
    if float(value) != row['decimal']:
        raise ArithmeticError('inconsistent rational/decimal representation')
    return value


def curveOrder(n):
    s = [2, -1]
    for j in range(2, n // 2 + 1):
        s.append(-s[j - 1] - 4 * s[j - 2])
    return 2 ** n + 1 - s[n // 2]


def thresholds():
    rows = []
    for n in (18, 30, 42, 60, 90, 126):
        nn = curveOrder(n)
        lo, hi = 2, 4
        while 16 * comb(hi, 3) < nn - 1:
            hi *= 2
        while lo + 1 < hi:
            middle = (lo + hi) // 2
            if 16 * comb(middle, 3) >= nn - 1:
                hi = middle
            else:
                lo = middle
        if not 16 * comb(hi - 1, 3) < nn - 1 <= 16 * comb(hi, 3):
            raise ArithmeticError('necessary base-size threshold failed')
        rows.append({'n': n, 'curve_order': nn, 'minimum_K_for_half_probability_ceiling': hi,
                     'necessary_binary_dimension': hi.bit_length(),
                     'minimum_current_s3_entries': hi * (hi - 1),
                     'minimum_current_s3_setup_multiplications': 7 * comb(hi, 2),
                     'kind': 'necessary counting threshold, not an achieved success rate'})
    return rows


def main():
    rawPath = HERE / 'raw.jsonl'
    rows = [json.loads(s) for s in rawPath.read_text().splitlines()]
    contract = json.loads((HERE / 'contract.json').read_text())
    publication = json.loads((HERE / 'publication.json').read_text())
    if rows[-1]['kind'] != 'completion' or rows[-1]['failures']:
        raise ArithmeticError('incomplete or failed bound campaign')
    if rows[0]['commit'] != publication['local_frozen_source_commit']:
        raise ArithmeticError('bound publication mapping mismatch')
    for name, digest in rows[0]['sha256'].items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != digest:
            raise ArithmeticError('frozen source changed: ' + name)
    oldPath = ROOT / contract['input']
    if hashlib.sha256(oldPath.read_bytes()).hexdigest() != contract['input_sha256']:
        raise ArithmeticError('prior raw evidence changed')
    old = [json.loads(s) for s in oldPath.read_text().splitlines()]
    validation = next(r for r in rows if r['kind'] == 'validation')
    if validation['tiny_target_space_checks'] != 265 or validation['failures']:
        raise ArithmeticError('small-field validation incomplete')
    identity = json.loads((HERE / 'pullback_identity_results.json').read_text())
    if identity['source_commit'] != publication['identity_local_source_commit']:
        raise ArithmeticError('identity publication mapping mismatch')
    if hashlib.sha256((HERE / 'pullback_identity.py').read_bytes()).hexdigest() != identity['source_sha256']:
        raise ArithmeticError('identity source changed')
    if identity['identity_checks'] != 2048 or identity['exhaustive_nonzero_GF64_target_checks'] != 756 or identity['failures']:
        raise ArithmeticError('identity checks incomplete')
    spaces = [r for r in rows if r['kind'] == 'space_bounds']
    descriptors = [r for r in old if r['kind'] == 'space']
    key = lambda r: (r['n'], r['base_kind'], r['d'])
    if sorted(map(key, spaces)) != sorted(map(key, descriptors)):
        raise ArithmeticError('frozen space coverage changed')
    targetChecks = 0
    fresh = 0
    for row in spaces:
        descriptor = next(r for r in descriptors if key(r) == key(row))
        if any(row[k] != descriptor[k] for k in ('basis', 'coefficients')):
            raise ArithmeticError('bound support changed')
        n, k = row['curve_order'], row['K_abscissas']
        if n != curveOrder(row['n']):
            raise ArithmeticError('curve order mismatch')
        total = row['valid_signed_triples']
        if 8 * comb(k, 3) - row['infinity_triples'] - row['target_x_collision_triples'] != total:
            raise ArithmeticError('exclusion sum mismatch')
        if sum(row['valid_mass_by_trace']) != total or min(row['valid_mass_by_trace']) < 0:
            raise ArithmeticError('trace mass mismatch')
        bound = Fraction(sum(min(a, b) for a, b in zip(row['fiber_capacity'], row['valid_mass_by_trace'])), n - 1)
        if decode(row['trace_probability_ceiling']) != bound:
            raise ArithmeticError('trace ceiling mismatch')
        if not bound <= decode(row['corrected_probability_ceiling']) <= decode(row['old_probability_ceiling']):
            raise ArithmeticError('bound ordering failure')
        if sum(row['fiber_capacity']) != n - 1:
            raise ArithmeticError('wrong target measure')
        batch = next(r for r in old if r['kind'] == 'batch' and key(r) == key(row))
        measured = batch['counters']['field_api']['totals']['multiplications']
        if row['batch_measured_multiplications'] != measured or row['batch_multiplication_floor'] > measured:
            raise ArithmeticError('S3 floor/measurement mismatch')
        batchChecks = [r for r in row['target_checks'] if r['cohort'].startswith('frozen-batch-')]
        if Counter((tuple(q['target']), 'frozen-batch-' + q['stratum']) for q in batch['queries']) != Counter((tuple(q['target']), q['cohort']) for q in batchChecks):
            raise ArithmeticError('batch input mismatch')
        for q in batchChecks:
            saved = next(r for r in batch['queries'] if r['target'] == q['target'])
            if q['exact_projected_relations'] != len(saved['solutions']):
                raise ArithmeticError('saved batch oracle disagreement')
        cold = [r for r in old if r['kind'] == 'instance' and key(r) == key(row)]
        coldChecks = [r for r in row['target_checks'] if r['cohort'] not in ('fresh-bound-holdout',) and not r['cohort'].startswith('frozen-batch-')]
        if Counter((tuple(q['target']), q['cohort']) for q in cold) != Counter((tuple(q['target']), q['cohort']) for q in coldChecks):
            raise ArithmeticError('cold input mismatch')
        for q in coldChecks:
            saved = next(r for r in cold if r['target'] == q['target'])
            if q['exact_projected_relations'] != len(saved['expected']):
                raise ArithmeticError('saved cold oracle disagreement')
        hybrid = sum(q['hybrid_complete_multiplication_floor'] for q in batchChecks)
        if hybrid != row['same_batch_hybrid_multiplication_floor']:
            raise ArithmeticError('hybrid complete-work floor mismatch')
        if decode(row['hybrid_floor_over_measured_s3_multiplications']) != Fraction(hybrid, measured):
            raise ArithmeticError('component ratio mismatch')
        targetChecks += len(row['target_checks'])
        fresh += sum(r['cohort'] == 'fresh-bound-holdout' for r in row['target_checks'])
    if (targetChecks, fresh) != (164, 60):
        raise ArithmeticError('missing larger/fresh bound checks')
    audit = {'spaces': len(spaces), 'tiny_target_space_checks': 265,
             'larger_target_checks': targetChecks, 'fresh_uniform_targets': fresh,
             'identity_checks': 2048, 'exhaustive_identity_target_checks': 756,
             'failures': 0, 'source_hashes_match': True,
             'raw_sha256': hashlib.sha256(rawPath.read_bytes()).hexdigest(),
             'same_batch_hybrid_floor_exceeds_measured_s3_multiplications': all(r['same_batch_hybrid_multiplication_floor'] > r['batch_measured_multiplications'] for r in spaces),
             'trace_fiber_improvement_beyond_chart_exclusion': any(decode(r['trace_probability_ceiling']) < decode(r['corrected_probability_ceiling']) for r in spaces),
             'larger_impossible_trace_classes': sum(len(r['impossible_trace_classes']) for r in spaces)}
    pruning = []
    for row in spaces:
        h, s = row['same_batch_hybrid_multiplication_floor'], row['batch_measured_multiplications']
        for costFraction in (Fraction(1), Fraction(4, 5)):
            removed = max(Fraction(0), 1 - costFraction * s / h)
            pruning.append({'n': row['n'], 'base_kind': row['base_kind'], 'd': row['d'],
                            'target_fraction_of_s3_multiplications': str(costFraction),
                            'necessary_removed_branch_fraction': str(removed),
                            'decimal': float(removed),
                            'assumptions': 'surviving branches retain the current seven-multiplication floor; zero pruning/setup/recovery overhead, so optimistic; complete enumeration only'})
    summary = {'audit': audit, 'spaces': spaces, 'thresholds': thresholds(), 'pruning_targets': pruning,
               'classification': 'accounting and conditional architecture bounds',
               'common_operation_speedup': None, 'full_dlp_S': None, 'rho_ratio': None}
    (HERE / 'audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    (HERE / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    render(summary)
    print(json.dumps(audit, indent=2))


def render(summary):
    rows = summary['spaces']
    lines = ['# Bounds identify coefficient branching as the next obstacle', '',
             'The existing hybrids cannot beat the measured direct S3 table in field-multiplication '
             'count on complete enumeration of any of the ten frozen eight-target batches: even '
             'their branch-only lower bounds exceed the table\'s measured total. The smallest gap '
             'is 15.73 times and the largest is 53.03 times. These are component bounds, not timed '
             'hybrid completions, total-operation speedups, or first-hit comparisons.', '',
             '[Proofs and assumptions](PROOFS.md) · [next experiment](NEXT.md) · [raw](raw.jsonl) · '
             '[summary](summary.json) · [audit](audit.json)', '',
             '## Exact matched-work comparison', '',
             'All columns below are field-multiplication calls for the same complete eight-target '
             'batch at the stated support. Hybrid entries are lower bounds excluding setup and '
             'recovery; S3 entries are measured totals including setup and verification. SAT '
             'work has no comparable completed multiplication measurement. A timeout has not '
             'been converted into a completed result.', '']
    htmlTables = []
    for n in (18, 30):
        selected = [r for r in rows if r['n'] == n]
        header = ['Variant'] + ['%s d%d' % (r['base_kind'], r['d']) for r in selected]
        table = []
        for variant in ('hybrid-reference', 'hybrid-filtered', 'pair-invariants-s3', 'chained-s3 SAT', 'symmetric-S4 SAT'):
            values = [variant]
            for r in selected:
                if variant.startswith('hybrid-'):
                    values.append('>= %s' % format(r['same_batch_hybrid_multiplication_floor'], ','))
                elif variant == 'pair-invariants-s3':
                    values.append(format(r['batch_measured_multiplications'], ','))
                else:
                    values.append('unmeasured')
            table.append(values)
        lines += ['%d-bit field:' % n, '', '| ' + ' | '.join(header) + ' |',
                  '|' + '|'.join(['---'] * len(header)) + '|']
        lines += ['| ' + ' | '.join(r) + ' |' for r in table] + ['']
        htmlTables.append('<div class="table-scroll"><table><caption>%d bits: field-multiplication calls for the same eight-target complete enumeration. Hybrid values are lower bounds; S3 values are measured.</caption><thead><tr>%s</tr></thead><tbody>%s</tbody></table></div>' %
                          (n, ''.join('<th>' + x + '</th>' for x in header),
                           ''.join('<tr>' + ''.join('<td>' + x.replace('>=', '&ge;') + '</td>' for x in r) + '</tr>' for r in table)))
    lines += ['| Bits | Base | d | Hybrid floor / measured S3 multiplication count | Measured S3 / its conservative floor |',
              '|---|---|---|---|---|']
    for r in rows:
        lines.append('| %d | %s | %d | %.4f | %.4f |' % (r['n'], r['base_kind'], r['d'],
                     decode(r['hybrid_floor_over_measured_s3_multiplications']), decode(r['batch_measured_over_floor'])))
    lines += ['', 'Every row is an accounting/architecture diagnostic. The multiplication bounds '
              'leave additions, squarings, binary work, allocation and control outside this one '
              'axis. They do not provide a calibrated common-operation speedup or rho ratio.', '',
              '## Yield ceilings', '',
              'The F4-defined curve orders are exactly 262,926 at 18 bits and 1,073,781,414 at '
              '30 bits. Both are composite; these probabilities are for uniform affine targets '
              'in the full curve group. No prime-subgroup probability is inferred.', '',
              '| Bits | Base | d | K abscissas | Signed triple mass after exclusions | Success ceiling | Uniform attempts lower bound |',
              '|---|---|---|---|---|---|---|']
    for r in rows:
        attempts = decode(r['mean_uniform_attempts_lower_bound']) if r['mean_uniform_attempts_lower_bound'] else None
        lines.append('| %d | %s | %d | %d | %d | %.6f%% | %s |' %
                     (r['n'], r['base_kind'], r['d'], r['K_abscissas'], r['valid_signed_triples'],
                      100 * decode(r['trace_probability_ceiling']), 'infinite' if attempts is None else '>= %.3f' % attempts))
    lines += ['', 'These are ceilings, not observed success rates. They correct the earlier signed '
              'triple count for infinity and target-abscissa collisions. Signed multiplicity '
              'differs from projected x-tuples; retaining multiplicity is essential for the proof.', '',
              '**Negative trace result:** every one of the 54 trace classes has positive '
              'admissible triple mass on every larger base. The trace-fiber ceiling equals the '
              'corrected global first-moment ceiling in all ten cases. Thus this quotient adds '
              'no stronger global success ceiling or whole-target rejection on this panel. '
              'It does not prove that finer coefficient constraints are useless.', '',
              '## Necessary scaling thresholds', '',
              'To make even the loose counting ceiling reach 50%, one needs 8*C(K,3) >= '
              '(N-1)/2. The dimensions below assume as optimistically as possible that every '
              'nonzero abscissa in V lifts. Real point support, exclusions and sum collisions '
              'can require more. Larger rows are exact derived thresholds, not solver experiments.', '',
              '| Field bits | Necessary K for 50% ceiling | Necessary binary dimension | Current S3 entries at that K |',
              '|---|---|---|---|']
    for r in summary['thresholds']:
        lines.append('| %d | %s | %d | %s |' % (r['n'], format(r['minimum_K_for_half_probability_ceiling'], ','),
                     r['necessary_binary_dimension'], format(r['minimum_current_s3_entries'], ',')))
    lines += ['', '## Architecture bound and next step', '',
              'Under independent uniform targets, complete base scanning per query, and a '
              'requirement for a constant multiple of M independent relation rows, the explicit '
              'pair-table architecture incurs Omega(M^2 + N/M), minimized at Omega(N^(2/3)). '
              'This is a conditional bound for that architecture. It is not a universal '
              'index-calculus lower bound, a measured exponent fit, or a full-DLP rho comparison.', '',
              'The next proposal is **coefficient-block rejection before branch solving**. '
              'The general-coefficient normalization in NEXT.md simplifies the target constant '
              'to B/r^2 and exposes a restricted linear equation. The existing fixed-branch '
              'pullback alone still visits a quadratic number of branches. A block relaxation must '
              'prove an entire block inconsistent at less cost than visiting it; whether that '
              'happens is unproved. No such block solver was implemented in this round.', '',
              'On the 30-bit prefix d10 batch, surviving branches with the present seven-call '
              'multiplication floor would require at least 97.6259% pruning to match the S3 '
              'multiplication total, or 98.1007% pruning to be 20% below it. These necessary '
              'thresholds optimistically assign zero cost to pruning, setup and recovery. '
              'A different per-branch formula requires a new bound. They are not a forecast '
              'that such pruning will succeed; exact rational thresholds are in summary.json.', '',
              '## Verification and reproduction', '',
              'Exact signed-triple enumeration and instrumented original branch counts agree '
              'on 265 tiny target/space cases, including three fresh spaces. The larger audit '
              'covers all 24 cold targets, 80 batch targets and 60 new uniform targets. The '
              'normalization passes 2,048 exact residual-identity checks plus 756 exhaustive '
              'nonzero-target constant checks across all F4 coefficient choices at GF64. '
              'All source/input hashes match and there are zero failures.', '',
              '`python research/nagao_relations/bounds_01/analyze.py` rechecks saved arithmetic, '
              'provenance, input coverage and prior-output agreement and regenerates these '
              'tables. To repeat expensive oracle/branch checks, run run.py from a separate '
              'checkout of its frozen source commit in publication.json; evidence files refuse '
              'overwrite. The identity checker has its own frozen source mapping. No solver '
              'performance candidate, timing comparison or full-DLP run was added here.', '']
    (HERE / 'RESULTS.md').write_text('\n'.join(lines))
    panel = ('  <!-- NAGAO BOUNDS BEGIN -->\n'
             '  <div class="panel" id="nagao-bound-audit"><div class="panel-head"><h2>Bounds require skipping coefficient branches</h2>'
             '<p>Source: <code>research/nagao_relations/bounds_01/summary.json</code>. '
             'On all ten matched eight-target batches, the hybrid branch-only multiplication '
             'floor exceeds the measured direct S3 total by 15.73–53.03 times. These are '
             'complete-enumeration component bounds, not timed completions or first-hit claims. '
             'Classification: accounting and conditional architecture bound. Common-cost '
             'speedup, S, rho and full-DLP floor ratios remain unmeasured.</p></div>' + ''.join(htmlTables) +
             '<div class="panel-head"><p>All 54 trace fibers have positive mass in every larger '
             'base: no whole-target rejection from this quotient. At 30 bits, the prefix d10 '
             'success ceiling is 15.896911%, not a measured success rate. Under uniform targets, '
             'full base scans and a linear number of required relation rows, explicit pair '
             'setup plus queries has the conditional lower bound Omega(M^2 + N/M). '
             'The next hypothesis is rejecting coefficient blocks before solving individual '
             'branches. Old timing and completion panels remain historical evidence.</p></div></div>\n'
             '  <!-- NAGAO BOUNDS END -->\n\n')
    scoreboard = ROOT / 'docs/index-calculus-scoreboard.html'
    text = scoreboard.read_text()
    if '  <!-- NAGAO BOUNDS BEGIN -->' in text:
        if text.count(panel) != 1:
            raise ArithmeticError('existing bound panel differs from evidence')
    else:
        marker = '  <!-- ============ NOTES ============ -->'
        if text.count(marker) != 1:
            raise ArithmeticError('scoreboard marker changed')
        scoreboard.write_text(text.replace(marker, panel + marker))


if __name__ == '__main__':
    main()
