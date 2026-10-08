"""Replay certificates and publish literal tables from frozen measurements."""
import hashlib
import json
from pathlib import Path
import sys

import run as experiment

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def check_certificates(row):
    if sorted(c['xs'] for c in row['certificates']) != sorted(row['solutions']):
        raise ArithmeticError('certificate coverage mismatch')
    if len({tuple(x) for x in row['solutions']}) != len(row['solutions']):
        raise ArithmeticError('duplicate output relation')
    if any(len(c['points']) != 3 for c in row['certificates']):
        raise ArithmeticError('certificate arity mismatch')


def main():
    rawPath = HERE / 'raw.jsonl'
    raw = [json.loads(s) for s in rawPath.read_text().splitlines()]
    contract = json.loads((HERE / 'contract.json').read_text())
    trials = [r for r in raw if r['kind'] == 'trial']
    batches = [r for r in raw if r['kind'] == 'batch']
    instances = {r['instance_sha256']: r for r in raw if r['kind'] == 'instance'}
    if len(trials) != 240 or len(batches) != 10 or len(instances) != 24:
        raise ArithmeticError('campaign incomplete')
    publication = json.loads((HERE / 'publication.json').read_text())
    if raw[0]['commit'] != publication['local_frozen_source_commit']:
        raise ArithmeticError('published source mapping mismatch')
    for name, digest in raw[0]['sha256'].items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != digest:
            raise ArithmeticError('frozen source changed: ' + name)
    oracles = {}
    for row in raw:
        if row['kind'] == 'space':
            f, c, v = experiment.context(row, counted=False)
            if row['base_kind'] == 'f4-stable':
                meta = v.metadata(c.a2)
                if not meta['f4_linear'] or not meta['frobenius4_stable'] or meta['proper_subfield_containments']:
                    raise ArithmeticError('structured support invariant failed')
            oracles[row['n'], row['base_kind'], row['d']] = experiment.Oracle(f, c, v)
    certificates = 0
    for digest, item in instances.items():
        instance = {k: v for k, v in item.items() if k not in ('kind', 'instance_sha256', 'expected')}
        if hashlib.sha256(json.dumps(instance, sort_keys=True).encode()).hexdigest() != digest:
            raise ArithmeticError('instance hash mismatch')
        rr = [r for r in trials if r['instance_sha256'] == digest]
        if sorted((r['variant'], r['mode']) for r in rr) != sorted(
                (v, mode) for v in contract['variants'] for mode in ('first', 'enumerate')):
            raise ArithmeticError('missing or duplicate matched variant')
        oracle = oracles[item['n'], item['base_kind'], item['d']]
        expected = oracle.expected(tuple(oracle.f.fromCoords(x) for x in item['target']))
        if expected != {tuple(x) for x in item['expected']}:
            raise ArithmeticError('recorded oracle set changed on replay')
        for row in rr:
            for key in ('n', 'd', 'base_kind', 'basis', 'coefficients', 'target', 'stratum', 'cohort'):
                if row[key] != item[key]:
                    raise ArithmeticError('cross-variant instance mismatch')
            experiment.check(row, item, oracle)
            check_certificates(row)
            certificates += len(row['certificates'])
            if abs(sum(row['phase_seconds'].values()) - row['all_phase_seconds']) > 1e-6:
                raise ArithmeticError('exclusive timer sum mismatch')
    old = [json.loads(s) for s in (HERE.parent / 'subfield_01/raw.jsonl').read_text().splitlines()]
    prior = {(r['n'], tuple(r['target'])) for r in old if r['kind'] == 'instance' and r['d'] == 8}
    replayed = {(r['n'], tuple(r['target'])) for r in instances.values() if r['cohort'].startswith('frozen-')}
    if prior != replayed:
        raise ArithmeticError('frozen d8 target coverage mismatch')
    counterexample = None
    batchSummary = []
    for row in batches:
        if row['batch_targets'] != 8 or len(row['queries']) != 8:
            raise ArithmeticError('incorrect amortization denominator')
        if row['all_phase_seconds'] > contract['batch_budget_seconds']:
            raise ArithmeticError('batch did not complete within its named budget')
        oracle = oracles[row['n'], row['base_kind'], row['d']]
        if row['table']['signed_base_size'] != 2 * len(oracle.base):
            raise ArithmeticError('table and independent oracle support differ')
        for query in row['queries']:
            experiment.check(query, {**row, 'target': query['target']}, oracle)
            check_certificates(query)
            certificates += len(query['certificates'])
            if counterexample is None and row['base_kind'] == 'f4-stable' and query['certificates']:
                f, c = oracle.f, oracle.c
                target = tuple(f.fromCoords(x) for x in query['target'])
                conjugate = tuple(f.frob(x, 2) for x in target)
                if conjugate != target:
                    total = None
                    vv = experiment.space.Space(f, row['basis'])
                    points = query['certificates'][0]['points']
                    transformed = []
                    for point in points:
                        pp = tuple(f.frob(f.fromCoords(x), 2) for x in point)
                        if pp[0] not in vv.indices or not c.onCurve(pp):
                            raise ArithmeticError('Frobenius support failure')
                        transformed.append([f.toCoords(x) for x in pp])
                        total = c.add(total, pp)
                    if total != conjugate or total == target:
                        raise ArithmeticError('fixed-target counterexample failed')
                    counterexample = {'n': row['n'], 'd': row['d'], 'basis': row['basis'],
                                      'target': query['target'], 'original_points': points,
                                      'conjugate_target': [f.toCoords(x) for x in conjugate],
                                      'conjugated_points': transformed,
                                      'conjugated_relation_is_for_original_target': False}
        querySeconds = sum(q['all_phase_seconds'] for q in row['queries'])
        relations = sum(len(q['solutions']) for q in row['queries'])
        strata = {s: {'targets': sum(q['stratum'] == s for q in row['queries']),
                      'targets_with_relations': sum(q['stratum'] == s and bool(q['solutions']) for q in row['queries']),
                      'verified_relations': sum(len(q['solutions']) for q in row['queries'] if q['stratum'] == s)}
                  for s in ('uniform', 'known_decomposable')}
        batchSummary.append({k: row[k] for k in ('n', 'd', 'base_kind', 'batch_targets', 'setup_seconds', 'all_phase_seconds', 'table', 'counters')} |
                            {'query_seconds': querySeconds, 'verified_relations': relations, 'strata': strata,
                             'setup_plus_queries_seconds': row['setup_seconds'] + querySeconds,
                             'external_validation_and_orchestration_seconds': row['all_phase_seconds'] - row['setup_seconds'] - querySeconds})
    if counterexample is None:
        raise ArithmeticError('no saved fixed-target Frobenius counterexample')
    groups = []
    keys = sorted({(r['n'], r['base_kind'], r['d'], r['variant'], r['mode']) for r in trials})
    for n, base, d, variant, mode in keys:
        rr = [r for r in trials if (r['n'], r['base_kind'], r['d'], r['variant'], r['mode']) == (n, base, d, variant, mode)]
        groups.append({'n': n, 'base_kind': base, 'd': d, 'variant': variant, 'mode': mode,
                       'attempts': len(rr), 'resolved': sum(r['status'] != 'timeout' and r['within_budget'] for r in rr),
                       'verified_relations': sum(r['verified_unique_relations'] for r in rr),
                       'all_phase_seconds': sum(r['all_phase_seconds'] for r in rr),
                       'common_operations': None, 'full_dlp_S': None, 'rho_ratio': None})
    paired = []
    for digest, item in instances.items():
        for mode in ('first', 'enumerate'):
            rr = {r['variant']: r for r in trials if r['instance_sha256'] == digest and r['mode'] == mode}
            ref = rr['hybrid-reference']
            for variant, candidate in rr.items():
                if variant == 'hybrid-reference':
                    continue
                admissible = ref['status'] != 'timeout' and candidate['status'] != 'timeout' and ref['within_budget'] and candidate['within_budget']
                if mode == 'enumerate' and admissible and ref['solutions'] != candidate['solutions']:
                    raise ArithmeticError('paired complete outputs differ')
                a, b = ref['counters']['field_api']['totals'], candidate['counters']['field_api']['totals']
                paired.append({'instance_sha256': digest, 'n': item['n'], 'd': item['d'], 'base_kind': item['base_kind'],
                               'variant': variant, 'mode': mode, 'equal_completed_work': admissible,
                               'reference_status': ref['status'], 'candidate_status': candidate['status'],
                               'reference_field_primitives': a, 'candidate_field_primitives': b,
                               'candidate_binary_word_operations': candidate['counters']['binary_word_operations'],
                               'field_api_sum_ratio_candidate_over_reference': b['fieldOperations'] / a['fieldOperations']
                               if admissible and variant in ('hybrid-filtered', 'pair-invariants-s3') else None,
                               'common_operation_speedup': None})
    profile = {}
    for n in (18, 30):
        rr = [r for r in old if r['kind'] == 'supplemental' and r['n'] == n]
        total = sum(r['all_phase_seconds'] for r in rr)
        profile[n] = {'setup_fraction': sum(r['phase_seconds']['setup'] for r in rr) / total,
                      'support_verification_fraction': sum(r['phase_seconds']['support_extract_verify'] for r in rr) / total}
    audit = {'matched_trials': len(trials), 'instances': len(instances), 'batches': len(batches),
             'batch_queries': sum(len(r['queries']) for r in batches), 'certificates_replayed': certificates,
             'source_hashes_match': True, 'frozen_targets_replayed': len(replayed),
             'validation_failures': 0, 'frobenius_counterexample': counterexample,
             'raw_sha256': hashlib.sha256(rawPath.read_bytes()).hexdigest()}
    (HERE / 'audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    summary = {'groups': groups, 'batches': batchSummary, 'audit': audit, 'prior_profile': profile,
               'classification': 'engineering diagnostic', 'common_operation_speedup': None,
               'full_dlp_S': None, 'rho_ratio': None}
    (HERE / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (HERE / 'comparison.json').write_text(json.dumps({'pairs': paired, 'limits': 'API ratios are descriptive only. Word/SAT work and complete ECDLP lack common calibration; timed-out or unequal pairs have null ratios.'}, indent=2) + '\n')
    render(raw, summary)
    print(json.dumps({'audit': {k: v for k, v in audit.items() if k != 'frobenius_counterexample'},
                      'totals': {variant: {'resolved': sum(r['status'] != 'timeout' and r['within_budget'] for r in trials if r['variant'] == variant),
                                          'attempted': sum(r['variant'] == variant for r in trials)} for variant in contract['variants']}}, indent=2))


def render(raw, summary):
    trials = [r for r in raw if r['kind'] == 'trial']
    overview = ['| Variant | First-mode resolved / 24 | Exact enumeration / 24 | Common-cost speedup | Correctness failures |',
                '|---|---|---|---|---|']
    for variant in ('hybrid-reference', 'hybrid-filtered', 'pair-invariants-s3', 'chained-s3', 's4-symmetric'):
        counts = [sum(r['variant'] == variant and r['mode'] == mode and r['status'] != 'timeout' and r['within_budget']
                      for r in trials) for mode in ('first', 'enumerate')]
        overview.append('| %s | %d/24 | %d/24 | unmeasured | 0 |' % (variant, *counts))
    tables = []
    htmlTables = []
    for n in (18, 30):
        for mode in ('first', 'enumerate'):
            columns = [('prefix', 8), ('prefix', 9), ('prefix', 10), ('f4-stable', 8), ('f4-stable', 10)]
            header = ['Variant'] + ['%s d%d' % p for p in columns]
            rows = []
            for variant in ('hybrid-reference', 'hybrid-filtered', 'pair-invariants-s3', 'chained-s3', 's4-symmetric'):
                row = [variant]
                for base, d in columns:
                    g = next(g for g in summary['groups'] if (g['n'], g['base_kind'], g['d'], g['variant'], g['mode']) == (n, base, d, variant, mode))
                    row.append('%d/%d' % (g['resolved'], g['attempts']))
                rows.append(row)
            tables += ['%d-bit field, %s; resolved within three seconds / attempts:' % (n, mode), '',
                       '| ' + ' | '.join(header) + ' |', '|' + '|'.join(['---'] * len(header)) + '|']
            tables += ['| ' + ' | '.join(r) + ' |' for r in rows] + ['']
            htmlTables.append('<div class="table-scroll"><table><caption>%d bits, %s: resolved within three seconds / attempts.</caption><thead><tr>%s</tr></thead><tbody>%s</tbody></table></div>' %
                              (n, mode, ''.join('<th>' + h + '</th>' for h in header),
                               ''.join('<tr>' + ''.join('<td>' + x + '</td>' for x in r) + '</tr>' for r in rows)))
    lines = ['# Structured-base scaling: the direct S3 table changes the comparison', '',
             'Executed 240 matched cold trials, including all eight frozen d8 targets and sixteen '
             'fresh holdouts, plus ten charged eight-target enumeration batches. All saved answers '
             'and complete sets pass independent replay. The pair-invariant method is a direct '
             'Semaev solver and is therefore a stronger comparator to the hybrid; a win for it '
             'cannot be claimed as a function-first advantage.', '',
             '[Proof and contract](README.md) · [raw evidence](raw.jsonl) · [summary](summary.json) · '
             '[paired counter comparison](comparison.json) · [audit](audit.json)', '',
             'The matched three-second cold results are:', ''] + overview + ['',
             'The table is an engineering diagnostic. First-mode resolutions include proofs that '
             'a target has no admissible relation. At 18 bits the hybrid resolves more first-mode '
             'trials (11/12) than the table (8/12); the table advantage in complete enumeration '
             'does not establish universal first-hit superiority. All d10 cold table trials '
             'time out while building the table.', '',
             '## Matched cold trials', '',
             'First-relation resolutions include proved empty targets; enumeration resolves only '
             'after the exact whole set is returned. Timeout is unknown and partial enumeration '
             'remains incomplete. The three-second budget includes supplied-basis expansion, '
             'curve and solver setup, search, extraction and verification. SAT soft overruns are retained.', ''] + tables
    lines += ['## Charged pair-table batches', '',
              'Each row builds a new S3 table and enumerates eight new targets: four uniform and '
              'four sampled from signed factor triples. No hybrid or SAT batch comparison was run, '
              'so these are capacity/amortization diagnostics. Setup remains quadratic. The batch '
              'wall clock also includes the independent harness oracle checks between queries; '
              'their time plus orchestration is exposed separately in summary.json. The table '
              'below reports the measured setup and sum of timed queries, which both include solver-side '
              'point verification. Counts are projected relations per target, not matrix-independent rows.', '',
              '| Bits | Base | d | Signed points | S3 table entries | Setup s | Eight queries s | Batch wall s | Uniform relations | Supported relations |',
              '|---|---|---|---|---|---|---|---|---|---|']
    for b in summary['batches']:
        lines.append('| %d | %s | %d | %d | %d | %.4f | %.4f | %.4f | %d | %d |' %
                     (b['n'], b['base_kind'], b['d'], b['table']['signed_base_size'], b['table']['table_entries'],
                      b['setup_seconds'], b['query_seconds'], b['all_phase_seconds'],
                      b['strata']['uniform']['verified_relations'], b['strata']['known_decomposable']['verified_relations']))
    a = summary['audit']
    lines += ['', '## What the experiments establish', '',
              '- Both structured dimensions are real F4-linear, Frobenius4-stable subspaces; they '
              'are not contained in a proper subfield. Structured d9 is impossible because F4-linearity '
              'forces even binary dimension. Measured signed-point support varies substantially: '
              'at 18 bits and d8 the structured base has 148 points versus 254 for the prefix base. '
              'Their timing difference cannot be interpreted as a fixed-support speedup.',
              '- The early filter preserves exact function acceptance on 10,638 candidate checks. '
              'All 1,920 parity-membership checks agree with Gaussian image recovery. All five '
              'solvers agree with exhaustive signed triples on 106 tiny target/space pairs.',
              '- Independent replay verifies %d signed-point certificates and every larger complete '
              'oracle set. Frozen targets and source hashes match. No correctness failures were found.' % a['certificates_replayed'],
              '- The saved Frobenius counterexample shows that conjugated factors remain in the '
              'structured base but sum to R^4, not the original R. Stable support does not justify '
              'fixed-target orbit collapse.',
              '- Ordinary setup caching was overestimated in the proposed plan: previous d8 setup '
              'was only 0.24455%% and 0.20975%% of wall time. The new pair table is different: it deliberately '
              'pays quadratic precomputation once per named batch.' .replace('%%', '%'), '',
              'The counting ceiling min(1,8*C(M/2,3)/(#E-1)) changes when the base changes. No '
              'subquadratic cold-search theorem or exponent fit follows from these cells. The direct '
              'S3 table uses quadratic setup; the hybrid still enumerates O(2^(2d)) branches. '
              'Any lower field-API sum for the parity filter must also account for its new binary '
              'word operations. comparison.json retains both vectors and makes incomplete-pair '
              'ratios null.', '',
              '## Hypothesis decisions', '',
              '| Proposal | Observation or proof | Decision |',
              '|---|---|---|',
              '| Cache ordinary image setup for a 20% gain | Even removing all prior setup saves under 0.25% | Disproved for this measured workload |',
              '| Early parity rejection improves completion | Same 13/24 first resolutions and 0/24 complete enumerations as the reference | Correct filter; completion hypothesis not supported |',
              '| Stable support permits fixed-target Frobenius pruning | Saved conjugated relation sums to a different target | Disproved without a target stabilizer |',
              '| Pair-invariant S3 strengthens enumeration | 14/24 cold completions versus 0/24 for both hybrids | Retain as a stronger Semaev baseline |',
              '| Charged reuse reaches larger bases | All ten eight-target batches complete, including d10 at 30 bits | Capacity established; setup remains quadratic |',
              '| Larger structured bases guarantee useful random-target yield | Only one of twenty uniform 30-bit batch targets has a relation | Yield remains an obstacle; supported targets do not estimate it |', '',
              'The filter uses fewer field-API calls on ten of the thirteen matched completed '
              'first-relation searches, and more on three. Its candidate/reference API-sum ratios '
              'range from 0.5037 to 1.1799; these exclude its separately recorded binary word work '
              'and are not common-operation speedups. Regressions remain in comparison.json.', '',
              'The next function-first experiment must reduce coefficient branches and compare '
              'against the direct S3 table on identical cold and named-batch workloads. Improving '
              'the early filter alone leaves O(2^(2d)) coefficient branching intact. Random-target '
              'yield and first-relation latency must remain separate from complete enumeration.', '',
              'Classification: engineering diagnostic. Calibrated common-operation speedup, '
              'full-DLP S, rho/floor ratios and the original 20% goal remain unproved. These small '
              'family-specific trials include frozen replays and fresh holdouts, but do not pass '
              'the broader 60-input/repetition and full-pipeline promotion gates. No full ECDLP '
              'pipeline or scalar recovery was run in this campaign.', '',
              '## Replay', '',
              'Run `python research/nagao_relations/structured_01/analyze.py` to rebuild the '
              'independent group-law oracles, check every saved result and certificate, and '
              'regenerate the report. An identical existing scoreboard panel is accepted; a '
              'different panel is rejected. This analysis is outside the measured budget.', '',
              'For a timing rerun, use a separate checkout of the published frozen source commit '
              'in publication.json, with Python 3.12 and pycryptosat 5.14.7. The measured runner '
              'refuses to overwrite raw.jsonl. The analysis script and reports were added after '
              'the campaign and are not part of its frozen measured source.', '']
    (HERE / 'RESULTS.md').write_text('\n'.join(lines))
    panel = ('  <div class="panel" id="nagao-structured-scaling"><div class="panel-head">'
             '<h2>Structured bases and a stronger direct S3 comparator</h2>'
             '<p>Frozen evidence: <code>research/nagao_relations/structured_01/summary.json</code> and '
             '<code>raw.jsonl</code>. 240 matched cold trials and 80 complete batched queries; '
             '%d certificates replayed, zero correctness failures. The pair-invariant table is a '
             'Semaev solver. Its gains cannot be relabelled as a function-first advantage. '
             'All rows are engineering diagnostics; common-cost speedup, S, rho and floor ratios '
             'remain unmeasured.</p></div>%s<div class="panel-head"><p>Structured d9 is impossible '
             'over F4. Ordinary prior setup was under 0.25%% of wall time. The direct S3 table '
             'amortizes quadratic setup over eight named targets; those batches are unmatched '
             'capacity checks. No subquadratic cold scaling or full-DLP improvement is established. '
             'Prior panels remain historical measurements.</p></div></div>\n\n' %
             (a['certificates_replayed'], ''.join(htmlTables)))
    scoreboard = ROOT / 'docs/index-calculus-scoreboard.html'
    text = scoreboard.read_text()
    if 'id="nagao-structured-scaling"' in text:
        if text.count(panel) == 1:
            return
        raise ValueError('scoreboard result panel differs from frozen evidence')
    marker = '  <!-- ============ NOTES ============ -->'
    if text.count(marker) != 1:
        raise ArithmeticError('scoreboard insertion marker changed')
    scoreboard.write_text(text.replace(marker, panel + marker))


if __name__ == '__main__':
    main()
