"""Independent replay and literal report for the conditional batch."""
from collections import Counter
import importlib.util
import json
import experiment as exp

HERE, ROOT = exp.HERE, exp.ROOT
spec = importlib.util.spec_from_file_location('block_analysis_helpers', HERE / 'analyze.py')
helpers = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helpers)
raw_bytes, dump, field, words, table = [getattr(helpers, name) for name in ('raw_bytes', 'dump', 'field', 'words', 'table')]


def run():
    raw = [json.loads(s) for s in raw_bytes('batch_raw.jsonl').splitlines()]
    contract = json.loads((HERE / 'batch_contract.json').read_text())
    items = json.loads((HERE / 'batch_targets.json').read_text())
    rows = [r for r in raw if r['kind'] == 'batch']
    if Counter(r['variant'] for r in rows) != Counter(contract['variants']):
        raise ArithmeticError('batch missing a variant')
    if next(r for r in raw if r['kind'] == 'targets')['items'] != items:
        raise ArithmeticError('batch target mismatch')
    publication = json.loads((HERE / 'publication.json').read_text())
    if raw[0]['commit'] != publication['batch_local_source_commit']:
        raise ArithmeticError('batch source publication mismatch')
    if raw[0]['cold_raw_sha256'] != exp.hashlib.sha256(raw_bytes()).hexdigest():
        raise ArithmeticError('cold selection evidence differs')
    for path, digest in raw[0]['sha256'].items():
        if exp.digest(ROOT / path) != digest:
            raise ArithmeticError('batch source/input changed')
    # Also enforce inherited cold dependency hashes when replaying this audit alone.
    cold_provenance = json.loads(raw_bytes().splitlines()[0])
    for path, digest in cold_provenance['sha256'].items():
        if exp.digest(ROOT / path) != digest:
            raise ArithmeticError('inherited dependency changed')
    f, c, v = exp.old.context(items[0], False)
    oracle = exp.old.Oracle(f, c, v)
    checkers = [exp.Checker(item) for item in items]
    relations = rejections = 0
    for row in rows:
        if len(row['queries']) != len(items):
            raise ArithmeticError('wrong amortization denominator')
        for item, query, checker in zip(items, row['queries'], checkers):
            if query['target'] != item['target'] or query['stratum'] != item['stratum']:
                raise ArithmeticError('batch input order changed')
            if query['expected_count'] != exp.check(query, item, oracle, checker):
                raise ArithmeticError('batch expected count differs')
            if len(query['solutions']) != len({tuple(x) for x in query['solutions']}):
                raise ArithmeticError('duplicate batch output')
            if any(len(cert['points']) != 3 for cert in query['certificates']):
                raise ArithmeticError('batch relation arity')
            relations += len(query['certificates'])
            rejections += len(query['rejection_certificates'])
        total = row['setup_seconds'] + row['query_seconds'] + row['orchestration_seconds']
        if abs(total - row['all_phase_seconds']) > 1e-6 or row['orchestration_seconds'] < -1e-6:
            raise ArithmeticError('batch timing not exclusive')
        counts = row['counters']['field_api']
        for name, value in counts['totals'].items():
            if value != sum(phase[name] for phase in counts['phases'].values()):
                raise ArithmeticError('batch counter sum differs')
    ref = next(r for r in rows if r['variant'] == 'pair-invariants-s3')
    summary_rows = []
    for row in rows:
        equal = row['status'] == ref['status'] == 'complete' and row['within_budget'] and ref['within_budget']
        if equal and [q['solutions'] for q in row['queries']] != [q['solutions'] for q in ref['queries']]:
            raise ArithmeticError('complete batches differ')
        count = sum(len(q['solutions']) for q in row['queries'])
        summary_rows.append({'variant': row['variant'], 'status': row['status'], 'within_budget': row['within_budget'],
                             'complete_queries': sum(q['status'] == 'complete' for q in row['queries']),
                             'verified_relations': count, 'field_primitives': field(row), 'selected_binary_words': words(row),
                             'all_phase_seconds': row['all_phase_seconds'], 'setup_seconds': row['setup_seconds'],
                             'query_seconds': row['query_seconds'], 'orchestration_seconds': row['orchestration_seconds'],
                             'equal_complete_work_to_s3': equal,
                             'field_api_ratio_to_s3': field(row)['fieldOperations'] / field(ref)['fieldOperations'] if equal else None,
                             'field_api_per_verified_relation': field(row)['fieldOperations'] / count if equal and count else None,
                             'common_operation_speedup': None, 'rho_ratio': None, 'floor_ratio': None})
    audit = {'batches': len(rows), 'matched_queries': len(rows) * len(items), 'fresh_targets': len(items),
             'relation_certificates_replayed': relations, 'rejection_certificates_replayed': rejections,
             'correctness_failures': 0, 'hashes_match': True,
             'raw_sha256': exp.hashlib.sha256(raw_bytes('batch_raw.jsonl')).hexdigest()}
    # Constructor-specific multiplication floor, independent of pruning success.
    d = contract['d']
    per_b = 3 * d * d + 8 * d + 3
    counted, curve, support = exp.old.context(items[0])
    r = counted.fromCoords(items[0]['target'][0])
    h, b = next((h, b) for h in support.values for b in curve.asRoots(counted.add(h, r)) if b)
    powers = [(x, counted.sqr(x), counted.frob(x, 2)) for x in support.basis]
    before = counted.report()['totals']['multiplications']
    exp.block.Circuit(counted, curve, support, r, b, h, powers)
    if counted.report()['totals']['multiplications'] - before != per_b:
        raise ArithmeticError('constructor multiplication formula differs')
    pruned = next(row for row in rows if row['variant'] == 'bilinear-blocks')
    coefficients = sum(q['stats']['coefficient_values'] for q in pruned['queries'])
    floor = coefficients * per_b
    if floor > field(pruned)['multiplications']:
        raise ArithmeticError('constructor floor exceeds measured work')
    s3_muls = field(ref)['multiplications']
    bounds = {'scope': 'Frozen bilinear Circuit constructor only, complete matched batch; not a universal bound.',
              'per_b_multiplications_formula': '3*d^2+8*d+3', 'per_b_multiplications': per_b,
              'coefficient_values': coefficients, 'constructor_multiplication_floor': floor,
              'measured_s3_total_multiplications': s3_muls, 'floor_over_s3_multiplications': floor / s3_muls,
              'necessary_constructor_saving_for_20pct_fewer_muls': 1 - .8 * s3_muls / floor,
              'constructor_counter_check_passed': True}
    summary = {'rows': summary_rows, 'audit': audit,
               'constructor_bound': bounds,
               'strata': {s: {'targets': sum(i['stratum'] == s for i in items),
                             'targets_with_relations': sum(q['stratum'] == s and q['expected_count'] > 0 for q in ref['queries'])}
                          for s in ('uniform', 'known_decomposable')},
               'classification': 'engineering diagnostic; common total cost uncalibrated'}
    dump('batch_summary.json', summary)
    dump('batch_audit.json', audit)
    field_rows, time_rows, word_rows = [], [], []
    for variant in contract['variants']:
        r = next(r for r in summary_rows if r['variant'] == variant)
        ff = r['field_primitives']
        ratio = r['field_api_ratio_to_s3']
        field_rows.append([variant, *[f'{ff[k]:,}' for k in ('additions', 'multiplications', 'squarings', 'fieldOperations')],
                           f'{ratio:.4f}' if ratio is not None else 'unmeasured', r['complete_queries'], r['verified_relations']])
        time_rows.append([variant, f'{r["setup_seconds"]:.4f}', f'{r["query_seconds"]:.4f}', f'{r["all_phase_seconds"]:.4f}'])
        word_rows.append([variant, f'{r["selected_binary_words"]:,}', 'unmeasured', 'unmeasured'])
    result = '\n'.join(['<!-- block-batch-results -->', '', '## Fresh eight-target batch with shared setup', '',
                        '[Batch contract](batch_contract.json) · [targets](batch_targets.json) · '
                        '[compressed raw](batch_raw.jsonl.gz) · [summary](batch_summary.json)', '',
                        'The cold screen selected n30, prefix d8 for this follow-up. Every variant gets the '
                        'same four fresh uniform and four fresh supported targets, with 60 seconds for the '
                        'entire eight-target enumeration batch. Field/curve/support setup is charged once '
                        'per batch; S3 pair tables and hybrid support images are also charged once. '
                        'Target-dependent coefficient circuits are rebuilt and charged per target.', '',
                        table(['Variant', 'Add calls', 'Mul calls', 'Square calls', 'Field API sum', 'API ratio to S3', 'Complete / 8', 'Verified relations'], field_rows), '',
                        'This table is in field API calls, not a calibrated machine-operation unit. '
                        'Complete workloads include empty targets and all verification. The following '
                        'separate counter tracks only the explicitly instrumented binary operations; '
                        'zero does not mean zero machine bit work. These components cannot be added '
                        'without a measured conversion.', '',
                        table(['Variant', 'Selected binary word operations', 'Calibrated total-cost ratio', 'Ratio to full-DLP floor'], word_rows), '',
                        'Wall time is a single observation, in seconds; no confidence interval or '
                        'unqualified runtime improvement is claimed:', '',
                        table(['Variant', 'Setup seconds', 'Queries seconds', 'All phases seconds'], time_rows), '',
                        f'The replay checks {relations} relation certificates and {rejections:,} block '
                        'rejection certificates, with zero failures. Uniform and supported outcomes: '
                        + json.dumps(summary['strata'], sort_keys=True) + '.', '',
                        'The original three-size calibrated cost goal remains unmet. Both the cold '
                        'regressions and this selected follow-up must accompany any claim about pruning.',
                        '', 'The frozen constructor alone requires ' + f'{floor:,}' +
                        ' multiplications, versus ' + f'{s3_muls:,}' +
                        ' measured total S3 multiplications on this batch. Even free perfect pruning '
                        'cannot remove that cost. See [the new bound and next hypothesis](BOUNDS_AND_NEXT.md).',
                        '<!-- block-batch-results-end -->', ''])
    (HERE / 'BATCH_RESULTS.md').write_text(result)
    path = HERE / 'RESULTS.md'
    text = path.read_text().split('<!-- block-batch-results -->')[0].rstrip()
    path.write_text(text + '\n\n' + result)
    # Literal batch values are also part of the canonical scoreboard deliverable.
    panel = '<!-- block-batch-scoreboard -->\n<div class="panel" id="nagao-block-batch"><div class="panel-head"><h2>Coefficient blocks versus a reused S3 table</h2><p>Source: <code>research/nagao_relations/blocks_01/batch_summary.json</code>. Eight identical fresh n30 d8 targets, setup and unsuccessful targets charged. API sums below are uncalibrated component counts; selected binary work is separate. Classification: engineering diagnostic, no full-DLP advance.</p></div><div class="table-scroll"><table><thead><tr><th>Variant</th><th>Field API calls</th><th>API ratio to S3</th><th>Complete queries / 8</th><th>Verified relations</th></tr></thead><tbody>'
    for r in field_rows:
        panel += '<tr>' + ''.join('<td>' + str(x) + '</td>' for x in [r[0], r[4], r[5], r[6], r[7]]) + '</tr>'
    panel += '</tbody></table></div><div class="panel-head"><p>Common-cost, rho and full-DLP floor ratios remain unmeasured. A passing pruning screen is not a Semaev crossover or the original three-size goal. Earlier cold regressions remain visible.</p></div></div>\n<!-- block-batch-scoreboard-end -->\n'
    (HERE / 'batch_scoreboard.html').write_text(panel)
    path = ROOT / 'docs/index-calculus-scoreboard.html'
    page = path.read_text()
    if '<!-- block-batch-scoreboard -->' in page:
        a, b = page.index('<!-- block-batch-scoreboard -->'), page.index('<!-- block-batch-scoreboard-end -->') + len('<!-- block-batch-scoreboard-end -->')
        page = page[:a] + panel.rstrip() + page[b:]
    else:
        page = page.replace('<!-- nagao-blocks-01-end -->', '<!-- nagao-blocks-01-end -->\n' + panel)
    path.write_text(page)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    run()
