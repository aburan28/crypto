"""Replay original, recovered and batch evidence; never price unresolved work."""
from collections import Counter
import gzip
import hashlib
import html
import json
import trial

HERE, ROOT, base = trial.HERE, trial.ROOT, trial.base


def raw(name):
    plain, packed = HERE / name, HERE / (name + '.gz')
    data = gzip.decompress(packed.read_bytes()) if packed.exists() else plain.read_bytes()
    if plain.exists() and plain.read_bytes() != data:
        raise ArithmeticError('compressed evidence differs: ' + name)
    return data


def dump(name, value):
    (HERE / name).write_text(json.dumps(value, indent=2) + '\n')


def field(row):
    return row['counters']['field_api']['totals'] if row['counters'] else None


def words(row):
    return sum(sum(p.values()) for p in row['counters']['binary_word_operations'].values()) if row['counters'] else None


def resolved(row):
    return row['status'] in ('first', 'complete') and row['within_budget'] and not row.get('prior_unresolved_runner_incident')


def counters(row):
    if row['counters'] is None:
        assert row['status'] == 'watchdog-timeout'
        return
    counts = row['counters']['field_api']
    for name, value in counts['totals'].items():
        assert value == sum(p[name] for p in counts['phases'].values()), name
    assert counts['totals']['fieldOperations'] == sum(counts['totals'][k] for k in ('additions', 'multiplications', 'squarings'))


def table(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |', '| ' + ' | '.join(['---'] * len(headers)) + ' |'] +
                     ['| ' + ' | '.join(map(str, row)) + ' |' for row in rows])


def main():
    names = ['raw.jsonl', 'continuation_raw.jsonl', 'batch_raw.jsonl']
    streams = {name: [json.loads(x) for x in raw(name).splitlines()] for name in names}
    original, continuation, batch = [streams[name] for name in names]
    contract = json.loads((HERE / 'contract.json').read_text())
    publication = json.loads((HERE / 'publication.json').read_text())
    assert original[0]['commit'] == publication['local_frozen_source_commit']
    for stream in (continuation, batch):
        assert stream[0]['commit'] == publication['recovered_run_source_commit']
    for stream in streams.values():
        for path, digest in stream[0]['sha256'].items():
            assert base.digest(ROOT / path) == digest, path
    assert continuation[0]['original_raw_sha256'] == hashlib.sha256(raw(names[0])).hexdigest()
    predecessor = gzip.decompress((ROOT / contract['predecessor']).read_bytes())
    assert hashlib.sha256(predecessor).hexdigest() == contract['predecessor_uncompressed_sha256']
    prior_items = [r for r in map(json.loads, predecessor.splitlines()) if r['kind'] == 'instance']
    items = json.loads((HERE / 'targets.json').read_text())
    signature = lambda x: (x['n'], x['d'], tuple(x['basis']), tuple(x['coefficients']), tuple(x['target']), x['stratum'])
    assert {signature(x) for x in prior_items} == {signature(x) for x in items if x['cohort'].startswith('frozen-')}
    trials = [r for r in original + continuation if r['kind'] == 'trial']
    assert len([r for r in original if r['kind'] == 'trial']) == 198
    assert len(trials) == contract['expected_trials'] == 288
    instances = {r['instance_sha256']: r for r in original + continuation if r['kind'] == 'instance'}
    assert len(instances) == len(items) == 24
    oracles, relations, rejections = {}, 0, 0

    def oracle_for(item):
        key = item['n'], tuple(item['basis']), tuple(item['coefficients'])
        if key not in oracles:
            f, c, v = base.old.context(item, False)
            oracles[key] = base.old.Oracle(f, c, v)
        return oracles[key]

    def audit(row, item, oracle, checker):
        nonlocal relations, rejections
        assert row['expected_count'] == base.check(row, item, oracle, checker)
        assert len(row['solutions']) == len({tuple(x) for x in row['solutions']})
        assert all(len(c['points']) == 3 for c in row['certificates'])
        relations += len(row['certificates'])
        rejections += len(row.get('rejection_certificates', []))

    for item in items:
        digest = item['instance_sha256']
        assert hashlib.sha256(json.dumps({k: v for k, v in item.items() if k != 'instance_sha256'}, sort_keys=True).encode()).hexdigest() == digest
        assert all(instances[digest][k] == v for k, v in item.items())
        oracle, checker = oracle_for(item), base.Checker(item)
        expected = oracle.expected(tuple(oracle.f.fromCoords(x) for x in item['target']))
        assert expected == {tuple(x) for x in instances[digest]['expected']}
        rr = [r for r in trials if r['instance_sha256'] == digest]
        assert Counter((r['variant'], r['mode']) for r in rr) == Counter((v, m) for v in contract['variants'] for m in contract['modes'])
        for row in rr:
            assert all(row[k] == v for k, v in item.items())
            audit(row, item, oracle, checker)
            counters(row)
            if row['phase_seconds'] is not None:
                assert abs(sum(row['phase_seconds'].values()) - row['all_phase_seconds']) < 1e-6
    cold_audit = {'trials': len(trials), 'targets': len(items), 'relation_certificates': relations,
                  'rejection_certificates': rejections, 'correctness_failures': 0,
                  'incident_excluded_cells': sum(bool(r.get('prior_unresolved_runner_incident')) for r in trials)}
    pairs = []
    for item in items:
        for mode in contract['modes']:
            rr = {r['variant']: r for r in trials if r['instance_sha256'] == item['instance_sha256'] and r['mode'] == mode}
            for candidate in ('reuse-circuit', 'reuse-batch-inverse'):
                for reference in ('bilinear-blocks', 'reuse-circuit', 'pair-invariants-s3'):
                    if candidate == reference:
                        continue
                    a, b = rr[reference], rr[candidate]
                    equal = resolved(a) and resolved(b)
                    if equal and mode == 'enumerate':
                        assert a['solutions'] == b['solutions']
                    pairs.append({'instance_sha256': item['instance_sha256'], 'n': item['n'], 'd': item['d'],
                                  'cohort': item['cohort'], 'stratum': item['stratum'], 'mode': mode,
                                  'candidate': candidate, 'reference': reference, 'equal_completed_work': equal,
                                  'candidate_field_primitives': field(b), 'reference_field_primitives': field(a),
                                  'candidate_selected_words': words(b), 'reference_selected_words': words(a),
                                  'field_api_ratio': field(b)['fieldOperations'] / field(a)['fieldOperations'] if equal else None,
                                  'common_operation_speedup': None})
    cold = [{'variant': v, **{m: {'resolved': sum(resolved(r) for r in trials if r['variant'] == v and r['mode'] == m),
                                 'attempts': len(items)} for m in contract['modes']}} for v in contract['variants']]
    groups = []
    axes = ('n', 'd', 'cohort', 'stratum', 'variant', 'mode')
    for key in sorted({tuple(r[k] for k in axes) for r in trials}):
        rr = [r for r in trials if tuple(r[k] for k in axes) == key]
        groups.append({**dict(zip(axes, key)), 'attempts': len(rr), 'resolved': sum(map(resolved, rr)),
                       'performed_field_api_calls': sum(field(r)['fieldOperations'] for r in rr if field(r)),
                       'cost_interpretation': 'Performed work only; not an equal-completion ratio.'})
    batch_items = json.loads((ROOT / contract['batch_targets']).read_text())
    assert base.digest(ROOT / contract['batch_targets']) == contract['batch_targets_sha256']
    assert next(r for r in batch if r['kind'] == 'targets')['items'] == batch_items
    rows = [r for r in batch if r['kind'] == 'batch']
    assert Counter(r['variant'] for r in rows) == Counter(contract['batch_variants'])
    checkers = [base.Checker(i) for i in batch_items]
    for row in rows:
        assert len(row['queries']) == len(batch_items)
        for q, item, checker in zip(row['queries'], batch_items, checkers):
            assert q['target'] == item['target'] and q['stratum'] == item['stratum']
            audit(q, item, oracle_for(item), checker)
        counters(row)
        assert abs(row['setup_seconds'] + row['query_seconds'] + row['orchestration_seconds'] - row['all_phase_seconds']) < 1e-6
        assert row['orchestration_seconds'] >= -1e-6
    ref = next(r for r in rows if r['variant'] == 'pair-invariants-s3')
    old = next(r for r in rows if r['variant'] == 'bilinear-blocks')
    summaries = []
    for row in rows:
        equal = row['status'] == ref['status'] == 'complete' and row['within_budget'] and ref['within_budget']
        if equal:
            assert [q['solutions'] for q in row['queries']] == [q['solutions'] for q in ref['queries']]
        if row['variant'].startswith('reuse-') and row['status'] == old['status'] == 'complete':
            assert [q['rejection_certificates'] for q in row['queries']] == [q['rejection_certificates'] for q in old['queries']]
            assert row['counters']['field_api']['phases']['block_rejection'] == old['counters']['field_api']['phases']['block_rejection']
        summaries.append({'variant': row['variant'], 'status': row['status'], 'complete_queries': sum(q['status'] == 'complete' for q in row['queries']),
                          'verified_relations': sum(len(q['solutions']) for q in row['queries']),
                          'field_primitives': field(row), 'selected_binary_words': words(row),
                          'field_api_phases': row['counters']['field_api']['phases'],
                          'peak_cached_field_slots': max((q.get('stats') or {}).get('cached_field_slots', 0) for q in row['queries']),
                          'seconds': row['all_phase_seconds'], 'equal_complete_work_to_s3': equal,
                          'field_api_ratio_to_s3': field(row)['fieldOperations'] / field(ref)['fieldOperations'] if equal else None,
                          'mul_ratio_to_s3': field(row)['multiplications'] / field(ref)['multiplications'] if equal else None,
                          'common_operation_speedup': None, 'full_dlp_S': None, 'rho_ratio': None, 'floor_ratio': None})
    rejection_floor = old['counters']['field_api']['phases']['block_rejection']['fieldOperations']
    bound = {'scope': 'Complete frozen eight-target batch, changes only to coefficient preparation with identical rejection traversal.',
             'unchanged_rejection_field_api_floor': rejection_floor,
             's3_total_field_api_calls': field(ref)['fieldOperations'],
             'floor_over_s3_field_api_calls': rejection_floor / field(ref)['fieldOperations'],
             'not_a_calibrated_cost_bound': True}
    audit_report = {'cold': cold_audit, 'batch_queries': len(rows)*len(batch_items),
                    'all_relation_certificates': relations, 'all_rejection_certificates': rejections,
                    'correctness_failures': 0, 'source_and_input_hashes_match': True,
                    'raw_sha256': {name: hashlib.sha256(raw(name)).hexdigest() for name in names}}
    summary = {'cold': cold, 'groups': groups, 'batch': summaries, 'bound': bound, 'audit': audit_report,
               'classification': 'engineering diagnostic; total operation unit uncalibrated', 'three_size_goal_met': False}
    dump('comparison.json', pairs)
    dump('summary.json', summary)
    dump('audit.json', audit_report)
    dump('batch_summary.json', {'rows': summaries, 'bound': bound})
    cold_table = table(['Variant', 'First / 24', 'Enumeration / 24'], [[r['variant'], r['first']['resolved'], r['enumerate']['resolved']] for r in cold])
    batch_table = table(['Variant', 'Add', 'Mul', 'Square', 'Field API sum', 'API ratio to S3', 'Complete / 8', 'Relations'],
                        [[r['variant'], *[f"{r['field_primitives'][k]:,}" for k in ('additions', 'multiplications', 'squarings', 'fieldOperations')],
                          f"{r['field_api_ratio_to_s3']:.4f}" if r['field_api_ratio_to_s3'] is not None else 'null', r['complete_queries'], r['verified_relations']] for r in summaries])
    diagnostic = table(['Variant', 'Selected binary word operations', 'Cached field slots', 'Observed seconds'],
                       [[r['variant'], f"{r['selected_binary_words']:,}", r['peak_cached_field_slots'], f"{r['seconds']:.4f}"] for r in summaries])
    (HERE / 'RESULTS.md').write_text('# Coefficient reuse: complete audit\n\n'
        'Engineering diagnostic. The original calibrated three-size goal remains unmet. '
        'All variants use identical F4-defined curves, prefix supports, targets and arity three.\n\n'
        '[Contract](contract.json) · [raw](raw.jsonl.gz) · [continuation](continuation_raw.jsonl.gz) · '
        '[batch raw](batch_raw.jsonl.gz) · [audit](audit.json) · [paired comparisons](comparison.json)\n\n'
        'Cold 3-second screen, 16 frozen and eight fresh targets at n18/n30, d8/d9. '
        'The unresolved original S3 runner incident is excluded from completion counts and cost pairs, '
        'even when replay succeeds. These counts are conservative and depend on the budget.\n\n' + cold_table + '\n\n'
        'Complete enumeration of the same eight n30 d8 targets as the predecessor batch: four uniform '
        'and four supported. Setup, failed targets, extraction and exact verification are charged.\n\n' + batch_table + '\n\n'
        'Field API calls are an uncalibrated equal-weight component sum. The next table keeps selected '
        'binary work separate; zero is not zero machine work. Cached slots are logical storage, not RSS. '
        'Timing is one observation, with no confidence interval or runtime claim.\n\n' + diagnostic + '\n\n'
        f"The audit replays {relations:,} relation certificates and {rejections:,} rejection certificates with zero failures. "
        'Complete reuse batches preserve the predecessor rejection certificates and rejection-phase counter vector exactly.\n\n'
        f"Changing coefficient preparation alone leaves a floor of {rejection_floor:,} rejection field API calls, "
        f"{bound['floor_over_s3_field_api_calls']:.4f} times the measured total S3 API count. "
        'This is a bound on this traversal and this unit, not a universal or calibrated lower bound. '
        'Cached fixed column spans are the next falsifiable change. See [execution notes](EXECUTION_NOTES.md).\n')
    panel = '<!-- nagao-reuse-01-start -->\n<div class="panel" id="nagao-reuse-01"><div class="panel-head"><h2>Coefficient reuse: rank tests remain dominant</h2><p>Engineering diagnostic. Source: <code>research/nagao_relations/reuse_01/summary.json</code>. Complete matched eight-target n30 d8 batch; all setup, empty targets and verification charged. Field API calls are uncalibrated; common-cost, full-DLP S and rho ratios remain null.</p></div><div class="table-scroll"><table><thead><tr><th>Variant</th><th>Field API calls</th><th>API ratio to S3</th><th>Complete / 8</th><th>Relations</th></tr></thead><tbody>'
    for r in summaries:
        vals = [r['variant'], f"{r['field_primitives']['fieldOperations']:,}", f"{r['field_api_ratio_to_s3']:.4f}" if r['field_api_ratio_to_s3'] is not None else 'null', r['complete_queries'], r['verified_relations']]
        panel += '<tr>' + ''.join('<td>' + html.escape(str(x)) + '</td>' for x in vals) + '</tr>'
    panel += '</tbody></table></div><div class="panel-head"><p>The three-size calibrated goal is unmet. Cold first-hit regressions, incomplete trials and the excluded runner incident remain in RESULTS.md. All certificates replay with zero failures.</p></div></div>\n<!-- nagao-reuse-01-end -->\n'
    (HERE / 'scoreboard.html').write_text(panel)
    path = ROOT / 'docs/index-calculus-scoreboard.html'
    page = path.read_text()
    start, end = '<!-- nagao-reuse-01-start -->', '<!-- nagao-reuse-01-end -->'
    if start in page:
        a, b = page.index(start), page.index(end) + len(end)
        page = page[:a] + panel.rstrip() + page[b:]
    else:
        page = page.replace('<!-- nagao-blocks-01-start -->', panel + '\n<!-- nagao-blocks-01-start -->')
    path.write_text(page)
    print(json.dumps({'cold': cold, 'batch': [{k: r[k] for k in ('variant', 'field_primitives', 'field_api_ratio_to_s3', 'seconds')} for r in summaries], 'bound': bound, 'audit': audit_report}, indent=2))


if __name__ == '__main__':
    main()
