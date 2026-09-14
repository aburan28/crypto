"""Audit immutable evidence, compare equal work and render literal tables."""
from collections import Counter
import html
import gzip
import json
from pathlib import Path
import experiment as exp

HERE, ROOT = exp.HERE, exp.ROOT


def raw_bytes(name='raw.jsonl'):
    plain, packed = HERE / name, HERE / (name + '.gz')
    if packed.exists():
        data = gzip.decompress(packed.read_bytes())
        if plain.exists() and plain.read_bytes() != data:
            raise ArithmeticError('compressed evidence differs from original')
        return data
    return plain.read_bytes()


def dump(name, data):
    (HERE / name).write_text(json.dumps(data, indent=2) + '\n')


def field(row):
    return row['counters']['field_api']['totals']


def words(row):
    return sum(sum(phase.values()) for phase in row['counters']['binary_word_operations'].values())


def resolved(row):
    return row['status'] != 'timeout' and row['within_budget']


def main():
    raw = [json.loads(s) for s in raw_bytes().splitlines()]
    contract = json.loads((HERE / 'contract.json').read_text())
    targets = json.loads((HERE / 'targets.json').read_text())
    trials = [r for r in raw if r['kind'] == 'trial']
    instances = {r['instance_sha256']: r for r in raw if r['kind'] == 'instance'}
    validation = next(r for r in raw if r['kind'] == 'validation')
    publication = json.loads((HERE / 'publication.json').read_text())
    if raw[0]['commit'] != publication['local_frozen_source_commit']:
        raise ArithmeticError('published source mapping mismatch')
    if len(trials) != contract['expected_trials'] or len(instances) != len(targets):
        raise ArithmeticError('incomplete campaign')
    for path, digest in raw[0]['sha256'].items():
        if exp.digest(ROOT / path) != digest:
            raise ArithmeticError('changed source/input: ' + path)
    prior = [json.loads(s) for s in (ROOT / contract['predecessor']).read_text().splitlines()]
    if exp.digest(ROOT / contract['predecessor']) != contract['predecessor_sha256']:
        raise ArithmeticError('predecessor evidence changed')
    expected_frozen = {(r['n'], tuple(r['target'])) for r in prior if r['kind'] == 'instance'
                       and r['base_kind'] == 'prefix' and r['d'] == 8}
    got_frozen = {(r['n'], tuple(r['target'])) for r in targets if r['cohort'].startswith('frozen-')}
    if expected_frozen != got_frozen:
        raise ArithmeticError('predecessor coverage mismatch')
    oracles = {}
    relations = rejections = 0
    for item in targets:
        digest = item['instance_sha256']
        content = {k: v for k, v in item.items() if k != 'instance_sha256'}
        if exp.hashlib.sha256(json.dumps(content, sort_keys=True).encode()).hexdigest() != digest:
            raise ArithmeticError('input hash mismatch')
        inst = instances[digest]
        if any(inst[k] != value for k, value in item.items()):
            raise ArithmeticError('raw instance differs from frozen input')
        key = item['n'], tuple(item['basis']), tuple(item['coefficients'])
        if key not in oracles:
            f, c, v = exp.old.context(item, False)
            oracles[key] = exp.old.Oracle(f, c, v)
        oracle = oracles[key]
        expected = oracle.expected(tuple(oracle.f.fromCoords(x) for x in item['target']))
        if expected != {tuple(x) for x in inst['expected']}:
            raise ArithmeticError('stored oracle result differs')
        rr = [r for r in trials if r['instance_sha256'] == digest]
        if Counter((r['variant'], r['mode']) for r in rr) != Counter(
                (variant, mode) for variant in contract['variants'] for mode in contract['modes']):
            raise ArithmeticError('missing or duplicate comparison cell')
        checker = exp.Checker(item)
        for row in rr:
            if any(row[k] != value for k, value in item.items()):
                raise ArithmeticError('trial inputs differ')
            if row['expected_count'] != exp.check(row, item, oracle, checker):
                raise ArithmeticError('wrong expected count')
            if len(row['solutions']) != len({tuple(x) for x in row['solutions']}):
                raise ArithmeticError('duplicate projected relation')
            if any(len(cert['points']) != 3 for cert in row['certificates']):
                raise ArithmeticError('wrong relation arity')
            if abs(sum(row['phase_seconds'].values()) - row['all_phase_seconds']) > 1e-6:
                raise ArithmeticError('phase timers overlap or omit work')
            counts = row['counters']['field_api']
            for name, value in counts['totals'].items():
                if value != sum(p[name] for p in counts['phases'].values()):
                    raise ArithmeticError('phase operation mismatch')
            if counts['totals']['fieldOperations'] != sum(counts['totals'][k] for k in ('additions', 'multiplications', 'squarings')):
                raise ArithmeticError('double-counted operation sum')
            relations += len(row['certificates'])
            rejections += len(row.get('rejection_certificates', []))
    pairs = []
    for item in targets:
        for mode in contract['modes']:
            rr = {r['variant']: r for r in trials if r['instance_sha256'] == item['instance_sha256'] and r['mode'] == mode}
            for reference in ('bilinear-reference', 'hybrid-filtered', 'pair-invariants-s3'):
                for candidate in ('bilinear-reference', 'bilinear-blocks'):
                    if reference == candidate:
                        continue
                    a, b = rr[reference], rr[candidate]
                    equal = resolved(a) and resolved(b)
                    if equal and mode == 'enumerate' and a['solutions'] != b['solutions']:
                        raise ArithmeticError('matched complete outputs differ')
                    pairs.append({'instance_sha256': item['instance_sha256'], 'n': item['n'], 'd': item['d'],
                                  'cohort': item['cohort'], 'stratum': item['stratum'], 'mode': mode,
                                  'reference': reference, 'candidate': candidate, 'equal_completed_work': equal,
                                  'reference_status': a['status'], 'candidate_status': b['status'],
                                  'reference_field_primitives': field(a), 'candidate_field_primitives': field(b),
                                  'reference_binary_words': words(a), 'candidate_binary_words': words(b),
                                  'field_api_ratio_candidate_over_reference': field(b)['fieldOperations'] / field(a)['fieldOperations'] if equal else None,
                                  'binary_word_ratio_candidate_over_reference': words(b) / words(a) if equal and words(a) else None,
                                  'common_operation_speedup': None})
    gate_pairs = [p for p in pairs if p['mode'] == 'enumerate' and p['reference'] == 'bilinear-reference'
                  and p['candidate'] == 'bilinear-blocks' and p['equal_completed_work']
                  and instances[p['instance_sha256']]['target'][0]]
    a = sum(p['reference_field_primitives']['fieldOperations'] for p in gate_pairs)
    b = sum(p['candidate_field_primitives']['fieldOperations'] for p in gate_pairs)
    aw = sum(p['reference_binary_words'] for p in gate_pairs)
    bw = sum(p['candidate_binary_words'] for p in gate_pairs)
    gate = {'complete_equal_enumeration_pairs': len(gate_pairs), 'reference_field_api_sum': a,
            'candidate_field_api_sum': b, 'candidate_over_reference_field_api': b / a if a else None,
            'reference_binary_words': aw, 'candidate_binary_words': bw,
            'passed': len(gate_pairs) >= 4 and 5 * b <= 4 * a and bw <= aw,
            'interpretation': 'Exploratory screening only; no calibrated total-cost claim.'}
    groups = []
    for n, d, cohort, stratum, variant, mode in sorted({(r['n'], r['d'], r['cohort'], r['stratum'], r['variant'], r['mode']) for r in trials}):
        rr = [r for r in trials if (r['n'], r['d'], r['cohort'], r['stratum'], r['variant'], r['mode']) == (n, d, cohort, stratum, variant, mode)]
        groups.append({'n': n, 'd': d, 'cohort': cohort, 'stratum': stratum, 'variant': variant, 'mode': mode,
                       'attempts': len(rr), 'resolved': sum(map(resolved, rr)),
                       'verified_relations': sum(len(r['solutions']) for r in rr),
                       'seconds': sum(r['all_phase_seconds'] for r in rr)})
    profiles = []
    for variant in ('bilinear-reference', 'bilinear-blocks'):
        for n in (18, 30):
            for mode in contract['modes']:
                rr = [r for r in trials if (r['variant'], r['n'], r['mode']) == (variant, n, mode)]
                stats, times, field_counts, binary_counts, sizes = Counter(), Counter(), Counter(), Counter(), Counter()
                for r in rr:
                    stats.update({k: value for k, value in r.get('stats', {}).items() if isinstance(value, int)})
                    times.update(r['phase_seconds'])
                    field_counts.update(field(r))
                    for phase in r['counters']['binary_word_operations'].values():
                        binary_counts.update(phase)
                    for k, value in r.get('stats', {}).get('rejections_by_free_bits', {}).items():
                        sizes[k] += value['blocks']
                profiles.append({'variant': variant, 'n': n, 'mode': mode, 'stats': dict(stats),
                                 'phase_seconds': dict(times), 'field_primitives': dict(field_counts),
                                 'binary_word_operations': dict(binary_counts), 'rejections_by_free_bits': dict(sizes),
                                 'scope': 'Work actually performed, including timeouts; unequal completed workloads.'})
    audit = {'matched_trials': len(trials), 'frozen_inputs': len(got_frozen), 'fresh_inputs': len(targets) - len(got_frozen),
             'relation_certificates_replayed': relations, 'rejection_certificates_replayed': rejections,
             'correctness_failures': 0, 'source_hashes_match': True,
             'raw_sha256': exp.hashlib.sha256(raw_bytes()).hexdigest(), 'source_commit': raw[0]['commit'],
             'validation': validation}
    summary = {'audit': audit, 'groups': groups, 'profiles': profiles, 'batch_gate': gate,
               'classification': 'engineering diagnostic; no advance established',
               'common_operation_speedup': None, 'full_dlp_S': None, 'rho_ratio': None, 'floor_ratio': None}
    dump('audit.json', audit)
    dump('summary.json', summary)
    dump('comparison.json', {'pairs': pairs, 'batch_gate': gate,
                            'limits': 'Ratios only for resolved matched tasks. API sums and word diagnostics are not a calibrated total cost.'})
    render(trials, targets, contract, summary, pairs)
    if (HERE / 'batch_raw.jsonl').exists() or (HERE / 'batch_raw.jsonl.gz').exists():
        import analyze_batch
        analyze_batch.run()
    print(json.dumps({'audit': audit, 'batch_gate': gate}, indent=2))


def table(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |', '|' + '|'.join(['---'] * len(headers)) + '|'] +
                     ['| ' + ' | '.join(map(str, row)) + ' |' for row in rows])


def render(trials, targets, contract, summary, pairs):
    count_rows, details, cost_rows = [], [], []
    for variant in contract['variants']:
        counts = [sum(resolved(r) for r in trials if r['variant'] == variant and r['mode'] == mode) for mode in contract['modes']]
        count_rows.append([variant, f'{counts[0]}/16', f'{counts[1]}/16', '0', 'unmeasured'])
    for n, d in sorted({(r['n'], r['d']) for r in targets}):
        for variant in contract['variants']:
            rr = [r for r in trials if (r['n'], r['d'], r['variant']) == (n, d, variant)]
            counts = [f'{sum(resolved(r) for r in rr if r["mode"] == m)}/{sum(r["mode"] == m for r in rr)}' for m in contract['modes']]
            details.append([n, d, variant, *counts])
    # Comparisons to the strongest Semaev baseline, only on jointly resolved tasks.
    for variant in contract['variants']:
        for mode in contract['modes']:
            selected = []
            for item in targets:
                rr = {r['variant']: r for r in trials if r['instance_sha256'] == item['instance_sha256'] and r['mode'] == mode}
                if resolved(rr[variant]) and resolved(rr['pair-invariants-s3']):
                    selected.append((rr[variant], rr['pair-invariants-s3']))
            counts = sum(field(a)['fieldOperations'] for a, _ in selected)
            ref_counts = sum(field(b)['fieldOperations'] for _, b in selected)
            cost_rows.append([variant, mode, len(selected), f'{counts:,}' if selected else '—',
                              f'{counts/ref_counts:.4f}' if ref_counts else '—',
                              'unmeasured', 'unmeasured'])
    gate = summary['batch_gate']
    audit = summary['audit']
    overview = table(['Variant', 'First resolved', 'Enumeration complete', 'Correctness failures', 'Common-cost speedup'], count_rows)
    breakdown = table(['Bits', 'd', 'Variant', 'First resolved', 'Enumeration complete'], details)
    costs = table(['Variant', 'Mode', 'Joint resolved tasks', 'Field API sum', 'API ratio to matched S3', 'Total-cost ratio', 'Ratio to floor'], cost_rows)
    lines = ['# Block rejection: measured solver-stage outcome', '',
             '[Derivation and rejection proof](README.md) · [contract](contract.json) · [raw](raw.jsonl.gz) · '
             '[summary](summary.json) · [paired comparison](comparison.json) · [audit](audit.json)', '',
             'The 192 three-second cold trials cover 16 identical inputs per variant: eight frozen '
             'and eight fresh, at 18 and 30 bits with dimensions 8 and 9. All six variants are rerun. '
             'First-mode resolutions include proved empty targets. Timeouts remain unknown and '
             'over-budget completions do not count as resolved.', '', overview, '',
             'Classification: engineering diagnostic; no advance or calibrated total-cost speedup '
             'established. One timing repetition cannot establish a runtime improvement. The '
             'separate field and binary counters do not include a measured conversion for SAT, '
             'allocation or control. The broad regression, full-DLP and three-size goal gates remain open.', '',
             '## Matched cost diagnostics', '',
             'Every numeric cost cell below is an equal-weight sum of field additions, multiplications '
             'and squarings, with expanded inversion internals. Ratios use the direct S3 costs on the '
             'same jointly resolved tasks. Each row can have a different subset, so compare its ratio '
             'within the row, not totals across rows. API ratios are component diagnostics. '
             'SAT has no resolved pair here unless shown. Complete primitive vectors, binary work '
             'and every regression are retained in comparison.json.', '', costs, '',
             '## Predeclared pruning screen', '',
             f'The batch gate **{"passes" if gate["passed"] else "fails"}**. There are '
             f'{gate["complete_equal_enumeration_pairs"]} complete equal-output nonzero-target enumeration pairs '
             'between the pruned and unpruned circuits. The gate requires at least four, 20% fewer '
             'field API calls and no increase in separately counted binary word work.', '',
             '```json', json.dumps(gate, indent=2), '```', '',
             'A failed screen ends this candidate before optional matched batches. Pruned-branch '
             'counts alone cannot reverse the verdict. Profiles include work actually performed '
             'on timed-out searches and must not be read as complete-work speedups.', '',
             '## Matched coverage by size', '', breakdown, '',
             'Cohort and uniform/supported strata remain separate in summary.json. No random-target '
             'yield estimate is inferred from supported targets.', '',
             '## Correctness and boundaries', '',
             f'The independent replay checked {audit["relation_certificates_replayed"]:,} relation certificates '
             f'and {audit["rejection_certificates_replayed"]:,} rejection certificates in the larger trials, '
             'with zero failures. Tiny validation covers all 53 finite GF64 targets on three d4 spaces: '
             '159 target/space pairs, 318 new-solver runs, 7,018 rejection certificates and '
             '7,168 exact circuit identity evaluations including fresh 18/30-bit samples. '
             'Both new solvers agree with exhaustive signed triples and an independent pair oracle.', '',
             'The signed-triple success ceiling from bounds_01 is unchanged. On this prefix support '
             'it is 100% at 18 bits (a vacuous ceiling), and 0.333461% / 2.233163% at 30 bits '
             'for d8 / d9. These are counting ceilings over uniform affine targets in the full curve '
             'group, not measured success rates or prime-subgroup claims. The old seven-multiplication '
             'branch floor does not apply after this circuit change. No ratio to an invented new '
             'full-DLP floor is reported.', '',
             f'Raw SHA256: `{audit["raw_sha256"]}`. All recorded source and input hashes match. '
             'Run analyze.py to reproduce the audit and these literal tables. Frozen solver source '
             'is mapped to the published Git tree in publication.json.', '']
    (HERE / 'RESULTS.md').write_text('\n'.join(lines))
    headers = ['Variant', 'First resolved / 16', 'Enumeration complete / 16', 'Correctness failures', 'Common-cost speedup']
    htable = '<div class="table-scroll"><table><caption>Three-second cold matched trials; first includes proved empty targets.</caption><thead><tr>' + ''.join('<th>' + html.escape(x) + '</th>' for x in headers) + '</tr></thead><tbody>'
    htable += ''.join('<tr>' + ''.join('<td>' + html.escape(str(x)) + '</td>' for x in row) + '</tr>' for row in count_rows)
    htable += '</tbody></table></div>'
    panel = '<!-- nagao-blocks-01-start -->\n<div class="panel" id="nagao-blocks-01"><div class="panel-head"><h2>Certified coefficient blocks: pruning must pay for its rank tests</h2><p>Frozen source: <code>research/nagao_relations/blocks_01/summary.json</code>. 192 matched cold trials, eight frozen inputs and eight fresh holdouts. The bilinear circuit and every saved rejection certificate pass independent replay. The direct S3 table remains a Semaev baseline. Classification: engineering diagnostic, no advance established.</p></div>' + htable + '<div class="panel-head"><p>Predeclared batch screen: ' + ('passes' if gate['passed'] else 'fails') + '. Complete equal-output enumeration pairs: ' + str(gate['complete_equal_enumeration_pairs']) + '. Field API counts and binary word operations remain separate; calibrated total cost, S, rho and floor ratios are unmeasured. The signed-triple yield ceiling is unchanged. One timing repetition and two field sizes do not meet the original goal. Earlier panels remain historical evidence.</p></div></div>\n<!-- nagao-blocks-01-end -->\n'
    (HERE / 'scoreboard.html').write_text(panel)
    path = ROOT / 'docs/index-calculus-scoreboard.html'
    page = path.read_text()
    start, end = '<!-- nagao-blocks-01-start -->', '<!-- nagao-blocks-01-end -->'
    if start in page:
        a, b = page.index(start), page.index(end) + len(end)
        page = page[:a] + panel.rstrip() + page[b:]
    else:
        anchor = '  <div class="panel" id="nagao-bound-audit">'
        if page.count(anchor) != 1:
            raise ArithmeticError('scoreboard insertion anchor missing')
        page = page.replace(anchor, panel + '\n' + anchor)
    path.write_text(page)


if __name__ == '__main__':
    main()
