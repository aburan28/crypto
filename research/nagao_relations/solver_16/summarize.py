"""Aggregate raw.jsonl into summary.json and the tables of RESULTS.md.
Cites; never re-solves.  Fits an exponent only where the contract allows
(>= 4 complete sizes, every uniform cell decided with a conflict count)."""
import json
import math
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load():
    """raw.jsonl holds the CryptoMiniSat arm; raw_wdsat.jsonl the WDSat arms,
    which ran in a second process (see run.py, SOLVER16_RAW).  Both frozen."""
    rows = []
    for name in ('raw.jsonl', 'raw_wdsat.jsonl'):
        path = HERE / name
        if path.exists():
            for line in path.read_text().splitlines():
                if line.strip():
                    rows.append(json.loads(line))
    return rows


def fit(xs, ys):
    """Least squares y = a x + b; returns slope, intercept, residual SE."""
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    a = sxy / sxx
    b = my - a * mx
    res = [y - (a * x + b) for x, y in zip(xs, ys)]
    se = math.sqrt(sum(r * r for r in res) / max(1, n - 2))
    return a, b, se


def main():
    rows = load()
    contract = json.loads((HERE / 'contract.json').read_text())
    targets = json.loads((HERE / 'targets.json').read_text())['targets']
    out = {'cells': len(rows), 'errors': [r['cell'] for r in rows if r['status'] == 'ERROR'], 'tables': {}, 'fits': {}, 'matched_rr_over_s4': [],
           'wdsat_plain_vs_full_tree': [], 'xg_audit': {}, 'correctness': {}}

    # correctness gates
    bad_uniform = [r['cell'] for r in rows if r['stratum'] == 'uniform' and r.get('genuine_relations', 0) > 0]
    planted_missed = [r['cell'] for r in rows if r['stratum'] == 'planted' and r['solver'] in ('cms', 'wdsat_plain')
                      and r['status'] not in ('TIMEOUT', 'FALSE_SAT_PARTIAL_ASSIGNMENT') and not r.get('planted_found')]
    unverified = [r['cell'] for r in rows if any(not v['relation'] for v in r.get('verification', []))]
    out['correctness'] = {'uniform_cells_with_a_relation': bad_uniform, 'planted_cells_decided_without_the_planted_triple': planted_missed,
                          'cells_with_an_unverified_witness': unverified,
                          'total_verified_relations': sum(r.get('genuine_relations', 0) for r in rows if r['solver'] != 'wdsat_xg_audit')}

    # CMS tables
    for enc in ('s4', 'rr'):
        table = []
        for d in sorted({r['d'] for r in rows if r['solver'] == 'cms' and r['encoding'] == enc}):
            cells = [r for r in rows if r['solver'] == 'cms' and r['encoding'] == enc and r['d'] == d]
            uni = [r for r in cells if r['stratum'] == 'uniform']
            pla = [r for r in cells if r['stratum'] == 'planted']
            uconf = [r['conflicts'] for r in uni if r['conflicts'] is not None]
            entry = {'d': d, 'admissible': targets[str(d)]['admissible'], 'pairs': targets[str(d)]['pairs'],
                     'uniform_decided': sum(1 for r in uni if r['status'] == 'UNSAT'), 'uniform_attempted': len(uni),
                     'uniform_timeouts': sum(1 for r in uni if r['status'] == 'TIMEOUT'),
                     'uniform_conflicts': uconf, 'uniform_conflicts_mean': statistics.mean(uconf) if uconf and len(uconf) == len(uni) else None,
                     'uniform_wall_mean': statistics.mean(r['wall'] for r in uni) if uni else None,
                     'planted_first_found': sum(1 for r in pla if r.get('planted_found')), 'planted_attempted': len(pla),
                     'planted_enumeration_complete': sum(1 for r in pla if r['status'] == 'UNSAT' and r.get('planted_found')),
                     'planted_conflicts': [r['conflicts'] for r in pla], 'planted_wall_mean': statistics.mean(r['wall'] for r in pla) if pla else None,
                     'planted_genuine_relations': sum(r.get('genuine_relations', 0) for r in pla)}
            entry['conflicts_over_pairs'] = (entry['uniform_conflicts_mean'] / entry['pairs']) if entry['uniform_conflicts_mean'] else None
            table.append(entry)
        out['tables']['cms_' + enc] = table
        complete = [e for e in table if e['uniform_conflicts_mean'] and e['uniform_decided'] == e['uniform_attempted'] and e['uniform_attempted'] > 0]
        if len(complete) >= 4:
            xs = [e['d'] for e in complete]
            ys = [math.log2(e['uniform_conflicts_mean']) for e in complete]
            a, b, se = fit(xs, ys)
            # the null object's slope over the same d
            pa, pb, pse = fit(xs, [math.log2(e['pairs']) for e in complete])
            out['fits']['cms_' + enc] = {'d_values': xs, 'log2_mean_conflicts': ys, 'slope_per_d': a, 'intercept': b, 'residual_se': se,
                                         'pairs_slope_per_d': pa, 'passes_contract_slope_below_1.8_and_se_below_0.5': (a < 1.8 and se < 0.5),
                                         'note': 'slope in d; |F| ~ 2^(d-1), so an exponent e in |F| is a slope e per d'}
            # also the fit excluding the linear-certificate regime (d <= 6), where it exists
            tail = [e for e in complete if e['d'] >= 7]
            if len(tail) >= 3:
                a2, b2, se2 = fit([e['d'] for e in tail], [math.log2(e['uniform_conflicts_mean']) for e in tail])
                out['fits']['cms_' + enc + '_d_ge_7'] = {'d_values': [e['d'] for e in tail], 'slope_per_d': a2, 'residual_se': se2,
                                                        'sizes': len(tail), 'admissible_as_fit': len(tail) >= 4}
        else:
            out['fits']['cms_' + enc] = {'admissible': False, 'complete_sizes': [e['d'] for e in complete],
                                         'reason': 'fewer than four sizes with every uniform cell decided and counted'}

    # matched rr / s4 under CMS
    by = {(r['encoding'], r['d'], r['stratum'], r['target_index']): r for r in rows if r['solver'] == 'cms'}
    for (enc, d, st, i), r in sorted(by.items()):
        if enc != 'rr':
            continue
        s = by.get(('s4', d, st, i))
        if not s:
            continue
        ratio = None
        kind = None
        if r['conflicts'] is not None and s['conflicts'] is not None:
            ratio, kind = r['conflicts'] / max(1, s['conflicts']), 'conflicts'
        elif r['status'] != 'TIMEOUT' and s['status'] != 'TIMEOUT':
            ratio, kind = r['wall'] / max(1e-3, s['wall']), 'wall (both complete, conflicts unavailable)'
        out['matched_rr_over_s4'].append({'d': d, 'stratum': st, 'index': i, 'rr_status': r['status'], 's4_status': s['status'],
                                          'rr_conflicts': r['conflicts'], 's4_conflicts': s['conflicts'], 'rr_wall': r['wall'], 's4_wall': s['wall'],
                                          'ratio': ratio, 'ratio_kind': kind,
                                          'rr_censored_while_s4_complete': r['status'] == 'TIMEOUT' and s['status'] != 'TIMEOUT'})

    # WDSat plain vs the full tree.  In FIND_ALL + symmetry mode the solver
    # prints assignments from a table the plain search does not maintain, so
    # an 'enum' cell's witnesses are unreadable by construction; its decision
    # count is the exhaustion cost, and the witness gate is the paired 'first'
    # cell on the same instance.
    firstCells = {(r['d'], r['stratum'], r['target_index']): r for r in rows if r['solver'] == 'wdsat_plain' and r['mode'] == 'first'}
    for r in rows:
        if r['solver'] == 'wdsat_plain':
            full = int(r['full_tree_over_sym'])
            status = r['status']
            plantedFound = r.get('planted_found')
            if r['mode'] == 'enum' and status == 'FALSE_SAT_PARTIAL_ASSIGNMENT':
                status = 'ENUMERATED (witness unreadable in FIND_ALL symmetry mode; see paired first cell)'
                pf = firstCells.get((r['d'], r['stratum'], r['target_index']))
                plantedFound = pf.get('planted_found') if pf else None
            out['wdsat_plain_vs_full_tree'].append({'cell': r['cell'], 'd': r['d'], 'stratum': r['stratum'], 'mode': r['mode'], 'status': status,
                                                    'decisions': r['decisions'], 'full_tree_2^3d_over_6': full,
                                                    'decisions_over_full_tree': (r['decisions'] / full) if r['decisions'] else None,
                                                    'pairs': r['pairs'], 'wall': r['wall'], 'planted_found': plantedFound})
    # -x audit
    audit = [r for r in rows if r['solver'] == 'wdsat_xg_audit']
    out['xg_audit'] = {'cells': len(audit), 'complete': sum(1 for r in audit if r['completeness'] == 'ok'),
                       'failed': [{'cell': r['cell'], 'status': r['status'], 'decisions_not_a_cost': r['decisions_not_a_cost']} for r in audit if r['completeness'] != 'ok'],
                       'by_args': {}}
    for r in audit:
        k = ' '.join(r['args'])
        out['xg_audit']['by_args'].setdefault(k, {'ok': 0, 'failed': 0})
        out['xg_audit']['by_args'][k]['ok' if r['completeness'] == 'ok' else 'failed'] += 1
    (HERE / 'summary.json').write_text(json.dumps(out, indent=1, sort_keys=True) + '\n')
    print(json.dumps({'cells': out['cells'], 'errors': out['errors'], 'correctness': out['correctness'], 'fits': out['fits'],
                      'xg_audit': {k: v for k, v in out['xg_audit'].items() if k != 'failed'}}, indent=1))
    for enc in ('s4', 'rr'):
        print('\ncms', enc)
        for e in out['tables']['cms_' + enc]:
            print(' d=%2d pairs=%9d uniform %d/%d decided (%d timeouts) mean conflicts %s  ratio/pairs %s | planted found %d/%d enum complete %d, genuine %d, mean wall %.2fs' % (
                e['d'], e['pairs'], e['uniform_decided'], e['uniform_attempted'], e['uniform_timeouts'],
                ('%.0f' % e['uniform_conflicts_mean']) if e['uniform_conflicts_mean'] else 'n/a',
                ('%.2f' % e['conflicts_over_pairs']) if e['conflicts_over_pairs'] else 'n/a',
                e['planted_first_found'], e['planted_attempted'], e['planted_enumeration_complete'], e['planted_genuine_relations'], e['planted_wall_mean'] or 0))
    for m in out['matched_rr_over_s4']:
        print(' rr/s4 d=%d %s%d: %s vs %s, ratio %s (%s)' % (m['d'], m['stratum'], m['index'], m['rr_status'], m['s4_status'],
                                                           ('%.1f' % m['ratio']) if m['ratio'] else 'censored', m['ratio_kind']))
    for w in out['wdsat_plain_vs_full_tree']:
        print(' wdsat %s: %s decisions=%s full=%d ratio=%s' % (w['cell'], w['status'], w['decisions'], w['full_tree_2^3d_over_6'],
                                                             ('%.3f' % w['decisions_over_full_tree']) if w['decisions_over_full_tree'] else 'n/a'))


if __name__ == '__main__':
    main()
