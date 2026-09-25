"""Development-only diversity retention and measured-combination proposals.

This module never sees confirmation/replay measurements and cannot promote an
algorithm. Complete implementations, not sums of the best stage times, compete.
"""
import copy
import itertools
import math

from oracle import require


AXES = ('factor_base', 'solver', 'linear_algebra', 'batch_trials', 'collection_window')


def family(arm):
    cfg = arm['config']
    base = cfg.get('factor_base', {})
    return (base.get('kind', 'subgroup_orbits'), cfg.get('solver', 'pair_table'),
            cfg.get('linear_algebra', 'sparse'))


def retain(comparisons, arms, *, width=6, exploration=1, seed=0):
    """Keep total-cost leaders, cell specialists, and a reproducible outsider.

    Only complete eligible development rows enter. Unselected candidates stay
    in the archive. Width limits evaluation, never the evidence retained.
    """
    import random
    require(width >= 2 and 0 <= exploration < width, 'invalid portfolio budget')
    by_id = {a['id']: a for a in arms}
    eligible = []
    for row in comparisons:
        if not row.get('eligible') or row['candidate'] not in by_id:
            continue
        vector = [row['candidate_over_baseline'], row['native_wall_candidate_over_baseline']]
        vector += list(row['per_cell'].values()) + list(row['native_wall_per_cell'].values())
        require(vector and all(math.isfinite(v) and v > 0 for v in vector), 'invalid portfolio cost')
        eligible.append(row)
    eligible.sort(key=lambda r: (r['candidate_over_baseline'], r['candidate']))
    if not eligible:
        return []
    keys = sorted(eligible[0]['per_cell'])
    require(all(sorted(r['per_cell']) == keys and sorted(r['native_wall_per_cell']) == keys
                for r in eligible), 'unmatched portfolio cells')
    def vector(r):
        return [r['candidate_over_baseline'], r['native_wall_candidate_over_baseline']] + [
            r[k][cell] for k in ('per_cell', 'native_wall_per_cell') for cell in keys]
    def dominates(a, b):
        av, bv = vector(a), vector(b)
        return all(x <= y for x, y in zip(av, bv)) and any(x < y for x, y in zip(av, bv))
    frontier = [r for r in eligible if not any(dominates(s, r) for s in eligible)]
    selected = []
    limit = min(width - exploration, len(eligible))
    def add(row, reason):
        if len(selected) < limit and row['candidate'] not in {s['candidate'] for s in selected}:
            selected.append({'candidate': row['candidate'], 'reason': reason})
    add(eligible[0], 'complete instruction-cost leader')
    add(min(eligible, key=lambda r: (r['native_wall_candidate_over_baseline'], r['candidate'])),
        'complete native-time leader')
    # Reserve diversity before filling with close variants of the leader.
    for row in frontier:
        known = {family(by_id[s['candidate']]) for s in selected}
        if family(by_id[row['candidate']]) not in known:
            add(row, 'non-dominated implementation family')
    for cell in keys:
        add(min(eligible, key=lambda r: (r['per_cell'][cell], r['candidate'])),
            'cell specialist: ' + cell)
    for row in frontier + eligible:
        add(row, 'Pareto frontier' if row in frontier else 'development reserve')
    remaining = [r for r in eligible if r['candidate'] not in {s['candidate'] for s in selected}]
    random.Random(seed).shuffle(remaining)
    selected.extend({'candidate': r['candidate'], 'reason': 'predeclared exploration slot'}
                    for r in remaining[:min(exploration, width-len(selected))])
    return selected


def factorial_candidates(base, *, panel='implementation', max_candidates=16):
    """Predeclare interactions before measuring; every tuple is a real job."""
    require(2 <= max_candidates <= 16, 'candidate budget must be 2..16')
    out = [{'id': 'incumbent', 'parent': None, 'hypothesis': 'Unmodified declared baseline.',
            'config': copy.deepcopy(base)}]
    if panel == 'implementation':
        axes = {'batch_trials': [1, 8, 64], 'linear_algebra': ['sparse', 'dense'],
                'collection_window': [None, 8]}
    elif panel == 'algebra':
        axes = {'solver': ['pair_table', 'enumerate', 'f4', 'f5', 'inherited_f4', 'sat_xor', 'sat_cnf'],
                'linear_algebra': ['sparse', 'dense']}
    elif panel == 'factor-base':
        axes = {'factor_base': [
            {'kind': 'subgroup_orbits', 'seed': 43, 'points': 52},
            {'kind': 'subgroup_orbits', 'seed': 71, 'points': 78},
            {'kind': 'subgroup_orbits', 'seed': 97, 'points': 104},
            {'kind': 'factor', 'index': 0},
            {'kind': 'frobenius_union', 'seed_masks': [1, 2, 4, 8]},
            {'kind': 'frobenius_union', 'seed_masks': [1, 2, 8]},
        ], 'linear_algebra': ['sparse', 'dense']}
    else:
        raise ValueError('unknown panel')
    for values in itertools.product(*axes.values()):
        cfg = dict(copy.deepcopy(base), **dict(zip(axes, values)))
        if cfg.get('collection_window') is None:
            cfg.pop('collection_window', None)
        if panel == 'factor-base':
            cfg.pop('factor_base_orbits', None)
            cfg.pop('factor_base_cube_root', None)
        if any(cfg == a['config'] for a in out):
            continue
        changed = {k: v for k, v in cfg.items() if base.get(k) != v}
        out.append({'id': f'{panel.replace("-", "_")}_{len(out):02d}', 'parent': 'incumbent',
                    'hypothesis': 'Measure complete-pipeline interaction: ' + ', '.join(changed),
                    'changed_parameters': changed, 'config': cfg,
                    'falsification': 'Incorrect certificate, incomplete workload or no total-cost improvement.'})
        if len(out) == max_candidates:
            break
    return out


def recombine(arms, retained_ids, *, max_candidates=16):
    """Pair disjoint configuration deltas, including individually slower parents.

    Source edits cannot be merged as JSON: those combinations need an explicitly
    built source candidate. Conflicting parameter edits are left separate.
    """
    require(2 <= max_candidates <= 16, 'candidate budget must be 2..16')
    incumbent = next(a for a in arms if a['id'] == 'incumbent')
    baseline = incumbent['config']
    parents = [a for a in arms if a['id'] in retained_ids and a['id'] != 'incumbent']
    require(len(parents) == len(set(retained_ids)-{'incumbent'}), 'unknown retained parent')
    out = [copy.deepcopy(incumbent)]
    for a, b in itertools.combinations(parents, 2):
        # Require the same executable/source; never discard a source mutation.
        if any(a.get(k) != incumbent.get(k) or b.get(k) != incumbent.get(k)
               for k in ('source_root', 'source_manifest_sha256', 'binary_relative')):
            continue
        if not set(baseline) <= set(a['config']) or not set(baseline) <= set(b['config']):
            continue
        da = {k:v for k,v in a['config'].items() if baseline.get(k) != v}
        db = {k:v for k,v in b['config'].items() if baseline.get(k) != v}
        if not da or not db or set(da) & set(db):
            continue
        cfg = dict(copy.deepcopy(baseline), **da, **db)
        if any(cfg == x['config'] for x in arms+out):
            continue
        out.append({'id': f'combination_{len(out):02d}', 'parents': [a['id'], b['id']],
                    'config': cfg, 'hypothesis': 'Test interaction of retained mechanisms in a complete run.',
                    'falsification': 'No gain on fresh development fixtures; confirmation remains unused.'})
        if len(out) == max_candidates:
            break
    return out
