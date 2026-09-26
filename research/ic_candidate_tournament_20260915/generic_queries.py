"""Replay bounded generic query records; this alone is not scientific admission.

The query RNG, exclusive timing, source binding and matrix rank audit remain
separate requirements. No performance totals or promotion decision are produced.
"""
from collections import Counter
import hashlib
import json
from oracle import Curve, require

OUTCOMES = {'witness', 'proved_unsat', 'incomplete', 'unsupported',
            'invalid_model', 'unresolved', 'identity'}


def integer(value, name):
    require(type(value) is int and value >= 0, f'invalid {name}')
    return value


def verify_queries(report, fixture, summands):
    require(report.get('query_schema_version') == 1, 'missing query schema')
    require(report.get('mode') == 'ic' and report.get('fixture') == fixture,
            'query fixture mismatch')
    require(2 <= summands <= 4, 'unsupported summand count')
    curve = Curve(fixture)
    base = [curve.decode(p) for p in report['factor_base']]
    require(base and all(p is not None for p in base), 'invalid query base')
    # Exact negative controls are bounded to the two-summand toy integration
    # inputs. Larger negative claims need a separate proof adapter.
    pair_sums = None
    counts = Counter()

    def check(attempt, trial, q, collecting):
        nonlocal pair_sums
        require(integer(attempt['trial'], 'trial') == trial, 'query chronology gap')
        a, b = integer(attempt['a'], 'a'), integer(attempt['b'], 'b')
        require(a < curve.r and b < curve.r, 'query coefficient outside subgroup')
        require(b == 0 if collecting else b > 0, 'query coefficient role mismatch')
        target = curve.add(curve.mul(curve.g, a), curve.mul(q, b))
        pdp = attempt['pdp']
        outcome, points, stats = pdp['outcome'], pdp['points'], pdp['stats']
        require(outcome in OUTCOMES and outcome != 'invalid_model', 'invalid PDP verdict')
        family = stats['family']
        require(family in {'none', 'groebner', 'sat', 'crossbred', 'weil_chart'},
                'unknown counter family')
        if family != 'none':
            native = stats['stats']
            if family == 'weil_chart':
                native = native['solver']
            require(type(native['unsupported']) is bool and type(native['exhausted']) is bool,
                    'invalid completion flags')
            for key, value in native.items():
                if key in {'unsupported', 'exhausted', 'refuted'}:
                    require(type(value) is bool, 'invalid boolean counter')
                else:
                    integer(value, key)
            require(not native.get('spurious', 0), 'invalid model retained')
            if outcome == 'unsupported':
                require(native['unsupported'] and native['exhausted'], 'unsupported flags lost')
            if outcome == 'incomplete':
                require(native['exhausted'] and not native['unsupported'], 'incomplete flags lost')
            if outcome == 'proved_unsat':
                require(not native['unsupported'] and not native['exhausted'], 'incomplete refutation')
                if family == 'sat':
                    require(native['refuted'], 'missing SAT refutation')
        if outcome == 'identity':
            require(target is None and points is None and family == 'none', 'false identity query')
        else:
            require(target is not None, 'identity query sent to PDP')
        if outcome == 'witness':
            require(isinstance(points, list) and len(points) == summands, 'wrong witness size')
            total = None
            for index in points:
                require(type(index) is int and 0 <= index < len(base), 'invalid witness index')
                total = curve.add(total, base[index])
            require(total == target, 'false query witness')
        else:
            require(points is None, 'nonwitness carries points')
        if outcome == 'proved_unsat':
            require(summands == 2 and curve.n <= 13 and len(base) <= 256,
                    'negative claim needs an independent proof adapter')
            if pair_sums is None:
                pair_sums = {curve.add(p, q) for p in base for q in base}
            require(target not in pair_sums, 'false negative: decomposition exists')
        counts[outcome] += 1
        return points

    trial = 0
    relations = []
    for batch in report['collection_reports']:
        attempts = batch['attempts']
        require(len(attempts) == integer(batch['trials'], 'batch trials'), 'dropped collection query')
        successes = 0
        for attempt in attempts:
            points = check(attempt, trial, None, True)
            if points is not None:
                successes += 1
                relations.append(dict(trial=trial, a=attempt['a'], points=points))
            trial += 1
        require(successes == integer(batch['relations'], 'batch relations'), 'batch yield mismatch')
    require(trial == integer(report['trials'], 'collection trials'), 'collection total mismatch')
    require(relations == report['relations'], 'relations do not match query witnesses')
    collection_counts = dict(sorted(counts.items()))
    descent_trials = 0
    seen = set()
    for solution in report.get('solutions', []):
        index = integer(solution['index'], 'solution index')
        require(index < len(fixture['targets']) and index not in seen, 'invalid target index')
        seen.add(index)
        q = curve.decode(fixture['targets'][index])
        attempts = solution['attempts']
        require(len(attempts) == integer(solution['trials'], 'descent trials'), 'dropped descent query')
        for t, attempt in enumerate(attempts):
            check(attempt, t, q, False)
        relation = solution['relation']
        if solution['recovered'] is not None:
            require(attempts and relation is not None, 'missing final descent query')
            last = attempts[-1]
            require(relation == dict(a=last['a'], b=last['b'], points=last['pdp']['points'])
                    and last['pdp']['outcome'] == 'witness', 'final witness mismatch')
            scalar = int(solution['recovered'])
            require(0 <= scalar < curve.r and curve.mul(curve.g, scalar) == q, 'false recovered scalar')
        else:
            require(relation is None, 'failed descent has a success certificate')
        descent_trials += len(attempts)
    payload = dict(collection=[batch['attempts'] for batch in report['collection_reports']],
                   descent=[solution['attempts'] for solution in report.get('solutions', [])])
    query_digest = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':'),
                                            allow_nan=False).encode()).hexdigest()
    return dict(schema_version=1, query_sha256=query_digest,
                collection_queries=trial, descent_queries=descent_trials,
                collection_outcomes=collection_counts, all_outcomes=dict(sorted(counts.items())),
                scope='query accounting and group replay only', promotion_eligible=False)
