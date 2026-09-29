#!/usr/bin/env python3
"""Tests for the round-0025 checker amendment (round25-oracle-mixed-summands.patch).

    python3 campaign_20260916/test_round25_oracle.py          # plain python
    python3 -m pytest -q campaign_20260916/test_round25_oracle.py

`OLD` is round 0023's frozen `evaluator/oracle.py`, the checker rounds 0023 and
0024 ran with. `NEW` is `evaluator-r25/oracle.py`, the same file plus the
patch. What is checked:

1. Hand-built certificates, from the checker's own group arithmetic on
   round 0024's `n13a0` fixture (r = 2003, so every logarithm is found by
   table lookup): a valid 3-term report under a 4-summand configuration
   verifies; a valid 4-term report verifies; a report declaring more
   summands than configured is refused, and so is one declaring fewer than
   three; mixing lengths within one report is refused; a tampered relation,
   descent or logarithm is refused at either length; and whenever the
   declared count equals the configured one the two checkers return the same
   certificate or the same refusal.
2. Real worker reports, when round 0024's frozen workers are present: the
   pair worker's 3-term report and the triple worker's 4-term report verify
   under `summands=4`, and tampering with either is refused.
3. Stored receipts, when the restored rounds are present: every VERIFIED
   report of rounds 0014-0020 that the old checker accepts, the new one accepts
   with the identical certificate, which is also the certificate the receipt
   recorded. By default one trial per (round, stage, cell, arm), profile
   report; `R25_RECEIPTS=rep0` checks repetition 0 of every trial;
   `R25_FULL=1` checks every receipt, profile and native report both.
"""
import copy
import importlib.util
import json
import os
import random
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


OLD = _load('oracle_r23', ROOT / 'runs/round-0023/evaluator/oracle.py')
NEW = _load('oracle_r25', HERE / 'evaluator-r25/oracle.py')
FIXTURES = json.loads((ROOT / 'runs/round-0024/fixtures.json').read_text())


def outcome(module, report, fixture, **kw):
    """A certificate, or the refusal's message: what a checker decided."""
    try:
        return ('ok', module.verify(report, fixture, **kw))
    except module.InvalidEvidence as exc:
        return ('refused', str(exc))


def refused(report, fixture, summands=4, module=NEW):
    kind, value = outcome(module, report, fixture, summands=summands)
    assert kind == 'refused', value
    return value


# ── Hand-built certificates ────────────────────────────────────────────────
class Toy:
    """A complete, honest certificate on a fixture small enough to take every
    logarithm by lookup. Built with the checker's own arithmetic, so what it
    exercises is the checker's acceptance logic, not a solver."""

    def __init__(self, fixture, orbits=2, seed=1):
        self.fixture = fixture
        self.c = c = NEW.Curve(fixture)
        self.log = {}
        p = None
        for d in range(c.r):
            self.log[p] = d
            p = c.add(p, c.g)
        self.target = c.decode(fixture['targets'][0])
        self.d = self.log[self.target]
        rng = random.Random(seed)
        reps, used = [], set()
        while len(reps) < orbits:
            rep = c.mul(c.g, rng.randrange(1, c.r))
            orbit, q = [], rep
            for _ in range(c.n):
                orbit += [q, c.neg(q)]
                q = c.frob(q)
            if used.isdisjoint(orbit) and len(set(orbit)) == len(orbit):
                used.update(orbit)
                reps.append(rep)
        self.reps = reps
        self.base = []
        for rep in reps:
            q = rep
            for _ in range(c.n):
                self.base += [q, c.neg(q)]
                q = c.frob(q)
        self.rng = rng

    def total(self, ids):
        q = None
        for i in ids:
            q = self.c.add(q, self.base[i])
        return q

    def relation(self, k):
        while True:
            ids = [self.rng.randrange(len(self.base)) for _ in range(k)]
            q = self.total(ids)
            if q is not None:
                return {'a': self.log[q], 'points': ids}

    def descent(self, k):
        while True:
            ids = [self.rng.randrange(len(self.base)) for _ in range(k)]
            q = self.total(ids)
            if q is None:
                continue
            b = self.rng.randrange(1, self.c.r)
            return {'a': (self.log[q] - b * self.d) % self.c.r, 'b': b, 'points': ids}

    def report(self, k, relations=6, lengths=None, descent_length=None):
        lengths = lengths or [k] * relations
        rels = [self.relation(m) for m in lengths]
        keys = {(r['a'], tuple(sorted(r['points']))) for r in rels}
        return {
            'schema_version': 1, 'mode': 'ic', 'status': 'complete',
            'fixture': self.fixture,
            'solutions': [{'index': 0, 'recovered': str(self.d),
                           'relation': self.descent(descent_length or k)}],
            'factor_base_orbits': [[str(x), str(y)] for x, y in self.reps],
            'column_convention': 'representative',
            'relations': rels,
            'column_logs': [{'point': [str(x), str(y)], 'log': str(self.log[(x, y)])} for x, y in self.reps],
            'columns': len(self.reps),
            'accepted_relations': len(keys), 'duplicate_relations': len(rels) - len(keys),
            'rejected_relations': 0, 'summands': k,
        }


def toy():
    fixture = next(c['fixture'] for c in FIXTURES['development'] if c['cell'] == 'n13a0')
    return Toy(fixture)


def test_three_term_report_verifies_under_a_four_summand_arm():
    t = toy()
    report = t.report(3)
    kind, proof = outcome(NEW, report, t.fixture, summands=4)
    assert kind == 'ok', proof
    assert proof['relation_summands'] == 3
    assert proof['solutions'] == [str(t.d)]
    # The unamended checker refuses exactly this, and nothing else about it.
    assert outcome(OLD, report, t.fixture, summands=4) == ('refused', 'changed summand count')
    # The same report at its own configured count: both checkers, one certificate,
    # which is the amended certificate without the new key.
    same = outcome(OLD, report, t.fixture, summands=3)
    assert same == outcome(NEW, report, t.fixture, summands=3)
    assert same[0] == 'ok' and 'relation_summands' not in same[1]
    assert {k: v for k, v in proof.items() if k != 'relation_summands'} == same[1]


def test_four_term_report_verifies_identically():
    t = toy()
    report = t.report(4)
    old = outcome(OLD, report, t.fixture, summands=4)
    assert old[0] == 'ok', old
    assert outcome(NEW, report, t.fixture, summands=4) == old
    assert 'relation_summands' not in old[1]


def test_more_summands_than_configured_is_refused():
    t = toy()
    five = t.report(5)
    assert outcome(NEW, five, t.fixture, summands=5)[0] == 'ok'  # an honest 5-term report
    assert refused(five, t.fixture, summands=4) == 'changed summand count'
    four = t.report(4)
    assert refused(four, t.fixture, summands=3) == 'changed summand count'


def test_fewer_than_three_summands_is_refused():
    t = toy()
    two = t.report(2)
    assert outcome(NEW, two, t.fixture, summands=2)[0] == 'ok'  # configured at two, as before
    assert refused(two, t.fixture, summands=4) == 'changed summand count'
    assert refused(two, t.fixture, summands=3) == 'changed summand count'


def test_declared_count_must_be_an_integer_below_the_configured_one():
    t = toy()
    report = t.report(3)
    for bad in ('3', 3.0, True, None, -3):
        r = copy.deepcopy(report)
        r['summands'] = bad
        assert refused(r, t.fixture, summands=4) == 'changed summand count', bad
        # At the configured count the old rule is kept to the letter, `3.0 == 3` included.
        assert outcome(OLD, r, t.fixture, summands=3) == outcome(NEW, r, t.fixture, summands=3), bad


def test_lengths_do_not_mix_within_a_report():
    t = toy()
    mixed = t.report(4, lengths=[4, 4, 3, 4, 4, 4])
    assert refused(mixed, t.fixture) == 'bad relation indices'
    mixed = t.report(3, lengths=[3, 3, 4, 3, 3, 3])
    assert refused(mixed, t.fixture) == 'bad relation indices'
    mixed = t.report(3, descent_length=4)
    assert refused(mixed, t.fixture) == 'bad descent relation indices'
    mixed = t.report(4, descent_length=3)
    assert refused(mixed, t.fixture) == 'bad descent relation indices'


def tamperings(report):
    """(label, tampered copy, expected refusal)."""
    r = report['relations'][0]
    out = []
    x = copy.deepcopy(report)
    x['relations'][0]['a'] = r['a'] % (int(report['fixture']['subgroup_order']) - 1) + 1
    out.append(('relation scalar', x, 'incorrect point relation'))
    x = copy.deepcopy(report)
    x['relations'][0]['points'][0] = (r['points'][0] + 2) % (len(report['factor_base_orbits']) * 2 * report['fixture']['degree'])
    out.append(('relation index', x, 'incorrect point relation'))
    x = copy.deepcopy(report)
    x['relations'][0]['points'] = r['points'][:-1]
    out.append(('relation shortened', x, 'bad relation indices'))
    x = copy.deepcopy(report)
    x['relations'][0]['points'] = r['points'] + [r['points'][0]]
    out.append(('relation lengthened', x, 'bad relation indices'))
    x = copy.deepcopy(report)
    rel = x['solutions'][0]['relation']
    rel['a'] = (rel['a'] + 1) % int(report['fixture']['subgroup_order'])
    out.append(('descent scalar', x, 'descent relation does not hold in the group'))
    x = copy.deepcopy(report)
    x['column_logs'][0]['log'] = str((int(x['column_logs'][0]['log']) + 1) % int(report['fixture']['subgroup_order']))
    out.append(('column log', x, 'bad column log'))
    x = copy.deepcopy(report)
    order = int(report['fixture']['subgroup_order'])
    x['solutions'][0]['recovered'] = str((int(x['solutions'][0]['recovered']) + 1) % order)
    out.append(('logarithm', x, 'incorrect scalar'))
    return out


def test_tampered_reports_are_refused_at_either_length():
    t = toy()
    for k in (3, 4):
        report = t.report(k)
        assert outcome(NEW, report, t.fixture, summands=4)[0] == 'ok'
        for label, bad, message in tamperings(report):
            got = refused(bad, t.fixture)
            # A shifted index can land on another valid point sum only by accident;
            # the refusal must then still come from the group, not the length.
            if label == 'relation index' and got != message:
                assert got in ('incorrect point relation', 'incorrect scalar-field row'), (k, label, got)
                continue
            assert got == message, (k, label, got)
            if k == 4:
                assert outcome(OLD, bad, t.fixture, summands=4) == ('refused', message), (k, label)


def test_equal_counts_decide_identically_on_random_certificates():
    rng = random.Random(25)
    for trial in range(40):
        t = Toy(toy().fixture, orbits=rng.choice((2, 3)), seed=trial)
        k = rng.choice((3, 4))
        report = t.report(k, relations=rng.randrange(3, 9))
        if rng.random() < 0.5:
            label, report, _ = rng.choice(tamperings(report))
        assert outcome(OLD, report, t.fixture, summands=k) == outcome(NEW, report, t.fixture, summands=k)


# ── Real worker reports (round 0024's frozen workers) ──────────────────────
PAIR_WORKER = ROOT / 'runs/round-0024/worker'
TRIPLE_WORKER = ROOT / 'runs/round-0024/source_candidates/counted/worker'
BASE = {'batch_trials': 1, 'linear_algebra': 'sparse', 'max_trials': 65536}


def _worker_report(worker, case, config):
    job = dict(case['job'], mode='ic', config=config)
    env = {'PATH': os.environ.get('PATH', ''), 'RAYON_NUM_THREADS': '1',
           'IC_ARTIFACT_CACHE': 'off', 'IC_F2_BACKEND': 'cpu'}
    done = subprocess.run([str(worker)], input=json.dumps(job), text=True,
                          capture_output=True, env=env, timeout=600)
    return json.loads(done.stdout)


def test_real_worker_reports_under_a_four_summand_arm():
    if not (PAIR_WORKER.is_file() and TRIPLE_WORKER.is_file()):
        print('  (skipped: round-0024 workers not present)')
        return
    for cell in ('n19a0', 'n37a0'):
        case = next(c for c in FIXTURES['development'] if c['cell'] == cell)
        pair = _worker_report(PAIR_WORKER, case, dict(BASE, solver='pair_table', summands=3))
        triple = _worker_report(TRIPLE_WORKER, case, dict(BASE, solver='triple_counted', summands=4))
        p = outcome(NEW, pair, case['fixture'], summands=4)
        assert p[0] == 'ok' and p[1]['relation_summands'] == 3, p
        assert outcome(OLD, pair, case['fixture'], summands=3)[1] == {
            k: v for k, v in p[1].items() if k != 'relation_summands'}
        q = outcome(NEW, triple, case['fixture'], summands=4)
        assert q[0] == 'ok' and q == outcome(OLD, triple, case['fixture'], summands=4), q
        assert p[1]['solutions'] == q[1]['solutions']
        for report in (pair, triple):
            for label, bad, message in tamperings(report):
                got = refused(bad, case['fixture'])
                if label == 'relation index' and got != message:
                    assert got in ('incorrect point relation', 'incorrect scalar-field row'), (cell, label, got)
                    continue
                assert got == message, (cell, label, got)


# ── Stored receipts: what the old checker accepted ─────────────────────────
RECEIPT_ROUNDS = ('round-0014', 'round-0015', 'round-0016', 'round-0017',
                  'round-0018', 'round-0018b', 'round-0019', 'round-0020')


def test_stored_receipts_verify_identically():
    full = os.environ.get('R25_FULL') == '1'
    every_trial = full or os.environ.get('R25_RECEIPTS') == 'rep0'
    checked = 0
    for name in RECEIPT_ROUNDS:
        rnd = ROOT / 'runs' / name
        if not (rnd / 'runs').is_dir():
            continue
        fixtures = json.loads((rnd / 'fixtures.json').read_text())
        cases = {(stage, c['id']): c for stage in fixtures for c in fixtures[stage]}
        arms = {a['id']: a for a in json.loads((rnd / 'candidates.json').read_text())}
        pattern = '*/*/*/rep-*/receipt.json' if full else '*/*/*/rep-0/receipt.json'
        sampled = set()
        for path in sorted((rnd / 'runs').glob(pattern)):
            receipt = json.loads(path.read_text())
            if receipt['status'] != 'VERIFIED':
                continue
            key = (receipt['stage'], receipt['cell'], receipt['arm'])
            if not every_trial and key in sampled:
                continue
            sampled.add(key)
            case = cases[(receipt['stage'], receipt['case'])]
            summands = arms[receipt['arm']]['config']['summands'] if receipt['arm'] in arms else 3
            for report_path in (['profile', 'native'] if full else ['profile']):
                report = json.loads((path.parent / report_path / 'stdout.json').read_text())
                kw = dict(expected_mode=receipt['mode'], summands=summands)
                old = outcome(OLD, report, case['fixture'], **kw)
                new = outcome(NEW, report, case['fixture'], **kw)
                assert old == new, (str(path), old, new)
                assert old[0] == 'ok', (str(path), old)
                if report_path == 'profile':
                    assert new[1] == receipt['certificate'], str(path)
                checked += 1
    mode = 'every receipt, profile and native' if full else (
        'repetition 0 of every trial, profile' if every_trial else 'one trial per round/stage/cell/arm, profile')
    print(f'  stored reports checked: {checked} ({mode})')


def main():
    tests = [(name, fn) for name, fn in globals().items() if name.startswith('test_') and callable(fn)]
    failed = 0
    for name, fn in tests:
        try:
            fn()
            print(f'PASS {name}')
        except Exception as exc:  # noqa: BLE001 -- report every failure, then exit non-zero
            failed += 1
            print(f'FAIL {name}: {exc!r}')
    print(f'{len(tests) - failed} passed, {failed} failed')
    return 1 if failed else 0


if __name__ == '__main__':
    raise SystemExit(main())
