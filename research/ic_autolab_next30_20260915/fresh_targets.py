import json
import random
import secrets
import subprocess
from pathlib import Path

def target_key(fixture):
    return (int(fixture['degree']), int(fixture['curve_a']),
            tuple(tuple(int(x) for x in point) for point in fixture['targets']))

def choose_seed(t, baseline_binary, out):
    prior = t.read(Path('/opt/previous-targets.json'))
    excluded = {(int(d), int(a), tuple(tuple(p) for p in points))
                for d, a, points in prior['target_keys']}
    attempts = []
    for _ in range(8):
        seed = secrets.randbits(63)
        rng = random.Random(seed)
        seen = set(excluded)
        confirmation = []
        collisions = 0
        for stage in t.STAGES:
            if stage == 'replay':
                continue
            cells = [(13, 0), (17, 1), (19, 0), (23, 0)]
            if stage == 'confirmation':
                cells.append((19, 1))
            count = {'development': 3, 'selection': 3, 'confirmation': 12}.get(stage, 1)
            for degree, a in cells:
                for _index in range(count):
                    public_seed = rng.getrandbits(64)
                    job = {'mode': 'fixture', 'degree': degree, 'curve_a': a,
                           'target_seeds': [public_seed], 'algorithm_seed': rng.getrandbits(64),
                           'factor_base': {'kind': 'subgroup_orbits', 'seed': 43, 'points': 6*degree},
                           'config': t.BASE_CONFIG}
                    raw = subprocess.run([str(baseline_binary)], input=json.dumps(job),
                        text=True, capture_output=True, check=True, timeout=30, env=t.child_env())
                    fixture = json.loads(raw.stdout)['fixture']
                    key = target_key(fixture)
                    if stage == 'confirmation':
                        collisions += int(key in seen)
                        confirmation.append(key)
                    seen.add(key)
        attempts.append({'seed': seed, 'confirmation_collisions': collisions})
        if collisions == 0:
            record = {'seed': seed, 'generated_after_candidate_locked': True,
                      'prior_target_keys': len(excluded), 'attempts': attempts,
                      'confirmation_target_keys': confirmation,
                      'previous_targets_sha256': t.digest(Path('/opt/previous-targets.json')),
                      'selection_uses_only_freshness': True}
            t.write(out/'holdout-seed.json', record, exclusive=True)
            return seed
    raise RuntimeError('No disjoint confirmation set in eight preparation-only seed attempts')

def verify_prepared_targets(t, round_dir, out):
    recorded = t.read(out/'holdout-seed.json')['confirmation_target_keys']
    expected = [(d, a, tuple(tuple(p) for p in points)) for d, a, points in recorded]
    actual = [target_key(c['fixture']) for c in t.read(round_dir/'fixtures.json')['confirmation']]
    t.require(actual == expected, 'prepared confirmation differs from freshness precheck')
