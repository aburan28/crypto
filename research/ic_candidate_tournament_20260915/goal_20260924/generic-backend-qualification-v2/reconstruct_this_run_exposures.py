"""Reconstruct every public point seed 2026092902 could have exposed.

The measured v2 campaign lost its archive. Later registrations must exclude the
25 points this seed's prepare path would generate after the sealed history,
supplemental corpora and the first censored run's points. This command never
reruns measurement; it only freezes the deterministic fixture schedule.
"""
import argparse
import hashlib
import json
from pathlib import Path
import random
import subprocess
import sys
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent
TOURNAMENT = HERE.parents[1]
ROOT = TOURNAMENT.parents[1]
sys.path.insert(0, str(TOURNAMENT))

from oracle import Curve, require
from run_improvement_v3 import EXPOSED
from target_history import extend, freeze_exposures, history_sets, key_for

ARCHIVES = ('ic-improvement-round1-20260925',
            'ic-improvement-round2-20260928',
            'ic-improvement-round3-20260928')
SOURCE = '765c3c5f19032bd852163805f257c56babef2040'
SEED = 2026092902
PRIOR = HERE / 'lost-campaign-exposures.json'
PRIOR_SHA256 = 'a728677b199eac02800d8338204d5306f391ec5da757c910bec1e51955fe7b41'


def extract_census(archive_path, root):
    decoder = subprocess.Popen(['zstd', '-dc', str(archive_path)],
                               stdout=subprocess.PIPE)
    count = 0
    try:
        with tarfile.open(fileobj=decoder.stdout, mode='r|') as archive:
            for member in archive:
                if not member.isfile() or '/round/tournament/' not in member.name:
                    continue
                suffix = member.name.split('/round/tournament/', 1)[1]
                if suffix not in ('fixtures.json', 'target-history.json', 'contract.json') and not (
                    suffix.startswith('fixture_generation/') and suffix.endswith('/stdout.json')):
                    continue
                path = root / suffix
                require(path.resolve().is_relative_to(root.resolve()), 'unsafe archive member')
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(archive.extractfile(member).read())
                count += 1
    finally:
        decoder.stdout.close()
    require(decoder.wait() == 0 and count >= 3, 'historical census extraction failed')


def fixture(worker, job):
    process = subprocess.run([str(worker)], input=json.dumps(job), text=True,
                             capture_output=True, check=True)
    result = json.loads(process.stdout)
    require(result['status'] == 'fixture', 'worker did not generate a fixture')
    Curve(result['fixture'])
    return result['fixture']


def reconstruct(worker, temporary):
    inventory = json.loads((TOURNAMENT / 'evidence/manifest.json').read_text())
    hashes = {entry['file']: entry['sha256'] for entry in inventory['archives']}
    rounds = []
    for label in ARCHIVES:
        path = TOURNAMENT / 'evidence' / (label + '.tar.zst')
        with path.open('rb') as stream:
            require(hashlib.file_digest(stream, 'sha256').hexdigest() == hashes[path.name],
                    'sealed historical archive changed')
        destination = temporary / label
        extract_census(path, destination)
        rounds.append(destination)
    history = json.loads((TOURNAMENT / 'goal_20260924/improvement/target-history.json').read_text())
    known = extend(history, rounds)
    supplemental = list(EXPOSED) + [PRIOR]
    require(hashlib.sha256(PRIOR.read_bytes()).hexdigest() == PRIOR_SHA256,
            'first-run exposure corpus changed')
    known = freeze_exposures(known, supplemental, temporary / 'known-exposures')

    for case in json.loads((rounds[2] / 'fixtures.json').read_text())['aa']:
        job = dict(case['job'])
        job.pop('public_targets', None)
        require(fixture(worker, job) == case['fixture'],
                'worker disagrees with a retained prior fixture: ' + case['id'])

    panel = json.loads((HERE / 'panel.json').read_text())
    require(panel['seed'] == SEED and panel['repetitions'] == 1
            and panel['cells'] == ['n17a1', 'n19a0', 'n23a0', 'n23a1', 'n31a0'],
            'second registration panel changed')
    used = history_sets(known)
    rng = random.Random(SEED)
    attempts = []
    for stage, count in (('aa', 1), ('smoke', 1), ('development', 3)):
        for cell in panel['cells']:
            degree, curve_a = map(int, cell[1:].split('a'))
            for index in range(count):
                case_id = f'{cell}-{index:03d}'
                algorithm_seed = rng.getrandbits(64)
                for number in range(1000):
                    seed = rng.getrandbits(64)
                    job = dict(mode='fixture', degree=degree, curve_a=curve_a,
                               target_seeds=[seed], algorithm_seed=algorithm_seed,
                               factor_base=dict(kind='subgroup_orbits', seed=43,
                                                points=6 * degree),
                               config=panel['candidates'][0]['config'])
                    generated = fixture(worker, job)
                    point = tuple(map(int, generated['targets'][0]))
                    occupied = used.setdefault(key_for(generated), set())
                    accepted = point not in occupied
                    attempts.append(dict(stage=stage, case=case_id, attempt=number,
                                         target_seed=seed, accepted=accepted,
                                         fixture=generated))
                    if accepted:
                        occupied.add(point)
                        break
                else:
                    raise ValueError('point generation exhausted for ' + case_id)
    require(len([row for row in attempts if row['accepted']]) == 25,
            'changed second-run distinct-point schedule')
    first_run = {tuple(map(int, row['fixture']['targets'][0]))
                 for row in json.loads(PRIOR.read_text())['attempts'] if row['accepted']}
    second = {tuple(map(int, row['fixture']['targets'][0]))
              for row in attempts if row['accepted']}
    require(not (first_run & second), 'second-run schedule collided with first-run exposures')
    return dict(schema_version=1, source_commit=SOURCE, seed=SEED,
                source_run_id=36580669479,
                prior_censored_exposures_sha256=PRIOR_SHA256,
                reconstruction=('deterministic frozen fixture generator with sealed history, '
                                'supplemental exposures and first-run censored points'),
                attempts=attempts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as name:
        result = reconstruct(args.worker.resolve(), Path(name))
    args.out.write_text(json.dumps(result, sort_keys=True, separators=(',', ':')) + '\n')
    print(hashlib.sha256(args.out.read_bytes()).hexdigest())


if __name__ == '__main__':
    main()
