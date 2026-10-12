"""Reconstruct every point the censored 2026092902 run could have exposed.

Run with an ic_tournament_worker built from Git commit 765c3c5f. The result
is lost-v2-campaign-exposures.json. This does not redispatch the timed-out
campaign. Exclusions include sealed improvement history, supplemental EXPOSED
corpora, and the reconstructed first-run corpus (lost-campaign-exposures.json).
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
SOURCE_RUN_ID = 36580669479
PRIOR_LOST = HERE / 'lost-campaign-exposures.json'
PRIOR_LOST_SHA256 = 'a728677b199eac02800d8338204d5306f391ec5da757c910bec1e51955fe7b41'
PANEL = HERE / 'panel.json'
PANEL_SHA256 = 'd283a869b0412228d1c66260fdfd8f387d7243bd15456c7febf3c46ee5da27a8'


def extract_census(archive_path, root):
    """Read only historical target evidence from a checked archive stream."""
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
                path = root/suffix
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


def file_sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def reconstruct(worker, temporary):
    inventory = json.loads((TOURNAMENT/'evidence/manifest.json').read_text())
    hashes = {entry['file']:entry['sha256'] for entry in inventory['archives']}
    rounds = []
    for label in ARCHIVES:
        path = TOURNAMENT/'evidence'/(label+'.tar.zst')
        require(file_sha256(path) == hashes[path.name], 'sealed historical archive changed')
        destination = temporary/label
        extract_census(path, destination)
        rounds.append(destination)
    history = json.loads((TOURNAMENT/'goal_20260924/improvement/target-history.json').read_text())
    known = extend(history, rounds)
    require(file_sha256(PRIOR_LOST) == PRIOR_LOST_SHA256,
            'changed frozen first-run target exposures')
    prior = json.loads(PRIOR_LOST.read_text())
    require(prior['seed'] == 2026092901 and prior['source_run_id'] == 36532455386
            and len(prior['attempts']) == 25
            and all(row['accepted'] for row in prior['attempts']),
            'incomplete first-run target reconstruction')
    # Match run_generic_backend_qualification_v2.py: EXPOSED then prior lost.
    known = freeze_exposures(known, list(EXPOSED) + [PRIOR_LOST],
                             temporary/'known-exposures')

    for case in json.loads((rounds[2]/'fixtures.json').read_text())['aa']:
        job = dict(case['job'])
        job.pop('public_targets', None)
        require(fixture(worker, job) == case['fixture'],
                'worker disagrees with a retained prior fixture: '+case['id'])

    require(file_sha256(PANEL) == PANEL_SHA256, 'changed premeasurement panel bytes')
    panel = json.loads(PANEL.read_text())
    require(panel['seed'] == SEED and panel['repetitions'] == 1
            and panel['cells'] == ['n17a1', 'n19a0', 'n23a0', 'n23a1', 'n31a0'],
            'second registration changed')
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
                        factor_base=dict(kind='subgroup_orbits', seed=43, points=6*degree),
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
                    raise ValueError('point generation exhausted for '+case_id)
    accepted = [row for row in attempts if row['accepted']]
    require(len(accepted) == 25, 'changed second-run distinct-point schedule')
    # No accepted v2 point may equal an accepted first-run point on the same curve.
    prior_points = {(row['fixture']['degree'], row['fixture']['curve_a'],
                     tuple(map(int, row['fixture']['targets'][0])))
                    for row in prior['attempts'] if row['accepted']}
    for row in accepted:
        key = (row['fixture']['degree'], row['fixture']['curve_a'],
               tuple(map(int, row['fixture']['targets'][0])))
        require(key not in prior_points, 'v2 reconstruction reused a first-run point')
    return dict(schema_version=1, source_commit=SOURCE, seed=SEED,
                source_run_id=SOURCE_RUN_ID,
                prior_censored_exposures_sha256=PRIOR_LOST_SHA256,
                registration_panel_sha256=PANEL_SHA256,
                reconstruction=('deterministic frozen fixture generator with sealed '
                                'history, supplemental exposures, and first-run '
                                'lost-campaign exclusions'),
                attempts=attempts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as name:
        result = reconstruct(args.worker.resolve(), Path(name))
    args.out.write_text(json.dumps(result, sort_keys=True, separators=(',', ':'))+'\n')
    print(hashlib.sha256(args.out.read_bytes()).hexdigest())


if __name__ == '__main__':
    main()
