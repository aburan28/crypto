"""Exact public-point exclusions with retained source provenance."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile

from identity import curve_record
from oracle import require


def file_hash(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def object_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def fixture_values(value):
    if isinstance(value, dict):
        required = {'targets', 'subgroup_order', 'irreducible', 'generator', 'lambda',
                    'cofactor', 'group_order', 'degree', 'curve_a'}
        if required <= value.keys():
            yield value
        else:
            for child in value.values():
                yield from fixture_values(child)
    elif isinstance(value, list):
        for child in value:
            yield from fixture_values(child)


def history_sets(history):
    require(history.get('schema_version') == 1, 'unknown target history schema')
    curves = {}
    for entry in history['curves']:
        key = entry['curve_id']
        require(key not in curves, 'duplicate history curve')
        points = entry['points']
        require(all(isinstance(p, list) and len(p) == 2 and
                    all(type(v) is int and v >= 0 for v in p) for p in points), 'invalid history point')
        frozen = {tuple(p) for p in points}
        require(len(frozen) == len(points), 'duplicate history point')
        curves[key] = frozen
    return curves


def key_for(fixture):
    return curve_record(fixture)['curve']['curve_id']


def reserve_fresh(fixture, used):
    key = key_for(fixture)
    points = [tuple(map(int, p)) for p in fixture['targets']]
    occupied = used.setdefault(key, set())
    if len(set(points)) != len(points) or any(p in occupied for p in points):
        return False
    occupied.update(points)
    return True


def validate_fresh(fixtures, history):
    used = history_sets(history)
    count = 0
    for stage in ('aa', 'smoke', 'development', 'selection', 'confirmation'):
        for case in fixtures.get(stage, []):
            require(reserve_fresh(case['fixture'], used), 'reused public point in '+stage)
            count += len(case['fixture']['targets'])
    if 'replay' in fixtures:
        require(fixtures['replay'] == fixtures.get('confirmation'), 'replay changed confirmation inputs')
    return count


def extend(history, rounds):
    """Reserve all exposed preparation points, including unfinished preparation."""
    result = json.loads(json.dumps(history))
    entries = {row['curve_id']: row for row in result['curves']}
    used = history_sets(result)
    additions = []
    for root in rounds:
        paths = sorted(set(root.glob('fixtures.json')) |
                       set(root.glob('fixture_generation/**/stdout.json')))
        require(paths, 'prior round has no retained fixture evidence')
        hashes = {}
        for path in paths:
            hashes[str(path.relative_to(root))] = file_hash(path)
            for fixture in fixture_values(json.loads(path.read_text())):
                key = key_for(fixture)
                if key not in entries:
                    entries[key] = dict(curve_id=key, cell=f"n{fixture['degree']}a{fixture['curve_a']}",
                                        subgroup_order=int(fixture['subgroup_order']), points=[])
                    used[key] = set()
                used[key].update(tuple(map(int, p)) for p in fixture['targets'])
        additions.append(dict(files=hashes, contract_sha256=file_hash(root/'contract.json')
                              if (root/'contract.json').exists() else None))
    for key, points in used.items():
        entries[key]['points'] = [list(p) for p in sorted(points)]
    result['curves'] = [row for _, row in sorted(entries.items())]
    result['parent_history_sha256'] = object_hash(history)
    result['round_additions'] = additions
    return result


def verify_sources(history, repository):
    """Reconstruct the original corpus from the pinned archive/file inventories."""
    reconstructed = {}
    identities = {}
    def collect(value):
        hits = 0
        for fixture in fixture_values(value):
            identity = {key:fixture[key] for key in ('degree', 'curve_a', 'irreducible',
                'group_order', 'subgroup_order', 'cofactor', 'generator', 'lambda')}
            digest = object_hash(identity)
            if digest not in identities:
                identities[digest] = key_for(fixture)
            key = identities[digest]
            reconstructed.setdefault(key, set()).update(tuple(map(int, p)) for p in fixture['targets'])
            hits += 1
        return hits
    for source in history['sources']:
        if 'archive' in source:
            path = repository/'research/ic_candidate_tournament_20260915/evidence'/source['archive']
            require(file_hash(path) == source['sha256'], 'changed historical archive')
            process = subprocess.Popen(['zstd', '-dq', '-c', str(path)], stdout=subprocess.PIPE)
            members = []
            try:
                with tarfile.open(fileobj=process.stdout, mode='r|') as archive:
                    for member in archive:
                        if not member.isfile() or not member.name.endswith('.json'):
                            continue
                        data = archive.extractfile(member).read()
                        if collect(json.loads(data)):
                            members.append(dict(member=member.name, sha256=hashlib.sha256(data).hexdigest()))
                require(process.wait() == 0, 'history archive decode failed')
            finally:
                process.stdout.close()
                if process.poll() is None:
                    process.terminate()
                process.wait()
            require(len(members) == source['matched_members'] and
                    object_hash(members) == source['member_list_sha256'], 'history member inventory changed')
        else:
            path = repository/source['path']
            require(file_hash(path) == source['sha256'], 'changed historical fixture source')
            collect(json.loads(path.read_text()))
    require(reconstructed == history_sets(history), 'history omitted or changed an exposed target')
    return dict(status='VERIFIED', curves=len(reconstructed), points=sum(map(len, reconstructed.values())))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--history', type=Path, required=True)
    parser.add_argument('--repository', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(verify_sources(json.loads(args.history.read_text()), args.repository.resolve())))


if __name__ == '__main__':
    main()
