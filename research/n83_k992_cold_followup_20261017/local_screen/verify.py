#!/usr/bin/env python3
"""Independent arithmetic replay of exposed-point lower-K N83 controls."""

import hashlib
import importlib.util
import json
import math
import tarfile
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE_SHA256 = '8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12'
BINARY_SHA256 = 'd71a9662c9d4e2dde58443bfc646a05737c6b9ca6b3d5e831b632ce09384c906'
ARCHIVE_SHA256 = '9269fde6cd09fffd8866396be119c05356e75e3975f660d4d2c88a54cf8f0bc0'


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


archive = HERE / 'raw.tar.gz'
if digest(archive) != ARCHIVE_SHA256:
    raise ValueError('raw archive differs from frozen digest')
expected_files = json.loads((HERE / 'manifest.json').read_text())['files']
temporary = tempfile.TemporaryDirectory(prefix='n83-lower-k-replay-')
run_root = Path(temporary.name)
with tarfile.open(archive, 'r:gz') as saved:
    names = set()
    for member in saved:
        name = member.name
        relative = Path(name)
        if (name in names or relative.is_absolute() or '..' in relative.parts
                or not member.isfile() or
                name not in expected_files and name != 'manifest.json'):
            raise ValueError(f'unsafe or unlisted archive member: {name}')
        names.add(name)
        content = saved.extractfile(member).read()
        if name == 'manifest.json':
            if content != (HERE / 'manifest.json').read_bytes():
                raise ValueError('archive manifest differs')
        elif hashlib.sha256(content).hexdigest() != expected_files[name]:
            raise ValueError(f'archive member digest differs: {name}')
        out = run_root / relative
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(content)
    if names != set(expected_files) | {'manifest.json'}:
        raise ValueError('archive member inventory incomplete')


verifier = module('lower_k_group', ROOT / 'research/n83_root11_dynamic_20261016/verify.py')
group = verifier.group
target_path = ROOT / 'research/n83_cold_factor_base_20261016/targets.jsonl'
answer_path = ROOT / 'research/n83_cold_factor_base_20261016/answers.jsonl'
target = json.loads(target_path.read_text())
answer = json.loads(answer_path.read_text())
point = group.point_from_json(target)
scalar = int(answer['scalar'])
assert digest(ROOT / 'examples/koblitz_orbit_dlp_fast.rs') == SOURCE_SHA256
assert group.multiply(group.G, scalar) == point

rows = []
for k in (768, 896, 960, 992, 1008, 1024, 1056, 1107):
    folder = run_root / f'k{k}root11'
    receipt = json.loads((folder / 'receipt.json').read_text())
    summary = json.loads((folder / 'stdout.jsonl').read_text())
    result = json.loads((folder / 'targets_result.jsonl').read_text())
    base_dump = json.loads((folder / 'base_dump.jsonl').read_text())
    if (receipt['k'] != k or receipt['status'] != 'EXIT_0' or
            receipt['exit_code'] != 0 or receipt['source_sha256'] != SOURCE_SHA256 or
            receipt['binary_sha256'] != BINARY_SHA256 or
            receipt['target_sha256'] != digest(target_path) or
            summary['orbit_columns'] != k or summary['rank'] != k or
            summary['factor_base_points'] != 166 * k or
            summary['targets_solved'] != 1 or summary['targets_failed'] != 0 or
            result['target'] != target or result['group_verified'] is not True or
            int(result['recovered_scalar']) != scalar or
            (folder / 'stderr.log').read_bytes()):
        raise ValueError(f'K{k}: producer, source, target or scalar mismatch')
    base = group.replay_base(base_dump, summary)
    indices = result['point_indices']
    if len(base) != 166 * k or len(indices) != 4 or any(
            not isinstance(i, int) or not 0 <= i < len(base) for i in indices):
        raise ValueError(f'K{k}: base or witness mismatch')
    relation = None
    for index in indices:
        relation = group.add(relation, base[index])
    if relation != point:
        raise ValueError(f'K{k}: relation mismatch')
    phases = [float(result[key]) for key in (
        'target_query_ms', 'target_pdp_ms', 'target_relation_check_ms',
        'target_descent_ms', 'target_recovery_check_ms')]
    online = float(result['target_ms'])
    if (not math.isfinite(online) or online <= 0 or
            any(not math.isfinite(v) or v < 0 for v in phases) or
            abs(math.fsum(phases) - online) > 1):
        raise ValueError(f'K{k}: online phase mismatch')
    rows.append({
        'k': k, 'status': 'PASS_INDEPENDENT_ARITHMETIC_REPLAY',
        'actual_usable_points': len(base),
        'effective_columns': k,
        'root_table_entries': summary['root_table_entries'],
        'base_hash': summary['base_hash'],
        'index_build_ms': summary['timing_ms']['index_build'],
        'rank_stage_ms': summary['timing_ms']['rank_stage'],
        'cold_process_ms': summary['timing_ms']['process_total'],
        'target_online_ms': online,
        'target_probes': result['probes'],
        'receipt_sha256': digest(folder / 'receipt.json'),
        'stdout_sha256': digest(folder / 'stdout.jsonl'),
        'base_dump_sha256': digest(folder / 'base_dump.jsonl'),
        'target_result_sha256': digest(folder / 'targets_result.jsonl'),
    })

output = {'schema': 'n83.local-lower-k-replay/v1',
          'source_sha256': SOURCE_SHA256, 'binary_sha256': BINARY_SHA256,
          'target_sha256': digest(target_path), 'answer_sha256': digest(answer_path),
          'status': 'PASS_INDEPENDENT_ARITHMETIC_REPLAY', 'rows': rows}
if output != json.loads((HERE / 'replay.json').read_text()):
    raise ValueError('independent replay differs from versioned result')
temporary.cleanup()
print(json.dumps({'status': output['status'], 'rows': len(rows)}, sort_keys=True))
