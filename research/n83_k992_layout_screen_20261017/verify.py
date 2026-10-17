#!/usr/bin/env python3
"""Replay archived N83 K992 layout bases, relations, scalars and file hashes."""

import hashlib
import importlib.util
import json
import math
import tarfile
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE_SHA256 = '8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12'
BINARY_SHA256 = 'd71a9662c9d4e2dde58443bfc646a05737c6b9ca6b3d5e831b632ce09384c906'
TARGET_SHA256 = '2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d'
ARCHIVE_SHA256 = '02993c7b3b7b7bf85af28bf999acfc48b41df2d61afa9408dfde3d9f440d5f55'
RUNS = ('k928root11', 'k960root11', 'k992root11', 'k992root22')


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


archive = HERE / 'raw.tar.gz'
if digest(archive) != ARCHIVE_SHA256:
    raise ValueError('raw archive differs from frozen digest')
manifest_blob = (HERE / 'manifest.json').read_bytes()
manifest = json.loads(manifest_blob)
if (manifest['schema'] != 'n83.local-k992-layout-archive/v1' or
        manifest['runs'] != list(RUNS)):
    raise ValueError('archive manifest differs')
files = manifest['files']
with tempfile.TemporaryDirectory(prefix='n83-k992-layout-replay-') as temporary:
    extracted = Path(temporary)
    seen = set()
    with tarfile.open(archive, 'r:gz') as saved:
        for member in saved:
            name = member.name
            relative = Path(name)
            if (name in seen or relative.is_absolute() or '..' in relative.parts or
                    not member.isfile() or
                    name not in files and name != 'manifest.json'):
                raise ValueError(f'unsafe or unlisted archive member: {name}')
            seen.add(name)
            content = saved.extractfile(member).read()
            if name == 'manifest.json':
                if content != manifest_blob:
                    raise ValueError('manifest member differs')
            elif hashlib.sha256(content).hexdigest() != files[name]:
                raise ValueError(f'archive member digest differs: {name}')
            if name in ('replay.json', 'screen.svg', 'screen.png', 'screen.pdf'):
                if content != (HERE / name).read_bytes():
                    raise ValueError(f'versioned result differs: {name}')
            destination = extracted / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(content)
    if seen != set(files) | {'manifest.json'}:
        raise ValueError('archive inventory incomplete')

    arithmetic = module('n83_k992_layout_arithmetic',
                        ROOT / 'research/n83_root11_dynamic_20261016/verify.py')
    flags = module('n83_k992_layout_flags',
                   ROOT / 'research/n83_root11_dynamic_20261016/run.py')
    target_path = ROOT / 'research/n83_cold_factor_base_20261016/targets.jsonl'
    answer_path = ROOT / 'research/n83_cold_factor_base_20261016/answers.jsonl'
    target = json.loads(target_path.read_text())
    answer = json.loads(answer_path.read_text())
    scalar = int(answer['scalar'])
    point = arithmetic.group.point_from_json(target)
    if (digest(ROOT / 'examples/koblitz_orbit_dlp_fast.rs') != SOURCE_SHA256 or
            digest(target_path) != TARGET_SHA256 or
            answer['point'] != target or
            arithmetic.group.multiply(arithmetic.group.G, scalar) != point):
        raise ValueError('source or exposed point differs')

    rows = []
    for name in RUNS:
        folder = extracted / name
        receipt = json.loads((folder / 'receipt.json').read_text())
        summary = json.loads((folder / 'stdout.jsonl').read_text())
        result = json.loads((folder / 'targets_result.jsonl').read_text())
        base_dump = json.loads((folder / 'base_dump.jsonl').read_text())
        k, layout = receipt['k'], receipt['layout']
        original_folder = Path(receipt['command'][4]).parent
        expected_flags = {key: value for key, value in
                          flags.environment(original_folder, layout).items()
                          if key.startswith('KIC_')}
        if (name != f'k{k}{layout.split("_")[0]}' or
                receipt['schema'] != 'n83.local-k992-layout-screen/v1' or
                receipt['status'] != 'EXIT_0' or receipt['exit_code'] != 0 or
                receipt['source_sha256'] != SOURCE_SHA256 or
                receipt['binary_sha256'] != BINARY_SHA256 or
                receipt['target_sha256'] != TARGET_SHA256 or
                receipt['environment_flags'] != expected_flags or
                Path(receipt['command'][0]).name != 'koblitz_orbit_dlp_fast' or
                receipt['command'][1] != f'construct:83:1:{k}' or
                Path(receipt['command'][2]).parts[-3:] !=
                ('research', 'n83_cold_factor_base_20261016', 'targets.jsonl') or
                receipt['command'][3] != '29' or
                original_folder.name != name or
                Path(receipt['command'][4]).name != 'targets_result.jsonl' or
                summary['orbit_columns'] != k or summary['rank'] != k or
                summary['factor_base_points'] != 166 * k or
                summary['targets_solved'] != 1 or summary['targets_failed'] != 0 or
                result['target'] != target or result['group_verified'] is not True or
                int(result['recovered_scalar']) != scalar or
                (folder / 'stderr.log').read_bytes()):
            raise ValueError(f'{name}: source, command, rank or scalar differs')
        base = arithmetic.group.replay_base(base_dump, summary)
        indices = result['point_indices']
        if (len(base) != 166 * k or len(indices) != 4 or any(
                not isinstance(i, int) or not 0 <= i < len(base) for i in indices)):
            raise ValueError(f'{name}: base or witness differs')
        relation = None
        for index in indices:
            relation = arithmetic.group.add(relation, base[index])
        if relation != point:
            raise ValueError(f'{name}: target relation differs')
        phases = [float(result[key]) for key in (
            'target_query_ms', 'target_pdp_ms', 'target_relation_check_ms',
            'target_descent_ms', 'target_recovery_check_ms')]
        online = float(result['target_ms'])
        if (not math.isfinite(online) or online <= 0 or
                any(not math.isfinite(v) or v < 0 for v in phases) or
                abs(math.fsum(phases) - online) > 1):
            raise ValueError(f'{name}: online phase sum differs')
        rows.append({
            'name': name, 'k': k, 'layout': layout,
            'status': 'PASS_INDEPENDENT_ARITHMETIC_REPLAY',
            'base_hash': summary['base_hash'],
            'actual_usable_points': len(base), 'effective_columns': k,
            'index_build_ms': summary['timing_ms']['index_build'],
            'rank_stage_ms': summary['timing_ms']['rank_stage'],
            'cold_process_ms': summary['timing_ms']['process_total'],
            'target_online_ms': online, 'target_probes': result['probes'],
            'receipt_sha256': digest(folder / 'receipt.json'),
            'stdout_sha256': digest(folder / 'stdout.jsonl'),
            'base_dump_sha256': digest(folder / 'base_dump.jsonl'),
            'target_result_sha256': digest(folder / 'targets_result.jsonl'),
        })
    replay = {'schema': 'n83.local-k992-layout-replay/v1',
              'source_sha256': SOURCE_SHA256,
              'binary_sha256': BINARY_SHA256,
              'target_sha256': TARGET_SHA256,
              'answer_sha256': digest(answer_path),
              'status': 'PASS_INDEPENDENT_ARITHMETIC_REPLAY',
              'rows': rows}
    if replay != json.loads((HERE / 'replay.json').read_text()):
        raise ValueError('independent replay differs from versioned result')

print(json.dumps({'status': replay['status'], 'rows': len(rows),
                  'archive_sha256': ARCHIVE_SHA256}, sort_keys=True))
