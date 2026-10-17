#!/usr/bin/env python3
"""Materialize exact IC1 run rows for three current-source local controls."""
import hashlib
import importlib.util
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE_SHA256 = '8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12'
BINARY_SHA256 = 'd71a9662c9d4e2dde58443bfc646a05737c6b9ca6b3d5e831b632ce09384c906'
CASES = (
    (1248, HERE, 'receipt_K1248.json', 'replay_K1248.json',
     'IC1N83Ckb1a1fb207168PDP4rootRCguidedLAgaussTDdirectISO0hc9eda193920f',
     'f2b9cbc454705add4b8691785fefc8a4498764cc8066e0e6e294ce931f16f685'),
    (1280, HERE, 'receipt_K1280.json', 'replay_K1280.json',
     'IC1N83Ckb1a1fb212480PDP4rootRCguidedLAgaussTDdirectISO0h69261fd42dac',
     'f2b9cbc454705add4b8691785fefc8a4498764cc8066e0e6e294ce931f16f685'),
    (2816, ROOT / 'research/n83_k2816_current_source_20261017',
     'receipt_K2816.json', 'replay_K2816.json',
     'IC1N83Ckb1a1fb467456PDP4rootRCguidedLAgaussTDdirectISO0hf693b1b33c8b',
     '61daed9917be546f46bf72cc8d634525fefbce5a27a01cf72fd56fc5f6cfde63'),
)


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


finegrid = module('smallk_local_workload', ROOT / 'research/n83_e2e_finegrid_20261016/materialize.py')


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def main():
    if digest(ROOT / 'examples/koblitz_orbit_dlp_fast.rs') != SOURCE_SHA256:
        raise ValueError('current source differs')
    rows, workloads = [], []
    for k, study, receipt_name, replay_name, candidate_id, archive_sha in CASES:
        receipt_path, replay_path = study / receipt_name, study / replay_name
        receipt, replay = json.loads(receipt_path.read_text()), json.loads(replay_path.read_text())
        candidate_path = study / f'candidate_{candidate_id}.json'
        candidate = json.loads(candidate_path.read_text())
        archive = study / 'raw.tar.gz'
        workload = finegrid.workload({'hashes': receipt['frozen_input_hashes']})
        phases = replay['online_phases_ms']
        if (receipt['status'] != 'PASS_PRODUCER' or
                replay['status'] != 'PASS_INDEPENDENT_ARITHMETIC_REPLAY' or
                receipt['source_sha256'] != SOURCE_SHA256 or
                receipt['binary_sha256'] != BINARY_SHA256 or
                replay['source_sha256'] != SOURCE_SHA256 or
                replay['binary_sha256'] != BINARY_SHA256 or
                replay['candidate_id'] != candidate_id or
                candidate['candidate_id'] != candidate_id or
                candidate['factor_base']['actual_usable_points'] != 166 * k or
                replay['actual_usable_points'] != 166 * k or
                replay['effective_columns'] != k or
                replay['receipt_sha256'] != digest(receipt_path) or
                digest(archive) != archive_sha or
                not math.isclose(math.fsum(phases.values()),
                                 replay['one_target_online_ms'], abs_tol=1.0)):
            raise ValueError(f'K{k}: source, replay, identity or archive differs')
        rows.append({
            'candidate_id': candidate_id,
            'workload_id': workload['workload_id'],
            'run_id': f"{candidate_id}W{workload['workload_id']}R1",
            'status': 'PASS_INDEPENDENT_ARITHMETIC_REPLAY',
            'hardware_tier': 'shared_host_exploratory_l0',
            'actual_usable_points': 166 * k,
            'effective_columns': k,
            'source_sha256': SOURCE_SHA256,
            'binary_sha256': BINARY_SHA256,
            'target_sha256': receipt['target_sha256'],
            'one_target_online_ms': replay['one_target_online_ms'],
            'online_phases_ms': phases,
            'cold_receipt_wall_ms': receipt['wall_seconds'] * 1000,
            'receipt_sha256': digest(receipt_path),
            'replay_sha256': digest(replay_path),
            'candidate_manifest_sha256': digest(candidate_path),
            'raw_archive_sha256': archive_sha,
            'study': str(study.relative_to(ROOT)),
        })
        workloads.append(workload)
    if len({json.dumps(item, sort_keys=True) for item in workloads}) != 1:
        raise ValueError('exposed one-target workload differs across controls')
    output = {'schema': 'n83.current-source-local-ic1-rows/v1',
              'workload': workloads[0], 'rows': rows}
    path = HERE / 'local_ic1_rows.json'
    blob = json.dumps(output, sort_keys=True, indent=2) + '\n'
    if path.exists() and path.read_text() != blob:
        raise ValueError('versioned local run rows differ')
    path.write_text(blob)
    print(json.dumps({'status': 'PASS', 'workload_id': workloads[0]['workload_id'],
                      'run_ids': [row['run_id'] for row in rows]}, sort_keys=True))


if __name__ == '__main__':
    main()
