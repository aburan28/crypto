#!/usr/bin/env python3
"""Deterministically freeze positive holdout Qs; no rank or log outcome."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import resource
import signal
import sys
import tarfile
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent / 'rotated_row_certificate_20260925'
ARCHIVE = PRIOR / 'evidence/raw.tar.gz'
ARCHIVE_SHA = '863f50ea872ce01b1999d0ef67bc0341fa0c764f2ad4187d16736768b0ab9b3d'
OLD_GATE = HERE.parent / 'rotated_subspace_support_20260925/gate.py'
DOMAIN = 'ECC2K130-ROTATED-RANK-20260925-v1'
ARMS = ((3, 5), (3, 6), (7, 5), (7, 6))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prior():
    spec = importlib.util.spec_from_file_location('frozen_old_gate_rank_input', OLD_GATE)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def peak_rss():
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return raw if sys.platform == 'darwin' else raw * 1024


def alarm(_signum, _frame):
    raise TimeoutError('holdout input selection exceeded 30-second wall cap')


def save(path, value):
    path.write_text(json.dumps(value, sort_keys=True, separators=(',', ':')) + '\n')


def score(beta, m, k):
    return hashlib.sha256(f'{DOMAIN}/{beta}/{m}/{k}'.encode('ascii')).digest()


def generate():
    assert digest(ARCHIVE) == ARCHIVE_SHA
    old = prior()
    f = old.Field(13, old.MODELS[13]['low'])
    curve = old.Curve(f)
    assert old.source_group_order(13) == 8012
    arms_manifest, point_arms, label_arms = [], [], []
    with tarfile.open(ARCHIVE, 'r:gz') as tar:
        for beta, m in ARMS:
            name = f'raw/n13-b{beta}-m{m}-roster.jsonl'
            member = tar.getmember(name)
            assert member.isfile() and member.size < 10_000_000
            handle = tar.extractfile(member)
            assert handle is not None
            roster = [json.loads(line) for line in handle.read().decode().splitlines()]
            assert len(roster) == 2003 and [r['k'] for r in roster] == list(range(2003))
            H = tuple(roster[1]['point'])
            assert H == (4793, 2429) and curve.on_curve(H)
            assert curve.scalar(H, 2003) is None
            selected, scanned = [], []
            for k in sorted(range(1, 2003), key=lambda x: (score(beta, m, x), x)):
                record = roster[k]
                point = record['point']
                assert tuple(point) == curve.scalar(H, k)
                cosets = record['coset_witness_indices']
                assert len(cosets) == 4
                positive = any(index is not None for index in cosets)
                item = {'k': k, 'score_sha256': score(beta, m, k).hex(),
                        'positive': positive, 'positive_cosets': sum(x is not None for x in cosets)}
                scanned.append(item)
                if positive:
                    case_id = f'b{beta}m{m}-{len(selected):02d}'
                    selected.append({'case_id': case_id, 'k': k, 'point': point})
                    if len(selected) == 16:
                        break
            assert len(selected) == 16
            arms_manifest.append({'beta': beta, 'm': m, 'scanned_candidates': scanned,
                                  'selected_case_ids': [x['case_id'] for x in selected],
                                  'selection_positive_count': 16,
                                  'selection_negative_count': sum(not x['positive'] for x in scanned)})
            point_arms.append({'beta': beta, 'm': m,
                               'cases': [{'case_id': x['case_id'], 'point': x['point']} for x in selected]})
            label_arms.append({'beta': beta, 'm': m,
                               'labels': [{'case_id': x['case_id'], 'k': x['k']} for x in selected]})
    manifest = {'schema': 'rotated_rank_holdout_v1', 'domain': DOMAIN,
                'source_archive_sha256': ARCHIVE_SHA, 'q': 2003,
                'generator': [4793, 2429], 'arms': arms_manifest}
    point_only = {'schema': 'rotated_rank_point_only_v1', 'arms': point_arms}
    sealed = {'schema': 'rotated_rank_sealed_labels_v1', 'arms': label_arms}
    return manifest, point_only, sealed, {'field_operations': dict(f.operations),
                                         'curve_operations': dict(curve.operations)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    assert not args.out.exists()
    args.out.mkdir(parents=True)
    begin, cpu = time.monotonic(), time.process_time()
    receipt = {'kind': 'holdout_input_construction', 'status': 'failed',
               'source_archive_sha256': ARCHIVE_SHA}
    signal.signal(signal.SIGALRM, alarm)
    signal.setitimer(signal.ITIMER_REAL, 30)
    try:
        manifest, points, sealed, operations = generate()
        for name, data in [('manifest.json', manifest), ('point_only.json', points),
                           ('sealed_labels.json', sealed)]:
            save(args.out / name, data)
        receipt['operations'] = operations
        receipt['input_sha256'] = {name: digest(args.out / name)
                                   for name in ('manifest.json', 'point_only.json', 'sealed_labels.json')}
        receipt['status'] = 'success'
    except Exception as exc:
        receipt['error'] = repr(exc)
        receipt['traceback'] = traceback.format_exc()
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        receipt['wall_seconds'] = time.monotonic() - begin
        receipt['cpu_seconds'] = time.process_time() - cpu
        receipt['peak_rss_bytes'] = peak_rss()
        save(args.out / 'construction_receipt.json', receipt)
        assert receipt['wall_seconds'] <= 30 and receipt['peak_rss_bytes'] <= 512 * 1024 * 1024


if __name__ == '__main__':
    main()
