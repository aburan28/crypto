#!/usr/bin/env python3
"""Charged lookup of saved witness rows, stripped of target logarithm labels."""
from __future__ import annotations

import argparse
import hashlib
import json
import resource
import sys
import tarfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent / 'rotated_row_certificate_20260925/evidence/raw.tar.gz'
SHA = '863f50ea872ce01b1999d0ef67bc0341fa0c764f2ad4187d16736768b0ab9b3d'
ARMS = ((3, 5), (3, 6), (7, 5), (7, 6))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rss():
    v = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return v if sys.platform == 'darwin' else v * 1024


def save(path, value):
    path.write_text(json.dumps(value, sort_keys=True, separators=(',', ':')) + '\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--points', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    start, cpu = time.monotonic(), time.process_time()
    assert digest(ARCHIVE) == SHA
    points = json.loads(args.points.read_text())
    assert [(a['beta'], a['m']) for a in points['arms']] == list(ARMS)
    output = []
    scans = []
    with tarfile.open(ARCHIVE, 'r:gz') as tar:
        for (beta, m), arm in zip(ARMS, points['arms']):
            stem = f'raw/n13-b{beta}-m{m}-rows.jsonl'
            member = tar.getmember(stem)
            assert member.isfile() and member.size < 50_000_000
            handle = tar.extractfile(member)
            assert handle is not None
            source = handle.read()
            wanted = {tuple(case['point']): case['case_id'] for case in arm['cases']}
            assert len(wanted) == len(arm['cases']) == 16
            found = {}
            scanned = 0
            for line in source.splitlines():
                row = json.loads(line)
                scanned += 1
                target = None if row['target_Q'] is None else tuple(row['target_Q'])
                if target in wanted and target not in found:
                    # Do not copy k, source indices, witness points, or signed term telemetry.
                    found[target] = {'point': list(target), 'torsion_index': row['torsion_index'],
                                     'row': row['row']}
            responses = []
            for case in arm['cases']:
                target = tuple(case['point'])
                responses.append({'case_id': case['case_id'], 'point': case['point'],
                                  'lookup_status': 'hit' if target in found else 'miss',
                                  'response': found.get(target)})
            output.append({'beta': beta, 'm': m, 'responses': responses})
            scans.append({'beta': beta, 'm': m, 'archive_rows_scanned': scanned,
                          'archive_bytes_read': len(source), 'hits': len(found),
                          'misses': 16 - len(found)})
    result = {'schema': 'rotated_archive_oracle_v1', 'source_archive_sha256': SHA,
              'point_only_sha256': digest(args.points), 'arms': output,
              'scan_cost': {'arms': scans, 'total_wall_seconds': time.monotonic() - start,
                            'total_cpu_seconds': time.process_time() - cpu,
                            'peak_rss_bytes': rss()}}
    save(args.out, result)
    assert time.monotonic() - start <= 120 and rss() <= 512 * 1024 * 1024


if __name__ == '__main__':
    main()
