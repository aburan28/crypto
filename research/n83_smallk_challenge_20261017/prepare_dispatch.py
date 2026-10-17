#!/usr/bin/env python3
"""Seal the finalist result binding that arms the five-arm cloud workflow."""

import argparse
import json
from pathlib import Path

import importlib.util

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('n83_k3072_dispatch_run', HERE / 'run.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('finalist_result_dir', type=Path)
    args = parser.parse_args()
    directory = args.finalist_result_dir.resolve()
    validated = directory / 'validated_rows.json'
    selection = directory / 'selection.json'
    _, choice = runner.choose_second(validated, selection)
    record = json.loads(validated.read_text())
    if record.get('run_id') is None or record.get('head_sha') is None:
        raise ValueError('finalist provider run or source head is absent')
    rows_path = directory / f"ic1_rows_{record['run_id']}.json"
    rows = json.loads(rows_path.read_text())
    if (rows.get('run_id') != record['run_id'] or
            rows.get('run_head_sha') != record['head_sha'] or
            rows.get('validated_rows_sha256') != runner.parent.digest(validated) or
            rows.get('selection_sha256') != runner.parent.digest(selection) or
            rows.get('host_strata') != record.get('host_strata') or
            len(rows.get('rows', [])) != 132 or
            len(rows.get('rho_rows', [])) != 128 or
            len(rows.get('workloads', [])) != 16):
        raise ValueError('finalist materialized provider binding differs')
    dispatch = {
        'schema': 'n83.smallk-challenge-dispatch/v1',
        'finalist_run_id': record['run_id'],
        'finalist_head_sha': record['head_sha'],
        'finalist_validated_sha256': runner.parent.digest(validated),
        'finalist_selection_sha256': runner.parent.digest(selection),
        'finalist_rows_sha256': runner.parent.digest(rows_path),
        'comparison_candidate_id': choice['comparison_candidate_id'],
        'k1248_candidate_id': runner.K1248[3],
        'k1280_candidate_id': runner.K1280[3],
        'k2543_candidate_id': runner.K2543[3],
        'k3072_candidate_id': runner.K3072[3],
        'fallback_used': choice['fallback_used'],
    }
    path = HERE / 'dispatch.json'
    blob = json.dumps(dispatch, sort_keys=True, indent=2) + '\n'
    if path.exists() and path.read_text() != blob:
        raise ValueError('sealed dispatch differs')
    path.write_text(blob)
    print(json.dumps(dispatch, sort_keys=True))


if __name__ == '__main__':
    main()
