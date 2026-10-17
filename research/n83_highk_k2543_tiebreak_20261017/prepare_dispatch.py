#!/usr/bin/env python3
"""Seal the predecessor result binding that arms the high-K cloud workflow."""

import argparse
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('n83_highk_dispatch_run', HERE / 'run.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('predecessor_result_dir', type=Path)
    args = parser.parse_args()
    directory = args.predecessor_result_dir.resolve()
    validated = directory / 'validated_rows.json'
    selection = directory / 'selection.json'
    arms, choice = runner.choose_anchor(validated, selection)
    record = json.loads(validated.read_text())
    if record.get('run_id') is None or record.get('head_sha') is None:
        raise ValueError('predecessor provider run or source head is absent')
    rows_path = directory / f"ic1_rows_{record['run_id']}.json"
    rows = json.loads(rows_path.read_text())
    if (rows.get('run_id') != record['run_id'] or
            rows.get('run_head_sha') != record['head_sha'] or
            rows.get('validated_rows_sha256') != runner.parent.digest(validated) or
            rows.get('selection_sha256') != runner.parent.digest(selection) or
            rows.get('host_strata') != record.get('host_strata') or
            len(rows.get('rows', [])) != 30 or
            len(rows.get('rho_rows', [])) != 30 or
            len(rows.get('workloads', [])) != 10):
        raise ValueError('predecessor materialized provider binding differs')
    dispatch = {
        'schema': 'n83.highk-k2543-tiebreak-dispatch/v1',
        'predecessor_run_id': record['run_id'],
        'predecessor_head_sha': record['head_sha'],
        'predecessor_validated_sha256': runner.parent.digest(validated),
        'predecessor_selection_sha256': runner.parent.digest(selection),
        'predecessor_rows_sha256': runner.parent.digest(rows_path),
        'candidate_ids': [arm[3] for arm in arms],
        'anchor_candidate_id': choice['anchor_candidate_id'],
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
