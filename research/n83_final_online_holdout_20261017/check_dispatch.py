#!/usr/bin/env python3
"""Check sealed final holdout dispatch before allocating measurement work."""

import json
from pathlib import Path

import prepare_dispatch

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def check():
    path = HERE / 'dispatch.json'
    sealed = json.loads(path.read_text())
    if sealed.get('schema') != 'n83.final-online-holdout-dispatch/v1':
        raise ValueError('final holdout dispatch schema differs')
    paths = tuple((ROOT / name).resolve() for name in sealed['panel_paths'])
    if len(paths) != 4 or any(not path.is_relative_to(ROOT) for path in paths):
        raise ValueError('upstream panel paths escape repository')
    expected = prepare_dispatch.record(paths)
    if sealed != expected:
        raise ValueError('sealed dispatch, upstream evidence or source differs')
    for half in (0, 1):
        import run
        events = run.schedule(half, tuple(tuple(arm) for arm in sealed['ordered_arms']))
        if len(events) != 6 * len(sealed['ordered_arms']) * 4:
            raise ValueError('frozen event count differs')
    return sealed, len(events)


def main():
    sealed, event_count = check()
    print(json.dumps({'status': 'PASS',
                      'candidate_ids': sealed['choice']['ordered_candidate_ids'],
                      'panel_sha256': sealed['panel_sha256'],
                      'events_per_host': event_count}, sort_keys=True))


if __name__ == '__main__':
    main()
