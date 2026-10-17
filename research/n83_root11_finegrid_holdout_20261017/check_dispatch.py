#!/usr/bin/env python3
"""Check every sealed input and event before cloud measurement starts."""

import json
from pathlib import Path

import prepare_dispatch
import run

HERE = Path(__file__).resolve().parent


def check():
    sealed = json.loads((HERE / 'dispatch.json').read_text())
    if sealed != prepare_dispatch.record():
        raise ValueError('dispatch, arm choice, source or fixture differs')
    arms = tuple(tuple(item) for item in sealed['ordered_arms'])
    seen = []
    for stratum in range(4):
        events = run.schedule(stratum, arms)
        if len(events) != 6 * len(arms) * 4:
            raise ValueError('stratum event count differs')
        seen += [event['point_index'] for event in events
                 if event['kind'] == 'ic']
    if sorted(seen) != sorted(index for index in range(24)
                              for _ in arms):
        raise ValueError('point by arm coverage differs')
    return sealed, len(events)


def main():
    sealed, event_count = check()
    print(json.dumps({'status': 'PASS',
                      'candidate_ids': sealed['choice']['ordered_candidate_ids'],
                      'events_per_host': event_count}, sort_keys=True))


if __name__ == '__main__':
    main()
