#!/usr/bin/env python3
"""Compare deterministic replay evidence while excluding host timings and RSS."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def without_clock(value):
    if isinstance(value, dict):
        return {key: without_clock(item) for key, item in value.items()
                if key not in ('cpu_ns', 'wall_seconds', 'peak_rss_bytes',
                               'timestamp_utc', 'platform', 'python')}
    if isinstance(value, list):
        return [without_clock(item) for item in value]
    return value


def compare(left, right, fields):
    for name in fields:
        assert without_clock(left[name]) == without_clock(right[name]), name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--frozen-toy', type=Path, required=True)
    ap.add_argument('--fresh-toy', type=Path, required=True)
    ap.add_argument('--frozen-exact', type=Path, required=True)
    ap.add_argument('--fresh-exact', type=Path, required=True)
    args = ap.parse_args()
    old, new = (json.loads(p.read_text()) for p in (args.frozen_toy, args.fresh_toy))
    compare(old, new, ('schema', 'source_sha256', 'parameters', 'geometry',
                       'challenge', 'bases', 'base_meta', 'occupancy',
                       'target_streams', 'phase_costs', 'target_prefix_costs',
                       'codomain_target_prefix_costs', 'variants',
                       'cold_cost_to_rank_or_512', 'controls'))
    old, new = (json.loads(p.read_text()) for p in (args.frozen_exact, args.fresh_exact))
    compare(old, new, ('schema', 'status', 'preflight_sha256',
                       'public_challenge_input', 'caps', 'lines'))
    assert len(new['lines']) == 2
    assert all(line['status'] == 'PASS' for line in new['lines'])
    print('PASS: all deterministic toy and exact full-point outcomes and operation counters match')


if __name__ == '__main__':
    main()
