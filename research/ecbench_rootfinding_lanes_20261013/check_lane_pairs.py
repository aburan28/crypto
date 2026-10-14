#!/usr/bin/env python3
"""Check sealed ecbench serial/batch pairs under an identical lane schedule."""
import json
import sys
from pathlib import Path

session = Path(sys.argv[1])
rows = [json.loads(line) for line in (session / 'records.jsonl').read_text().splitlines()]
rows = [row for row in rows if not row['warmup']]
by_key = {(row['workload']['workload_id'], row['round'], row['arm']): row for row in rows}
keys = sorted({(row['workload']['workload_id'], row['round']) for row in rows})


def phase_native(row, name):
    return next(p for p in row['phases'] if p['name'] == name)['native']


checks = {}
for width in (4, 8, 16):
    serial_name, batch_name = f'b48-serial-{width}', f'b48-batch-{width}'
    outcomes = []
    for key in keys:
        serial = by_key[key + (serial_name,)]
        batch = by_key[key + (batch_name,)]
        fields = ('targets_tried', 'relations_found', 'matrix_rank', 'matrix_rows', 'matrix_work')
        counters_equal = all(serial['counters'][f] == batch['counters'][f] for f in fields)
        rel_s, rel_b = phase_native(serial, 'relations'), phase_native(batch, 'relations')
        native_fields = ('walk_steps', 'walk_restarts', 'repeated_targets_skipped',
                         'lookups', 'target_guard_probes', 'lane_rounds', 'lane_prepared',
                         'lane_unused_prepared')
        native_equal = all(rel_s[f] == rel_b[f] for f in native_fields)
        outcome_equal = (serial['outcome']['recovered'] == batch['outcome']['recovered']
                         and serial['outcome']['status'] == 'verified' and batch['outcome']['status'] == 'verified')
        cost_equal = serial['cost']['total_gae'] == batch['cost']['total_gae']
        # The current ecbench record omits GroupOps.field_ops.  Report that
        # separately instead of treating absent inversion evidence as zero.
        field_available = 'field_ops' in serial and 'field_ops' in batch
        inversion_reduced = (batch['field_ops']['invs'] <= serial['field_ops']['invs']) if field_available else None
        outcomes.append((counters_equal and native_equal and outcome_equal and cost_equal,
                         field_available, inversion_reduced))
    checks[str(width)] = {'pairs': len(outcomes), 'all_exact': all(item[0] for item in outcomes),
                          'field_ops_available': all(item[1] for item in outcomes),
                          'inversion_reduced_when_available': all(item[2] for item in outcomes) if all(item[1] for item in outcomes) else None,
                          'failed_pairs': [i for i, item in enumerate(outcomes) if not item[0]]}
result = {'session': str(session), 'checks': checks,
          'all_pass': all(value['all_exact'] for value in checks.values())}
print(json.dumps(result, indent=2))
if not result['all_pass']:
    raise SystemExit(1)
