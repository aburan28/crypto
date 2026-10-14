#!/usr/bin/env python3
"""Bind fresh fixed-prefix inputs to the frozen ecbench workload plan."""
import hashlib
import json
from pathlib import Path

study = Path(__file__).parent
plan_file = study / 'fixed-prefix-fresh-plan.json'
plan = json.loads(plan_file.read_text())
spec = json.loads((study / 'specs/fixed-prefix-fresh.json').read_text())
inputs = json.loads((study / 'fixed-prefix-fresh-input.json').read_text())
old = [json.loads(line) for line in (study / 'sessions/coverage-holdout-local-20261013/records.jsonl').read_text().splitlines()]
old_points = {tuple(row['workload']['target']) for row in old}

assert spec['workloads']['target_seed'] == 20261025
assert spec['measurement']['warmup'] == 0 and spec['measurement']['rounds'] == 1
assert inputs['plan_sha256'] == hashlib.sha256(plan_file.read_bytes()).hexdigest()
assert inputs['spec_sha256'] == plan['spec_sha256']
assert len(plan['workloads']) == len(inputs['targets']) == 8
for index, (workload, target) in enumerate(zip(plan['workloads'], inputs['targets'])):
    assert workload['curve']['r'] == inputs['subgroup_r'] == 562427
    assert workload['target_index'] == target['target_index'] == index
    assert [target['x'], target['y']] == workload['target']
    assert tuple(workload['target']) not in old_points
    assert workload['planted'] is None
    assert target['workload_sha256'] == workload['workload_sha256']
    workload_number = int(workload['workload_id'][1:], 16)
    canonical = json.dumps({'domain': 'ecbench.algorithm_seed',
                            'parts': [spec['measurement']['seed'], workload_number, 0]},
                           sort_keys=True, separators=(',', ':'), ensure_ascii=True).encode()
    assert target['algorithm_seed'] == int.from_bytes(hashlib.sha256(canonical).digest()[:8], 'big')

print(json.dumps({'ok': True, 'targets': len(inputs['targets']),
                  'fresh_vs_previous_holdout': True,
                  'plan_sha256': inputs['plan_sha256'],
                  'input_sha256': hashlib.sha256((study / 'fixed-prefix-fresh-input.json').read_bytes()).hexdigest()},
                 indent=2))
