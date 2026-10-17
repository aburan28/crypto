#!/usr/bin/env python3
"""Check the million-assignment exact test and non-admission behavior."""

import json
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

spec = spec_from_file_location('n83_intermediate_selector_tested',
                               Path(__file__).with_name('select.py'))
select = module_from_spec(spec)
spec.loader.exec_module(select)


assert select.sign_flip_p([0.0] * 20) == 1.0
assert select.sign_flip_p([1.0] * 20) == 2 / (1 << 20)
assert select.sign_flip_p([1.0, -1.0] * 10) == 1.0
rejected = select.select({'schema': 'n83.intermediate-holdout-validated/v1',
                          'source_sha256': select.SOURCE_SHA256,
                          'targets_sha256': select.TARGETS_SHA256,
                          'evidence_status': 'VERIFICATION_FAILED'})
assert rejected['online']['selected_candidate_id'] is None
assert rejected['cold']['selected_candidate_id'] is None
print(json.dumps({'status': 'PASS', 'cases': 4}))
