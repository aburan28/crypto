"""A deterministic stand-in producer for the command adapter's tests.

It burns a fixed amount of CPU, then writes an ICMS metrics file whose counts
depend only on the seed and the declared base size, the way a real producer's
operation counts depend only on its inputs.
"""
import json
import os
import sys

seed, size = int(sys.argv[1]), int(sys.argv[2])
acc = seed
for i in range(20000 * size):
    acc = (acc * 6364136223846793005 + i) % (1 << 64)
ops = 1000 * size + seed % 7
with open(os.environ["ICMS_METRICS_PATH"], "w") as fh:
    json.dump({
        "outcome": {"status": "complete", "verified": True},
        "units": {"count.group_additions": {"total": ops, "deterministic": True, "host_dependent_because": []}},
        "metrics": {"factor_base": {"usable_points": 2 * size, "columns": size, "producer_specific": acc % 97},
                    "system": {"n_vars": size, "cnf_clauses": 3 * size, "conflicts": 5}},
        "phases": {"factor_base": {"ops": size, "ops_unit": "count.group_additions", "wall_ns": None},
                   "target_pdp": {"ops": ops - size, "ops_unit": "count.group_additions", "wall_ns": None}},
        "windows": {"cold_end_to_end": {"conformance": "exact_operations", "ops": ops, "ops_unit": "count.group_additions",
                                        "wall_ns": None}},
    }, fh)
print(json.dumps({"ops": ops}))
