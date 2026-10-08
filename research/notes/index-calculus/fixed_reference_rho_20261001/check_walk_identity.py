#!/usr/bin/env python3
"""Check that `koblitz_rho_fixture ... strong` (library koblitz_strong_rho) walks the same
trajectory as `koblitz_rho_batch_ks_strong` rung 3 (32 lanes, dp_bits 4) for the same seeds.

Run from the repository root after
`cargo build --release --example koblitz_rho_fixture --example koblitz_rho_batch_ks_strong`.
Writes walk_identity.json next to this script and exits non-zero on any mismatch.
"""
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.join("target", "release", "examples")
CASES = [(n, seed, None) for n in (37, 41, 53) for seed in (531310, 1, 2, 3)]
CASES.append((53, 531310, 476811900269))  # the selected-panel target (PR #1090)


def rows(cmd, env=None):
    full_env = dict(os.environ)
    full_env.update(env or {})
    out = subprocess.run(cmd, capture_output=True, text=True, env=full_env, check=True).stdout
    return [json.loads(line) for line in out.splitlines() if line.startswith("{")]


results = []
for n, seed, explicit in CASES:
    env = {"KIC_RHO_RUNG": "3", "KIC_RHO_LANES": "32", "KIC_RHO_DP_BITS": "4"}
    if explicit is not None:
        env["KIC_RHO_EXPLICIT_SCALAR"] = str(explicit)
    example = [r for r in rows([os.path.join(EX, "koblitz_rho_batch_ks_strong"), str(n), "0",
                                "signed_frobenius", "1", str(seed)], env)
               if r.get("kind") == "rho_ks_batch_fixture"][0]
    cmd = [os.path.join(EX, "koblitz_rho_fixture"), str(n), "0", "signed_frobenius", "1", "strong", str(seed)]
    if explicit is not None:
        cmd.append(str(explicit))
    fixture = rows(cmd)[0]
    same = all(example[k] == fixture[k] for k in ("walk_steps", "recovered_fixture_scalar", "published_fixture_scalar"))
    results.append({"n": n, "batch_seed": seed, "explicit_scalar": explicit,
                    "example_walk_steps": example["walk_steps"], "fixture_walk_steps": fixture["walk_steps"],
                    "planted": fixture["published_fixture_scalar"], "recovered": fixture["recovered_fixture_scalar"],
                    "identical": same})
    print(results[-1], flush=True)
json.dump({"cases": results, "all_identical": all(r["identical"] for r in results)},
          open(os.path.join(HERE, "walk_identity.json"), "w"), indent=2)
sys.exit(0 if all(r["identical"] for r in results) else 1)
