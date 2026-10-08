# N85 (a=0) single-target IC/rho boundary run

Seventh primary single-target online rung — the **a=0 arm past n=73**
(`y^2 + xy = x^3 + 1` over GF(2^85), modulus `x^85 + x^8 + x^2 + x + 1`,
subgroup `r = 14351744671810121 ≈ 2^53.7`, cofactor 2,695,534,732), using
the parallel guided rank precompute. One compact-orbit index-calculus
solve and one Pollard-rho solve per paired run on the same frozen public
point; a toy boundary rung, not ECC2K-130 evidence or a key-recovery
claim.

| Run | Seeds (rho / ic) | IC online | rho online | rho / IC | Replay |
|---|---|---:|---:|---:|---|
| R1 | 202610071 / 31 | 499.3 ms | 113,566 ms | 227.4× | PASS |
| R2 | 202610072 / 32 | 400.5 ms | 277,792 ms | 693.4× | PASS |
| R3 | 202610073 / 33 | 215.5 ms | 239,175 ms | 1,109.7× | PASS |

The workload freezes the single public point
`Q = ["12599156368032549807479514","36570408432718250497127020"]` (standalone
Python GF(2^85) scalar multiplication; no Rust producer involved in
freezing). IC receives only that point; the known-answer scalar
`9876543210987654` lives in `frozen/fixture.json` as a validation-only
sidecar. The IC online clock starts after target-independent preparation
(base load, S3 root index build, guided rank-600 log table) and includes
target relation extraction, the group-lift check, log recovery, and the
in-process `[d]G=Q` check; rho's online clock includes per-target jump
setup, the walk, and in-process validation. Point construction is
excluded from both arms. The target relation is deterministic (same
probe count and decomposition on every run).

**Precompute (excluded, fully logged):** parallel guided rank
(`KIC_RANK_THREADS=12`), 600 relations, mean ~4.2M probes/relation,
~9 min wall, 0 failures, 0 rows without gain; peak RSS ~4.35 GB (root
table 30,598,167 entries over 30.6M states; no pair table, no edge
selectors).

Operational note: the original R2 launch was interrupted by an
unrelated process-group kill after its rho arm completed; the recovery
driver (`recover_n85.py`) reran only the missing IC arm and launched R3.
All receipts (wait4(2) peak RSS, stdout hashes, heartbeats) are
retained per arm.

## Files

- `base_n85_K600.jsonl` — retained factor base (600 orbit columns,
  102,000 points, base hash `e1b8e4a61202962d…`).
- `frozen/` — fixture (sidecar scalar, curve validation) and the frozen
  target point.
- `runs/<run>/` — raw producer rows, timing receipts, and
  `independent_replay.json` receipts.
- `freeze_target.py`, `run_n85_pairs.py`, `recover_n85.py`,
  `verify_n85_runs.py`, `assemble_claim_report.py`.

## Reproduce

```bash
python3 freeze_target.py
python3 run_n85_pairs.py                   # 3 paired runs
python3 verify_n85_runs.py                  # standalone replay (PASS ×3)
python3 assemble_claim_report.py            # fail-closed claim
```

Claim boundary: public synthetic Koblitz n=85 a=0, one unseen online
target, identical frozen public point in both arms, constant-factor win
only. Not ECC2K-130 evidence (the challenge is n=131, same curve
family — see `RESEARCH_ECC2K130_IC_FEASIBILITY.md`), not an asymptotic
sub-rho claim, not key recovery, no deployed-curve impact.
