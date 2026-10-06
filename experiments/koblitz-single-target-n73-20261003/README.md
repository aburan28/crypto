# N73 single-target IC/rho boundary run

Fifth primary single-target online rung. One compact-orbit index-calculus
solve and one Pollard-rho solve per paired run on the same frozen public
point. The point is a synthetic target on the Koblitz curve over
\(\mathbb{F}_{2^{73}}\) (a=0, b=1, modulus \(x^{73}+x^4+x^3+x^2+1\),
subgroup \(r = 86020738150056119 \approx 2^{56.3}\), cofactor 109,796).
This is a toy boundary rung, not ECC2K-130 evidence or a key-recovery
claim.

| Run | Targets | IC online | rho online | rho / IC | Replay |
|---|---:|---:|---:|---:|---|
| R1 (seed 202610041 / 7) | 1 | 249.024 ms | 46,034.255 ms | 184.9× | PASS |
| R2 (seed 202610042 / 8) | 1 | 326.159 ms | 683,096.952 ms | 2,094.7× | PASS |
| R3 (seed 202610043 / 9) | 1 | 543.772 ms | 646,986 ms | 1,189.7× | PASS |

**Median online speedup 1,189.7× (range 184.9×–2,094.7×).**  Both arms
verify every run; the recovered scalar is identical across arms and runs
(validation-only sidecar, never a solver input).  The target relation is
deterministic: every run reports the same 457,561-probe decomposition.

The workload freezes the single public point
`Q = ["1852858978593996699777","7821257062942871962177"]` in
polynomial-basis decimal encoding.  IC receives only that point.  The
known-answer scalar `12345678901234567` lives in
`frozen/fixture.json` as a validation-only sidecar.

The IC online clock starts after target-independent preparation (base
load, S3 root index build, guided rank-600 log table) and includes
target relation extraction, the group-lift check, log recovery, and the
in-process `[d]G=Q` check.  rho's online clock includes per-target jump
setup, the walk, and in-process validation; point construction is
excluded from both arms.

The precompute is fully charged and logged, excluded from the online
claim per the frozen contract: the guided rank acquired 600 relations
(mean 57–70M probes per relation, heavy-tailed; R1 5.1 h single-thread;
R2/R3 ran while an unrelated verification process shared the host, so
their walls are conservative).  The parallel guided rank
(`KIC_RANK_THREADS`, added 2026-10-05) reproduces the frozen target row
exactly and cuts this stage 11.5× at 12 threads uncontended — see
`research/sat_factor_base_review_20260908/autolab/runs_manual/koblitz_parallel_rank_n73_20261005/`.

## Files

- `claim_report_vs_rho.json` — assembled claim (fail-closed; requires
  three PRODUCERS_COMPLETE runs with PASSing replays).
- `single-target-results.csv` — per-run measurements.
- `frozen/` — fixture, frozen target, curve validation.
- `runs/<run>/` — raw producer rows, timing receipts, per-process peak
  RSS (wait4(2)), and `independent_replay.json` receipts.
- `run_single_target.py` (R1), `run_single_target_r2_r3.py` (R2/R3) —
  frozen launch drivers.
- `verify_single_target.py`, `verify_single_target_r2_r3.py` —
  standalone standard-library GF(2^73) replays.
- `assemble_claim_report.py`, `promote_to_ledger.py` — claim assembly
  and ledger promotion.

## Reproduce

```bash
python3 run_single_target_r2_r3.py R2          # rho then IC, ~6 h
python3 verify_single_target_r2_r3.py R2        # standalone replay
python3 assemble_claim_report.py                 # fail-closed claim
python3 promote_to_ledger.py                    # ledger promotion
```

Claim boundary: public synthetic Koblitz n=73, one unseen online
target, identical frozen public point in both arms, constant-factor win
only.  Not ECC2K-130 evidence (the challenge is n=131, same curve
family — see `RESEARCH_ECC2K130_IC_FEASIBILITY.md`), not an asymptotic
sub-rho claim, not key recovery, no deployed-curve impact.
