# N83 (a=1) paired single-target IC/rho runs

Three paired runs recovered and independently replayed the same frozen
public point. The recorded rho/IC wall ratios are 8.2×, 5.5×, and 12.7×
(median 8.2×); IC target extraction required 8,845,441 probes in each repeat.
All three runs shared the host with unrelated verification processes, and
there is no qualifying host-isolation receipt. The retained rho interval
includes setup whose target dependence is unresolved, so these ratios are
exploratory CPU measurements.

This is the a=1 arm: `y^2 + xy = x^3 + x^2 + 1` over GF(2^83), modulus
`x^83 + x^7 + x^4 + x^2 + 1`, subgroup `r = 8569786107849059 ≈ 2^52.9`,
cofactor 1,128,547,018. The IC precompute uses parallel guided rank.

| Run | Seeds (rho / ic) | IC target interval | rho recorded interval | rho / IC recorded wall | Replay |
|---|---|---:|---:|---:|---|
| R1 | 202610061 / 21 | 7,939.1 ms | 65,042 ms | 8.2× | PASS |
| R2 | 202610062 / 22 | 7,746.0 ms | 43,003 ms | 5.5× | PASS |
| R3 | 202610063 / 23 | 5,623.3 ms | 71,324 ms | 12.7× | PASS |

Both arms verify every run; the recovered scalar is identical across arms
and runs (validation-only sidecar `3141592653589793`, never a solver input).

The workload freezes the single public point
`Q = ["2655898345537166180532960","3881526307991673265292615"]` (standalone
Python GF(2^83) scalar multiplication; no Rust producer involved in
freezing). IC receives only that point. The IC online clock starts after
target-independent preparation (base load, S3 root index build, guided
rank-600 log table) and includes target relation extraction, the
group-lift check, log recovery, and the in-process `[d]G=Q` check. The
retained rho interval includes jump setup, the walk, and in-process
validation; the receipts do not establish that all setup is target dependent.
Point construction is outside both intervals.

The n=83 subgroup (2^52.9) is smaller than n=73's (2^56.3). The retained
rho walks used 6–12M steps here versus 10–44M at n=73; IC target extraction
used 8.8M probes here versus 0.46M at n=73. Probes and walk steps have
different costs, and host contention prevents attributing the recorded wall
ratio difference quantitatively to these counters.

**Precompute (excluded, fully logged):** the parallel guided rank
(`KIC_RANK_THREADS=12`; per-column scheduling-independent scalars; solved
logs identical to the sequential path — verified byte-identical against
the frozen n=71/n=73 target rows) acquired all 600 relations in
**4.8 min wall** (mean 2.48M probes/relation, 0 failures, 0 rows
without gain) vs the 5.1 h sequential rank at n=73. Peak RSS 5.75 GB
(root table 29,878,182 entries; no pair table, no edge selectors).

## Files

- `base_n83_K600.jsonl` — retained factor base (600 orbit columns,
  99,600 points, base hash `d4a4da26d6e7cb09…`).
- `frozen/` — fixture (sidecar scalar, curve validation) and the frozen
  target point.
- `runs/<run>/` — raw producer rows, timing receipts, per-process peak
  RSS (wait4(2)), and `independent_replay.json` receipts.
- `freeze_target.py`, `run_n83_pairs.py`, `verify_n83_runs.py`,
  `assemble_claim_report.py`, `promote_to_ledger.py`.

## Reproduce

```bash
python3 freeze_target.py                    # standalone GF(2^83) freeze
python3 run_n83_pairs.py                   # 3 paired runs, ~35 min
python3 verify_n83_runs.py                  # standalone replay (PASS ×3)
python3 assemble_claim_report.py            # fail-closed claim
python3 promote_to_ledger.py                # ledger promotion
```

`derived_exploratory_claim_report.json` links the raw receipts and keeps
the controlled online speedup unknown pending a qualifying host-isolation
receipt, matched resource envelope, canonical candidate and workload
manifests, and a frozen rho online boundary.
