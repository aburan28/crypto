# Koblitz decomposition n=31 dim-16 m=2 (2026-10-14)

Control-plane launch of `koblitz.decomposition.n31_m2_dim16` on `origin/main`
(`6ea104b4`). Public synthetic fixtures only. Status
**`PENDING_INDEPENDENT_VALIDATION`**. No ledger promotion.

| Field | Value |
|---|---|
| Run | `20261014T043448Z-f2b0ba0fd2` |
| Beat | `koblitz.decomposition.n31_m2_dim16` |
| Host | macOS arm64 (`Mac.attlocal.net`), local |
| Producer | `symmetrised_oracle_bench --only 0 31 2 --targets 1 --no-sat --node-budget 50000` |
| `n` / `a` / `m` / dim | 31 / 0 / 2 / 16 |
| Unknowns | 32 (x-chained); 31 (sym) |
| System degree | quadratic |
| Oracle class | groebner (F4) |
| FFD | min=3 max=4 mean=3.5 (x-chained FFD=3, sym FFD=4; 1 draw) |
| Median ms/target | 13915.0 (x-chained found); sym refuted median 7826.0 |
| Disagreements | 0 (gates ok) |
| Verdict mix | x-chained: found=1; sym: refuted=1; inconclusive=0 |
| Whole-process wall | 158918.43 ms (timed sample); cold-start wall 2376343.37 ms |
| Peak RSS (children) | 3,738,796,032 (~3.48 GiB); cap 16 GiB; under_cap |
| Claim-check `decomposition` | PASS |
| `verify` | PASS (13 files) |
| Status | `PENDING_INDEPENDENT_VALIDATION` |
| Ledger promotion | **none** |

Producer exit 0. Empty-class medians in `claim_draft.json` are IEEE NaN
literals (Python `json` allow_nan); they correspond to unused found/refuted
arms on this single-target draw.

## Prior related archives

- `runs_manual/n31_decomp_m2_dim16_20260911/` (manual drafts)
- `autolab_groebner_hyperopt_20260921/runs/selected_default_n31_*` (policy selection)
