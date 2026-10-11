# n=71 single-target online `vs_rho` record — ledger evidence (2026-10-02)

The `u128` compact-orbit implementation recovered the same frozen n=71 public
point in three paired runs. The median recorded IC and automorphism-discounted
Pollard-rho online wall ratio is 2,864.1x; independent replay verifies both
scalars in every pair.

## Result

One frozen public target `Q` (see `ledger-freeze/`), three paired observations
(`runs/ledger-R1`, `runs/ledger-R2`, `runs/ledger-R3`), one process per arm,
sequential, same 16 GiB envelope:

| Run | IC online (ms) | rho online (ms) | rho / IC |
|---|---:|---:|---:|
| ledger-R1 | 14.457 | 70,747.3 | 4,893.5x |
| ledger-R2 (median) | 10.208 | 29,237.4 | 2,864.1x |
| ledger-R3 | 11.874 | 18,200.9 | 1,532.8x |

Median paired wall ratio **2,864.1x** (range 1,532.8x–4,893.5x). The IC relation is
deterministic: the same 14,554-probe relation on all three runs. Rank 600/600
with 0 failures and 0 rows without gain on every run (guided policy).

Curve: `K_0` over `GF(2^71)` (`x^71+x^5+x^3+x+1`), subgroup
`r = 5513228015079457` (~2^52.3), cofactor 428,276. Base: 600
signed-Frobenius orbit columns, 85,200 points (hash `f3acde0c…`), root index
25,558,188 entries, 25,560,000 regular states. **Zero pair-table entries,
zero edge selectors at every stage.** Peak RSS ~5.6 GB (IC) / ~2.1 GB (rho).
Automorphism discount √(2n), A=142 signed Frobenius.

Online clocks start after reusable base/S3-root-index/guided rank-600 log
preparation (~39–48 s excluded setup) and at rho's first target-dependent
walk step. Both arms recover the identical scalar on every run and verify
`[d]G = Q`. Standalone Python GF(2^71) replay passes 15/15 checks on all
three runs (`independent-replay-ledger.json`, via `validate_ledger_runs.py`).
A current autolab `claim-check --stage vs_rho` returns **FAIL** for the retained
`claim_report_vs_rho.json`: canonical identity and replay-digest fields, exact
five-phase IC timing, the complete resource envelope, and rho policy remain
unrecorded there. The raw pairs and independent replay remain the result
receipts.

## Evidence scope

The field is GF(2^71) and the subgroup has order 5,513,228,015,079,457.
The recorded wall ratios are exploratory because these runs have no auditable
host-level CPU isolation receipt. IC probes and rho walk steps use different
native units; calibrated total work `S` and controlled online speedup remain
unknown. The three repetitions time one point, so they measure timing
variation for that point.

## Provenance notes

- The `u128` field extension (library twins `Gf2_128`/`FastBinaryCurve128`,
  wide sparse search, producer wide paths) keeps every `n ≤ 63` code path
  untouched: the n=61 R4 reproduction and the n=61 rho command rerun byte-identical
  after the change (same base hash, scalar, relation, walk).
- The decimal-string wire format for wide field values fixes a serde `u128`
  range panic; legacy `u64` fixtures are unchanged.
- `runs/pilot-R1-concurrent-session/` preserves a concurrent session's pilot
  attempt (its rho observation with the wide backend is valid; its IC arm ran
  against a pre-fix binary and emitted no rows). The top-level `README.md`
  distinguishes this ledger series from the earlier 16.78x pilot pair;
  `candidate-*.json` and `single-target-results.csv` remain pilot artifacts.
- `base_n71_K600.jsonl` is the retained base (deterministic x-scan;
  identical hash across independent constructions).
- `dump_out.jsonl`, `smoke_out.jsonl` are construction smoke receipts
  (scalar 314159265358979, since solved; not evidence targets).

## Files

- `claim_report_vs_rho.json` — retained claim report and native counter ledger
  (R2 median row; current claim-check status FAIL).
- `ledger-freeze/` — frozen target, sidecar scalar, freeze receipt.
- `runs/ledger-R{1,2,3}/` — paired raw rows (`ic.jsonl`, `rho.jsonl`) plus
  `/usr/bin/time -l` receipts (`ic_time.txt`, `rho_time.txt`).
- `validate_ledger_runs.py`, `independent-replay-ledger.json` — replay.
