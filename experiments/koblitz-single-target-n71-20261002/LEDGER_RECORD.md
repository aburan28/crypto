# n=71 single-target online `vs_rho` record — ledger evidence (2026-10-02)

Fourth rung of the primary single-target ladder, and the first past the old
u64-packing ceiling: compact-orbit index calculus on `u128` field words beats
automorphism-discounted Pollard rho on the identical frozen public point.

## Result

One frozen public target `Q` (see `ledger-freeze/`), three paired observations
(`runs/ledger-R1`, `runs/ledger-R2`, `runs/ledger-R3`), one process per arm,
sequential, same 16 GiB envelope:

| Run | IC online (ms) | rho online (ms) | rho / IC |
|---|---:|---:|---:|
| ledger-R1 | 14.457 | 70,747.3 | 4,893.5x |
| ledger-R2 (median) | 10.208 | 29,237.4 | 2,864.1x |
| ledger-R3 | 11.874 | 18,200.9 | 1,532.8x |

Median paired ratio **2,864.1x** (range 1,532.8x–4,893.5x). The IC relation is
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
A schema-complete claim report (`claim_report_vs_rho.json`) passes the
autolab `claim-check` for stage `vs_rho`.

## Explicit non-claims

Public synthetic fixture, 52-bit subgroup; constant-factor win only — no
asymptotic sub-rho claim (total operation-count boundary not comparable, S
unknown); not ECC2K-130 evidence; no key recovery, external targets, or
deployed-curve security impact. Multi-target amortized results stay secondary.

## Provenance notes

- The `u128` field extension (library twins `Gf2_128`/`FastBinaryCurve128`,
  wide sparse search, producer wide paths) keeps every `n ≤ 63` code path
  untouched: the n=61 R4 reproduction and the n=61 rho command rerun byte-identical
  after the change (same base hash, scalar, relation, walk).
- The decimal-string wire format for wide field values fixes a serde `u128`
  range panic; legacy `u64` fixtures are unchanged.
- `runs/pilot-R1-concurrent-session/` preserves a concurrent session's pilot
  attempt (its rho observation with the wide backend is valid; its IC arm ran
  against a pre-fix binary and emitted no rows). The top-level `README.md`,
  `candidate-*.json`, and `single-target-results.csv` are that session's
  pilot documents (16.78x observation with Sage replay); the ledger record
  above is independent of them.
- `base_n71_K600.jsonl` is the retained base (deterministic x-scan;
  identical hash across independent constructions).
- `dump_out.jsonl`, `smoke_out.jsonl` are construction smoke receipts
  (scalar 314159265358979, since solved; not evidence targets).

## Files

- `claim_report_vs_rho.json` — schema-complete claim (R2 median primary).
- `ledger-freeze/` — frozen target, sidecar scalar, freeze receipt.
- `runs/ledger-R{1,2,3}/` — paired raw rows (`ic.jsonl`, `rho.jsonl`) plus
  `/usr/bin/time -l` receipts (`ic_time.txt`, `rho_time.txt`).
- `validate_ledger_runs.py`, `independent-replay-ledger.json` — replay.
