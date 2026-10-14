# binary.decomposition.l8_pairs_baseline control-plane evidence (2026-10-12)

Status **`PENDING_INDEPENDENT_VALIDATION`**. Public synthetic binary Semaev
pairs-and-solve baseline only. No ledger promotion. This is **not** a
sub-2^(2ℓ) oracle win — pairs-and-solve is Θ(2^(2ℓ)).

Control-plane launch of `binary.decomposition.l8_pairs_baseline` after wiring
`semaev_decomp_bench` into `boundary_autolab.py` (producer
`examples/semaev_decomp_bench`, args `8`).

| Field | Value |
|---|---|
| Run id | `20261012T191735Z-29ca2b1d16` |
| Host | Mac.attlocal.net / macOS-26.6-arm64-arm-64bit |
| ℓ / n / m | 8 / 24 / 3 |
| Pairs wall at ℓ=8 | **0.007 s** |
| Median ms/target | **7.0** |
| Speedup vs triples (gen) | 34.7× |
| Triples gen wall at ℓ=8 | 2.191 s |
| Oracle class | `pairs_and_solve` |
| FFD on pairs arm | inapplicable (no Boolean Macaulay/GB system) |
| Whole-process wall | 31.3 s |
| Peak RSS | 56000512 bytes (~53 MiB) |
| Claim-check | `PASS` |
| Autolab status | `PENDING_INDEPENDENT_VALIDATION` |
| Verify | PASS (13 files; ledger JSON omitted from archive, pinned by `inputs/ledger_pin.json`) |

## Ladder (pairs vs triples)

| ℓ | pairs s | triples gen s | speedup vs triples |
|---:|---:|---:|---:|
| 5 | 0.000 | 0.004 | 4.1× |
| 6 | 0.000 | 0.025 | 8.2× |
| 7 | 0.002 | 0.259 | 18.0× |
| 8 | 0.007 | 2.191 | 34.7× |

Tiny ℓ pair times round to 0.000 s in the producer table (exact rows in
`claim_draft.json`).

## Ledger gate assessment

- `sub_2_to_2l_oracle=False` — this measures the pairs baseline, not a faster oracle
- `ffd_logged=False` on the pairs arm (status `inapplicable`)
- Related S4 subspace FFD at ℓ=8 remains under
  `runs_manual/binary_s4_ffd_l8_20260920` (not this beat)

## Reproduce

```bash
python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py preflight
python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py \
  launch --beat binary.decomposition.l8_pairs_baseline
python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py claim-check \
  --report research/sat_factor_base_review_20260908/autolab/runs/<run-id>/artifacts/claim_draft.json \
  --stage decomposition
```

## Files

`autolab_runs/20261012T191735Z-29ca2b1d16/` holds `state.json`, `artifacts/`,
`receipts/`, `inputs/` (protocol + ledger pin; full `boundary_targets.json`
omitted, SHA-256 in `ledger_pin.json`), and producer logs.
