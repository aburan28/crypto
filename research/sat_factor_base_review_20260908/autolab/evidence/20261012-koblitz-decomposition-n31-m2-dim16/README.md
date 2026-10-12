# koblitz.decomposition.n31_m2_dim16 control-plane evidence (2026-10-12)

Status **`PENDING_INDEPENDENT_VALIDATION`**. Public synthetic Semaev
decomposition oracle panel only. No ledger promotion.

Control-plane launch of
`koblitz.decomposition.n31_m2_dim16` with `--fixtures 8` after wiring the
decomposition stage into `boundary_autolab.py` (producer
`examples/symmetrised_oracle_bench`, args
`--only 0 31 2 --targets 8 --no-sat --node-budget 50000`).

| Field | Value |
|---|---|
| Run id | `20261012T173525Z-27b5dd9f1d` |
| Host | Mac.attlocal.net / macOS-26.6-arm64-arm-64bit |
| Targets | 8 |
| Dim / m / unknowns | 16 / 2 / 32 |
| Median ms/target (honesty max class) | **12212.0** |
| FFD min/mean/max | 3 / 3.5 / 4 |
| Disagreements | 0 |
| Producer wall (paired_bench) | 175.1 s |
| Timed whole-process wall | 175.1 s |
| Peak RSS | 5235032064 bytes (4.88 GiB) |
| Under 16 GiB cap | True |
| Claim-check | `PASS` |
| Autolab status | `PENDING_INDEPENDENT_VALIDATION` |
| Verify | PASS (13 files; ledger JSON omitted from archive, pinned by `inputs/ledger_pin.json`) |

## Arms

| Arm | found | refuted | inconc | found ms med | refuted ms med | ffd | gate |
|---|---:|---:|---:|---:|---:|---:|---|
| x-chained F4 | 5 | 3 | 0 | 3274.0 | 12212.0 | 3 | ok |
| sym F4 | 4 | 4 | 0 | 1109.0 | 6114.0 | 4 | ok |

Acceptance vs ledger gates:
`dim_16_m_2_quadratic=True`,
`disagreements_zero=True`,
`ffd_logged=True`,
`median_under_1h=True`.

## Reproduce

```bash
python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py preflight
python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py \
  launch --beat koblitz.decomposition.n31_m2_dim16 --fixtures 8
python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py claim-check \
  --report research/sat_factor_base_review_20260908/autolab/runs/<run-id>/artifacts/claim_draft.json \
  --stage decomposition
```

## Files

`autolab_runs/20261012T173525Z-27b5dd9f1d/` holds `state.json`, `artifacts/`,
`receipts/`, `inputs/` (protocol + ledger pin; full `boundary_targets.json`
omitted, SHA-256 in `ledger_pin.json`), and producer logs.
