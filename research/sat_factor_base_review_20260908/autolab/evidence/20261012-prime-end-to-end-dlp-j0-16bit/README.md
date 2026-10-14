# prime.end_to_end_dlp.j0_16bit control-plane evidence (2026-10-12)

Status **`PENDING_INDEPENDENT_VALIDATION`**. Public synthetic j=0 known-answer
IC vs Pollard ρ only. No ledger promotion. Not a vs_ρ crossover (IC wall ≫ ρ).

Control-plane launch of `prime.end_to_end_dlp.j0_16bit` after wiring
`j0_known_answer_bits` into `boundary_autolab.py` (producer
`examples/j0_known_answer_bits`, args `{bits}` → `16`). Protocol-registered
rung is **16 bits** (agent_priorities prose may mention ≥20; this beat follows
`protocol.json`).

| Field | Value |
|---|---|
| Run id | `20261012T200056Z-0fa4798ad8` |
| Host | Mac.attlocal.net / macOS-26.6-arm64-arm-64bit |
| Bits / curve | 16 / `j0-bench-16bit` (n=65521) |
| IC wall | **2623.62 ms** |
| ρ wall | **4.13 ms** |
| IC recovers d / [d]G=Q | true |
| IC matches truth | true |
| IC agrees ρ | true |
| Target orbits / trials / ρ steps | 40 / 200000 / 65536 |
| Whole-process wall (timed) | 1.807 s (3 repeats; stdout_stable=false) |
| Peak RSS | 55508992 bytes (~53 MiB) |
| Claim-check | `PASS` (`end_to_end_dlp`) |
| Autolab status | `PENDING_INDEPENDENT_VALIDATION` |
| Verify | PASS (13 files; ledger JSON omitted from archive, pinned by `inputs/ledger_pin.json`) |

## Claim boundary

- Public synthetic known-answer only
- Not key recovery on deployed curves
- Not asymptotic sub-ρ
- Not ledger promotion until independent validation
- IC slower than ρ on this instance (expected)

## Reproduce

```bash
python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py preflight
python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py \
  launch --beat prime.end_to_end_dlp.j0_16bit
python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py claim-check \
  --report research/sat_factor_base_review_20260908/autolab/runs/<run-id>/artifacts/claim_draft.json \
  --stage end_to_end_dlp
```

## Files

`autolab_runs/20261012T200056Z-0fa4798ad8/` holds `state.json`, `artifacts/`,
`inputs/` (protocol + ledger pin), `logs/`, and `receipts/`.
