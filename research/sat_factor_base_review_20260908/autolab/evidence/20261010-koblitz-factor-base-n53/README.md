# Koblitz factor_base n=53 (2026-10-10)

Control-plane launch of `koblitz.factor_base.n53` on the IC accounting-repair tip
(with merge-corrupted `koblitz_rank_fixture.rs` restored from `fe07f6088^1` and
`run_timed` emitting `children_peak_rss_bytes` again).

| Field | Value |
|---|---|
| Run | `20261010T044040Z-2b5ad63dff` |
| Beat | `koblitz.factor_base.n53` |
| Host | macOS arm64 (local), not hosted-isolated |
| Claim-check `factor_base` | PASS |
| Status | `PENDING_INDEPENDENT_VALIDATION` |
| Ledger promotion | **none** |

Whole-process wall was recovered from the producer `charged_total_ms` after the
autolab crashed mid-receipt write (pre-fix). Re-run through a fixed
`boundary_autolab.py` for a full resource receipt if an independent host needs
the children peak RSS from `getrusage`.

Does not touch the n=61 L=16384/65536 K=700 panels from PR #1419.
