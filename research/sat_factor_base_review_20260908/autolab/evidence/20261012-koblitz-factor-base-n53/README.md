# Koblitz factor_base n=53 (2026-10-12)

Fresh control-plane launch of `koblitz.factor_base.n53` on `origin/main`
(`c62555b7a`) after the unmerged 2026-10-10 draft in closed PR #1621. Public
synthetic fixtures only. Status **`PENDING_INDEPENDENT_VALIDATION`**. No ledger
promotion.

| Field | Value |
|---|---|
| Run | `20261012T161605Z-ac2e7ab560` |
| Beat | `koblitz.factor_base.n53` |
| Host | macOS arm64 (`Mac.attlocal.net`), local (not hosted-isolated) |
| `n` / `a` / η | 53 / 0 / 1/16 |
| Pair mode | `signed_quotient` (protocol beat; SUMMARY_ONLY) |
| Factor-base points \|F\| | 19,928 |
| Orbit columns K | 188 |
| Construction wall | 6,231.55 ms |
| Retained bytes | 89,887,680 |
| Support table allocated | 92,274,688 |
| Observed peak RSS (children) | 3,189,948,416 (~2.97 GiB) |
| Resource cap | 16 GiB; under_cap = true |
| Base hash | `e2c0276a010a4787cdd617043ad9833a528c919e9e3882b7de377a0e9e4f7015` |
| Whole-process wall | 608,525.14 ms (includes cold start / process setup) |
| Claim-check `factor_base` | PASS |
| `verify` | PASS (13 files) |
| Status | `PENDING_INDEPENDENT_VALIDATION` |
| Ledger promotion | **none** |

Producer exit 0. Resource receipt includes `children_peak_rss_bytes` from
`getrusage`. Small producer log archived as `logs/direct.stdout.jsonl`
(41 KiB) plus empty `direct.stderr.txt`.

This beat is the protocol `eta=1/16` signed-quotient SUMMARY_ONLY construction
under the memory gate. It is not the ledger's rank-guided η=1/10 220-orbit
record (`d8593190…`), which remains the published factor_base frontier.

## Prior closed draft

Closed PR #1621 archived a mid-receipt-crash recovery of the same beat as
`20261010T044040Z-2b5ad63dff` and never merged. This run replaces that draft
with a clean receipt on current main.

## Later independent replay

Control-plane replay archived at `../20261014-koblitz-factor-base-n53/` (run `20261014T014101Z-9264d0fee7`).
