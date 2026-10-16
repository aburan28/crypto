# Koblitz factor_base n=53 (2026-10-14)

Independent control-plane replay of `koblitz.factor_base.n53` on `origin/main`
(`38459ee115`) after the 2026-10-12 archive (`20261012-koblitz-factor-base-n53`,
PR #1727). Public synthetic fixtures only. Status
**`PENDING_INDEPENDENT_VALIDATION`**. No ledger promotion.

| Field | Value |
|---|---|
| Run | `20261014T014101Z-9264d0fee7` |
| Beat | `koblitz.factor_base.n53` |
| Host | macOS arm64 (`Mac.attlocal.net`), local (not hosted-isolated); high host load during measurement |
| `n` / `a` / η | 53 / 0 / 1/16 |
| Pair mode | `signed_quotient` (protocol beat; SUMMARY_ONLY) |
| Factor-base points \|F\| | 19,928 |
| Orbit columns K | 188 |
| Construction wall | 90,942.12 ms |
| Retained bytes | 89,887,680 |
| Observed peak RSS (children) | 2,881,667,072 (~2.68 GiB) |
| Resource cap | 16 GiB; under_cap = true |
| Base / fixture hash | `25354dec50ff2d6062dc7e699b1b6b3ddb0ddfffd7cbe3e01834aff162c1b0a4` |
| Whole-process wall | 4,995,823.64 ms (includes cold start; host-contended) |
| Claim-check `factor_base` | PASS |
| `verify` | PASS (13 files) |
| Status | `PENDING_INDEPENDENT_VALIDATION` |
| Ledger promotion | **none** |

Producer exit 0. Resource receipt includes `children_peak_rss_bytes` from
`getrusage`. Producer log archived as `logs/direct.stdout.jsonl`
plus empty `direct.stderr.txt`.

Same protocol knobs as 2026-10-12: η=1/16 signed-quotient SUMMARY_ONLY under
the memory gate. This is not the ledger's rank-guided η=1/10 220-orbit record
(`d8593190…`). Base size / orbit counts match the prior archive (|F|=19,928,
K=188). Construction wall is higher under host contention and is not a new
frontier claim.

## Prior archive

- `evidence/20261012-koblitz-factor-base-n53/` (PR #1727)
