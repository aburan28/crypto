# K3300 root4 versus root11 on 4,000 frozen N83 points

The selected K3300 factor base has 547,800 distinct usable subgroup points,
3,300 folded columns, and complete base digest
`a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b`.
This experiment keeps that exact base and compares its four-byte `2^31`
root table with reconstructed 11-byte states against the selected five-byte
`2^30` root table with the same states. The latter's sealed candidate is
`IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0h0f8ab8bfcff7`.
The four-byte candidate receives its IC1 identity after its full base and
algorithm manifest are replayed.

Both arms use the public 4,000-point confirmation fixture frozen in
`research/n83_factor_base_confirm_20261015/targets.jsonl`, SHA-256
`451df8524b2c356dd4a523260620c4dff406975a122762b98ddddb7987b75c2c`.
Its independent scalar sidecar SHA-256 is
`8e18c9c97ba9f57e5eb43f1b39b9e9575b5cfab87840c25c428dac150c868101`.
The points are disjoint from earlier N83 selection and fresh validation
fixtures. The source has SHA-256
`5e9fd826d22816b429a82e658c517d9387916aabe7660f708e4c29901cce647e`.

| Order | Layout | Exact slots | Retained state bytes | Child cap |
| ---: | --- | ---: | ---: | ---: |
| 1 | K3300/root4 | `2^31 × 4` | 11 | 1,200 s |
| 2 | K3300/root11 | `2^30 × 5` | 11 | 1,200 s |

Each arm uses a fresh process with rank seed 29, four rank threads, one
target thread, an empty factor-log/index setup, and a 14-GiB observed peak
RSS guard. The verifier checks the complete base and subgroup orbit, 4,000
four-point relations and recovered scalars, exclusive target-online phase
sums, both base digests, and equality of exact root-key counts. The
predeclared contrast is
`Δ = mean_i(online_ms(root11, i) − online_ms(root4, i))`.
The 95% interval resamples 80 contiguous 50-target block means 10,000
times with seed `2026101599`. A positive lower endpoint confirms root4's
online advantage under this allocation; a negative upper endpoint confirms
root11's. A crossing interval is inconclusive. The allocation, setup,
matrix solve, memory, process wall, failures and timeouts remain in the raw
receipt. This secondary layout comparison uses no rho arm; the separate
fresh one-target Tier 1 IC/rho pairs are in PR #1855.
