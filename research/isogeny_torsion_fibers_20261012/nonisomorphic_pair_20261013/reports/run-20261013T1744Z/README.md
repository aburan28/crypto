# V2 run: clean exit and exact replay

Start: 2026-10-13 17:44 UTC. Host: Darwin 25.6.0 arm64, shared, L0.
Dependencies: [protocol](../../PROTOCOL.md) committed at `70b119bfab` and
[wrapper remediation](../../REMEDIATION_V2.md) committed at `fa1a6d4dee`.
The release source SHA-256 is
`72b25721734346d82e4268c54622251312b050aa2f23e573ec52616e165d60b2`;
binary SHA-256 is
`cc28892bcaa8830ddc6d1c643088c715c7ad4e0c6d097a5c40d618952c277cb6`.

Run and replay from repository root:

```sh
target/release/examples/isogeny_fiber_nonisomorphic_20261013 > research/isogeny_torsion_fibers_20261012/nonisomorphic_pair_20261013/raw-v2.tsv 2> research/isogeny_torsion_fibers_20261012/nonisomorphic_pair_20261013/run-v2.log
target/release/examples/isogeny_fiber_nonisomorphic_20261013 > research/isogeny_torsion_fibers_20261012/nonisomorphic_pair_20261013/replay-v2.tsv 2> research/isogeny_torsion_fibers_20261012/nonisomorphic_pair_20261013/replay-v2.log
```

Both executions returned **exit 0**. The two TSVs and two logs compare byte
for byte; the v1 raw TSV also equals the v2 raw TSV. The final
[raw records](../../raw-v2.tsv) contain 24/24 verified arm rows for eight
frozen source targets. The [log](../../run-v2.log) contains the map and target
preflight. The [summary](../../summary-v2.tsv) and
[phase table](../../phase_summary-v2.tsv) are derived by the native
`summarize.rs`. Totals: rho 664, source relations 3,312, destination
relations 3,794 counted group operations. Source/destination collection
trials totaled 231/290. The [independent checker](../../independent_v2.log)
verified every scalar, reconstructed all 16 relation solves, checked the
source and destination models and eight rational target fibers, and
reconciled phase sums.

![V2 counted group operations](group_ops-v2.svg)

The [run PDF](report.pdf) embeds this graph. The result
[report](../../REPORT.md) states the complete interpretation and limits.
