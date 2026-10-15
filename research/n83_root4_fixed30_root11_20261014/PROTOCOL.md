# N83 four-byte fixed30 slots with 11-byte states

The four-byte exact-root slot encodes a 29-bit state index, one root selector
bit, and two fingerprint bits. When `2*S(K) < 2^30`, the `2^30`-slot table can
hold every canonical root of the symmetric S3 state set. This experiment
allows the existing 11-byte first-root/reconstruct-second-root state layout
under `KIC_INDEXED4_FIXED30=1`, retaining the exact-key check and first
witness. It changes the early table request from `2^30` to `2^29` entries,
which the constructor doubles to `2^30` four-byte slots. The mode keeps the
29-bit state-ID and complete-root capacity guards.

The corrected seven-arm execution uses IC Rust source SHA-256
`3112f5617548d68dab023d8dd9133cd68a1d07fd89e05311c3a2f4b532990c74`.
The first execution used source SHA-256
`5281bfc1c370255d8c03851fe85b4b35499ee55c066ba460b5de414b59d202a3`.
Its direct 21-byte control has the indexed-insertion error recorded in the
preceding fixed30 campaign; retain that execution's raw arm statuses and
verification receipt separately from the corrected run.
The same 1,000 public N83 points frozen for the preceding factor-base
refinement are used in all arms. Their SHA-256 is
`35e34151dcc1ad6048939560adf28fc7c1b1017833fe005a656c3f4e2bc07de8`;
the independently checked scalar sidecar SHA-256 is
`f8cbc90f7f592aae5014ca4f2693a3327a72e935970f22e8f02868e4aeed7f71`.
Producers receive only public points. These previously frozen points make
this a secondary layout and factor-base selection comparison.

| Order | Arm | Usable points | Index slots | State bytes |
| ---: | --- | ---: | ---: | ---: |
| 1 | K3300 root4f30r11 | 547,800 | 2^30 × 4 | 11 |
| 2 | K3300 root11 | 547,800 | 2^30 × 5 | 11 |
| 3 | K3300 root4f30 | 547,800 | 2^30 × 4 | 21 |
| 4 | K3350 root4f30r11 | 556,100 | 2^30 × 4 | 11 |
| 5 | K3400 root4f30r11 | 564,400 | 2^30 × 4 | 11 |
| 6 | K3450 root4f30r11 | 572,700 | 2^30 × 4 | 11 |
| 7 | K3500 root4f30r11 | 581,000 | 2^30 × 4 | 11 |

Each arm starts a fresh process with rank seed 29, four rank threads, one
target thread, a 1,200-second child cap, and a 14-GiB observed peak-RSS
guard. The order above is fixed before measurement. K3300 controls isolate
the table width and state reconstruction on the identical base and target
points. The verifier replays all full factor-base points, 1,000 group
relations and scalars per completed arm, the target-online phase sums, and
same-K base hash, key count and probe parity. Failed and over-cap arms stay
in the result. Among verified arms within the RSS guard, the descriptive
selection rule is the lowest arithmetic mean of 1,000 target-online
intervals; paired target-block differences and bootstrap intervals remain
alongside it. This secondary panel has no rho arm.

Five further public points were frozen before this new mode was run at
`research/n83_root4_fixed30_20261014/reserve_targets.jsonl`, SHA-256
`fef03b92aef0352e2e21ff7362878b8bed8f4826fe8d90503d70174508d631f1`.
Their answer sidecar SHA-256 is
`a1979e2bbcd3d82d47e29f38d66056228e3733faf125b0d28f8e4c32831b5665`.
A changed selected candidate uses these five unseen points for alternating
same-point Tier 1 IC/rho validation with complete cost and noise receipts.
