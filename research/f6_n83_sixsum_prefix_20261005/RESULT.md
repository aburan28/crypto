# One fixed n83 source summand does not unlock the root reduction

The [preregistered protocol](PROTOCOL.md) froze 16 source-point indices,
the public T001 target, and a structural success rule before this
run. The native process completed all 16 rows under the monitored
120-second and 7-GiB RSS envelope. Its [raw JSONL](panel_1.jsonl)
records every exact source abscissa, code, rank, affine row, operation
count, and resource observation; [status](panel_1.status) and
[stderr](panel_1.stderr.txt) are retained with [SHA-256
checksums](SHA256SUMS). This is a finite deterministic
panel, not a natural-yield estimate or a one-target IC solve.

The existing direct system on K0 curve
`icv1-f2m83-tm6151469093347-debefd74` has 428 variables and 415
cubic equations. The dimension-16 standard source subspace has 64,907
geometric points. In the setup control, the six planted source points
from indices `[0,2,4,6,8,10]` satisfied every original equation and
every equation after the first source code was substituted. The
ordinary input was the fixed public subgroup point T001,
`(355fb5df7a905f16921eb,5900a390f42d290f1bbe)`, lifted to source
preimage `[4^{-1} mod r]T001` with torsion offset 0. The same subgroup
point and source preimage were used in all 16 rows.

| Own-degree gate | Prefixes | Distinct columns | Rank | Contradictions | Rows involving only source bits | Outcome |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| No source prefix, preceding [direct-system result](../f6_n83_sixsum_wide_20261005/RESULT.md) | 1 | 384,895 | 415 | 0 | 0 | One affine row with two free intermediate bits |
| One complete valid first source summand, this frozen panel | 16 | 362,063–362,162 | 415 in every row | 0 | 0 | 1 affine row for index 0; 2 for each other index |

All 16 source abscissae were distinct. The common affine row is
`v80 ⊕ v345 ⊕ v426 = 0`; the second row when present has support
`[v16,v96,v177]` and a constant bit determined by the fixed source
code. Both include free intermediate variables. Each row alone can be
satisfied for either value of its source bit, so these rows do not
exclude a valid first source point or constrain the remaining source
codes without additional equations. The success criterion of at least
one contradiction or source-only affine row was **not met**. This
statement is limited to the frozen 16 prefixes and one own-degree
reduction. It is not a proof that no stronger algebraic method could
prune or that T001 has no decomposition.

The 16 reductions used 248,567,323 counted word XOR operations in the
internal elimination phase, with 1,370,810–1,395,530 source-system term
occurrences after substitution. These are stage counters, excluding
base setup, symbolic construction, and all IC phases. The process
exited 0 in 10 elapsed seconds; peak sampled RSS was 490,720 KiB and
`getrusage` peak RSS was 504,397,824 bytes. Another Rust build was
active on this host, so the elapsed times are contended feasibility
observations only. No CPU speedup or isolated performance result is
reported.

Source commit for the run:
`c71298b99a579a712a1bef30412f3cbe36a3c107`. SHA-256 of the
512-bit solver source:
`a4fc551dc19bfebd1be1c0cf1b85a0aa304d64e01bf1391e7d99f44eb12690e8`;
probe source:
`6d3d9b80fb1323b6e5280329ca68639e2e4f6ee7976a28de2cc3ee03e92a4ac6`;
protocol:
`14164499e76db3041c0f36d07a4b2ffd07b7581c153675bb90492b3de794772b`;
native executable:
`1ca58295c78227cbdd34c56bd15bcc77c9271d7b5a9ac5603343dff7e8b68273`;
raw JSONL:
`b0051619234e91de80f113a79cfa80a907d23809d571f4f9fe32b51ea293c4a9`.
The focused substitution tests passed 3/3 and the example compiled.
The host was an Apple M4 Pro, Darwin arm64, Rust 1.93.1; it had no
auditable exclusive CPU partition.

Rebuild and replay with the accepted native Rust dependency set:

```sh
CARGO_INCREMENTAL=0 CARGO_TARGET_DIR=/Volumes/SSD990/crypto/worktrees/f6-n83-geometric-20261004/target cargo build --offline --release --example f6_n83_sixsum_prefix_probe
./research/f6_n83_sixsum_prefix_20261005/run.sh panel_1
```

Decision: stop the one-prefix root reduction as this F6 pruning path.
Fixing two full source points would couple the first S3 link more
strongly, but enumerating all unordered source pairs would already
require `C(64,907+1,2) = 2,106,491,778` pair instances. A follow-on
method needs a justified way to convey that two-point constraint without
materialising or querying the full pair space. There is still no
ordinary verified six-summand relation, factor-base relation yield or
rank, target scalar recovery, complete IC candidate ID, matched
one-target rho run, or measured F4/F5/F6/IC speedup.
