# n83 five-summand dimension-18 F6 system: exact, root pruning fails

The [preregistered structural gate](PROTOCOL.md) passes; its
search-pruning gate fails. On the registered K0 curve
`icv1-f2m83-tm6151469093347-debefd74`, the standard dimension-18
source subspace enumerates **261,447 geometric points**. Cofactor-four
projection gives **261,444 distinct nonidentity subgroup-usable points**
and 130,722 sign-folded columns. The direct five-summand S3 chain has
339 Boolean variables, 332 cubic coordinate equations, and four S3
links. It is a smaller algebraic system than the preceding dimension-16
six-summand chain, which had 428 variables and 415 equations. This is an
algebraic feasibility comparison, not a complete F6 candidate.

The planted source indices `[0,2,4,6,8]` made every equation vanish.
Their full curve sum replayed to
`(26c1d2a9079c67a42c102,4bd614a90524edd0b4e9b)`; its cofactor-four
projection was
`(28c50c9a125ad41cfb39f,298b6d142b459ceeff33f)`.
The subgroup preimage plus exactly one verified rational 4-torsion
offset (2) reconstructed that source sum. The [raw control](planted.jsonl)
and [exit status](status.tsv) preserve the check. This is a planted
correctness control and contributes no natural relation-yield estimate.

The ordinary input was the same public subgroup point T001,
`(355fb5df7a905f16921eb,5900a390f42d290f1bbe)`. For each of its
four checked torsion lifts, the native own-degree Macaulay reduction
completed under the 120-second and observed 7-GiB limits:

| Torsion offset | Term occurrences | Columns | Rank | Affine row constant | Source-only rows | Contradictions | XOR ops |
| ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: |
| 0 | 1,104,702 | 299,558 | 332 | 0 | 0 | 0 | 17,264,647 |
| 1 | 1,108,576 | 299,558 | 332 | 0 | 0 | 0 | 17,259,692 |
| 2 | 1,097,593 | 299,558 | 332 | 1 | 0 | 0 | 17,260,678 |
| 3 | 1,099,513 | 299,558 | 332 | 1 | 0 | 0 | 17,259,670 |

Every offset has exactly one affine row with support
`[v72,v256,v337]`. Only `v72` is a source-coordinate bit; `v256` and
`v337` are free bits of the third intermediate x-coordinate. The row
can be satisfied for either value of `v72`, so it does not prune a
source summand. There is no contradiction. The four raw files
([offset 0](ordinary_0.jsonl), [1](ordinary_1.jsonl),
[2](ordinary_2.jsonl), [3](ordinary_3.jsonl)) retain all affine
supports, construction times, reduction times, memory observations,
and the source-preimage coordinates; the [status table](status.tsv)
records exit zero for the planted control and all four ordinary runs.
Each process had an empty stderr file. Maximum reported process RSS was
389,283,840 bytes. The macOS kernel memory limit was not accepted in
the earlier six-summand probe, so the 7-GiB check here is observed at
phase boundaries, not a hard kernel allocation limit.

Using the **actual** usable base count, the nominal unordered
five-multiset count divided by the subgroup order is
`C(261444+4,5)/2417851639230796216685689 = 4.21016`.
This is a counting capacity before duplicate sums and solver cost;
no ordinary decomposition, relation yield, or rank from natural
queries was measured. The positive capacity cannot turn a root-only
matrix row into a witness.

The exact native release build is in [`build.log`](build.log). Binary
SHA-256:
`953b9b089a9f9aad86a5685ccdbffa2a62da1eca565ff053bd6cdf159111fd32`.
Probe source SHA-256:
`e788e6541777b58d99001dd4a0f02a0398457576721dd317ec43ac7b57b02a8c`.
The protocol and [runner](run.sh) hashes are
`a944763cacb84b7db95947ad08983db3e6dd974266985d3e96386d9bffc7842e`
and `742c5de3a2fe3dbbff49fcebb95c4f56949188f3e1004e2801b09df78b522946`.
The [SHA-256 receipt](SHA256SUMS) verifies every raw artifact.
Rebuild with the shared offline release target and run
`sh research/f6_n83_fivesum_d18_20261005/run.sh`.

Physical host: Apple M4 Pro, arm64 macOS 26.6, Rust 1.93.1. It had no
auditable exclusive CPU partition, so the observed elapsed times are
feasibility diagnostics only. The target-independent base setup was
repeated in every process and is not an IC online interval. A complete
IC candidate manifest and ID remain unset; `candidate_id: null`.
There is no verified ordinary relation, F4/F5/F6 complete-call timing,
target scalar recovery, matched one-target rho run, or IC speedup.

Decision: retain the generic five-summand system as exact algebraic
admission evidence, and stop the own-degree root-only search path for
this arity and base. A follow-on method must produce constraints on
source-coordinate bits or a full-group-verified ordinary relation before
any F6 throughput comparison.
