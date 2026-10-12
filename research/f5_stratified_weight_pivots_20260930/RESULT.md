# Stratified current-weight F5 pivots: rejected at the structural gate

The [frozen protocol](PROTOCOL.md) preceded the opt-in candidate in commit
`dc953a352c4752b9c251bda7ee4c4252f4bc0289`. The release GF(2)
elimination tests passed 7/7 and the F5 tests passed 10/10. The candidate was
then run under the repository's shared benchmark lock on Apple ARM64.

The [raw structural receipt](structural_dc953a352.json) contains all eight
complete processes: reference and candidate for each of 16 and 32 bands on
seed XORs `0` and `badc0de1`. Every process exited successfully. All seven
cases matched the reference on rank, canonical row-space fingerprint, F5
criterion work, row and column counts, and pruning. The route indicator was
active only on the primary n24 degree-4 case. The two repeated reference
processes per seed agreed exactly.

| Seed XOR | Bands | Reference terms | Candidate terms | Candidate/reference |
| --- | ---: | ---: | ---: | ---: |
| `0` | 16 | 13,734,979 | 15,614,091 | 1.1368 |
| `0` | 32 | 13,734,979 | 14,364,207 | 1.0458 |
| `badc0de1` | 16 | 13,722,549 | 15,533,367 | 1.1320 |
| `badc0de1` | 32 | 13,722,549 | 14,343,495 | 1.0453 |

Neither candidate met the preregistered term ratio of at most 0.75 on both
seeds; both produced *more* terms. Per protocol, no paired speed measurement
was run and no speedup is claimed. The runtime option was removed. The exact
measured candidate source diff is retained as
[`measured_candidate.patch.gz`](measured_candidate.patch.gz), SHA-256
`6bab061f6bce4577eb9cd28da27c6b841cdb53c25edcbb5d5996bfb069b220b8`.
The raw receipt SHA-256 is
`26c878e24778c3ef8635dbac679fb40ad43c225439480622c3086e1a4a3940d7`.
This is a Boolean matrix-F5 solver-stage result, not an IC single-target or
Pollard-rho comparison.
