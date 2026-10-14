# Run report: cubic equivariant coordinate

UTC run start: 2026-10-14 05:51. Frozen protocol commit `b9f3259739`.
Host: macOS Darwin arm64, shared/L0. The full rational-point sweep is
deterministic; no random seed, wall-time comparison, or confidence
interval applies.
The [source](../../equivariant.rs) was compiled with native Rust and run
once on all rational points of the four source curves. The [raw TSV](../../equivariant-v1.tsv)
has four rows; the [run log](../../equivariant-v1.log) records kernel
values and nonfactorization checks. Every requested identity assertion
passed in the original run.

The verifier checked two translation identities on respectively 618,
2,634, 9,678 and 38,562 rational points; `f³=id` was checked on the
same number of projective-line points. The [diagram](../../equivariant.svg)
and [run PDF](report.pdf) explain the finite result. The
[result report](../../REPORT.md) gives the algebraic derivation, proof
status and next experiment. No solver or performance measurement was made.

Source SHA-256:
`54345e8f02f704ca2992171a21b5caa417f40a7ddb944e7455b79acd0e618b79`.
Binary SHA-256:
`3eead58c31b8b8a882af9e71bf033e5b594481051dd4e1014694b950005effd4`.
Raw TSV SHA-256:
`812b41e58d6d4608d6f92e5bdc4041ca37d0fd708b03cbdd61ff5fea73de2e74`.

The first replay attempt exited 127 because its binary had disappeared
from `/private/tmp` between invocations. Its empty
[TSV](../../equivariant-replay-v1.tsv) and [error log](../../equivariant-replay-v1.log)
are retained as an infrastructure failure, not a coordinate failure.
Recompiling the unchanged source reproduced the binary SHA-256 above.
The corrected [replay TSV](../../equivariant-replay-v2.tsv) and
[replay log](../../equivariant-replay-v2.log) each match the original
byte-for-byte (`cmp` exit 0); the replay TSV has the same SHA-256.
This is a same-host replay, not an independent-host check.

Exact commands, from repository root:

```sh
rustc --edition=2021 -O research/isogeny_torsion_fibers_20261012/torsion_orbit_20261013/equivariant.rs -o /private/tmp/isogeny-torsion-equivariant-20261013
/private/tmp/isogeny-torsion-equivariant-20261013 > research/isogeny_torsion_fibers_20261012/torsion_orbit_20261013/equivariant-v1.tsv 2> research/isogeny_torsion_fibers_20261012/torsion_orbit_20261013/equivariant-v1.log
```
