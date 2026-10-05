# n83 K0 mixed-base six-summand capacity screen

The follow-on [Frobenius-closure result](ORBIT_RESULT.md) found a more
column-efficient high-arity proposal on the same curve: 332,166 actual
usable points in 2,001 signed Frobenius columns. It remains a structural
screen with no ordinary relation or complete solver measurement.

The [registered protocol](PROTOCOL.md) used the exact K0 confidence-gate
curve from PR #1341. One native Rust process constructed standard polynomial
subspaces, enumerated their curve points, projected every point by the
public cofactor, removed identity and duplicates, and checked that the
dimension-12 projected set is contained in each larger set. Every row
completed. The dimension-12 control reproduced **4,054 usable points and
2,027 signed columns**. The [raw JSONL](raw.jsonl) contains exact integer
counts and BLAKE3 digests of the projected point sets.

| Larger base dimension | Geometric points | Usable points `B` | Signed columns | Mixed `4B+2B₁₂` combination count / subgroup order | Six points all from larger base, same ratio |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 14 | 16,125 | 16,122 | 8,061 | 0.009572844 | 0.010096 |
| 15 | 32,657 | 32,654 | 16,327 | 0.161075343 | 0.696716 |
| 16 | 64,907 | 64,904 | 32,452 | 2.513794177 | 42.950379 |

The fixed terminal base has `B₁₂ = 4,054`. The mixed numerator is exactly
`C(B+3,4) × C(B₁₂+1,2)`; the all-large numerator is `C(B+5,6)`.
Each table ratio is an **upper bound before capping at one** on coverage of
a uniformly selected subgroup target. The mixed roles overlap because the
smaller base is a subset of the larger base, and different multisets may
have the same group sum. These effects can only lower actual coverage.
The exact subgroup order is `2,417,851,639,230,796,216,685,689`.

The registered 1% *capacity* gate fails at dimension 14 (0.9573%) and
passes at dimensions 15 (16.11%) and 16 (count above the subgroup order).
**No ordinary relation was sought or found.** A count above the subgroup
order gives no positive lower bound on coverage or solver success. This is
a proposal, with no complete candidate ID, online time, cold cost, or
F4/F5/F6/rho speedup. The point-set digest uses lexicographically sorted
33-byte encodings: one affine tag byte, followed by two little-endian
`u64` words each for `x` and `y`.

For comparison of *problem sizes only*, an ordinary six-summand S3 chain
with four dimension-16 coordinates and two dimension-12 coordinates has
`4×16 + 2×12 + 4×83 = 420` Boolean variables. The existing eight-summand,
all-dimension-12 chain has 594. The mixed design raises the potential
signed relation columns from 2,027 to 32,452, and the current `u64`/
`u128` monomial backends cannot directly represent the 420-variable
system. Thus this screen identifies a smaller high-arity algebraic system,
not an implemented faster solver or a cheaper end-to-end IC method.

Run command, with the protocol and source committed before execution:

```sh
RAYON_NUM_THREADS=1 CARGO_TARGET_DIR=/Volumes/SSD990/crypto/worktrees/f6-n83-geometric-20261004/target timeout 900s cargo run --offline --release --example f6_wide_n83_k0_mixed_base_inventory > research/f6_n83_mixed_base_20261004/raw.jsonl 2> research/f6_n83_mixed_base_20261004/run.stderr
```

The first online Cargo attempt failed before execution because the crates.io
index was unreachable; the same frozen source then built and completed with
`--offline`. Its stderr was overwritten by the successful rerun, so
[startup-failure.stderr](startup-failure.stderr) is an exact transcription
of the captured tool output rather than an original raw log. The successful build log is
[run.stderr](run.stderr); the inventory process exited 0. Source SHA-256:
`87bfa6c28f2dace302a4a9fdac9be74676a9da77f5422d6c1e7d5f2ac75b3efa`.
Protocol SHA-256:
`fd878d8fa90c14ded5de12a87012c06480ed181e8406cc247a349b065d257d72`.
Raw JSONL SHA-256:
`2158d64cb5c00f266399898581ca403edeb3c4f7591b10ee7207d7ea08f896d1`.
No CPU performance ratio is reported because this was an unisolated
structural run.

The next solver gate remains an exact ordinary, independently verified
six-summand K0 relation under a frozen budget, including failed attempts
and memory. A compositional solver must communicate proved constraints
between small algebraic blocks and the exact terminal pair index; a
budget stop must have an explicit inconclusive outcome. The larger base's
relation rank and final sparse linear algebra cost must then be measured
before a full one-target online IC comparison.
