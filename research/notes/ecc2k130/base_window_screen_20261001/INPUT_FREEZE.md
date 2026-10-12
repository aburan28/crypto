# Frozen held-out input set

**Status: input freeze; no scored arm has run.** This is the five-block public-Q
workload for the [preregistered four-window screen](PROTOCOL.md), anchored to
the [source lock](SOURCE_LOCK.md) merged as [PR #1129](https://github.com/aburan28/crypto/pull/1129).
It is a batch pre-selector diagnostic on `icv1-f2m41-tm2308219-7f48b14a`,
not the repository's primary one-target IC-versus-rho comparison.

The immutable preparation tree and source-lock merge are both
`b97d77a0932abd57feaa2d3734fcdf68193a694d`. Nineteen source files are
byte-pinned to that commit. The inventory contains all 121 tracked prior
point corpora at that tree, with digest
`56647a01fd54be86b0c4ece88d2e891ad6839262685cb97ffdc4f0b30b2490b2`.
The independent replay checked 18,446 earlier n=41 point rows and orbit
keys. Future point corpora do not alter this historical inventory because the
verifier reads the committed blobs from the recorded tree.

| Block | Public Q | Candidates tried | Rejected candidates | Point-only SHA-256 |
| ---: | ---: | ---: | ---: | --- |
| 0 | 1,024 | 1,024 | 0 | `13d8f931a7f5cc9d33e71930049f01b9a19b751fe4324a7d2b44579801c7f340` |
| 1 | 1,024 | 1,024 | 0 | `952f885f2258508eb486f7978d869e55d6d2523b8141f44e406905b3724f2707` |
| 2 | 1,024 | 1,024 | 0 | `2193c5d208f002d96833021b4eb038aa38f99d87de8115bdc36fb7740477d5a6` |
| 3 | 1,024 | 1,024 | 0 | `5ce44c0d7129ac92f47e08e16ddb41af956307412d21bda9977e2fd08492b6ab` |
| 4 | 1,024 | 1,024 | 0 | `4757d3fdd76a8c2eaa4795d362ec7f2d85a04c522c295bbcff780e45e65aff69` |

`FROZEN.json` SHA-256 is
`6487e1c30ef1ae1db5cb99fd39297150b1a327baaff727e3966d3115bc88e2a7`;
`INPUT_RECEIPT.json` SHA-256 is
`570f80c9bc5939b06b4cde81d3ee2c8bf2d93a568814d3717606bc483516d0e4`.
The latter records a PASS after independently replaying every deterministic
candidate, all 5,120 `[d]G = Q` equations, subgroup membership, and disjoint
signed-Frobenius orbit keys. Each block's known scalars are in a separate
`.fixture.jsonl`; the scored runner opens only `.points.jsonl`.

The next outcome PR must retain a fresh frozen release build, every raw
generator/compact/rho child receipt from the one-shot 45-arm schedule, quiet
machine and A/A evidence, independent rank/witness/scalar replay, and the
frozen decision. No timing, method speedup, or n=131 transfer is asserted by
this input freeze.
