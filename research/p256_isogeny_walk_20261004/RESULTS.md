# P-256 degree-11 isogeny-walk feasibility result

Status: **gate passed**

Protocol: [`PROTOCOL.md`](PROTOCOL.md)

Date executed: 2026-10-04

Result class: **correctness and feasibility diagnostic, not an ECDLP
speedup**

## Result

The native Rust walker emitted exactly 4,096 deterministic,
non-backtracking degree-11 edges from registered P-256. The in-command replay
and a second standalone replay verified all 4,096 edges. The path contained
4,097 distinct `j`-invariants including the starting vertex, so no cycle was
encountered in this prefix.

This passes the frozen bounded-walk gate. It does not establish an anomalous
curve, improve an index-calculus stage, compare a complete solver against rho,
or authorize the proposed `2^40` campaign.

| quantity | result | boundary / interpretation |
|---|---:|---|
| requested emitted edges | 4,096 | frozen exact count |
| emitted edges | 4,096 | pass: ratio to requested = 1.000 |
| replay-verified edges | 4,096 | pass: ratio to emitted = 1.000 |
| failed edges | 0 | pass: required zero |
| distinct `j`-invariants | 4,097 | includes the starting P-256 vertex |
| first detected cycle | none | no cycle within the measured prefix |
| generation time | 955.011 s | stage practicality diagnostic |
| generation throughput | 4.28896 edges/s | not an attack speedup |
| isolated command wall time | 961.029 s | includes generation, replay and output |
| standalone replay wall time | 5.857 s | all edges checked again |
| serialized edge bytes | 3,230,565 | 788.712 bytes/emitted edge |
| compressed certificate bytes | 744,997 | SHA-256 in `MANIFEST.json` |

Both isolated runs report `contended: false`. They ran on CPU 4 of the
five-vCPU environment, an AMD EPYC 9V74 host (`x86_64`, Linux 6.18.44), using
`rustc 1.98.0 (88d9e12ae 2026-08-18)`. The generation command used 16,528 KiB
maximum RSS; the standalone replay used 17,588 KiB.

## Degree census

The native Frobenius census corrected the earlier idea that the first small
edge automatically supplies a walk. Degrees 3 and 5 are ramified: each has one
`F_p`-defined direction and therefore forces the dual edge back. Degree 11 is
the first split prime and the first usable non-backtracking degree.

| `ell` | `t^2 - 4p mod ell` | class | `F_p` edges | walkable |
|---:|---:|---|---:|:---:|
| 2 | 1 | no rational 2-torsion | 0 | no |
| 3 | 0 | ramified | 1 | no |
| 5 | 0 | ramified | 1 | no |
| 7 | 6 | inert | 0 | no |
| 11 | 3 | split | 2 | **yes** |
| 13 | 10 | split | 2 | yes |
| 17 | 15 | split | 2 | yes |
| 19 | 15 | inert | 0 | no |
| 23 | 1 | split | 2 | yes |
| 29 | 23 | split | 2 | yes |
| 31 | 15 | inert | 0 | no |

## Verification boundary

Each edge carries the target short-Weierstrass model, `j`-invariant,
predecessor, degree, modular-polynomial identity, and a nonidentity target
point `Q`. Replay performs all of the following:

1. recompute the target `j` and reject singular or noncanonical models;
2. evaluate the pinned classical `Phi_11(source_j, target_j)`;
3. check that `[n]Q = O` for P-256's prime order `n`;
4. use the Hasse upper bound `< 2n` to conclude that the target curve's order
   is exactly `n`, rejecting the opposite-trace quadratic twist;
5. enforce the exact start, edge continuity, first-root rule, no immediate
   backtracking, record count and edge-stream digest.

The edge-stream SHA-256 is
`a2b82f718c7436368268d3d5cbbb5f9b105504a951eddd61c4f7b11ffdb2af2c`.
A modular-polynomial root without the exact-order point is not accepted.

## Independent start-edge fixture

The slower cross-check constructs `psi_11`, obtains the two Frobenius
eigenspace kernel polynomials from the x-coordinate equations for eigenvalues
2 and 8 (represented by scalars 2 and 3 modulo sign), verifies closure under
doubling, and applies Vélu's coefficient sums. It produces the same two
neighbour `j`-invariants as the published classical modular polynomial:

| eigenspace scalar | degree-five kernel constant | codomain `j` |
|---:|---|---|
| 3 | `f136e7e36906a7ed31e3dc2ef549d5e0001c8c6237f7205b0d86d62ddfacb549` | `1e83ce42e311c91e26e534220f6c8af76ef81a2c22246611462ef9b5b044ff88` |
| 2 | `c066675b2ac76605dbb2aa7c9a20755a5906ef4eb6babd4be3d19dcf720b014d` | `f17b32fb5796946d1d2b09fd642c7661428bf217216ae83fb6354a5c87fefbc0` |

These values are frozen in the Rust test rather than learned from the
4,096-edge certificate.

## Reproduction

```bash
cargo run --release --bin p256_isogeny_walk -- fixture
cargo run --release --bin p256_isogeny_walk -- \
  generate --steps 4096 \
  --output research/p256_isogeny_walk_20261004/certificate.jsonl.gz
cargo run --release --bin p256_isogeny_walk -- \
  verify \
  --input research/p256_isogeny_walk_20261004/certificate.jsonl.gz
```

The measured commands were wrapped with `tools/isolated_bench.py run`; exact
commands and host conditions are in `isolated-run.jsonl` and
`isolated-verify.jsonl`. `MANIFEST.json` pins every committed artifact.

## Decision

Gate 1's bounded native walk is feasible and reproducibly verified. The next
permitted step is a separately preregistered bounded storage pilot; this result
does not satisfy the native payable-work witness requirement for a Cairn
objective, and it does not justify production scale.

No scoreboard or leaderboard panel changes: this experiment measured a walk
construction stage and made no complete one-target ECDLP cost or speedup claim.
