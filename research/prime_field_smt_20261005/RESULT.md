# Prime-field SMT result

## Outcome

The exact frontend and all three encodings are implemented. Planted SAT and
exhaustive UNSAT controls pass, but every frozen 20-bit SMT formulation
returned `unknown` at the fixed 120,000 ms cvc5 limit. The functional
hypotheses that the frozen queries would complete were therefore falsified.

The native exact reference found and independently verified the frozen
decomposition

```text
(137, 533365) + (3709, 528803) = (36271, 187981)
```

with slope `500576`. It enumerated 3,992 signed factor points and found the
witness after 146 membership probes.

| Frozen 20-bit method | Query variables | Query bytes | cvc5 wall time | Outcome | Pollard-rho ratio | End-to-end S |
|---|---:|---:|---:|---|---:|---:|
| native exact `R-P` reference | N/A | N/A | not timed | verified SAT, 146 probes | null | null |
| widened `QF_BV` coordinates | 5 bit-vectors | 1,976 | 120,018 ms | `unknown` | null | null |
| native `QF_FF` coordinates | 29 field variables | 4,747 | 120,011 ms | `unknown` | null | null |
| native `QF_FF` Semaev `S3` | 26 field variables | 125,520 | 120,028 ms | `unknown` | null | null |

The S3 query scanned all 4,096 candidate abscissae and admitted exactly
1,996 liftable values, recorded by digest in its receipt. Its larger query
is the explicit membership trie; it removes `y1`, `y2`, and `lambda` from
the algebraic relation but did not make cvc5 finish this instance.

## Correctness controls

- Eight focused Rust tests pass, covering exact integer/model parsing,
  overflow-safe export, full-coordinate and S3 export, signed S3 lifting,
  independent group-law replay, SAT/UNSAT exhaustive agreement, and invalid
  input rejection.
- The `QF_BV` toy SAT model verified in 51 ms; its toy UNSAT answer agreed
  with the complete native oracle in 35 ms.
- The corrected native-field coordinate controls produced a verified SAT
  model in 102 ms and the expected UNSAT result in 20 ms.
- The S3 controls produced a verified SAT model in 10 ms and the expected
  UNSAT result in 15 ms.

Development compatibility evidence is retained rather than erased: the
non-CoCoA cvc5 artifact rejected `QF_FF`; the first CoCoA SAT control exposed
a model-literal parser form; and the first one-bit UNSAT control exposed an
invalid unary `ff.bitsum`. Those issues were corrected before the frozen
native-field trials, and the failed control directories remain archived.

## Interpretation

This answers the prime-field WDSat question negatively in its literal form:
`F_p` has no nontrivial extension basis to descend. A direct modular
SAT/SMT formulation is possible and is now available, but general-purpose
cvc5 was already inconclusive at 20 bits while the specialized exact
two-summand lookup completed immediately. The present solver is useful as a
sound research oracle and as a base for higher-degree decomposition work;
it is not evidence for a method faster than Pollard rho at 51 or 83 bits.

No relation collection, sparse linear algebra, or discrete-log recovery was
performed. Consequently speedup, Pollard-rho ratio, and end-to-end `S` stay
`null`, and no scoreboard or leaderboard entry changes.

## Reproducibility

The cross-repository input is the exact decimal translation of
`aburan28/cryptanalysis` commit
`a14ca487a5659fc6e6572133f7c89a3326b3e1b5`, path
`challenges/ecc/curves/fp-random-b20.json`, Git blob
`40b380492c7a250babaa2b5ae86a4cb5592fe598`. The bound `B=4096` was frozen
before solving.

Native-field runs used the official cvc5 1.4.1 CoCoA-enabled static GPL
archive with SHA-256
`d0b54324ec2129697975da8753767fd255309947b87832917d32a78da7d16666`.
Every run receipt records the executable BLAKE3 before and after execution,
the exact arguments, hashes of all byte streams, and the independent
verification decision. Raw artifacts are under `results/` and are never
overwritten.
