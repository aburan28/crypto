# Decision: finite-fibre recursive S3 exactness passes

The local-fibre and affine-chain claims in [PROOF.md](PROOF.md) hold for
`E: y²+xy=x³+1` over any finite characteristic-two field when every
admitted rational factor-x fibre contains **both** point signs (or the
singleton at x=0). The four exhaustive local cases prove that finite roots
of `S3(a,b,c)` are exactly x-coordinates of rational finite sums. Induction
by negating all earlier factors proves that every all-affine recursive-S3
path has a globally consistent signed tuple, and either terminal target
sign can be attained. This is a mathematical implication of the proof;
the bounded controls below independently check its implementation and
edge cases rather than extrapolating from sample counts.

The frozen source commit is `bf29e854772f8c990dbaf581c86cc1a58c15bb75`
([PR #777](https://github.com/aburan28/crypto/pull/777)); the pre-outcome
[FROZEN.json](FROZEN.json) SHA-256 is
`d87bfa84b6d6ccffe859cf044b22cbf0fc901ac71e3d0d01b35469a7760855e3`.
The pair and chain domains, stop rules, source hashes, and independent
implementation were fixed in [PROTOCOL.md](PROTOCOL.md) before the first
small-field run. The local hash-only CI and the draft PR preflight passed
before the run. One producer and one independent verifier then passed; no
failed or censored attempt occurred.

| Field degree | Rational x fibres | Ordered pairs | S3 evaluations | Distinct root occurrences | Pair mismatches |
|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 4 | 8 | 3 | 0 |
| 2 | 4 | 16 | 64 | 21 | 0 |
| 3 | 2 | 4 | 32 | 3 | 0 |
| 4 | 8 | 64 | 1,024 | 105 | 0 |
| 5 | 22 | 484 | 15,488 | 903 | 0 |
| 6 | 28 | 784 | 50,176 | 1,485 | 0 |
| 7 | 58 | 3,364 | 430,592 | 6,555 | 0 |
| **Total** | — | **4,720** | **497,384** | **9,075** | **0** |

Every irreducible modulus and Koblitz group order was checked independently.
The complete ordered-pair rows include the `a=b=0`, `a=b≠0`,
`a≠b,ab=0`, and `a≠b,ab≠0` cases. The independent verifier rebuilt the
entire root list with polynomial-product/reduction arithmetic and Euclid
inversion, then compared all 4,720 rows with the producer's bit-serial/
Fermat implementation.

| Complete chain panel | Rational x tuples | Signed point tuples | Distinct affine candidate paths | All-affine signed tuples | Signed tuples with O prefix | Signed tuples ending at O | Candidate-only paths / missing target signs |
|---|---:|---:|---:|---:|---:|---:|---:|
| n=3, m=4 | 16 | 81 | 8 | 24 | 45 | 21 | 0 / 0 |
| n=4, m=4 | 4,096 | 50,625 | 19,208 | 41,160 | 6,525 | 3,165 | 0 / 0 |
| n=5, m=3 | 10,648 | 79,507 | 37,485 | 75,852 | 1,849 | 1,806 | 0 / 0 |
| **Total** | **14,760** | **130,213** | **56,701** | **117,036** | **8,419** | **4,992** | **0 / 0** |

The O-prefix and O-terminal columns are separate, possibly overlapping
properties of signed tuples; they must not be added as disjoint classes.
For every recorded candidate path, the independently rebuilt signed sums
covered the **entire** rational terminal x fibre, not merely one sign.
The first prerequisite controls were also found exactly as frozen:
nonliftable `n=3,a=b=2` has formal S3 root `c=3` but no rational factor
lift; restricting `n=2,a=1,b=2` to points `(1,0),(2,0)` retains sum x=3
but loses formal S3 root x=2. Neither is a counterexample to the
sign-complete rational-fibre theorem.

| Cold child | Wall s | CPU s | Peak RSS bytes | Native operation counters |
|---|---:|---:|---:|---|
| Producer | 2.878 | 2.783 | 94,912,512 | 4,011,414 field mul; 1,858,424 square; 241 inverse; 329,528 point adds |
| Independent verifier | 3.222 | 3.114 | 134,807,552 | 3,146,762 field mul; 861,697 square; 307,457 inverse; 329,528 point adds |

Both children are below their frozen 180-second and 512-MiB caps. The
verifier's Euclid inversion is uncached; the producer caches Fermat
inverses, so these counters are audit data, not a solver benchmark.
[The evidence archive](evidence/README.md) retains full compressed pair
and chain rows, result and verifier JSON, stdout/stderr, process receipt,
SHA-256/byte manifest, and an archive-only `ci_replay.py --evidence` check.
The receipt SHA-256 is
`56c0ec80d605bec68411fb76a02647c271bd058799705e1eae2978dcb83cf76b`.

**Decision.** The affine recursive-S3 equations are semantically exact for
finite point sums on complete rational x fibres. They can still miss a
valid PDP witness whose prefix equals O, as [PR #774](https://github.com/aburan28/crypto/pull/774)
demonstrated for a fixed n=13 target; an O target also needs explicit
handling. A factor-x mask with an O witness might have another all-affine
witness, so that exceptional witness alone does not prove a mask is missed.
An exporter must encode the infinity/inverse branches and be checked
against the full point oracle before treating absence of an affine model
as full PDP UNSAT. This theorem does not validate a direct higher S_m
resultant, denominator-cleared Boolean equations, a sign-incomplete base,
solver timing, n=131 relation yield, full ECDLP cost S, or a matched-rho
speed ratio. Those quantities remain unset.
