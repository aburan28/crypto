# Frozen design: native complete-support comparison

Status: design only, blocked. No executable/source/input hash has been frozen
and no task may run until the readiness checks below pass. Changes after
readiness require a new version. No measurement is reported.

## Hypothesis and arms
The cached semilinear operator lowers complete support cost after setup,
rather than only squarefree-classification cost.
A: optimized native reference using sparse coefficient squaring and monic
reduction, including three binary nilpotence squares (not eight repeated
multiplies or generic dense multiplication).
B: same arithmetic, cached C6/C8 and unchanged support decisions.
C: B plus r=0 early acceptance and proved squarefree/nilpotence dispatch.
Odd-characteristic controls retain a proved odd-field implementation; do not
apply the characteristic-two square formula.

## Inputs, independent controls and output contract
Fresh fixed seeds2026101701/02/03 are assigned to the three shards. Before
timing, generate and hash every field, modulus/basis, domain and polynomial.
Bind the generator revision; never select cases after comparative timing.
No repeated identical H may inflate sample size. Repeated roots and repeated
queries are different regimes and must remain separate.

Oracle shard: exhaustive monic quintics over GF4/F3/F5; literal allowed-root
products; independent degree1/2 square-divisor oracle; GF9/GF25 nonrational
repeated factors; compare every intermediate remainder vector, not merely a
Boolean support result. Include the GF4 and odd-characteristic counterexamples.

Binary shard: degrees4,83,163,233,283,409,571; prefix dimension min(18,n-1)
and full-width affine-scaled variants. In particular preserve the actual
GF83 modulus X^83+X^45+X^2+X+1 and d18 case. Choose distinct H across16cases
per field, with six multiplicity partitions and unsupported roots outside
the domain. The d18 domain is implicit, never enumerated.

Odd shard: the five NIST primes plus GF9/GF25,8-element explicit distinct
domains and16distinct quintics per field, supported and unsupported cases.
No claim that a prime field has a proper Fp-linear subspace.

For every H return exact support and, when extraction is in scope, the complete
canonical root/multiplicity list under the declared domain. Original
multivariate root counts remain null. This is not a complete generalized solve.

## Cost and limits
At most three experiments; each worker one CPU,4GiB,90seconds. Native code,
oracles, harness and analysis only. No network or compute rental in a worker.
No simultaneous compilation or other benchmark. Preflight work/storage caps.
Three alternating cold pairs, plus A/A at the same batch size, must charge input,
field/domain setup, cache construction, all recurrence steps, root extraction,
independent verification, serialization and complete child lifetime.
Preserve censored runs; clear exact counts on timeout or truncated output.
A90second budget is not a prediction that the whole panel fits.

Run the same cases against the fastest available FLINT/NTL/Magma implementation
of the same arithmetic/support computation; record exact versions and semantics.
Missing external comparison leaves superiority unknown.
A predicate-only ratio cannot be multiplied into a full-solver ratio.

## Decisions and readiness
Any mismatch refutes the implementation. Any field regression rejects blanket
adoption; retain the complete panel. If gain disappears with full accounting,
return to upstream reuse/search instead of raising an unchanged timeout.
A/A spread above10percent or unqualified isolation makes wall evidence
incomplete. The end-to-end research goal still requires the immutable GATE
plus requested8x against the fastest matched genuine all-roots F4/F5/F6,
20hard completed matching cases with fastest baseline>=1second, no more than
1.10worst slowdown and the requested field-family coverage.

Before ready: implement native runner/checker, pass controls, record code and
binary digests, hash all inputs, pin a merged code commit, retain protocol Git
blob, validate this card, and confirm the executor is reachable. No task is
posted to cairn and no funding or infrastructure is authorized by this card.
