# Frozen symbolic full-point addition control for K0

Status: protocol and source freeze before any outcome run. This is a bounded
engineering gate for the ECC2K-130 index-calculus PDP program, not an ECDLP
experiment. The upstream width audit in [PR #796](https://github.com/aburan28/crypto/pull/796)
found that the current monomial/Boolean exporter has a `u64` variable ceiling;
it has **not** admitted a symbolic n131 rotated PDP. This PR tests whether a
different expression representation can at least state a complete full-point
addition edge without that ceiling.

## Hypothesis, relation, and reference

For `E: y²+xy=x³+1` over `GF(2^n)`, represent a point as `(o,x,y)` where
`o=1` is only the canonical `(1,0,0)` infinity and `o=0` requires the curve
equation. The DAG has `7n+3` Boolean inputs: three point flags, six n-bit
field coordinates, and one n-bit existential slope `λ`. It uses only constants,
inputs, XOR, AND, and hash-consing. Input IDs and node IDs are unrestricted
Python integers; a model is an exact-length vector of 64-bit limbs, rejecting
extra/omitted limbs, nonzero padding, and out-of-range coordinates. This is a
solver-neutral expression graph; it contains no CNF/ANF exporter or solver.

For valid points P,Q,R, the output predicate must hold for **some** λ exactly
when `R=P+Q`. The five derived, disjoint cases are P=O (copy Q), P≠O and
Q=O (copy P), finite same-x inverses (result O), finite same-x equal points
with nonzero x (doubling), and finite different-x points (ordinary addition).
The relation has an explicit branch-cover constraint and validates P, Q, and
R, including canonical O. It imposes the following denominator-free equations
inside the slope branches:

| Branch | Slope equation | Output equations |
|:--|:--|:--|
| double | `λx=x²+y` | `xR=λ²+λ`; `yR=x²+(λ+1)xR` |
| different x | `λ(xP+xQ)=yP+yQ` | `xR=λ²+λ+xP+xQ`; `yR=λ(xP+xR)+xR+yP` |

Copy and inverse branches leave λ unconstrained; slope branches must have one
λ each. The relation explicitly handles `x=0` self-inverse 2-torsion through
the inverse branch. The independent reference is the standard affine binary
curve law in `src/binary_ecc/curve.rs`, implemented separately in the verifier
with polynomial long multiplication/reduction and Euclidean inversion. The
producer uses a different shift/reduce product and brute-force inversion.

## Frozen inputs and falsification

The exhaustive fields are exactly `(n,modulus)=(2,0x7),(3,0xb)`, with
polynomial basis bit 0 as constant term. Check irreducibility before building
the DAG. Enumerate every canonical rational curve point P,Q,R and every
`λ in GF(2^n)`, in lexicographic point order with O first. Store one raw row
for every `(n,P,Q,R)` with the number of satisfying λ and the reference
answer. The exact expected witness count is `2^n` for the correct R in copy
or inverse branches, `1` for the correct R in slope branches, and `0` for all
other R. Every case must occur in at least one frozen field.

Enumerate every other raw `(o,x,y)` triple in each field and place it in each
of P,Q,R while the other points are O. Every λ must reject it. Reject seven
malformed-model/coordinate probes per field. Separately check a synthetic
920-input, 15-limb Boolean DAG: the highest bit is read correctly, and a
missing limb or nonzero padding is rejected. This is a **width-decoding**
control only; it is not an n131 curve circuit or n131 PDP.

A single wrong truth-table entry, missing branch, invalid point accepted,
malformed model accepted, irreducibility failure, reference discrepancy,
replay mismatch, timeout, or budget overrun fails the gate. Preserve the first
counterexample and every raw row. No field, pair, slope, or bad point may be
removed after seeing outcomes. The exact row count is
`Σ_n #E(GF(2^n))³`, and valid model evaluations are
`Σ_n #E(GF(2^n))³·2^n`; these counts are outcomes, not frozen constants.

## Execution and accounting

Before the first producer invocation, SHA-256 pin this protocol, all Python
sources, exact field list, domain, and caps in `FROZEN.json`; commit them and
open a draft PR. Hash-only CI must pass on that exact PR head and receive peer
review before an outcome run. Then invoke `run.py` once. It starts separate
cold producer and independent verifier processes, each capped at 180 seconds
reported child wall time, 195 seconds external wall time, and 512 MiB peak
RSS. Record UTC, command, exit/timeout, wall/CPU/RSS, stdout/stderr hashes,
source/freeze hashes, compressed full truth-table rows, first mismatch, and
independent replay receipt. Never overwrite a failed archive. Committed CI
checks artifact hashes and reruns the verifier from the raw archive in a fresh
process. A failed result is a reportable negative/censored control, not an
invitation to remove the failing case.

The success condition is exact truth-table agreement, all invalid controls
rejected, all five branches exercised, and matching independent replay within
caps. Report the actual XOR/AND node counts, model limbs, cold resources, and
raw case counts; do not quote a SAT runtime or an end-to-end speedup. The
reference boundary for eventual ECDLP claims remains matched automorphism-
aware Pollard rho with **all** index-calculus phases charged, per `AGENTS.md`.
This toy gate advances representation correctness only. The next admission
decision is whether this relation can be exported, model-lifted, and costed
under a frozen n13/n19 solver panel; n131 remains `NOT_ADMITTED` until its
factor-bit, prefix-coordinate, memory, and full-cost gates are met.
