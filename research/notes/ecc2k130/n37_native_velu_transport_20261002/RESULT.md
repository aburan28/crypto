# Native signed transport prerequisite for the archived degree-73 descent

Status: **correctness preflight passed; no PDP or full-cost result**.

The hypothesis was that the archived degree-73 kernel polynomial could be used
in Rust to map complete points, preserving the sign of an unknown target. The
existing `velu_x_map` yielded only an abscissa. Picking a separate lift for a
base point and a target would leave their relative sign unknown, so a signed
relation or recovered logarithm could not be transported honestly.

The frozen input is
[`raw.json.gz`](../../../koblitz_isogeny_descent_37_results_20260925/raw.json.gz),
SHA-256 `eb4773d556886672e8b80735fc486c16b73765fd5069ff81133100bfc3eafa90`,
with its original
[`contract.json`](../../../koblitz_isogeny_descent_37_20260925/contract.json),
SHA-256 `ba5e9844517ca00d1a798b2caf4552ac37ffab73e875dfcba1e04d5fa89ed8c0`.
The test reads the archived kernel and normalized codomain directly from that
gzip JSON. The field polynomial is `x^37+x^5+x^4+x^3+x^2+x+1`; the source
group order is 137439487532 and the designated subgroup order is 230603167.
The certificate records endomorphism-order discriminants `-7` on the source
and `-37303 = -7·73²` on the descendant. This is a real order change, but
by itself it neither increases useful relation yield nor gives the descendant
a cheap same-curve Frobenius action. That remains a costed policy question.
The old Sage/Python archive remains historical evidence; this implementation
and both new correctness tests run in Rust.

`velu_point_map` evaluates the ordinate without splitting the kernel. For
`h(x)=∏(x+u_i)`, it evaluates `h'/h` and the second and third Hasse
derivatives. Newton identities recover the cubic reciprocal sum needed by
Vélu's ordinate formula. The map sends a kernel point to infinity and returns
a fully signed point on the normalized codomain.

The acceptance checks for this correctness prerequisite were any mismatch in
archived codomain coefficient, a point off the codomain, failure to preserve
negation/addition, or a transported scalar mismatch. These did not occur:

| Native check | Outcome |
|:--|:--|
| Four degree-3 kernels over `GF(2^8)` | Every rational input point maps on-curve, has the expected abscissa and sign, and matches a direct pointwise Vélu ordinate sum; sampled pairs preserve addition. |
| Archived degree-73 kernel over `GF(2^37)` | Computed normalized `a6` and `y_shift` match the archive; eight deterministically selected subgroup points map on-curve with preserved order and sign; seven adjacent pairs preserve addition. |
| Planted scalar correctness probe | For `d=12345678`, `φ([d]P)=[d]φ(P)` as complete points. This is a correctness test only; the scalar is not a PDP oracle input or attack result. |

Replay with `cargo test --lib cryptanalysis::binary_velu::tests::full_point_map_preserves_sign_and_addition_on_small_curve` and `cargo test --lib cryptanalysis::binary_velu::tests::archived_degree_73_descent_has_oriented_native_point_map`. No speed or stage-cost number was measured; setup, target transport, relation search, solver, rank, and matched rho costs remain **unknown** for this path.

The next scored comparison has an additional identity gate. The archived map
uses low field-polynomial bits `{0,1,2,3,4,5}`; the frozen n37/L1024 source
control in `disjoint_cold_v2_20261001/FROZEN.json` uses `{0,1,4,6}` and its
registered source model is `icv1-f2m37-tm534059-32aad96b`. Although both
fields have 2^37 elements and the same abstract group order, their bit-packed
points are **not interchangeable**. A native, tested field-basis isomorphism
must transport the frozen generator and held-out targets before a fair
original/descendant-native/transported/pullback factor-base comparison can
begin. Using a planted logarithm to recreate a held-out target would leak the
answer and is excluded.

The next bounded experiment should build that basis bridge, verify both
directions on random field elements and subgroup points, and freeze a
target-blind four-policy base at equal useful projected sign-class size.
Because the current descendant-native `K=42` base has a counting support
ceiling of 16.10084% at arity five, the first viable native trial is arity
six; that counting threshold does not predict yield or solver cost. Admit a
policy only after independent signed relation checks and rank progress on
held-out targets, then charge every setup and online phase against the
existing matched automorphism-aware rho reference. A failed bridge or failed
rank gate is a reported result, not an inferred speedup.
