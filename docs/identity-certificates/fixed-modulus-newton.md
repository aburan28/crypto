# Fixed-modulus Newton reciprocal: proof and controls

TASK-20261016-de577e; cached-fixed-modulus-reduction is off by default.
No performance result or expanded-map score follows from this source change.
The native identity control issues and rechecks the adjacent IDC1 record.
Native identifier: IDC1h41e9973f839c85bf. Issuance and32 extra domain-separated checks passed in implementation-02/newton-identity.json; independent source review remains pending.

## Exact arithmetic argument

Let K be a field, m a nonzero trimmed polynomial of degree d, L=d+1,
and B(X)=X^d m(1/X). The constant coefficient B_0 is nonzero. There is a
unique formal inverse: G_0=B_0^-1 and, for k>0,

    G_k = -B_0^-1 * sum(i=1..k, B_i * G_(k-i)).

This recurrence proves existence and uniqueness inductively. Newton's
update H=G(2-BG) satisfies

    B*G*(2-B*G) = 1 - (1-B*G)^2.

If BG=1 modulo X^k, every squared-residual coefficient below 2k vanishes:
for i+j<2k at least one of i,j is less than k. Thus H agrees with the
unique inverse through precision 2k. Truncate the last doubling to exactly
L. Zero coefficients are part of formal precision and must not be removed.

For a trimmed dividend a of degree n>=d and count=n-d+1, let A be the
degree-n reverse of a. The reversed quotient is

    Qrev = A * (B^-1 mod X^count) mod X^count.

Reverse exactly count coefficients of Qrev to obtain q, retaining zeros
at the reversal boundary. The top count coefficients of a-qm cancel;
deg(a-qm)<d. Quotient and remainder are unique by division over K.

By uniqueness, the precision-L inverse has exactly the same first count
coefficients as an independently computed precision-count inverse.
The preserved predicate is L>=128 and L<=a.len()<2L, so count ranges
from 1 through L. Reusing the stored prefix therefore gives the same q
and remainder, including a.len()=2L-1. Every other length uses unchanged
schoolbook divrem. A nonzero constant modulus takes that original path;
zero or trailing-zero moduli are rejected explicitly at construction.

The context borrows its exact field and modulus; methods accept only
operands. There is no modulus clone or field/modulus/inverse replacement
API. Construction is charged once per verification, after degree and
squarefree checks, with its own durable trace stage.

## Proof map and limits

| Obligation | Evidence and status |
|---|---|
| Scalar residual identity for all integer b,g | Checked Lean FixedModulusNewton.newton_scalar |
| Residual vanishing implies doubled square precision | Checked Lean residual_square_precision, integer sequences and explicit finite convolution |
| Formal inverse existence/uniqueness over arbitrary K | Written recurrence above; complete Lean formal-series bridge open |
| Quotient reversal equals division | Written cancellation/uniqueness argument; complete Lean polynomial bridge open |
| Actual field coefficients, inverse precision and quotient | Native residual, corrupt-inverse, threshold and schoolbook controls |
| Selected recurrence and modular composition | Preserved reference plus dense recurrence and Horner controls |
| Public admission and all20 transports | P192/13, P521/7, P224/1471 controls; full A/B mathematical byte comparison |

Lean source and pinned toolchain are in
research/isogeny_world_frontier_20261016/replay-optimization-v1/proof/.
Lean4.35.0, commit eb6cafb7290e8bd5f2937522d8d804f58f02bafc, checks
both declarations with only standard propext and Quot.sound axioms.
No admitted facts, custom claim axioms, kernel bypasses or native evaluation
are used. The environment lacks lean4checker; a separate fresh checker was
not run. The first failed source and log remain retained.

The IDC1 evaluation is over the module's fixed field p=2^61-1. Deterministic
evaluation points are not a proof over integer coefficients; no random
false-accept probability is claimed. Its proof pointer stays empty because
the adjacent Lean results have distinct scopes. Exact coefficient tests
cover the declared finite grid, not all possible Rust inputs. Full
formal-series and reversal formalization remains open.

Scientific timing, sealed50341 replay, default activation and board
promotion require later frozen approval and independent review.

