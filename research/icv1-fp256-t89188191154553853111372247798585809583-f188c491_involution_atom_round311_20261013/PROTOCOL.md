# P-256 involution-eigenspace atom factorization, round 311: preregistration

Date preregistered: 2026-10-13

The first local preregistration attempt was made on 2026-10-09, but its
managed workspace was replaced before the branch was pushed and before any
experiment ran.  This durable preregistration preserves the same hypotheses
and gates; it does not claim evidence from that interrupted attempt.

## Objective and fixed boundary

Continue the search for a factor-base mechanism at Pollard-rho parity on
`icv1-fp256-t89188191154553853111372247798585809583-f188c491`, starting
from registered base `FB1h2f8621cda105` and Round-309 half-turn base
`FB1h6255ce9746fe`.

Round 309 constructed a real same-base involution with two useful rows but
priced fixed S17 against a product target.  Round 310 proved that the complete
fixed-sum linear image is one-dimensional only for affine-log transport.  This
round instead diagonalizes the involution in the coefficient domain and asks
whether its actual point atoms reduce the unknown-log dimension or only
relabel it.

Pollard rho is fixed at `S=operations/sqrt(n)=1.3`.  Promotion still requires:
zero errors on complete checked instances; structured residual degree at most
five; 138,031 independent usable rows below `2^120` operations; cost below
`2^103` per usable row; peak materialized storage below `2^50` bytes; exact
recovery; and no discarded probabilistic branch described as exhaustive.

## Frozen inputs

- Round 294 variable-arity result, SHA-256
  `9716fc1ba7d24a023627e8701094ca870b2121a84f7a0148784f70565640815a`.
- Round 309 result, SHA-256
  `1fa17da4d11345f2dce55d6704355d193e081eae6912d3b45a9605e2efa32bc1`.
- Round 309 complete factor-base dump, SHA-256
  `677f5ddbe0f2ea1e0c2bd29272d76461acc0d0197d11843d3ebf77d91eaf15ac`.
- Round 310 result, SHA-256
  `f8c39500b97394e633dcbcd36a92bf2cb896176aff7e4ffcb357e201b4dbfdc8`.

Reject every hash, schema, curve, factor-base, column-count, replay, or gate
mismatch.  Do not inherit a degree claim from either parent base.

## Construction under test

For each actual two-cycle `{i,pi(i)}`, oriented by smaller column index, set

```text
A_i = P_i + P_pi(i),       B_i = P_i - P_pi(i).
```

For symmetric coefficients `a` and antisymmetric coefficients `b`, the
product target separates exactly:

```text
sum a_i A_i = (Q+G)/2,     sum b_i B_i = (Q-G)/2.
```

The corresponding original-log rows are `(2*c_plus,-1,-1)` and
`(2*c_minus,-1,+1)`, whose last two columns have determinant `-2`.  Locally,

```text
A_i+B_i=2*P_i,             A_i-B_i=2*P_pi(i),
```

so atomization alone is invertible, not a variable reduction.  Promotion
therefore requires actual atom collisions, cross-channel overlap, or another
verified known-log transport.

## Actual-base census

Parse the complete Round-309 `ecbench.factor_base_dump/v1-wide` file and
independently rebuild all 32,862 cycles.  For every cycle:

1. validate parent and atom points on P-256;
2. replay both group identities above;
3. for `P=(x,y)` and `P_pi=(-x,z)`, verify without division
   `x(A)+x(B)=b/x^2` and `x(A)x(B)=(x^2+a)^2/(4x^2)`;
4. count negation-folded duplicates inside A, inside B, across A/B, and
   against the parent base;
5. screen complete-set transports for scalars `{2,3,4,8,16}`, preserving all
   partial matches without promotion credit; and
6. emit the deduplicated atom union in the repository's wide factor-base
   format with a content-derived FB1 identity.

Do not preregister collision counts or the new FB1 identifier.

## Complete controls

Exhaust `F_7`, two involution cycles, all 400 projective nonzero log profiles,
and all `7^4=2,401` atom coefficient tuples: 960,400 states.  Require direct
and factorized outputs to agree, inverse output transport to recover both atom
sums, augmented row rank two in both row orders, and zero discrepancies,
false positives, and false negatives.

On P-256, use a deterministic scalar-labelled target only to replay
`X+Y=Q` and `X-Y=G`.  It receives no actual-relation credit.

## Cost accounting

For each atom channel of width `N`, split sorted indices into deterministic
halves and use exactly ten signed distinct atoms from each half:

```text
L0=2^10*C(floor(N/2),10),  L1=2^10*C(ceil(N/2),10).
```

Price two targets as `2*(L0+L1)` entries, 82 bytes each, in `S` and relative
to rho.  This target-only diagnostic assumes atom logs and is not S17: two
arity-20 atom equations can contain up to 80 original point terms.

For 138,031 usable rows, allocate 69,016/69,015 between channels and import
Round 294's exact negation-folded collision law

```text
T_K=sqrt(2*n)*Gamma(K+1/2)/Gamma(K).
```

Report both an unrealistically free paired sample and the charge for two group
sums per sample.  Using `T_69016` for both channels omits the extra stopping
time for both quotas, so it is a favourable lower projection, never an
exhaustive run.  Add `2*K*nnz+K^2` sparse-algebra operations and keep setup,
collection, target decomposition, verification, and recovery separate.

## Hypotheses and stop condition

- H1: atomization factorizes the target exactly with zero replay errors.
- H2: actual atom collision or transport may reduce independent log dimension.
- H3: the target-only stage may approach rho while complete collection does not.
- H4: injective atoms plus a complete floor above rho are an accounting or
  representation change, not an end-to-end advance.

Commit the Rust runner before execution.  Preserve deterministic JSON, the
complete atom dump, isolation receipts, exact commands, tests, a requirement
ledger, one-unit boundary table, dashboard update, and a stacked PR.  Attempt
no unplanted full-depth relation unless every promotion gate passes.
