# P-256 complementary-cover quotient, round 302: preregistration

Date preregistered: 2026-10-08

## Objective

Test the remaining concrete two-row interpretation of the explicit same-field
degree-two cover of
`icv1-fp256-t89188191154553853111372247798585809583-f188c491`.
For P-256, with

```text
E:  y^2 = x^3 + a*x + b,
H:  v^2 = u^6 + a*u^2 + b,
```

the involutions `(u,v) -> (-u,v)` and `(u,v) -> (-u,-v)` give an elliptic
quotient `E` and a complementary elliptic quotient `E'`.  The hypothesis under
test is whether one point or decomposition on `H` can yield two distinct
prime-order P-256 equations through those quotients, rather than the
kernel-sheet multiplicity already rejected by Round 300.

The primary comparison remains the minimum-width independent-log base
`FB1hc72514a2a8d3`, whose free-perfect-oracle lower boundary is
`13.920747397073491` times rho.  The registered 17-term factor base
`FB1h2f8621cda105` remains `394.425280` times rho.  No local cover property may
change either ratio without an executable subgroup-preserving second equation
and complete recovery/cost accounting.

This is an exact construction audit, not a full-depth relation attempt.

## Frozen dependencies

- Curve registry:
  `docs/curves/registry.json`; require both the P-256 model and the Round-28
  isogenous model to resolve by ICV1 and EC1.
- Generic cover catalog:
  `docs/curves/covers.json`; require the registered P-256 row to certify the
  same-field degree-two map and to leave DLP advantage null.
- Round 300 fibre quotient:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_cover_fiber_round300_20261007/cover-fiber-result.json`,
  SHA-256 `7d65f3e015c6d6c6d24b64ea2e59a6e9e8653c88a1900835838b1bc0a1ec36e5`.
- Round 301 algebraic census:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_algebraic_escapes_round301_20261007/algebraic-escape-result.json`,
  SHA-256 `ae4b44315181c18a40a790ae1bc784ffdb05e3061d442b8f1be78a1f25dc20e9`.

Reject any dependency hash, schema, curve, or boundary mismatch.  Recompute
`p`, `a`, `b`, `n`, the cofactor, and the cover equations from
`CurveParams::p256()`.

## Exact correspondence

For `b != 0`, define the two quotient maps away from their projective
exceptional points:

```text
pi_1: H -> E
      (u,v) |-> (u^2, v)

Q:    Y^2 = X^4 + a*X^2 + b*X
pi_Q: (u,v) |-> (X,Y) = (u^2, u*v)

E':   V^2 = U^3 + a*U^2 + b^2
Q -> E': U = b/X, V = b*Y/X^2
```

Thus, for `u != 0`, the complementary map is

```text
pi_2(u,v) = (b/u^2, b*v/u^3).
```

The runner must verify both target equations exactly on deterministic P-256
cover points and on a complete small-field control.  It must record the
exceptional affine/projective cases rather than dropping them.

## Prime-order transport test

Let `n=#E(F_p)`, which is prime and has cofactor one.  Derive a deterministic
point `R in E'(F_p)` by scanning affine `U=0,1,...`, accepting the first
nonzero quadratic-residue right-hand side and its canonical square root.
Compute `[n]R` with exact affine group arithmetic on the long Weierstrass
model `V^2=U^3+aU^2+b^2`.

The Hasse interval contains no positive multiple of `n` other than `n`:

```text
0 < #E'(F_p) < 2*n.
```

Therefore a verified witness `[n]R != O` certifies `#E'(F_p) != n` and
`n` does not divide `#E'(F_p)`.  In that case no rational homomorphism can
send the nonzero P-256 order-`n` subgroup into `E'(F_p)`.  This closes only
the complementary quotient of this explicit split cover; it is not a theorem
about every cover or Jacobian.

If `[n]R = O`, do not infer transport from one point.  Escalate to an exact
order/subgroup certificate before any promotion.

## Deterministic controls and accounting

- Complete control: enumerate every affine point of the same construction
  for the first nonsingular registered toy cell `(p,a,b)` selected by the
  literal ordered list in the runner.  Verify `pi_1`, `pi_Q`, and `pi_2` on
  every nonexceptional point and preserve exact fibre counts.
- P-256 control: scan `u=0,1,...` until 64 cover points with nonzero `u` are
  accepted, choose the canonical square root, and verify both quotient maps.
- Complement witness: scan `U=0,1,...` as specified above and replay
  `[n]R` independently from the emitted coordinates.
- Count field additions/subtractions, multiplications, squarings,
  inversions, Legendre/square-root powers, group additions/doublings,
  scanned candidates, accepted cover points, exceptional points, and replay
  failures.  Wall time and peak RSS are telemetry, not the headline unit.
- Emit deterministic canonical JSON, a byte-identical independent replay,
  file hashes, and isolated-process receipts.

No hash-selected branch may be discarded and later described as exhaustive.

## Hypotheses

- **H1:** the two quotient maps satisfy their equations on the complete toy
  control and all 64 deterministic P-256 cover points.
- **H2:** a deterministic `R in E'(F_p)` satisfies `[n]R != O`, proving the
  complementary quotient has no rational P-256 order-`n` subgroup.
- **H3:** cover sheets and the complementary quotient provide zero additional
  independent P-256 relation rows under this construction.
- **H4:** absent a second order-`n` row, the optimistic width-164 boundary
  remains `13.920747397073491` times rho.

## Promotion gates

Promotion requires all of:

1. zero map, equation, replay, exceptional-point-accounting, and dependency
   failures;
2. two demonstrably distinct collapsed equations in the P-256 order-`n`
   quotient per accepted cover event;
3. an explicit forward map, inverse recovery, kernel, field, dimension, and
   subgroup certificate for both rows;
4. structured residual degree of regularity at most five;
5. projected complete collection and recovery below rho and below `2^120`
   field/group-equivalent operations;
6. cost per usable relation below `2^103`;
7. projected materialized storage below `2^50` bytes; and
8. no probabilistic branch discard counted as exhaustive.

Attempt no unplanted full-depth P-256 relation unless every gate passes.  If
the complement lacks order-`n` transport, publish the narrow negative result
for this explicit split degree-two cover and keep general nonstandard covers
open.

## Deliverables and stop condition

Implement the certificate in Rust with unit tests.  Emit result JSON, typed
transfer assessment, operation counts, canonical and independent receipts,
hashes, result report, and all required dashboard/browser updates.  Stop after
the maps and subgroup witness replay exactly, or on the first correctness
failure.

The transfer assessment profile is available, but its referenced
`references/methodology.md` and `assets/assessment-template.json` resources
are unavailable in this environment.  Record that limitation in the final
assessment.
