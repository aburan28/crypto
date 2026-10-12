# Additive screen: Frobenius closure of the K0 terminal base

Registered before the orbit-closure run on 2026-10-04. This extends the
same high-arity factor-base capacity PR without changing or replacing the
earlier mixed-base rows. The curve, subgroup, and cofactor-projected
dimension-12 base are frozen by `PROTOCOL.md` and PR #1341.

Hypothesis: closing that base under all 83 Frobenius powers and sign gives
a much larger *actual usable* subgroup point set while preserving at most
the original 2,027 signed relation columns. The subgroup eigenvalue check
on the generator must pass, and each seed point must return after 83
Frobenius steps. Iteration uses `KoblitzCurve::frobenius` and exact point
negation; deduplicate encoded points, canonically choose the smallest
encoded member of each signed Frobenius orbit, and record both the point
count and orbit count. The source uses the fixed 33-byte affine encoding
defined in `RESULT.md`. A cofactor image commutes with Frobenius, so the
closure of projected points is exactly the projected closure of the
geometric seed base after removing identity.

Compute exact `C(B+3,4)` and `C(B+4,5)` numerators for the resulting
actual point count `B`, divided by the pinned prime subgroup order.
These are uniform-target coverage *ceilings*, never measured yields.
Admit a five-summand *capacity* proposal only if all checks pass, the
projected closure has no more than 2,027 signed Frobenius columns, and
its five-summand ceiling reaches 1%. Preserve the raw row and errors.

The available Boolean solver does not encode the union-of-charts law at
n83. There is no ordinary relation, candidate manifest, one-target IC
measurement, or speedup claim from this screen. A separate solver must
handle chart choices and prove exactness or report budget exhaustion.
