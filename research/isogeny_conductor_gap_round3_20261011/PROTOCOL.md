# Protocol: round 3, the two deferred leads of the conductor-gap study

Frozen 2026-10-11, before any E4 system was solved and before the E5
sessions ran.  Rounds 1 and 2 are
[`../isogeny_conductor_gap_ic_20261010/`](../isogeny_conductor_gap_ic_20261010/README.md) and
[`../isogeny_conductor_gap_round2_20261011/`](../isogeny_conductor_gap_round2_20261011/README.md).
Their certified prime pairs are reused unchanged (`pairs/prime_pairs.json`
of round 1); every curve is already registered.

Result class:
- E4 is a **stage diagnostic** of the algebraic decomposition step: the
  F4 solving degree and work of one summation-polynomial system per target.
  It is not a whole-method cost, and no rho comparison is made or implied.
- E5 is **engineering**, a negative-result audit: round 1's combined Weil
  descent, which round 1 reported only from unaudited warm-up records.

## E4. Two-torsion-symmetrised summation polynomials, crater against floor

### What round 2 deferred, and what it reduces to

Round 2 deferred this lead for want of an `F_p` Gröbner solver.  The
repository has one: `cryptanalysis::f4_fp`, a degree-bounded F4 over
`F_p`, `p < 2³²`.  It reports the solving degree, the widest matrix up to
the solution, and the `F_p` multiplications of row reduction.

**Derivation, stated before measuring.**  Let `H ⊆ E(F_p)` be a rational
subgroup of 2-power order, so `H ⊆ E(F_p)[h]`.  Translation by `T ∈ H`
preserves a two-term decomposition: `(P₁ + T) + (P₂ − T) = R`.  If the
factor base `F` is `H`-stable, the decomposition variety is a union of
`H`-orbits.

The quotient by that action is the Vélu isogeny `ψ: E → E/H`.  A
decomposition `P₁ + P₂ = R` maps to `ψ(P₁) + ψ(P₂) = ψ(R)` over
`ψ(F)`.  Conversely, a decomposition on `E/H` lifts to
`P₁ + P₂ ∈ R + H`.  That is a valid row of the relation matrix, which
every pipeline here writes over `[h]P` (round 1, §A.1).

So the symmetrised system is Semaev's `S₃` on `E/H`.  Its factor-base
polynomial has degree `|F|/|H|` instead of `|F|`, in the same two unknowns.
This is the same quotient as round 2's partial fiber, read algebraically.

**Consequences.**
- **P-E4a (gain).**  At a matched factor base of `B` abscissae on `E`,
  the quotient system's solving degree falls about `|H|`-fold against the
  direct system.  Its row-reduction work falls by at least `|H|²`.  Both
  are monotone in `|H|`.
- **P-E4b (gap).**  At matched `B` and matched `|H|`, crater, floor and
  isomorphic floor model have the same solving degree and widest matrix,
  and row-reduction work within 10 percent.
  - The crater of an ℓ = 2 volcano can have rational `E[2] ≅ (Z/2)²`,
    which its floor (rank-1 2-torsion) cannot.
  - But wherever 4 divides the cofactor, the floor has a cyclic `C4` of the
    same order.
  - So the crater's richer 2-torsion changes which `H` exist, not what a
    given `|H|` costs.
- **What would refute P-E4b:** `E[2]` on a rank-2 crater solving at a
  lower degree, or with ≥ 25 percent less work, than `C4` on its floor or
  the floor's model at the same `B`.  That would be a structure effect
  beyond the order of `H`.
- **P-E4c (correctness).**  On the planted target, the planted solution
  annihilates the basis in every run.  Inconsistency is decided the same
  way by F4 and by direct enumeration of the base pairs.  The Vélu map is
  checked on 32 random sums per subgroup before any system is built.  Any
  failure stops the lead.

### Inputs

- Every round-1 prime pair with even cofactor (24 pairs, 22 and 26 bits),
  each with crater, floor and isomorphic floor model.
- `B ∈ {32, 64}` abscissae, the smallest on the curve, taken with their
  `H`-orbits so the base is `H`-stable.
- `H ∈ {1 (direct), C2, C4, E[2]}`, as the curve has them.  Where the curve
  has several, `C2` is its first rational 2-torsion point and `C4` its first
  element of order 4.
- Three targets per (curve, `B`, `H`): two pseudo-random points and one
  planted sum of two base points.
- F4 under grevlex, degree bound `2·deg + 8`, 20 s budget per system.
- The runner is `examples/isogeny_gap_two_torsion_f4.rs`, appending one
  JSON line per system to `e4/records.jsonl`.  It can restart
  curve by curve after a host reboot.
- Field-operation counts and degrees are exact.  Wall times are host
  observations, not admitted measurements.

## E5. Round 1's combined Weil descent, audited

Round 1's arm A3 (`descent-algebraic-fiber:fiber=combined`, the target
abscissa free and the fiber polynomial adjoined) was cut short by host
reboots four times.  It was reported from warm-up records:
- verified in 335 to 578 s on the a = 1 Koblitz curve over F_2^17;
- 629 to 654× the blind F4 work;
- maximum degree 4.

This host still reboots every 20 minutes to 2 hours, so the rerun is cut to
fit.

- **Specs.**  Round 1's arms unchanged: rho, blind F4, F4 over all lifts,
  combined F4, all on `binary-subspace:dimension=8`.
  - `e5_combined_n17_t0`: one target, no warm-up, two measured rounds.
  - `e5_combined_n17_t1`: the same with another target seed.
- **What is predicted.**  The audited records reproduce round 1's warm-up
  figures:
  - the combined arm verifies;
  - its F4 work is within 600 to 700× the blind arm's;
  - its maximum solving degree is 4.

  Solver work is deterministic, so a repeat of round 1's per-target counts
  is expected exactly where the same target recurs.
- **Stop.**  A session cut by a reboot is kept under `interrupted/` with
  its log and rerun once.  A second cut is reported and not rerun.

## Visual record

Each run gets a run report, a graph and a PDF under `reports/run-<UTC>/`
(the `research-visuals` skill).  The round's README carries the results
graphs and a PDF.
