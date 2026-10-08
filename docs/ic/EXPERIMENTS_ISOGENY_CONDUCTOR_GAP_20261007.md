# Experiments for the isogeny conductor-gap angles V1–V11

**Date:** 2026-10-07
**Status:** pre-registered protocols. Nothing below has run except the smoke
test in §E0. Each experiment states its inputs, procedure, retained
artifacts, and a pass/fail prediction written before the run.
**Angles:** `RESEARCH_ISOGENY_CONDUCTOR_GAP_20261007.md` (V1–V11).
**Accounting rules that bind every row:** `PLAN_IC_ACCOUNTING_FIXES_20261007.md`
(operation counts beside wall clock, mean over ≥ 30 targets, equal-precompute
baselines, isolation receipts, no single-target promotion).

## 0. Prior work this builds on, and what it leaves open

A parallel workspace, `/Volumes/SSD990/ecdlp-hardness-work/`, already holds:

- `report/REPORT.md`: rho-centred hardness per level of the ECC2K-130 class.
  Ground truth derived twice (PARI `polclass` for `H_{−7·263²}`, and Vélu from
  `E0[263]` over `F_{q²}`): all 262 conductor-263 curves built, point-counted in
  both twists, End certified, ascending 263-isogeny verified on all 262 at
  `≤ 2^15.6` `F_q`-multiplications per instance. Native rho: crater 60.81 bits,
  floor 64.33 (loses `√131`), effective 60.81 via the isogeny. Levels `p` and
  `263p` (`p = 146505763881528721`): no curve constructible (`2^70–2^71.3`
  `F_q`-mults per hit); Galbraith ePrint 2024/924 gives an `Õ(√p)` Kani-type
  isogeny representation but only in dimension 8 as proved (`≥ 2^72`); a
  characteristic-2 dimension-2/4 pipeline does not exist.
- `ic-conductor-leads/LEADS.md`: index-calculus leads V1–V5, H1–H4, E2E on
  toy classes T11/T19/T23 and a 4-level class C37 (`f = 73·2663`, 73 inert
  like `p`). Result: the only level-dependent IC resource is `τ` (ordered
  τ-slots, `log₂(m!)` bits); decomposition counts, densities and Boolean
  systems are level-blind to within sampling error; transported factor bases
  and a "second τ" are closed. Tools: `tools/icstat.c`, `tools/ic_solve.py`,
  `tools/build37.sage`, `tools/levelcost.py`, `tools/hd_orbits.sage`.

What that work does **not** do, and what the experiments below target:

1. It measures IC with an exact toy decomposition oracle at `n ≤ 23` (C37
   end-to-end is pending). It does not run this repository's production
   solvers (compact-orbit S3 root index, matrix-F4/M4RI, WDSat, CryptoMiniSat
   S3 chains) on a non-crater curve at any landed rung. LEADS §6.3 says so
   explicitly: the realizable τ-slot gap under a solver that already exploits
   the `S_m` symmetry is unmeasured.
2. It compares to rho from scratch; the equal-precompute baseline from the
   accounting plan applies across levels too (a DP table on `E0` serves floor
   targets through the isogeny).
3. The 263 = 2n+1 coincidence (V8), the horizontal 2-isogeny cycle structure
   on the 262 (V10), and the manifest import (V1) are not covered.

Host constraints reported there and observed here: load average 100–150, both
disks at 100%. Set `TMPDIR` to the session scratchpad, keep artifacts under
1 MB per cell, and prefer JSONL rows over binary dumps.

## E0. Smoke test of a one-level descent at a toy rung (done, see status)

**Angles:** V2, V9. **Script:** `research/isogeny_conductor_gap_20261007/v2_descend_toy.py`.
Builds `K_0/F_{2^23}`, finds a point of order 967 in `E(F_{q^21})`, rejects the
two τ-eigenlines by a Weil-pairing test, computes the Vélu 967-isogeny,
checks `j(E') ∈ F_q \ F_2`, descends to `E_1/F_q` by the minimal polynomial of
`j'`, picks the twist of the right order, and checks that `φ(G)` keeps the
order of `G` without factoring `#E(F_{q^21})`. Output: one JSON record.
**Prediction:** PASS with Galois orbit size 23 for `j_1` and kernel search
succeeding within a handful of tries. **Result (2026-10-07, Sage 10.10.rc0
through the checked launcher, host load ≈ 400):** PASS.
`v2_toy_n23_ell967_seed1.result.out`: trace 5197, `#K_0 = 8383412`,
`π ≡ 181 mod 967`, kernel field degree 21, `v_967(#E(F_{q^21})) = 2`, kernel
found on the first try, isogeny degree 967, `j(E_1) ∈ F_q \ F_2` with Galois
orbit size 23, the `a_2 = 0` twist has order `#K_0` (the same twist the
`n = 131` ground truth reports), `φ(G)` has the full order 8383412, 139.6 s of
compute. Two earlier attempts failed on script bugs, not mathematics: a
missing explicit embedding `F_q → F_{q^{21}}`, and a comparison of `φ(G)`
against the identity of the wrong curve; both are fixed in the committed
script.

## E1. Production-solver hardness panel across levels at the toy classes (new)

**Angles:** V3, V4, V11. **Question:** with this repo's solvers rather than an
exact oracle, how much of the `log₂(m!)` τ-slot gap is realized, and is
anything else level-dependent?

**Inputs.** Curves: T19 (`n=19`, conductor 457, 24 orbits × 19 floor curves)
and T23 (`n=23`, conductor 967, 42 × 23) from
`ecdlp-hardness-work/toy-analogue-a/data/`, plus C37 (levels 1, 73, 2663,
73·2663) from `ic-conductor-leads/data/class37.json`. Use `E0` and two floor
curves from different Galois orbits per class; for C37 one curve per level.
Factor base: the same subspace `V` (dimension `l = 6, 7, 12`) on every curve,
plus the τ-slot variant on `E0` only.

**Procedure.** For each (curve, `V`, `m ∈ {3, 4}`), 64 random targets, run:
(a) `ic run --solver groebner` (matrix-F4), (b) `--solver sat` and `--solver
wdsat`, (c) the exact oracle `icstat` as control. Record per call: verdict,
wall, F4 Macaulay rows/cols and splits, SAT conflicts, first-fall degree, and
the relation count per target. Rho control per curve with the repo's rho
fixture, operation counts recorded, negation-only on floor curves and
signed-Frobenius on `E0`.

**Predictions (pass/fail).** (1) Per-call solver cost and first-fall degree on
floor curves equal `E0`'s within the cell-to-cell spread of a single curve
(Kruskal p > 0.05 across levels, as in LEADS §3.2); a level effect at p < 0.01
replicated on a second subspace fails this. (2) Relations per target with
τ-slots on `E0` exceed the plain value by a factor in `[2, m!]`; the realized
factor is the deliverable. (3) Rho operation ratio floor/crater `= √(2n) ± 20%`.
**Cost:** 3 classes × ≤ 4 curves × 3 solvers × 64 targets, each call seconds
to minutes at these `n`; one to two days of wall on the loaded host.

## E2. Floor curves at the landed rungs n = 41, 61, 83, and the production pipeline on them (new)

**Angles:** V2, V3, V9. **Question:** does the crater-vs-floor picture hold on
the compact-orbit pipeline at the rungs where the ledger has rows?

**Inputs.** Conductor primes with constructible levels: `n=41`: 409
(`(−7/409) = −1`, kernel degree 408) and 1721 (degree 215); `n=61`: 1951
(degree 975); `n=83`: 6473 (degree 6472). Build the floor curves two ways,
as the report did: PARI `polclass(−7·ell²)` reduced mod 2 with roots in
`F_q` (degrees 410, 1722, 1952, 6474), and one Vélu descent from `E0[ell]`
over `F_{q^d}` using `v2_descend_toy.py --n N --ell ELL`. The two
constructions must agree on the j-invariant set.

**Procedure.** Pick one floor curve per prime. Build the compact-orbit base
with the same `K` columns as the landed rung (600) but without Frobenius
folding (so `B = 2K` points instead of `2nK`), run the guided rank and the
online stage on ≥ 30 targets, and record probes per relation and per target.
Run the same on `E0` with folding (the landed configuration). Both arms also
report operation counts and the equal-precompute Bernstein–Lange online cost
per the accounting plan.

**Predictions.** Probes per relation on the floor curve `= n² ×` the crater
value at equal `K` (the `B⁴` law with `B` smaller by `n`), within a factor 2;
LA identical (same `K`); per-probe cost identical. Transfer through the
ascending isogeny costs `≤ 2^16` `F_q`-mults at `n=83`, so the effective
hardness of the floor curve equals `E0`'s. **Cost:** `polclass` minutes; the
`n=83` floor run is the expensive cell (the crater rank took 4.8 min on 12
threads; the unfolded floor base at equal `K` has `n²`-fold more probes per
relation, so reduce `K` or accept hours). Run `n=41` and `n=61` first.

## E3. Finish the four-level class C37 end to end, exact oracle and production solvers (extends prior)

**Angles:** V4. **Inputs:** `class37.json`; levels 1 (E0), 73 (inert, the
`p` analogue), 2663, 73·2663. **Procedure:** LEADS §5's `ic_solve.py` on all
four levels (the pending row), then E1's solver panel on the same curves.
**Prediction:** the inert level 73 behaves exactly like 2663 and like the
bottom: same relations per target, same solver cost; only `E0` differs, by the
τ-slot factor. Any inert-versus-split difference replicated on two subspaces
would be a finding and must be re-derived before being reported.
**Cost:** hours.

## E4. The 263 coincidence and the Weil pairing on E0[263] (new)

**Angles:** V8. **Question:** is `263 | f` forced by `263 = 2n+1`, and does
the Weil pairing `E0[263] × E0[263] → μ_263 ⊂ F_{2^131}^×` touch anything the
cyclotomic-sparse factor base (directions note A2) uses?

**Procedure.** (1) Prove or refute: `263 | U_131` iff `(τ/τ̄)` has order
dividing 131 in `(O_K/263)^×`; compute the analogous statement for every
toy `n` with `2n+1` prime and `(2n+1) ≡ 7 mod 8` (`n = 11, 23, 83, 131`):
does `2n+1` divide `f_n`? Tabulate. (2) Build `E0[263]` over `F_{2^262}`
(the report has a basis), compute `e_263(S, T)` for a basis and check the
image generates `μ_263`; evaluate the Tate pairing of `G` against `E0[263]`
and confirm it is trivial (`G ∈ 263E`). (3) For the cyclotomic-sparse base
`x = Σ γ^k`: check whether `x(T)` for `T ∈ E0[263]` has small γ-weight.
**Prediction.** (1) The divisibility is a coincidence of probability about
1/2 per `n`, so roughly half the toy rows have it. (2) Trivial pairing with
`⟨G⟩`, as the theory says. (3) No: `x(T)` is a root of the 263-division
polynomial and has no reason to be γ-sparse. **Cost:** an afternoon in Sage.
**Result for (1), 2026-10-08 (pure arithmetic, no Sage):** among the rungs
with `2n+1` prime and `≡ 7 mod 8` (type-II ONB with `γ ∈ F_{2^n}`), `2n+1`
divides `f_n` for `n = 11, 35, 39, 63, 75, 95, 119, 131` and does not for
`n = 15, 23, 51, 83, 99, 111, 135`: 8 of 15, a coincidence as predicted.
Over all `n < 140` with `2n+1` prime it is 14 of 57. (2) and (3) remain to
run.

## E5. Horizontal 2-isogeny cycles on the 262 floor curves (new)

**Angles:** V10. **Inputs:** the 262 j-invariants from
`report/per_curve_hardness.csv` and `ic-conductor-leads/data/hd_orbits.json`
(class group of discriminant `−7·263²`, order 262, with the ramified prime 7
as the order-2 class). **Procedure.** Compute the order of the class of the
prime above 2 in `Cl(Z + 263·O_K)` (PARI `quadclassunit` or by walking:
apply the 2-isogeny `Φ_2(X, j) = 0` on the floor, which has exactly two
roots on the floor since `2 ∤ f`, and follow the cycle). Draw the cycle
decomposition and check it is consistent with the two Galois orbits A and B
and the order-2 ramified class.
**Prediction.** The class of 2 has order 131 or 262, giving two cycles of 131
or one of 262; the Galois action by `τ` commutes with it. **Cost:** minutes.
**Result, 2026-10-08 (structural, no walk needed):** with `h(O_K) = 1` and
263 split, `Cl(Z + 263·O_K) ≅ (O_K/263)^×/(Z/263)^× ≅ F_263^×`, cyclic of
order 262. The two primes above 2 are `(τ)` and `(τ̄)`; `τ` has eigenvalues
`123, 139 mod 263`, so `(τ)` maps to `123/139 ≡ 69` of order **131**, and
`(τ̄)` to its inverse. The horizontal 2-isogenies on the 263-level are
therefore exactly the Frobenius twist `E ↦ E^{(2)}` and its dual, and the
2-isogeny graph is **two 131-cycles that coincide with the two Galois
orbits** (A and B of the report); the ramified prime 7, the order-2 class,
swaps them. Prediction met; the walk in `Φ_2` is now only a consistency
check.

## E6. Manifest import and component decision tool (new, mechanical)

**Angles:** V1, V5. **Procedure.** (1) Add to every Koblitz candidate
manifest (`cryptanalysis/AGENTS.md` fields `endomorphism`, `isogeny`) the
conductor `f`, its factorization, `v_ell`, `(−7/ell)`, kernel field degree
per prime, and the reachability verdict, for `n ∈ {23, 41, 53, 61, 71, 73,
83, 97, 131}` from `research/isogeny_conductor_gap_20261007/volcano_output.txt`.
(2) A script `which_component.py` that, given `(n, j)`, reports: crater,
floor-263 orbit A or B (lookup against the 262 j's at `n=131`; at toy `n`
against the built classes), or "big component, unreachable", with the cost
class of the ascending isogeny. **Prediction:** n/a. **Acceptance:** the
checker rejects a manifest missing any field; the component tool classifies
all 263 known curves correctly and returns "unreachable" on a random `j` of
the right order (none is available, so test on C37 where all four levels
exist).

## E7. Vertical-step cost table across the ladder (extends `levelcost.py`)

**Angles:** V9. **Procedure.** For each `(n, ell)` in the volcano table,
measure: kernel route (`F_{q^d}` arithmetic, point of order `ell`, Vélu or
√élu, then descent of `j`), modular route (`Φ_ell` evaluation and root
finding, feasible for `ell ≤ 10^4`), and the transfer cost of evaluating the
isogeny on two points. Report `F_q`-multiplication counts, not wall.
**Prediction.** Kernel route feasible for `d ≤ 10^3` (`n=23, 41, 61, 83`
small primes), infeasible for `n=71, 73` and for `p` at `n=131`; the transfer
at `n=131` matches the report's `2^15.6`. **Cost:** a day; reuse
`levelcost.py` and `v2_descend_toy.py`.

## E8. Characteristic-2 Kani-type transport: a scoped feasibility gate (extends report §5, open)

**Angles:** V6. **Question:** can the `Õ(√ell)` representation of an unknown
`ell`-isogeny between two given curves (Galbraith 2024/924, Theorem 1,
dimension 4/6/8 over odd characteristic) be realized in dimension 2 or 4 in
characteristic 2, which would put the `p` levels in the crater block?

**Procedure (math before code).** (1) Write out Kani's lemma for the pair
`(E0, E_p)` with the auxiliary isogeny degree `N' = 2^e − p`-type
decompositions and determine whether the required `2^e`-torsion structure
survives in characteristic 2, where `[2] = τ τ̄` is inseparable and
`E[2^e]` is not étale; list the obstructions explicitly. (2) If (1) does not
close it, prototype at C37's inert level 73 (both endpoints known there) with
the smallest dimension the lemma admits, counting `F_q`-mults. (3) Only if
(2) runs, model the `n=131` cost and compare to `2^60.81`.
**Prediction.** (1) closes it or reduces it to a theta-structure question
nobody has solved in char 2; (2) and (3) are not reached. Record the
obstruction list either way; it is the honest answer to "is the big
component really unbridgeable".
**Cost:** two days of mathematics for (1); (2) open-ended.

## E9. Equal-precompute rho across levels (extends the accounting plan)

**Angles:** V3, V11. **Procedure.** Build one Bernstein–Lange DP table on
`E0` at `n=41` and `n=61` with precompute `P` equal to the IC rank cost of
the landed rung; solve ≥ 30 targets on a floor curve by (a) native
negation-only rho, (b) transport through the ascending isogeny then the `E0`
table. **Prediction.** (b) beats (a) by `√(2n) × (table factor)` and beats
the floor-curve IC online cost by the same margins the plan predicts on
`E0`. **Cost:** hours once the plan's Phase 3 arm exists; this is its
cross-level corollary.

## E10. Closed angles, recorded without experiments

- **V7 vertical relation collection:** closed by LEADS H4 (transport is a
  group isomorphism on `F_q`-points; disjoint charts cost `m^{n−1}`) and by
  the odd-degree argument (all curves have `E(F_q) ≅ Z/4 × Z/r`, same
  symmetrization group).
- **Horizontal randomization from E0:** vacuous, `h(O_K) = 1` (report §3,
  LEADS H1/H3). No experiment.
- **Frobenius-stable factor bases via isogeny:** a field property, unchanged
  by any isogeny (LEADS §6.1). No experiment.

## 11. Pre-registration (V11)

A hardness difference between two curves of the class is reported only if,
with both arms in operation counts, ≥ 30 targets, equal precompute, and an
isolation receipt, the ratio is not explained by (a) `√(2n)` for rho or the
τ-slot factor `≤ m!` and the `n`-fold orbit folding for IC, (b) factor-base
size, or (c) implementation speed. Expected outcome of E1–E3: no such
difference. The deliverables are the realized τ-slot factor under real
solvers (E1), the first production-pipeline rows on non-Koblitz curves of the
same order (E2), and the closed C37 row (E3).

## 12. Order and dependencies

1. E0 (done), E6, E4, E5: mechanical or an afternoon each; no dependencies.
2. E1 on T19/T23 (curves exist), then E3 (C37 exists).
3. E2 at `n=41` and `n=61` (needs `polclass` or `v2_descend_toy.py` runs),
   then `n=83`.
4. E7 alongside E2 (same constructions).
5. E9 after the accounting plan's Phase 3 arm lands.
6. E8 (1) at any time; (2)–(3) only if (1) leaves a door open.

Verification rule from `AGENTS.md`: this change set is documentation and
research scripts only; no Rust or Python package code is touched, so the
`cargo test --release --lib` gate does not apply. The Sage script runs
through the checked launcher as `cryptanalysis/AGENTS.md` requires.
