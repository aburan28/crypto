# Protocol: fiber-aware relation generation, and index calculus across a large conductor gap

Frozen 2026-10-10, before any run.  Two bounded experiments on public
synthetic curves.  Result class for both: **stage diagnostic and boundary
control**; `S`, end-to-end cost and speedup are reported from `ecbench`
sessions and classified by AGENTS.md §3, never read as a speedup unless
the whole-pipeline ratio to the matched rho moves.  No deployed-curve
claim; no key recovery.

Companion protocols this one instruments: I-4 (index-calculus yield
across a class) and I-5 (endomorphism ring, vertical edges) in
[`../isogeny_class_difficulty_20261008/`](../isogeny_class_difficulty_20261008/),
both frozen 2026-10-08 and unrun until now.

## A. Fiber-aware relation generation

### A.1 What "fiber" means here

The IC pipelines of this repository build the factor base from *raw*
curve points `P` with `x(P)` in a subspace (binary) or interval (prime),
and write every relation over the cofactor projection `[h]P` (the
relation loop multiplies the target coefficients by `h`,
`ic_boundary::collect_and_solve_with_completion`).  A subgroup target `R`
therefore has a **cofactor fiber**: the `h` raw points `T + K`,
`K ∈ E(F_q)[h]`, with `[h](T + K) = [h]T`.  This is the fiber the n = 83
five-sum notes name ("attempt all four raw target fibers under
multiplication by four", cryptanalysis `pdp-scaling/five-sum-next-20260929`).

The existing oracles are **fiber-blind**: a target `T` is decomposed only
as `Σ P_i = T` exactly.  A decomposition `Σ P_i = T + K` for any `K` in
the fiber is an equally valid relation row, because `[h]K = O`.

### A.2 Variants (the arms)

| arm | oracle | what changes |
|:--|:--|:--|
| A0 blind | `mitm`, `mitm-frobenius`, `descent-algebraic` as landed | one lift per target |
| A1 lifts | `…:fiber=lifts` | every lift `T + K` is probed / descended: `h` oracle calls per target |
| A2 closed | `mitm-fiber:fiber=closed`, `mitm-frobenius-fiber:fiber=closed` | the pair table holds `P_i + P_j + K` for every `K`: one probe per target against a table `h` times larger |
| A3 combined | `descent-algebraic:fiber=combined` | one Weil-descended system per target with the target abscissa `X` free (`n` more unknowns) and the fiber polynomial `g_T(X) = Π_K (X + x(T + K))` adjoined |

The rational `h`-torsion is computed by the oracle from `[r]P` over base
points, charged to oracle set-up; the fiber size it found is recorded
beside `h` (they agree when `gcd(h, r) = 1`, which holds on every curve
here).

### A.3 Derivation, stated before measuring

Let `B = |F|` (raw points), `N = #E`, `m = 2`.  Per target, the blind
MITM probes once and hits with probability about `B²/N`; the closed
table hits with probability `h·B²/N` at one probe; the lifts variant hits
with the same `h·B²/N` but at `h` probes and `h − 1` extra additions.
Hence:

- **P-A2.**  Relations per relation-phase group operation rise by a
  factor `h` (closed) against blind; the pair table costs `h` times the
  additions and memory.  At fixed `B` the relation phase shrinks by `h`
  and oracle set-up grows by `h`; re-optimising `B ∝ (N/h²)^{1/3}` the
  cold total falls by about `h^{1/3}` in the crude model.
- **P-A1.**  Relations per group operation are unchanged against blind
  (the `h` lifts cost `h` probes); relations per *target* rise by `h`.
- **P-A3.**  The combined system has `m·n' + n` unknowns and `2n`
  equations; `S₃` acquires a degree-3 boolean term (`x₁x₂X`), `g_T`
  descends to degree ≤ 2 for `h = 4`.  Prediction left open: the solving
  degree of F4 rises by at most one, and the solver cost per relation is
  compared with `h` blind solves.  Either outcome is a result.
- **P-A0 control.**  On a prime-order curve (`h = 1`) every variant
  reduces to blind and all counters agree.

### A.4 Frozen inputs

| item | value |
|:--|:--|
| Koblitz | `K_0 / F_{2^17}`, `K_0 / F_{2^19}`, `K_0 / F_{2^23}` (`h = 4`); `K_1 / F_{2^19}` (`h = 2`) |
| binary random | `binary_random` n = 17, 19 with seeds 1..4, `max_cofactor = 8`; cofactors as found |
| prime | `prime_search` 20-bit seed 59297 (`h = 1` control); the CM crater and floor curves of part B (`h = ℓ²·c`) |
| factor bases | Koblitz: `koblitz-orbit:divisor=0;1` (MITM) and `binary-subspace:dimension=8` (algebraic); binary random: `binary-subspace:dimension=8`; prime: `prime-abscissa:size=32` |
| oracles | `m = 2` throughout; solver `f4-f2` for the algebraic arms |
| targets per curve | 6; `target_seed = 20261010`; rounds 3, warmup 1, `alternate` |
| unit | ecbench group operations (`S`), relation-phase ops per relation, solver monomial operations, solving degree; wall only as a practicality note |

### A.5 Pass / fail and stop

- Q-A1: every arm verifies every target (`[d]G = Q`); any wrong answer
  halts the arm.
- Q-A2: A2's relation-phase ops per relation is below A0's by a factor
  within [0.5h, 2h] on every binary curve with `h = 4`; outside that
  band the model in A.3 is wrong and is corrected before any reading.
- Q-A3: A3's solving-degree max is reported next to A0's; no pass line,
  the number is the result.
- Stop: the arms above, one session per family, replays verified.
  Inadmissible: changing the base dimension, `m`, or trial caps between
  arms; reporting a relation-phase ratio as a speedup.

## B. The same index calculus across a large conductor gap

### B.1 Pairs

Two `F_q`-isogenous curves whose endomorphism-ring conductors differ by
a large factor: the **crater** (`End = O_K`) and the **floor**
(`End = Z[π]`) of an `ℓ`-volcano, joined by an explicit, kernel-certified
chain of `ℓ`-isogenies of total degree `ℓ^h` (the gap).

| family | construction | gap | certification |
|:--|:--|:--|:--|
| prime | CM by `D_K ∈ {−7, −8, −11, −19, −3, −4}`, `f_π = ℓ^h`, trace chosen so that `E[ℓ]` is rational above the floor; descend by rational Vélu kernels, excluding every `j` already visited | `ℓ^h ∈ {3^6, 5^4, 7^3, 2^8, 3^4, 11^2, 13^2}` at `p ≈ 2^22`, `2^26` | crater: `j = j(O_K)` (class number one); floor: cyclic `ℓ`-part of `E(F_p)`; levels: rank-2 `E[ℓ]`; chain: kernel points of order `ℓ`, `#E′ = #E`, image of the generator of order `r`, planted scalar transported |
| Koblitz | `K_0 / F_{2^n}`, the crater of its class; descend by one `ℓ`-isogeny for the prime `ℓ = f_π` (`n` prime) or by every `ℓ ∣ f_π` (`n` even) with kernels built in `F_{q^m}`, `m = ord_ℓ(t/2)` (`binary_torsion_walk`) | `271` (n=17), `457` (n=19), `967` (n=23), `627 = 3·11·19` (n=20), `7917` (n=28) | floor: rank-1 `E[ℓ]` over `F_{q^m}`; kernel trace and image abscissae rational; `#E′ = #E`; image of the generator of order `r` |
| binary random | random `b` over `F_{2^n}`, `n ∈ {17, 19, 23}`, classes whose `f_π` has a prime factor `ℓ ≥ 11`; the vertex found is placed by the rank of `E[ℓ]` over `F_{q^m}` and the other end of the height-one volcano is reached by one descending or ascending edge | the largest `ℓ ∥ f_π` found per `n` | as Koblitz |
| extension `F_{p^k}` | **not attempted** in this round (no explicit-curve constructor in `ecbench`) | — | — |

### B.2 Arms

Every pair is one `ecbench` spec with both curves as explicit workloads
and identical arms on both:

- reference: `rho.negation` (and `rho.signed_frobenius_strong` on the
  Koblitz crater only, where it applies);
- IC: prime `prime-abscissa:size=32` + `mitm:negation_folded=1`;
  binary `binary-subspace:dimension=8` + `descent-algebraic:m=2` with
  `f4-f2` (solving degree) and `mitm` (yield);
- Koblitz crater only: `koblitz-orbit:divisor=0;1` + `mitm-frobenius:m=2`
  (the τ-orbit fold, which the floor curve cannot run: it is not defined
  over `F_2`);
- control: `rho.negation` repeated (A/A).

### B.3 Metrics and predictions

Per curve: relation yield `γ = relations/trials` (hit rate), solving
degree mean and max (algebraic arms), `S` cold, online wall, floor
ratio.  Per pair: `γ_floor/γ_crater`, `S_floor/S_crater` with the paired
bootstrap interval.

- **P-B1 (I-4 null).**  `γ` and `S_IC` agree between crater and floor
  within the A/A interval for the same factor-base family and oracle;
  the Macaulay/F4 solving degree is identical (binary `S₃` sees `b` only
  as a constant; prime `S₃` is model-independent up to `|F|`).
- **P-B2 (I-5).**  Rho `S` is flat across the gap.  The only structural
  difference is the endomorphism available for folding: the τ-orbit
  (Koblitz) and the `j ∈ {0, 1728}` automorphism folds (prime `D_K = −3,
  −4`) exist on the crater and not on the floor; their effect is the
  known fold factor, not a new one.
- Kill line: a `γ` or `S_IC` ratio outside the A/A interval on a pair
  with `j ∉ {0, 1728}` and the same oracle, surviving the isomorphic-model
  control, is a reproducible unexplained anomaly (I-4 decision rule).

### B.4 Stop and inadmissible moves

Bounded: the pairs the constructor certifies at the sizes above; one
session per family; `ecbench verify --replay`.  Inadmissible: comparing
curves with different factor-base dimensions; reading a stage ratio as a
speedup; reporting a curve without its certificate; a claim about any
registered target curve.

## Amendments (recorded before the affected sessions ran)

1. **Binary factor-base dimension** (A.4, B.2): `dimension = 8` was
   frozen for every `n`; the sessions use 8 at `n = 17`, 9 at `n = 19`
   and 11 at `n = 23`, so that an algebraic run costs about
   `2^{n − dim}` solver calls (columns `≈ 2^{dim−1}`, hit rate
   `≈ 2^{2·dim−n−1}`).  Within one spec every arm and both curves of a
   pair use the same dimension, which is the comparison the protocol
   protects.
2. **A3 at `n = 17` only.**  The smoke run of the combined descent on
   `K_1 / F_{2^17}` took 67 s per target against 0.35 s for the blind
   descent (solving degree 3.0 against 1.95; 240× the solver
   operations), so A3 runs in its own spec (`fiber_combined_n17.json`,
   two targets, two rounds) on `K_1 / F_{2^17}` and two random curves
   over `F_{2^17}`, and not at `n = 19, 23`.  Its result at `n = 17` is
   the result.
3. **Targets and rounds.**  Four targets per curve and three rounds
   (one warm-up) for every spec but A3, instead of the six targets of
   A.4, so that the prime sessions (114 curves) finish in the session.
4. **`D_K = −7` chains.**  `(t² − f²D)/4` with `t, f` odd is even when
   `D ≡ 1 (mod 8)`, so for `D_K = −7` the conductor is `2ℓ^h`, the chain
   is the `ℓ`-descent followed by one 2-descent, and the gap is `2ℓ^h`.
   `ℓ = 2, h = 8` fails for `D_K = −7` at both sizes (no prime pair in
   the window) and is recorded as not constructed.
5. **Binary random pairs.**  The first scan re-found the Koblitz class
   (random `b` with the Koblitz trace); those three pairs are kept as
   `pairs/binary_pairs_koblitz_class_rediscovered.json` and not run.
   The scan excludes `|t| = |t(K_a)|`; the largest prime conductor
   factors found among 600 curves per degree are 17, 41 and 23.  A
   longer scan for `ℓ ≥ 60` (3000 curves per degree, `m·n ≤ 5000`)
   found none at `n = 17, 19`: `pairs/binary_pairs_large_ell.json`.
6. **Prime part B runs the fiber-closed MITM, not the blind one.**  The
   first `prime_gap_22` session (kept, `sessions/prime_gap_22`, status
   `interrupted`) exhausted its 100-million-trial cap on the cofactor-784
   pair without a relation: a fiber-blind pair table of `|F|²/2 ≈ 512`
   sums holds about `512/h` sums inside the order-`r` subgroup, below the
   33 independent rows the elimination needs once `h > 16`.  Every CM
   pair has `h ≥ 4` and most have `h > 50`, so the IC arm across the gap
   is `mitm-fiber:fiber=closed` on both curves (`prime_gap_22_fiber`,
   `prime_gap_26_fiber`).  The blind oracle's failure on large cofactors
   is itself a part-A result and is measured, capped at two million
   trials, in `fiber_prime_large_h`.
