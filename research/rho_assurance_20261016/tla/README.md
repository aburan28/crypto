# TLA+ model of the strong distinguished-point rho

`StrongRho.tla` models the protocol of `StrongRhoG::solve` in
`src/cryptanalysis/koblitz_strong_rho.rs` as committed in `448be6146`
(git blob `c24377d08`; every line number below refers to that blob). That
commit changed the table hasher, the linear maps of `canonicalize` and added
`debug_assert!`s; the walk itself is the same as in `680f3cf9d`. Modelled:
lockstep lanes, walk starts from a strided cursor, signed-Frobenius
canonicalization, distinguished points, the
distinguished-point table, collisions, fruitless cycles in both modes, and
the walk and step caps. TLC checks it exhaustively on small instances.

Run everything with `./run_tlc.sh` (or `./run_tlc.sh a.cfg b.cfg` for a
subset). The script downloads `tla2tools.jar` (the v1.8.0 release URL) into
`/tmp/claude-0/tla/` if it is missing, checks its pinned SHA-256, runs TLC
with 2 workers on every `*.cfg`, writes each output to `out/<config>.txt`
and a table to `out/SUMMARY.txt`, and exits nonzero when any config misses
its expected verdict (the `\* EXPECT:` line at the top of each cfg). It
puts TLC's `-metadir` on `/dev/shm` because the shared disk filled up during
one development run.

## The model

The prime-order subgroup `<G>` becomes the additive group Z_R: `G = 1`,
`Q = d`, a point is a residue, `0` is the point at infinity, and `k·P` is
`k*P mod R`. The automorphism group is `{±LAM^k : 0 ≤ k < NF}` with `LAM` of
multiplicative order `NF` mod `R` (the Frobenius eigenvalue λ and the field
degree n of the Rust), so it has `A = 2·NF` elements.

| Rust | TLA+ | What the model keeps, what it replaces |
|---|---|---|
| `from_parts` 912–943 (assert 930): `lam_pow`, `assert!(cur == S::one())` | `LamPow`, `ASSUME Md(LamPow[NF-1]*LAM) = 1`, `Mults`, `A` | the multiplier table and the order check |
| `canonicalize` 1143–1171 (with `least_rotation` 678, `choice_exhaustive` 1111) | `CanonMul`, `Canon`, `CanonState` | canonical point = least element of the orbit `{m·x : m ∈ Mults}`; `a` and `b` are multiplied by the multiplier `m`; infinity is returned unchanged (multiplier 1). The normal-basis rotation that computes the least orbit member and its `±λ^k` is **not** modelled: the model takes the least residue directly. Whether rotation by `k` really corresponds to `λ^k` is the subject of the Rust unit test `canonical_form_is_orbit_invariant_and_multiplier_is_correct` (line 1567), not of this model. |
| `partition` 790, `dp_hash` 849, `dp_mask` 1211 | `Partition(x) = (x*C1 + C2) % J`, `IsDP(x) = (x*C3 + C4) % D = 0` | the 64-bit mixers are replaced by fixed affine functions of the canonical point; `D` plays `2^dp_bits` |
| `jumps` 1174–1186 (`JUMPS = 32`, `s` uniform in `[1, r)`) | `NJ` ∈ {2, 3} constants `J1..J3` | a fixed jump table, two or three nonzero scalars; six different tables across the configs |
| `solve` 1204–1207: `stride_a`, `stride`, `cursor_a`, `cursor = cursor_a·G + Q` | `Init`: `d`, `strideA`, `cursorA` chosen from ranges, `cursor = cursorA + d` | the seeded RNG draws become nondeterministic choices, so TLC explores every target and every start sequence in the ranges |
| `walk_cap = 8 << dp_bits` 1212, `step_cap` 1214–1216 | constants `WalkCap`, `StepCap` | small numbers instead of `8·2^dp_bits` and `≥ 10^6` |
| `table: HashMap` 1218, `lanes` 1220–1232 | `table` (function point ↦ ⟨a, b⟩), `lane` (record per lane: `pt, a, b, prev, len, live`, plus the ghost `nprev`) | `previous` is a 4-tuple initialised with infinity, exactly as in the Rust |
| `while steps < step_cap` 1243 | `Begin` (phase `begin`) | when `steps ≥ StepCap` the state becomes `capped` (`None`) |
| phase 1, start loop 1244–1267 | `StartLane(l)` (one turn of the `loop`), `EndStart` | takes `cursor`, advances `cursor += stride`, `cursor_a += stride_a`, `walks += 1`; if the taken point is infinity the same lane loops again (`continue`, 1253); otherwise the lane gets `canonicalize(start, start_a, b = 1)`, `previous = [∞; 4]`, `length = 0`, `live` |
| phase 2, distinguished points 1268–1332 | `DpLane(l)`, `EndDp` | lane leaves (`live = false`); table miss → insert (1277–1289); hit with `b = b'` (no inverse) → `wasted` (1290–1294); otherwise `candidate = (hit_a − a)·(b − hit_b)^{-1}` (1295–1296), returned as `solved` if `candidate·G = Q` (1298), else `failed_collisions += 1` (1331) |
| phase 3, one step per live lane 1333–1405 | `StepLane(l)`, `EndStep` | `next = canonicalize(P + s_j·G, a + s_j, b)`; `steps += 1`, `length += 1` (1370–1371); fruitless if `next == current ∨ next ∈ previous` (1372): abandon (1373–1378) or escape to `canonicalize(2·next, 2a, 2b)` with `previous` reset (1379–1390); else `length > walk_cap` → abandon (1392–1397); else shift `previous` and move (1398–1404) |
| `return None` 1407 / `return Some` 1299 | `status ∈ {"capped", "solved"}`, `Done` | explicit terminal states that stutter, so TLC's deadlock check stays on and means "a running state with no move" |

**Lane order.** With `LaneOrderND = FALSE` each phase takes its pending
lanes in index order, as the Rust's `for` loops do; the whole system is then
deterministic after `Init` (outdegree 1). With `LaneOrderND = TRUE` any
pending lane may go next in every phase, so TLC also explores every
interleaving of phase 1 (which lane receives which start), phase 2 (which of
two lanes on the same distinguished point inserts and which collides) and
phase 3.

**Batch inversion.** Phase 3 in the Rust first computes all partitions and
slope denominators, inverts them together (`batch_inv`, 1349–1350) and then
applies the additions, falling back to the full group law when a
denominator is zero (1354–1362). A lane's step reads only that lane's state,
so stepping lanes one at a time is equivalent; field arithmetic itself is
abstracted to exact addition in Z_R.

### Ghost variables (not in the Rust)

* `sentinel`: counts fruitless detections where `next` is infinity, is not
  the current point and is not one of the `previous` slots actually
  written since the start or last escape (`nprev`) — that is, detections
  caused only by the `[∞; 4]` initial fill.
* `longest`: the longest genuine cycle a fruitless detection has caught.

They only feed the reachability witnesses below.

## Properties

* `TypeOK`.
* `RepInv`: every lane (live or not) satisfies `pt = a + b·d (mod R)`, every
  table entry `k ↦ ⟨a, b⟩` satisfies `k = a + b·d`, and
  `cursor = cursor_a + d`.
* `CanonInv`: every live lane stands on the canonical representative of its
  orbit; every table key is canonical and distinguished.
* `NoFailedCollision`: `failed = 0`.
* `SolvedCorrect`: `status = "solved" ⇒ result = d`.
* Deadlock freedom (TLC's default check): every running state has a
  successor; the only terminal states are `solved` and `capped`.
* `Terminates` (liveness, under `WF_vars(Next)`): every behaviour reaches
  `solved` or `capped`.
* `WalkCapRespected` (diagnostic): a live lane never has `len > WalkCap`.

## Parameter sets

Every config explores all of `stride_a ∈ [StrideLo, StrideHi]`,
`cursor_a ∈ [CursorLo, CursorHi]` and `d ∈ [DLo, DHi]` (full ranges
`[1, R)`, `[0, R)`, `[0, R)` for R = 31 and 61; for R = 113, stride and
cursor in `[1, 16]` and `[0, 16]`, d over all of `[0, 113)`).
`StepCap = 14` everywhere. The arrow lists the walk graph on canonical
points (`*` = distinguished), which is what makes each set interesting:

| name in cfg | R, LAM (order, A) | jumps | partition, DP | WalkCap | canonical walk graph |
|---|---|---|---|---|---|
| `r31_l5_T2_*` | 31, 2 (5, 10) | {2, 9} | `(3x+1)%2`, `(5x+2)%3=0` | 4 | 0→5\*, 1→3, 3→5\*, 5\*: every walk reaches a DP |
| `r31_l5_T3_*` | 31, 2 (5, 10) | {5, 28, 20} | `(2x+1)%3`, `(2x+1)%3=0` | 4 | 0→3, 1\*, 3→0, 5→3: a 2-cycle through infinity |
| `r31_l3_*` (except T4) | 31, 5 (3, 6) | {28, 21} | `(3x+1)%2`, `(2x+1)%3=0` | 4 | 0→2, 1\*, 2→8, 3→0, 4\*, 8→2 |
| `r31_l3_T4_*` | 31, 5 (3, 6) | {4, 17} | `x%2`, `x%5=0` | 6 | contains a DP-free 5-cycle |
| `r61_l3_*` | 61, 13 (3, 6) | {2, 59, 9} | `(3x+1)%3`, `(x+1)%2=0` | 3 | 12→8→6→4→2→0→2: walks longer than walk_cap |
| `r113_l7_*` | 113, 16 (7, 14) | {108, 21} | `(3x+1)%2`, `(2x+1)%3=0` | 4 | 0→3, 2→3, 3→2, 5→0, 6→5, 9→1\*, 10\*, 13\*, 1\* |

R = 113 with λ of order 7 (A = 14) is the closest analogue of a Koblitz
curve: λ has the order of the field degree, and 2 is not in `{±λ^k}`. In
the two R = 31, LAM = 2 sets, 2 is an automorphism, so the escape's doubling
lands back on the same canonical point (only `a, b` change); the escape
configs on those sets check `RepInv` for the doubled coefficients but do
not exercise an escape that moves the walk. The R = 31 LAM = 5, R = 61 and
R = 113 sets do.

## Results

Final run (`out/SUMMARY.txt`, TLC 2026.10.15 from the v1.8.0 release URL,
2 workers). "Distinct" and "depth" are TLC's distinct states and the depth
of the state graph; for a config that stops at a violation they are the
numbers at the stop, which vary by a few states between runs because the
two workers race to the first violation (the depth, i.e. the length of the
shortest counterexample, does not). Pass verdicts are exhaustive and their
counts are exact. All 41 configs met their expectation (`failures: 0`).

| config | setting | expected | TLC verdict | distinct states | depth | time |
|---|---|---|---|---|---|---|
| `diag_walkcap_abandon` | Diagnostic: with escape off, a live walk never exceeds walk_cap (r=31, lambda=5) | pass | pass | 1006818 | 76 | 13s |
| `diag_walkcap_escape` | Diagnostic: with escape on, a live walk can exceed walk_cap because escape steps skip the length check (r=113, lambda=16, 1 lane) | violate WalkCapRespected | violate WalkCapRespected | 855314 | 31 | 10s |
| `live_r31_l3_L2_escape` | Liveness: under weak fairness every behaviour ends solved or step-capped (r=31, lambda=5, 2 lanes, escape) | pass | pass | 1036206 | 71 | 47s |
| `live_r31_l5_T2_L2_abandon_nd` | Liveness, lanes in any order within a phase (r=31, lambda=2, jumps {2,9}, 2 lanes, abandon) | pass | pass | 1194430 | 39 | 54s |
| `mut_canon_b` | MUTANT canon_b: canonicalize multiplies a but not b by the orbit multiplier (r=31, lambda=5) | violate RepInv | violate RepInv | 58596 | 3 | 2s |
| `mut_canon_b_downstream` | MUTANT canon_b with RepInv not checked, escape on (doubling makes b vary): the broken representation reaches a failed collision | violate NoFailedCollision | violate NoFailedCollision | 606380 | 22 | 10s |
| `mut_canon_b_never_solves` | MUTANT canon_b, escape off: b stays 1 in every state, so every table hit has b = b' (wasted) and the solver can never succeed; W_Solved holds | pass | pass | 2638410 | 116 | 28s |
| `mut_escape_b` | MUTANT escape_b: the fruitless-cycle escape doubles the point and a but not b (r=31, lambda=5) | violate RepInv | violate RepInv | 229650 | 9 | 5s |
| `mut_sign` | MUTANT sign: candidate = (a - hit_a)(b - hit_b)^-1; the [d]G = Q check rejects it | violate NoFailedCollision | violate NoFailedCollision | 173933 | 7 | 5s |
| `mut_sign_noverify` | MUTANT sign_noverify: wrong sign and the candidate returned without [d]G = Q | violate SolvedCorrect | violate SolvedCorrect | 173936 | 7 | 4s |
| `mut_sign_solved` | MUTANT sign with NoFailedCollision not checked: the [d]G = Q check keeps SolvedCorrect | pass | pass | 2585778 | 116 | 36s |
| `r113_l7_L1_abandon` | r=113, lambda=16 (order 7, A=14), jumps {108,21}, 1 lane, abandon; stride in 1..16, cursor in 0..16, d in 0..112 | pass | pass | 1586627 | 113 | 18s |
| `r113_l7_L2_escape` | r=113, lambda=16 (order 7, A=14), jumps {108,21}, 2 lanes, escape; stride in 1..16, cursor in 0..16, d in 0..112 | pass | pass | 1280525 | 78 | 21s |
| `r31_l3_L1_abandon` | r=31, lambda=5 (order 3, A=6), jumps {28,21}, 1 lane, abandon | pass | pass | 1464750 | 107 | 18s |
| `r31_l3_L2_abandon` | r=31, lambda=5 (order 3, A=6), jumps {28,21}, 2 lanes, abandon | pass | pass | 1006818 | 76 | 15s |
| `r31_l3_L2_abandon_nd` | r=31, lambda=5 (order 3, A=6), jumps {28,21}, 2 lanes, abandon, lanes in any order within a phase | pass | pass | 2422650 | 76 | 51s |
| `r31_l3_L2_escape` | r=31, lambda=5 (order 3, A=6), jumps {28,21}, 2 lanes, escape | pass | pass | 1036206 | 71 | 18s |
| `r31_l3_L3_escape` | r=31, lambda=5 (order 3, A=6), jumps {28,21}, 3 lanes, escape | pass | pass | 867318 | 68 | 18s |
| `r31_l3_T4_L1_abandon` | r=31, lambda=5 (order 3, A=6), jumps {4,17}, DP 1 in 5, walk_cap 6, 1 lane, abandon | pass | pass | 2568660 | 90 | 26s |
| `r31_l3_T4_L2_escape` | r=31, lambda=5 (order 3, A=6), jumps {4,17}, DP 1 in 5, walk_cap 6, 2 lanes, escape | pass | pass | 1731660 | 61 | 26s |
| `r31_l5_T2_L1_abandon` | r=31, lambda=2 (order 5, A=10), jumps {2,9}, 1 lane, abandon | pass | pass | 696570 | 53 | 9s |
| `r31_l5_T2_L2_abandon` | r=31, lambda=2 (order 5, A=10), jumps {2,9}, 2 lanes, abandon | pass | pass | 519250 | 39 | 8s |
| `r31_l5_T2_L2_escape` | r=31, lambda=2 (order 5, A=10), jumps {2,9}, 2 lanes, escape | pass | pass | 519250 | 39 | 9s |
| `r31_l5_T3_L1_escape` | r=31, lambda=2 (order 5, A=10), jumps {5,28,20}, 1 lane, escape | pass | pass | 2166280 | 96 | 27s |
| `r31_l5_T3_L2_abandon` | r=31, lambda=2 (order 5, A=10), jumps {5,28,20}, 2 lanes, abandon | pass | pass | 980530 | 78 | 16s |
| `r31_l5_T3_L2_escape` | r=31, lambda=2 (order 5, A=10), jumps {5,28,20}, 2 lanes, escape | pass | pass | 1457310 | 71 | 25s |
| `r31_l5_T3_L2_escape_nd` | r=31, lambda=2 (order 5, A=10), jumps {5,28,20}, 2 lanes, escape, lanes in any order within a phase | pass | pass | 3516950 | 71 | 109s |
| `r61_l3_L2_abandon` | r=61, lambda=13 (order 3, A=6), jumps {2,59,9}, 2 lanes, abandon | pass | pass | 9959592 | 85 | 126s |
| `r61_l3_L2_escape` | r=61, lambda=13 (order 3, A=6), jumps {2,59,9}, 2 lanes, escape | pass | pass | 11323674 | 79 | 159s |
| `wit_cycle5` | Reachability witness: TLC must violate W_Cycle5, i.e. a genuine cycle of length 5 (not 4) is caught as fruitless (r=31, lambda=5, jumps {4,17}) | violate W_Cycle5 | violate W_Cycle5 | 864899 | 31 | 9s |
| `wit_escape` | Reachability witness: TLC must violate W_Escape, i.e. a fruitless cycle is escaped by doubling | violate W_Escape | violate W_Escape | 228791 | 9 | 5s |
| `wit_fruitless` | Reachability witness: TLC must violate W_Fruitless, i.e. a walk is abandoned in a fruitless cycle | violate W_Fruitless | violate W_Fruitless | 228794 | 9 | 5s |
| `wit_infinity_start` | Reachability witness: TLC must violate W_InfinityStart, i.e. a walk start lands on infinity and the start loop continues | violate W_InfinityStart | violate W_InfinityStart | 28862 | 2 | 1s |
| `wit_lane_at_infinity` | Reachability witness: TLC must violate W_LaneAtInfinity, i.e. a live lane stands on infinity | violate W_LaneAtInfinity | violate W_LaneAtInfinity | 227851 | 9 | 5s |
| `wit_sentinel` | Reachability witness: TLC must violate W_Sentinel, i.e. a step to infinity is taken for a fruitless cycle because previous starts filled with infinity | violate W_Sentinel | violate W_Sentinel | 228791 | 9 | 4s |
| `wit_solved` | Reachability witness: TLC must violate W_Solved, i.e. a collision with b != b' is resolved and verified | violate W_Solved | violate W_Solved | 172988 | 7 | 4s |
| `wit_step_capped` | Reachability witness: TLC must violate W_StepCapped, i.e. the global step cap ends the solve with None | violate W_StepCapped | violate W_StepCapped | 1004970 | 75 | 13s |
| `wit_table_two` | Reachability witness: TLC must violate W_TableTwo, i.e. the table holds two distinguished points | violate W_TableTwo | violate W_TableTwo | 172987 | 7 | 5s |
| `wit_walk_capped` | Reachability witness: TLC must violate W_WalkCapped, i.e. a walk is abandoned at walk_cap (r=31, lambda=5, 1 lane, escape) | violate W_WalkCapped | violate W_WalkCapped | 797945 | 31 | 10s |
| `wit_walk_capped_abandon` | Reachability witness: TLC must violate W_WalkCapped, i.e. a walk is abandoned at walk_cap with escape off (r=61, lambda=13) | violate W_WalkCapped | violate W_WalkCapped | 6230452 | 33 | 53s |
| `wit_wasted` | Reachability witness: TLC must violate W_Wasted, i.e. a table hit with b = b' is counted as wasted | violate W_Wasted | violate W_Wasted | 354981 | 14 | 6s |

### Expected-pass models

The 18 model configs (`r*`) pass with no invariant violated and no
deadlock: lanes 1, 2 and 3; escape on and off; lanes in index order and
in every order; R = 31 with two different λ, R = 61 and R = 113; six
jump tables. Both liveness runs (`live_*`, `FairSpec`, property
`Terminates`) pass: every behaviour ends `solved` or `capped`.

### Reachability witnesses (`wit_*`)

Each `W_*` is an invariant TLC must violate; the counterexample is a
behaviour that takes the branch named, so the safety results above are not
vacuous. Every one is violated as expected: a verified solve, a wasted
`b = b'` hit, the start loop skipping infinity, the global step cap, a
fruitless abandon, an escape, the walk cap (escape on and off), a live lane
standing on infinity, two table entries, a detection caused only by the
infinity sentinel, and a genuine 5-cycle.

### Mutants (`mut_*`)

| mutant | first property TLC reports violated | shortest trace |
|---|---|---|
| `canon_b` (`mut_canon_b`): `canonicalize` does not multiply `b` | `RepInv` | depth 3: the first walk start, `b` stays 1 while the point is scaled |
| `canon_b`, `RepInv` not checked, escape on (`mut_canon_b_downstream`) | `NoFailedCollision` | depth 22: a collision of two states that do not represent their point |
| `canon_b`, escape off (`mut_canon_b_never_solves`) | nothing: `SolvedCorrect`, `NoFailedCollision` and `W_Solved` all hold | `b = 1` in every state, so every table hit has `b = b'` and is wasted; the solver can never succeed and every behaviour ends step-capped |
| `sign` (`mut_sign`): candidate `(a − hit_a)(b − hit_b)^{-1}` | `NoFailedCollision` | depth 7: `d = 1`, the candidate is `−d = 30`, which `[d]G = Q` rejects |
| `sign`, `NoFailedCollision` not checked (`mut_sign_solved`) | nothing: `SolvedCorrect` holds | the `[d]G = Q` check stops every wrong candidate |
| `sign_noverify` (`mut_sign_noverify`): wrong sign and no `[d]G = Q` check | `SolvedCorrect` | depth 7: returns 30 for `d = 1` |
| `escape_b` (`mut_escape_b`): the escape doubles the point and `a`, not `b` | `RepInv` | depth 9 |

So `SolvedCorrect` rests on the `[d]G = Q` check alone: a scalar is only
returned after it passes, so no mutant but the one that drops the check
can return a wrong `d` (`mut_sign_solved` and `mut_canon_b_never_solves`
confirm it exhaustively for those two bugs). A wrong representation
shows up only as failed collisions or as a solver that never succeeds,
which is why `RepInv` and `NoFailedCollision` are the properties worth
watching (the `debug_assert!`s added in `448be6146` at 1282–1285 and
1315–1330 check exactly these).

## Findings about the Rust

No property failed on the unmutated model: within these instances the
walk keeps `P = a·G + b·Q` in every lane and table entry, every live lane
and table key is canonical, no collision fails, a returned scalar is always
`d`, nothing deadlocks and every behaviour terminates. No bug that yields a
wrong answer was found. What the model did surface:

1. **The infinity sentinel in `previous`** (`[RawPointG::Infinity; 4]` at
   1227, 1262 and 1389; the test at 1372). A walk whose step lands on the
   point at infinity within four steps of its start or of an escape is
   reported as a fruitless cycle although it has not repeated a point
   (witness `wit_sentinel`). With `escape_fruitless` on, the escape doubles
   infinity to infinity, so the lane stays on infinity with `a, b` doubled
   (`diag_walkcap_escape` trace, states 13–19). Correctness is unaffected
   (`RepInv` holds at infinity) and at real sizes a step hits infinity with
   probability about `1/r`, so this only matters for toy groups.
2. **Escape steps skip the walk cap.** The escape branch ends with
   `continue` (1390) before the `length > walk_cap` test (1392), so a live
   walk can be longer than `walk_cap` (`diag_walkcap_escape`: `WalkCapRespected`
   fails; it holds with escape off, `diag_walkcap_abandon`). The walk is
   still capped at its next ordinary step, so the overshoot is bounded by
   the run of consecutive escapes, but `abandoned_trail_lengths` can hold
   lengths above `walk_cap + 1` when escape is on.
3. **Cycles of length 5 are caught, not "≤ 4".** `previous` holds the four
   states before the current one, so `next ∈ previous ∨ next == current`
   catches cycles of length 1 to 5 (`wit_cycle5`). The field doc of
   `StrongRhoOutcomeG::fruitless` (line 499) says "length ≤ 4". Documentation only.
4. **"Two merged walks escape alike" holds for the escape map but not for
   the detection** (reasoned, not model-checked). The escaped state is a
   function of `next`, as the `escape_fruitless` doc says, but whether a
   step is flagged depends on each walk's own `previous`, which still holds
   pre-merge points (and sentinel slots) for four steps after a merge. Two
   walks that merge just before entering a short cycle can flag it at
   different points and escape to different states. This costs efficiency
   (a lost merge), never correctness.

## What is abstracted, and what this does not show

* The curve group is replaced by Z_R. Field arithmetic, the group law,
  `batch_inv`, `add_with_slope` and its fallback, and `scalar_mul` are
  exact operations in Z_R. Their correctness is not checked here.
* `canonicalize` is modelled by its specification (least orbit member and
  its multiplier), not by the normal-basis least-rotation algorithm or the
  y-based sign choice; there are no ties in Z_R except at 0.
* The 64-bit hash mixers (`partition`, `dp_hash`) are replaced by fixed
  affine maps; 32 jumps by 2 or 3; `walk_cap = 8·2^dp_bits` and
  `step_cap ≥ 10^6` by small constants.
* The seeded `StdRng` draws (`stride_a`, `cursor_a`) and the target `d` are
  nondeterministic choices in `Init` over the ranges above; the jump table
  is fixed per config. The Rust solves one target per call and restarts
  nothing else, so there are no other random restarts to model.
* `d = 0` (Q = infinity) is included; the Rust does not exclude it.
* Counters that only feed statistics (`charges`, trail-length vectors) are
  not modelled, except `walks`, `fruitless`, `escapes`, `capped`, `wasted`
  and `failed`.
* Scalars are reduced mod R everywhere; the `u64`/`u128` widths, `MulMod`
  and overflow are not modelled.

TLC enumerates every reachable state of each instance, so each verdict is
exhaustive for that instance: those R, λ, jump tables, hashes, caps and
ranges. It is not a proof for every r, every jump table or the real hash
functions; small-instance model checking finds protocol mistakes that
appear in small groups (as all four mutants do within 3–9 steps), and the
Lean work in `../lean/` is where general statements belong.
