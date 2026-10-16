# TLA+ model of the strong distinguished-point rho

`StrongRho.tla` models the protocol of `StrongRhoG::solve` in
`src/cryptanalysis/koblitz_strong_rho.rs` (file SHA-256
`54978d00bbb2c50f090df4f6a649deb830ace9b7b356d2b0feaefa30017f8a16`, last
changed in commit `680f3cf9`): lockstep lanes, walk starts from a strided
cursor, signed-Frobenius canonicalization, distinguished points, the
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
| `from_parts` 841–872: `lam_pow`, `assert!(cur == S::one())` | `LamPow`, `ASSUME Md(LamPow[NF-1]*LAM) = 1`, `Mults`, `A` | the multiplier table and the order check |
| `canonicalize` 1072–1103 (with `least_rotation` 658, `choice_exhaustive` 1040) | `CanonMul`, `Canon`, `CanonState` | canonical point = least element of the orbit `{m·x : m ∈ Mults}`; `a` and `b` are multiplied by the multiplier `m`; infinity is returned unchanged (multiplier 1). The normal-basis rotation that computes the least orbit member and its `±λ^k` is **not** modelled: the model takes the least residue directly. Whether rotation by `k` really corresponds to `λ^k` is the subject of the Rust unit test `canonical_form_is_orbit_invariant_and_multiplier_is_correct` (line 1498), not of this model. |
| `partition` 765, `dp_hash` 778, `dp_mask` 1143 | `Partition(x) = (x*C1 + C2) % J`, `IsDP(x) = (x*C3 + C4) % D = 0` | the 64-bit mixers are replaced by fixed affine functions of the canonical point; `D` plays `2^dp_bits` |
| `jumps` 1106–1118 (`JUMPS = 32`, `s` uniform in `[1, r)`) | `NJ` ∈ {2, 3} constants `J1..J3` | a fixed jump table, two or three nonzero scalars; four different tables across the configs |
| `solve` 1136–1141: `stride_a`, `stride`, `cursor_a`, `cursor = cursor_a·G + Q` | `Init`: `d`, `strideA`, `cursorA` chosen from ranges, `cursor = cursorA + d` | the seeded RNG draws become nondeterministic choices, so TLC explores every target and every start sequence in the ranges |
| `walk_cap = 8 << dp_bits` 1144, `step_cap` 1145–1148 | constants `WalkCap`, `StepCap` | small numbers instead of `8·2^dp_bits` and `≥ 10^6` |
| `table: HashMap` 1150, `lanes` 1151–1163 | `table` (function point ↦ ⟨a, b⟩), `lane` (record per lane: `pt, a, b, prev, len, live`, plus the ghost `nprev`) | `previous` is a 4-tuple initialised with infinity, exactly as in the Rust |
| `while steps < step_cap` 1174 | `Begin` (phase `begin`) | when `steps ≥ StepCap` the state becomes `capped` (`None`) |
| phase 1, start loop 1175–1197 | `StartLane(l)` (one turn of the `loop`), `EndStart` | takes `cursor`, advances `cursor += stride`, `cursor_a += stride_a`, `walks += 1`; if the taken point is infinity the same lane loops again (`continue`, 1184); otherwise the lane gets `canonicalize(start, start_a, b = 1)`, `previous = [∞; 4]`, `length = 0`, `live` |
| phase 2, distinguished points 1199–1263 | `DpLane(l)`, `EndDp` | lane leaves (`live = false`); table miss → insert (1214–1219); hit with `b = b'` (no inverse) → `wasted` (1221–1225); otherwise `candidate = (hit_a − a)·(b − hit_b)^{-1}` (1226–1227), returned as `solved` if `candidate·G = Q` (1229), else `failed_collisions += 1` (1262) |
| phase 3, one step per live lane 1264–1337 | `StepLane(l)`, `EndStep` | `next = canonicalize(P + s_j·G, a + s_j, b)`; `steps += 1`, `length += 1` (1301–1302); fruitless if `next == current ∨ next ∈ previous` (1303): abandon (1304–1309) or escape to `canonicalize(2·next, 2a, 2b)` with `previous` reset (1310–1321); else `length > walk_cap` → abandon (1323–1327); else shift `previous` and move (1329–1335) |
| `return None` 1338 / `return Some` 1230 | `status ∈ {"capped", "solved"}`, `Done` | explicit terminal states that stutter, so TLC's deadlock check stays on and means "a running state with no move" |

**Lane order.** With `LaneOrderND = FALSE` each phase takes its pending
lanes in index order, as the Rust's `for` loops do; the whole system is then
deterministic after `Init` (outdegree 1). With `LaneOrderND = TRUE` any
pending lane may go next in every phase, so TLC also explores every
interleaving of phase 1 (which lane receives which start), phase 2 (which of
two lanes on the same distinguished point inserts and which collides) and
phase 3.

**Batch inversion.** Phase 3 in the Rust first computes all partitions and
slope denominators, inverts them together (`batch_inv`, 1280) and then
applies the additions, falling back to the full group law when a
denominator is zero (1285–1294). A lane's step reads only that lane's state,
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

RESULTS_PLACEHOLDER
