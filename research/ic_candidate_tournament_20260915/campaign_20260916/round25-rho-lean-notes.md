# Round 25 candidate: a lean matched rho

`round25-rho-lean.patch` (`-p1`, against `round24-sources/matched`) makes the
round-0024 matched rho cheaper per step without changing its walk. On all 72
round-0024 fixtures measured, the lean and matched rho took **the same walk**.
Iterations, walk additions and restarts were identical on every fixture, and
both recovered the same, oracle-verified logarithm. The patch lowered the
`rho_solve` phase to **0.28–0.58×** of the matched rho's instructions, by cell.
At `n43a1` the cost per walk addition fell from **2,620 to 1,497**
instructions, pooled over six fixtures.

Everything here is callgrind instruction counts (`Ir`) in one container. Native
time was not measured.

## 1. Where the matched rho's instructions went

Profiles are of `rho_solve`, the worker's own client-request dump around the rho
call. Each fixture is `development/<cell>-000`, run on a line-table build of the
same source. The categories are inclusive costs by caller; the tooling is
described in §6.

| `n43a1`, 7,096 walk additions | matched Ir/addition | lean Ir/addition |
|---|--:|--:|
| setup: `[a]G + [b]Q` for 16 jumps + 32 walks | **869** | 143 |
| walk addition arithmetic | **855** (of which the Fermat inverse: 264) | ≈ 600 (batched Euclid: 291; inlined λ, x₃, y₃: 166; bounds-checked reduction/loop code, mostly: 145) |
| canonicalisation | 704 (of which `NormalBasis::canon`, IC code: 259) | 616 (same `canon`: 259) |
| collision candidate check | 53 | 6 |
| rest: examine loop, caches, stored table, hashes, basis setup | 203 | 139 |
| **total `rho_solve` / walk additions** | **2,684** | **1,507** |

| `n23a1`, 398 walk additions | matched | lean |
|---|--:|--:|
| setup | **3,143** | 654 |
| walk addition arithmetic | 1,099 (of which the Fermat inverse: 647; only 5 walks share it) | ≈ 596 |
| canonicalisation | 553 | 491 |
| collision candidate check | **429** | 31 |
| rest | 344 | 231 |
| **total** | **5,568** | **2,003** |

The "~3,000 Ir per walk addition" figure divides `rho_solve` by walk additions,
so it includes setup and the collision check amortised over the walk. Most of
the gap to an IC probe was not in the step at all:

1. **Setup was a third of `rho_solve` at `n43a1` and over half at `n23a1`.**
   Every jump and every walk start was built with two López–Dahab
   double-and-add products and three Fermat inversions. There are two
   inversions in `ld_to_affine` and one in `FastCurve::add`, each `n`
   squarings and `n` multiplications.
2. **The batched addition inverted by Fermat** through the library `Gf2`:
   264 Ir/addition at `n43a1` with 32 walks sharing an inversion, and 647 at
   `n23a1` with 5. `Gf2` also squares by bit-spreading, and its reduction
   indexes a `Vec` with bounds checks.
3. **The candidate check ran in the general `BinaryPoint` arithmetic**, which
   allocates on every field operation: 53 Ir/addition at `n43a1` and 429 at
   `n23a1`, for one check.
4. **Smaller items in the step:**
   - the partition hash was computed twice per step, and the mix twice when
     examining a state;
   - two 128-bit `%` reductions (`__umodti3`) scaled the coefficients;
   - two `%` reductions summed them;
   - the basis changes were done as `coords`, rotate and `from_coords`,
     including when `k = 0`.

What is not overhead, because the IC probe pays it too or it is intrinsic to rho:

- `NormalBasis::canon` (259 Ir/addition at `n43a1`) is the IC arm's own code,
  and it is unchanged.
- One more field multiplication than a probe, for the ordinate y₃.
- The ordinate's basis change and the abscissa's inverse basis change (about
  170–220 Ir/addition at `n43a1`, with the IC arm's nibble tables).
- The coefficient bookkeeping.
- The recent-state cache and the stored-point table.

## 2. What the patch changes

Only rho code paths change. The walk is unchanged:

- the same iteration function: partition, jumps, `+J` and signed-Frobenius
  canonical form;
- the same class names (normal-basis least rotation, with the loop fallback
  for `x ∈ {0, 1}`) and the same sign rule;
- the same recent cache and stored-point rule;
- the same fruitless-cycle escape;
- the same restart and `max_iterations_per_restart` (`max_trials`) semantics;
- the same RNG draw order and the same charge ledger;
- the same worker report and final verification (the worker is untouched).

Each change computes the same value more cheaply.

1. **Field arithmetic is the IC arm's.** The walk uses `koblitz_tiny_ic::Arith`
   (made `pub(crate)`, no behaviour change): carry-less squaring, a Euclidean
   inverse, and its `batch_inv`, `add`, `add_with_lambda` and `double`. The
   matched walk used `FastCurve` and `Gf2` for these. This is the same parity
   argument by which round 0024 gave rho the IC arm's `NormalBasis`.
2. **Setup uses the IC arm's fixed-base product.** Each restart draws the
   same `(a, b)` scalars in the same order, jumps first and then walks. It
   forms every `[a]G` and every `[b]Q` with `Arith::fixed_base_mul_many` from
   doubling tables of G and Q, built once per run. It then adds the pairs with
   one batched inversion. The points are the same group elements, so the
   jumps and starts are identical.
3. **The candidate check** `[d]G = Q` runs in the same single-word arithmetic
   (`Arith::fixed_base_mul` on G's table), as the IC arm checks its own
   candidates. The boolean is the same. The worker still verifies the reported
   log in the general arithmetic, and the oracle checks it independently.
4. **The advance hashes the partition once per walk** and reuses the index for
   the coefficient update. The examine loop computes one mix for both the
   cache slot and the storage rule.
5. **Canonicalisation.**
   - When `k = 0` there is no basis change, since the point is already the
     representative up to sign.
   - Otherwise the ordinate's round trip is one call,
     `NormalBasis::frobenius_poly`. It uses the same nibble tables and returns
     the same value as `from_coords(frobenius_coords(coords(y), k))`.
   - `from_coords`, which only rho calls, consumes the word a nibble at a
     time.
   - `coords`, which the IC arm's `canon` uses, is untouched.
6. **Coefficient arithmetic.** Scaling by ±λ^k uses a Shoup precomputed
   quotient per (k, sign), which is exactly equal to `mulmod_u64`. The jump
   sums use a conditional subtraction instead of `%`.

The matched path is kept: `lean = false`, or any curve outside the lean path's
conditions, runs the old code. The lean path requires `k = 1`, `n ≤ 61`, a
normal basis, `r < 2^62` and G, Q ≠ O. Two new unit tests use this:

- **`lean_rho_walk_matches_the_matched_walk_exactly`** covers 9 curves shaped
  like the tournament's cells from `n = 13` to `n = 43`, with 1, 4 and 32 walks
  each. It requires identical iterations, restarts, jump-table rebuilds, walk
  counts, the full charge ledger and the progress-event sequence, and the
  right log.
- **`lean_rho_setup_points_match_double_and_add`** checks the batched setup
  points against `fc.add(fc.mul_u64(G, a), fc.mul_u64(Q, b))`, the Shoup
  products against `mulmod_u64`, and the fast candidate check.

`normal_basis_inverts_and_rotates_as_squaring` now also checks
`frobenius_poly` against repeated squaring.

**Why this is fair.** The walk's point sequence is not merely "the same
distribution". It is the same sequence: every fixture's iterations, additions
and restarts matched, and so did the unit test's full ledger. Walk-length
fairness therefore holds by identity rather than by sampling.

The field and curve arithmetic, the setup product and the candidate check are
the IC arm's own functions, so rho gets no field operation cheaper than the IC
arm's. Both use the same reduction tables and the same non-inlined carry-less
multiply call.

Two things were tightened only on rho's side:

- the rho-only operations, which have no IC counterpart: coefficient scaling,
  and the inverse basis change `from_coords`;
- the forward basis change inside `frobenius_poly`. This is the same
  nibble-table walk as the IC arm's `coords`, written in the
  nibble-consuming form. Whether that form costs fewer instructions than
  `coords`' shift-by-4i loop was not measured separately; at 11 nibbles
  (`n43a1`) any difference is a few tens of Ir per addition at most.

## 3. Tests

On the lean tree (the raw output is appended to `round25-rho-lean-check.txt`):

- **`cargo test --release --offline --locked --lib rho`:** 41 passed, 1
  failed, 4 ignored. The one failure is
  `single_word_rho_walk_matches_the_reference_step_for_step`, at the
  iterations assertion (6 vs 5 on its first curve). That test pins the fast
  walk to the reference walk's lexicographic representative, which the
  normal-basis naming of round 0024 already changed. Since the lean walk is
  the matched walk, the failure is the matched tree's own. **Nothing else
  fails.** The 41 passes include the two new tests.
- **`cargo test ... --lib cryptanalysis::koblitz_tiny_ic`:** 10 passed, 0
  failed.

The matched tree's tests were not re-run here.

## 4. Measurement (`round25_rho_lean_check.py`, output in `round25-rho-lean-check.txt`)

The fixtures are round 0024's development and selection cases, 3 + 3 per cell.
For `n19a1`, `n29a1`, `n59a0` and `n61a1` they are the first 6 confirmation
cases instead. The config is `batch_trials 1, sparse, max_trials 65536,
pair_table, 3 summands`, in `mode: rho`.

- Every report passed `oracle.verify(..., expected_mode='rho')`. The oracle
  accepts only d < r with [d]G = Q, which is the fixture's logarithm.
- Lean and matched agreed on every logarithm.

The bands are 95% intervals from the t distribution over each cell's six log
ratios. Ir per addition is `rho_solve` Ir over walk additions, pooled over the
cell.

| cell | fixtures | Ir lean/matched (whole run) | `rho_solve` lean/matched | median iterations (matched = lean) | iteration ratio | Ir/addition matched | Ir/addition lean |
|---|--:|--:|--:|--:|--:|--:|--:|
| `n13a0` | 6 | 0.485 [0.477, 0.492] | **0.284 [0.268, 0.302]** | 9 | 1.000 [1.000, 1.000] | 56,719 | 16,141 |
| `n17a1` | 6 | 0.456 [0.449, 0.463] | **0.314 [0.289, 0.341]** | 70.5 | 1.000 [1.000, 1.000] | 16,747 | 5,313 |
| `n19a0` | 6 | 0.467 [0.463, 0.471] | **0.319 [0.296, 0.343]** | 88.5 | 1.000 [1.000, 1.000] | 16,048 | 5,168 |
| `n19a1` | 6 | 0.448 [0.445, 0.451] | **0.311 [0.302, 0.319]** | 97.5 | 1.000 [1.000, 1.000] | 14,570 | 4,535 |
| `n23a0` | 6 | 0.483 [0.465, 0.501] | **0.373 [0.335, 0.415]** | 339 | 1.000 [1.000, 1.000] | 6,463 | 2,446 |
| `n23a1` | 6 | 0.451 [0.434, 0.468] | **0.339 [0.305, 0.375]** | 298.5 | 1.000 [1.000, 1.000] | 6,352 | 2,184 |
| `n29a1` | 6 | 0.488 [0.480, 0.496] | **0.276 [0.271, 0.281]** | 31 | 1.000 [1.000, 1.000] | 43,341 | 11,973 |
| `n31a0` | 6 | 0.440 [0.435, 0.444] | **0.310 [0.284, 0.340]** | 205 | 1.000 [1.000, 1.000] | 12,307 | 3,884 |
| `n37a0` | 6 | 0.465 [0.400, 0.539] | **0.402 [0.313, 0.516]** | 2108 | 1.000 [1.000, 1.000] | 4,098 | 1,737 |
| `n43a1` | 6 | 0.603 [0.570, 0.638] | **0.559 [0.515, 0.607]** | 6808.5 | 1.000 [1.000, 1.000] | 2,620 | 1,497 |
| `n59a0` | 6 | 0.951 [0.944, 0.959] | **0.579 [0.544, 0.616]** | 11878.5 | 1.000 [1.000, 1.000] | 2,869 | 1,683 |
| `n61a1` | 6 | 0.609 [0.584, 0.635] | **0.576 [0.540, 0.614]** | 12169.5 | 1.000 [1.000, 1.000] | 2,839 | 1,662 |

Over all 72 fixtures, the whole-run ratio is 0.515 [0.490, 0.542].

- **The walk-length ratio is exactly 1** on every fixture, so each band is
  degenerate. Median iterations are therefore equal.
- **At `n43a1`, per fixture:** matched 2,215–3,098 Ir/addition, lean
  1,422–1,593. The spread comes from setup amortised over walks of different
  lengths.
- **The whole-run ratio includes phases outside rho**: startup, curve and
  target construction, and the worker's final verification, which are the
  same in both arms. At `n59a0`, curve construction alone is about 259M Ir of
  a roughly 283M run, which is why that cell's whole-run ratio is 0.95 while
  its `rho_solve` ratio is 0.58.
- **Small cells stay dominated by fixed setup.** The per-addition figure there
  measures setup more than the step. At `n13a0` a run is 7–23 iterations with
  17 setup points, doubling tables, the normal basis and a 4,096-slot cache.

## 5. Caveats

- **The IC arm computes exactly the same thing, but its instruction count
  moved slightly.**
  - `--ic-check 1` ran one fixture per cell in `ic` mode through both builds.
    The reports were identical, excluding `elapsed_seconds`.
  - The IC hot phase, `collection_and_decomposition`, had identical Ir at
    every cell.
  - Other IC phases moved by a few hundred to about ten thousand Ir:
    `factor_base_and_tables` +99 to +10,254, while `log_certification` and
    `verify_filter_and_linear_algebra` went down by −116 to −1,552 each.
  - The net change is between −0.053% and +0.013% of the IC run. This is
    consistent with code-generation shifts from rho now also calling
    `Arith::doubling_table`, `fixed_base_mul_many` and the others. Callgrind
    counts of a single build also vary by a few Ir from run to run.
  - If the patch were adopted, the tree the IC arm runs from would carry these
    shifts. They were measured on one fixture per cell only.
- **This is still not the leanest possible rho.**
  - Per addition, the lean rho is about 1,500 Ir at `n43a1`, or about 1,360
    excluding setup. The IC arm's own primitives set that level: the table
    reduction with bounds-checked indexing and a non-inlined carry-less
    multiply, nibble-table basis changes, and `canon`. At `n43a1`, the batched
    inverse's three multiplications plus its share of one Euclidean inverse
    cost 291 Ir/addition.
  - Changing those primitives would change the IC arm too, or give rho
    arithmetic the IC arm lacks, so it was not done.
  - The "~1,100 Ir per IC probe" figure is the other lane's. It was not
    re-measured here.
- **Rare paths still run the matched code:**
  - the canonical-form fallback loop for `x ∈ {0, 1}` (still `Gf2` squarings);
  - the escape's `2·coefficient % r`;
  - the matched walk itself when the lean conditions fail. No tournament cell
    hits this (`n ≤ 61`, `k = 1`).
- **Timing fields:** the lean path builds its doubling tables and Shoup
  constants once, before the restart loop. They are therefore outside the
  report's `setup_ns` wall-clock field. The worker does not output that field,
  and the Ir above count the tables because they sit inside `rho_solve`.
- **Sample size:** six fixtures per cell, one configuration, and one fixture
  set (round 0024's, with no draw from a new seed). Instructions only.

## 6. Reproducing

```sh
S=<scratch>; C=research/ic_candidate_tournament_20260915/campaign_20260916
cp -a $C/round24-sources/matched $S/matched && cp -a $C/round24-sources/matched $S/lean
(cd $S/lean && patch -p1 -i $OLDPWD/$C/round25-rho-lean.patch)
for t in matched lean; do (cd $S/$t && cargo build --release --offline --locked --jobs 2 \
    --example ic_tournament_worker --target-dir $S/target) && \
  cp $S/target/x86_64-unknown-linux-musl/release/examples/ic_tournament_worker $S/$t.bin; done
python3 $C/round25_rho_lean_check.py $S/matched.bin $S/lean.bin --per-cell 6 --ic-check 1
(cd $S/lean && cargo test --release --offline --locked --lib rho; \
              cargo test --release --offline --locked --lib cryptanalysis::koblitz_tiny_ic)
```

The §1 profiles came from builds with `CARGO_PROFILE_RELEASE_DEBUG=line-tables-only`.
Those builds' `rho_solve` Ir matched the plain builds' to within about 200 Ir
at `n43a1`. The builds were run as
`valgrind --tool=callgrind --dump-instr=yes` on `development/n43a1-000` and
`development/n23a1-000`, and the `rho_solve` dump (`*.3`) was read with
`callgrind_annotate --inclusive=yes --tree=caller`.

The binaries measured have these sha256 values (also in the check output):

- matched `b2c291e1…732c3`
- lean `1ca7f379…90eda`
