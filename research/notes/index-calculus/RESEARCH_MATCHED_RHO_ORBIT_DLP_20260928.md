# Matched-arithmetic recheck: compact-orbit IC vs batched Kuhn–Struik rho (n=53, L=1,024)

Written and committed **before the matched-arithmetic binary is built or run**,
per `AGENTS.md` "state the boundary before measuring" / "declare the
falsification target in advance". Nothing below this line is a result.

## Why this note exists

PR #830 (merged 2026-09-27, "Compact-orbit DLP beats frozen Kuhn-Struik batched
rho at n=53") reports wall-clock ratios of 0.246 (n=53) and 0.339 (n=41) for
`examples/koblitz_orbit_dlp_fast.rs` (compact-orbit index calculus, "IC")
against `examples/koblitz_rho_batch_ks.rs` (batched Pollard rho, "KS") at
`L = 1,024` targets — see `docs/ic/BOUNDARY_TARGETS.md` "Compact-orbit
shared-log DLP" and
`research/sat_factor_base_review_20260908/autolab_orbit_extract_20260924/RESULT.md`
§3, "Growing n against batched rho at fixed L = 1,024", n=53 row: K=440, IC
14,926 ms vs rho 59,513 ms, ratio 0.246. The claim is correctly marked
`PENDING_INDEPENDENT_VALIDATION` and was never promoted to the scoreboard's
headline `vs_rho` row, so there is no live overclaim — but the comparison uses
wall clock (AGENTS.md §6, "what does not count": "Wall-clock time as the
headline"), and, more importantly, the two binaries do not spend the same
arithmetic per group-law step.

**The asymmetry, read from the current source:**

- `koblitz_rho_batch_ks.rs`'s `raw_canonicalize`, in `SignedFrobenius` mode
  (the mode this exact comparison uses), walks the signed-Frobenius orbit by
  calling `raw_square` — one full carryless-multiply-and-reduce over
  `GF(2^53)` — up to `2·curve.n` times per canonicalization, and
  canonicalization runs once per `raw_add` in the main walk loop.
- The same file's `raw_inverse` is Fermat's-little-theorem square-and-multiply
  (`value^(2^n − 2)`): `n − 1` squarings plus up to `n − 1` multiplications,
  called once inside every `raw_add`/`raw_double`, i.e. once per rho step.
- `koblitz_orbit_dlp_fast.rs` (the IC side of the same comparison) instead uses
  a `NormalBasis` (`Gf2::sqr`-power-of-2 conjugates) with Frobenius as an O(1)
  bit rotation (`NormalBasis::rotate`) and inversion via Itoh–Tsujii
  (`invert()`) using the same rotation for the `a^(2^k)` steps — about 7–8 real
  field multiplications total, against Fermat's up to ~104 (52 sqr + 52 mul at
  n=53).

This is the same failure mode PR #818 found and fixed in a different code path
(`research/ic_triple_counted_20260923/rho-normal-basis.patch`, applied to
`src/cryptanalysis/koblitz_tiny_ic.rs`'s internal rho) and that PR #869 is
independently redoing for a third code path right now. Neither touches
`examples/koblitz_rho_batch_ks.rs`. This note freezes a matched-arithmetic
rerun of exactly the n=53, L=1,024 cell from RESULT.md §3.

## The exact frozen protocol reproduced

From
`research/sat_factor_base_review_20260908/autolab_orbit_extract_20260924/growing_n_vs_batched_rho/growing_n_vs_batched_rho.sh`
and `chosen_K_n53.txt` (K=440):

```sh
KS=examples/koblitz_rho_batch_ks.rs   # built in release; see reproducibility gap below
IC=target/release/examples/koblitz_orbit_dlp_fast
KIC_RHO_DP_BITS=4 KIC_RHO_BATCH_CORPUS=n53-ks-growing-1024-v1 \
  "$KS" 53 0 signed_frobenius 1024 531310 > ks_n53.jsonl
# scalars_n53-ks-growing-1024-v1.txt = published_fixture_scalar from ks_n53.jsonl
"$IC" construct:53:0:440 scalars_n53-ks-growing-1024-v1.txt 7 ic_n53.jsonl
```

n = 53, a = 0, quotient mode `signed_frobenius`, L (fixtures) = 1,024, K
(orbit columns) = 440, batch_seed = 531310, dp_bits = 4, rank_seed = 7.
Target/scalar generation is a deterministic blake3 function of
`batch_seed` + corpus name (`n53-ks-growing-1024-v1`) + index (see `main()` in
`koblitz_rho_batch_ks.rs`), so a run with the identical corpus name and
batch_seed regenerates bit-identical targets; the run stays a true paired
comparison run to run.

**Reproducibility gap, disclosed:** the frozen binary path the original
script referenced,
`research/sat_factor_base_review_20260908/autolab_batched_rho_n53_20260922/frozen/v2/koblitz_rho_batch_ks`,
is not present in the repository (an uncommitted local build artifact from the
original lab run). This rerun instead builds the current, unmodified
`examples/koblitz_rho_batch_ks.rs` from this branch's `main`-derived source —
the same CLI/env-var interface and JSON schema, confirmed by reading the file
— as the baseline. If that historical binary differed from the current source
in some undocumented way, this rerun cannot detect it; it can only confirm or
refute the claim against the source as it stands today.

## What this rerun changes, and what it does not

**New file:** `examples/koblitz_rho_batch_ks_matched_arith.rs`, a copy of
`examples/koblitz_rho_batch_ks.rs` with exactly two functional changes (CLI
args, JSON schema, walk logic, DP table, jump table and negation handling are
otherwise identical, so the pair remains a controlled comparison):

1. `raw_canonicalize`'s `SignedFrobenius` branch reads the orbit off a normal
   basis (`NormalBasis`, built once from `crypto_lib::cryptanalysis::semaev_decomp::Gf2`
   over the curve's own irreducible polynomial) instead of walking it with
   `raw_square`.
2. `raw_inverse` is Itoh–Tsujii inversion with the `a^(2^k)` steps done as the
   same normal-basis rotation, instead of Fermat's square-and-multiply.

Both are copied from `koblitz_orbit_dlp_fast.rs`'s own `NormalBasis`/`invert()`
machinery (duplicated locally in the new file rather than factored into a
shared module, to keep the change minimal and scoped, and to leave both
protected files — `koblitz_rho_batch_ks.rs` and `koblitz_orbit_dlp_fast.rs` —
untouched).

**Correctness gates** (run as `cargo test --example
koblitz_rho_batch_ks_matched_arith`, required to pass before this file's
output is trusted for anything below):

1. `Gf2::mul`/`Gf2::sqr` agree with this file's own `raw_mul_field`/`raw_square`
   on thousands of random field elements, at every admitted n (7…53). This is
   the check that `Gf2`'s field representation and reduction polynomial
   actually match `KoblitzCurve`'s raw arithmetic: both are built from the
   identical `IrreduciblePoly` (`curve.curve.irreducible`, same `degree` and
   `low_terms`) returned by `KoblitzCurve::new`, which is a deterministic
   function of `(a, n)` with no randomness, so the two binaries independently
   construct the same field with no cross-process synchronization needed.
2. The normal-basis rotation reproduces `k` chained `raw_square` calls for
   every `k` in `0..n`, at every admitted n — the exact value the rewritten
   canonicalization loop now reads off in O(1).
3. The rewritten `raw_canonicalize` (`SignedFrobenius` mode) picks the exact
   same canonical representative point and the exact same coefficient
   multiplier as the original repeated-squaring implementation, replayed
   side by side on a batch of random points at every admitted n. The
   algorithm's output must not change; only its cost may.
4. The Itoh–Tsujii `raw_inverse` is a genuine multiplicative inverse
   (`raw_mul_field(x, inv(x)) == 1`) for hundreds of random nonzero field
   elements, at every admitted n.

## Operation-count unit: what it charges, and a discovered limitation

**Primary/headline metric: valgrind `--tool=callgrind` retired-instruction
counts**, on the whole process, for both the matched-arithmetic KS binary and
the unmodified `koblitz_orbit_dlp_fast` binary (AGENTS.md §10 explicitly
endorses this as the deterministic measure beside wall time). This charges
every instruction either process retires — process/library startup, curve and
factor-base construction, jump-table/index build, the rank stage, linear
algebra, every target, and final verification — identically and completely on
both sides, with no risk of a missed call site.

**Why not exact `Gf2::mul`/`Gf2::sqr` call counting on both sides**, which was
attempted first and is the more direct measure of the specific claim (real
field multiplications/squarings, since normal-basis rotations and
`Linear`-table applies are O(1) word ops by comparison): the matched-arithmetic
KS file is fully self-contained (`raw_mul_field`/`raw_square` wrapped with
global counters at their own two call sites, at
`examples/koblitz_rho_batch_ks_matched_arith.rs`), so its count is exact and
complete. `koblitz_orbit_dlp_fast.rs`'s own direct `gf.mul`/`gf.sqr` calls
(inside `S3Solver`, `NormalBasis`, `invert()`, `build_index()`, `extract()`)
could be wrapped the same way in a second, instrumented copy — but a real
fraction of that binary's group-law cost runs inside
`crypto_lib::cryptanalysis::koblitz_fast_arith::FastBinaryCurve` (`add`,
`double`, `scalar_mul`, `points_with_x`), which is shared-library code that
calls its own internal `Gf2` field directly (`pub gf: Gf2` on
`FastBinaryCurve`, called as `self.gf.mul(...)`/`self.gf.sqr(...)`), not
through any counting wrapper this task's own copy could intercept without
either editing the shared library (which this task and AGENTS.md's "prefer a
counting wrapper over modifying the shared library" both direct against) or
duplicating `FastBinaryCurve`'s point-addition/doubling/scalar-multiplication
formulas inside the new file — a correctness risk for no clear benefit, given
callgrind already gives a complete answer.

This is not a negligible gap: at n=53 (44-bit subgroup order) a single
`FastBinaryCurve::scalar_mul` is ~44 doublings plus up to ~44 additions, each
of which is one `Gf2::inv` (Itoh–Tsujii, ~7–8 multiplications) plus 2–3 more
multiplications/squarings — call it ~10–12 real field operations per point
op, ~900 per scalar multiplication. The rank stage (K=440 attempts) and the
target stage (L=1,024 targets, each generated from its published scalar via
one `fast.scalar_mul`, exactly mirroring how the KS side also generates each
target from its scalar via one `raw_scalar_mul`, so this is charged on both
sides) both call `FastBinaryCurve::scalar_mul` at least once per unit, which
on its own is on the order of `(440 + 1,024) × 900 ≈ 1.3` million field
operations — large enough that omitting it would understate the IC side's
real cost and bias the ratio in the IC arm's favor, which is the wrong
direction to be wrong in for this specific recheck. So: exact call counting is
reported as a **supplementary diagnostic for the KS arm only** (it shows
concretely where the savings come from — canonicalization goes from up to
`2n` squarings to zero field operations, inversion goes from ~104 calls to
~8), and callgrind is the metric the pre-registered rule below is scored on.

**What callgrind's count does and does not charge:** it is a whole-process
instruction count on one host, one build (`--release`, this note's recorded
`rustc --version` and CPU features), single run per binary (not an A/A/B/A/B
interleave — see the host-contention note below for why a single paired run
is used here rather than the five-round interleave AGENTS.md §10 asks for
wall-clock comparisons, and why that is an accepted limitation for this
one-off recheck rather than a full campaign). It does not distinguish
field-arithmetic instructions from everything else (hashing, JSON
serialization, table lookups) the way the exact call count would for the KS
arm alone — that is the price of a metric that is complete on both sides.

## Pre-registered pass/fail rule (verbatim, fixed now, not renegotiated after seeing results)

> The compact-orbit IC result survives if matched-arithmetic rho's operation
> count divided by IC's operation count, at n=53 L=1,024, is below 0.8x; it
> dies if that ratio is 1.0x or above; report the exact ratio explicitly if it
> falls in [0.8, 1.0) or elsewhere, honestly, without forcing it into
> survive/die.

"Operation count" here is the callgrind retired-instruction count from the
section above. "Survives"/"dies" describes only this one n=53, L=1,024,
K=440 cell measured this way, on this host, in this one rerun — not the
n=41 cell, not any other L, and not independent replication (see Scope,
below).

**Inadmissible, stated in advance:** changing K, L, dp_bits, batch_seed or the
corpus name from the values above; choosing a favorable seed; skipping the
correctness gates; counting only one phase (rank stage only, or targets only)
instead of the whole process; and reporting wall clock as the metric the rule
is scored on (it is recorded for reference only, per AGENTS.md §6).

## Host and contention (recorded before running, AGENTS.md §10)

To be filled in immediately before the run, in this same commit's follow-up:
`uname -a`, `rustc --version`, CPU model/features, logical core count,
memory, `uptime`/`top` before and after each measured run, and the exact
build/run commands. If callgrind at L=1,024 does not finish in a reasonable
wall-clock budget, this note will be amended (not rewritten) to name the
smaller L actually used for the callgrind leg specifically, with the
Kuhn–Struik batching caveat: batched rho shares its jump table and
distinguished-point table across the whole batch, so its per-target cost
falls as L grows, and a smaller L is not a free stand-in for L=1,024 on that
arm. The wall-clock leg (§ below) is still run at the full frozen L=1,024
regardless.

## What this note does and does not authorize

This is an **accounting** check under AGENTS.md §3: correcting an unmatched
baseline, not proposing a new algorithm. If the matched rerun changes the
verdict, this note's classification is "accounting", not "advance" — no
index-calculus count changed, what the KS reference was measured against did.
This one rerun, on one host, does not itself promote or demote
`PENDING_INDEPENDENT_VALIDATION`; it states plainly that it is one rerun,
unreplicated on a second host/agent, same limitation the original measurement
had.
