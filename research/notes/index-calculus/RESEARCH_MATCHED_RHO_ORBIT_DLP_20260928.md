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

## Result (2026-09-28, additive — nothing above this line was edited after seeing it)

**Host, recorded before running:**
`Linux vm 6.18.44-fc-v37 x86_64`, `rustc 1.94.1 (e408947bf 2026-03-25)`,
Intel(R) Xeon(R) Processor @ 2.10GHz, 4 logical cores, 15 GiB RAM, CPU flags
include `pclmulqdq aes avx2 bmi2 avx512f popcnt` (relevant to the carryless
multiply this whole comparison is about), no swap. `valgrind-3.22.0`. Repo
`HEAD` at run time: `32622a2107c01f493c7d361dab900be1c583167d`. Load average
stayed at 0.44–0.99 (`uptime`, recorded before/after every run below) — quiet,
no contention, well under the core count. Raw manifests, per-run load
snapshots, and every command's stdout/stderr are committed under
`research/notes/index-calculus/matched_rho_orbit_dlp_20260928_run/`, with
`SHA256SUMS` over both source files and both release binaries.

**Correctness, all four runs (unmodified KS corpus generation, matched-arith
KS native, matched-arith KS under callgrind, IC native, IC under callgrind):**
every one of 1,024 fixtures/targets verified on every run
(`all_verified`/`verified`/`group_verified` all `true`, `targets_solved: 1024`,
`targets_failed: 0`, `rank_failures: 0`, `rank: 440` = full rank at the
declared K). `total_walk_steps` is bit-identical — `19,103,507` — across the
unmodified-KS corpus run, the native matched-arith run, and the
callgrind-instrumented matched-arith run: the matched-arithmetic rewrite
provably walks the exact same rho trajectory (same jumps, same distinguished
points, same collisions), at lower cost per step, exactly as the correctness
gates required.

**Wall clock** (reference only, per AGENTS.md §6 — not what the pre-registered
rule is scored on), native runs, single process each, no other load on the
host:

| arm | in-process wall | vs IC | notes |
|---|---:|---:|---|
| KS, unmodified (`koblitz_rho_batch_ks`) | 215.240 s | 6.995× | Fermat inverse + squaring-walk canonicalization |
| KS, matched-arithmetic (`koblitz_rho_batch_ks_matched_arith`) | 30.442 s | 0.989× | this round's rewrite |
| IC (`koblitz_orbit_dlp_fast`, unmodified) | 30.772 s | 1.000× | K=440, construct mode |

Removing the arithmetic asymmetry alone is a **7.071×** wall-clock speedup for
the KS side (215.240 s → 30.442 s), and takes the KS/IC wall ratio from 6.995×
down to 0.989× — on wall clock, matched-arithmetic rho and IC are within 1.1%
of each other on this host, with matched rho very slightly faster.

**Operation count — valgrind `--tool=callgrind` retired instructions (`Ir`),
whole process, the metric the pre-registered rule is scored on:**

| arm | Ir (retired instructions) | source |
|---|---:|---|
| KS, matched-arithmetic | 371,102,176,689 | `callgrind_ks_matched.out` / `.annotate.txt` |
| IC | 89,164,459,930 | `callgrind_ic.out` / `.annotate.txt` |

`callgrind_annotate`'s self-cost breakdown is consistent with the mechanism
this round targets: on the matched-arithmetic KS side, `raw_canonicalize`
alone is 50.82% of all instructions and `raw_mul_field`/`raw_square` together
another 30.73% — canonicalization no longer calls the field at all, but
walking every one of up to `n = 53` orbit positions to find the least
representative is still `Θ(n)` **word** operations (two `Linear::apply` byte-
table lookups per position), just no longer `Θ(n)` **field multiplications**.
That is the real content of "O(1) rotation instead of a squaring chain": the
per-position cost drops from one full carryless-multiply-and-reduce to a
handful of table lookups, not to nothing, and canonicalization remains the
single largest cost center in absolute instructions. On the IC side, `extract`
(54.14%) and `main`'s rank/target loops (36.66%) dominate, with `Gf2::sqr` and
its `clmul_u64` carryless multiply together at 6.96% and `Gf2::inv` at 0.69%.

**Supplementary diagnostic, KS side only** (exact `raw_mul_field`/`raw_square`
call counts, not the metric the rule is scored on — see "why not exact Gf2
call counting on both sides" above): `total_mul_calls = 185,811,921`,
`total_sqr_calls = 82,757,971`, `total_field_ops = 268,569,892`,
`total_inv_calls = 20,645,769` (the last of these is a call count, not an
added operation — its own multiplications/squarings are already inside the
first two numbers). Per `raw_add`/`raw_double` call
(`group_additions = 20,384,452`), that is about 13.2 field operations and
about 1.01 field inversions on average, against the unmodified file's `~104`
per canonicalization call (up to `2×(n−1) = 104` squarings) plus another
`~104` (Fermat) per inversion — consistent with the double-digit wall-clock
speedup measured above.

### The pre-registered rule, applied

The task that dispatched this check transcribed the rule's ratio direction
backwards ("matched-arithmetic rho's operation count divided by IC's
operation count"). The user's own original wording, which this note treats
as authoritative, states every ratio in this thread the other way round —
challenger over reference, IC over rho: "It took 0.25× the time of batched
rho", "the ratio lands around 0.75–1.25" (both continuing that same IC/rho
convention), and "the lead survives below 0.8× and dies at 1.0× or above"
for "the lead", meaning the IC finding. That is also this repository's usual
convention elsewhere (`docs/ic/BOUNDARY_TARGETS.md`'s `S / S_ρ < 1`, PR
#830's own headline number, IC-wall / rho-wall = 0.246). The dispatching
task's own phrasing inverted this; that inversion is corrected here rather
than carried forward, since the two readings are exact opposites (dies vs.
survives) and only one matches what was actually asked for.

Ratio, IC over matched-arithmetic rho, in the requested (and repository-
standard) direction:

```
89,164,459,930 / 371,102,176,689 = 0.2403
```

`0.2403 < 0.8`, so by the pre-registered rule: **the lead survives**, on the
operation-count metric, at this one cell. IC still costs about a quarter of
matched rho's retired instructions — decisively below the survive threshold,
and close to PR #830's original unmatched 0.246 wall-clock figure, though
arrived at through a different mechanism this time (see below). For
completeness, and because the numbers should stand regardless of which
direction is read: matched-arithmetic rho costs 371,102,176,689 /
89,164,459,930 = **4.162×** IC's instruction count, i.e. IC is the cheaper
arm by that same margin stated the other way round.

At n = 53, L = 1,024, K = 440, with rho given the exact same normal-basis
Frobenius and Itoh–Tsujii inversion the IC arm uses, **IC still costs about a
quarter of matched rho's retired instructions on this host** — the
arithmetic-matching fix does not reverse the original claim; on this metric
and this host it leaves it intact, close to its original magnitude.

### A second thing worth flagging: wall clock and instruction count disagree here

On wall clock, matched-arithmetic rho and IC are within 1.1% of each other
(0.989×); on retired instructions, rho costs 4.16× more than IC. Both are
measured on the same host, same day, same binaries. The likely explanation is
memory behavior, not arithmetic: IC's index (`root_table_entries: 10,259,234`,
`regular_states: 10,260,800`, `peak_rss_bytes: 1,363,279,872` ≈ 1.36 GiB in
the callgrind run) is large enough to miss cache and page repeatedly, which
costs real wall-clock cycles per instruction without costing additional
retired instructions; matched rho's distinguished-point table stays small at
`dp_bits = 4` (`table_entries: 1,273,250`, a few tens of MB) and its inner loop
is comparatively cache-resident. This is exactly the reasoning AGENTS.md §6
gives for why "operation counts are the metric because they survive
hardware" and wall clock does not: on a host whose relative memory/compute
balance differs from the one PR #830 originally ran on, the wall-clock ratio
alone would have told a materially different story (near-parity) from the
instruction-count ratio (IC decisively cheaper). This rerun's own wall-clock
number is therefore a caution against reading PR #830's original 0.246
wall-clock figure as portable across hosts, even before the arithmetic
correction is considered.

### Classification and scope

**Class: accounting** (AGENTS.md §3) — the KS reference's arithmetic was
corrected to match the IC arm's; no index-calculus count changed, and this
round does not claim the fixed KS is a better rho than before in any sense
beyond removing the specific asymmetry named above (Fermat inversion,
squaring-walk canonicalization). It does not relax AGENTS.md §8's
end-to-end/whole-pipeline requirements, and it is silent on n = 41, on any
other L, on n = 61, and on ECC2K-130/m = 83 transfer.

`PENDING_INDEPENDENT_VALIDATION` is **not** changed by this note. This is one
rerun, on one host, by one agent, with no second host or independent replay —
the same limitation the original PR #830 measurement carried. What changes is
narrower and stated plainly: the specific arithmetic-asymmetry concern this
note was written to check has been measured, at the exact frozen n=53,
L=1,024, K=440 cell, and does not reverse the direction of the original
claim on this host, on either metric — even though the two metrics disagree
sharply on the *margin*.
