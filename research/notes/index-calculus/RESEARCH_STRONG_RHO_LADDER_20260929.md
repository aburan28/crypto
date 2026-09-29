# Strongest-available batched rho vs compact-orbit IC (n=53, L=1,024): reference-strengthening ladder

Written and committed **before any of rungs R1–R3 is built or run**, per
`AGENTS.md` §1 ("state the boundary before measuring") and §4 ("declare the
falsification target in advance"). Nothing below the line "Results" is a result.

## Erratum on PR #955 (this note supersedes its verdict where they conflict)

PR #955 (`RESEARCH_MATCHED_RHO_ORBIT_DLP_20260928.md`) gave rho the IC arm's
normal-basis rotation and Itoh–Tsujii inversion, called that rho "matched", and
reported IC/rho = 0.2403 in retired instructions, "the lead survives". Two
things were wrong with calling that reference matched, both visible in the
callgrind annotation PR #955 itself committed
(`matched_rho_orbit_dlp_20260928_run/callgrind_ks_matched.annotate.txt`) and
both missed by its author:

1. **Hardware-class asymmetry.** The IC arm multiplies with the library
   `Gf2` (`src/cryptanalysis/semaev_decomp.rs`): one hardware `pclmulqdq` plus a
   byte-table reduction. The rho file's `carryless_product` uses hardware
   carry-less multiply **only on aarch64**; on x86-64 it falls back to a
   bit-serial software loop, followed by a bit-serial `raw_reduce`. On the x86-64
   host of PR #955 that is ~429 retired instructions per `raw_mul_field` and ~415
   per `raw_square` (79.7 B and 34.4 B over 185.8 M and 82.8 M calls), against a
   handful for IC. PR #955's instruction ratio therefore favoured IC by a
   host-specific handicap on the rho side that the original PR #830 host
   (aarch64, where rho had PMULL) did not have. PR #955's *wall-clock* figures
   inherit the same problem.
2. **"Matched" is not "best".** `AGENTS.md` §1 defines the reference as "the
   best algorithm that already solves the same problem". After PR #955 the rho
   still spent 50.82 % of all instructions in `raw_canonicalize` (Θ(n) byte-table
   basis changes and a `u128 %` — `__umodti3`, 4.15 % — per orbit position, for
   every step), and one full Itoh–Tsujii inversion per step, although the
   repository's own `Gf2::batch_inv` (Montgomery's trick) is exactly the
   primitive high-performance rho implementations use to amortise inversions.

PR #955's arithmetic-asymmetry finding (rho's squaring-walk and Fermat inversion
were a large unmatched cost) stands. Its *verdict* — "IC/rho = 0.2403, the lead
survives" — is a statement about IC versus a rho that is neither
hardware-matched nor best-effort, and must not be quoted as a statement about
IC versus rho. This note replaces it with a ladder of progressively stronger
references.

## What was already seen before writing this (disclosure)

The ladder below was motivated by PR #955's committed callgrind breakdown of the
R0 rho (371,102,176,689 Ir): `raw_canonicalize` 50.82 %, `raw_mul_field` 21.47 %,
`raw_square` 9.26 %, `main` 9.25 %, `raw_inverse` 4.30 %, `__umodti3` 4.15 %. The
rungs are therefore not blind to the profile. What has **not** been seen: the
instruction count, step count or wall time of any rung R1–R3. The decision rule
below is fixed before any of them runs.

## Frozen cell (unchanged from PR #830 / PR #955)

n = 53, a = 0, quotient mode `signed_frobenius`, L = 1,024 fixtures, K = 440
(IC only), batch_seed = 531310, `KIC_RHO_DP_BITS=4`, corpus name
`n53-ks-growing-1024-v1`, rank_seed 7 (IC). Targets are a deterministic blake3
function of (corpus, batch_seed, index) and regenerate bit-identically for every
rung. Host: x86-64, `pclmulqdq`, 4 logical cores, single-threaded runs.

## The ladder (cumulative; each rung contains the ones before it)

| rung | change to the rho reference | algorithm changed? |
|:--|:--|:--|
| R0 | PR #955's `koblitz_rho_batch_ks_matched_arith` (normal-basis rotation, Itoh–Tsujii; software clmul on x86-64) | — (measured: 371,102,176,689 Ir) |
| R1 | every field multiply/square is the library `Gf2::mul`/`Gf2::sqr` — the same code IC runs (hardware clmul, table reduction) | no: identical field values, so the walk must be **bit-identical** to R0 |
| R2 | canonicalize by the least rotation of the x normal coordinates (integer rotate/compare), choose the sign by comparing y, convert the winner back once; orbit multiplier from a precomputed `λ^k` table; `mul_mod` by float-quotient estimate instead of `u128 %` | yes: a different (still class-invariant) orbit representative, so the walk changes; statistically equivalent |
| R3 | `W` walks advance in lockstep; their `x₁+x₂` denominators are inverted together with the library `Gf2::batch_inv` (one inversion + 3(W−1) multiplications per round) | no change to the per-walk map; scheduling only |

All rungs live in one new file, `examples/koblitz_rho_batch_ks_strong.rs`,
selected by `KIC_RHO_RUNG=0..3`. Rung 0 in that file is a regression check that
it still reproduces PR #955's walk; PR #955's own binary remains the R0 source of
truth for instruction counts. `koblitz_rho_batch_ks.rs`,
`koblitz_rho_batch_ks_matched_arith.rs` and `koblitz_orbit_dlp_fast.rs` are not
modified.

**Admissibility.** A rung may use only techniques from the public literature on
binary-field rho (normal-basis Frobenius, canonical orbit representatives,
Montgomery batch inversion); must solve the same 1,024 targets; must not use any
planted scalar except to verify at the end; and must charge everything it does
(setup, target generation, walks, table, verification) in the same whole-process
instruction count. Wall clock is recorded for reference only (`AGENTS.md` §6).

## Gates a rung must pass before its number is used

- **G1 (correctness).** All 1,024 targets recovered, each verified `[d]G = Q`
  and equal to the planted scalar; zero failures.
- **G2 (R1 only).** `total_walk_steps == 19,103,507` and every per-fixture
  `walk_steps` equal to R0's, since R1 changes no field value.
- **G3 (R2, R3).** Total walk steps within `[0.85, 1.25] ×` R0's 19,103,507 (a
  different representative or schedule may change the trajectory but must not
  change the walk's health), and mean steps per target consistent with
  `√(πr/(2A))` shrinking under batching as in R0.
- **G4 (unit tests).** R1: `Gf2` equals the raw field ops (existing test).
  R2: the new canonicalization is class-invariant (`canon(φᵏP) = canon(±P)` for
  all k and both signs), its multiplier `m` satisfies `canon.point = [m]P` by
  group-law replay, and `mul_mod` equals the `u128` reference on random and edge
  inputs. R3: end-to-end at n ∈ {17, 19, 23} with every target verified.

A rung that fails a gate is reported with the failure and is not eligible to be
"the best rho".

## Metric and decision rule (fixed now)

Metric: `valgrind --tool=callgrind` retired instructions (`Ir`), whole process,
for both arms, on the same host. The IC figure is **re-measured** on the current
source of `examples/koblitz_orbit_dlp_fast.rs` (sha256 `8dea682a…`, which
differs from the `2afd019a…` measured in PR #955 by an optional rank-trace dump
and one `row.clone()`), not reused.

Let `rho_best` be the smallest `Ir` among rungs R1, R2, R3 that pass their gates,
and `ρ* = Ir(IC) / Ir(rho_best)` (challenger over reference, the convention of
PR #830 and `docs/ic/BOUNDARY_TARGETS.md`; `< 1` means IC is cheaper).

- The lead **survives** if `ρ* < 0.8`.
- The lead **dies** at this cell if `ρ* ≥ 1.0`.
- In `[0.8, 1.0)` the result is reported as inconclusive, without being forced
  into survive/die.

**Inadmissible, stated in advance:** changing K, L, dp_bits, batch_seed or the
corpus; dropping a cost from either side; scoring on wall clock; selecting the
rung after seeing the numbers other than by the rule above (the smallest passing
`Ir`); and reporting a rung whose gates failed.

## Scope, asymmetry and classification

- This tests **IC as merged** against **the best rho the rungs above produce**.
  IC receives no comparable micro-optimisation (for example its S3 root solver
  inverts once per candidate pair). That is deliberate: this note is the
  adversarial test of the lead. If the lead dies here, the follow-up question is
  whether equally optimising IC restores it; that is *not* answered by this note
  and no "IC cannot win" conclusion may be drawn from it.
- One cell (n = 53, L = 1,024, K = 440), one host class (x86-64 with
  `pclmulqdq`), one implementation of each arm. Nothing here speaks to n = 41 or
  61, to other L, to aarch64, to m = 83 or to ECC2K-130 transfer
  (`AGENTS.md` §8a/§8b); those are separate cells.
- Class (`AGENTS.md` §3): **accounting** for R1 (same values, same walk, cost of
  the same operation on the same hardware); **engineering** for R2 and R3
  (reference algorithm improved; not an advance for either arm).

---

## Results (2026-09-29, additive — nothing above this line was edited after seeing it)

**Host and build.** `Linux vm 6.18.44-fc-v37 x86_64`, Intel Xeon @ 2.10 GHz, 4
logical cores, 15 GiB, flags include `pclmulqdq popcnt aes avx2 bmi2 avx512f`;
`rustc 1.94.1`, `valgrind-3.22.0`; repo `HEAD` `d7b3957b8`. Single-threaded runs.
`examples/koblitz_rho_batch_ks_strong.rs` sha256 `3a4eb7b4…`, its release binary
`988a3657…`; `examples/koblitz_orbit_dlp_fast.rs` sha256 `8dea682a…`, binary
`e7466f49…` (full hashes in `strong_rho_ladder_20260929_run/SHA256SUMS_*`,
`host_manifest.txt`). Instruction counts do not depend on host load; other
callgrind jobs were running on other cores during several runs (load 1–3 of 4).

**Gates.** All four rungs pass G1 (1,024/1,024 targets recovered, verified and
equal to the planted scalar). G2: rungs 0 and 1 walk the *bit-identical*
trajectory to PR #955's run — every one of the 1,024 per-fixture `walk_steps`
equal, total 19,103,507, same table size (1,273,250) and cross-target solves
(1,020) — so switching to the library `Gf2` changed no field value. G3: rung 2
19,247,253 steps (× 1.0075), rung 3 19,569,835 (× 1.0244; the excess is the
predicted work abandoned in the other lanes when a target completes), both
inside `[0.85, 1.25]`. G4: the six unit tests in the new file pass, including
the group-law replay of the rung-2 multiplier and the exhaustive-definition
check of the fast canonical form. The four rungs solved the identical 1,024
targets as PR #955. **Independent replay** (`independent_rho_replay.py`, the
repository's pure-Python GF(2ⁿ)/Koblitz arithmetic, no Rust): for every rung,
1,024/1,024 records have generator on the curve with order r, target on the
curve, `[d]G = Q`, recovered = published scalar, and `[r]Q = ∞` on a 64-record
sample (`independent_replay_rung{0..3}.json`).

**The ladder (whole-process retired instructions, `valgrind --tool=callgrind`;
generated by `strong_rho_ladder_20260929_run/ladder_table.py` from the raw logs).**

| reference rho | Ir | IC / rho | steps | verified |
|:--|--:|--:|--:|--:|
| R0 as merged in PR #955 (its own binary) | 371,102,176,689 | 0.240 | 19,103,507 | 1,024/1,024 |
| R0 regression: this binary, rho-local software field | 406,208,773,598 | 0.220 | 19,103,507 | 1,024/1,024 |
| R1 library `Gf2` field (hardware clmul, table reduction) | 281,027,207,019 | 0.317 | 19,103,507 | 1,024/1,024 |
| R2 + normal-coordinate canonicalization, fast `mul_mod` | 68,926,289,783 | 1.294 | 19,247,253 | 1,024/1,024 |
| R3 + 32 lockstep lanes, `Gf2::batch_inv` | 46,384,789,571 | **1.923** | 19,569,835 | 1,024/1,024 |
| IC `koblitz_orbit_dlp_fast` (current source) | 89,190,346,806 | — | — | 1,024/1,024 |

IC re-measured on the current source: 89,190,346,806 Ir, 1,024/1,024 solved,
rank 440, 0 failures — 0.03 % above PR #955's figure for the earlier source, so
the rank-trace change was indeed benign.

**Decision under the pre-registered rule.** `rho_best` is R3 (the smallest `Ir`
among gate-passing R1–R3), and `ρ* = 89,190,346,806 / 46,384,789,571 = 1.923 ≥ 1.0`.
**The lead dies at this cell**: against the strongest rho this repository can
build from the same field code, compact-orbit IC costs 1.92× the retired
instructions at n = 53, L = 1,024, K = 440. PR #955's "the lead survives,
IC/rho = 0.2403" is superseded; it compared IC with a rho that was neither
hardware-matched nor best-effort.

**Which change moved the verdict.** R1 alone (same field code as IC) takes the
reference from 371 B (PR #955's build) or 406 B (this binary's rung 0) to 281 B:
between a quarter and a third of the reference's cost, depending on which build
of rung 0 is the base, was the x86-64 software multiply; IC/rho is still 0.317
there. The flip happens at R2 (−212 B):
PR #955's "matched" rho still found the canonical orbit representative by a
Θ(n) scan of polynomial-basis conversions with a `u128 %` per orbit position,
which the IC arm does not do (it takes a least rotation of integers). R3 removes
a further 22.5 B by amortising inversions. So the survival in PR #955 rested on
an implementation choice in the reference, not on the algorithms.

**Where the remaining cost is.** R3 (2,370 Ir/step): `FastCanon::apply` 23.87 B
(51.5 %, 1,142 Ir per call over 20.9 M calls — still a Θ(n) least-rotation scan),
`run` 8.85 B (19.1 %), `Gf2::batch_inv` 4.36 B (9.4 %), `Field::inv` 2.18 B and
`Gf2::inv` 2.06 B (walk starts, verification, setup), `clmul_u64` 1.28 B, default
`SipHash` on the DP table 1.4 B (3.0 %). IC (89.19 B): `extract` 48.27 B (54.1 %),
`main` 32.71 B (36.7 %), `Gf2` operations ~7 B. **R3 is therefore not the best
rho that could be written** — a faster least-rotation, a cheaper DP-table hash
and a rotation-invariant step selector would each cut it further — and every such
cut raises ρ* above 1.92. The reported 1.92 is a lower bound on the constant IC
loses by, for this reference family.

**Build-to-build scale.** Rung 0 of the new binary costs 9.5 % more instructions
than PR #955's binary (406.2 B vs 371.1 B) for the *identical* walk. `run`'s
self-cost (223.93 B) is the same in R0 and R1 of the new binary and equals PR
#955's `raw_canonicalize` + `main` (188.6 + 34.3 B), so the difference is
inlining and layout of the software multiply, not work. About 10 % is thus the
scale of uncontrolled code-generation variation between two builds of the same
algorithm; R2's ρ* = 1.29 clears it, R3's 1.92 clears it by a wide margin.

**Wall clock (reference only, `AGENTS.md` §6).** Native in-process wall: R0
35.4 s, R1 22.9 s, R2 7.3 s, R3 5.35 s (some concurrent load). IC's native wall
at this cell was 30.8 s in PR #955's run. The Ir/second of the two arms differs
threefold (IC ≈ 2.9 G, R3 ≈ 8.7 G instructions per second): IC's 1.36 GiB index
runs at a lower IPC, which is consistent with the memory-bound reading in PR
#955 but was **not** tested here (no cache simulation was run), so it stays a
hypothesis. Wall clock and instructions point the same way at this cell: the IC
arm loses on both, by 1.9× in instructions and about 5.8× in wall time.

**Scope, limits, classification.** As pre-registered: one cell, one host class
(x86-64 with `pclmulqdq`), one implementation of each arm, IC as merged with no
comparable micro-optimisation, so this measures IC-as-merged against a
strengthened reference and does **not** show that no IC variant could win. It
does not address n = 41 or 61, other L, aarch64 (where PR #830 ran; there the
rho file already had PMULL, so the hardware-class artefact would be absent but
the canonicalization gap would remain), m = 83 or ECC2K-130 transfer. R1 is
class **accounting**; R2 and R3 are **engineering** on the reference. The
scaling question this cell cannot answer is taken up in
`RESEARCH_STRONG_RHO_SWEEP_PROTOCOL_20260929.md`.

**Superseded statements.** (i) `RESEARCH_MATCHED_RHO_ORBIT_DLP_20260928.md`'s
verdict ("the lead survives … 0.2403") and its use of the word "matched" for
the x86-64 rho; its arithmetic-asymmetry finding stands. (ii) Any reading of
PR #830's wall-clock 0.246 / 0.339 ratios as evidence that compact-orbit IC
beats batched rho: those compared IC with the unmodified rho; the strengthened
reference is cheaper than IC at n = 53, L = 1,024 in both instructions and wall
time. PR #830's own status (`PENDING_INDEPENDENT_VALIDATION`, not promoted)
already withheld the promotion this result now argues against.
