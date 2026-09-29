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
