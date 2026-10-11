# lit2026 W11: generic-method engineering from five 2025–26 papers

**Workstream.** W11 of the lit2026 speedup program
([`../lit2026_speedup_program/`](../lit2026_speedup_program/README.md)), stacked on
aburan28/crypto#1681. **Date.** 2026-10-11. **Run.**
[`reports/run-20261011T2049Z/`](reports/run-20261011T2049Z/) (run report
[`report.pdf`](reports/run-20261011T2049Z/report.pdf)); results PDF
[`results.pdf`](results.pdf).

**Question.** Do five generic-method engineering ideas from the 2025–26 literature buy
anything on this repository's own code paths, measured against the program baseline
`4b6e09b78`?

| | paper | what we built | verdict (this run) |
|---|---|---|---|
| (a) | Cenk–ElSheikh–Keskinkurt Paksoy–Hasan, ePrint 2026/2102 | `bsgs_fast` step kernels `paired` and `coz` behind `BSGS_STEP` | **speedup** (paired): 1.15–1.67× wall vs original at 40–56 bits (L0); counted cost per point 16.75 → 3.71 M-equivalents |
| (b) | Banegas–Sarkis–Smith, ePrint 2026/704 | x-only binary isogeny evaluators (affine fraction, product, inversion-free projective) + one-inversion pair evaluation in the torsion walk behind `BTW_XEVAL` | step: 1.50× fewer `F_{q^8}` mults, half the inversions, 1.29× wall; **end-to-end class walk: no change** |
| (c) | Zheng–Lin–Zhao, ePrint 2025/1652 | cubical (biextension) ladder for the `k = 2` Tate pairing on `eds_mov`'s path | **speedup in counted ops**: 17–19% fewer `F_p` M-equivalents than denominator-free Miller at 32–52-bit `r`; identical reduced values |
| (d) | Panny et al., ePrint 2026/1823 | DP trail-length recording in `koblitz_strong_rho` | **diagnostic**: DP trails follow a competing-risk geometric law; fruitless-cycle abandonment truncates trails at dp = 12 (mean 0.43–0.48 of geometric) and costs 1.37–1.58× the ideal steps |
| (e) | Dallaire-Demers–Doyle–Foo, arXiv 2508.14011 | the `secpkk1` challenge ladder (Algorithm 3) in Rust, solved by `prime_fast::rho_parallel` | **diagnostic**: 720 solves, 30–60 bits, iterations/√(πn/4) = 1.077 ± 0.020 (mean ± s.e.); the paper's 8/12/16-bit cards use the opposite prime-search direction from its own Algorithm 3 |

None of this is an exponent claim, a deployed-curve claim, or an IC-versus-rho result.
The BSGS speedup is a constant-factor change to one solver; no ecbench method-level
comparison was run (the ecbench BSGS arms are the counted `generic.rs` walks, which do
not call `bsgs_fast`; see "Limits").

## Papers and access

- ePrint 2026/2102 and 2026/704: abstracts read; the PDFs returned HTTP 403 from this
  host, so formulas are our own derivations, cited as "after" the paper. The co-Z arm is
  our reading of the 2026/2102 abstract (Jacobian giant steps, co-Z paired giant steps,
  windowed batch inversion).
- ePrint 2025/1652: abstract read (ladder with y-recovery; 17.53% fewer `F_p` mults for
  the reduced Tate pairing on SS-1536). Our ladder is a three-lane cubical ladder with an
  `unsquare` step; the PDF was not read.
- ePrint 2026/1823: abstract read (distinguished-subgraph analysis under exploration
  modes). Item (d) applies its question to the repo's DP walk; it does not reproduce the
  paper's analysis.
- arXiv 2508.14011v2: PDF read (Algorithm 3, p. 10; Appendix A cards, pp. 25–28).

## Host, arms, isolation

- Host: Apple M4 Pro, 14 cores, macOS 26.6, shared; 1-minute load 72–155 during the run
  (recorded per sample). **Isolation L0** for every wall-time number.
- `original`: program baseline `4b6e09b78`, built from the read-only worktree
  `/Volumes/SSD990-2/crypto/wt/lit2026-baseline` (the BSGS timing example is not in the
  baseline tree, so it was compiled against the baseline library from a scratch crate;
  default release profile, same as the repo).
- `previous` / `candidate`: one binary from this branch with the switch off/on
  (`BSGS_STEP=chain` vs `paired`/`coz`, `BTW_XEVAL` unset vs `fraction`).
- Arms rotated per round (ABCD, BCDA, …); ≥ 5 rounds; intervals are min–max of per-round
  ratios. Deterministic counters (field operations, `F_{q^8}` operations) are the gate.

## Results

### (a) BSGS step kernels (`src/cryptanalysis/bsgs_fast.rs`)

Kernels (`BSGS_STEP`, window `BSGS_PAIR_K`):

- `chain` (original): one affine addition per point, `W = 8` chains share an inversion.
- `paired`: affine centre `C`, offsets `jA` (`j ≤ K`), `C ± jA` share
  `1/(x(jA) − x(C))`, non-centre points are x-only; one batch inversion over
  `W·(K+1)` denominators (windowed, memory `O(W·K)`).
- `coz`: Jacobian centre advanced by mixed additions, co-Z conjugate pairs after lifting
  `jA` to the centre's `Z`, the window's `Z`s inverted together.

Counted per visited point (giant phase, `W = 8`, single-threaded counting field;
inversion = Fermat, counted as `log₂p` S + popcount M; identical `(index, x)` visits to
the chain in every case):

| curve | kernel | K | M | S | I | M-equivalent |
|---|---|---|---|---|---|---|
| secp48k1 | chain (original/previous) | – | 4.625 | 1.000 | 0.1250 | 16.75 |
| secp48k1 | paired | 32 | 2.539 | 1.000 | 0.00193 | 3.71 |
| secp48k1 | coz | 32 | 5.117 | 2.054 | 0.00193 | 7.34 |
| secp56k1 | chain | – | 4.625 | 1.000 | 0.1250 | 18.88 |
| secp56k1 | paired | 32 | 2.539 | 1.000 | 0.00193 | 3.74 |

Without the inversion share the paired step is 3.54 vs 5.63 M+S (1.59×); the rest of the
4.5× in M-equivalents is the larger inversion window. The co-Z arm costs about twice the
paired arm's multiplications: with affine inputs, co-Z conjugate addition plus
normalisation does strictly more work than the affine shared-denominator pair.

Wall time, full-group solve with negation map (`BsgsFast::build` + `solve_many`),
`rayon` 14 threads, curves from item (e)'s ladder (generator = "Challenge" point),
targets 64/32/16/8/4; 5 interleaved rounds; seconds per batch:

| curve (`n`) | original | previous | paired K=32 | coz K=32 | original/paired | previous/paired | original/coz | original/previous |
|---|---|---|---|---|---|---|---|---|
| secp40k1 (< 2^40) | 4.67 | 4.28 | 3.20 | 3.49 | 1.459 [1.054, 1.917] | 1.337 [1.223, 2.037] | 1.337 [0.801, 1.526] | 1.091 [0.655, 1.091] |
| secp44k1 (< 2^45) | 5.63 | 6.12 | 4.27 | 4.22 | 1.316 [1.015, 2.002] | 1.430 [1.262, 1.737] | 1.333 [1.151, 1.834] | 0.920 [0.804, 1.399] |
| secp48k1 (< 2^48) | 8.63 | 9.55 | 6.64 | 8.16 | 1.301 [1.299, 1.603] | 1.439 [1.203, 1.849] | 1.058 [1.056, 1.814] | 0.904 [0.847, 1.157] |
| secp52k1 (< 2^53) | 15.55 | 19.40 | 13.58 | 14.77 | 1.145 [0.996, 1.976] | 1.429 [1.204, 2.466] | 1.052 [0.969, 1.917] | 0.801 [0.726, 1.002] |
| secp56k1 (< 2^57) | 45.84 | 48.99 | 27.44 | 35.32 | 1.670 [1.349, 1.670] | 1.785 [1.190, 1.849] | 1.298 [1.219, 1.751] | 0.936 [0.860, 1.327] |

Ratios are ratio-of-medians with [min, max] of per-round ratios. Verified answers and
table sizes are identical across all four arms at every size. The `previous` arm (this
branch, `BSGS_STEP=chain`) is slower than `original` at four of five sizes (0.80–0.94×
ratio of medians, intervals reaching above 1 at all of them); the chain loop was made
generic over a field trait, which should compile to the same code, and the difference
is within this host's L0 spread, but it is recorded rather than explained away.

### (b) x-only binary isogeny evaluation

**sect163-sized walk** (`F_{2^163}`, `z^163 + z^7 + z^6 + z^3 + 1`,
`y² + xy = x³ + x² + b`, `b` the first scan value with an x-rational ℓ-kernel; 8 steps
per ℓ, two carried abscissae; every evaluator equals `velu_x_map`, the Vélu point map's
abscissa, and the image lies on the `velu_codomain` curve). Counted per step (two
points), `F_{2^163}` mul/inv:

| ℓ (s) | batch inverse | kernel poly (h given) | fraction | fraction, 1 inversion | product (normalised) | projective |
|---|---|---|---|---|---|---|
| 3 (1) | 2 / 2 | 6 / 2 | 4 / 2 | 7 / 1 | 8 / 2 | 14 / 0 |
| 5 (2) | 8 / 2 | 8 / 2 | 8 / 2 | 11 / 1 | 14 / 2 | 20 / 0 |
| 7 (3) | 14 / 2 | 14 / 2 | 12 / 2 | 15 / 1 | 20 / 2 | 26 / 0 |

At these `s` the inversion dominates (one `flt_inverse` ≈ 162 squarings + ~10 mults), so
time follows inversions: at ℓ = 7 the median per step was 132.8 µs for the original
`velu_x_map`, 72.7 µs for the one-inversion fraction form (1.98× [0.91, 2.65]) and
1.27 µs for the projective form (48× [9.7, 222]), whose outputs stay projective across
steps. The product form `x·Π((x x_i + √b)/(x + x_i))²` is `(Π x_i)²·φ(x)` (a scaled
model; checked on every kernel) and costs `3s` M against the fraction form's `2s` M.
The paper's 5s/4s counts could not be compared formula-for-formula (PDF unavailable).

**Torsion-walk step** (`(n, a₂) = (19, 1)`, ℓ = 457, `E[ℓ] ⊂ E(F_{q^8})`,
s = 228, all 458 kernels, x(P) and x(Q) per kernel): `velu_x_from_abscissae` twice vs
`velu_x_pair_fraction` (`BTW_XEVAL=fraction`): 627 460 → 419 070 `F_{q^8}` mults
(1.497×), 916 → 458 inversions, identical outputs on all 458 kernels; wall per kernel
1564 → 1211 µs (ratio of medians 1.29, per-round 0.94–1.93, 9 rounds).

**End-to-end class walk** (`koblitz_isogeny_class_walk 19:1`, 12 interleaved rounds,
original binary from `4b6e09b78`): snapshots identical across all 36 runs (with the
`seconds` field removed). Walk seconds (reported by the example): original 1.280, previous
1.164, candidate 1.198 (medians); original/candidate 1.069 [0.714, 1.831],
previous/candidate 0.972 [0.561, 1.498]. **No change** detectable at L0: the evaluation
step is not where this walk spends its time.

### (c) k = 2 pairing: cubical ladder vs Miller (`src/cryptanalysis/biext_k2.rs`)

The repo's `k = 2` path is `eds_mov` (supersingular `y² = x³ + x`, `p ≡ 3 mod 4`,
distortion map). Three loops on one counting field, final-exponentiated: the affine Miller
loop of `eds_mov::miller_fp2` (ported), a denominator-free Jacobian Miller loop, and a
three-lane cubical ladder `(nP, (n+1)P, nP+Q)` with the cubical lift of `P − Q`. All three
equal `eds_mov::modified_tate` at `(p, r) = (1283, 107)`, and each other at every size;
**the cubical value is the square of the pairing** (measured on every instance), so
`unsquare` raises it to `(r+1)/2` (counted). Instances: `r` the largest prime below
`2^k`, `p = c·r − 1` the first prime with `4 | c`.

| r bits | p bits | affine Miller M+S / I | Jacobian Miller M+S | cubical M+S / I | unsquare M+S | (cubical+unsquare)/Jacobian, M+S | same, M-equivalent |
|---|---|---|---|---|---|---|---|
| 16 | 20 | 438 / 52 | 506 | 363 / 1 | 64 | 0.844 | 0.909 |
| 32 | 38 | 996 / 120 | 1148 | 737 / 1 | 147 | 0.770 | 0.830 |
| 40 | 46 | 1215 / 146 | 1401 | 912 / 1 | 181 | 0.780 | 0.840 |
| 52 | 58 | 1611 / 194 | 1857 | 1188 / 1 | 241 | 0.770 | 0.828 |

So 16–19% fewer `F_p` M-equivalents than the Jacobian Miller loop at 28–52-bit `r`
(the paper reports 17.53% on SS-1536). Wall times (7 rounds) are dominated by load
spikes and are recorded in `raw/biext_k2.jsonl` without a claim.

### (d) Distinguished-point statistics (`koblitz_strong_rho`, signed Frobenius)

Strong reference parameters (32 lanes, fruitless walks abandoned), `K₀/F_{2^53}`
(`r` = 21044858204113, 64 solves per row) and `K₀/F_{2^61}` (`r` = 162888033982417,
32 solves per row); all logarithms verified.

| curve | dp | steps | DPs | hits/step ÷ θ | mean trail / geometric | mean trail / competing-risk | abandon hazard /step | KS geo | KS c-r | KS 95% | steps / ideal |
|---|---|---|---|---|---|---|---|---|---|---|---|
| K0/F2^53 | 4 | 29.2 M | 1 945 590 | 1.066 | 0.993 | 0.997 | 2.97e-4 | 0.0023 | 0.0006 | 0.0010 | 0.82 |
| K0/F2^53 | 8 | 42.3 M | 165 190 | 1.001 | 0.921 | 0.991 | 2.96e-4 | 0.031 | 0.0042 | 0.0033 | 1.18 |
| K0/F2^53 | 12 | 56.4 M | 13 518 | 0.982 | 0.433 | 0.950 | 2.91e-4 | 0.301 | 0.022 | 0.012 | 1.58 |
| K0/F2^61 | 4 | 37.0 M | 2 465 006 | 1.066 | 0.994 | 0.998 | 2.55e-4 | 0.0019 | 0.0005 | 0.0009 | 0.80 |
| K0/F2^61 | 8 | 46.8 M | 183 269 | 1.003 | 0.932 | 0.993 | 2.56e-4 | 0.026 | 0.0025 | 0.0032 | 1.01 |
| K0/F2^61 | 12 | 63.4 M | 15 276 | 0.986 | 0.482 | 0.983 | 2.54e-4 | 0.264 | 0.0088 | 0.011 | 1.37 |

Deviations from the plain geometric model: (1) at dp = 4 the hit rate is 6.6% above θ
because a walk whose start point is distinguished counts a DP with zero steps; (2) the
abandonment hazard (fruitless cycles, ~2.5–3.0e-4 per step) truncates trails, so the
DP-trail law is the competing-risk geometric `(1−θ)(1−f)` (KS distance within or near the
95% bound) rather than `1−θ` (KS 0.26–0.30 at dp = 12); (3) the steps spent on abandoned
walks show up as steps/ideal rising from 0.80 (dp = 4) to 1.37–1.58 (dp = 12) — with
abandonment on, larger `dp_bits` costs work; the walk's `escape_fruitless` option exists
for that regime.

### (e) The arXiv 2508.14011 challenge ladder

`src/cryptanalysis/challenge_ladder.rs` derives `secpkk1` for `k ≤ 62` (CM order by
Cornacchia, twist chosen by annihilation, SHA-256 NUMS points). Checks against the
paper's Appendix A: the 6-, 24-, 32- and 48-bit cards are reproduced exactly by
Algorithm 3 as printed (descending prime search); **the 8-, 12- and 16-bit cards
(`p` = 163, 2089, 32803) are the first hits of an ascending search from `2^{k−1}`**,
which the descending search does not give (`p` = 211, …, 65419). Both are pinned in the
tests. The 30–60-bit levels used here are Algorithm 3 as printed (`raw/ladder_derive.jsonl`).

`prime_fast::rho_parallel` (negation map, r-adding walk, DP table), 64 seeds per level at
30–48 bits, 16 at 50–56, 8 at 58–60 (720 solves, every `κ` verified by `κQ = P`):
iterations/√(πn/4) pooled mean 1.077 (s.d. 0.546, s.e. 0.020), median 0.998; per-level
means 0.80–1.24 with per-level s.e. 0.06–0.19 (`raw/summary.json`). The configuration's
DP start-up overhead estimate `walkers·2^dp / √(πn/4)` is 0.12–0.18.

## Code (all behind switches; defaults unchanged)

- `src/cryptanalysis/bsgs_fast.rs`: `StepKernel` (`BSGS_STEP`, `BSGS_PAIR_K`),
  `FieldOps`, `CountingField`, `count_kernel_ops`; tests pin identical visited points,
  identical solves at interval seams, degenerate rounds, and the counted costs.
- `src/cryptanalysis/binary_velu.rs` Part E: `xeval_*` evaluators, `F2mCtx` counting,
  `rational_kernels_wide`, `y_for_x_odd`; test against Vélu at 163 bits, ℓ = 3, 5, 7.
- `src/cryptanalysis/binary_torsion_walk.rs`: `velu_x_pair_fraction`, `BTW_XEVAL`,
  per-thread `F_{q^m}` op counters; test of identity and inversion count.
- `src/cryptanalysis/biext_k2.rs`: counted `F_p`/`F_{p²}`, three pairing loops,
  `unsquare`; test against `eds_mov::modified_tate`.
- `src/cryptanalysis/koblitz_strong_rho.rs`: `dp_trail_lengths`,
  `abandoned_trail_lengths` in the outcome.
- `src/cryptanalysis/challenge_ladder.rs`: Algorithm 3 and the published-card test.
- Examples `examples/lit2026_w11_*.rs` (timing, op counts, binary, biext, DP stats,
  ladder, summary statistics).

## Commands

```sh
export CARGO_TARGET_DIR=/Volumes/SSD990-2/crypto/targets/lit2026-w11
cargo test --release --lib -- bsgs_fast challenge_ladder biext_k2 binary_velu::tests::xonly \
  binary_torsion_walk::tests::pair_fraction koblitz_strong_rho        # 31 passed
cargo build --release --example lit2026_w11_{bsgs_time,bsgs_opcount,binary_xonly,biext_k2,dp_stats,challenge_ladder,summarize} \
  --example koblitz_isogeny_class_walk
R=research/lit2026_w11_generic_engineering_20261011; RUN=$R/reports/run-20261011T2049Z/raw
lit2026_w11_challenge_ladder derive 30 60 2 > $RUN/ladder_derive.jsonl
sh $R/scripts/run_opcount.sh <bin>/lit2026_w11_bsgs_opcount $RUN/bsgs_curves.jsonl $RUN/bsgs_opcount.jsonl 40 48 56
sh $R/scripts/run_bsgs_arms.sh $RUN/bsgs_time.jsonl 5 <baseline>/lit2026_w11_bsgs_time <bin>/lit2026_w11_bsgs_time 32 $RUN/bsgs_curves.jsonl
sh $R/scripts/run_classwalk_arms.sh <scratch>/classwalk 12 <baseline>/koblitz_isogeny_class_walk <bin>/koblitz_isogeny_class_walk .
sh $R/scripts/run_rest.sh <bin> $RUN          # binary, biext, DP statistics
sh $R/scripts/run_ladder.sh <bin> $RUN/ladder_solve.jsonl
lit2026_w11_summarize $R/reports/run-20261011T2049Z <scratch>/classwalk
python -I $R/scripts/plot.py $R/reports/run-20261011T2049Z   # plots recorded values only
```

The baseline BSGS binary: a scratch crate depending on the baseline worktree
(`crypto = { path = "/Volumes/SSD990-2/crypto/wt/lit2026-baseline" }`, its
`Cargo.lock`), containing `examples/lit2026_w11_bsgs_time.rs` unchanged; the baseline
class-walk binary: `cargo build --release --example koblitz_isogeny_class_walk` in that
worktree with `CARGO_TARGET_DIR=/Volumes/SSD990-2/crypto/targets/lit2026-baseline`.

## Limits

- L0 only; load 72–155 on 14 cores. Wall ratios carry wide min–max intervals.
- (a) is not an ecbench method-level comparison: ecbench's BSGS arms are the counted
  `ecbench/generic.rs` walks, deliberately not `bsgs_fast`, so the reference-arm session
  was not run. The rho walk (`prime_fast`) does not share `bsgs_fast`'s step, so the
  kernels were not applied to it. K = 32 was fixed in advance, not tuned on the timing
  session.
- (b) The sect163-field kernels are x-rational kernels on curves found by scanning `b`
  (the first hits were scan values 1, 3 and 1 for ℓ = 3, 5, 7), not a registered sect163
  curve (whose rational ℓ-torsion was not examined); `s ≤ 3`. The paper's own formulas
  were not available.
- (c) The cubical ladder returns the squared pairing; the 17–19% includes the
  `unsquare` correction. Toy-to-58-bit `p` only (`u64` field).
- (d) Two Koblitz curves, one walk configuration (32 lanes, abandonment on).
- (e) 720 solves; per-level means at 50–60 bits rest on 8–16 seeds.
- Preliminary runs (fewer seeds/reps, and a biext run that stopped at an over-range
  `r_bits = 56` level) are kept in `raw/prelim/` and not used in the tables.

## Graphs checked and left unchanged

`docs/ic/progress-timeline.json`, `docs/index-calculus-scoreboard.html` and
`docs/performance-gains/` track index-calculus and Koblitz IC throughput; none of the five
items changes an IC stage, so they were left unchanged. The program's
`journey/journey.svg` and `board.svg` are re-rendered with this workstream's entries.
