# Isogeny-loop GLV, round 3: constructed curves, point-level validation, timings, binary GLS

**Date.** 2026-10-12.  **Question.** Round 2 (`../isogeny_loop_glv_20261011/`, PR #1698)
left one constructive lead and five open items: the 4-dimensional GLV-GLS split
(ψ plus a moderate-D loop endomorphism on a quadratic GLS twist) was measured
only on lattices and in an operation-count model. Does the derived eigenvalue
hold on points of a real curve? What does the split measure in time at
D = −619? Does a binary-field version exist? And what do the native sweeper
port and an x-only re-run of the searches say?
**PDF.** [`REPORT.pdf`](REPORT.pdf) (source [`REPORT.typ`](REPORT.typ)). The
figures are in [`figures/`](figures/), as SVGs with the plotting-only script
`plot_figures.py`, which reads the JSON records below.
**Run reports.** See [`reports/`](reports/), one directory per run (§7).

Status words: **measured** (a native run recorded it), **derived** (follows by
the stated argument), **proposed** (untested). All computation is native Rust in
`isogeny_algos/`. Python only plots recorded values and renders the PDFs.

## 0. Summary

| round-2 item | round-3 status |
|:--|:--|
| point-level validation of the GLV-GLS eigenvalue | **measured, all pass.** Curves were built from H_D at D ∈ {−619, −23, −71, −431, −3019}, with p of 127 bits. ψ acts as +μ, and each closing loop acts as ±λ(±s), on 24 points per check (§2). |
| benchmark at D = −619 | **measured.** GLV-GLS-4 with a projective loop takes 0.872× the time of GLS-2 (median/median, 1.147× faster) and 0.498× the time of plain wNAF. The A/A control reads 0.987. The measured gain is below the round-2 model at every D tested (§3). |
| binary-field GLS version | **constructed and validated at toy size.** The construction uses Frobenius-assisted loops 𝔩·F^k. Scans at m = 83, 127 and 131 find loops on prime-order GLS twists. An x-only point check at m = 31 over GF(2^62) passes on 2 discriminants (§4). Not timed. |
| native sweeper port | §5 (subagent run; see there for status) |
| x-only re-run of the searches | §5 |

## 1. What was built (code)

- `isogeny_algos/src/cm.rs`: the Hilbert class polynomial H_D by the complex-analytic
  method. It evaluates j(τ) = E4³/Δ for each reduced form, computing 1/q directly
  as exp(π√|D|/a) and the η-product by the pentagonal series. Arithmetic is
  fixed point on the crate's `Int`. Precision grows until every product
  coefficient rounds with error below 2^-16, and that bound is asserted.
  Tests check the class numbers and H_D for D = −3, −4, −7, −8, −11, −15, −23 and −163.
- `isogeny_algos/src/gls_loop.rs`:
  - certified cheapest loops (exhaustive up to the norm bound 2^{C/κ});
  - instance search, 4p = t² + |D|v² with r = #E′(F_{p²})/h prime;
  - CM curve construction from the roots of H_D mod p;
  - crater walks with Vélu isogenies (`api::isogenies`), keeping only crater
    j-invariants and forbidding backtracking;
  - the closing isomorphism, and evaluation of the twisted chain on E′;
  - 2/4-dimensional lattice reduction and Babai rounding;
  - 2-adic helpers for the binary case.

  Tests cover the cheapest loops for −23, −71, −619 and −3019, and a loop that
  closes and maps points correctly on a 40-bit D = −23 instance.
- `isogeny_algos/examples/gls_loop_validate.rs`: §2.
- `isogeny_algos/examples/gls_loop_bench.rs`: §3.
- `isogeny_algos/examples/binary_gls_scan.rs` and `binary_gls_validate.rs`: §4.

## 2. Point-level validation (measured)

`gls_validate.{md,json,log}`; command in the JSON. p is 127 bits and the seed is 20261012.

| D | loop (steps) | h(D) | H_D precision / worst rounding | F_p check (r_E bits) | ψ | twisted loops act as | verdict |
|--:|:--|--:|:--|:--|:--|:--|:--|
| −619 | 5·5·7 | 5 | 303 bits / 2^-90 | −λ(s), −λ(−s) (61) | +μ 24/24 | −λ(s); −λ(−s) | PASS |
| −23 | 2·3 | 3 | 150 / 2^-86 | (74) | +μ 24/24 | −λ(−s); −λ(s) | PASS |
| −71 | 2·3·3 | 7 | 238 / 2^-92 | (62) | +μ 24/24 | −λ(−s); −λ(s) | PASS |
| −431 | 2·2·3·3·3 | 21 | 631 / 2^-118 | (66) | +μ 24/24 | +λ(−s); +λ(s) | PASS |
| −3019 | 5·5·5·13 | 7 | 523 / 2^-88 | (63) | +μ 24/24 | −λ(−s); +λ(s) | PASS |

What each check means:
- E/F_p is built from a root of H_D mod p, then twisted if its trace is −t.
- The walk follows the certified degree sequence on the crater and returns to
  j₀. The walk found exactly 2 closing walks for each D, consistent with one
  per orientation (𝔞 and 𝔞̄).
- On E(F_p)[r_E], where r_E is a large prime factor of #E(F_p), each closing
  loop acts as one of ±λ(±s) mod r_E.
- On E′(F_{p²})[r] with r of 248–253 bits:
  - ψ = φ∘π_p∘φ⁻¹ equals [+μ] with μ = (p−1)/t mod r. This is measured, not just
    up to sign.
  - The twisted loop φ∘α∘φ⁻¹ equals one of [±λ(s)], [±λ(−s)], with
    λ(s) = (A + Bs)/2 and s = (2μ − t)/v. The two closing walks always land on
    *conjugate* embeddings, as the ring map π ↦ μ predicts.
- The sign ± is the [−1] automorphism and does not matter for decomposition.
  The benchmark fixes μ and λ, sign included, on its first sample point and asserts them on every pair.

This closes round 2's "not established": the derived eigenvalue holds on points.

## 3. Timing (measured)

`bench_run1.*` (D = −619) and `bench_D{23,71,431,3019}.*`, each with its `_env.txt`
(load averages before and after, CPU) and `_start_utc.txt`.

**Setup.**
- A single binary times every arm on the same arithmetic: `Fp2<FpM<2>>`
  (Montgomery, 2×64-bit limbs, Karatsuba), Jacobian `dbl-2007-bl`, mixed
  addition, w = 5 tables normalised with one batched inversion per call.
- The prime p is generic (127 bits, no special form). All absolute times
  therefore belong to this crate's arithmetic, and only ratios are reported.
- 256 (k, P) pairs, 30 rounds, arm order alternating per round.
- `taskset -c 3`, host Intel Xeon 2.10 GHz (cloud container, not isolated).
  The load average at start was 0.02–0.39.
- Before timing, every arm is checked on every pair: all arms agree, ψ(P) = [μ]P,
  α′(P) = [λ]P, and each decomposition ≡ k mod r. All checks passed in every run.

D = −619 (r 251 bits). Times are ns per scalar multiplication over 30 rounds:

| arm | median | min | wNAF / arm | GLS-2 / arm |
|:--|--:|--:|--:|--:|
| wNAF (w = 5) | 142,773 | 124,949 | 1.000 | 0.572 |
| GLS-2 | 81,665 | 71,901 | 1.748 | 1.000 |
| GLV-GLS-4, affine loop | 73,066 | 67,012 | 1.954 | 1.118 |
| GLV-GLS-4, projective loop | 71,172 | 64,788 | 2.006 | 1.147 |
| wNAF again (A/A control) | 144,657 | 125,440 | 0.987 | — |

Stage diagnostics at D = −619 (median ns):
- 4-dimensional decomposition: 2,850
- α′(P), affine: 7,776
- α′(P), projective: 5,602

Maximum coefficient sizes: 125 bits for GLS-2 and 63 bits for GLV-GLS-4.

Across D (GLS-2 median / GLV-GLS-4 median; figure `figures/bench_trend.svg`):

| D | loop cost (model M) | measured, projective loop | measured, affine loop | round-2 model (derived) |
|--:|--:|--:|--:|--:|
| −23 | 37.5 | 1.214 | 1.155 | 1.316 |
| −71 | 60 | 1.120 | 1.049 | 1.292 |
| −431 | 97.5 | 1.115 | 1.019 | 1.248 |
| −619 | 127.5 | 1.147 | 1.118 | 1.220 |
| −3019 | 210 | 1.062 | 1.043 | 1.158 |

**Observations.**
- GLV-GLS-4 was faster than GLS-2 at every D tested. The A/A control ranges
  0.987–1.019, so these are single-host, single-build measurements, and a
  difference under about 2% is within noise.
- The measured gain is 0.05–0.17 below the round-2 operation-count model.
  Two costs the model omitted are visible in the stage timings:
  - the decomposition (about 2.8 µs, roughly 3.5% of a GLS-2 call);
  - α′(P): 1.7–8.5 µs projective, which is 2–10.5% of a GLS-2 call. It is
    charged in full per call, together with its table and the ψ images. The
    model counted the loop as Σ 7.5ℓ M and charged nothing for the
    decomposition or the closing isomorphism.
- The measured α′(P) stage time rises with the model loop cost: 1.7, 2.7, 4.2, 5.6
  and 8.5 µs for D = −23, −71, −431, −619 and −3019. The end-to-end ratio is
  not monotone, because −619 reads above −71 and −431. That reversal is within
  run-to-run spread: the GLS-2 medians alone range 79.1–82.2 µs across the five
  runs. Pairwise orderings between D under about 3% are not resolved.
- Scope:
  - Variable-time code.
  - Generic p: a special p would cut every arm's field cost and leave the
    loop's fixed costs relatively larger.
  - One build on one host.
  - No constant-time variant.
  - No comparison with an optimised library.

## 4. Binary-field GLS with Frobenius-assisted loops

**Setting (derived).** Take E/F_q, q = 2^m, ordinary, with Frobenius
π = (t + v√D)/2 and t odd. 2 splits in O_K exactly when D ≡ 1 (mod 8). π
generates 𝔭₂^m, where 𝔭₂ is the prime above 2 that F (the 2-power Frobenius
isogeny) corresponds to on the crater. The quadratic twist E′/F_{q²} carries
ψ with ψ² = −1, as in the prime case, and #E′ = (q−1)² + t² = 2r.

In characteristic 2 an F step costs one squaring per coordinate. A loop may
therefore combine an odd ℓ-step 𝔩 with k F-steps whenever [𝔩][𝔭₂]^k = 1. Its
degree is ℓ·2^kk with kk = min(k, m−k), and its cost is c(ℓ) plus about 2kk
squarings. The direction is fixed 2-adically: pick s₂ with
v₂(t + v·s₂) > m, so that α = (A + B√D)/2 is F-oriented if v₂(A + B·s₂) > kk.
Round 2's erratum still applies: a V step costs a doubling. Only F steps are
cheap, which is why the scan keeps F-oriented α only.

**Scans (measured; `binary_scan_m{83,127,131}.{md,json,log}`).** The scan
requires all of the following:
- fundamental D ≡ 1 mod 8;
- ord[𝔭₂] = m;
- an odd t from Cornacchia (4q = t² + |D|v²);
- r = ((q−1)² + t²)/2 prime;
- the cheapest F-oriented loop with ℓ ≤ 61, ranked by model cost 7.5ℓ + 2kk M.

Each discriminant is then also measured on 5,000 scalars for its GLS-2 and
GLV-GLS-4 lattice bits.

| m | |D| range | scanned | ord[𝔭₂] = m | prime r | with a loop, ℓ ≤ 61 | cheapest (D: ℓ, F steps, model M) | GLS-2 / GLV-GLS-4 max bits |
|--:|:--|--:|--:|--:|--:|:--|:--|
| 83 | ≤ 8·10⁶ | 810,576 | 440 | 4 | 0 | — | — |
| 127 | ≤ 4·10⁶ | 405,292 | 319 | 3 | 2 | −1532327: 3, 53, 128.5; −2854663: 43, 21, 364.5 | 126 / 64 |
| 131 | ≤ 8·10⁶ | 810,576 | 459 | 4 | 3 | −283583: 13, 16, 129.5; −359383: 7, 60, 172.5; −1148567: 29, 59, 335.5 | 130 / 66 |

m = 127 is the field of the GLS254 curves of Oliveira–López–Aranha–Rodríguez-Henríquez
(F_{2^254}, GLS-2). m = 131 is ECC2K-130's base field, but these are *CM curves of
moderate D twisted over F_{2^262}*, not ECC2K-130. The challenge curve itself is
a D = −7 Koblitz curve, and round 1 covered it.

**Point-level validation (measured; `binary_validate.{md,json,log}`).** H_D has
degree h ≥ 127 at m = 127, and 1397 for the m = 127 hit, so this run is a toy
check at m = 31: curves over F_{2^31}, twisted over GF(2^62).
- E is built from a root of H_D mod 2, with H_D computed natively at
  precision 4,162 bits for h = 62 and 14,081 bits for h = 186.
- The odd ℓ-isogenies come from kernel polynomials (Kohel), and the walk is the
  ℓ-step then kk F-steps back to j₀.
- On points of prime order r (61 bits) of E′, the check compares
  x(ψ(P)) = x^q with x([μ]P), and x(α′(P)) = (R(x))^{2^kk} with x([λ(±s)]P).
- Twist and isomorphisms leave x unchanged in characteristic 2, so this is an
  x-only check: it validates λ up to sign, which is all a decomposition needs.

| D | h | loop | ψ = [μ] (x-only) | loop acts as | verdict |
|--:|--:|:--|:--|:--|:--|
| −39583 | 62 | ℓ = 31, 11 F steps, α = (−463 − √D)/2 | 24/24 | ±λ(s) 24/24 | PASS |
| −482623 | 186 | ℓ = 31, 12 F steps, α = (−159 + √D)/2 | 24/24 | ±λ(s) 24/24 | PASS |

The recorded run searches down to D = −500007 (`--dmax 500000`), in 1 min 19 s, verdict ALL PASS. Two discriminants met
every condition with h ≤ 320 and ℓ ≤ 61. Other attempts:
- **attempt 1** (hmax 160, ℓ ≤ 41, |D| < 10⁷): found only the first, PASS.
  Kept as `binary_validate_run1_hmax160.*`.
- **attempts 2 and 3** (hmax 320, ℓ ≤ 61, |D| < 10⁷): both printed the same two
  PASS rows, then the container's memory-cgroup OOM killer stopped them at 11.5
  and 13.9 GB resident, while they worked on a later candidate. Attempt 3's
  partial output is kept as `binary_validate_attempt3_killed.*`.

`--dmax` was added after attempt 3. Which stage exhausts memory was not
isolated. All commands and the kernel log lines are in `binary_validate_cmd.txt`.

**Binary cost model (derived, not measured).**
- Over F_{2^254}, a 4-dimensional split replaces about 127 shared doublings
  with about 64.
- At the López–Dahab or λ-coordinate doubling cost of about 4M (squarings
  free), that saves about 250 M. Against it stands a once-per-call loop of
  128.5 M for D = −1532327: one 3-isogeny (22.5 M) and 53 F steps charged at
  2 M each.
- The 53 F steps compose to x ↦ x^{2^53}, a linear map that a multi-squaring
  table can apply for about 1–2 M. The model's charge is therefore an upper bound.
- So the binary lead looks *more* favourable than the prime one: doublings
  cost less, but the loop's F steps cost almost nothing.
- It is untimed. The crate has no fast GF(2^254) arithmetic, and the prime-case
  experience (§3) is that the model overstates the gain by 0.05–0.17.
- Security status of such curves was not examined beyond prime order. That
  includes GHS/Weil descent for F_{2^{2m}} with m prime, which Oliveira et al.
  discuss for GLS254, and the moderate class number h = 1397.

## 5. Native sweeper port and x-only re-run

*Pending.* The native sweeper port and x-only re-run are in progress in a separate worktree. This section and `sweeper_native/` are filled in when that run lands.

## 6. Open after round 3

- A constant-time, special-prime implementation of GLV-GLS-4 with a moderate-D
  loop. Example: Mersenne p = 2^127 − 1 needs D with 4p = t² + |D|v², a
  condition on D, not on p. The question is whether the 0.05–0.17 shortfall
  against the model holds once the field is fast.
- A binary GLV-GLS-4 benchmark over F_{2^254} with carry-less multiplication
  (§4 model).
- Merging the x-only loop costs (§5) into the GLV-GLS cost model.

## 7. Evidence index

| file | what | run report |
|:--|:--|:--|
| `gls_validate.*`, `validate_run_start_utc.txt` | §2 | `reports/run-…-gls-validate/` |
| `bench_run1.*`, `bench_D*.*` | §3 | `reports/run-…-bench/` |
| `binary_scan_m{83,127,131}.*` | §4 scans | `reports/run-…-binary-scan/` |
| `binary_validate.*`, `binary_validate_run1_hmax160.*`, `binary_validate_attempt3_killed.*`, `binary_validate_cmd.txt` | §4 validation | `reports/run-…-binary-validate/` |
| `sweeper_native/` | §5 | see §5 |
