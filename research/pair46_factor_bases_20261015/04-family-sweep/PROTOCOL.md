# 04-family-sweep: PROTOCOL (frozen before data)

Pair-46 factor-base campaign, item 04. Written 2026-10-15 before the real run.
It was preceded by two dry runs on scratch output (T=10 targets, 20..30 bits and
30 bits only), which were used only to measure timing and to check the pipeline.
Dry-run outcome disclosed: every primary comparison was null at T=10; at that T the
M1/M2 positive control was not detected, and the M3 control was (p=0.006).
This file is not edited after the run starts. Any change goes in a dated amendment
at the bottom.

## Question

Pair 46 lies in the one-parameter family p = 5t^2+19t+19, which gives exactly
t^2 - 4p = -19 (t+2)^2. The CM field is Q(sqrt(-19)) and the Frobenius conductor is v = t+2.
Does any small-x factor-base property of the family's j = -96^3 (crater) curves
differ from (N1) random D=-19 curves of prime order that are not in the family, or
from (N2) random prime-order curves, at matched bit size?

## Claim classes

- Exact: the identity t^2-4p = -19(t+2)^2. The j and order checks on the pair-46 curves.
  The S3 self-test. The fact that every j=-96^3 curve over F_p (D=-19, units {±1})
  has Aut = {±1}, so the factor base has no extra automorphism orbit compression.
- Empirical comparative: "the family's metric M is within / outside the null
  distribution at matched size". Verdict wording comes from the decision rule below.
- Out of scope: any claim about floor curves (End = Z[pi]) of family members.
  Producing them needs an isogeny-volcano descent, which this item does not run.
  It is cut, and it is listed as cut.

## Arms (one instrument for all)

- **F**: every t >= 1 with p = 5t^2+19t+19 of 20..30 bits, p prime and N = p+1-t prime.
  An enumeration (dry run) gives 100 members. Curve: j=-96^3, with the twist chosen so
  that #E = N (verified by N·P = O on 3 random points), then a uniformly random model
  (a c^4, b c^6).
- **N1** (null, CM-matched): for each F member, 2 curves. A random prime p with the same bit
  length, with 4p = a^2+19y^2, *excluding* the family (y = |±a+2|). Trace ±a is chosen
  among the signs that give prime N (random when both do). The j=-96^3 curve with
  #E = N, random model.
- **N2** (null, generic): for each F member, 2 curves. A random prime p with the same bit
  length, random (a,b), prime order certified by BSGS over the Hasse interval
  (a unique multiple of a random point, which is prime, and confirmed on 3 points).
- **PC** (positive control): for each F member, the first of its N2 curves, with a known
  planted effect. The model is best-of-2 by small-x count (raises M3). The membership bound
  is B·1.03 while the normalization stays at B (expected M1 ratio 1.03^2 = 1.061, and M2 rises
  from ~0.632 to ~0.655).

## Instrument (identical on every arm)

B = ceil(sqrt(2N)). Factor base F_x = { x in [0,B) : x^3+ax+b is a nonzero square } (small-x).
T = 120 targets per curve, each a uniformly random affine point R (a random x until it lies on the curve).
For each target x3, every x1 in F_x with x1 != x3: solve S3(x1, X, x3) = 0 (the Semaev
3rd summation polynomial, which is quadratic in X) and count roots x2 with x1 <= x2, x2 in F_x.
- **M1** (normalized Semaev m=2 decomposition yield) = total count / (T · 2|F_x|^2/N).
- **M2** (relation rate) = fraction of targets with >= 1 decomposition.
- **M3** (small-x yield) = |F_x| / (B/2).
The S3 formula is self-tested on 200 random sums before any data.

## Statistics and decision rule

Statistic: difference of arm means of the per-curve metric, X - Y. Stratified permutation test
within bit length, 9999 permutations, two-sided p = (1+#{|d*|>=|d|})/10000. 95% CI from a
stratified bootstrap (2000 resamples), with the ratio of means and its CI.
Primary family: F vs N1 and F vs N2 for M1, M2, M3 (6 tests), Bonferroni α = 0.05/6 = 0.00833.
Power check: 2/(9999+1) = 2e-4 < α. Positive control: PC vs N2 on M1, M2, M3, same test and α.
If PC is not detected for a metric, that metric's primary verdict is "positive control
failed — inconclusive". Secondary (descriptive, not used for verdicts): N1 vs N2, and every test
per bit band 20-24, 25-27, 28-30.
Verdict wording: "within the control distribution (p=…, diff …[CI])" or "differs from the null
(p=…, diff …[CI])".
Expected sensitivity at T=120: per-curve SE of M1 is about 0.09, the difference SE about 0.011,
and the MDE (80% power) about +-0.04 in M1. For M2 the difference SE is about 0.005, and for M3 it is
about 0.0007 at 30 bits.

## pair-46-scale add-on (M3 only)

The given crater and floor models at p = 5947698030651343, with B = 2^16, are compared with
200 random models of the crater class and with 200 random (a,b) over the same p (rank reported).
Semaev M1/M2 at 53 bits is cut, because about 5e7 root-finds per target are infeasible.

## Seeds, sizes, stop rule

Master seed 46004, with a per-curve RNG derived from it (see src/main.rs). Bits 20..30, T=120.
Pre-flight estimate: the 30-bit-only dry run took about 22 ms per target of user time, including
generation, and the weight is about 54 thirty-bit member-equivalents. That gives about 54·6·120·0.022 ≈ 850 s of
one core (wall time is about 1.4x on the shared machine). Stop rule: the run is fixed-N and is not
stopped early on results. If the wall time exceeds 35 min, the 30-bit stratum is dropped first, then the
29-bit stratum, and the cut is reported.

Command: `cargo build --release && taskset -c 3 target/release/family_sweep . 120 20 30 9999 2000`
