# 13-xy-windows: PROTOCOL (frozen before data)

Pair-46 factor-base campaign, item 13. Written 2026-10-15 before any data was generated.
Not edited after the run starts. The sha256 is recorded in `PROTOCOL.sha256` before the run.

## Question

Can a factor base cut by **joint** windows on x and y, i.e. by y-sign, y-bit or
|y|-window constraints together with an x-window, beat an x-only window of matched
size? The comparison is made either in relation yield or in a way an index-calculus
solver could use.

## Step 0: coverage check (documentation, done before this protocol)

- `research/factor_base_lab.md` (crypto-autoresearcher) lists `bit_window` as "fixed high or low
  radix bits of x": the predicate uses x only. `src/crypto_autoresearcher/index_calculus/factor_base.py`
  implements small_x, random, subgroup, dickson, known_log, planted, small_x_offset. None of these
  constrains y, and there is no bit_window implementation. The ecbench plug-in families in
  `docs/ecbench/FACTOR_BASE_PROTOCOL.md` are also x-only.
- The proposal IDEA-20261015-5888a7 (ledger, "joint x-y windows", notes say "nothing measured")
  states the same gap. It predicts conserved matched-size yield, an x-scan price of p/B2 lifts per
  element, and one open measurement: box-count dispersion below the Weil scale.
- Conclusion: bit_window does **not** cover joint windows. This item implements them.

## Claim classes

**Exact** (proved in RESULTS.md; the computation only checks them):
- X1. For any Y subset of F_p, the factor base {P in E : x(P) in X, y(P) in Y}, taken modulo +-,
  has x-projection {x in X : f(x) is a nonzero square, and +-sqrt(f(x)) meets Y}. Semaev
  relations R = +-P1 +- P2 (and S_{m+1} in general) depend only on that x-set. Every joint
  window is therefore an x-predicate.
- X2. If Y and -Y together cover F_p^* (the half-plane y < p/2, or y even, since p is odd and
  y -> p-y flips parity), the x-projection equals that of the x-only window. Every relation is
  identical, target by target.
- X3. Membership degree: K distinct abscissae have a univariate vanishing polynomial of degree
  exactly K, and a nonzero function on E vanishing at K affine points has pole order >= K.
  A y-constraint cannot lower either bound below the base size. This is the campaign's existing
  RR/bounded-degree obstruction, specialised here.

**Empirical:**
- E1. Toy matched-size relation yield (M2, D) of the joint-window arms against the x-only window.
- E2. Construction price: x-candidates scanned (sqrt lifts) per base element, per arm.
- E3. Pair 46 (crater, floor) against 4 random curves on the same p: box-count dispersion and the
  small-x/small-y corner density.

## Toy arms (one instrument; all arms of size exactly K = ceil(sqrt(N/2)) +- classes, except X_half)

Each arm scans x = 0,1,2,... and keeps the first K liftable x (f(x) nonzero square) that pass its
predicate. The exception is R. cy = min(y, p-y) is the canonical |y|.
- X: x-only window (the matched-size baseline; the comparator).
- XY_sign: keep x if some y with y^2=f(x) has y < p/2 (X2 predicts that this equals X).
- XY_par: keep x if some root y is even (X2 predicts that this equals X).
- XY_box2, XY_box4, XY_box6: keep x if cy < ((p+1)/2) >> k, k = 2,4,6. These keep a fraction 2^-k
  of the liftable x.
- XY_bits4: keep x if bits 1..4 of cy are zero. This is a 1/16 y-residue constraint, matched to box4.
- Xbits4: keep x if x = 0 mod 16. This is the x-only bit_window family at the same 1/16 thinning,
  the x-only comparator for XY_bits4 and XY_box4.
- R (negative control): K uniformly random liftable x in F_p (seeded).
- PC (positive control): the K classes {jG : j = 1..K}, an arithmetic progression in the log.
  The yield concentrates on targets with log in [-2K, 2K], so M2 must collapse (< 0.5 of X)
  at matched size.
- X_half (size control): the first ceil(K/2) elements of X. D must fall to ~1/4.

Toy curves: `input_curves_from_04.csv` (copied verbatim from item 05/04, sha256 88c26329...7c12),
all 100 F rows (family p=5t^2+19t+19, j=-96^3, prime N) and the first 100 N2 rows (random prime
order), 20..30 bits. Per curve: check N·P = O on 3 random points. T = 40 random targets R = kG,
shared by all arms. Seed 46013 xor curve index.

## Metrics

- M2: fraction of targets with >= 1 decomposition R = +-P1 +- P2, x(P1), x(P2) in V, counting
  unordered pairs.
- D: mean number of decompositions per target.
- Effect: ratio of the pooled sums sum_c(arm)/sum_c(X) for M2 and D, with a 95% percentile
  bootstrap over curves (paired by curve, 2000 resamples, seed 46013).
- Scan cost c_scan = number of x candidates scanned to fill V, divided by K.
- Exact checks per target: D(XY_sign) = D(XY_par) = D(X), and the V-sets are identical.

## Pair-46 arms (E3)

Curves: crater, floor (both with N·P = O checked on 3 points), plus 4 random curves
y^2=x^3+ax+b over the same p (a, b from seed 46013; any order). The random curves are the
negative-control population. For each curve, scan x in [0, 2^23) and record cy for liftable x.
- Dispersion: boxes [x0, x0+B1) x [c, c+H) in (x, cy) space, with B1 in {2^8, 2^12, 2^16} and
  H = floor(16 p / B1), so lambda ~ 16. That gives 2^23/B1 x-blocks, and the y-bands tile
  [0, p/2) (the last partial band is dropped). Index I = sample variance / mean of the box counts.
  95% interval: bootstrap over x-blocks (2000 resamples). Null: I = 1.
- Planted positive control: the crater counts at each shape with +8 added to a hashed 10% of
  boxes. The expected index is ~1.34. The detector passes iff the interval excludes 1 and
  lies above 1.2.
- Corner density: among liftable x in [0, 2^23), the fraction with cy < ((p+1)/2) >> k for
  k = 0..16, as obs/expected (expected 2^-k) with a 95% Wilson interval. This is the
  realized density of the XY_box_k base at pair-46 size, and it fixes the x-scan price
  (1/(density · 1/2) scans per element).

## Decision rule (fixed in advance)

- Coverage: answered in Step 0 (not covered).
- **Null (yield conserved)** iff, for XY_box2/4/6, XY_bits4 and Xbits4, the D-ratio and M2-ratio
  intervals each lie inside [0.85, 1.15], AND R's intervals lie inside [0.85, 1.15], AND PC's
  M2 ratio is < 0.5 with its interval excluding 1, AND X_half's D ratio interval contains 0.25.
  (If a positive control fails, the instrument is declared blind and no null is claimed.)
- **Advance** would require some joint arm with a D or M2 interval entirely above 1.15. Even then,
  the 2^k scan price must be paid.
- Exact checks X2 must hold on every curve and target. A failure is an instrument bug.
- E3 null iff every real-curve dispersion interval contains 1 or lies inside [0.97, 1.03], and the
  planted control passes. If pair-46 values fall outside the range of the 4 random curves,
  that is reported as a difference.

## Stop rule

Fixed sizes, one run, no adaptive extension. Budget guard: toy curves are processed in file
order, and processing stops when wall time exceeds 25 minutes. Any cut is recorded in
results.json and RESULTS.md. Single thread, `taskset -c 3`.

## Implementation

Rust, `src/main.rs`, single thread. The arithmetic helpers are copied in `src/common.rs` from item 05.
Python only plots recorded values (`plot.py`) and lays out the PDF (`make_pdf.py`, copied from item 05).
