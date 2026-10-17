# 09-cm-norm — preregistered protocol (frozen before any data)

Item 09 of the pair-46 factor-base campaign: norm-form / CM-class-ideal factor bases.

## Step 0: what is a factor-base ELEMENT? (settled before data)

Two candidate readings, both fixed here:

- **(I) Endomorphism / principal-ideal points.** A factor-base element is a
  point alpha(P) with alpha in End(E) of small norm, i.e. alpha = a + b*w with
  Nm(alpha) = a^2 + a*b + 5*b^2 <= X (crater, O_K = Z[w], h = 1, so every ideal
  is principal and "ideal of small norm" = "endomorphism of small degree" up to
  units +-1). On E(F_p) ~ Z/N, alpha acts as the scalar a + b*lambda mod N.
  On the floor, End = Z[pi] and pi acts as the identity on E(F_p), so
  a + b*pi acts as the integer a + b.
- **(II) Non-principal class-group ideals (floor, h(Z[pi]) ~ v).** An ideal
  class [a] acts on the set of curves with End = Z[pi]; its element is an
  isogeny E -> E_a, not a point of E(F_p). This reading is assessed as an
  exact (theoretical + exactly computed) statement, not measured.

## Claim classes

- C1 (exact, computed): pair-46 consistency: N = p + 1 - t; lambda^2 - lambda + 5 = 0 mod N;
  v*lambda = 2 mod N (from pi = v*w - 1 and pi = 1 on rational points);
  h(Z[pi]) = v - (-19/v) (class-number formula, units of Z[pi] are +-1);
  minimal norm of a non-integer element of Z[pi] is >= 19 v^2 / 4.
- C2 (exact, proof in RESULTS): reading (I) factor base F_X is the image of a
  lattice-point set in an ellipse, so |F_X + F_X| <= #{Nm <= 4X} ~ 4|F_X|
  (small doubling); on the floor, reading (I) degenerates to an arithmetic
  progression for X < 19 v^2/4.
- C3 (empirical, toy): relation coverage of reading-(I) factor bases versus a
  random set of the same size, on toy prime-order curves with CM by -19.

## Metric (same instrument on every arm)

For a symmetric factor base S subset of Z/N (logs of the points, S = -S, 0 not in S)
of size 2B, coverage_m(S) = |{s_1 + ... + s_k : k <= m, s_i in S}| / N,
computed exactly with a bitset over Z/N. The relation probability of a uniform
random target R = [r]P + [c]Q under "R is a sum of <= m factor-base points" is
exactly coverage_m. Effect size: ratio coverage_m(arm) / mean coverage_m(random arm)
on the same curve and B, with a 95% interval from the 2.5/97.5 percentiles of
the random arm (null distribution) and a percentile bootstrap over curves.

## Arms

- A_ideal_crater: B smallest-norm alpha up to sign (ordered by norm, then a, then b), S = {+-(a+b*lambda)}; run for both roots lambda of x^2-x+5 mod N.
- A_floor: reading (I) on the floor: S = {+-1..+-B} (integers, by C2 exact degeneration). This is identical to the positive control and is reported as such.
- NULL (negative control): uniformly random symmetric S of size 2B, 20 seeds (seed = 1000*curve_index + rep).
- SMALLX (prior baseline): the B points with smallest x-coordinate (x = 0,1,2,...), logs found by enumerating all multiples of the generator.
- POS (positive control): arithmetic progression S = {+-1,...,+-B}. The detector must show coverage_2 ratio <= 0.1 for B >= 256; otherwise the instrument is declared blind and C3 is not reported.

## Sizes, curves, seeds

- Toy curves: j = -96^3 over F_p, choose the twist of prime order N; p is the
  smallest primes above 2^20 (6 curves) and above 2^22 (2 curves) with 4p = t^2 + 19u^2
  and a prime-order twist. Generator: the point with smallest x of order N.
- B in {32, 64, 128, 256, 512, 1024}; m in {1, 2}; m = 3 for B <= 128.
- 20 NULL seeds per (curve, B).

## Stop rule

Fixed design, run once. Hard cap 40 min wall on 1 core; if exceeded, drop
the 2^22 curves first, then m = 3, and say so. No reruns with changed parameters.

## Decision rule

- "advances" only if A_ideal coverage_2 ratio's lower interval endpoint > 1 at B >= 256.
- "negative" if its upper endpoint < 1 (worse than random) for B >= 256.
- Otherwise null.
