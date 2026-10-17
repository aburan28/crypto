# Residual zero classes, model V2: exact 8th-power mass, random thinning (registered 2026-10-16)

V1 ([PROTOCOL.md](PROTOCOL.md)) failed on its training primes 37-47
(overdispersion about 5, zero counts about twice the prediction). Its
diagnosis used those primes, so they are not pristine for V2. V2 was
derived only from exhaustive small-prime censuses (`iso1_lambda_census`,
every lambda of the relevant subgroup with a certified order):

* p = 11 full census: membership in T is a random thinning of
  Omega = (F_Q^*)^(2^e) inside each class (dispersion 1.012, 31 zero
  classes against 32.4 expected).
* p = 7, 11, 13, 17 (e = 4, 3, 3, 5): the number N_8(t) of 8th-power
  Legendre parameters with |trace| = t is exactly a local 2-adic factor
  times the odd-conductor Hurwitz number. 32nd powers at p = 17 are not.

## Registered formula

For an ordinary positive trace t with k = v2(f_pi) >= 2, s = (D_K / 2)
and H_odd(t) = (h_K / (w_K/2)) * sum_{d | f_odd} d prod_{l | d} (1 - (D_K/l)/l):

    L_2 = (2 - s) * 2^(k-2)                      if t != +-(Q+1) (mod 32)
    L_2 = 5 (2 - s) * 2^(k-2) - 6 (1 - s)        if t == +-(Q+1) (mod 32)
    N_8(t) = H_odd(t) * L_2

(N_8 = 0 when k = 1, by the theorem.) The weak count is modelled as a
thinning of the N_8(t)/12 absolute-Frobenius orbits with rate
r = |T| / |(F_Q^*)^8| = 8 / (p^2 - 1): weak orbits ~ Binomial(N_8/12, r),
P(zero) = (1 - r)^(N_8 / 12). Nothing is fitted.

## Tests

A. Exactness of N_8: exhaustive 8th-power censuses at p = 19 and 23
   (not used above). Pass: formula equals the census count on every row.
B. Holdouts p = 53, 59 (pristine) and every prime from the certified
   queue (61 onward); p = 37-47 reported separately as non-pristine.
   Per prime: observed zero rows against the Poisson-binomial 95% interval
   (pass if inside); dispersion of weak orbit counts against the binomial
   variance for rows with mean >= 5 (reported; near 1 expected); ROC AUC.
C. A pass is a statistical description beyond the proved obstruction, not
   a criterion that decides an individual class.

## Follow-up registration (2026-10-16, after p = 97 was scored; before p = 101 completes)

Observation. Over the ten scored primes 53-97 the observed number of zero
rows exceeds the V2 prediction by +95.3 (sd 49.8, z = 1.91), about 1.9% of
the expected total, spread over the bulk of the predicted-probability range
rather than the tails. Splitting by e = v2(p^2 - 1) does not isolate it
(`iso1_residual_model pool`: e = 3 z = 1.62, e = 4 z = -0.39, e = 5 z = 1.19,
e = 6 z = 1.82, e >= 4 z = 1.12), which rules out the idea that it comes
from the 8th-power model being inexact at e > 3. This hypothesis was
suggested by data already scored, so those primes cannot test it.

Registered test, fixed before any of these censuses completes. Hypothesis H:
V2 under-predicts the number of zero rows. Fresh primes: 101, 103, 107, 109,
113 (their class data and predictions are committed before each census
finishes). Statistic: the pooled `obs - exp` over those five primes divided
by `sqrt(sum(var))`, with each prime's variance `((hi - lo) / 3.92)^2` as in
`iso1_residual_model pool`. Decision: if the one-sided z exceeds 2.0, V2 is
rejected as a parameter-free description (it is miscalibrated by about 2%)
and a V3 is registered with its own holdouts; otherwise V2 stands, and H is
reported as not confirmed. If fewer than five of these primes complete, the
test is run on those that do and reported as underpowered.
