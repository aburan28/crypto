# Protocol: cubic GHS to genus 3, priced end to end against rho (item 7)

Frozen 2026-10-14, before any data. Campaign rules: `../RULES.md`.
Design: `reports/Index calculus leads next steps.md`, "### 7.".

## Check first (done before freezing)

- `research/ghs_cubic_quartic_ic_20261012/`: plane-quartic IC, 90/90 correct at
  q in {5, 7, 11}; rho arm records step counts only, so no IC-vs-rho ratio exists.
- `research/ghs_cubic_hyperelliptic_genus3_20261013/`: the hyperelliptic
  `m = 3, r = 5` cover is built and the transfer is checked (4800/4800
  homomorphism checks) at q <= 13 with Zech tables of `F_{q^6}`; it runs no
  IC and no DLP.
- No directory, benchmark row or commit prices this transfer against rho.
  Item 7 has not been run.

## Hypothesis

H1 (from the report): on the hyperelliptic branch, genus-3 index calculus
with double large primes after the cubic GHS transfer costs about `q^{4/3}`,
below rho on `E0(F_{q^3})` at about `q^{3/2}`, with a crossover at some `q`.

Kill criterion (from the report): a fitted IC exponent >= 1.35 in `q` (not
separated from rho's 1.5), or an `ell`-part that line relations with large
primes cannot reach at `q >= 2^9`, closes the route.

## Branches

- **Hyperelliptic branch (`m = 3, r = 5`, split arm): run.** `E0 : y^2 =
  Q(x)(x - r0)(x - r1)` over `K = F_{q^3}`, `Q = (x-a)(x-b)` in `F_q[x]`,
  `r_i = p^{q^i}`. `C : w^2 = h(z)`, `deg h = 7`, over `F_q`, built from the
  conic parametrisation of `D = C/iota` with an `F_q`-rational Weierstrass
  point moved to infinity. The transfer is
  `phi(P) = sum_{j=0..2} Fr^j(pi0^*[P] - 4 inf)`, computed entirely over `K`
  with Cantor arithmetic (`pi0^*[P]` is the Mumford pair
  `(xn - x_P xd, y_P R(z) mod U)`).
- **Plane-quartic branch (`m = 2`, Diem `O~(q)`): BLOCKED up front.** The
  repository has no class arithmetic on the plane quartic
  (`problems/plane-quartic-jacobian-class-arithmetic-missing.md`); building
  it is the weeks-long engineering the report names, outside a 3-hour item.

## Arms (one instrument, all in Python 3 with explicit operation counters)

1. **IC-2LP (test arm).** Transfer `G`, `T = kG` to `Jac(C)(F_q)`; r-adding
   walk (r = 20) on `alpha D1 + beta D2`; smoothness test `x^q = x mod u`;
   factor base = the `N_FB = ceil(kappa q^{2/3})` degree-1 primes (up to sign)
   with smallest `x`; other degree-1 primes are large primes; relations with
   <= 2 large primes kept; large primes eliminated on the spanning tree of the
   component of the large-prime graph that contains the factor-base vertex;
   linear algebra by Lanczos on the normal equations mod `ell` with the
   unknown `k` as an extra column. The walk restarts from a fresh random point
   when it revisits a divisor (cost counted).
2. **IC-basic (positive control for the exponent detector).** Same code,
   factor base = all degree-1 primes, no large primes. Theory: relation phase
   `~ q`, linear algebra `~ q^2`; the fit must separate it from IC-2LP.
3. **Rho (reference).** Single-walk r-adding rho (r = 20) with distinguished
   points on the short Weierstrass model of `E0` over `K`, affine
   coordinates, no negation map. Sensitivity: the ideal negation-map
   rho (`sqrt(pi ell / 4)` steps at the measured per-step price).
4. **Null control for the instrument.** The rho fit itself: its exponent in
   `ell` must be 0.5 within its interval, and in `q` about 1.5.

`kappa` is chosen once, by a pilot at the smallest `q` only, from
{1, 2, 4} (lowest mean total cost); the pilot is recorded and then frozen.

## Sizes

Primes `q` nearest above 2^8, 2^9, 2^10, 2^11, 2^12, 2^13 (2^14 if the
budget allows), three curves per size (`Q`, `p` drawn from seed 20261014,
accepted when `#E0(K) = c ell` with `ell` prime, `c <= 64`). One planted
`k` per curve for IC-2LP; rho runs 3 targets per curve (1 at the largest
size if wall time forces it). IC-basic at the sizes it finishes in 45 min.
Order of `Jac(C)(F_q)` from PARI `hyperellcharpoly`; checked by
`#Jac * D = 0` and `ell phi(G) = 0`, `phi(G) != 0`, and 3 random
homomorphism checks `phi(P+Q) = phi(P) + phi(Q)` per curve.

## Metric

Cost in **`F_q`-multiplication equivalents**:
`cost = M_q + c_I I_q + c_ell M_ell`, with `M_q` counted `F_q`
multiplications and squarings (including inside `K` arithmetic), `I_q`
counted `F_q` inversions, `M_ell` counted multiplications mod `ell` in the
linear algebra. Primary: `c_I = ceil(log2 q)` (Fermat inversion),
`c_ell = 1` (one machine word at these sizes). Sensitivity: `c_ell = 9`.
IC cost includes the transfer, the walk, smoothness tests, root finding,
large-prime processing (no field ops, not counted) and linear algebra.
Group operations are recorded as a second unit.

Fit: least squares of `log cost` on `log q` per arm; 95% interval by
bootstrap over curves (2000 resamples). Crossover: the `q` where the fitted
lines meet, reported as an extrapolation when outside the measured range,
and the measured ratio IC/rho at every size.

Correctness: every solve is compared with the planted `k`; a wrong or
missing answer is recorded as a failure, never dropped.

## Decision rule

- POSITIVE: IC-2LP solves all instances, fitted exponent < 1.35 with the
  upper end of its interval < 1.5, and IC/rho decreasing in `q`.
- NEGATIVE: fitted IC-2LP exponent >= 1.35 (kill criterion), or IC-2LP
  cannot reach the `ell`-part at `q >= 2^9`.
- INCONCLUSIVE otherwise (for example exponent below 1.35 with an interval
  that reaches 1.5).
The plane-quartic branch is reported BLOCKED regardless.

## Amendment A1 (2026-10-14, after the kappa pilot and one timing run at 2^12, before any recorded data)

- Pilot (`results/pilot_kappa_k{1,2,4}.jsonl`, q = 257, curves 0-2):
  mean IC-2LP cost 541,757 / 501,231 / 576,169 for kappa = 1 / 2 / 4.
  **kappa = 2 frozen.**
- Point counting moved from PARI `hyperellcharpoly` (33 s per try at 2^14)
  to SEA `ellcard` on `E0` over `F_{q^3}`; equality with
  `hyperellcharpoly` is recorded as a check at q <= 2^11.
- Wall times are far below the plan (2^12: IC 9 s, rho 2 s), so the sizes are
  extended to q = next prime above 2^b for b = 8..16, and b = 17 if the
  budget allows. Rho: 3 targets per curve for b <= 14, 2 for b >= 15.
  IC-basic for b <= 13. IC and rho run as separate processes on the same
  instance and the same G (drawn from the same seed before any target), so
  their planted k differ.
- The timing run at b = 12 (curve 0) is not part of the record; the
  recorded run regenerates it from its seed.

## Deviations recorded after data (2026-10-14, not part of the frozen plan)

- D1. The recorded run stopped at b = 15. At b = 15 one IC-2LP solve took
  622 s and one rho target 512 s of wall time on the shared machine; a b = 16
  rho job (2 targets) would reach the 45-minute per-job cap, and the item
  budget is 3 hours. The runner was stopped before any b = 16 job started; the
  b = 15 curve-2 jobs were launched by hand with the same arguments.
- D2. The rho arm builds its 20 multipliers from two random combinations by a
  Fibonacci chain (to keep set-up cost at 4 scalar multiplications). Its
  measured constant came out at about 1.76 sqrt(ell), above the random-mapping
  1.25 sqrt(ell). A control with 20 independent random multipliers
  (`RHO_RANDOM_MULT=1`, b = 9..12, 6 targets per curve) was added after seeing
  this. The verdict therefore uses, besides the measured arm, a **rho model**:
  measured cost per rho step times sqrt(pi ell / 2) (and sqrt(pi ell / 4) for
  an ideal negation map), which is the strongest rho the instrument supports.
