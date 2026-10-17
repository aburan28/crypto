# Frozen protocol: 61-bit prime-conductor-gap set (2026-10-17)

Written before the CM construction ran. The user asked for a separate set of
61-bit curves with a **large prime conductor gap**: the Frobenius conductor
is a single large prime `f = l`, so the class has exactly two strata, the
crater (`End = O_K`) and the floor (`End = Z[pi]`), with crater-to-floor
index `l`. The user approved a CM construction, which gives up seeded
coefficients and `p = 2^61 - 1`.

## Prime field: class-number-one CM

- `d` ranges over every class-number-one discriminant:
  `-3, -4, -7, -8, -11, -19, -43, -67, -163`, with `j(d)` known exactly.
- For each `d`, the construction walks `l` over primes descending from
  `L = isqrt(2^63 / |d|)`, and then `t = 1, 2, ...` with the parity that
  makes `p = (t^2 + |d| l^2) / 4` an integer. A candidate is accepted when:
  - `2^60 < p < 2^61` and `p` is prime;
  - `n = p + 1 - t` is prime (prime order, cofactor 1; `t > 0` so `n < p`);
  - a model is found:
    - `d = -3`: `y^2 = x^3 + b` (secp256k1 form), with the least `b >= 1`
      whose curve has order `n`;
    - `d = -4`: `y^2 = x^3 + a x` (`j = 1728`), with the least `a >= 1`
      whose curve has order `n`;
    - otherwise: `p = 3 (mod 4)` is required, and the `j(d)` curve is put
      in the twist with order `n`, then into `a = -3` form (P-curve shape)
      when `-3/A` is a fourth power; if it is not, the candidate is
      rejected and the search moves on.
- Orders are proved by BSGS uniqueness plus `n` prime with `n > 4 sqrt(p)`.
  `Delta = t^2 - 4p = d l^2` is checked to factor as `d_K = d`, `f = l`.
- The first accepted candidate per `d` is recorded. The featured
  prime-field gap curve is the `a = -3` candidate with the largest `l`; the
  `d = -3` candidate is recorded as the secp256k1-form gap curve. Every
  other `d` is recorded too.

## Binary field

Over `GF(2^m)` with `m` prime, `(pi) = P^m` for a prime `P` above 2, so
either `P` is principal (then `8 = a^2 + |d| b^2`, so `d = -7`: the Koblitz
classes) or `m | h(d_K)`. Class-number-one CM therefore gives no new binary
class at `m = 61`. The binary gap curve is the Koblitz curve
`icv1-f2m61-t158598901-ab42b6c5`, whose conductor `1951 * 587551` has two
large prime gaps; its `O_1951` stratum has its own prime gap 587551 to the
floor. A seeded binary search is not expected to find a large `l`: the
largest conductor prime in 32,768 seeded B-61 curves was 79
(`data/tally-B61.json`). CM with `61 | h(d)` needs a Hilbert class
polynomial of degree at least 61, which this study does not build.

## Not claimed

No ECDLP cost, weakness or security statement.

## Note added during the run (before any `d != -3` result was seen)

For `d = -4` and `d = -8`, `t^2 + |d| l^2 = 4p` forces `t` even, so
`n = p + 1 - t` is even and prime order is impossible. The first run spent
10 minutes on `d = -4` without noticing this. These two discriminants are
recorded as impossible, and the `l` loop is capped at `(L/2, L]`.
The second run stalled on `d = -7` for the same kind of reason: for
`d = 1 (mod 8)` odd `t` and odd `l` give `t^2 + |d| l^2 = 0 (mod 8)`, so
`p` would be even; a prime conductor `l` is impossible there. `d = -7` is
recorded as impossible. The rule is otherwise unchanged.
The third run stalled on `d = -11`: there `p = t^2 + 2 l^2 (mod 3)`, so `p`
prime forces `t = 0 (mod 3)` and then `3 | n`; prime order is impossible.
Each `d` is now capped at the first 2000 primes `l` below `L`, and an
exhausted cap is recorded as "none found within the cap".
