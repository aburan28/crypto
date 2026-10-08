# E5 result: Φ_ℓ mod p construction, measured 2026-10-08

Status: **implemented and measured; one arm of E5 as registered, plus one
arm not registered (the transform), reported separately**

Result class: **engineering speedup of one walker stage; no ECDLP
speedup; `S` unset.**

## What changed

`src/cryptanalysis/isogeny_walk/modpoly.rs`, same `ModPoly` output and the
same vanishing checks on every exponent of `Φ_ℓ(j(q), j(q^ℓ))`:

- the `ℓ + 2` powers of `j` by a length-`n` transform over `F_{p²}` when
  `n | p² − 1` (`n` the power of two above `2ℓ(ℓ+1)`), two transforms per
  power with the series of `j` kept transformed; Kronecker substitution
  into one big-integer product when the field lacks the roots of unity;
- `σ₃, σ₅` by a divisor sieve and the series inverse by Newton iteration;
- the triangular solve by blocks `b = ℓ+1, …, 0`: `O(1)` lazy terms per
  unknown inside a block, then the block row and its mirror applied to
  the residual as two `O(ℓ³)` products; `O(ℓ⁴)` in all.

Correctness: the old schoolbook construction is kept as a test-only
reference and the fast one reproduces its coefficient tables and check
counts for `ℓ ∈ {3, 5, 7, 11, 13}` over `F_1009` and the P-256 field; the
literature tables for `Φ_2, Φ_3` still pass; `checks_passed` is unchanged
at every `ℓ` measured.

## Timings, P-256 field, one core, Apple Silicon

| ℓ | before (ms) | after (ms) | j-series | powers | solve | checks |
|--:|--:|--:|--:|--:|--:|--:|
| 11 | 8 | 10 | | | | 64 |
| 31 | 1,338 | 316 | 36 | 233 | 47 | 474 |
| 47 | 8,465 | | | | | 1,090 |
| 61 | 33,540 | 1,858 | 205 | 1,279 | 374 | 1,839 |
| 97 | not run | 12,701 | 1,177 | 9,268 | 2,256 | 4,665 |
| 127 | not run | 19,215 | 1,108 | 11,042 | 7,064 | 8,010 |
| 151 | not run | 43,006 | 2,066 | 28,098 | 12,636 | 11,334 |
| 199 | not run | 181,399 | 4,079 | 75,849 | 101,469 | 19,710 |
| 251 | not run | 242,507 | 10,598 | 145,868 | 86,039 | 31,384 |

Phase columns are cumulative differences from the trace and are
approximate.  The 47 row was timed before the change only.

Reading against E5's registered windows: the old construction's fitted
exponent from 31 to 61 is 4.75, inside the registered `≥ 4.5`.  The new
construction fits 3.6 from 61 to 251 (powers 3.3, solve 3.9); the solve is
`O(ℓ⁴)` by design and the powers are `O(ℓ · ℓ² log ℓ)`.  E5.1's Kronecker
arm (`≤ 4.2`) and BLS arm (`≤ 3.5`) were not run as registered; the
transform arm was not registered and is reported here without a verdict
on E5.1.  E5.3 (`Φ_1009` under one hour) is not reached: extrapolating the
solve at `ℓ⁴` gives about 6 hours at `ℓ = 1009`, so BLS or a transform-based
solve is still needed for that target.

## Walk above the old cap

`isogeny_walk walk --curve p256 --primes 97,101,103 --max-curves 8
--threads 4 --no-class-audits`, with the binary before the transform-domain
refinement: `Φ_ℓ` built in 18.9, 19.7 and 20.1 s; 8 curves, 2 expanded,
**8 verified edges at degrees 97, 101 and 103, 0 failures**, 3.2 s of walk.
These are the first certified P-256 edges above degree 61 in this
repository.  Output kept in the session scratchpad only; a registered run
with replay belongs to the protocol, not to this note.

## Memory

The stored powers cost `(ℓ + 2) · ℓ² · 32` bytes: about 0.5 GB at
`ℓ = 251`, 4 GB at `ℓ = 509`.  A streamed variant that recomputes rows
would trade this for time.
