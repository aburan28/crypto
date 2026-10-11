# Conductor-gap CM pairs: crater vs floor curves with a certified large gap

**Status: construction and verification, 2026-10-10. No ECDLP measurement
has been run on these pairs yet; `PROTOCOL.md` is frozen for that.**

## 1. Why these pairs

The isogeny-volcano probe in this repository compared curves at adjacent
levels (conductor gaps of 2 and 3) and found no structural method that
cares. The open question was whether a *large* gap would change that, and
the obvious route, walking ℓ-isogenies up a volcano for ℓ in the hundreds
or thousands, is blocked by the modular polynomial Φ_ℓ.

The gap does not need the walk. With a class-number-one CM field it can be
certified from the Frobenius alone:

- Pick `D = −19` (`h(D) = 1`, `j_D = −884736 = −96³`), a prime `ℓ`, and odd
  `t` with `q = (t² + 19 ℓ²)/4` prime and `N = q + 1 − t` prime.
- Then `t² − 4q = −19 ℓ²`, so every ordinary curve over `GF(q)` with trace
  `t` has endomorphism ring of conductor `f_E ∈ {1, ℓ}` (the conductor
  divides `f_π = ℓ` and `ℓ` is prime).
- The curve with `j = j_D` is the crater (`End = O_K`, `f_E = 1`); the
  twist with order `N` is the one used. Every other curve of order `N` is
  a floor curve (`f_E = ℓ`), because the class group of `O_K` is trivial
  so the crater has exactly one `j`.
- Floor curves are found by sampling `(a, b)`, one scalar multiplication
  `[N]P` per sample, and keeping curves where it vanishes (choosing the
  twist when `[N]P` vanishes on it instead). `N` is prime and larger than
  the Hasse interval, so one point of order `N` already pins `#E = N`.

Zero isogenies are computed. The two levels are joined only by degree-ℓ
isogenies, so for `ℓ ≥ 263` nothing in the repository can even exhibit the
isogeny; the pairs sit exactly in the regime the question is about. At a
fixed field size `ℓ` can be traded against `t` up to about `√q/2`, so a
ladder of gaps at one `q` is cheap.

### Discriminants that do not work

Found while prototyping, each verified by exhaustive search over `t`:

- `D = −7`, and any `D ≡ 1 mod 8`: `t` and `f` are forced even, so the
  conductor is never an odd prime.
- `D = −11`: `3 | q·N` for every odd `t` (`q ≡ 0 mod 3` unless `3 | t`,
  `N ≡ 0 mod 3` unless `3 | t − 2`), so `q` and `N` cannot both be prime.
- `D ∈ {−19, −43, −67, −163}` work; `−19` is used here.

## 2. What was built

| file | ℓ | q bits | curves | floor samples per hit | wall per floor hit |
|---|---|---|---|---|---|
| `results/ell_263_q19.txt` | 263 | 19 | 1 crater + 4 floor | ≈1.3k (predicted ≈1.2k) | 0.5 s |
| `results/ell_8209_q29.txt` | 8209 | 29 | 1 crater + 3 floor | ≈96k (predicted ≈39k, 3 hits) | ≈50 s |
| `results/ell_131101_q37.txt` | 131101 | 37 | 1 crater + 1 floor | ≈970k (predicted ≈620k, 1 hit) | ≈1000 s |

Prediction: with `h(O_K) = 1` there are `ℓ − (−19/ℓ) ≈ ℓ` isomorphism
classes of conductor `ℓ` at trace `t`, and as many at trace `−t`. A random
`(a, b)` lands in a given class with probability `1/(2q)` (each class has
`(q − 1)/2` Weierstrass models among `q²` pairs), so accepting both twists
gives about `q/ℓ ≈ 19ℓ/4` samples per floor hit. Observed rates agree
within the Poisson noise of one to three hits.

`gap_pairs.rs` is the single-file prototype (no dependencies; build with
`rustc -O gap_pairs.rs`; run as `gap_pairs <k> <floor_curves>` for
`ℓ ≈ 2^k`). It is affine, single-threaded, and was run on a machine
saturated by other builds. Projective arithmetic and `rayon` make
`ℓ ≈ 2^25` at 53 bits a matter of minutes per floor curve.

## 3. Verification

`results/sage_verification.txt` is the output of `verify_sage.py` over all
three files. For every curve, independently of the prototype: `#E = N`
exactly (SEA), the j-invariant matches, the crater has `j = −884736` and no
floor curve does, the Frobenius order has discriminant `−19 ℓ²`, and the
embedding degree exceeds 20. `q`, `N`, `ℓ` prime and `t² − 4q = −19 ℓ²` hold
in every file.

## 4. What Lean can and cannot prove about this

Provable now, with Mathlib and no curve-specific library:

- **Generic invariance.** A generic-group algorithm's cost and output
  distribution is identical on isomorphic groups. The crater and floor
  curves are both cyclic of order `N`, so rho, BSGS and the kangaroo
  cannot see the conductor. This is the null hypothesis of
  `PROTOCOL.md`, as a theorem rather than a guess.
- **The transfer bound** (the formal shape of Jao–Miller–Venkatesan):
  an ECDLP oracle on one curve in an isogeny class solves the ECDLP on
  any other at the cost of the isogeny walk.

Not formalizable today: anything that mentions endomorphism rings,
conductors, or volcanoes. Mathlib has no isogenies and no `End(E)`. So
Lean proves the generic methods are blind, and the experiment asks whether
the implemented non-generic ones (`ic.pipeline` and any structural
candidate) are too.

## 5. ecbench

`specs/*.json` are `ecbench.spec/v1` specs, one per pair, crater first,
as `prime_explicit` constructions (`make_spec.py` derives them from
`results/`). They are not yet run. The intended `CurveSpec` variant that
builds a pair from `(ℓ, t, seed)` inside ecbench, so the search is
replayed from the record instead of pasted, lives in
`src/cryptanalysis`, which is held by T-237 at the time of writing; it is
the follow-up, not part of this directory.

## 6. Falsifiers

- A generic arm with `log₂ S` excluding 0: harness fault (check the
  control arm first).
- `ic.pipeline` with `log₂ S` excluding 0 at both ℓ = 263 and ℓ = 8209
  with the same sign: the first implemented method that sees the
  conductor. Build a third pair at a new ℓ before claiming anything.
- Everything at `S = 1`: the honest negative, extending the small-gap
  negative of the volcano probe to gaps of 2^8 to 2^17 at no extra cost,
  and a reason to stop looking for conductor effects in the current
  methods.
