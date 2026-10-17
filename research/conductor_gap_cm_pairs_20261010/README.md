# Conductor-gap CM pairs: crater vs floor curves with a certified large gap

**Status: constructed, Sage-verified, measured (five ecbench sessions,
2026-10-10). Result: no implemented method sees the gap; see §5.**

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

- **Generic invariance**, in `lean/GenericInvariance.lean` (no Mathlib).
  An algorithm whose control flow depends only on group-operation outcomes
  (equality answers) and its random tape has the same cost and output
  distribution on isomorphic groups. That covers BSGS-type algorithms. The
  crater and floor curves are both cyclic of order `N`, so BSGS cannot see
  the conductor. **It does not cover the rho and kangaroo that ecbench
  runs.** They choose each step by hashing the point's representation
  (`mix(g.key(&pos))` in the kangaroo, `src/cryptanalysis/ecbench/generic.rs:275`;
  `mix(self.g.key(x))` in the rho walk, `TunedWalker::index` in
  `src/cryptanalysis/ic_boundary.rs`). That is outside the theorem's model.
  Sheriff case SHF-20261016-3caa6f (finding L5) gives a toy demonstration,
  not run on ecbench: on `Z_1019` and its isomorphic image in `F_2039^*`, a
  representation-hashing r-adding walk had a different cost on 197 of 200
  transported instances, while BSGS had a different transcript on 0 of 200.
  For rho and kangaroo, `S = 1` is the expected null, not a theorem.
  (Corrected 2026-10-16. Earlier text said the theorem covers rho and the
  kangaroo.)
- **The transfer bound** (the formal shape of Jao–Miller–Venkatesan), in
  `lean/TransferBound.lean`: any solver on one curve in the isogeny class,
  generic or not, solves the ECDLP on any other at two evaluations of the
  isogeny, returning the same logarithm; with maps both ways the best costs
  on the two levels differ by at most that. For these pairs the isogeny has
  degree ℓ, which is what the construction makes large.

Not formalizable today: anything that mentions endomorphism rings,
conductors, or volcanoes. Mathlib has no isogenies and no `End(E)`. So
Lean proves that BSGS-type methods are blind. The experiment asks whether
the implemented methods outside the theorem are blind too: the
representation-hashing rho and kangaroo, `ic.pipeline`, and any structural
candidate.

## 5. ecbench: the measurement (2026-10-10)

Five sessions, all on this Mac at L0 (operation counts only; wall time is
never read), every record verified, every audit `OK` with 12 exact replays
(`sessions/*.audit.json`). `crater_vs_floor.py` pairs each arm's mean `S`
on the crater's workloads against its mean on the floor workloads, with a
bootstrap over workloads; `results/crater_vs_floor.md` has every row.

| ℓ | q bits | session | arms | targets per curve | floor/crater (95%) |
|---|---|---|---|---|---|
| 263 | 19 | `ell263_q19` | rho, bsgs, kangaroo, ic.pipeline, control | 4 | rho 0.97 [0.78, 1.28], bsgs 0.84 [0.61, 1.32], kangaroo 0.80 [0.58, 1.20], ic 0.974 [0.953, 0.992] |
| 8209 | 29 | `ell8209_q29` | rho, bsgs, kangaroo, control | 4 | rho 1.22 [1.005, 1.46], bsgs 0.98 [0.65, 1.50], kangaroo 0.95 [0.68, 1.49] |
| 131101 | 37 | `ell131101_q37` | rho, bsgs, kangaroo, control | 4 | rho 1.28 [0.81, 2.16], bsgs 1.25 [0.82, 1.79], kangaroo 1.06 [0.68, 1.46] |
| 263 | 19 | `ell263_q19_ic16` | rho, ic.pipeline, control | 16 | ic 0.998 [0.970, 1.026], rho 1.05 [0.96, 1.16] |
| 8209 | 29 | `ell8209_q29_ic16` | rho, ic.pipeline, control | 16 | ic 0.96 [0.87, 1.07], rho 1.04 [0.91, 1.19] |

Reading, by the decision rule of `PROTOCOL.md`:

- Two intervals excluded 1 in the first three sessions: ic at ℓ = 263
  and rho at ℓ = 8209, which was 1.22 [1.005, 1.46]. Each rested on four
  crater workloads, and they had opposite signs across ℓ. **The rho
  interval is unexplained.** Correction (2026-10-16, SHF-20261016-3caa6f
  L5): this bullet used to call it "noise by theorem". It is not. The
  generic-invariance theorem does not cover ecbench's rho, which hashes the
  point representation (§4). The `ell8209_q29_ic16` session put rho at
  1.04 [0.91, 1.19], with a fresh target seed (20261011) and 16 targets
  per curve, so the deviation did not reproduce there. All five sessions,
  though, used the same walk seed (`measurement.seed = 1`). Still open: a
  rho-only re-test at ℓ = 8209 under a fresh measurement seed. The A/A
  arm `rho-neg-aa` matches `rho-neg` row for row in
  `results/crater_vs_floor.md`, so it does not act as an independent
  control for this question. The ic lead was followed up first.
- The follow-up with 16 targets per curve and a fresh target seed puts
  `ic.pipeline` at 0.998 [0.970, 1.026] at ℓ = 263 and 0.96 [0.87, 1.07] at
  ℓ = 8209, with the control arm at 1 in both. **No implemented method sees
  the conductor gap**, at gaps of 2^8, 2^13 and 2^17, on curves that share
  q, N, t and the Frobenius discriminant and differ only in `End(E)`.
- This is the honest negative the falsifiers section predicted: it extends
  the small-gap negative of the volcano probe to the largest gaps the
  construction gives at these sizes, at no isogeny cost.

`specs/gap-ell263-q19-constructor.json` is the same pair expressed through
the in-harness `conductor_gap` constructor (`src/cryptanalysis/ecbench/
conductor_gap.rs`, `CurveSpec::ConductorGap`), which rebuilds the pair from
`(ℓ, t, level, seed)` so a record replays the search; it passes `ecbench
plan` and is not run here because the explicit-form sessions above already
measure the same curves.

### Relation to the ECC2K-130 conductor-gap study on main

`research/isogeny_conductor_gap_20261007/` and
`docs/ic/EXPERIMENTS_ISOGENY_CONDUCTOR_GAP_20261007.md` study the same
question inside the binary Koblitz class of ECC2K-130 (`D = −7`, gaps at
primes such as 263 and 967 reached by Vélu descent from the crater). Its E1
panel ran the production index-calculus pipeline on the crater and on floor
curves at the toy rung `n = 23` and found the per-call statistics
level-blind (trials per relation 2.41–2.52 on the crater, 2.50–2.78 on the
floor). Its only level-dependent effect is orbit folding: the crater's
Frobenius endomorphism `τ` is defined over `F_2`, so its factor base folds
by about `n`, and floor curves lose that.

This directory is the prime-field complement. Here `D = −19` gives a crater
whose extra endomorphisms have norm at least 5 and are not used by any
method in ecbench, so there is no folding to lose. Every method then
measures as level-blind. Together the two say the same thing from both
characteristics: in the implemented methods, the conductor matters only
through an endomorphism the method actually exploits, never through the
ring's index as such. That is the Lean pair of §4 in practice. Methods
whose control flow depends only on group-operation outcomes (BSGS-type) are
blind by theorem. The measured rho and kangaroo hash the representation, so
for them blindness is measured, not proved (the ℓ = 8209 rho interval in §5
is still unexplained). A structural method is no better than the
endomorphisms it uses.

## 6. Falsifiers

- A BSGS-type arm with `log₂ S` excluding 0: harness fault (check the
  control arm first). A rho or kangaroo arm with `log₂ S` excluding 0:
  unexplained until re-tested under a fresh measurement seed. These
  methods hash the representation, so generic invariance does not rule out
  a real effect (§4, correction of 2026-10-16).
- `ic.pipeline` with `log₂ S` excluding 0 at both ℓ = 263 and ℓ = 8209
  with the same sign: the first implemented method that sees the
  conductor. Build a third pair at a new ℓ before claiming anything.
- Everything at `S = 1`: the honest negative, extending the small-gap
  negative of the volcano probe to gaps of 2^8 to 2^17 at no extra cost,
  and a reason to stop looking for conductor effects in the current
  methods.
