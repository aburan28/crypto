# Binary fields of three words and more (B4's design)

**Written 2026-10-01, before any B4 code.** The plan is
`research/notes/index-calculus/IC_TOOL_PROGRAM.md` (§9, B4). B4's
protocol, [`../rounds/B4-multi-word/PROTOCOL.md`](../rounds/B4-multi-word/PROTOCOL.md),
declares its measurement.

## 1. Where the tool stands

With Track B's stack (B0–B3b, on record in [`../track-b/`](../track-b/README.md)):
- **one word** (`n ≤ 62`): every pipeline, on `FastCurve`;
- **two words** (`63 ≤ n ≤ 126`): `rho-koblitz` (B3) and `kic` (B3b), on
  `GfWide`, a `u128`;
- **beyond 126**: refused as `field-wider-than-two-words` (C055, C086).

So the curves this workstream is about are still refused:
- ECC2K-130, `n = 131` (AGENTS.md §8b);
- the standard Koblitz curves, `sect163k1` to `sect571k1`.

`ic fixed` reaches `n = 131`, in Python, outside the schema.

## 2. What B4 adds

**One more pipeline, generic in the word count `W`**, from 3 to 9. A
class key packs to `2(x + 1) + s`, so `W` words take `n ≤ 64W − 2`, and
nine words take `sect571k1`. The plan's B4 is three words (plan §9's
`n ≤ 191` is the field's width; a class key needs `n ≤ 190`); the
pipeline is generic, so nine cost no more code, and the
standard curves from `sect233k1` to `sect571k1` pass the field gate.
Their subgroups then meet the scalar gate (§2's table), and none is
affordable at F0 anyway.

| layer | B4's kernel |
|:--|:--|
| field | `GfMulti<W>`: elements in `[u64; W]`. Products from 64-bit carry-less limb products (PCLMULQDQ or PMULL where present, portable otherwise), Karatsuba at `W = 3`. Reduction by the modulus's tail when it is a trinomial or pentanomial, as every registered modulus is, and Barrett otherwise. Squaring, Itoh–Tsujii inversion, the trace, the half-trace, and the quadratic solve |
| curve | `MultiCurve<W>`: the curves `kic` takes (Koblitz, and curves over `F_{2^k}`, `k ≤ 8`), affine, with batched inversion. `MultiCanon<W>`: the normal basis and the least rotation that name a signed-Frobenius class |
| rho | `MultiRho<W>`: `WideRho`'s walk, step for step |
| kic | select, build, scan and descent at `W` words. The pair table is keyed by a 64-bit hash of the class key, in the one-word table's layout |
| scalars and logs | as at two words: `r < 2^127`, in `u128` scalars and B3b's two-limb residues. A larger `r` is refused as `scalar-wider-than-127-bits`, which B7b lifts for F1 at `n = 131`, where `r ≈ 2^129` |
| import | v2 documents for `K_0` at `n = 131` through the same importer, in place of `ic fixed`'s files: the challenge's points are polynomial-basis (§4) |
| targets | v1's hashed target (`public_hash_seed`) past one word, by the rule below |

**What earlier steps left to B4.** B3 and B3b refuse three of v1's rules
past one word as `not-yet-supported`, naming B4 (schema v2 §5.2's
`not-yet-supported`; `src/bin/ic/v2.rs` on B3b's branch):
- v1's hashed target (`public_hash_seed`);
- v1's random target (`random_seed`);
- v1's generator rule (`koblitz_search_v1`), past one word or under a
  modulus other than the repository's.

B4 defines each past one word so that it is v1's rule at one word, below.

**Three earlier cases wait for B4**, by their `until`:
- **C055**: `rho-koblitz` refuses `n = 131` as `field-wider-than-two-words`
  (B3). From B4 it admits it.
- **C086**: `kic` refuses `n = 131` the same way (B3b). From B4 it admits
  it.
- **C056**: a v1 hashed target on a field past one word is refused as
  `not-yet-supported` (B3). From B4 it is derived.

**The hashed target past one word.** v1's rule (`public_hash_target` in
`src/bin/ic/workflow.rs`) hashes a fixed domain, the curve and the seed
with a counter, using BLAKE3:
- `x` is the digest's first 8 bytes as a little-endian word, masked to
  `n` bits;
- byte 8 chooses between the two points with that `x`;
- the point is multiplied by the cofactor, and the counter moves on if
  `x` has no point or the multiple is the identity.

B4 reads `8W + 1` bytes of BLAKE3's extended output over the same input:
- `x` is the first `W` little-endian words, masked to `n` bits;
- byte `8W` makes the choice.

At one word this is v1's rule byte for byte, since the extended output
begins with the digest. A test checks that at every one-word degree, and
the rule serves the two-word pipelines too (C056's successor runs at
`n = 83`).

**The random target past one word.** v1 seeds `StdRng` with
`random_seed ⊕ 0x534f4c5645525447` and draws `k` uniformly from
`[1, r)`, reading `r` as its low word (`known_scalar` in
`src/bin/ic/workflow.rs`). That reading is right only for `r < 2^64`.
- For `r < 2^64`, B4 draws exactly as v1 does.
- For `r ≥ 2^64`, it draws `k` uniformly from `[1, r)` with
  `num-bigint`'s `gen_biguint_range`, from the same seeded generator.

**v1's generator rule past one word, and under any modulus.** v1 takes
the order's largest prime `r`, then steps a 64-bit LCG from a seed of
`n` and `a`, reads each state, masked to `n` bits, as an abscissa, and
keeps the first lift, in the order `points_with_x` returns them, whose
cofactor multiple is not the identity (`KoblitzCurve::subfield`).
- **The subgroup is the document's.** Past one word the order's
  largest prime cannot be found in general, so v2 documents state `r`
  and `h`. The checks verify them as before.
- **An abscissa takes `W` successive LCG states**, as its little-endian
  words, masked to `n` bits. At one word that is one state, as in v1.
- **The modulus is the document's.** The abscissa is read in the
  document's polynomial basis. Under the repository's modulus at one
  word that is v1's curve.

Each rule has a test at every one-word degree against v1's own code.

## 3. Three pipelines, one algorithm

**The one- and two-word pipelines are left as they are.** Every round's
F0 timing depends on the first, and B3b's measurement on the second.

**The multi-word pipeline makes the same choices everywhere**, as B3b's
did: the draws, the lift, the orbit and row order, the fold, the
canonical key's normal element and tie rule, the table's placement, and
the walks' steps.

**Sameness is checked three ways.** The pipeline is generic, so it also
compiles at `W = 1` and `W = 2`, in tests only:
- at one-word degrees, `W = 1` must equal the one-word pipeline;
- at two-word degrees, `W = 2` must equal B3b's;
- at `W = 3`, on the F0 instances below, the answers must verify.

Each comparison is counter for counter: the same base and table, the same
relations in the same order, the same logarithms and counters.

**Two rules where no narrower pipeline has one**, settled while writing
the first code and before this declaration merged:
- **A base's abscissae past `n = 127`.** The narrower pipelines draw an
  abscissa from `[1, 2^n)` with `rng.gen_range` in a `u64` up to
  `n = 63` and in a `u128` up to `n = 127`; past that the multi-word
  pipeline draws it with `gen_biguint_range(1, 2^n)` from the same
  generator, the random target's rule (§2).
- **The pair table's hash of a key past two words.** The two-word hash
  mixes the high word into the low one through `fmix64` before the
  one-word hash. The multi-word hash folds every word above the lowest
  in from the top the same way, so on every key below `2^128` it is the
  two-word hash, and below `2^64` the one-word hash.

Both make the pipeline at `W` words choose what the narrower pipelines
choose wherever they run, so the sameness tests also run it at `W = 3`
and `W = 9` on one- and two-word degrees.

## 4. Instances

**F0 at three words: six Koblitz curves of prime degree.** No proper
intermediate subfield is involved, as with B3b's instances. Each has a
prime `r` dividing its order once, small enough to solve end to end,
and the signed Frobenius group of order `2n` acts freely on the
subgroup: `π` acts there as `λ mod r`, of order exactly `n`.

| `n` | `a` | `log₂ r` | `r` | cofactor bits | slug |
|--:|--:|--:|--:|--:|:--|
| 127 | 0 | 30.06 | 1118452171 | 97 | `icv1-f2m127-t24589614856193402413-a00e5890` |
| 137 | 0 | 37.48 | 191818802977 | 100 | `icv1-f2m137-t574016314927011818501-11ceb4d6` |
| 151 | 0 | 34.71 | 28016550571 | 117 | `icv1-f2m151-tm97937651335354183120307-ee31df7e` |
| 157 | 1 | 47.14 | 154781804431543 | 110 | `icv1-f2m157-t158060339213695215877259-816eafb2` |
| 173 | 1 | 27.63 | 208110697 | 146 | `icv1-f2m173-tm67870783603944754053042229-ddf2de63` |
| 179 | 0 | 39.58 | 820651535909 | 140 | `icv1-f2m179-t1681527843948629186391379613-d1b73c24` |

- **How they were found** ([`../rounds/B4-multi-word/instances.py`](../rounds/B4-multi-word/instances.py),
  output in `instances.json`): every `K_a` with `127 ≤ n ≤ 190`, its
  order's prime factors in `[2^25, 2^52)` found by ECM within a time
  budget, kept when the action is free.
- **`r` is not the order's largest prime.** At prime degree the order's
  large part is almost never `2^52`-smooth, so a smaller prime factor
  names the subgroup. The cofactor, 97 to 146 bits, only clears points
  into it.
- **Why the action matters.** The folded table counts `2n` points in
  every signed-Frobenius class. The class is that large exactly when
  `λ` has order `n` and, at even `n`, `λ^{n/2} ≠ −1`.
  - At even `n` the order contains the half-degree twist's, and on a
    subgroup of the twist `π^{n/2}` acts as `−1`, so its classes have `n`
    points.
  - `kic` and `rho-koblitz` already refuse even `n` (`even-extension-degree`,
    schema v2 §5.2), so that case cannot reach them.
  - At odd composite `n`, `r` can divide a proper subfield curve's order.
    Its points are then defined over that subfield, and their classes
    are smaller.
  - The search keeps only instances where the action is free. B4's six
    are of prime degree, where it is free for every `r > 4`.

**F1 at ECC2K-130** (`icv1-f2m131-tm22283658519494248867-115e0dc5`),
`r ≈ 2^129`: B7b, with these kernels.

**The challenge's points are polynomial-basis, and B4 imports no normal
basis.**
- Certicom's generator and target satisfy `y² + xy = x³ + 1` in the
  polynomial basis of `x^131 + x^13 + x^2 + x + 1` (checked for this
  design).
- The repository already holds them in that basis:
  `docs/ic/params/ecc2k130-fixed.json` (`"basis": "polynomial"`), and
  B1's conformance document `ecc2k130-challenge.json`.
- Three places say otherwise:
  - schema v2 §10–11, which defers a "normal-basis import" to B4;
  - `ic params ecc2k-130`'s profile ("point coordinates not imported");
  - the curve registry builder's reason for leaving ECC2K-130 without an
    EC1 identity.

  B4's code corrects the last two. This design corrects the first, and
  the schema's deferral lapses.
- `ic fixed`'s `onb_root` converts the other way, from the polynomial
  basis to the permuted type-II normal basis its kernels use. Nothing
  here needs it.

**The standard Koblitz curves** (`sect163k1` to `sect571k1`) and
ECC2K-130 pass the field gate from B4 and meet the scalar gate, since
their subgroups are 129 to 570 bits (C096, C099). None is affordable at
F0. A three-word subgroup below `2^127` that is too large for F0 is
admitted and refused as over budget (C102: the 95-bit prime of
`#K_0(GF(2^127))`, whose order is `4 · 1118452171 · r`).

## 5. Tests the code must carry

- The field against the general arithmetic (`binary_ecc`): every degree
  from 127 to 574, with the repository's modulus for each, and every
  registered curve's modulus.
- The curve and `MultiCanon` against `binary_ecc`'s points and a
  Frobenius walk.
- Sameness at `W = 1` and `W = 2`, as §3 defines it.
- v1's three rules past one word (§2), each equal to v1's own code at
  every one-word degree.
- Each three-word F0 instance at a reduced size, within the test budget.
- ECC2K-130's published generator and target, read from B1's v2
  document at three words: on the curve, in the subgroup, and the
  generator's multiple by `r` the identity.

## 6. What it does not do

- **No speed claim.** The multi-word pipeline is a new path, measured as
  new rows.
- **Not F1 at `n = 131`.** Wiring F1's samples to these kernels is B7b.
- **No subfield methods.** An instance whose field has a proper
  intermediate subfield says so (AGENTS.md §8b), and nothing here uses
  the subfield.
