# Extension fields `GF(p^k)` (B5a's design)

**Written 2026-10-01, before any B5a code.** The plan is
`research/notes/index-calculus/IC_TOOL_PROGRAM.md` (§9, B5). B5a's
protocol, [`../rounds/B5a-extension-fields/PROTOCOL.md`](../rounds/B5a-extension-fields/PROTOCOL.md),
declares its measurement. It came to `main` with B5a's native
re-declaration (2026-10-06), which changed §4's generator and nothing
else.

## 1. Where the tool stands

The plan's B5 has two halves, and they share no code:
- **prime fields past one word**, which `rho-bignum` already runs, slowly;
- **extension fields** `GF(p^k)`, which nothing runs.

This design is the second, **B5a**. The first, a multi-limb `rho-negation`
for `p ≥ 2^63`, is **B5b**, declared on its own later.

With Track B's stack (B0–B4, on record in [`../track-b/`](../track-b/README.md)):
- **Parsing and validation exist.** Schema v2 reads a `prime_extension`
  field (§3.1), and B2 validates it (§4.2–§4.5):
  - Rabin's test that the modulus is irreducible;
  - the curve, in short or general Weierstrass form;
  - `#E` by §4.5's method 3 (`r > 4√q`) or method 4 (`q ≤ 2^32`).
- **Routing goes nowhere.** Every valid extension instance is refused as
  `no-pipeline-for-field` (C050, whose `until` is B5).
- **No name.** ICV1 had no extension-field part. Schema v2 §6 asks for one
  before the tool writes an extension-field result, so this declaration
  adds it (§2.1).
- **The code that exists:**
  - B2's `curves.rs`: `Fpk`, big-integer arithmetic for any `p`, any `k`
    and the document's own modulus, and `Short<F>`, the curve over it.
    Validation uses them.
  - `rho_bignum.rs`: rho with the negation map over any `BigGroup`.
  - `ic_boundary.rs`: `rho_reference_negation`, the matched rho over any
    `CountedGroup`, with 64-bit keys and scalars.
  - `gaudry_cubic.rs`: Gaudry's index calculus on `E(GF(p³))`, with the
    factor base `{P : x(P) ∈ GF(p)}`, from the residual-walk thread
    (`research/notes/index-calculus/RESEARCH_RESIDUAL_WALKS.md` §11). It
    has three oracles, dense or Wiedemann linear algebra, and large
    primes. Its field is `GF(p)[t]/(t³ − c)`, and its group has prime
    order.

## 2. What B5a adds

### 2.1 ICV1 for `GF(p^k)`

`docs/curves/ICV1.md` gains an extension kind, which changes no
existing identity.
- **The field part** is `fpk-<p>-<k>-<modhash8>`, hashed over the
  modulus's low coefficients `c_0, …, c_{k−1}`. The modulus is part of
  the model, as it is for binary fields.
- **The trace** is `p^k + 1 − #E`.
- **`j`** is written as its `k` coefficients.
- **The slug's field tag** is `fp<bits of p>k<k>`.

`scripts/curve_id.py` implements it now, with self-test vectors: brute
force over `GF(5²)` and `GF(7³)`, checked against the trace recurrence
from `GF(p)`. `scripts/check_curve_names.py` reads the new tag. The Rust
port, `curve_id::extension`, comes with B5a's code: it is pinned to the
reference's vectors in `tests/curve_id.rs`, and it gives a report its
`curve_id`.

### 2.2 `rho-negation` on `GF(p^k)`

**The same walk on a third field kind.**
- It is `rho_reference_negation` with:
  - the look-ahead step;
  - the short-cycle escape by doubling;
  - distinguished points;
  - stride starts.
- `ExtCurve` implements its `CountedGroup` over the document's own field
  and basis:
  - **Coefficients:** one-word, so `p < 2^31` whenever `k ≥ 2` and
    `q ≤ 2^62`. Products are summed in `u128` and reduced modulo `p`.
  - **The modulus:** the general monic one, as given.
  - **Inversion:** the extended Euclidean algorithm in `GF(p)[t]`.
  - **Storage:** elements sit in a fixed-capacity array, eight
    coefficients for `k ≤ 8` and 32 above. Since `5^k ≤ 2^62`,
    `k ≤ 26`.
- **The point key** is `2(x̂ + 1) + s`.
  - `x̂ = Σ x_i p^i < q` packs the abscissa.
  - `s` says whether the packed ordinate `ŷ` exceeds the packing of
    `−y`.
  - The key is injective for `q ≤ 2^62`, as binary's `2(x + 1) + s` is
    for `n ≤ 62`.
  - `{P, −P}` is represented by its point with the smaller key, as on
    the other kinds (`NegationClasses`).
- **Gates:**
  - `field-wider-than-one-word` when `q > 2^62`, since the key needs
    `⌈log₂ q⌉ + 2` bits;
  - `scalar-wider-than-63-bits` when `r ≥ 2^63`.
- **Counting:** one affine addition a step, as on prime fields.

### 2.3 `rho-bignum` on `GF(p^k)`

A `BigGroup` over B2's `Fpk` makes every valid extension instance run
`solve: rho`, within its budget.
- **Any `p` and any `k`.**
- **The jump and distinguished-point hash** is `hash_limbs` over the limbs
  of the packed abscissa `x̂`.
- **The class representative** is the sign of the packed ordinate.
- **The table's key** is the bytes of `(x̂, ŷ)`.

### 2.4 `ic-gaudry-cubic`: Gaudry's index calculus on `E(GF(p³))`

**Imported, as B2 imported `ic boundary`'s pipelines.**

**What it admits.** The first three gates are new:

| condition | why | gate when it fails |
|:--|:--|:--|
| `k = 3` | the module is the cubic case | `extension-degree-not-three` |
| the modulus is `t³ − c`: `c₂ = c₁ = 0`, so the document's `c₀` is `−c` | the module's field is `GF(p)[t]/(t³ − c)`, and its Weil restriction is written in that basis | `modulus-not-binomial` |
| `h = 1` | the module's relations are taken modulo `#E`, with every point in `⟨G⟩` | `cofactor-not-one` |
| `r < 2^63` | its linear algebra is modulo `n` in `u64`, so `p < 2^21` | `scalar-wider-than-63-bits` |
| a recipe object has only this pipeline's keys (below) | the module has its own knobs | `recipe-not-taken` |

`p ≡ 1 (mod 3)` follows from validation: `t³ − c` is irreducible only
for a non-cube `c`, and non-cubes exist only then.

**The field is the document's.** The basis `1, t, t²` is the module's,
with `c = −c₀ mod p`, so elements need no map.
- `Fp3::with_c(p, c)` is a constructor added to the module. Its `new`
  draws `c` at random, and it is otherwise unchanged.
- `Curve3` takes the document's `a`, `b` and `G`, with `n = r`.
- `Instance3` takes the target as `q`. Its `d` is the known logarithm
  when the document gives one; the module reads `d` only for its own
  report. Every arm's answer is verified by the tool, as everywhere:
  `[d]G = Q`, then replayed in `Short<Fpk>`.

**The recipe** is the module's options:

| key | values | the module's |
|:--|:--|:--|
| `oracle` | `mitm`, `groebner`, `pair-only` | `Solver` |
| `sparse_la` | true or false | Wiedemann, or dense elimination |
| `small_base`, `max_large_primes`, `max_merge_level` | integers | the large-prime variation |
| `max_residuals` | an integer | the residual cap |

- `recipe: auto` is the module's `GaudryOptions::default()`: the
  meet-in-the-middle oracle, dense elimination, no large primes and
  `10⁶` residuals. This follows design §5.3's rule for imported
  pipelines, and the report names the default.
- A recipe object is this pipeline's only when every key is in the
  table. Any other key, `kic`'s `summands` say, is refused at the gate
  as `recipe-not-taken`, as B2's amendment 1 refuses a recipe on the
  other imported pipelines. A key whose value is out of range is
  `value-invalid`.
- The relation seed is fixed, as `ic-prime-s3`'s is.

**Set-up and online.**
- **Set-up:** the subspace factor base (`SubspaceBase::build`, about
  `p/2` points).
- **Online:** relation collection on random `R = [a]G + [b]Q`, the
  linear algebra and the logarithm. With the `groebner` oracle, the
  symmetrised `S₄`'s once-per-curve precomputation is inside this
  interval too: the module runs it inside its run, and the import does
  not change the module.
- As for B2's imported pipelines, the online work is not split into
  AGENTS.md's five exclusive phases, so the report marks the pair
  ineligible for a speedup claim.

**Units.** The module's:
- affine additions in `E(GF(p³))`;
- `GF(p)` multiplications, converted at the ratio it measures on the
  instance (`fp_muls_per_add`);
- `S = total / √r`.

Rho's steps are affine additions in the same group, so both arms' `S`
are in one unit.

**The matched reference is `rho-negation` on `GF(p³)`** (§2.2), with
`A = 2`. The module's own `run_rho3` has no negation map, so it is not
the reference.

**Disclosures:**
- `study-pipeline`: the residual-walk thread's implementation. At the
  sizes it was measured, it costs `528×` to `4,000×` rho (§11.7 there).
  Asymptotically it is `Õ(p^{4/3})` against rho's `p^{3/2}`, with `p`
  the base field's size (Gaudry 2009, `Õ(p^{2−2/k})` at `k = 3`).
- `factor-base`: the points whose abscissa is in the prime field
  `GF(p)`. That is the base field, not a proper intermediate subfield
  (AGENTS.md §8b), so `method-uses-subfield` is not set.

### 2.5 Routing

**For an extension instance:**
- `kic`, `ic-binary-s4`, `ic-prime-s3` and `rho-koblitz` keep
  `no-pipeline-for-field`.
- `ic-gaudry-cubic` and `rho-negation` have their gates (§2.2, §2.4).
- `rho-bignum` has no gate.

**For a binary or prime instance:** `ic-gaudry-cubic` is listed with
`no-pipeline-for-field`, and nothing else changes.

**Arm choice.** Under `paired`:
- the index calculus arm is `ic-gaudry-cubic` if it admits the instance;
- the rho arm is `rho-negation`, or `rho-bignum` when `rho-negation`
  does not admit it.

An instance with no admitted index calculus is `no-ic-route`, as on the
other kinds.

### 2.6 Estimates

Each estimate is a model with its source, as in B2's. The constants are
in [`../rounds/B5a-extension-fields/estimates.json`](../rounds/B5a-extension-fields/estimates.json).

| pipeline | estimate | constants |
|:--|:--|:--|
| `rho-negation` | `√(πr/4)` steps at a step cost on `GF(p^k)` | provisional, then set by B5a's measurement 6 |
| `rho-bignum` | the same, at its step cost by the width of `q` | provisional, then set by measurement 6 |
| `ic-gaudry-cubic` | its measured cost at its largest rung, carried along its fitted exponent, as `imported_ic` does for `ic-prime-s3` | the anchor is §11.2's: `p = 2083`, `n = 2^33.1`, `361.7·10⁶` operations, `ops ∝ n^{0.69}`. Measurement 6 re-anchors it on the reference host |

### 2.7 Targets and F1

**Targets.** `known_log` and `random_seed` are B2's rules, as on every
kind. `public_hash_seed` is the prime field's rule read on the packing:
- `x` is the element whose packing `Σ x_i p^i` is `H mod q`, where `H` is
  the first `bits(q) + 64` bits of B2's counter-mode output under its
  label;
- `y` is the root with the smaller packing;
- the target is that point times `h`, the first that is not the identity.

**F1.** F1 has no model on extension fields: B7a samples `kic` and the
rho walks on binary and prime fields. `fidelity: F1` on an extension
instance is refused as `no-f1-model`; F0 runs there.

### 2.8 Disclosures for every extension instance

`intermediate-subfields` (B2) states the proper divisors `d` of `k`.
B5a adds `curve-subfield` for extension fields:
- it states the least `d` dividing `k` with `a, b ∈ GF(p^d)`, tested as
  `a^{p^d} = a`;
- when `d < k`, it also states that the `p^d`-power Frobenius is an
  endomorphism `rho-negation` leaves unused, as the binary disclosure
  does.

## 3. Sameness

B5a adds field kinds and changes no binary or prime path.
- **The pin:** every suite row's output is unchanged (protocol,
  measurement 3).
- **The module:** on `gaudry_cubic`'s own instances
  (`generate_instance3`), written as v2 documents, the tool's
  `ic-gaudry-cubic` gives the module's own run, counter for counter and
  logarithm for logarithm. The tool adds the routing, the reading of the
  document and the verification, and nothing in between.

## 4. Instances

Found by `icprog b5a instances` (`src/bin/icprog/b5a.rs`), in arithmetic
of its own that shares nothing with the tool. Its output is
[`instances.json`](../rounds/B5a-extension-fields/instances.json). #1178's
Python generator found the same records, which the native one writes
again byte for byte.
- **Group orders:** baby steps and giant steps over the Hasse interval,
  narrowed by several points until one order is left. This is exact, and
  it is B2's generator's method.
- **Subgroups:** factored by trial division and Brent's method, each
  prime factor certified by B1's exact test.
- **Labels:** every choice comes from SHAKE-256 of a public label, with
  `:1`, `:2`, … appended until it works.

| id | field | modulus | `h` | what it exercises |
|:--|:--|:--|:--|:--|
| `G1`, `G2`, `G3` | `GF(271³)`, `GF(523³)`, `GF(1039³)` | `t³ − c`, `c` the least non-cube | 1 | `ic-gaudry-cubic` paired with `rho-negation`, at §11.2's first three primes, `n ≈ 2^24`, `2^27`, `2^30` |
| `H1` | `GF(271³)` | `t³ − c` | `> 1` | `cofactor-not-one` |
| `E2` | `GF(p²)`, `p = 2^31 − 1` | a labelled irreducible quadratic | `> 1` | the one-word arithmetic at its widest, `q = (2^31 − 1)² < 2^62` |
| `E5` | `GF(2053⁵)` | a labelled irreducible quintic | `> 1` | a general degree, `q ≈ 2^55` |
| `E11` | `GF(37^11)` | a labelled irreducible of degree 11 | `> 1` | past eight coefficients, the 32-coefficient capacity, `q ≈ 2^57` |
| `B2` | `GF(p²)`, `p` the least prime above `2^35` | a labelled irreducible quadratic | `> 1` | `q ≈ 2^70`, past one word: `rho-bignum` |

**Subgroups for rho:**
- each is the largest prime factor `r` of `#E` with `4√q < r < 2^44`,
  so §4.5's method 3 establishes `#E`, and F0 takes seconds;
- `r²` does not divide `#E`.

**Targets.** Each instance has two, as B4's do:
- a known logarithm's multiple;
- a public point `T001`, whose logarithm the tool is not given.

**C050's document,** `GF(1009³)` with a general modulus and
`h = 1524`, is run as it is (its successor case, §5).

## 5. Tests the code must carry

- `ExtCurve`'s field against B2's `Fpk` on random elements:
  - `k` from 2 to 8, 11 and 26;
  - `p` from 5 to `2^31 − 1`;
  - random irreducible moduli, binomial and not.
- `ExtCurve`'s group law against `Short<Fpk>`, and keys injective on
  every point of a small curve, enumerated.
- `rho-negation` and `rho-bignum` recover known logarithms on small
  instances at `k = 2, 3, 5`.
- The `BigGroup` over `Fpk`: canonical classes, and agreement with
  `Short<Fpk>`.
- §3's sameness with the module, on `generate_instance3`'s instances at
  `p = 271` and `523`.
- `curve_id::extension` equals the reference's vectors: the self-test's,
  and each instance's.

## 6. What it does not do

- **Prime fields past one word.** That is B5b.
- **A basis change.** `ic-gaudry-cubic` refuses these, though rho runs
  them:
  - a cubic modulus that is not `t³ − c`;
  - `p ≡ 2 (mod 3)`, where no `t³ − c` is irreducible.

  A field isomorphism would lift the first: a root of the document's `f`
  in `GF(p)[s]/(s³ − c)`, found by equal-degree splitting, maps the
  document's elements to the module's. The second needs the module's
  arithmetic over a general cubic.
- **Cofactors in the index calculus.** The relations would be taken on
  `[h]P_i` in `⟨G⟩`.
- **Other degrees.** `gaudry_quartic.rs` (`k = 4`, §11.16–§11.19 of the
  residual-walk note) is the next import. Gaudry's method at `k = 2` does
  not beat rho.
- **Characteristic 3.** It is still refused at validation.
- **The `p^d`-power Frobenius** of a curve over a subfield, which rho
  could fold. It is disclosed (§2.8), not used.
- **No speed claim.** These are new rows, measured as new rows.
