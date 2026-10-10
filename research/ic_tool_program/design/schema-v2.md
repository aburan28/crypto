# Track B design: parameter schema v2, validation and routing

**Declared 2026-10-01, before any B1 or B2 code.** The plan's B1 and B2
steps (`research/notes/index-calculus/IC_TOOL_PROGRAM.md` §9) are built
and judged against this design. It fixes:
- the schema;
- the checks and their codes;
- the routing rules;
- the conformance cases.

Each step still declares its own `PROTOCOL.md`; B1's is
[`../rounds/B1-schema-v2/PROTOCOL.md`](../rounds/B1-schema-v2/PROTOCOL.md).
Once a step has run its first measurement, this design changes only by
a dated amendment at the end.

## 1. Where the tool stands

Five entry points read curve parameters, in four formats. None of them
solves, at speed, a curve the user names:

| entry point | input | fields | curves | what it does |
|:--|:--|:--|:--|:--|
| `ic inspect` | inspection schema v1, or a built-in name (`src/bin/ic/params.rs`) | binary of degree 2..=571 with an explicit polynomial; prime up to 512 bits | any `a`, `b` | checks the parameters; never solves |
| `ic workflow`, `price`, `run` | workflow schema v1 (`src/bin/ic/workflow.rs`): `curve: {degree, curve_a, subfield, curve_b}` | binary, `n ≤ 63`, with fast arithmetic to 62; the tool picks the modulus | `K_a`, and curves over `GF(2^k)`, `k ≤ 8`, with `n/k` odd and at least 3 | builds its own instance and solves it: the modulus by `find_irreducible_sparse`, the generator by a seeded search, and the largest prime subgroup |
| `ic rho`, `boundary`, `bench` | sizes on the command line | prime `p < 2^63`; random binary of degree 5..=32; Koblitz | generated | generated instances only |
| `ic fixed` | a fixed-parameter file (`docs/ic/FIXED_PARAMETERS.md`) | binary to degree 131, in Python | `K_0` only | solves the named instance, slowly |
| `icx` | a catalog name | 52 curves | named | runs a smaller same-family substitute, never the named curve |

So a curve the user names by its own modulus, generator and target goes
either to a checker or to a slow path. That includes the m = 83 gate
curve. G1 and G3 both need two things: one format that names the
instance exactly, and a way to hand that instance to the fast pipelines.

Three facts make that cheap:
- **The fast field takes any modulus.** `semaev_decomp::Gf2::new` builds
  its reduction table from any irreducible polynomial of degree at most
  63. Only the AVX-512 batched addition depends on the modulus's shape.
  It falls back to the portable path when the shape does not suit.
- **The curve can be filled in from outside.** `KoblitzCurve`'s fields
  are public. An importer can fill them from validated input, compute
  `λ` with `frobenius_eigenvalue_q`, and hand the curve to the pipeline
  unchanged.
- **The checks already run at every size.** `ic inspect` validates with
  `BigUint` arithmetic (`binary_ecc`, `ecc`). Validation needs no new
  arithmetic, and it need not wait for B3–B5.

## 2. Principles

1. **The named instance, or a refusal.** The tool runs exactly the
   input's:
   - field and modulus;
   - curve;
   - subgroup and generator;
   - target.

   It never substitutes another instance. `icx`'s scaled runs stay in
   `icx`.
2. **Validity and support are separate questions.**
   - Validation decides whether the input is a well-defined
     discrete-log instance.
   - Routing decides whether the tool can run it, and at which level
     (plan §3).

   A valid instance the tool cannot run is `unsupported`, with the
   reason and the step that would lift it. It is never reported as
   invalid.
3. **One way to say each thing.** Parsing is strict:
   - an unknown key is refused at every level;
   - integers have one syntax;
   - every default the tool fills in is written back into the report.
4. **Every refusal has a stable code.** Messages may change. A code
   keeps its meaning and is never reused.
5. **v1 keeps its outputs.** A v1 file runs exactly as before. Its v2
   translation (§8) gives identical outputs, and B1 checks that on all
   90 suite rows.
6. **Exact before probable.** Each check states whether it is exact or a
   screen, and the report carries that. A screen is never reported as a
   proof.

## 3. Schema v2

A document is JSON, at most 1 MiB, with these keys:

| key | required | what |
|:--|:--|:--|
| `schema_version` | yes | `2` |
| `name` | yes | a label of 1–120 non-control characters; not an identity |
| `field`, `curve`, `subgroup` | yes, unless `named` is given | §3.1–§3.3 |
| `named` | instead of `field`, `curve` and `subgroup` | a standard name, §3.7 |
| `target` | yes | exactly one, §3.4 |
| `method` | no | §3.5 |
| `budget` | no | §3.6 |

**Integers.** The schema has two kinds:
- **Mathematical integers** are strings: field elements, `p`, `r`,
  `h`, scalars and coordinates.
  - A string is decimal digits, or `0x` followed by hexadecimal digits.
  - It takes no sign and no spaces, and has at most 1,024 bits.
  - The report writes them back in one form: decimal for scalars and
    orders, lower-case hexadecimal for binary-field elements.
- **Configuration integers** are JSON numbers: degrees, seeds, counts
  and budgets.

### 3.1 The field

A binary field, in a polynomial basis:

```json
"field": {"kind": "binary", "degree": 83, "modulus": "0x800000000200000000007"}
```

- `degree` is 2..=571.
- `modulus` is the polynomial as an integer, bit `i` the coefficient of
  `z^i`, as in ICV1's model JSON (`docs/curves/ICV1.md`). The example is
  `z^83 + z^45 + z^2 + z + 1`, AGENTS.md §8a's polynomial.
- `modulus` may be omitted only when no coordinate appears in the
  document. The tool then uses the repository's rule, `ICV1.md`'s table
  (`find_irreducible_sparse`), and writes the polynomial into the
  report.
- Coordinates are in the polynomial basis. Normal-basis coordinates, in
  which ECC2K-130 publishes its points, wait for B4 (§10).

A prime field, with `p ≥ 5`:

```json
"field": {"kind": "prime", "p": "2147483647"}
```

An extension of a prime field. B2 parses and validates it, and it is
routed only from B5 (B5a routes it: [`extension-fields.md`](extension-fields.md)):

```json
"field": {"kind": "prime_extension", "p": "1000003", "degree": 3, "modulus": ["2", "0", "1"]}
```

`modulus` lists `c_0 … c_{k−1}` of the monic polynomial
`t^k + c_{k−1} t^{k−1} + … + c_0`. An element of the field is an array
of `k` integers in the basis `1, t, …, t^{k−1}`.

### 3.2 The curve

| form | field | equation | keys |
|:--|:--|:--|:--|
| `binary_weierstrass` | binary | `y² + xy = x³ + ax² + b` | `a`, `b` |
| `koblitz` | binary | `y² + xy = x³ + ax² + 1`, `a ∈ {0, 1}` | `a` (a number) |
| `general_weierstrass` | any | `y² + a₁xy + a₃y = x³ + a₂x² + a₄x + a₆` | `a1`, `a2`, `a3`, `a4`, `a6` |
| `short_weierstrass` | prime, extension | `y² = x³ + ax + b` | `a`, `b` |
| `montgomery` | prime | `By² = x³ + Ax² + x` | `A`, `B` |
| `twisted_edwards` | prime | `ax² + y² = 1 + dx²y²` | `a`, `d` |

The tool converts `general_weierstrass`, `montgomery` and
`twisted_edwards` to `binary_weierstrass` or `short_weierstrass`:
- It uses the standard isomorphisms and birational equivalences, and
  maps the input's points the same way.
- It records the map in the report.
- A group isomorphism leaves a discrete logarithm unchanged, so the
  answer needs no mapping back.
- ICV1 names the converted model.

When it converts a binary curve and `n` is odd, the tool also takes `a`
to 0 or 1, by `y ↦ y + cx`. A curve defined over a subfield then stays
recognisable after the conversion. A `binary_weierstrass` input is run
as given. If its `a` could be normalised that way, the `curve-subfield`
disclosure says so.

A binary `general_weierstrass` curve with `a₁ = 0` is supersingular. It
is valid but unsupported (§4.3).

v2 has no form for v1's subfield curves (`subfield`, `curve_a`,
`curve_b`). Those indices refer to a basis of `GF(2^k)` that the tool
computes, which is why ICV1 retired that notation. In v2, `a` and `b`
are field elements, and the router finds the subfield itself (§5).

### 3.3 The subgroup

```json
"subgroup": {"order": "2417851639230796216685689", "cofactor": "4",
             "generator": {"x": "0x…", "y": "0x…"}}
```

- `order` is the prime `r`, and `cofactor` is `h = #E / r`.
  - Both may be omitted when the tool can find `#E` exactly (§4.5).
    The tool then divides out `h` by trial division up to `2^20`, and
    needs what remains to be a prime greater than `h`. Otherwise the
    input is refused with `subgroup-underivable`.
- `generator` is either a point or `{"rule": "koblitz_search_v1"}`.
  - The rule is v1's seeded search (`KoblitzCurve::subfield`).
  - The v1 translation needs it, to give identical outputs.
- `order_certificate` is optional. It is a Pocklington certificate for a
  prime `r` or `p` beyond the exact test's range (§4.2).

### 3.4 The target

Exactly one of:

| key | the target | its scalar |
|:--|:--|:--|
| `point` | a point `Q`, or `"identity"` | unknown to the tool |
| `public_hash_seed` | v1's public hash-to-curve point | never constructed |
| `known_log` | `Q = [k]G` with `0 < k < r` | a known answer, marked as such |
| `random_seed` | `Q = [k]G` with `k` drawn from the seed | a known answer, marked as such |

A panel of targets is a set of files, one target each. This is
AGENTS.md's single-target rule.

### 3.5 The method

```json
"method": {"solve": "paired", "fidelity": "auto",
           "index_calculus": {"pipeline": "auto", "recipe": "auto"},
           "rho": {"pipeline": "auto", "seed": 2293760}}
```

- `solve` takes one of four values:
  - `paired`, the default: the index calculus and rho on the same point,
    the programme's primary comparison;
  - `index_calculus` alone;
  - `rho` alone;
  - `check`: validate and route, and run nothing.
- `fidelity` is one of:
  - `auto`, the default: F0 when the estimate fits the budget, a
    refusal otherwise;
  - `F0`;
  - `F1`, from B7a ([`f1-sampled.md`](f1-sampled.md)): `kic` and every
    rho pipeline, extrapolated from samples on the instance. It runs only
    when asked for by name;
  - `F2`: the estimate alone.
- `pipeline` is `auto` or a pipeline id from §5.1.
- `recipe` is `auto` (§5.3), or, on `ic-gaudry-cubic`, an object of that
  pipeline's keys (B5a, [`extension-fields.md`](extension-fields.md) §2.4), or an object with v1's knobs, checked as v1
  checks them:
  - `summands`, `descent_summands`;
  - `collection_window`, `collection_aim`, `collection`;
  - `pair_table_bytes`, `pair_table_tier`;
  - `solver`, `seed`, `max_trials`;
  - `linear_algebra`, `factor_base`.
- `rho` takes v1's `baseline` knobs: `seed` and `max_iterations`.

### 3.6 The budget

```json
"budget": {"wall_seconds": 3600, "memory_bytes": 8589934592, "threads": 1}
```

- **Before the run**, the router refuses a run whose estimate exceeds
  the budget (§5.4).
- **During the run**, a run that exceeds its wall budget stops, with
  exit status 1 and the budget named on stderr. It is kept as a failure
  (plan §10).
- **The defaults** are:
  - a wall budget of 86,400 seconds, written into the report, so that a
    run of `2^128` steps is refused rather than started;
  - memory bounded by the pair table's default budget;
  - one thread.

### 3.7 Named curves

`"named": "sect163k1"` stands for the field, curve and subgroup of a
curve that a standards body or a public challenge published.
- The definition comes from the curve registry
  (`docs/curves/registry.json`) and `params::named`.
- The report writes the full definition out.
- A named curve without a recorded generator is refused with
  `named-incomplete`. `"named": "ecc2k-130"` is one today. A v2 file can
  still carry the challenge's polynomial-basis points, which
  `docs/ic/params/ecc2k130-fixed.json` already holds.

## 4. Validation

The checks run in the order below. A check whose prerequisite failed is
reported as `not_checked`, naming that prerequisite, as `ic inspect`
does. Every check is exact unless it is marked as a screen.

### 4.1 Parsing

| code | refused when |
|:--|:--|
| `file-too-large` | the file exceeds 1 MiB |
| `json-invalid` | the file is not JSON |
| `schema-version` | `schema_version` is not 2. v1 files go through §8. |
| `unknown-key` | a key outside the schema appears, at any level |
| `missing-key` | a required key is absent |
| `value-invalid` | a value has the wrong JSON type, or is not one of the values its key allows (a `kind`, `form`, `solve` or pipeline name, say) |
| `integer-syntax` | an integer is empty, signed, spaced, neither decimal nor `0x` hexadecimal, or over 1,024 bits |
| `name-syntax` | `name` is empty, over 120 characters, or has a control character |
| `target-count` | there is not exactly one target form |
| `modulus-required` | a coordinate appears and the binary modulus is omitted |
| `form-field-mismatch` | the curve's form does not belong to the field's kind |

### 4.2 The field

| code | refused when | test |
|:--|:--|:--|
| `degree-range` | a binary degree outside 2..=571, or an extension degree outside 2..=64 | |
| `modulus-degree` | the modulus's degree is not the declared one, or its constant term is 0 | |
| `modulus-reducible` | the modulus is reducible | Rabin's test: `z^{2^n} ≡ z`, and `gcd(z^{2^{n/ℓ}} − z, f) = 1` for each prime `ℓ` dividing `n`; over `F_p`, `p` takes the place of 2 |
| `p-range` | `p < 5`, or `p` has over 1,024 bits | |
| `p-composite` | `p` is not prime | Miller–Rabin with the first 13 primes as bases. It is exact below `ψ₁₃ = 3,317,044,064,679,887,385,961,981` (`2^81.46`; Sorenson and Webster). Above that it is a Baillie–PSW screen, unless a Pocklington certificate makes it exact. |

The m = 83 gate's `r = 2417851639230796216685689` (`2^81.0`) lies below
`ψ₁₃`, so its primality is exact. ECC2K-130's `r` (`2^129`) needs a
certificate.

### 4.3 The curve

| code | refused when |
|:--|:--|
| `coefficient-range` | a coefficient is not reduced: at least `2^n` (binary), or at least `p` (prime) |
| `curve-singular` | binary: `b = 0`; short Weierstrass: `4a³ + 27b² ≡ 0`; Montgomery: `B(A² − 4) ≡ 0`; twisted Edwards: `ad(a − d) ≡ 0` |

A binary curve with `a₁ = 0` is valid but unsupported, under the code
`supersingular-binary`. Its embedding degree is at most 4, so the MOV
and Frey–Rück reductions apply. The tool has no route for them.

### 4.4 The subgroup and the target

| code | refused when |
|:--|:--|
| `order-composite` | `r` is not prime, by §4.2's test |
| `certificate-invalid` | a supplied certificate does not verify |
| `generator-off-curve` | `G` is not on the curve |
| `generator-identity` | `G = O` |
| `generator-order` | `[r]G ≠ O` |
| `cofactor-mismatch` | `h·r ≠ #E`, with `#E` found as in §4.5 |
| `order-squared` | `r` divides `h`. Then `r²` divides `#E`, and a point of order `r` need not lie in `⟨G⟩`. |
| `target-off-curve` | `Q` is not on the curve |
| `target-outside-subgroup` | `[r]Q ≠ O`. Since `r² ∤ #E`, `[r]Q = O` puts `Q` in `⟨G⟩`. |
| `known-log-range` | `k` is not in `[1, r)` |

The target `"identity"` is valid. Its logarithm is 0, and its route is
`trivial` (§5.1).

### 4.5 The group order, exactly

`cofactor-mismatch` needs `#E`. The tool finds it by the first method
below that applies, and the report names the method.

1. **A Koblitz curve** (`a ∈ {0, 1}`, `b = 1`): the trace recurrence,
   `koblitz_order` in `params.rs`.
2. **A curve defined over a subfield** `GF(2^k)`. Here
   `a^{2^k} = a`, `b^{2^k} = b`, `k` divides `n`, and `k ≤ 16`.
   - The tool counts `E(GF(2^k))` by enumeration.
   - It then applies the recurrence, as `KoblitzCurve::subfield` does.
3. **A large subgroup**, `r > 4√q`.
   - `#E` is then the only multiple of `r` in the Hasse interval, and
     `h = ⌊(√q + 1)²/r⌋`. This is SEC 1's check.
   - It needs `r` prime, `G ≠ O` and `[r]G = O`.
4. **A small field**, `q ≤ 2^32`, within B0's bound. It enumerates up to
   `q = 2^16`, and above that uses Mestre's method: the orders of points
   on the curve and its quadratic twist, intersected in the Hasse
   interval. Both are exact (B2's amendment 1).

When none of these applies, the instance is unsupported, under the code
`cardinality-unknown`.
- Point counting would lift it: Schoof or SEA for prime fields, Satoh or
  AGM for binary ones.
- It rarely binds. A curve chosen for cryptography has a small cofactor,
  so method 3 applies.

### 4.6 Disclosures

These never refuse. AGENTS.md §8b requires them.

| code | what it states |
|:--|:--|
| `intermediate-subfields` | the proper divisors `d` of `n` (binary) or of `k` (an extension). For `n = 51`: `GF(2^3)` and `GF(2^17)`. |
| `curve-subfield` | the least `d` with `a, b ∈ GF(2^d)` |
| `frobenius-module` | for a prime `n`: `ord_n(2)`, and the number of irreducible factors `(n − 1)/ord_n(2)` of the cyclotomic block. 31 gives 5, so six blocks; 53, 83 and 131 give one block. |
| `embedding-degree` | the least `e ≤ 20` with `r` dividing `q^e − 1`, or "> 20" |
| `anomalous` | a prime field with `#E = p`. Smart's attack solves it in polynomial time; the tool has no route for it. |
| `koblitz-endomorphism` | `End(E) ⊇ Z[τ]`, ICV1's `end = −7`, and the signed-Frobenius classes of size `2n` that rho and the index calculus use |
| `method-uses-subfield` | set by the router when the chosen pipeline uses a proper intermediate subfield |
| `primality-screen` | a prime whose test was a screen, not a proof |
| `study-pipeline` | set when the route uses a study pipeline, with the reason: `ic-prime-s3`, since no subexponential index calculus is known for prime fields (added at B2's declaration) |

## 5. Routing

### 5.1 The pipelines

| id | arm | admits | gate today | lifted by | imported by |
|:--|:--|:--|:--|:--|:--|
| `kic` | index calculus: the pair-table pipeline `ic price` runs | binary; `a, b ∈ GF(2^k)` with `k ≤ 8`; `n = k·e`, `e` odd and at least 3 | `n ≤ 62`, `r < 2^63`; `k = 1` until B2b | B3b (`n ≤ 126`), B4 (`n ≤ 191`) | B1 (`k = 1`), B2b (`k > 1`) |
| `rho-koblitz` | rho on signed-Frobenius classes, the matched reference | Koblitz curves (`k = 1`), `n` odd and at least 3 | `n ≤ 62` | B3 (`n ≤ 126`, `r < 2^127`), B4 | B1 |
| `rho-negation` | rho with the negation map | any ordinary binary curve; prime fields; from B5a, extension fields | binary `n ≤ 62`; prime `p < 2^63`; extension `q ≤ 2^62` | B3, B5b | B2, B5a |
| `ic-binary-s4` | index calculus: `ic boundary`'s generic binary pipeline | any ordinary binary curve | `n` in 5..=32 | — | B2 |
| `ic-prime-s3` | index calculus: `ic boundary`'s prime pipeline (see below) | prime fields | `p < 2^63` | B5 | B2 |
| `rho-bignum` | rho with the negation map, on `BigUint` arithmetic | any valid instance (extension fields from B5a) | the estimate must fit the budget | — | B2, B5a |
| `ic-gaudry-cubic` | index calculus: Gaudry's, on `E(GF(p³))` (the residual-walk thread's `gaudry_cubic`) | extension fields with `k = 3`, the modulus `t³ − c` and `h = 1` | `r < 2^63` | — | B5a |
| `trivial` | none | the target `"identity"` | — | — | B1 |

- `ic-prime-s3` is a study pipeline. No subexponential index calculus is
  known for prime fields, and the report says so.
- `rho-bignum` makes every valid instance with a small enough subgroup
  runnable, whatever the field's width. It is never the matched
  reference while a faster rho admits the instance (plan §8, A5).
- The signed-Frobenius classes (`SignedFrobeniusClasses`) exist only
  for `k = 1`, and so does `ic price`. A curve over a subfield with
  `k > 1` therefore pairs `kic` with `rho-negation`, from B2. That
  reference leaves the `q`-power Frobenius unused, and the report says
  so (plan §8, A5).

### 5.2 The decision

For a valid instance:
1. The target `"identity"` routes to `trivial`.
2. The router goes through every pipeline. For each one it records
   whether the pipeline admits the instance and, if not, its gate's
   code.
3. Under `paired`, each arm takes the admitted pipeline of its kind
   with the least estimate.
   - For the rho arm that means `rho-koblitz` for a Koblitz curve,
     otherwise `rho-negation`. `rho-bignum` is used only when neither
     admits the instance.
   - If an arm has no admitted pipeline, the instance is unsupported,
     under `no-ic-route` or `no-rho-route`.
4. A pipeline named in the input is used if it admits the instance.
   Otherwise the input is refused with that pipeline's gate code.
5. **The level.**
   - The run is F0 when every arm's estimate fits the budget.
   - Otherwise it is refused with `over-budget` and the estimates. From
     B7a, the refusal also gives F1's own estimated cost. `auto` never
     falls back to F1: a run asked for as a measurement never becomes an
     extrapolation unasked.
   - The arms of a paired run share one wall budget, so the router
     compares their estimates' sum with it (B2's amendment 1).
   - `fidelity: F2` runs nothing, so it skips the size gates: the five
     word-width codes `field-wider-than-one-word`,
     `field-wider-than-two-words`, `prime-wider-than-one-word`,
     `scalar-wider-than-63-bits` and `scalar-wider-than-127-bits`. It
     returns the estimate of every pipeline whose family admits the
     instance, with `status: estimated`.

Steps 3 and 4 come before step 5. An instance with no route is
therefore `unsupported`, whatever its budget.

The gate codes:

| code | gate |
|:--|:--|
| `field-wider-than-one-word` | binary `n > 62`: on `rho-koblitz` until B3, on `kic` until B3b. From B5a, also an extension field with `q > 2^62` on `rho-negation`: its point key packs to `2(x̂ + 1) + s`, as a binary one does |
| `field-wider-than-two-words` | binary `n > 126` on `rho-koblitz`, from B3, and on `kic`, from B3b. A class key packs to `2(x + 1) + s`, which needs 129 bits at `n = 127`. Lifted by B4. |
| `prime-wider-than-one-word` | `p ≥ 2^63` (B5) |
| `scalar-wider-than-63-bits` | `r ≥ 2^63`: on `rho-koblitz` until B3, on `kic` until B3b |
| `scalar-wider-than-127-bits` | `r ≥ 2^127` on `rho-koblitz`, from B3. No curve that passes the field gate reaches it, since `r ≤ #E < 2^127` there; it is stated so that the gate is total. |
| `even-extension-degree` | `kic` and `rho-koblitz` need `e` odd. Their one-word point lifting solves `z² + z = c` by the half-trace, and the orbit maps assume the distinct factors of `x^n − 1` that an odd `n` gives (`koblitz_index_calculus.rs`). |
| `subfield-curve-unsupported` | `k > 1` on `kic` before B2b, and on `rho-koblitz` always |
| `subfield-too-large` | `k > 8` |
| `not-a-subfield-curve` | no `k ≤ 8` with `a, b ∈ GF(2^k)` |
| `enumeration-bound` | `ic-binary-s4` at `n > 32` |
| `no-pipeline-for-field` | an extension field before B5a, or characteristic 3; and, per pipeline, a field kind the pipeline has no implementation for (B2's amendment 1) |
| `extension-degree-not-three` | `ic-gaudry-cubic` on `k ≠ 3` (B5a) |
| `modulus-not-binomial` | `ic-gaudry-cubic` on a cubic modulus other than `t³ − c`: the module's field is `GF(p)[t]/(t³ − c)` (B5a) |
| `cofactor-not-one` | `ic-gaudry-cubic` on `h > 1`: its relations are taken modulo `#E`, every point in `⟨G⟩` (B5a) |
| `recipe-not-taken` | a recipe object, which holds `kic`'s knobs, on any other index calculus pipeline (B2's amendment 1); on `ic-gaudry-cubic`, a recipe object with a key outside that pipeline's (B5a) |
| `no-recipe` | `recipe: auto` where §5.3 has no rule |
| `subgroup-smaller-than-cofactor` | `kic` and `rho-koblitz` need `r > h` at one word, as `KoblitzCurve`'s own construction does; past one word the refusal is lifted (B3b, amendment 1) |
| `not-the-identity` | `trivial` takes the identity only |
| `over-budget` | an estimate exceeds the budget |
| `no-f1-model` | `fidelity: F1` on a pipeline with no F1 model: `ic-binary-s4` and `ic-prime-s3` (B7a) |
| `not-yet-supported` | a feature of the schema that a later step implements. The message names the feature and the step. B1 uses it for prime and extension fields, the other curve forms, named curves, one arm alone, `recipe: auto`, F1, F2 and more than one thread. From B3, `solve: rho` runs. Under it, v1's hashed and random targets past one word wait for B4, and so does v1's generator rule. B2 lifts all of B1's list except F1 (B7) and more than one thread (plan §8, A6). B2's amendment 1 leaves `kic` alone, and `kic` beside any rho but `rho-koblitz`, to B2b. |

### 5.3 Recipes

`recipe: auto` gives the index calculus a recipe that the user need not
tune:
- **`kic` on a Koblitz curve:**
  - §20's rules (`research/ic_exponent_20260926/make_params.py`), at
    the column count and descent summands that §20's model gives for
    `(r, n)` (`predict.py`, `optimum`);
  - at the suite's sizes, the suite's recipe, so that `auto` there is
    the suite's row.
- **`kic` on a subfield curve with `k > 1`:** §20's model does not
  apply. `auto` is refused with `no-recipe`. The user gives a recipe, or
  asks for `factor_base: {"mode": "search"}`.
- **The pipelines B2 imports** (`ic-binary-s4`, `ic-prime-s3`): each
  pipeline's own default from `ic boundary`, named in the report.

### 5.4 Estimates

Every route carries an F2 estimate for each arm, in seconds on the
reference host, labelled as a model:
- **Rho.**
  - The expected number of steps is `√(πr/(2A))`, where `A` is the
    class size: `2n` for signed Frobenius, 2 for the negation map.
  - It is multiplied by that pipeline's step cost at the nearest
    measured size, from the newest baseline's rows (`baselines.json`).
  - `rho-bignum` uses its own step cost, which B2 measures.
- **The index calculus on `kic`.**
  - §20's phase model (`predict.py`, `ic_phases`), at the recipe's `|F|`
    and `m`, converted at the newest baseline's unit.
  - B2 reports this model's error against v0's measured cold times at
    the eleven suite sizes.
- **The other pipelines.** Each `ic boundary` pipeline's fitted
  exponent, from the boundary ledger.

The router compares these estimates with the budget. They are not
measurements, and the report keeps them apart from the run's own
figures.

**The keys** (fixed at B2's declaration, for the cases to read):
- `estimate.<arm>` is `{pipeline, seconds, model}`.
- For rho it also carries `steps`, `√(πr/2A)`, and `steps_bits`, its
  `⌈log₂⌉`.
- `model` names the step cost or the phase model the seconds come
  from.

## 6. The report

Every operation that reads a v2 document adds these keys to its report:

| key | what |
|:--|:--|
| `input` | the file's SHA-256 |
| `resolved` | the instance as the tool understood it, as a complete v2 document (below) |
| `curve_id` | the resolved model's ICV1 identity: `{icv1, slug, model_json}` (AGENTS.md §11). The EC1 alias of the representation is computed from `resolved` by `tools/curve_identity.py`, which has no Rust port yet. |
| `checks` | every check in §4, as `{code, status, exact, details}`. `status` is `pass`, `fail` or `not_checked`; `exact` is true or false. |
| `disclosures` | §4.6 |
| `route` | the chosen pipeline for each arm, the level, and every pipeline considered, with its gate code |
| `estimate` | §5.4, for each arm |
| `refusal` | when the run is refused: `{code, class, message}`, where `class` is `invalid`, `unsupported` or `over_budget` |
| `result` | after a run: `{scalar, verified, known_answer}`, the scalar both arms recovered, in decimal |
| `conversion` | for a converted curve: `{from, to, map}`. `from` is the input's form; `to` is `short_weierstrass` or `binary_weierstrass`; `map` states the substitution, points included (added at B2's declaration) |

The cases (§9) read these keys by path:
- `route` is `{ic, rho, level, considered}`.
  - `ic` and `rho` are each `{pipeline}`.
  - `considered` lists `{pipeline, arm, admitted, gate}`, one entry for
    every pipeline.
- Each disclosure is `{code, value}`.
- `curve_id` is `{icv1, slug, model_json}`.
- A run keeps v1's keys as well: `status`, `counts`, `certificates`,
  `rho_counts`, `all_verified` and `ic_and_rho_agree`. A pin therefore
  compares a v2 run with a v1 run directly.

In `resolved`:
- every default is filled in;
- every conversion is applied and recorded;
- the modulus is written out.

Run as an input, `resolved` gives the same result as the original.

**Exit statuses:**

| status | meaning |
|:--|:--|
| 0 | complete and verified; or `check` found the input valid |
| 1 | the run failed: a timeout, an unverified scalar, or out of memory |
| 2 | invalid; this is also clap's status for bad usage |
| 3 | unsupported |
| 4 | over budget |

A `prime_extension` instance had no ICV1 identity until B5a: ICV1
defined binary and prime models only. B5a's declaration extends ICV1 to
`GF(p^k)` and registers its instances (AGENTS.md §11), before the tool
writes any extension-field result.

## 7. Entry points

- **`ic check --params <file>`** validates, routes, estimates and stops.
  - It also reads inspection v1 files, so it supersedes `ic inspect`.
  - `inspect` stays, unchanged, for compatibility.
  - It exits 0 for a valid input even when no route admits it. The route
    and its codes are in the report.
- **`ic price --params <file>` and `ic workflow`** read v1 or v2,
  according to `schema_version`. Under v2, `price` is single-target by
  construction.
- **`--curve <name>`**, on `check` and `price`, builds a document with
  `named`.

## 8. v1

- **Workflow v1 files** keep their parser, their meaning and their
  outputs. B1 adds a translation from a v1 file to v2:

  | v2 part | taken from the v1 file |
  |:--|:--|
  | field | `binary`, with the repository's modulus written out |
  | curve | `a` and `b` in the polynomial basis, computed from v1's indices through the subfield basis the tool builds |
  | subgroup | `{"rule": "koblitz_search_v1"}`, with `r` and `h` written out |
  | target | the file's only target |
  | method | `solve: paired`; `index_calculus: {"pipeline": "kic", "recipe": …}` with v1's knobs copied; `rho: {"pipeline": "rho-koblitz", …}` with v1's rho seed and iteration cap |

  The pipelines are named, not `auto`, so that a later step adding a
  pipeline cannot change a translated row's route.

  The translation must give identical outputs, and B1 checks that on
  every suite row. A v1 file with several targets becomes one v2
  document per target.
- **Inspection v1 files** translate to `solve: check`.
- **`ic fixed`'s files** are left as they are. B4 replaces them with v2
  documents for `K_0` at `n = 131`, with polynomial-basis points.

## 9. Conformance cases

Conformance suite v1 (`../conformance/v1/`) holds C001–C008. The cases
below have their expected outcomes fixed here.
- **B1's cases (C009–C031) are frozen now**, in
  [`../conformance/v2/cases.json`](../conformance/v2/cases.json).
  - `../conformance/v2/make_cases.py` writes their parameter files and
    `SHA256SUMS`.
  - It checks each file's intended property with its own arithmetic,
    which shares nothing with the Rust.
- **B2's cases (C032–C051)** have their inputs described below. Their
  files were written at B2's declaration (2026-10-01), before any B2
  code, by `../conformance/v2-b2/make_cases.py`. B1's own generator and
  files are frozen, so B2's sit beside them. Their expected outcomes do
  not change, and the declaration made them exact as report keys
  (`conversion`, `estimate.<arm>`, `study-pipeline`).
- **Steps.** Each case names the step that must make it pass. A step is
  accepted only when its own cases and every earlier step's cases pass.
- **`until`.** A case that expects a refusal a later step will lift
  names that step in `until`. That step changes the case's expectation
  with a dated note in `cases.json`, and adds the newly possible run as
  a case of its own.
- **`supersedes`** (B3b's amendment 2). A step can move an expectation
  without being the case's `until` step. Its own case then names the
  earlier case it `supersedes`, and while that case is run the earlier
  one is not. B3b supersedes C031 and C053 this way, since it moves
  `kic`'s width gate. These two rules are the only ways an expectation
  changes.

The curves the cases use:
- **Curve A** is `icv1-f2m31-tm90707-c95f16f5`:
  `y² + xy = x³ + 1` over `GF(2^31)` with the modulus `z^31 + z^3 + 1`.
  It is the smoke curve: `r = 1439393`, `h = 1492`.
- **Curve A′** is the same equation with the modulus `z^31 + z^6 + 1`:
  another model of the same abstract curve.
- **The gate curve** is `icv1-f2m83-tm6151469093347-debefd74`, with
  AGENTS.md §8a's modulus.

| case | input | expected | step |
|:--|:--|:--|:--|
| C009 | the v2 translation of the smoke row `M1-T01` | exit 0, with the same outputs as the v1 file | B1 |
| C010 | curve A, written out, with another generator `G′` and a `known_log` target | exit 0, recovering `known_log` | B1 |
| C011 | curve A′, its own generator, a `known_log` target | exit 0, recovering `known_log`; `curve_id` is A′'s model, not A's | B1 |
| C012 | an unknown key inside `subgroup` | exit 2, `unknown-key` | B1 |
| C013 | `schema_version: 3` | exit 2, `schema-version` | B1 |
| C014 | two target forms | exit 2, `target-count` | B1 |
| C015 | a generator, and no modulus | exit 2, `modulus-required` | B1 |
| C016 | `"order": "+1439393"` | exit 2, `integer-syntax` | B1 |
| C017 | the modulus `z^31 + 1` | exit 2, `modulus-reducible` | B1 |
| C018 | `b = 0` | exit 2, `curve-singular` | B1 |
| C019 | a generator off the curve | exit 2, `generator-off-curve` | B1 |
| C020 | a generator on the curve, outside the order-`r` subgroup | exit 2, `generator-order` | B1 |
| C021 | `order` = `3r`, a composite | exit 2, `order-composite` | B1 |
| C022 | `cofactor` = 1491 | exit 2, `cofactor-mismatch` | B1 |
| C023 | a target on the curve, outside the order-`r` subgroup | exit 2, `target-outside-subgroup` | B1 |
| C024 | a target off the curve | exit 2, `target-off-curve` | B1 |
| C025 | `known_log` = `r` | exit 2, `known-log-range` | B1 |
| C026 | the target `"identity"` | exit 0, recovering 0 by the `trivial` route | B1 |
| C027 | the gate curve with a frozen generator and a public target, under `price` | exit 3, `no-ic-route`, with `field-wider-than-one-word` on `kic` and `rho-koblitz`; until B3 | B1 |
| C028 | the same file under `check` | exit 0; `order-composite` passes exactly; disclosures include `frobenius-module` (one block) | B1 |
| C029 | the suite's `icv1-f2m45-tm6236725-40939294` in v2, with a generator and a `known_log` target | exit 0, verified; disclosures include `intermediate-subfields` (3, 5, 9, 15) | B1 |
| C030 | `y² + xy = x³ + 1` over `GF(2^32)`, its modulus in the file | exit 3, `no-ic-route`, with `even-extension-degree` on `kic` | B1 |
| C031 | ECC2K-130, with the challenge's polynomial-basis generator and target, under `check` | exit 0; `primality-screen` for `r`; `field-wider-than-one-word` on `kic` until B4 | B1 |
| C032 | a prime curve, `p ≈ 2^24`, `known_log` target, `solve: rho` | exit 0, `rho-negation`, recovering `known_log` | B2 |
| C033 | the same, `solve: paired` | exit 0, `ic-prime-s3` and `rho-negation`, both verified, with the note that no subexponential index calculus is known | B2 |
| C034 | a Montgomery curve, `p ≈ 2^24`, `solve: rho` | exit 0, verified; the conversion recorded | B2 |
| C035 | a twisted Edwards curve, `p ≈ 2^24`, `solve: rho` | exit 0, verified; the conversion recorded | B2 |
| C036 | binary `general_weierstrass` with `a₁ = 0` | exit 3, `supersingular-binary` | B2 |
| C037 | binary `general_weierstrass` with `a₁ ≠ 0`, Koblitz after conversion | exit 0, `kic`, verified | B2 |
| C038 | a generic binary curve, `n = 23`, `paired` | exit 0, `ic-binary-s4` and `rho-negation`, verified | B2 |
| C039 | a generic binary curve, `n = 37`, `paired` | exit 3, `no-ic-route` with `enumeration-bound` | B2 |
| C040 | the same file, `solve: rho` | exit 0, `rho-negation`, verified | B2 |
| C041 | a composite `p` | exit 2, `p-composite` | B2 |
| C042 | `4a³ + 27b² ≡ 0` | exit 2, `curve-singular` | B2 |
| C043 | a prime curve with `r²` dividing `#E` | exit 2, `order-squared` | B2 |
| C044 | a prime curve, `p ≈ 2^34`, `r < 4√p` | exit 3, `cardinality-unknown` | B2 |
| C045 | an anomalous prime curve, `p ≈ 2^20`, `solve: rho` | exit 0, verified; disclosure `anomalous` | B2 |
| C046 | a curve over `GF(2^4)` taken over `GF(2^100)`, with a subgroup of about `2^20`, `solve: rho` | exit 0, `rho-bignum`, verified | B2 |
| C047 | `named: secp256k1`, with a target point, `solve: rho` | exit 4, `over-budget`, the estimate near `2^128` steps | B2 |
| C048 | the v2 translation of a suite row at `2^47.2`, with `wall_seconds: 1` | exit 4, `over-budget` | B2 |
| C049 | the gate file with `fidelity: F2` | exit 0, `status: estimated`, with estimates for `kic` and `rho-koblitz` | B2 |
| C050 | a valid `F_{p^3}` curve | exit 3, `no-pipeline-for-field` | B2 |
| C051 | an `F_{p^3}` modulus that is reducible | exit 2, `modulus-reducible` | B2 |

## 10. Steps and their acceptance

| step | delivers | accepted when |
|:--|:--|:--|
| **B1** | <ul><li>schema v2's parser (§3)</li><li>the checks a binary instance needs (§4.1–§4.2 for binary, §4.3, §4.4, §4.5 methods 1–2)</li><li>the disclosures for binary fields</li><li>`kic`, `rho-koblitz` and `trivial` reading imported instances</li><li>the v1 translation</li><li>`ic check`</li></ul> Estimates, `recipe: auto` and pre-run budget refusals come with B2. Under B1, `fidelity: auto` means F0, and the budget is enforced while the run goes. | <ul><li>C001–C031 pass</li><li>the pin holds on all 90 rows</li><li>every row's v2 translation gives identical outputs</li><li>no size regresses beyond its A/A band</li></ul> |
| **B2** | <ul><li>prime and extension fields</li><li>the other curve forms and their conversions</li><li>§4.5 methods 3–4</li><li>the importers for `rho-negation`, `ic-binary-s4`, `ic-prime-s3` and `rho-bignum`</li><li>estimates and budgets</li></ul> Declared 2026-10-01 in [`../rounds/B2-fields-forms-estimates/PROTOCOL.md`](../rounds/B2-fields-forms-estimates/PROTOCOL.md); its cases are in [`../conformance/v2-b2/`](../conformance/v2-b2/cases.json), since B1's files are frozen. | <ul><li>C001–C051 pass</li><li>the same pin, translation and timing checks as B1</li><li>the estimate's error reported at the suite's sizes</li></ul> |

B3–B7 are unchanged from the plan, with these refinements:
- **B3** (declared 2026-10-01 in
  [`../rounds/B3-two-word-rho/PROTOCOL.md`](../rounds/B3-two-word-rho/PROTOCOL.md))
  lifts `rho-koblitz` to `n ≤ 126` and `r < 2^127`, and adds
  `solve: rho`.
  - Its cases are [`../conformance/v2-b3/`](../conformance/v2-b3/cases.json),
    C052–C058.
  - C027 is retired by the `until` rule, and C054 succeeds it: `paired`
    on the gate file is refused with `no-ic-route`, and `rho-koblitz` is
    admitted.
  - The plan's gate run is B3's measurement 5: the gate file with
    `solve: rho`, at F0.
  - B1's files are frozen, so the dated note C027 needs is in B3's
    `cases.json`, and the programme's runner, `../conformance/run.py`,
    applies the rule.
- **B2b** (B2's amendment 1) takes three items from B2:
  - `kic` on curves over a subfield with `k > 1`, paired with
    `rho-negation`;
  - `kic` alone;
  - the verification of order certificates.

  It was declared 2026-10-01 in
  [`../rounds/B2b-subfield-kic-certificates/PROTOCOL.md`](../rounds/B2b-subfield-kic-certificates/PROTOCOL.md).
  Its cases are C059–C070, in
  [`../conformance/v2-b2b/`](../conformance/v2-b2b/cases.json). The
  declaration makes two points of this design exact:
  - §3.3's certificate for `p` is the key `field.p_certificate`, and
    `subgroup.order_certificate` certifies `r`. Both have the form and
    rules the protocol states.
  - `kic`'s estimate on a curve over a subfield is §20's model with `e`
    in the place of `n`, labelled as outside the range the model was
    fitted on.
- **B3b** lifts `kic` to `n ≤ 126` and runs the index calculus at F1.
  C054 then changes by the `until` rule: at `n = 83` the index
  calculus's estimate exceeds any day-long budget, so `paired` becomes
  `over-budget`. (This was B3's own refinement before B3 was split.)
  Its amendment 2 adds C086 and C087, which supersede C031 and C053:
  from B3b `kic` refuses `n = 131` as wider than two words and admits
  `n = 83`.
- **`rho-negation` on two-word fields** follows its importer, which is
  B2's.
- **B4** adds the normal-basis import with the challenge.
- **B5** is split in two (B5a's declaration, 2026-10-01).
  - **B5a** routes `F_{p^k}`: `rho-negation` and `rho-bignum` on it, and
    `ic-gaudry-cubic` on `E(GF(p³))`. Its declaration extends ICV1 first
    (`docs/curves/ICV1.md`). Design: [`extension-fields.md`](extension-fields.md);
    protocol: [`../rounds/B5a-extension-fields/PROTOCOL.md`](../rounds/B5a-extension-fields/PROTOCOL.md).
  - **B5b** lifts prime fields past one word, declared on its own later.

## 11. Deferred, and why

- **Normal-basis coordinates.** The only consumer is ECC2K-130, which
  is B4's. `ic fixed` already holds the conversion (`onb_root`).
- **Point counting beyond §4.5.** It is needed only for curves with a
  large cofactor that no other method covers.
- **Characteristic 3.** No pipeline exists. It is refused with
  `no-pipeline-for-field`.
- **Pairing transfer and Smart's attack.** They are disclosed, not
  routed. The tool's subject is the index calculus against rho.
