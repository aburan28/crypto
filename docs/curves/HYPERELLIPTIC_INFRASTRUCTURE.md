# Hyperelliptic cover infrastructure

`hyperelliptic-cover` is a bounded, certificate-oriented command-line tool.
It constructs and verifies the repository's two explicit cover families and
runs checked Jacobian arithmetic only when the infinity configuration matches
an implemented representation. It is not an index-calculus or key-recovery
solver.

The additive catalog/evidence schema is
[`hyperelliptic-infrastructure.schema.json`](hyperelliptic-infrastructure.schema.json).
That schema records EC1/ICV1 identity, source curves, Jacobian contexts,
morphisms, induced maps, evidence, and unsupported cases. It is not the CLI
request schema. The CLI wire contract is specified below.

## Direct packaged invocation

After unpacking a release artifact, invoke the executable directly:

```text
./hyperelliptic-cover --help
```

A local release build has the equivalent path:

```text
./target/release/hyperelliptic-cover --help
```

Users of a packaged artifact do not need `cargo run`. The four command shapes
are:

```text
./hyperelliptic-cover construct --input MODEL.json [--output RESPONSE.json]
./hyperelliptic-cover verify --input VERIFY.json [--output RESPONSE.json]
./hyperelliptic-cover arithmetic --input OPERATION.json [--output RESPONSE.json]
./hyperelliptic-cover catalog-export --registry REGISTRY.json [--output covers.json]
```

For the first three commands, `--input -` reads standard input. Omitting
`--output` writes JSON to standard output. Successful `catalog-export` writes
the existing deterministic `curve-covers/v1` catalog shape; it does not
silently claim enriched Jacobian or transfer evidence. `--registry` is
required: release archives contain the executable, not a registry or catalog.
All input paths, including the catalog registry, are capped at 32 MiB. A
larger file or standard-input stream returns `invalid_input` rather than being
partially parsed.

Command response status and process exit status are distinct and stable:

| JSON `status` | Exit | Meaning |
|---|---:|---|
| `ok` | 0 | Requested bounded operation completed |
| `unsupported` | 2 | Well-formed request is outside the implemented mathematical boundary |
| `invalid_input` | 3 | Malformed JSON, model, certificate, divisor, or operation arguments |
| `error` | 1 | Internal certificate replay or output failure, including a failed stdout write |

Every enveloped response has:

```json
{
  "schema": "hyperelliptic-infrastructure-command/v1",
  "command": "construct | verify | arithmetic",
  "status": "ok | unsupported | invalid_input | error",
  "result": {}
}
```

The strings separated by `|` above document alternatives; an actual response
contains one string.

## Construct and verify

`construct` accepts a supported elliptic-model JSON shape. It does not look
up a slug or consult a registry. A registry's `model_json` value is one way to
obtain such an object; the surrounding registry row is not accepted. Unknown
model fields are rejected. For example, the model represented by
`icv1-f2m7-t13-616700dd` uses this complete input:

```json
{
  "a": "0x0",
  "b": "0x1",
  "field": "f2m-7-e1a6f45a",
  "form": "y^2+xy=x^3+a*x^2+b",
  "modulus": "0x83",
  "v": "1"
}
```

Run:

```text
./hyperelliptic-cover construct --input f2m7-model.json --output f2m7-cover.json
```

The command constructs the fixed family and immediately replays its geometry
certificate. A verification request contains exactly `model` and
`certificate`. For the same binary fixture, this is a complete request:

```json
{
  "model": {
    "a": "0x0",
    "b": "0x1",
    "field": "f2m-7-e1a6f45a",
    "form": "y^2+xy=x^3+a*x^2+b",
    "modulus": "0x83",
    "v": "1"
  },
  "certificate": {
    "construction": "binary_cubic_pullback_v1",
    "degree": 3,
    "f": ["0x0", "0x1", "0x0", "0x0", "0x0", "0x0", "0x0", "0x1"],
    "genus": 3,
    "h": ["0x0", "0x0", "0x1"],
    "x": ["0x0", "0x0", "0x0", "0x1"],
    "y_0": ["0x1"],
    "y_v": ["0x0", "0x1"]
  }
}
```

```text
./hyperelliptic-cover verify --input f2m7-verify.json --output f2m7-verified.json
```

The supported constructors are deliberately separate:

- odd characteristic greater than three: the fixed quadratic pullback
  `prime_quadratic_pullback_v1`, producing a genus-two monic sextic with two
  rational infinity points; and
- characteristic two: the fixed ordinary binary cubic pullback
  `binary_cubic_pullback_v1`, producing a genus-three odd-degree model with
  one rational infinity point.

An unsuccessful bounded construction never proves that no cover exists.
The `construct`, `verify`, and arithmetic `model` paths reject additional
properties even when the required coefficients are otherwise valid.

## Arithmetic request contract

An arithmetic request is a JSON object with no unknown fields and the
following members:

| Member | Required | Contract |
|---|---|---|
| `schema` | yes | Exact string `hyperelliptic-jacobian-operation/v1` |
| `model` | exactly one curve selector | Exact binary or prime elliptic model object accepted by `construct` |
| `odd_curve` | exactly one curve selector | Separate odd-prime, odd-degree, one-infinity hyperelliptic curve object |
| `operation` | yes | `validate`, `identity`, `equal`, `negate`, `add`, or `scalar_multiply` |
| `left` | by operation | Required by `validate`, `equal`, `negate`, `add`, and `scalar_multiply` |
| `right` | by operation | Required by `equal` and `add` |
| `scalar` | by operation | Required by `scalar_multiply`; unsigned decimal or `0x` string |

`left` and `right` are canonical Mumford pairs:

```json
{
  "u": ["coefficient of x^0", "coefficient of x^1"],
  "v": ["coefficient of x^0", "coefficient of x^1"]
}
```

Coefficient arrays are in ascending-power order. The zero polynomial is the
empty array. Coefficients must be canonical field elements; trailing zeroes
are rejected. In particular, binary Mumford `u` and `v` polynomials must not
carry a redundant final zero coefficient. The identity divisor is
`{"u":["0x1"],"v":[]}`.

The two curve selectors are intentionally not interchangeable.

### Binary cover arithmetic

For `model`, the tool verifies the elliptic model, reconstructs its fixed
cover, verifies the cover certificate, and then initializes the checked
one-infinity binary Jacobian. This complete request returns the identity on
the F2^7 catalog cover:

```json
{
  "schema": "hyperelliptic-jacobian-operation/v1",
  "model": {
    "a": "0x0",
    "b": "0x1",
    "field": "f2m-7-e1a6f45a",
    "form": "y^2+xy=x^3+a*x^2+b",
    "modulus": "0x83",
    "v": "1"
  },
  "operation": "identity"
}
```

An addition request uses the same `schema` and `model`, plus:

```json
{
  "operation": "add",
  "left": {"u": ["0x1"], "v": []},
  "right": {"u": ["0x1"], "v": []}
}
```

Those members are illustrative fragments to merge into the complete request;
the first binary block is directly runnable JSON.

### Separate odd-prime one-infinity arithmetic

`odd_curve` has exactly these fields:

| Member | Contract |
|---|---|
| `family` | Exact string `odd_prime_one_infinity_v1` |
| `p` | Unsigned odd prime of at most 64 bits; primality is checked deterministically over the complete supported range |
| `genus` | Positive integer |
| `f` | Ascending coefficients of monic squarefree `f`, with exact degree `2g+1` |
| `infinity` | Exact string `one_rational` |

This complete request returns the identity for
`y^2=x^5+x^2+1` over F3:

```json
{
  "schema": "hyperelliptic-jacobian-operation/v1",
  "odd_curve": {
    "family": "odd_prime_one_infinity_v1",
    "p": "3",
    "genus": 2,
    "f": ["0x1", "0x0", "0x1", "0x0", "0x0", "0x1"],
    "infinity": "one_rational"
  },
  "operation": "identity"
}
```

Its complete successful response is:

```json
{
  "schema": "hyperelliptic-infrastructure-command/v1",
  "command": "arithmetic",
  "status": "ok",
  "result": {
    "operation": "identity",
    "context": {
      "family": "odd_prime_one_infinity_v1",
      "field": {
        "characteristic": "3",
        "representation": "prime_residue"
      },
      "curve": {
        "family": "odd_prime_one_infinity_v1",
        "p": "3",
        "genus": 2,
        "f": ["0x1", "0x0", "0x1", "0x0", "0x0", "0x1"],
        "infinity": "one_rational"
      },
      "infinity_configuration": "one_rational_point"
    },
    "result": {
      "divisor": {
        "u": ["0x1"],
        "v": []
      }
    }
  }
}
```

## Prime even-sextic boundary

Passing a supported prime elliptic `model` to `arithmetic` does not apply the
odd-degree code to its constructed cover. The cover is a monic sextic with
two rational infinity points and needs balanced/real-model arithmetic. The
command returns exit 2 and this structured result:

```json
{
  "schema": "hyperelliptic-infrastructure-command/v1",
  "command": "arithmetic",
  "status": "unsupported",
  "result": {
    "reason": "the explicit prime cover is an even-degree sextic with two rational points at infinity; the checked odd-characteristic implementation currently supports only odd-degree one-infinity models",
    "construction_status": "explicit_verified",
    "arithmetic_status": "unsupported_infinity_configuration"
  }
}
```

Generic odd-degree/one-infinity prime arithmetic is available only through
the separate `odd_curve` selector. That does not turn the prime cover into an
odd-degree model and is not a transfer certificate for that cover.

### Prime validation and diagnostic bounds

`HyperellipticCurveP` accepts only deterministically validated odd primes of
at most 64 bits. This replaces the earlier probable-prime wording for the
feature implementation. Its optional quadratic-extension reference backend
has the narrower `p <= 4096` boundary, matching the bounded quadratic point
counter: `Fp2Ctx` rechecks that `p` is an odd prime, constructs and validates
a canonical nonzero quadratic nonresidue, and every element operation
canonicalizes both coordinates before use. Forged public element fields
therefore cannot bypass reduction.

The exhaustive point counters are diagnostics, not general large-field
algorithms:

- linear `F_p` enumeration requires `p <= 1,000,000`; and
- quadratic `F_{p^2}` enumeration and `Fp2Ctx` require `p <= 4096`.

Requests outside a checked bound return a structured error where the API is
fallible; convenience wrappers may panic and are not public-input entry
points. These counting bounds do not restrict ordinary Jacobian arithmetic
within its declared 64-bit prime-field contract.

## Transfer boundary and evidence records

The library feature snapshot has a bounded binary transfer for rational
elliptic point classes, rational source-point divisors, and its own typed
pullback fibers. It checks infinity, the exceptional ramified two-torsion
fiber, the degree-three composition, and the exact F2^7/order-29 catalog
subgroup. It does not expose a CLI transfer subcommand and does not implement
the norm of an arbitrary Mumford class. That operation returns a typed
`Unsupported(ArbitraryJacobianPushforward)` result in the library API.

The generic subgroup classifier proves only the conditional statement
“pullback is injective if the supplied integer is independently known to be
the exact subgroup order and is coprime to 3.” The CLI accepts no subgroup
order and makes no unconditional subgroup claim. The F2^7 evidence fixture is
stronger because it independently checks the named generator is nonidentity,
is killed by the prime 29, and therefore has exact order 29 before applying
the degree-inverse argument.

The conforming machine-readable envelope is
[`../../research/hyperelliptic_cover_infrastructure_20261006/evidence.json`](../../research/hyperelliptic_cover_infrastructure_20261006/evidence.json).
Its F2^7 record binds all of the following:

- ICV1 `icv1-f2m7-t13-616700dd` and EC1
  `EC1N7Ce0hb6d297a2ca08`;
- the exact source `H`, `Jac(H)`, curve morphism, field identity, and bounded
  induced-map domains;
- the supported degree-three composition and exact order-29 subgroup check;
- content hashes and local test receipts; and
- a structured refutation of arbitrary-class pushforward support.

The same envelope records the prime F827 cover as an explicit verified
construction while marking use of one-infinity arithmetic on its even sextic
as `refuted`. The original PR #1383 artifact contains 115 rows (28 prime and
87 binary); the later audit-baseline catalog contains 121. These populations
must not be conflated.

## Catalog export

Rebuild the existing cover catalog directly from the exact registry:

```text
./hyperelliptic-cover catalog-export \
  --registry docs/curves/registry.json \
  --output covers.json
```

The exporter preserves exact coefficients and curve identifiers. Schema
validation alone is not a mathematical certificate: consumers must also
recompute EC1/ICV1 identities, resolve evidence and implementation references,
replay the cover equations and infinity checks, and enforce each declared
domain boundary.

The exported catalog identifies its wrapper with
`generated_by = "hyperelliptic-cover catalog-export"`. Its
`producer_source_sha256` is SHA-256 of the exact UTF-8 checker source followed
immediately by the exact UTF-8 `hyperelliptic-cover` wrapper source, in that
order. For the documented feature snapshot it is
`3927efaddd7e7e966bb284affc56ac22a823f98f498017509e266b603b244e37`.
Thus changing either producer changes the provenance digest. Both
`hyperelliptic-cover catalog-export` and `curve_cover_check` require callers
to supply an external registry; packaged executables do not embed one.
