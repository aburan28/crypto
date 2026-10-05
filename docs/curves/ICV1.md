# ICV1: how this repository names an elliptic curve

Every elliptic curve the repository names, in text, tables, reports, code
and new file names, is named by its **ICV1 identity**.  Before this, one
curve went by five names (`K_0 / GF(2^41)`, `K_0/2^41`, `K_0/F_2^41`,
`K₀/GF(2^41)`, `k0n41`), a prime curve by the generator that happened to
find it (`bench-20bit`, `generated-24bit-10935329`), and nothing said
which modulus, basis or model a name meant.  A number quoted against the
wrong curve is a wrong number, so the name is now computed from the curve
rather than chosen.

The identity was first written for the isogeny-volcano campaign
([`research/isogeny_volcano_ic_20260924/README.md`](../../research/isogeny_volcano_ic_20260924/README.md));
this document pins the parts that campaign left open, above all the model
JSON, and makes it the rule for the whole repository (`AGENTS.md` §11).

## The identity

```
ICV1:<field>:<trace>:<order>:<j>:<end>:<level>:<path>:<model12>
```

| part | binary field `GF(2^m)` | prime field `GF(p)` |
|:--|:--|:--|
| `field` | `f2m-<m>-<modhash8>`, `modhash8` the first 8 hex digits of `SHA-256("f2m-modulus:" + <modulus hex>)` | `fp-<p>`, decimal |
| `trace` | `2^m + 1 − #E`, signed decimal | `p + 1 − #E`, signed decimal |
| `order` | `#E(GF(2^m))`, the whole group, decimal | `#E(GF(p))`, decimal |
| `j` | `1/b` in the polynomial basis, lower-case hex (`y² + xy = x³ + ax² + b`) | `1728 · 4a³ / (4a³ + 27b²) mod p`, decimal |
| `end` | discriminant of `End(E)` when certified, else `unk` | the same |
| `level` | level in a certified ℓ-volcano, else `unk` | the same |
| `path` | root-relative isogeny path, `r` for a curve not placed in one | the same |
| `model12` | the first 12 hex digits of `SHA-256(model JSON)` | the same |

The **slug** is the name used in text:

```
icv1-<field tag>-t<trace>-<model8>        field tag: f2m<m> or fp<bits of p>
icv1-f2m41-tm2308219-7f48b14a             a negative trace is written tm<|t|>
```

The slug is display; the canonical record is the ICV1 string and the model
JSON beside it.  Never identify a curve by `j` or by its order alone:
twists share `j`, and an isogeny class shares the order.

### The model JSON

The model hash is taken over the curve's defining data, serialised with
sorted keys, no whitespace and ASCII escapes (`json.dumps(obj,
sort_keys=True, separators=(",", ":"), ensure_ascii=True)`):

```
binary: {"a":"0x0","b":"0x1","field":"f2m-41-39c74c32","form":"y^2+xy=x^3+a*x^2+b","modulus":"0x20000000009","v":"1"}
prime:  {"a":"5320418","b":"8535318","field":"fp-10935329","form":"y^2=x^3+a*x+b","p":"10935329","v":"1"}
```

Binary coefficients are polynomial-basis field elements in lower-case hex
(bit `i` is the coefficient of `x^i`), reduced modulo the modulus; prime
coefficients are decimal and reduced modulo `p`.  **The modulus is part of
the model.**  One abstract curve under two moduli, or in a normal basis,
is two models and two identities.  That is deliberate: a factor base, a
pair table or a log database is a set of field elements, and it means
nothing under another modulus.

### The fields every generator uses

| family | modulus |
|:--|:--|
| Koblitz `K_a`, `n < 64` | the least `x^n + low` with `low` odd of weight at most four that is irreducible (`koblitz_index_calculus::find_irreducible_sparse`) — what `KoblitzCurve::new` and `ic` build |
| random binary curves (`ic_boundary::random_binary_instance`) | the same rule |
| a standard or challenge curve | its published polynomial: ECC2K-95 `x^97+x^6+1`, ECC2K-130 `x^131+x^13+x^2+x+1`, sect163k1 `x^163+x^7+x^6+x^3+1`, and SEC 2's for the other `sect*k1` |
| the m = 83 confidence gate (`AGENTS.md` §8a) | `x^83+x^45+x^2+x+1`, the polynomial §8a designates |
| any other Koblitz degree `≥ 64` | the least sparse irreducible by the same rule, searched without the 64-bit cap |

### The endomorphism ring

A Koblitz curve `K_a : y² + xy = x³ + ax² + 1` over `GF(2^n)`, `a ∈ {0,1}`,
certifies `end = -7`: `End(E) ⊇ Z[τ]` with `τ² − t₁τ + 2 = 0`, `t₁ = ±1`,
whose discriminant `t₁² − 8 = −7` is fundamental, so `Z[τ]` is already the
maximal order of `Q(√−7)`.  Every other curve records `unk` until a proof
is attached.

## Names you may write

- **The slug**, anywhere: prose, tables, report fields, file names.
- **The full ICV1 string** in machine records (`curve_id.icv1`).
- **A standard name** for a curve a standards body or public challenge
  published (`secp256k1`, `P-256`, `sect163k1`, `ECC2K-130`, `ECC2K-95`,
  …).  Those are globally fixed and the registry maps each to its model.
- **Family notation with a free parameter** (`K_a / GF(2^n)`, `E(F_{p³})`)
  when a sentence is about a family, not a curve.

## Names that are retired

These denote one curve each but carry no model, so they may not be
written in new text or emitted by new code.  Frozen reports keep them, and
the registry resolves every occurrence:

| retired form | example | what it lacked |
|:--|:--|:--|
| `K_a / GF(2^n)`, `K_a/2^n`, `K_a/F_2^n`, `K₀/GF(2^n)` with a number for `n` | `K_0 / GF(2^41)` | the modulus; five spellings of one curve |
| `k<a>n<n>` as a curve name | `k0n41` | the same; still valid as the stem of a frozen parameter file |
| `bench-<bits>bit` | `bench-20bit` | anything but the roster slot |
| `generated-<bits>bit-<p>` | `generated-24bit-10935329` | `a` and `b` |
| `random-binary-n<n>-b<b>` | `random-binary-n27-b845462` | `a` and the modulus |
| `E_{a,b}/GF(2^k) over GF(2^n)` | `E_{0,2}/GF(4) over GF(2^14)` | the subfield basis the indices refer to |

## Where it lives

| file | what |
|:--|:--|
| [`scripts/curve_id.py`](../../scripts/curve_id.py) | the reference implementation: `python3 scripts/curve_id.py koblitz 0 41`, `… resolve 'K_0 / GF(2^41)'`, `… selftest` |
| [`src/cryptanalysis/curve_id.rs`](../../src/cryptanalysis/curve_id.rs) | the Rust port; `tests/curve_id.rs` pins it to the reference's vectors |
| [`registry.json`](registry.json) | every curve the repository names, with each spelling that denotes it; built by [`scripts/build_curve_registry.py`](../../scripts/build_curve_registry.py) |
| [`sources/generated.json`](sources/generated.json) | curves the repository generated from a seed and never recorded, rebuilt by `examples/curve_id_generated.rs` |
| [`scripts/check_curve_names.py`](../../scripts/check_curve_names.py) | the lint CI runs: no retired form in Markdown or HTML, none added to code, every slug registered |

`KoblitzCurve::label`, the `ic` instances (`PrimeInstance`,
`BinaryInstance`) and `ic`'s `curve_label` emit the slug.  Code that
replays a frozen report, or looks an instance up in a table keyed the old
way (`docs/ic/calibration.json`), matches through
`curve_id::same_curve`, which accepts either name.

## Adding a curve

1. Generate it with code that emits its slug, or compute the slug with
   `scripts/curve_id.py`.
2. If the curve came from a seed and its coefficients are not in any
   tracked JSON, add its specification to `docs/curves/sources/generated.json`
   by rerunning the example (see [`README.md`](README.md)).
3. Run `python3 scripts/build_curve_registry.py` and commit the registry
   with the work that first names the curve.

## Additional prime model forms

The native standards importer preserves these original equations in model JSON:

| `form` | Additional keys (decimal strings) |
| --- | --- |
| `B*y^2=x^3+A*x^2+x` | `A`, `B` |
| `a*x^2+y^2=1+d*x^2*y^2` | `a`, `d` |
| `x^2+y^2=c^2*(1+d*x^2*y^2)` | `c`, `d` |

Each also has `v: "1"`, `p` and `field`, as for short Weierstrass models.
Coefficients are canonical residues. The model hash names the original equation,
not its normalized short Weierstrass equation. Trace and order refer to the
smooth projective model; j is computed using the verified birational change in
[COVERS.md](COVERS.md). Existing identities are unchanged. Native import is in
`src/bin/curve_standards`; [the standards inventory](standards/README.md) explains
provenance and [the YAML graph](cover-links.yaml) links ICV1, EC1 and covers.
