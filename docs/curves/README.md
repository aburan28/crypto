# Curves

The [weak-curve family catalogue](weak-families/README.md) organizes 21
published or verified families and structural conditions, with their recognition
rules, construction domains and evidence. Its [searchable browser](weak-families/index.html)
indexes 814 exact models and the complete 294-class ordinary p7 oracle.
Model membership, class support, subgroup compatibility and measured cost
have separate labels; [canonical JSON](weak-families/catalog.json) retains
the source attribution and completeness boundary for each family.

[CM Jacobian certificates](JACOBIAN_CERTIFICATES.md) check principal
polarizations on elliptic squares using exact Hermitian arithmetic and all
ideal classes. The [versioned JSON](jacobian-certificates.json) and
[SQLite import](jacobian-certificates.sql) support lookup by full curve UID
when bound, and by order or certificate UID when unbound. Geometric
endomorphism hypotheses remain explicitly conditional.

[Automatic hyperelliptic cover checks](COVERS.md) attach replayable
same-field cover certificates to each supported catalog model in
[`covers.json`](covers.json). Each curve page in the lab browser displays
the verified genus, map degree and field assumptions. Descent and DLP
advantage remain unmeasured.

The [hyperelliptic infrastructure guide](HYPERELLIPTIC_INFRASTRUCTURE.md)
documents the standalone construction, verification, arithmetic, and catalog
commands; the one-infinity arithmetic boundary; bounded binary transfer; and
the additive evidence schema. In particular, the prime even-sextic cover does
not use the odd-degree one-infinity Jacobian implementation.

[IC curve records and cross-repo links](ic/README.md), including
[typed links](ic/curve-links/README.md), retain the exact EC1 representations,
optional trait statuses, and the factor-base/isogeny storage contract beside
this ICV1 model registry.

Curves isogenous to a registered prime-field curve are found and recorded
by the native isogeny walker, `src/bin/isogeny_walk.rs`
([`ic/README.md`](ic/README.md#walking-an-isogeny-class)): it writes each
curve in the `ic/curves.yaml` format with its ICV1 slug, EC1 identity and
traits, and each kernel-certified edge as an `IW1` route.

Every curve this repository names is named by its **ICV1 slug**
([`ICV1.md`](ICV1.md), `AGENTS.md` §11).  Each curve in the registry also
carries the **EC1 identity** of each exact representation the repository
records (subgroup and generator included), for comparisons joined across
repositories ([`../curve-identities.md`](../curve-identities.md)). Those records use the
encoding of crypto's existing producer, `examples/koblitz_curve_records.rs`, so a UID in the
registry equals the one a thread's `curve_ids.json` carries; the builder fails if they ever differ.

| file | what |
|:--|:--|
| [`ICV1.md`](ICV1.md) | the naming rule: the identity, the model JSON, the fields, the retired forms |
| [`registry.json`](registry.json) | every curve named anywhere in the repository: slug, ICV1, parameters, every legacy spelling that denotes it, and its EC1 representations |
| [`sources/specs.txt`](sources/specs.txt) | the constructor calls whose curves no tracked record states in full |
| [`sources/generated.json`](sources/generated.json) | those curves rebuilt: model, subgroup, generator, legacy handle, ICV1 computed in Rust |
| [`TRAITS.md`](TRAITS.md) | size-independent traits of every registered curve (CM field, conductor, subfield of definition, volcano depths, cofactors, embedding degree) and how to group curves or find similar ones with `curve_traits` |
| [`traits.json`](traits.json) | those traits, one record per registry curve, each value with its status; built by `cargo run --release --bin curve_traits -- build` |

## Two identities, and which to use

| | ICV1 slug | EC1 alias and curve UID |
|:--|:--|:--|
| identifies | the curve model: field, modulus, coefficients | one representation: the model plus its subgroup and generator |
| written in | prose, tables, report `name` fields, file names | cross-repository comparisons, UI exports, candidate manifests |
| example | `icv1-f2m41-tm2308219-7f48b14a` | `EC1N41Ck0h…` with `urn:ec-record:1:sha256:…` |
| one per | model | representation; a model can have several |

A slug in text says which curve; an EC1 UID in a record says which
representation was measured.  Join on the EC1 UID, never on the slug's
digest prefix or on a degree.

## Regenerating

```bash
cargo build --release --example curve_id_generated
./target/release/examples/curve_id_generated $(grep -v '^#' docs/curves/sources/specs.txt) \
    > docs/curves/sources/generated.json
python3 scripts/build_curve_registry.py
python3 scripts/check_curve_names.py
```

The builder reads every tracked JSON file, so stage a new source file
before rebuilding.  It fails when a name used in Markdown or HTML does not
resolve to exactly one curve, when a Rust-computed ICV1 differs from the
reference, or when a pin in `docs/ic/calibration.json` no longer resolves.

## Resolving a name

```bash
python3 scripts/curve_id.py resolve 'K_0 / GF(2^41)'
python3 scripts/curve_id.py koblitz 0 41
```

A legacy spelling resolves through one of two routes, and a curve's
`sources` in the registry say which:

- **by record** — a tracked report stores the curve's coefficients next to
  the name it used, so the name denotes that model (`sources` lists the
  report);
- **by the constructor rule** — `K_a / GF(2^n)` written in prose denotes
  the curve `KoblitzCurve::new(a, n)` builds, on the modulus
  [`ICV1.md`](ICV1.md) states (`sources` reads `prose, by the constructor
  rule: <file>`).  This is a versioned adapter, not an inference from the
  degree: a note that meant another modulus (another repository's `K_0`
  over `GF(2^13)` is on `x^13+x^5+x^2+x+1`, this one's on
  `x^13+x^4+x^3+x+1`) is a different model, and must say so.
