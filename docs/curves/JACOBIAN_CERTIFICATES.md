# CM Jacobian existence certificates

`jacobian_certificate` checks a supplied positive unimodular Hermitian
matrix over an imaginary quadratic maximal order. It exports exact lattice
findings and the **conditional** consequence that a genus-two Jacobian has
underlying abelian variety E squared. This is a sufficient test, not a
complete decision algorithm for every matrix, order or elliptic curve.

```sh
cargo test --bin jacobian_certificate
cargo run --bin jacobian_certificate -- --sql docs/curves/jacobian-certificates.sql
cargo run --bin jacobian_certificate -- --check --sql docs/curves/jacobian-certificates.sql
sqlite3 curves.sqlite < docs/curves/jacobian-certificates.sql
sqlite3 curves.sqlite "SELECT discriminant,status,json_extract(record_json,'$.lattice.minimum') FROM jacobian_certificates;"
```

The committed input is [jacobian-candidates.json](jacobian-candidates.json).
The outputs are [JSON](jacobian-certificates.json) and an idempotent
[SQLite import](jacobian-certificates.sql). JSON is the portable interchange
format; SQLite is a rebuildable lookup index. No running database is contacted.
The [JSON Schema](jacobian-certificates.schema.json) describes the record boundary.
`--input`, `--output`, `--registry`, `--sql` and `--budget` customize the run.
`--check` independently recomputes the mathematics and compares exact output
bytes, including source and input digests. A hash alone is not a proof.

## Algorithm and mathematical contract

Write D for a negative fundamental discriminant, t=D mod 2, and
omega=(t+sqrt(D))/2. The input matrix is
`[[a, b0+b1*omega], [conjugate(b0+b1*omega), d]]`.
Every calculation uses bounded exact i128 integers; no floating point or
probable-prime tests occur. Inputs are checked before arithmetic:
`-1000000 <= D <= -3`, `1 <= a,d <= 1000`, `|b_i| <= 1000`.
These bounds put even the largest norm intermediate far below i128 limits.
Larger inputs receive `unsupported`, rather than being truncated.

1. Check that D is fundamental by squarefreeness and congruence, and that
   `a*d-N(b)=1`. Together with a>0 this proves positive definiteness and
   unimodularity.
2. Enumerate all reduced primitive positive binary quadratic forms
   `(A,B,C)` of discriminant D: `|B| <= A <= C`, boundary B>=0, and
   `A <= floor(sqrt(|D|/3))`. Set M to the largest leading coefficient.
3. Determine the exact Hermitian minimum using
   `a*h(x,y)=N(a*x+b*y)+N(y)`. Initially `h_min <= min(a,d)`.
   Enumerate both norm balls with bound `a*min(a,d)`, retaining
   `z=b*y mod a`. Every vector below the initial bound occurs this way;
   the coordinate-vector witness supplies the boundary value. The coordinate
   identity `4*N(u+v*omega)=(2*u+t*v)^2+|D|*v^2` gives exact finite bounds.
4. A decomposable unimodular Hermitian module is an orthogonal sum of
   rank-one ideal modules in inverse ideal classes. The normalized norm
   minimum in each class is its reduced form's A. Thus its minimum is at
   most M. A minimum greater than M certifies indecomposability. Minimum
   one certifies a free orthogonal splitting. Other cases remain inconclusive.

Norm-ball coordinate trials and vector-pair trials consume the reported
`work_units`; the default budget is 2,000,000 per candidate. Reduced-form
enumeration is separately bounded by the discriminant limit. Exhaustion
does not return a partial minimum as exact. Up to 100 candidates are accepted.

For D=-619, the forms are `(1,1,155)`, `(5,+/-1,31)`, `(7,+/-5,23)`.
The proposed matrix has determinant one, minimum 12 and M=7. In particular,
this checks nonprincipal ideal classes; minimum greater than one alone
would not be a sufficient indecomposability test.

## What remains conditional

The program does not identify End(E), prove ordinarity, or certify the
field of definition of endomorphisms. Required geometric hypotheses are
always recorded, and `exists_over_base_field` remains null. Evidence links
are unverified references and cannot promote the result. No caller-supplied
boolean can turn a lattice check into an unconditional curve theorem.

After an independent audit establishes that an ordinary E over the declared
finite field has this geometric maximal order acting over that field, the
Hermitian correspondence gives a principal polarization on E squared.
Geometric indecomposability then gives a genus-two Jacobian over that field
by Weil's theorem (Howe--Nart--Ritzenthaler, Theorem 1.3).

`decomposable_lattice` or a rejected candidate says nothing about whether
another polarization or another cover exists. Curve equations, cover maps,
subgroup transfers and computational advantage remain unconstructed/unknown.

## Identity and database integration

The proposed CryptoPro-B lattice is deliberately **unbound**: the registry
snapshot used for this implementation contains no CryptoPro-B record. No
EC1 alias, subgroup metadata or curve UID is fabricated from the name.

To associate a candidate with an existing catalog representation, set
`binding` to `{"slug":"<existing ICV1 slug>","curve_uid":"<full existing UID>"}`
and supply `--registry docs/curves/registry.json`. Exactly one matching
representation is required. The output retains its EC1 alias, field, curve
and model metadata plus the registry digest. This is a metadata association;
it does not authenticate the registry or certify its mathematics, and it
does not establish that the proposed order actually acts on that E.

The input uses `jacobian-candidates/v1`; records use
`jacobian-certificate/v1`; the envelope uses `jacobian-certificates/v1`.
`certificate_uid` hashes the complete record before that UID is inserted,
using BLAKE3 over sorted-key compact serde_json bytes. It includes candidate,
binding, result, assumptions and provenance. This is separate from the
existing SHA-256 curve identity and never modifies its preimage.

SQLite indexes `(curve_uid,status)` and `(discriminant,status)`. Join on the
full curve UID; unbound records have SQL NULL and remain visible by order
or certificate UID. Never join on the display label or discriminant alone.
Useful queries include:

```sql
SELECT certificate_uid, status, record_json
FROM jacobian_certificates WHERE curve_uid = :curve_uid;

SELECT certificate_uid, json_extract(record_json,'$.lattice.minimum') AS minimum
FROM jacobian_certificates
WHERE discriminant = -619 AND status = 'conditional_existence';

SELECT certificate_uid FROM jacobian_certificates WHERE curve_uid IS NULL;
```

Imports use a transaction and content-addressed insert-on-conflict-do-nothing.
Distinct evidence versions remain separate records. SQL literals are escaped.
These files provide the database ingestion boundary; no live database or UI
deployment is claimed by generating them.

## Validation and references

Native tests include the -619 certificate, an independent direct evaluation
of small Hermitian forms, even/odd fundamental discriminants, the -20
nonprincipal-class case, determinant mutation, extreme integer rejection,
budget exhaustion, strict inputs, exact catalog binding and deterministic
record identities. CI replays JSON/SQL and exercises a real SQLite import.

- Gélin, Howe and Ritzenthaler, *Principally polarized squares of elliptic
  curves with field of moduli equal to Q*, Propositions 3.1 and 3.3:
  https://arxiv.org/abs/1806.03826
- Ritzenthaler, *Optimal curves of genus 1, 2 and 3*, Section 3.1, for the
  ordinary finite-field Hermitian-module interpretation:
  https://pmb.centre-mersenne.org/item/10.5802/pmb.a-137.pdf
- Howe, Nart and Ritzenthaler, *Jacobians in isogeny classes of abelian
  surfaces over finite fields*, Theorem 1.3:
  https://www.numdam.org/item/10.5802/aif.2430.pdf
