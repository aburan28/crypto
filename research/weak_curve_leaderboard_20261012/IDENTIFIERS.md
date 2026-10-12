# WID1: an identifier for a curve weakness finding

The repository's [ICV1](../../docs/curves/ICV1.md) identifies an exact field
representation and curve model. WID1 identifies **one mechanism claim on one
ICV1 model or selected subgroup**. A curve may have several WID1 findings.
A named standard is an alias pointing to ICV1, never a second curve identity.
This is a catalog convention, not a security-bit estimate or a theorem that
the associated attack is practical.

The printable form is `WID1:<MECHANISM>:<12 lowercase SHA-256 hex digits>`.
The digest is SHA-256 of compact UTF-8 JSON with sorted keys and no whitespace:

```json
{"criterion":"ghs_binary","mechanism":"GHS","scope":"selected_subgroup","selector":"r=3805993847215893016155463826195386266397436443;tower=5x31;g=15","subject_id":"ICV1:f2m-155-dea39df2:-414891960790832521345347:45671926166590716193865565914344635196769237316:0x61a2976b585b70cc897ad03a7f74a4ecc23bf4a:unk:unk:r:b865721c660a","v":1}
```

The corresponding identifier is `WID1:GHS:d089b2bdf750`. Oakley Group 4's
distinct ICV1 model and subgroup have `WID1:GHS:419192b724ad`. If a 12-digit
prefix collision is ever found, extend the digest for colliding records and
retain aliases for the old display values.

The identity key contains the mechanism, criterion, exact subject ID, scope,
and selector. `scope=model` means a claim about a model or covering geometry;
`scope=selected_subgroup` fixes an actual group order or subgroup route in the
selector. A full curve order and a generator order are separate quantities.
The identity excludes evidence grade, source links, costs, timestamps and
free-form prose. Those can change as a finding is checked without renaming it.
Changing the curve, selected order or route produces another WID1.

| Code | Mechanism in this snapshot | Current claim scope |
| --- | --- | --- |
| `ANOM` | Trace-one anomalous prime-field condition | Selected subgroup |
| `PH` | Smooth selected cyclic order and Pohlig–Hellman | Selected subgroup |
| `MOV` | Pairing embedding-degree candidate | Selected subgroup |
| `GHS3C` | Cubic norm-one genus-three model label | Exact model |
| `GHS3Q` | Quadratic-torus genus-three model label | Exact model |
| `GHS` | Binary GHS tower invariant and cover candidate | Selected subgroup |

Each finding has a separate `grade`: `exact_condition`,
`covering_model`, `candidate`, or `exact_invariant_candidate`. The two Oakley
GHS findings have the last grade. For each, the recorded magic number five is
an exact computation on the model; genus 15 follows a GHS formula, but the
cover and selected-subgroup transport are not constructed in this snapshot.
The repository transfer file includes a formula cost projection whose chosen
solver's hypothesis is explicitly false. That projection is stored as
`formula_only` and is never promoted to a measured advantage.

Oakley Group 3 and Group 4 are the IKE EC2N groups over
`GF(2^155)` and `GF(2^185)`, respectively, in
[RFC 2409 §§6.3–6.4](https://www.rfc-editor.org/rfc/rfc2409.html).
[RFC 2412 Appendix E](https://www.rfc-editor.org/rfc/rfc2412.html) specifies
their full group orders as `12r` and `4r`; its published generators have
orders `4r` and `2r`. The factors `r` are large primes, so composite full
order alone does not mark the selected prime-order DLP as smooth. The checked
Sage replay in `validation/oakley.sage` verifies the defining polynomials,
RFC points and the exact orders of the derived `[4]P` and `[2]P` subgroup
generators, using the RFC group-order claim as source data.

This snapshot imports 103 named curves from the repository registry, including
101 curves with standards provenance and two named challenge curves. It does
not enumerate every published elliptic curve or every possible mechanism.
Rows without a finding carry no WID1 until a precise criterion is evaluated;
absence of a WID1 is not a classification as secure.
