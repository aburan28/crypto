# IC curve records and ICV1 crosswalk

[`curves.yaml`](curves.yaml) and its [schema](curves.schema.json) mirror the exact small curve registry in
`aburan28/cryptanalysis` at
`experiments/ic-candidate-catalog/curves.yaml`. The mirrors have the same
bytes and retain full `EC1` curve UIDs, exact field/curve records,
known/unknown trait statuses and isogeny route references. Evidence paths
inside it are relative to **cryptanalysis**. Change both copies in one
paired update; compare their SHA-256 before merging either update.

`ICV1` is crypto's model/display name. `EC1` plus the full
`urn:ec-record:1:sha256:...` UID identifies the field representation,
model, subgroup and generator used by a DLP experiment. The toy N13 and
ECC2K-130 source entries match existing ICV1 **models**, but their EC1
representations are not currently entries in
[`registry.json`](../registry.json). The degree-263 codomain has an exact
EC1 record and verified route in cryptanalysis; its general binary
Weierstrass model is not yet registered under ICV1. Its ICV1 fields are
therefore `null`. Do not synthesize an ICV1 slug from a translated or
isomorphic model without a recorded isomorphism and registry support.

The [factor-base browser](https://aburan28.github.io/crypto/browser/#factor-bases)
is generated from committed ecbench sessions and lists `FB1` summaries.
It is not a point store. Large bases remain in cryptanalysis's
[`fb-archive`](https://github.com/aburan28/cryptanalysis/tree/main/experiments/fb-archive):
compressed point sets, orbit representatives, or content-addressed shards
with counts and hashes. Only attach an `FB1` to an archive entry after
the curve UID, point-set encoding/digest, quotient rule and actual usable
point count have been reconciled. `fb<B>` cannot be formed from a recipe
with unknown `B`.

[Typed curve-link rules](curve-links/README.md) distinguish twists, same-field
model changes, and field extensions from isogenies. The mirrored curve YAML
keeps explicit `not_enumerated` inventories for links that have not been
constructed or checked. A shared j-invariant is not a base-field point map.

The directed isogeny graph and ordered `IW1` routes live in
cryptanalysis's
[`isogeny_routes.json`](https://github.com/aburan28/cryptanalysis/blob/main/experiments/ic-candidate-catalog/isogeny_routes.json).
An isogenous neighbor gets its own curve UID. Its conductor, prime-specific
volcano level, map, kernel and target-log transport proof are separate
evidence. An unknown value is YAML `null` with a status, never `0`
or a question mark. A complete IC claim still needs an exact candidate,
frozen workload, verified one-target run, and paired rho run under the
existing measurement contract.
