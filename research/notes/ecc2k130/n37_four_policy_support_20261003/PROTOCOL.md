# Frozen n37 four-policy useful-support gate

Status: **preregistered, unmeasured**. This gate is a prerequisite for a
natural-target PDP and complete-cost comparison, not a speed claim. The
registered source curve is `icv1-f2m37-tm534059-32aad96b`, with prime
subgroup order `r = 230603167`, cofactor 596 and degree-73 descending map to
the archived leaf. The matched attack reference in any later cost panel is
cold signed-Frobenius rho on the same public source targets.

## Question and frozen inputs

Can the original, transported, descendant-native and pullback policies each
provide **42 independent log columns and 1,554 distinct signed support
classes** at n37, while preserving a checked relation label under the
degree-73 map? The factor `37` is the source Frobenius orbit length in the
prime subgroup. A class contains a point and its negative, so the target is
3,108 physical points per policy. Equal column count alone is insufficient:
the earlier native K42 base without closure has only 42 signed classes.

- The original source base is the exact `b00_ic_a.base.jsonl` member of
  `research/notes/ecc2k130/disjoint_cold_v2_outcome_20261001/evidence_run_36803331080/raw/n37_L1024.tar.gz`.
  Its member SHA-256 is
  `0a32de24a5680ff46baf9543e8bf8e32447323491e2a2adedce2548c14a25f75`
  and its internal `base_hash` is
  `8423b135df3515b0e284126d9eeb95f71e3901ae12908b3cad8eb03d1f31d4bc`.
  It records 42 representatives, 3,108 full points and their orbit labels.
  The source selection algorithm is the ascending abscissa scan in
  `examples/koblitz_orbit_dlp_s3_batch.rs`; the archived member fixes the
  input for this gate. Later cold runs must repeat and charge that scan.
- The descendant-native seeds are exactly the 42 `leaf` rows of
  `research/notes/ecc2k130/n37_native_basis_bridge_20261002/NATIVE42.json`,
  SHA-256 `bd4bd8af982bcc65234ae1fae1ad52eb0ed41a1c203bf4d2556fd8e87f788a4c`.
  They were selected by ascending leaf abscissa before the held-out targets.
  Do not replace a seed after seeing orbit collisions or target outcomes.
- The archived degree-73 map is
  `research/koblitz_isogeny_descent_37_results_20260925/raw.json.gz`,
  SHA-256 `eb4773d556886672e8b80735fc486c16b73765fd5069ff81133100bfc3eafa90`.
  The source-to-archive field-basis generator image is `10156182909`.
  Rebuild the map from this input; reject any model or order mismatch.
- Public held-out target inputs for the subsequent PDP gate are the n37/L1024
  block-03 and block-04 point-only files from
  `research/notes/ecc2k130/disjoint_cold_v2_20261001/fixtures/`, SHA-256
  `84400a2914f06e4a001d0f113f0195952f692634d4ff8285e23ade2599e2bde2`
  and `72c5b5361a0b846c9626269ce3928820cb7d420aae479d31dc26f5cbf7e78aa6`.
  No target file or fixture scalar is read during this support gate.

## Policies and exact admission rule

Let `lambda` be the source Frobenius eigenvalue modulo `r`, checked by
`[lambda]G = Frobenius(G)`. For each seed, use exponents `k = 0..36` and
both signs, with coefficient `±lambda^k mod r` on that seed's log column.
The policies are fixed before measurement:

| Policy | Curve | Seed origin | How 37 signed classes per column are obtained |
| --- | --- | --- | --- |
| Original | source | archived compact source base | cheap coordinate Frobenius |
| Transported | leaf | image of each source seed | map each source orbit member through the degree-73 isogeny; negation supplies its pair |
| Descendant-native | leaf | frozen `NATIVE42.json` leaf seeds | scalar multiplication by `lambda^k`; every evaluation is a costed operation in later cold runs |
| Pullback | source | inverse image of each frozen leaf seed | exact rational pullback once per seed, then source Frobenius |

The four sets need not contain the same points. Within each paired policy,
however, the map must preserve every **full point**, sign and label:
original ↔ transported and descendant-native ↔ pullback. The native seeds
must occupy 42 distinct `[lambda]`-and-negation orbits; if they do not, the
fixed K42 gate fails and this run may not silently substitute new seeds.
Check all 3,108 points per arm for curve and order membership, all 1,554
signed classes for distinctness, all 42 labels on the seed representatives,
and the 37th action returning to each seed. A preimage is accepted only when
its complete forward image equals the given leaf seed. Infinity and the
2-torsion point are separate negative/exception controls and never factors.

Record the action calls, map calls, inverse calls, source and leaf group
operations, field operations where instrumented, retained bytes and elapsed
phases. A count left unpriced must be marked as such. This support gate
does **not** divide those counts by rho or report `S`: archived source
selection is an input here, and no PDP, rank or held-out target is run.
The full cold experiment must rebuild both selections and map in fresh
processes and include all those costs.

## Audit, decision and follow-on

The producer writes an immutable manifest with source/input hashes,
representative and point coordinates, labels, orbit-class counts and every
failure. A separate Rust replay recomputes the field bridge, map, inverse,
group action and set identities from the frozen inputs; it must not trust
the producer's point or label arrays. Commit the manifest, replay receipt,
source and input hashes, commands, host facts and decision in this PR. If a
limit or mismatch stops a policy, preserve the partial result and mark the
gate failed, without reporting the missing policy as equal-sized.

On a pass, the next PR must run a frozen natural-target three-summand PDP
gate on the two held-out blocks, independently check every witness and
rank transition, and report proved misses separately from timeouts or
unknown. Any solver choice or factor-base change needs a new protocol.
Only after a complete cold, one-target and batch comparison in a calibrated
common unit may a policy be compared to matched rho or the generic floor.
No n41/n53, n83, degree-263, n131, or ECC2K-130 performance conclusion
follows from this support gate.
