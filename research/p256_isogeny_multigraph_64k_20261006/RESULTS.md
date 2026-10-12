# P-256 multi-prime isogeny walk: 65,536-curve results

Status: **local execution and hardened replay complete; large-artifact
publication blocked**

Date: 2026-10-06

Result class: **bounded graph-enumeration and screening diagnostic, not an
ECDLP or index-calculus speedup**

## Outcome

The preregistered 65,536-curve gate passed locally.  Starting from registered
P-256, the native walker produced 65,536 distinct curve identities and 302,653
explicitly certified directed isogeny edges across the 12 walkable odd degrees
through 61.  Construction recorded zero failures and proved the P-256 prime
group order on every curve.  A second, hardened verifier replayed every node,
kernel, Velu codomain, target isomorphism, edge id and route and returned
`pass`.

No preregistered screen fired.  The largest `qr_prefix_64` value was 48, below
52; the smallest canonical `b` was 236 signed bits, above 224; and the smallest
canonical `a` among curves without an `a = -3` model was 240 signed bits, above
224.  Therefore this run selects no curve for a factor-base experiment and
sets index-calculus speedup to **unset**.

The full raw evidence is not yet durable.  The managed environment reported no
configured cloud identity, so the 280,242,575 bytes of deterministic gzip
artifacts were not published to S3.  `MANIFEST.json` records their raw and
compressed hashes and the publication blocker.  Chat-provided access keys were
not used.

## Boundary table

The unit is a unique curve identity in a verified breadth-first prefix.  Edge
count is shown separately and is not substituted for curve coverage.

| variant | unique curves | verified edges | fraction of `2^32` | ratio to `2^32` boundary | `qr_prefix_64` max | ECDLP speedup |
|:--|--:|--:|--:|--:|--:|:--|
| existing reference | 20,000 | 78,672 | 0.000465661% | 1 / 214,748.365 | 48 | unset |
| 2026-10-06 expansion | 65,536 | 302,653 | 0.001525879% | 1 / 65,536 | 48 | unset |

This adds 45,536 curve identities beyond the earlier deterministic cap, a
3.2768-times larger population.  It still covers only `2^-16` of the requested
`2^32` population and is not an exhaustive enumeration of the P-256 isogeny
class.

## Enumeration checks

| check | frozen requirement | observed | result |
|:--|:--|:--|:--|
| curve nodes | exactly 65,536, unique identities | 65,536; hardened uniqueness replay passed | pass |
| construction failures | 0 | 0 | pass |
| order audits | 65,536 `proved_prime`, 0 refuted | 65,536 / 0 | pass |
| certified edges | full replay | 302,653 / 302,653 | pass |
| rational neighbours | 1 for 3 and 5; 2 for each split degree | exact on all 13,757 expanded curves | pass |
| detector validity | all nonsingular, all generators valid | 65,536 / 65,536 for both | pass |
| class-audit sample | identical to P-256 root | 8 / 8 identical | pass |

Every expanded curve had 22 rational outgoing isogenies: one each for degrees
3 and 5 and two each for
`11,13,17,23,29,37,41,43,47,59`.  Thus
`13,757 * 22 = 302,654`; exactly one otherwise-valid edge was omitted because
its new target would have crossed the 65,536-curve cap, matching
`edges_beyond_cap = 1`.

The walk reached depth 6.  Its node counts by depth were 1, 22, 241, 1,760,
9,680, 42,944 and 10,888.  There are 51,779 unexpanded nodes with certified
root routes: 40,891 at depth 5 and 10,888 at depth 6.  They are a concrete
anchor pool for a later distributed frontier protocol, but this run does not
yet Merkle-commit, partition or globally deduplicate work rooted at those
anchors.

## Frozen screens

The complete native trait pass covered all 65,536 curves and was bound to the
route hash in `metrics.json` and `collect.json`.

| screen | threshold | observed | decision |
|:--|:--|:--|:--|
| `qr_prefix_64` | at least 52 | range 16--48; four curves at 48 | no candidate |
| canonical `b` signed bits | at most 224 | minimum 236; one curve | no candidate |
| canonical non-`a=-3` `a` signed bits | at most 224 | minimum 240; two curves | no candidate |

There were 33,156 curves with an `a = -3` model and 32,380 without one.  The
`qr_prefix_64` distribution is in `runs/p256-ell61-64k/collect.json`.  It is a
model-level prefix statistic, not a measured relation yield, and no solver was
run after all three preregistered screens failed to select a finalist.

## Verification discrepancy found and closed

Reviewing the acceptance path found a verifier weakness, not a curve anomaly.
The original replay recomputed each node and edge but inserted their external
identifiers into maps without explicitly rejecting duplicates, and it did not
compare a node's recorded `j_invariant` with the recomputed model value.  The
generator itself deduplicates on `j`, so this did not show that the run
contained duplicates; it showed that an independently supplied record could
exploit an acceptance gap.

The verifier now:

- recomputes and compares every recorded `j_invariant`;
- rejects duplicate curve refs, `j` values, ICV1 slugs, EC1 ids and full UIDs;
- rejects duplicate edge ids and route ids; and
- has tamper tests for a forged `j` and duplicate node, edge and route records.

All five native `cryptanalysis::isogeny_walk::tests` passed, library clippy
passed with `-D warnings`, and the hardened release verifier replayed the
unchanged route SHA-256
`5dfe974e55fb3b83f237b7ba6646966803ab6f4da27352f4295cb0dcb7aea892`.
It accepted 65,536 unique nodes and 302,653 unique edges with two audit points.

## Cost

| stage | wall seconds | user seconds | peak RSS | exit |
|:--|--:|--:|--:|:--|
| construction and serialization | 2,186.917 | 7,837.557 | 4.141 GiB | 0 |
| hardened full replay | 482.756 | 1,860.185 | 4.953 GiB | 0 |
| complete trait pass | 25.659 | 79.678 | 4.028 GiB | 0 |

The first launch attempt stopped before invoking the walker because this host
does not provide `/usr/bin/time`.  Its 78-byte stderr log was preserved; the
protocol was amended before execution to allow the standard-library resource
wrapper.  That wrapper later wrote two literal terminal bytes `\\n` after each
otherwise valid resource JSON.  The originals were preserved and hashed, and
the committed copies replace only those bytes with one newline after a
successful JSON parse.  `MANIFEST.json` records both forms.

## Evidence and durability

Compact evidence committed here includes:

- `walk.json`: configuration, class invariants, counts, distributions and raw
  artifact hashes;
- `VERIFY.json`: the hardened full-replay receipt;
- `metrics.json` and `collect.json`: full trait-pass binding and distributions;
- resource records for generation, replay and trait detection; and
- `MANIFEST.json`: byte counts and SHA-256 values for every large raw and
  deterministic-gzip artifact.

The write-once S3 marker is absent, so this result is labelled locally verified
rather than durably published.  Publication must use a configured environment
identity and must reproduce every manifest hash; adding a URI without that
check does not close the blocker.

## Decision

The bounded enumeration succeeds, but the anomaly screen is negative.  Do not
promote any curve from this run to an index-calculus or whole-solver benchmark,
and do not claim a speedup.  The next useful enumeration step is a separately
preregistered frontier protocol over the 51,779 path-certified unexpanded
anchors, with immutable anchor membership, deterministic ownership, global
curve-UID deduplication and a checker that replays the anchor route plus each
worker edge.  Merely launching more copies of the old degree-11 prefix would
not add curves.
