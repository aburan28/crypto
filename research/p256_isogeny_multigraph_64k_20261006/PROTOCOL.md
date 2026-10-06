# P-256 multi-prime isogeny walk: 65,536-curve expansion

Status: **preregistered; execution not started**

Date frozen: 2026-10-06

Result class: **bounded graph-enumeration and screening diagnostic, not an
ECDLP or index-calculus speedup**

## Question

Does the deterministic native multi-prime walker extend the existing verified
20,000-curve P-256 prefix to 65,536 distinct curve identities without an edge,
order, identity, or graph-structure discrepancy, and does that larger prefix
contain a preregistered model-level screening outlier worth a separate
factor-base experiment?

This is not an exhaustive walk of the P-256 isogeny class.  It is one bounded
breadth-first prefix.  In particular, 65,536 curves cover only `2^-16` of the
previously requested `2^32` screening population and a vastly smaller fraction
of the full ordinary isogeny class.

## Frozen implementation and inputs

- implementation tree: `f63704c65689eb910dbfdc21e24cbcf3fad5adc6`;
  the execution commit may add this protocol and result records, but must not
  change `src/`, `Cargo.toml`, or `Cargo.lock` before the run;
- root: registered P-256,
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- requested degrees: every odd prime through 61;
- walkable degrees determined from the frozen class calculation:
  `3,5,11,13,17,23,29,37,41,43,47,59`;
- skipped Atkin degrees: `7,19,31,53,61`;
- `max_curves = 65,536`, no depth cap, root-finding seed 1;
- one point per construction-time order audit and two points per independent
  replay audit;
- canonical model rule `icwalk-canon/v1` and generator rule
  `icwalk-gen/v1`;
- class-audit sample: eight non-root curves;
- CPU thread count is an execution detail and does not change the output;
- expected deterministic run id: `p256-1688961d5b12d121`.

The generation command is:

```bash
target/release/isogeny_walk walk \
  --curve p256 --max-ell 61 --max-curves 65536 \
  --seed 1 --audit-points 1 --class-audit-sample 8 \
  --out "$RUN"
```

The independent replay command is:

```bash
target/release/isogeny_walk verify \
  --curve p256 --audit-points 2 --dir "$RUN"
```

The run is local and network-free.  Publishing is a separate transport step
and must use a configured environment identity; chat-provided access keys are
not inputs to this protocol.

## Reference and boundary

The fixed reference is the earlier P-256 run
`p256-aed06c3fc90f7645`, whose committed summary records 20,000 distinct
curves, 78,672 verified directed edges, zero construction failures, 20,000
prime-order proofs, and a `qr_prefix_64` range of 16 through 48.

The unit is a unique P-256-isogenous curve identity in a verified breadth-first
prefix.  Coverage is reported against the requested `2^32` population:

| variant | unique curves | fraction of `2^32` | ratio to `2^32` boundary | ECDLP speedup |
|:--|--:|--:|--:|:--|
| existing reference | 20,000 | 0.000465661% | 1 / 214,748.365 | unset |
| frozen expansion cap | 65,536 | 0.001525879% | 1 / 65,536 | unset |

The expansion is 3.2768 times the reference population.  Edge count is
reported separately and never substituted for unique-curve coverage.

## Preregistered checks and screens

The run passes the enumeration gate only if all of these hold:

1. exactly 65,536 curve nodes are emitted, with unique `j` invariants, ICV1
   slugs, EC1 identities, and curve UIDs;
2. construction records zero failed modular-polynomial roots and zero refuted
   order audits;
3. all 65,536 order audits have status `proved_prime`;
4. independent `verify` accepts every node, explicit kernel polynomial, Velu
   codomain, target isomorphism, edge id, and ordered `IW1` route;
5. every expanded node has one rational neighbour for degrees 3 and 5 and two
   for each walkable degree from 11 through 59 listed above;
6. all detector records say `non_singular = true` and
   `generator_valid = true`; and
7. the sampled class audits are byte-equivalent in verdict to the P-256 root.

The following are frozen *screening* thresholds.  Crossing one selects a curve
for a new preregistered experiment; it does not establish an attack gain.

- `qr_prefix_64 >= 52`.  Under the deliberately simple independent
  `Binomial(64, 1/2)` reference, the one-curve upper-tail probability is about
  `2.2833e-7` and the 65,536-curve expected exceedance count is about
  `0.01496`.  This statistic is only a cheap local proxy for usable
  factor-base yield.
- signed bit length of canonical `b <= 224`, or signed bit length of canonical
  `a <= 224` when the curve has no `a = -3` model.  This flags unusually small
  coefficients for arithmetic follow-up, not weaker ECDLP structure.

No other threshold may be introduced after reading the outputs and described
as preregistered.  Exploratory observations must be labelled exploratory.

## Evidence and acceptance

Preserve without overwriting:

- `curves.yaml`, `isogeny_routes.json`, and `walk.json` from construction;
- `VERIFY.json` from independent replay;
- GNU `/usr/bin/time -v` logs for construction and replay, or, when that
  executable is absent, a Python standard-library `subprocess`/`resource`
  wrapper recording exit status, wall time, user time, system time, and
  `ru_maxrss` without inspecting or changing the Rust workload;
- SHA-256 and byte length for every file;
- the exact source commit, Rust version, host CPU count, and peak RSS; and
- a compact result note containing the boundary table, graph counts, detector
  distributions, any threshold-crossing curve identities, and the decision.

The large route and curve files do not belong in Git.  Before a positive result
is treated as durable evidence, publish their deterministic gzip encodings to
write-once, hash-checked storage and record the object URIs and both compressed
and uncompressed digests in `MANIFEST.json`.  If no configured storage identity
is available, record that as a publication blocker; do not substitute a local
absolute path.

## Cost and stop conditions

Budget one construction and one full replay.  Expected construction resources,
extrapolated from the committed 20,000-curve run, are roughly 5.7 GiB peak RAM,
1.1--1.4 GiB raw output, and under 30 minutes on a many-core host.  Stop and
record the partial failure if construction or replay exceeds two hours, the
filesystem has less than 1.5 GiB free before output serialization, the process
is killed for memory, any acceptance check fails, or the source tree differs
from the frozen implementation plus research-only records.

Do not tune degrees, seed, cap, thresholds, or canonicalization and silently
rerun.  Any changed experiment gets a new protocol.  Even if a screen crosses
its threshold, an index-calculus claim remains unset until a separate frozen
factor-base test and then a cold, end-to-end, matched rho comparison include
relation generation, verification, linear algebra, and target-log recovery.
