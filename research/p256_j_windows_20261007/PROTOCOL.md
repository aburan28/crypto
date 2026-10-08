# P-256 global-row j-invariant continuation

Status: **preregistered; execution pending**

Date frozen: 2026-10-07

Result class: **bounded isogeny enumeration and structural screening, not an
ECDLP or index-calculus speedup**

## Requested scope and evidence plan

| requirement | acceptance evidence | status before execution |
|:--|:--|:--|
| continue the P-256 j-invariant search | a disjoint global-coordinate strip after the verified million-grid prefix | pending |
| preserve exact isogeny evidence | one replayed kernel certificate and codomain link into every emitted curve | pending |
| count cumulative coverage exactly | full-width j union audit against the frozen million certificate, with zero duplicates | pending |
| search structural anomalies and discrepancies | unchanged frozen detectors, thresholds, identity gates, and retained failures | pending |
| establish an index-calculus speedup | cold full-pipeline candidate and matched-rho measurements | not part of this enumeration round |
| exhaust `2^32` or the full P-256 isogeny class | exact cumulative count equal to that boundary | open beyond this round |

This round does not redefine an exhaustive class search as a coordinate scan.
Different product-grid coordinates are candidate curve identities until the
exact j-union audit proves them distinct within the executed population.

## Frozen curve and construction

- root: registered P-256,
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- field, group order, trace, canonical-model rule, deterministic-generator
  rule, root seed, modular-polynomial construction, and detector definitions:
  unchanged from `research/p256_isogeny_million_20261006/PROTOCOL.md`;
- vertical action: non-backtracking degree-13 spine;
- horizontal action: non-backtracking degree-11 rows;
- first edge chooses the numerically smaller rational root, later edges choose
  the rational root unequal to the predecessor j-invariant;
- production x coordinates: `0 <= x < 1000`;
- production y coordinates: `1000 <= y < 1064`;
- target population: exactly 64,000 emitted curves;
- construction order audit: one deterministic point beginning at `x = 0`;
- independent replay: two deterministic points beginning at `x = 7`.

Generation recomputes the degree-13 boundary from P-256. The curve at
`(0,999)` is a chain-bound context record and is not counted in the new
population. The first emitted curve `(0,1000)` carries its degree-13 kernel
certificate from that context. Every other emitted spine or row curve carries
its ordinary parent edge. No prior curve is counted again.

The existing million certificate is immutable and identified by:

- gzip SHA-256
  `87522b6de08e64804dcee17bba058eb91922d6b8dd9853cd32523f67a3d0cbdc`;
- exact decompression SHA-256
  `2a1a27d06f48b22b85ce3c96ac4017f5af6a0256e9d1ad277e3fb1f516475454`;
- record-chain SHA-256
  `4de6e93e0d7a90ff53a26f8f1fd1cbd03ab55f701a25f8e49874cb8bead40302`;
- 1,000,000 unique j-invariants already independently replayed.

## Certificate and replay gates

The strip uses a new schema. It must not change the legacy million schema or
the bytes produced by its existing `generate` command. The strip passes only
if all of the following hold:

1. the header freezes width 1,000, row start 1,000, height 64, global
   coordinates, source commit, and construction rules;
2. exactly 64,000 emitted curves and 64,000 incoming kernel certificates are
   generated and replayed;
3. the boundary context is independently reconstructed from P-256 and is not
   counted as a new curve;
4. every j-invariant, canonical model, full ICV1 identity, EC1 alias, and full
   curve UID is unique within the strip;
5. replay recomputes every modular-polynomial choice, kernel, Velu codomain,
   target isomorphism, canonical model, identity, detector, and record chain;
6. all order audits report `proved_prime`, every model is nonsingular, and
   every deterministic generator is valid;
7. an exact streaming union audit checks both certificate chains and summaries,
   parses canonical 256-bit j values, and reports 1,064,000 unique values with
   no duplicate; and
8. a deliberately overlapping small-window fixture is rejected by the same
   union path.

The union audit establishes exact identity accounting for certificates already
subjected to mathematical replay. It is not a substitute for replay and makes
no claim about unexecuted coordinates.

## Frozen structural screens

The prior thresholds remain unchanged:

- primary: `qr_prefix_64 >= 54`;
- comparison only: `qr_prefix_64 >= 52`;
- canonical signed `b` bit length at most 224;
- canonical signed `a` bit length at most 224 when no `a = -3` model exists.

A threshold hit selects a curve for a separately preregistered factor-base
experiment. It does not establish relation yield, a discrete-log solution, or
an end-to-end speedup. A repeated full-width j value is a search discrepancy;
the run stops and retains both locations. It is not adaptively replaced.

## Cost boundary and stop rules

The unit is one globally unique j-invariant with an independently replayed
kernel-certified incoming edge. Coverage is measured against the requested
`2^32` screening boundary.

| population | unique curves | fraction of `2^32` | ratio to boundary | ECDLP speedup |
|:--|--:|--:|--:|:--|
| frozen million prefix | 1,000,000 | 0.023283064% | 1 / 4,294.967296 | unset |
| planned new strip | 64,000 | 0.001490116% | 1 / 67,108.864 | unset |
| planned cumulative union | 1,064,000 | 0.024773180% | 1 / 4,036.624060 | unset |

The last two rows are plans, not measurements, until receipts pass.

Stop and preserve the attempt if generation or replay exceeds 60 minutes,
available temporary storage falls below 512 MiB, a process is killed, a root
choice fails, any duplicate or identity gate appears, or the source tree is not
clean at the recorded commit. Record source and binary hashes, exact commands,
exit status, wall and CPU time, peak RSS, host, compiler, thread count,
certificate hashes, summary, and union receipt. Do not overwrite a failed run.

## Reporting and transfer boundary

The result report must include an editable lattice diagram, a quantitative
coverage graph, and a PDF containing both. It must distinguish proposed,
derived, replayed, and measured statements. Canonical scoreboards change only
if this enumeration changes a value they are defined to show; no ECDLP
scoreboard row is created from curve coverage alone.

Every degree-11 or degree-13 edge is required to be a separable cyclic
isogeny over `F_p`. Since both degrees are coprime to the prime P-256 group
order, the induced map on that rational subgroup is an isomorphism. This is a
subgroup-preservation obligation, not a destination-solver advantage. Relation
generation, matrix work, target transport, recovery, matched Pollard rho, and
end-to-end cost remain unset unless separately measured.
