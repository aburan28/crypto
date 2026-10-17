# Recovery-handoff validation

This is a documentation and metadata change, with no new research execution.

## Checked

- Retrieved the exact source shortlist at commit
  `ef18058588077ef805c55c5dc25140b6bde2046b` and checked the selected broad
  curve's parameter strings against the recovery metadata. Source blob:
  `bf361c12401f9d1768c9bd4a967fb7eb5714fdc5`.
- Checked absence of the original experiment, reference and temporary extracted
  checkpoint directories, searched the accessible workspace/shared and temporary
  locations, and read Docker container/volume state. Docker reads succeeded
  after enabling socket/network access and returned empty lists. These scoped
  checks do not prove no backup exists outside the checked storage.
- Parsed the recovery JSON; checked archive count, names, size types, lowercase
  64-digit SHA-256 fields, checksum-file agreement, workload arithmetic and
  decimal-string preservation for large curve values.
- Checked local documentation links, the problem/index link pair, whitespace,
  and that the patch is restricted to this recovery study and problem record.
  The publishing step reads back committed content to check exact equality.

## Unavailable or not performed

- The original archives were not available: neither archive-content hashing,
  extracted-file checksums, graph/map certificates, point-order checks nor
  solver replay could be rerun. Expected hashes are prior manifest metadata,
  not freshly computed archive hashes.
- No solver source, binary, raw timing data, report PDF or diagram skill was
  recovered or regenerated. No new runtime, relation-yield or regularity result
  is reported.
- `cargo test --release --lib` could not run because Cargo is unavailable in
  this runtime. The original checkout is also absent. No library test pass is
  claimed; there are no touched native examples in this patch.
- `conductor check` was attempted before editing the proposed repository
  paths; it exited 127 because Conductor is unavailable. The referenced
  `.conductor/WORKFLOW.md` also returned 404 at the pinned base revision. No
  reservation or knowledge of other workers' leases is claimed.
- CI/review status is supplied by the live PR, not this frozen note. Pending or
  unavailable checks are not passes, and this handoff does not authorize
  bypassing a merge gate.

The PR's metadata assertions and link checks validate this handoff only. They
are not experimental verification. Canonical measurement graphs, scoreboards
and benchmark mirrors are unchanged because no new result is being published.
