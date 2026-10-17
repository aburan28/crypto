# K3072 target cancellation on the N83 root11 factor base

Combine the independently checked `KIC_TARGET_CANCEL=1` lookup path with the
K3072/root11 factor base. The existing K992 diagnostic showed a paired online
gain from cooperative cancellation; K3072 has a larger exact root index and
fewer expected probes, so this experiment measures the combined algorithm on
the same 128 public N83 points. This is an ordinary-host candidate diagnostic.

Freeze the executed Rust source, release binary, Python runner, verifier,
environment flags, 128-point input and answer sidecar before timing. The
factor base has 509,952 distinct usable points and 3,072 folded columns. Run
serial `KIC_TARGET_THREADS=1,KIC_TARGET_CANCEL=0` first, then two-worker
cooperative cancellation `KIC_TARGET_THREADS=2,KIC_TARGET_CANCEL=1`, each in
a fresh process with an empty base, index and log cache. Cap each process at
1,200 seconds. Preserve full stdout, stderr, base dump, target rows, exact
command and run receipts, including timeout and failure status. Do not retry
an ambiguous producer under the same run identity.

Independently reconstruct the factor base, rank, four-point relation and
scalar on every target. Check the exclusive target phase sum for each target,
and require matching base hashes across policies. Compare 128 paired
target-online intervals by geometric ratio and report the cold process
costs separately. A verified local ratio above one is the predeclared gate
for the sealed 24-point Tier 1 same-point IC/rho holdout. The local screen
alone does not choose a hardware-qualified winner.
