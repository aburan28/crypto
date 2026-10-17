# Reverse-order replication of the N83 K3350 factor-base lead

This experiment tests run-order sensitivity of the selected sparse root11
factor base. It repeats the exact 20,000 public points and producer source of
the first K3300/K3350 comparison on a new cloud instance. Each arm starts
from empty base, root index, and log caches. **K3350 executes first and
K3300 second**, reversing the original order. The known-scalar answer
sidecar is excluded from producer inputs and used only for independent replay.

The target and answer SHA-256 digests are
`fc8be4c2c946bdbcb78a6845c8e2be3d1a51728c2ed4d812f6546ea5b000055e`
and `7121dcaf478061b87f54f28a0dc8e5af28e2637fd34a6837dd8f6434b6e29e8b`.
The source file hash is
`c9824acebda0360416c3c507594bb07c17ca798c32cbf8f84bba95ca2c44d5b1`.
Both arms use sparse canonical encoding, 11-byte root states, the exact
five-byte indexed root table, four CPU workers for rank preparation, one
target worker, and Bloom disabled. The 14 GiB RSS gate and 1,800-second
wall cap per arm remain fixed.

The verifier independently reconstructs both geometric bases, replays all
40,000 scalar answers and four-point relations, and checks base hashes,
root-table entries, exclusive online phase sums, resource status, and exact
source/target receipts. The paired statistic is K3350-minus-K3300 mean online
time over 400 contiguous blocks of 50 targets. Use 10,000 bootstrap
resamples with seed `2026101519`. The reversed-order run replicates the
K3350 lead only if both arms pass every replay/resource gate and the upper
95% endpoint is below zero. Preserve a failed, capped, or inconclusive run
with that status and retain the first run's result as its own record.

This is a cloud-instance factor-base stage replication. It tests one source,
base family, public panel, and reversed process order; setup and peak RSS
are recorded separately from the target online interval. The next
end-to-end Tier 1 run uses the chosen base with sparse plus Bloom on unseen
same-point IC/rho targets.
