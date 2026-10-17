# N83 current-source small-base controls

## Question and frozen comparison

On the previously exposed N83 public point, K1248 and K1280/root11 dynamic
were 3.344 and 4.119 ms online in the 2026-10-16 cloud grid. Replay them
under the current finalist IC source and one local release binary, sequentially
on the same host as K2543, K2816, and K3072. The hypothesis is that at least
one small base retains a shorter verified one-target online interval than
K2543's 18.983 ms in this diagnostic. This screen cannot select a fresh-point
winner because the point and earlier timing results are exposed.

The input is the exact public target and answer in
`research/n83_cold_factor_base_20261016`; the field, generator, subgroup, base
construction and root11 dynamic flags come from the frozen 2026-10-16 runner.
Source SHA-256 is
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`
and local release binary SHA-256 is
`d71a9662c9d4e2dde58443bfc646a05737c6b9ca6b3d5e831b632ce09384c906`.
Each K gets one fresh process and an immutable output directory. Order is
K1248 then K1280, after the K2816 process exits. Producer cap is 600 seconds
per K. Full base reconstruction, four-point witness, group scalar replay,
source/binary hashes, exact candidate identity, and five exclusive online
phases are checked separately. Preserve any timeout or failure as a recorded
status; do not insert an unverified time in a comparison.

Report actual usable points, folded columns, root entries and memory, cold
process seconds, the five online phases and total, and replay status. The
primary next comparison is alternating same-point IC/rho on fresh sealed
points in Tier 1 cloud allocations. This local screen is a compatibility and
candidate-selection diagnostic; no observed timing is promoted to that
primary result.
