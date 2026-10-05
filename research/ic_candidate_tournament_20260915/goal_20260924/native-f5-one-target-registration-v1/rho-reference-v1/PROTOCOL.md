# Same-point strong-rho diagnostic, version 1

This is a **one-target exploratory reference** for the already frozen F5
target `[61889,74818]` on the same n17 curve and generator
`[43693,23339]`, subgroup order 65587. It is registered after the F5 result
and cannot retroactively turn that disclosed target into a fresh, randomized
paired campaign. It does not reopen or rerun the F5 preparation or target
capsules. The point is supplied directly; no fixture scalar is constructed.

Build `examples/ic_n17_same_point_rho.rs` from committed source with the locked
offline Cargo dependencies and record the exact source, binary and Cargo.lock
hashes. Run this source once through the shared busy wrapper on the same
physical Apple M4 Pro. The algorithm is the repository's strongest existing
`StrongRho` reference: signed-Frobenius orbit quotient, normal-basis lockstep
batch inversion, 32 lanes, four distinguished-point bits and default
2000× step-cap factor. Freeze jump seed 2026100502 and walk-start seed
2026100503. It uses only one target, one worker, no cross-target collision
table and no reused target-dependent work. The jump table and curve
construction are target-independent preparation.

The online interval starts before target-point validation and normal-basis
conversion, includes every walk and collision, and stops after an independent
polynomial-basis scalar replay against the supplied public point. Record
exclusive validation/conversion, walk and recovery-check nanoseconds and
require exact closure to the online interval. Retain the raw result even on
step-cap exhaustion. The distinguished-point entry count is a final table
size, not a measured peak RSS; retain operation counters as diagnostics and
do not compute an `S` score without a fixed accounting conversion.

This diagnostic has no separate one-use frozen rho controller or host-level
isolation receipt. Its result must keep source-bound execution admission,
fresh-paired qualification, headline and speedup claims false or null. In
particular, an arithmetic ratio of two unisolated wall times is not a
controlled IC win/loss. A later paired fresh-target campaign must separately
freeze both arms, their common resource envelope and run order, same public
points, target-dependent boundaries, raw failures, and an isolated-host noise
receipt. Keep the existing F5 and SAT incomplete evidence unchanged.
