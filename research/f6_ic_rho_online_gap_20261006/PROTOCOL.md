# Same-target F6-IC versus rho online gap

Registered before new timing. This measures the default one-target
objective, not a batch average. The IC arm is the exact, verified
`IC1N17Ckb1fb62PDP3f6RCsampleLAgaussTDpdpISO0h9628dfd41b76`
candidate from the compact-refuted pilot with every experimental
optimization and profiler flag off. Its reusable certified factor-base
logs and symbolic preparation are ready before the online interval.
The rho arm is this same worker binary's signed-Frobenius rho solver on
the identical supplied public point, with one walk, one Rayon thread,
seed `20261004039`, at most 65,536 iterations per restart and the
source-default restart/collision policy. Both have a 180-second process
cap. Input loading, process launch, point fixture construction and
target-independent preparation are outside each online interval.
The IC interval charges target query, every target PDP attempt, relation
check, descent and scalar replay. Rho starts at the first prepared walk
computation and stops after scalar replay.

Freeze worker binary/source hashes, two exact public targets T1/T7,
their archived workloads and both input hashes before timing. Run two
fresh-process repetitions per target in IC/rho/rho/IC order on the
same unisolated Mac, `RAYON_NUM_THREADS=1`, no artifact cache and the
CPU field backend. Retain every raw output, failure, timeout and
timestamp. IC must retain its candidate/workload/run IDs; rho is a
separate exact reference record keyed by workload, target, binary,
policy and run repetition. Freeze both algorithms' resource limits,
walk/collision parameters and source identity in the run record.

Compare online wall only when both algorithms complete, recover the
same scalar and independently replay it against the public point.
Require the IC five exclusive phases to sum exactly to online wall and
rho's exclusive solve/recovery phases to do likewise. Retain rho
iterations, restarts, walk additions and peak memory if available, and
IC attempts, F4 counts and all failures. Report per-target paired
`rho_online_ms / IC_online_ms`, observed two-repetition ranges and
the complete timing boundaries. An unverified result makes that pair's
speedup unknown. The physical host is not isolated, so all CPU ratios
are exploratory and cannot be promoted under AGENTS.md's isolation
gate. This small-curve comparison establishes no n83 relation yield.
