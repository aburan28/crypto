# Test F6 geometric pruning before inherited-basis specialisation

Registered before implementation or timing. The prior T7 profile found
105.4–105.6 ms in inherited-basis specialisation and 321 geometric
refutations across its eleven target PDP attempts. In the current
solver, a child basis is specialised before the node oracle can refute
or witness that node. This opt-in pilot asks the same oracle once at
node entry, before specialisation, then skips the duplicate first
in-loop call. It keeps later oracle calls after propagation in the same
order. The baseline retains the old order. The hypothesis is that
refuted nodes avoid otherwise wasted basis work.

Freeze two new exact IC1 candidates differing only in `early_oracle`.
Use the prepared n17 Koblitz curve, 62 actual usable factor-base
points, 29 folded columns, archived public T1/T7, imported certified
logs, three summands, degree three, 8,192 node budget, 32 target
trials, one Rayon thread and default algorithm environment. Disable
prior opt-in row/column/pack and profiler flags. Freeze binary, source,
candidate and input hashes before timing. Run paired fresh-process
baseline/candidate repetitions in ABBA order per target; retain every
status, failure, timestamp and raw stdout/stderr. Pair by exact target,
limits and resource envelope.

Before timing, test oracle call count and decisions on a deterministic
small Boolean system, and run the prepared F6 geometric-closure
control. A successful arm must complete, replay the archived scalar,
and satisfy the five-phase exclusive online sum. Require identical
target, recovered scalar, attempt count and per-attempt result status;
allow F4 operation counts to fall because the optimization deliberately
skips work. Record oracle calls/refutations/witnesses, F4 calls and
word operations, complete online wall, target PDP and matrix phases,
including failures. A changed solver outcome or unverified scalar is
a correctness failure.

Retention gate: both T7 paired complete online ratios must be at least
1.10 (baseline/candidate), no T1 paired ratio below 0.95, and all
outputs correct. Otherwise keep the feature default-off. The Mac host
is unisolated; even a gate pass is exploratory until an isolated-host
receipt. No n83 ordinary relation or one-target IC-versus-rho speedup
can be inferred from this small-curve stage pilot.
