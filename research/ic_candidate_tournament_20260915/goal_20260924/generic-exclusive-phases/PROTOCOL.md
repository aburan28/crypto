# Generic exclusive-phase instrumentation

This follows the query ledger in PR 799 and the supplied-point/outer interval in
PR 803. Keep the archived optimized producer and its sealed phase schemas
unchanged. Add an opt-in generic-library tracer and a separately versioned
generic report. Missing stage evidence cannot become an end-to-end comparison.

Partition actual executed work into setup, factor-base construction, reusable
precomputation, query generation, PDP, relation verification, matrix construction,
final scalar-field LA, target query/PDP/relation checking, target descent and scalar
replay. Boolean Macaulay reduction stays inside PDP. Every failed query and LA
attempt remains charged. Internal witness verification may require scoped phase
boundaries; do not invent a nonzero stage or count the same work twice.

All phase boundaries belong to the worker's controlling thread. Instrumented
collection uses the identical deterministic probes in serial order, with the
same witnesses and outcomes; one-thread internal algebra may still use Rayon.
Uninstrumented library callers keep their existing parallel collection API.
The measured worker fixes one Rayon thread and the Linux executor fixes one CPU.
Thread-local clocks must never silently ignore work executed in another phase.
The worker's strict session rejects a phase boundary from another thread instead
of silently leaving that work in its parent's interval.
Measure tracer effects with paired enabled/disabled controls on the same source
before performance qualification; passing clock closure alone does not establish
low overhead. Preserve the original worker interval mode behind an explicit
opt-in while the new stage adapter is under construction.

Before qualification, require exact online phase closure and whole-process
instruction closure, with process launch/input/report tail explicitly assigned
to setup. Preserve the raw native intervals and Callgrind dumps. Trace-disabled
calls must not issue profiling requests. A scope guard must disable tracing on
ordinary return and errors. Source identity, actual arithmetic/solver dispatch,
factor-base census and query/matrix replay remain separate admission gates.

Freeze controls at K_0 and K_1 over F_(2^9), one supplied public point from fixture
seed 2026092556, seven built-in PDP backends, dense/sparse final LA, and one thread.
Reuse the frozen 47-job query-law control panel, including its n9/n13 orbit-base
windows, and add rho on both n9 curves. Run three fresh process repetitions of
each job in both legacy and exclusive modes, alternating pair order. These 147
pairs quantify whole-mode effects (including serial outer collection), not pure
timer overhead. Retain all raw clocks and failures; no performance winner is
selected from these uncalibrated controls. All reported ratio diagnostics must
name this limitation.
Compare observed and unobserved query/certificate sequences before using timings;
test incomplete preparation and terminal failures. Stop on failed closure or
correctness. This is instrumentation admission, not a speedup experiment, and
does not consume either of the two remaining improvement rounds.

The exclusive policy assigns query scalar/state construction to query phases;
encoding, search, Boolean Macaulay reduction, sign-lifting arithmetic and pair
lookup to PDP; terminal enumeration membership/equality, lifted-candidate group
equality and explicit pair-witness sum checks to relation checking. The latter
also includes collector re-verification and duplicate filtering. Matrix row
construction is separate from final subgroup-field solving and log-table group
certification. Descent modular recovery and scalar replay have distinct phases.
This partitions existing executed work; it adds no duplicate mathematical check.
The raw tracer leaves every unentered phase null, including phases absent from
rho or a failed preparation. Any later scientific ledger must justify each
intentionally absent stage rather than globally replacing null by zero.
