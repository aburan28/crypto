# Preregistered 60-target full-image setup-reuse comparison

**Status: frozen protocol, before batch-run implementation or measurement.**
This gate follows the [one-target full-image result](../different_j_full_image_20261015/RESULTS.md)
and the [pending reuse question](../different_j_full_image_20261015/NEXT_GATE.md).
It tests whether paying rational-image construction once over distinct
public targets changes the complete local cost comparison. It makes no
asymptotic or counted-`S` claim.

## Boundary and inputs

The reference is **direct source-group negation rho on the exact same
target prefix**, with a direct source BSGS arm that may be cheaper. All
arms report total current-thread CPU in microseconds, setup, peak RSS,
primitive counts, prefix total and total divided by distinct target
count. Generic collision methods have a square-root order floor in
subgroup size, but this single order-61 batch cannot estimate an
exponent. Portable `S=operations/√61` against the counted rho frontier
remains a separate gate.

The curve and map are the hash-pinned [public V2 fixture](../different_j_full_image_20261015/PUBLIC_V2.json)
over `F_(67²)`, source `j=7`, destination `j=11`, degree-three Vélu
map, source generator of exact order 61, and the four fixed eight-column
bases. The public target order is generated once from seed 13 and step
37: for `i=0,…,60`, compute `(13+37i) mod 61` and omit zero. This is a
permutation of the 60 nonzero scalars. A native freeze process writes
60 source/destination point pairs to a public JSON file and their
scalars to a separate private checker file. The solver may read only
the public file; reject any answer-bearing solver input. An independent
arithmetic checker verifies distinctness, order and all 60 map pairs.
The first `N=1,4,16,60` points in this fixed order form the prefixes.

## Arms, runs and cost custody

Run three complete rounds with a committed arm rotation over seven
arms: source projected, isogeny-image projected, destination-direct,
full-image, direct source BSGS, direct negation rho, and an independent
rho A/A stream. The four relation arms reuse their own full setup
inside each 60-target process: base construction, maps, pair table and
validation are charged once; each target then gets its own capped
relation collection, rank, linear solve, recovery and scalar check.
BSGS builds its baby table once and charges every giant lookup and
verification. Rho charges all retries and checks per target. Each
round starts from a new process and incurs its own input parsing,
field/curve construction and target validation inside the timer.
Process startup itself remains outside the native timer and is
reported separately if a wall wrapper can measure it.

The per-target relation cap is 1,220 trials, with deterministic
independent streams derived from round, target position and arm. The
rho walk cap is the existing frozen order-61 implementation's bound,
recorded in the V2 source and raw output. Set a 512 MiB RSS ceiling
and 120 s wall ceiling per arm process. Freeze source, binary, public
input, private checker, script and exact stream formula hashes in a
committed V2 spec before the measured run. Preserve every arm's raw
stdout, stderr and exit code, including failures; never overwrite a
failed version.

## Acceptance and decision

The solver receives no answer scalar. For each of 60 public targets
and every arm/round, require a recovered scalar, source and mapped
point checks, a separate native arithmetic replay, and exact target
distinctness. Replay relation witnesses and factor-base/map facts.
An exhaustion or timeout is a censored arm/prefix, not a zero-yield
cost. Report all setup, query, linear, verification, failed-attempt,
memory, field and group counts. One table uses the same CPU unit for
all seven arms and each prefix; ratios use paired round totals and the
cheaper direct BSGS/rho control. Report cold totals and genuine
amortized values separately. A relation-yield or amortized phase ratio
alone does not establish a complete DLP gain.

The result needs a native run report, quantitative graph, PDF,
`benchmarks/` mirror, independent replay and updated obstruction
record. Portable counted `S`, isolation, larger-order scaling and
second-host performance replay remain open even if one of these
prefixes has a favorable local CPU ratio.
