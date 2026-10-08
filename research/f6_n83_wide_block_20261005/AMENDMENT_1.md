# Bounded n83 root-search diagnostic

Registered after the eight exactness tests passed and before this search
is run. Use the same frozen K0 curve, dimension-12 **source** subspace
and planted source points `[0,2,4]` from `PROTOCOL.md`. Build the exact
119-variable m=3 wide system and perform one `wide_groebner_decompose`
call on that planted source sum with `node_budget=1`. This tests the
cost and outcome of the root linearisation; it is not a relation-yield
sample or a proof that the full high-arity solver is viable.

Run a native Rust release probe through `gtimeout 120s`, with a best
effort 8 GiB virtual-address limit if the host accepts `ulimit -v`.
Record build and solve wall intervals separately, source and binary
hashes, exact point replay if a witness is found, node/reduction/matrix
counts, RSS, stdout/stderr, and exit status. `found` must replay in the
curve group. `exhausted` and `unsupported` remain distinct from
`refuted`; a timeout (exit 124) is inconclusive. A result within the
120-second envelope with no OOM is the admission signal for further
algorithm design; it is not a performance win. A failure remains in
the PR and stops this diagnostic.

The source-point system is not the projected factor-base candidate.
The cofactor bridge tested in `PROTOCOL.md` will be required before an
ordinary subgroup relation can be counted. No `IC1` identity or online
IC cost is assigned to this block test.
