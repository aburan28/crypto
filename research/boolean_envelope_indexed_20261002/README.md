# Indexed support-envelope continuation

The retained parameterized envelope study proves exact matrix reuse when
Boolean-polynomial coefficients change, including cancellation, degree drops,
zero/constant generators, newly eligible multipliers and support escapes. Its
frozen performance result is negative. This continuation tests the concrete
remaining construction hypothesis in [PROTOCOL.md](PROTOCOL.md).

The original five constructors are ported into one Rust binary with two new
arms. `envelope_degree` precomputes degree-cutoff lists in the original numeric
multiplier order. `envelope_indexed` additionally compacts actual output
columns through a predeclared universe and per-application occupancy bitmap.
Every returned row and column is still checked against the independent
set-parity oracle. The generator and caps remain bounded at n<=12.

The constructor currently passes the inherited exhaustive 8,192 coefficient,
degree and active-mask cases, focused fallback/cancellation/cap tests, and a
deterministic plan-count check. At n=12, a quadratic generator visits 13
eligible multipliers instead of scanning all 299 envelope plans; the same
compiled support still admits all 299 when a generator becomes constant. This
is an operation-count fact, not a wall-time speedup. Native campaign production,
independent replay, resource-qualified timing, fresh holdouts and the dramatic
gate are pending. No production solver or IC/rho result is claimed.

The old `research/boolean_support_envelope_20260922/run_01` files are immutable
historical evidence. New performance work will use a native Rust producer and
verifier with the repository's required external isolation controller.
