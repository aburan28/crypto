# Same-field non-Koblitz follow-up (frozen before execution)

The first panel's signed-Frobenius arm exhausted its fixed budget on 13 of
40 one-target workloads. Its complete, audited session remains unchanged.
The generic rho and BSGS arms each verified all 40 targets. This follow-up
tests whether GHS structure and generic ECDLP costs differ when `b` changes
at the **same field degree**. It does not reuse the unsuccessful arm as a
reference for a new claim.

Inputs in `non_koblitz.spec.json` are five exact, already registered binary
models: the two Koblitz controls at degrees 13 and 15, plus three generated
models with `b` outside `F_2`: one at degree 13 and two at degree 15. The
registry records their complete polynomial-basis field, group order, prime
subgroup, generator and EC1 representation. All run as `binary_explicit`
instances to avoid curve-search work inside measured children. Each curve
gets eight independent one-target workloads at seed 20261009, one round,
no warmup, alternating order, L0 operation counts only, 60-second timeout.

The methods are `rho.negation` (generic matched reference on every model)
and `bsgs.negation` (time-memory comparison). The boundary is the harness's
generic square-root floor. The unit is `ecbench.gae` and the reported
normalisation is `S = total_gae/sqrt(r)`. Both methods include their charged
setup, search and internal scalar check. Hash-table operations and other
unpriced native work remain identified in the records; no physical speedup
or index-calculus cost is inferred. The 13-bit and 15-bit within-field
comparisons are descriptive with eight targets each and one round.

Run `ghs_screen` and `ghs_transport` for each new model and every proper
field tower. With `a,b` outside the chosen subfield, `ghs_transport` should
report the curve-model precondition rather than an image. A higher-magic
GHS row is a structural lead only. Success requires both ECDLP arms to
recover and verify all 40 target scalars, an exact replay audit, and saved
structural output. Failures, timeouts, or mismatch remain in the record.
Do not change this spec or reinterpret incomplete runs as successful ones.
