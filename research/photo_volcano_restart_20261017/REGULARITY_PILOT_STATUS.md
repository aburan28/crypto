# Exact-resolution capability probe: stopped

The [frozen four-system probe](REGULARITY_PILOT_PROTOCOL.md) produced no
regularity values. Singular's resolution and regularity routines reported
`not implemented for rings with rings as coeffients` for the large-modulus
coefficient backend. The native runner correctly returned 2 and audit.passed
is false, even though Singular itself returned 0. The retained launcher also
records a Node spawnSync EPERM error; no successful launcher claim is made.

[Raw stdout](sessions/regularity_probe_1820/stdout.txt),
[audit](sessions/regularity_probe_1820/audit.json),
[resource receipt](sessions/regularity_probe_1820/resources.json),
[exact executor relay](sessions/regularity_probe_1820/executor_relay.txt) and
[director assessment](cards/CARD-regularity-20261017-2.assessment.yaml) preserve
the complete stopped attempt. The executor performed no retry. Compile-only
card 1 caught an absent Node path before execution; its failed specification
and the frozen invocation amendment remain preserved.

The original 14720-system affine Groebner ideal-equality replay still passes.
Hilbert cutoff 13, temporary polynomial degree 23 and pair sugar 26 keep
their existing definitions; this probe does not turn them into CM regularity.
The [certificate alternative](REGULARITY_CERTIFICATE_PLAN.md) is an unexecuted
plan with explicit hypotheses, proposed deductions and independent checks.
Its conditional values are not entered into a results table or heatmap.

All original rho/BSGS cases and their separate PARI audit are complete.
Arithmetic replay is running; genuine random-corpus IC recovery remains open.
