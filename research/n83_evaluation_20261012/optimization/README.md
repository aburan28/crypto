# Fair rho optimization on the exact legacy N83 workload

This follow-up implements the user's requirement to optimize rho as well as IC.
The bounded search, single configuration selection and fresh-target evaluation
are specified in [PROTOCOL.md](PROTOCOL.md). The original study's five-pair
comparison remains frozen evidence against default strong rho; it does not
establish that default rho is optimal.

The experiment is being prepared. No optimization timing result is claimed by
this initial document. Final status and results will be filled from the native
campaign receipts after execution and independent verification.

| Requirement | Initial status | Evidence |
|---|---|---|
| Audit existing strong rho | Verified implementation audit; optimization not yet measured | [Measurement review](MEASUREMENT_REVIEW.md) |
| Share existing optimized field primitives fairly | Implemented; four worker tests passed | [Configured worker](../../../examples/n83_rho_optimization_eval.rs) |
| Freeze training and held-out target law | Native generation complete; no timings yet | [Protocol](PROTOCOL.md) |
| Train on every declared candidate/case | Not started | Native runner under preparation |
| Freeze selection before held-out timing | Not started | Native runner under preparation |
| Default/selected rho on five held-out targets | Not started | Native runner under preparation |
| Same-target prepared IC online comparison | Not started | Original immutable IC binary |
| Independent scalar and decomposition replay | Not started | Standalone native verifier |
| Updated PDF and source-linked graphs | Under preparation | Reports beside frozen records |

This is a local diagnostic campaign. The exact legacy case still lacks the
sealed wide IC/ecbench adapter, calibrated common operation unit and independent
host replication. Shared IC preparation is reported separately and cannot
replace a cold single-target interval.
