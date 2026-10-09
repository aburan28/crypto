# SAT/PDP matched-pipeline workstream

The first preregistered protocol is [`PROTOCOL.md`](PROTOCOL.md). Its
[`pilot`](sessions/pilot) closed with 28 records: 24 verified and four
`ic-buchberger` timeouts at the fixed 60-second per-run limit. The
strong-rho, meet-in-the-middle, Frobenius meet-in-the-middle, and both
SAT encodings returned verified answers on both public targets in
both rounds. The timeout stops the original eight-target main panel
under its stated rule; [`specs/main.json`](specs/main.json) remains
frozen and unrun.

The first pilot also exposed a cost-replay defect in the historical
`ic.pipeline` method. Its
[`failed audit`](pilot.audit.failed.json) checked every seal and identity,
then replayed 11 of 12 measured verified runs identically. One SAT
native-XOR replay had the same answer and integer counters but differed
in the low bits of `total_gae` and phase GAE. Receipt SHA-256:
`7a25f644ddc2fc8137d2a4e2c1fbd54993798eeb3b543c2d28a2effa3b6feef1`.
That failed receipt is retained; the original pilot is not used for a
speed ratio.

[`PROTOCOL-2.md`](PROTOCOL-2.md) is a new, narrower preregistration.
It drops the timed-out Buchberger arm and uses the versioned
`ic.pipeline_counted` method, which reconstructs charged relation
cost from counts instead of subtracting a wall-derived floating charge.
The method keeps SAT conflicts and other units without a pinned
conversion explicitly unpriced. Its pilot and holdout panel have
separate frozen specs. An audit and result note will be added as those
sessions complete; their absence here does not promote the first
pilot to a complete comparison.
