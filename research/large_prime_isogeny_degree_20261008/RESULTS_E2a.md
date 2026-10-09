# E2a result: partial; the registered instrument did not finish

Run 2026-10-08 with `relation_lattice.gp` on registered P-256 under the
registered six-hour wall cap.  PARI/GP 2.17.3, one core, machine shared
with the Φ_509 timing for part of the window (27% CPU over the six hours
as reported by `time`).

Outcome: **`bnfinit(y² − D_π)` did not return within the cap.**  The stack
grew to 2 GB and no class-group, factor-base or relation output was
written.  Per the protocol this is reported as partial, no prediction
E2.1 to E2.4 receives a verdict, and E2b is not run.

Reading.  `bnfinit` builds the full Buchmann machinery for the number
field; for an imaginary quadratic field of a 258-bit discriminant the
dedicated `quadclassunit` is the right tool for the class number and the
cyclic structure, and discrete logarithms against its generators need a
form-based routine the instrument does not yet have.  The instrument was
chosen for `bnfisprincipal`, which is why it stalled.

Not registered, reported separately: a `quadclassunit(D_π)` run for the
class number alone, with a two-hour cap, to see whether E2.1's window
`log₂ h ∈ [120, 136]` is even testable this way.  Its result, if any, is
appended below and carries no verdict on E2.1 since the instrument
differs from the registered one.

Next step for E2a as registered: replace `bnfinit`/`bnfisprincipal` by
`quadclassunit` plus a discrete logarithm in the form class group
(baby-step giant-step on each cyclic factor, or PARI's `qfbsolve`-free
route through `bnfinit` on a *reduced* precision setting), re-freeze the
instrument, and rerun under the same cap.

## Appendix: the unregistered `quadclassunit` run, 2026-10-09

`quadclassunit(D_π)` with `parisizemax = 2³³` did not return within its
two-hour cap either (exit 124, no output).  `time` reports 24 CPU-minutes
over the two wall hours, so the process was mostly not computing; the
likely cause is PARI stack growth and reallocation under the 8 GB
ceiling rather than the arithmetic itself, but this was not diagnosed.
No class number for P-256's Frobenius order is available from this
session, and E2.1 remains untested.  Before any rerun: fix `parisize` at
the start instead of letting the stack grow, run on an otherwise idle
machine, and record CPU time as well as wall time in the instrument.
