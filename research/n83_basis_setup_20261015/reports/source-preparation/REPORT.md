# Prospective source and base preparation

No new training or holdout timing has run. This report records preparation
checks, not an IC/rho comparison. See ../../PROTOCOL.md for the exact boundary,
selection rule, seeds, all55planned processes and118planned answers.

Five candidates compare wide versus dense layout, a disjoint600-orbit window,
and enlarged1271/1798-orbit prefixes. Native candidate auditing verified every
point/orbit label, exact prefix/window recipes and zero signed-orbit overlap
between retained K600 and the new600 window. The draft first600 candidate was
identical to K600 and was prospectively corrected; its dump is retained.

Native fixture generation independently verified old16/32target artifacts,
replayed all64 dense-panel scalars, regenerated the1798-orbit public prefix
using independent half-trace arithmetic, and excluded1911 distinct signed
orbits. It generated two training points with seed202610156501. No training
answer has been sent to a measured worker; fresh holdouts do not yet exist.

Instrumented copies export per-attempt rows, witnesses, rank gain, probes,
exclusive clocks and endpoint memory. One synthetic algebra fixture and14
negative transcript mutations pass. Matching query clocks include generic
Q-validation, solve, generic scalar replay and answer write/flush on both sides.
Whole-process wait4 wall/CPU/RSS remains the primary accounting boundary.

The first full library run had15 environmental failures: missing sparse
fixtures and sandbox-denied loopback sockets. These were not algorithm failures;
the missing committed fixtures were materialized and the complete suite was
rerun with the network capability required by its local socket tests and passed
4089tests with0failures and117ignored.
Compilation failures before that (missing compile-time sparse files and the
new enlarged JSON macro's recursion limit) and the first audit-test harness
failure are retained under ../../checks. No test or guard was disabled.

Pending: final source freeze, measured training, independent audits, selection
freeze, fresh holdouts, native statistical analysis, graphs and inspected PDFs.
No winner, calibrated floor, L2, independent-host, asymptotic or security claim.
