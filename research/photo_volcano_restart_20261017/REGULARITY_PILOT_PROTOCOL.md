# Generator-homogenized IC-system regularity convention probe

This is a separate diagnostic, not an IC solver or a replacement timing.
Use the exact existing audited field and polynomial systems. Homogenize the
three supplied generators separately in z, obtaining J in F_p[x,y,z]. Compute
the graded free resolution with Singular mres(J,0), then regularity(M).
This is reg(J), not the regularity of the saturation/homogenization of the
whole affine ideal, nor the earlier top-forms Hilbert cutoff. If reg(J)>0,
the corresponding quotient regularity is reg(J)-1; do not silently conflate.

The convention control is (x^2,y^2,z^2), whose ideal regularity is 4. Compare
both the resolution builtin and mregular.lib regIdeal on that control.
Then probe exact global system IDs 0,1,2,22 from singular_all.sing: two planted
controls, one random setup, one random target, all on curve 0/generator 0.
No research-case output value is prescribed. Preserve all four observed
values, per-case CPU timer ticks (1000 ticks/second), whole-process wall/CPU,
RSS, faults, switches, stdout and stderr. Run one CAS process with cpus,
threads and flint-threads all 1 and public CAS seed 20261017. This seed is a
reproducibility setting, not instance-generation randomness.

Freeze this protocol, source runner, exact input, original full audited input,
Singular and native resource wrapper before execution. New source does not
edit original IC journals, fixtures or bindings. Retain the earlier unexecuted
API-input variant. The output directory sessions/regularity_probe_1820 must
not exist; the runner claims it atomically. No retries or output overwrites.

Invoke /usr/bin/node source/run_regularity_probe.mjs with the absolute campaign
root as its sole argument. The command card uses a 600-second infrastructure
timeout for this four-system capability probe, not a campaign spending cap.
Stop on nonzero status, a convention mismatch, missing cases, nonempty stderr
or Singular error output. Record a bounded infrastructure/capability failure
as such; it is not a security result or a proof about other IC algorithms.
Only after inspecting this probe may a separately frozen all-system run be
specified. The pilot is not all-curve regularity coverage.

Timing class is concurrent CPU-only L0 with the existing solvers/replay.
No fastest, isolated timing or random-corpus IC recovery claim is made.
Conductor is unavailable (exit 127); that limitation does not require another
user confirmation for this already requested diagnostic.
