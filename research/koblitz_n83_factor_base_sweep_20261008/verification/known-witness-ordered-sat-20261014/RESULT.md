# Ordered N83 fixed-witness SAT control

This control fixed the five summand x-coordinates and three affine partial-sum
x-coordinates of the independently replayed known witness in the worker's
first (`identity_mask=0`) chained-S3 pattern. It used the full K1182 small
domain and the 866-orbit residual suffix of K2048, including the worker's
first-pair order. The exact signed-Frobenius domains came from the two retained
S3 archives; [domains.json](domains.json) has 1,182 and 866 distinct
representatives. A second extraction of both archives reproduced its SHA-256
`63accc5b9dc9de6327c5b07e93b3e62717f8376a96c3f91f7bb70b77095f23bb`.
The [current-content S3 audit](../s3-stream-integrity-20261014/RESULT.md)
independently streamed and hashed every retained object.

The release probe source is frozen at commit
`55d5ed6771ad53ff6dcfe6e537ac3ec4598984bb`; its SHA-256 is
`06730771fda845339fdd14908ded0e259e0294ebc37d62af39e5444bc6079088`.
The exact run source is retained as [worker_source.rs](worker_source.rs).
The current example is its `rustfmt` rendering; the auditor checks that
equivalence against the original Git blob.
The exact native arm64 [worker](worker) has SHA-256
`64ce185bfb1a4db69dcc576fb0749108b56e804efd57f14da561f07d41f51752`.
Its source requires the same 465,738 domain clauses as the prior v7 control
and unit clauses for the eight fixed x-blocks before solving. It was launched
once with 5,000 conflicts and an external 840-second wall watchdog. The launch contract
is [launch.json](launch.json); the raw [stdout](solver.stdout) and
[stderr](solver.stderr) are both empty. The watchdog reported exit 124.
Its [exit record](exit.json) preserves that terminal observation and a
pre-timeout `ps` diagnostic at 12:30 wall time, 2:03.47 CPU, and 67,088 KiB
RSS. The watchdog terminated the `/usr/bin/time` wrapper too, so there is no
final resource summary, native construction marker, or solver counters. The
terminal exit code is transcribed into the retained exit record; the raw files
alone cannot replay that code. Consequently this receipt cannot establish
whether construction completed or the solver began searching.

The independent [audit](audit.json) binds the source, binary, prepared
domains, prior algebra replay, current S3 content receipt, launch and exit
records, and raw output hashes. It classifies the call
`VALIDATED_CENSORED` / `UNKNOWN_wall_cap`. This establishes neither SAT nor
UNSAT for the fixed witness. In particular, the timeout supplies no natural
relation yield, useful rank, linear-algebra time, individual logarithm, or
complete one-target cold runtime. The local N83 ledger conservatively charges
the full additional **1,200 seconds**, raising the cumulative charge from
21,000 to **22,200 seconds** with none remaining under this grant. Measured
process wall time is separately capped at 840 seconds; the two quantities are
not interchangeable.

Recheck the receipt from the repository root with:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/known-witness-ordered-sat-20261014/audit.py --check
```

This arm64 control is an encoding-access probe. Its wall time is not a matched
Linux one-CPU/cgroup runtime comparison with the earlier unrestricted v7
control. The unchanged [ordered algebra replay](../known-witness-first-pattern-replay-20261014/RESULT.md)
still verifies the group/S3 witness, and the earlier unrestricted call remains
censored at its own conflict cap. A factor-base winner requires observed
natural relation yield, rank, linear algebra, a fresh ordinary target and
matched complete cost accounting.

[Validation](VALIDATION.md) records the current focused checks, the earlier
complete release-library result, and the host-contended full rerun that was
intentionally stopped without a terminal full-suite pass.
