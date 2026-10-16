# Native stage-probe watchdog and raw-receipt preflight

The earlier [stage-marked worker](../known-witness-stage-preflight-20261014/RESULT.md)
can identify whether a capped ordered-witness control reached SAT search,
but its first self-test did not exercise a process watchdog or replay a
complete launch/exit receipt. The new [launcher](launch.py) pins that worker,
its Rust source, the trace checker and the full K1182/K2048 domains by
SHA-256. It starts the worker in a new process group, writes raw stdout and
stderr as separate files, records child and watchdog exits, sends TERM and
then KILL on a wall cap, checks process-group cleanup, and writes an audit
that can be replayed from the saved bytes. If a trace is missing or malformed,
the audit is unverified rather than a SAT/UNSAT or yield result.

The retained [self-test evidence](selftest-evidence/audit.json) ran only the
worker's `self-test-markers` mode. It emitted two ordered markers, exited zero,
and passed independent replay with process-group cleanup verified. The raw
[stdout](selftest-evidence/stdout.raw), [stderr](selftest-evidence/stderr.raw),
[launch](selftest-evidence/launch.json) and [exit](selftest-evidence/exit.json)
receipts are preserved. This path does not read the full-domain input into
the worker or build an N83 SAT model. The conservative local N83 charge
remains **22,200 seconds**; natural relations, rank and total runtime remain
unmeasured.

The `--run` mode is deliberately gated by a separate authorization JSON with
schema `n83.fixed-witness-stage-authorization/v1`, the public-control scope,
a nonempty concise approval reference, prior charge 22,200 seconds, a new
authorized maximum of at least 870 seconds, an 840-second watchdog, 5,000
conflicts and the pinned domains SHA-256. The launcher copies those bytes into
the new evidence directory and conservatively charges the authorized maximum
in its receipt. **No such authorization receipt exists in this preflight and
`--run` was not invoked.** The JSON is an audit record of an actual user grant,
not a way for a local file to create one.

This launcher is pinned to the retained **native macOS arm64** binary for a
stage-localization diagnostic comparable to the earlier censored native
control. It does not impose a Linux cgroup, one CPU or zero swap; it cannot
provide a matched Linux runtime comparison. A later full cold-cost experiment
requires its own source-matched Linux worker, resource isolation, natural
relation and rank measurements, target replay and rho comparison. The host
had roughly 475 MiB free on this study volume during preflight; the launcher
refuses to start if the output volume has less than 256 MiB free, and a fresh
run still needs a headroom check.

Recheck the retained evidence and focused tests from the repository root:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/known-witness-stage-launch-preflight-20261014/launch.py --check research/koblitz_n83_factor_base_sweep_20261008/verification/known-witness-stage-launch-preflight-20261014/selftest-evidence
python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008/verification/known-witness-stage-launch-preflight-20261014 -p 'test_*.py' -v
```

The four focused tests cover the retained marker receipt, raw-output mutation,
a separate sleeping process killed by the watchdog, and refusal to launch
the full-domain path without a valid grant. A preliminary invocation with a
mistyped checker digest was refused before creating an evidence directory;
the pin was corrected before the retained self-test. The present code and
receipt do not claim a full-domain solver outcome or a factor-base winner.
