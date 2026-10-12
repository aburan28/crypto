# Amendment 2: stop m = 31 once its outcome is fixed

Status: written while the m = 31 and m = 23 sessions were running. At that
point m = 31 had produced one strong-rho run (verified) and one inherited-F4
run, which hit the protocol's 3,600-second cap and is recorded as a
timeout. No other m = 31 IC outcome was known. Nothing in the specs, caps
or decision rules changes.

## The rule

`PROTOCOL.md` makes a size decidable only if at least six of its eight
targets complete in **both** IC arms. A target with any IC run that did not
verify can never count. So once **three** m = 31 targets have such a run,
m = 31 can no longer reach six. Its outcome is then fixed: not decidable.
The protocol's fallback ladder, m ∈ {7, 13, 19, 23}, then decides.

When that happens, the m = 31 session is stopped with SIGINT. ecbench
records it as interrupted. Its executed runs, timeouts included, are kept
and reported. The executions that never ran are listed as not run. They
are not timeouts and not failures, and they are never counted as wins.

## Why

The stop condition is purely logical. It can end the session only after
the decision it feeds can no longer change, so it cannot favour any
outcome. It saves compute only: each further timed-out IC run at m = 31
costs an hour of the cap, about 16 hours in total, and changes nothing.

If m = 31 does reach six complete targets, the session runs to the end
and this amendment does nothing.

## The fallback size

The `SPEC-n7.json` written from the protocol's fallback definition was
committed (`a8c5f733`) before any m = 31 IC outcome was known. It runs now,
in parallel. Its result enters the decision only on the fallback branch.
Otherwise it is reported as an extra size outside the primary fit.
