# Mixed-orbit search supervisor accounting revision

**Source-only validation passed at commit `ec48686c9`.** The existing
mixed-orbit search supervisor now requires a run-specific v2 additional-budget
record anchored at the **12,600-second charged total**, including the earlier
900-second import gate. It byte-pins that gate's passing audit as
`76a8bf74e46feb0ef4fbe650395f39733f4ff95d31c35f8829cb506a9ad04c7d`
and rejects a missing, altered or stale authorization before Docker. Build
and run both recheck the committed audit. The [machine validation](validation.json)
binds the source commit, changed Python source and tests, revised protocol,
passing [46-test study-suite log](study-python.log), and earlier import audit.

This revision launches no retained-N83 job and charges no new experimental
seconds. The [900-second import result](../mixed-orbit-import-900s-20261012/RESULT.md)
remains the latest measured mixed-orbit stage. The earlier
[source-validation receipt](../mixed-orbit-supervisor-source-20261012/RESULT.md)
is preserved as a historical check of the pre-import supervisor; its
11,700-second protocol is superseded by the [current protocol](PROTOCOL.md).
The updated two-residual search still needs at least 2,400 newly authorized
local seconds and a clean committed source before its first launch.

Partial yield, graph completion, independent rank, sparse linear algebra,
individual logarithm, complete cold single-target runtime and matched rho
remain unmeasured for this mixed-orbit arm. There is no evidence here for a
minimum-runtime factor base.
