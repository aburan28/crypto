# Corrected F5 v3 disclosed control registration

[PROTOCOL.md](PROTOCOL.md) and [panel.json](panel.json) were committed at
`2c492b402e537af250ec3f666c4d29798d75a767` before this new registration.
The actual accepted source is `7d483204c9f0359d1d8150c3e1a2e21908b409ba`:
runtime v3, input v2 and frozen transport v2. [REGISTRATION.json](REGISTRATION.json)
is the before-execution ledger. No native solver has run or claim been acquired
at this preregistration stage. Do not treat this historical stage label as a
later authorization: inspect the claim and published terminal result first.

| Frozen identity | Value |
| --- | --- |
| Candidate | `IC1N17Ckb1fb62PDP3f5RCsampleLAgaussTDpdpISO0h5cc365241aec` |
| Workload | `8b3d9065d872` |
| Run | `IC1N17Ckb1fb62PDP3f5RCsampleLAgaussTDpdpISO0h5cc365241aecW8b3d9065d872R2026100104` |
| External execution SHA-256 | `5c0a1473b7d2637a7e97666db68262d038c628f937f3b0ea0338fd83f949f98d` |
| Supplied point | `[52411,72106]`, already disclosed |
| Attempt cap / algorithm seed | 8 / `2026093032`, known correctness fixture |
| Native mathematical state | 63 geometric points, 62 usable points, 29 columns |

The complete [registration](registration/execution.json) retains 206 bounded
Python source components, the native code/build/assets and whole externally
sealed preparation before candidate/workload/run derivation. Native archive and
every nested root/dependency source byte are retained. The interpreter executable
and stdlib contents are bound by hashes; replay requires that exact local
interpreter, not a portable archived Python installation. Before-execution replay
reconstructed the mathematical arguments and verified runtime/asset seals,
current source/interpreter identity and absence of a claim. No timing or yield
was measured by registration.

Only after this preregistration PR is accepted with its exact-head applicable
checks, execute once through the busy wrapper, passing the external SHA above.
Use a new external output tree. The exclusive claim is created before launch;
failure, timeout or source rejection consumes the registration. Retain the claim
and all original output. Audit using that execution's preexecution-frozen
transport v2 in a separate output tree; never edit its raw native report or
restart its native job after an audit failure.

This point and seed are deliberately selected from known interface controls.
A successful native admission would establish only complete disclosed prepared
source execution, including its failed attempts. The existing natural ordinary
query evidence stays attached to the accepted preparation; there are zero new
ordinary queries in this job. Headline admissibility, fresh target qualification,
matched incumbent/rho speedup, promotion and full goal completion remain false
or unknown. All older registrations and all three confirmation sets stay closed.

## Terminal execution

The preregistration was accepted in PR #1159. Its one invocation is now consumed
and closed: [RESULT.md](RESULT.md) records complete source-bound recovery, the
original failed target attempts and exact frozen archive replay. Do not execute
this registration again. Fresh comparison and the full goal remain pending.
