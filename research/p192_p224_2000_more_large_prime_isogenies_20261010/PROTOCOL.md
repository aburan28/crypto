# At least 2000 additional large prime-degree maps

The requested deliverable is at least 2000 additional distinct independently
verified maps in the public P192 and P224 isogeny classes. Prime degrees exceed
1009. Repeated prime degrees are allowed; the distinct degree count is separate.
The earlier hundred-map campaign remains frozen and its maps do not count again.

| Requirement | Execution and acceptance evidence |
|---|---|
| At least 2000 additional maps | Two fresh 1050-step walks, subtracting every previously counted normalized map identity |
| Both public classes | P192 degree 2341 and P224 degree 1117; expected 1000 additional maps per class |
| Exact maps and subgroup action | Existing independent limb arithmetic, exact kernel and rational identity checks, and 20 public scalar checks per step |
| Distinctness | Source ICV1, prime degree, and canonical monic kernel coefficients; raw certificate hashes are provenance rather than the novelty key |
| Reproducibility | Clean source revision, binary and protocol hashes, command/exit/resource receipts, complete raw maps and replay certificates |
| Delivery | Extend ready PR 1606 with code, protocol, evidence, report, diagrams and PDF; update canonical models after verified findings |

Both walks start on the embedded published model and generator, seed 1, and
advance only after exact independent replay transports the subgroup generator.
The first 50 steps of each fresh walk are expected to reproduce the earlier
walks. They incur construction and replay costs and are explicitly excluded
from additional-map yield. Every later step is compared with all prior and
current normalized map identities. Merely requesting 2100 steps is insufficient:
the final unique additional count must reach 2000. A source-model revisit,
failed check, failed child, timeout or incomplete request remains unresolved.

The algorithms and acceptance gates are unchanged. The CLI request bound is
raised from 1000 to 10000 so these fresh 1050-step walks can be executed without
trusting an external source-model handoff. Coverage stays PARTIAL, one of two
Frobenius eigenlines. These are consecutive edges at two prime degrees; they are
not 2000 different prime degrees or 2000 direct maps from the published models.

Use the dedicated installed native executable, copied into a fresh campaign
directory. No Cargo, checkout or companion program is needed by that executable
at runtime. A native supervisor launches at most two walks, one per class,
retains stdout/stderr and per-walk exit receipts, and prevents system sleep
while executing. Constructor limits remain 1800 seconds and 8 GiB RSS per
invocation. Independent replay is outside those constructor limits. There is
no cloud allocation. Actual disk exhaustion stops work without accepting an
unsealed step. Raw output uses the system disk because the research SSD has
limited free space; final raw evidence must be archived in the PR or another
durable accessible location before completion is claimed.

Freeze the protocol, supervisor and CLI source in a commit before measured
execution. Preserve old unsuccessful attempts and old reports. Record separate
constructor, independent replay and complete invocation sums, including the
excluded repeated prefix costs. Timings are operational measurements on a
shared host and do not establish a matched performance comparison.

An integrity snapshot validates all linked file hashes and required replay
checks and reports the actual verified prefix. It does not redo the algebraic
replay. Exact independent replay already occurs inside each walk before its
manifest advances. Mark the campaign complete only when both child exits are
successful, both 1050-step requests are sealed, at least 2000 additional unique
maps remain after deduplication, and evidence and canonical reports are saved.

This is map-construction evidence with same-host independent arithmetic replay.
The IC performance ledger is unchanged. Runtime and release publication gates
remain separate from mathematical map acceptance.
