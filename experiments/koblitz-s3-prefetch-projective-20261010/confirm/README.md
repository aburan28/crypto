# Independent N53 C1200 combined-kernel confirmation

This panel uses T197–T204, separate from the T189–T196 selection screen.
T197–T200 were already frozen in the original public-point law. Checked Sage
extends that exact SHA-256 scalar law to T201–T204 and verifies the existing
T200 fixture before writing the four new public points and workload records.
`GENERATED_TARGETS.json` and `sage_runtime_info.json` bind that step.

Each point receives two independent processes per arm in ABBA or BAAB order,
alternating by point. Both arms use 14 preparation workers, one target worker,
the same 32 GiB RSS envelope, and identical one-target timing boundaries.
The reference and combined candidate keep their already screened source,
binary, factor-base, and algorithm identities. The runner stops on the first
invalid or over-budget run and retains all completed rows and raw reports.

The point, rather than a process repeat, is the independent sampling unit.
For each point, the online ratio is the geometric mean of the two reference
times divided by the geometric mean of the two combined times. The frozen
local confirmation gate requires all 32 native and independent Sage scalar
replays, the same first relation in every pair, a geometric online ratio of
at least **2×** across eight points, and a target-bootstrap 95% lower
endpoint above 1×. A 95% lower endpoint above 2× is reported separately as
evidence that the full 2× threshold is resolved rather than just reached by
the point estimate. Per-arm A/A spread and setup-time drift are independent
host-noise checks. A controlled wall-time claim still requires host-level
CPU isolation, affinity, frequency, NUMA, IRQ, and contention receipts.
