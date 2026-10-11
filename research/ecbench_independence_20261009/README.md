# Same-machine replay can change the environment class

Status: correctness regression, **not an independent validation**. This
note preserves the observed failure that motivated the hardware-class gate.

The [budget-qualified SAT/PDP session](../ecbench_sat_pdp_20261009/README.md)
was measured on this Mac with Homebrew `rustc 1.93.1` on `PATH`; its
`host.json` records `ECBENV2hd82681268e96`. On 9 October 2026, the same
Mac audited one of its runs with Rust 1.98 on `PATH`. The new receipt
[`same-machine-replay.audit.json`](same-machine-replay.audit.json) reports
`ok=true`, one reproduced replay, and
`auditor_env_class_id=ECBENV2hbb6f16e5c58e`. Under the old claim check,
the two different `ECBENV2` values would pass the host-class part of
`independent_validation`, despite both executions using the same physical
machine.

The source and auditor host capsules' `stable` objects are byte-identical
after deleting only `rustc_vv` and sorting JSON keys. The new, separate
`ECBHW1hb360e0ff59e1` hardware class excludes toolchain, kernel and
policy fields, so a receipt from this same Mac is rejected even though its
measurement environment class differs. The receipt's SHA-256 is
`430e1de5c5473f9627eb31f481436d5e58ca41aa00fcb6c94e82d5ab3b0bfa95`.

The auditor binary was built on an unpublished local validation merge of
Rust CI repair PR #1605 (`f18877bb6d9101b46a8503e8480d1282567cd247`)
and SAT/PDP port PR #1611 (`a148d750f0c3da0088ec4568df48f580394cabe7`),
plus the hardware-class change. The resulting auditor binary's SHA-256 was
`74a0259f5a478e0cd82ae89a65e10e32144ba98698dca50d3a10a3cc83e2d5ed`.
The receipt records that hash. This audit is a local
cross-revision check, not a second-host certificate. The full saved panel
also replayed 200/200 measured runs with zero audit problems on the local
merge before the hardware-class change; that receipt remains a local
validation artifact until a source commit on the accepted branch exists.

Validation on the local merge: release `ecbench` build passed;
`hardware_class_does_not_change_with_toolchain_or_kernel` and
`compiler_change_on_one_machine_cannot_supply_independence` each passed.
An external runner's custody record remains necessary: these fingerprints
are self-reported and do not cryptographically prove distinct physical hosts.
