# Prepared audit transport and one-use execution claims

This implementation follow-up supplies an isolated, preexecution-frozen audit
path for `prepared_sat_runtime_v1.py` and `prepared_f5_runtime_v2.py`. It closes
the shared launcher's ability to run one registration into different output
directories. It does not register or execute a new production SAT/F5 control,
generate a target, or qualify a family for the tournament.

The scientific goal remains a complete source-bound F4/F5 pipeline and a SAT
pipeline, natural ordinary-query yield with failed attempts retained, followed
by a newly frozen paired one-target comparison with the incumbent and rho.
The three historical confirmation sets and all consumed native registrations
remain closed. Work in this follow-up concerns the disclosed toy n17 control
and portable source-transport tests; it provides no larger-curve claim.

## Hypothesis and controls

The transport hypothesis is that a new retained prepared execution can be
audited with exactly the Python, interpreter and native-input bindings frozen
before execution, even after the caller's live checkout changes. The scheduling
hypothesis is that two launch attempts against the same retained registration
cannot both acquire its claim, including when they choose different outputs.

Use portable benign Python controls, with a synthetic known-name family auditor
for actual isolated CLI transport. These synthetic admissions explicitly say
that they contain no native solver, mathematical IC result or performance
evidence. The real family auditors' existing mathematical tests remain separate.
Run the existing prepared SAT/F5 and runtime tests to check compatibility.

Acceptance requires frozen-source transport, full before/after loaded-module
gates, exact external execution SHA-256 and candidate/workload/run agreement,
and rejection of changed artifacts, live imports, changed sources and output
reuse. Preserve rejected transport output and timeouts. A race for a claim must
have exactly one winner. A successful, failed, timed-out or setup-interrupted
claimed launch cannot retry the registration into a new directory. Invalid
preflight inputs must fail before acquiring the claim or starting a child.

These controls are correctness tests. Their whole-process test duration is
outside every IC online interval and is not a measured solver cost. No rho
comparison, calibrated timing, rate uncertainty or speedup is produced here.
Missing scientific costs remain unknown; there is no promotion gate to pass.

## Execution contract

New `sat_runtime_execution_v3.register` records include a versioned
`execution_policy` declaring exclusive creation of `execution-claim.json` in
the registration directory. The launcher validates the complete invocation,
sources, assets, interpreter, watchdog and separate output tree first. It then
atomically creates and flushes that claim **before** output setup or launch.
The claim is copied to the retained output and its canonical digest is checked
against the process and both source gates. A setup failure leaves a consumed
claim even if it cannot produce a complete process receipt.

The claim's original absolute directories are execution metadata, not candidate
identity. Archive copies may be audited at another location; retained claim
bytes must still match the process and gates. New launch code rejects legacy
registrations that lack this policy; their unmodified frozen code and existing
audit-only archive contracts remain historical evidence. Do not add a claim to
an old registration, delete a claim, or copy a consumed registration to retry.
This guard arbitrates cooperating local launchers per retained directory. It
does not attest against a malicious editor or registration copy; the external
registration/exposure ledger remains necessary.

## Audit command and boundaries

Freeze `prepared_runtime_transport_v1.py` with the complete runtime surface
before a future prepared control. After its one permitted execution, use a new
audit output directory outside the execution tree:

```sh
python3.12 research/ic_candidate_tournament_20260915/prepared_runtime_transport_v1.py audit \
  --execution /absolute/path/to/retained-execution \
  --expected-execution-sha256 EXTERNALLY_RECORDED_64_HEX_DIGEST \
  --out /absolute/path/to/new-independent-audit
```

The public CLI first checks the external seal, retained source gates and local
interpreter. It launches only the helper from `execution/extracted` under
`-I -S -B`, with a 180-second audit cap. That helper dispatches only the
registered prepared SAT v1 or F5 v2 mathematical auditor, checks its identities
and development claim boundary, and retains full import gates and admission
digests. Auditors inspect retained native receipts and mathematics; they do not
rerun a native solver. The accepted F5 v2 audit CLI remains available unchanged.
An archive lacking this helper cannot retroactively use it as a preexecution
auditor. A rejected audit is retained independently of the original execution
and never authorizes another solver attempt.

`PASS_FROZEN_PREPARED_TRANSPORT` means the independent transport passed.
Inspect `admission.json` for the family's complete/incomplete/failure status.
An incomplete or failed native control stays incomplete or failed even when
its transport passes. Every result keeps `promotion_eligible: false` and
`online_speedup: null`; disclosed controls remain ineligible for a headline
fresh single-target result.

## Remaining research gates

After this implementation is accepted, reconcile upstream executions and
freeze separate new one-shot disclosed SAT/F5 development controls with their
full sources, certificates, settings and external seals before invocation.
Retain complete, incomplete and failure outcomes. Then review the complete
historical/current exposure census, strong matched one-target references,
source-bound target-host builds, observer/calibration controls, resource limits
and arm order before any fresh paired protocol or target sampling. Neither
implementation acceptance nor synthetic transport tests closes those gates.
