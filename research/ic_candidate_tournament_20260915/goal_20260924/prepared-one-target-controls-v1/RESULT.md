# Closed prepared one-target controls

The two preregistered native invocations are consumed and closed. Neither
recovered a scalar within its frozen eight-attempt cap. Preserve both failures;
do not retry, resume, extend, regenerate the registration or select a new seed
to replace either outcome. The complete family/comparison goal remains open.

[PR #1110](https://github.com/aburan28/crypto/pull/1110) retained the protocol,
full invocation seals and identities at commit
`c51a0ca12d3fbe7d22af2c0c9dfa4acbd198e3f7` before native execution. Its upstream
merge `b9367b65cc1febe1edd3ecf5bfc3fc78b29968e9` on 2026-10-01 at 03:33 UTC
accepted the preregistration, not these later results. Executed sources remain
the separately frozen accepted snapshot
`78d7d74af200275c516f77bc5cde86f6b5f62637`. This results follow-up changes no
executed Python or Rust solver, preparation certificate or raw report.

## Outcomes

Both arms received the already disclosed public point `[52411,72106]` on
`EC1N17Ckb1hbbe2b5b6b1e6`, without its scalar. They used their distinct certified
preparations of the same mathematical state: 63 geometric PDP points, 62 usable
subgroup points, 29 folded columns and rank 29. Target-independent preparation
replay precedes the target interval. Native execution was sequential on the
uncalibrated macOS ARM64 development host, under the repository heavy-work lock,
one native worker and explicit single-thread library settings. No memory cap,
matched rho, incumbent arm or fresh-target comparison was registered.

| Arm | Frozen candidate suffix / workload / run number | Retained target attempts | Exact geometric feasibility of those queries | Original frozen audit | Verified scalar / verified online time / speedup |
| --- | --- | --- | --- | --- | --- |
| F5 | `h6c8a5463ba52` / `6bb7d623c8ab` / `2026100101` | 8 native `proved_unsat` | 0 of 8; all negatives independently checked | Rejected: `missing query schema` | none / unknown / unknown |
| SAT | `h5f7f1076e104` / `b79558c422c7` / `2026100102` | 8 `CONFLICT_BUDGET_INCONCLUSIVE` | 1 of 8 feasible in the finite group | Pass for an **incomplete** source-bound control | none / unknown / unknown |

The full `(candidate_id, workload_id, run_id)` keys, terminal controller and
claim records, raw status mixes and unset scientific metrics are in
[outcome/result.json](outcome/result.json). Both controllers exited zero without
a watchdog timeout; this means orderly termination, not successful IC recovery.
Both one-use claims were acquired and retained. No native retry occurred.
F5's original frozen mathematical admission failed; SAT's original admission
explicitly records `scalar_verified: false` and `online_wall_ns: null`.

The raw attempted target intervals and phase ledgers are preserved in the
native reports and SAT summary. They end at the failed-attempt cap rather than
verified recovery. The F5 report also leaves relation-check/descent costs null.
Do not replace unknown verified cost with an attempted interval, treat missing
phases as zero, divide these uncalibrated clocks, or infer a winner. There is no
population-yield estimate or paired uncertainty estimate from this fixed
eight-query control. No new ordinary queries or relation-LA solves ran; earlier
natural collection evidence belongs to the two preparation certificates.

## Diagnosis and decision

F5's accepted prepared native report builder in
`examples/ic_tournament_worker.rs` emits no top-level `query_schema_version`.
The preexecution-frozen `generic_queries_exact_v1.verify_queries` requires that
header, so the original audit failed at `oracle.InvalidEvidence: missing query
schema`. Complete execution source gates passing does not repair that contract.
The previous Python tests used a mocked cold report, which already carried the
header; they did not establish the actual prepared native-to-auditor interface.

[DIAGNOSIS-PROTOCOL.md](DIAGNOSIS-PROTOCOL.md) and
[diagnose_controls.py](diagnose_controls.py) describe separate postexecution
analysis restricted by both external invocation seals to this disclosed n17
fixture. Complete geometric three-sum enumeration includes repeated indices,
killed torsion and identity pair sums. All eight F5 negatives are genuinely
absent from that domain. A labelled in-memory copy adding only
`query_schema_version: 1` replays the mathematical/query/negative-proof checks.
It changes no raw file and cannot replace the rejected original frozen audit
or make the incomplete invocation a complete solver. This establishes a
specific report-contract defect, with no false refutation observed on these
eight queries and no claim of general F5 soundness.

SAT's trial 1 query `[62577,27783]` has a group-readded three-sum witness with
geometric base indices `[29,51,2]`. All eight native SAT calls nevertheless
exhausted the registered conflict budget. This group witness does not establish
a Boolean assignment for the exported CNF; independently checking its source
constraints is the next diagnostic gate. The other seven group queries are
absent. Budget-inconclusive results made no UNSAT claims, so none is a false
refutation. Retain the feasible case for solver/encoding investigation rather
than replacing it with favourable random queries.

The first diagnosis script assumed `fixture['field']['n']` and failed with
`KeyError: 'field'` before any mathematical analysis or output creation. Its
[original error](diagnosis-initial-fixture-error.stderr) is retained. The corrected
script checks the independently constructed `Curve.n`; its digest is in
[outcome/diagnosis.json](outcome/diagnosis.json). This analysis-script correction
changed no native invocation, input, report or original audit.

The next implementation work must repair the native F5 header and add an actual
native-report-to-Python-auditor regression, then retain a newly source-bound
build under a new protocol. SAT needs a source-constraint check and bounded
solver/encoding diagnostics on the retained feasible query. No corrected build
or additional native execution is part of this evidence closeout. Both current
registrations remain closed even after a later fix. F5 is bounded Macaulay
elimination with an F5 row criterion; no general incremental Gröbner-basis
implementation or globally fastest IC is established.

## Durable evidence and replay

[publication/evidence.tar.gz](publication/evidence.tar.gz) retains all original
registration, claim, execution, native/Python assets, source gates, raw output,
query/export/CMS files and original audit files, plus the diagnosis. It has 703
regular-file members and 61,601,692 compressed bytes; SHA-256 is
`c94aba5c67afbe60b5109169d2d57c1c44ce89b218c63d7b023ed023bdb64b24`.
The [publication receipt](publication/receipt.json) gives every member's digest,
byte count and mode. The Git-tracked archive is the durable location; local
absolute paths in historical receipts remain provenance metadata.

[publication-source-v1.py](publication-source-v1.py) is the exact packer source
used to create that archive. The current [pack_evidence.py](pack_evidence.py)
additionally checks both external invocation seals before future packing and
rejects special permission bits on restoration. Packing/restoration executes
no solver and does not replace or mutate original evidence. Restore into a new
directory from this repository:

```sh
PYTHONDONTWRITEBYTECODE=1 python3.12 tools/isolated_bench.py busy -- \
  python3.12 research/ic_candidate_tournament_20260915/goal_20260924/prepared-one-target-controls-v1/pack_evidence.py restore \
  --bundle research/ic_candidate_tournament_20260915/goal_20260924/prepared-one-target-controls-v1/publication \
  --out /private/tmp/ic-prepared-controls-restored-new
```

For frozen audit replay on the original bound macOS interpreter, use
`prepared_runtime_transport_v1.py audit` with `<restored>/sat-execution` and
external seal `43539f7d440289dae1ad4867ba1ca951bd4664bf01147b9070eaa68f3bad07cb`,
or `<restored>/f5-execution` and
`8c5afe8d3f355010b4b293405f69e76a2362db25ae163490a1685955883b0737`.
Always use a new separate `--out`. F5 is expected to exit 1 with the original
missing-header rejection; that is the reproduced outcome. SAT should admit
only an incomplete control. These commands execute isolated Python auditors,
never a native worker. Do not invoke either `execute` registrar on restored
evidence. An archive copy is not an additional experiment authorization.

Actual restoration verified every member against the receipt. Local relocated
SAT admission and transport are byte-identical to the originals. F5 reproduces
the same status, source binding, exit code and missing-header failure; traceback
paths and thus stderr digests change on relocation. Both initial loaded-module
gates match their originals, and the independent diagnosis is byte-identical.
[outcome/restoration-replay.json](outcome/restoration-replay.json) and its compact
raw files retain those checks. This proves local relocation on the same bound
interpreter, not replay on Linux or another hardware build.

[outcome/new-exposures.json](outcome/new-exposures.json) retains the public target,
all 16 recorded query points and their sign/Frobenius closure: 578 distinct
points in total. It is a partial
current-control census, not the complete historical exclusion union. Merge it
with the preparation-only and all historical/current exposure records before
any new fresh-target protocol. Fresh references, calibrated observation, source-
bound comparison-host builds, resource/order freezes and a new paired protocol
remain outstanding. All three old confirmation sets stay closed.
