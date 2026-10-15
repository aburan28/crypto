# Authorized basis/setup study — continuation

User approved the four-part bounded study and said "Going." This is execution,
not a reminder or permission request. Original matched study and PR1828 are
completed historical evidence; never rerun or rewrite them.

## First and isolation

Use cloud_environment.environment_status, current attached environment, then
read this file, /workspace/repos/crypto-n83-basis/AGENTS.md and
research/n83_basis_setup_20261015/PROTOCOL.md. Worktree:
/workspace/repos/crypto-n83-basis, branch codex/n83-basis-setup-20261015.
Sparse checkout: off-disk tracked files are materialized, not missing. About
3.8GiB disk free before training. Original /workspace/repos/crypto is untouched.

If any global /tmp/crypto-bench.lock reservation or measured worker/controller
is active, do no heavy work, builds/tests/replay/rendering/another solve. Do
not interrupt, do not foreground-poll and do not send redundant status. Keep
the completion watch enabled. New training exit file appears only after the
native reserve wrapper ends. If workspace unavailable, report concrete access
gap; do not invent results or start duplicate campaigns elsewhere.

## Prospective experiment

Exact legacy a=1/tail149 n83 r8569786107849059, G and lambda in protocol. Five
candidates: C0 retained K600 wide, C1 identical K600 dense, C2 next600 public
orbits after skipping first600 (disjoint from K600), C3 prefix1271 dense,
C4 prefix1798 dense. The initial C2 first600 draft equalled K600 exactly; it
was corrected BEFORE timing. Both draft and corrected preparation artifacts
are retained. Native audit verified full points/labels and prefix/window
relationships. Dense is layout, not mathematical base. C2 construction charges
all1200 generated orbits then discards/reindexes; all new candidates reconstruct
logs. C0/C1 start materialized K600; original construction is unpriced/unknown.

Training: two public targets generated outside1911 signed orbits including
original target, old16+32+64 targets, K600 and all1798prefix orbits. Generator
independently regenerated1798prefix with half-trace and verified old artifacts.
Training rank seeds29/43, same seed within each block across candidates.
Actual cold singleton5IC+rho on each training point, alternating candidate
order; five rho A/A pairs distributed2beforetraining,1aftertraining,
1betweencold/cohorts,1aftercohorts. First stage18processes/18answers includes
first3A/A pairs. Native write-once training matrix is frozen before launch.
All55processes/118answers across final study; no downscale/retry partial dirs.

Selection: minimum arithmetic mean actual complete-process wall of two
complete fully audited training runs; lexicographic candidate tie. Every failed
arm/cost retained. No eligible arm -> explicit no selection, duplicated C0
control only. Freeze selection receipt+identities+Git commit BEFORE creating
8fresh holdouts (seed202610156701). Holdout C0/selected/rho triplets on all8
points:24actualcoldprocesses, matched rank seeds101/103/107/109/113/127/131/137.
Then3genuinelyreused8targetcohorttriplets:9processes, rank seeds149/151/157,
rotations0/3/5, reversecohort2, sixmethodorderpermutations for cold/cyclecohorts.
Even selected=C0 gets independent duplicated process. Last2A/A pairs add4
processes. Total32ICprocesses/74ICanswers;23rhoprocesses/44rhoanswers.
Each CPU3/oneOSthread/12GiB/2100s. Rho Fold/reference-u64/256lanes/DP6/escape-off/
capfactor2000 fixed; per-seed G jump setup charged; no cross-Q DP sharing.
No1685targetexperiment, extrapolated crossover deferred.

## Paths and exact native tools

Training:
- /workspace/n83_basis_setup/training-plan/matrix.json (native frozen matrix)
- /workspace/n83_basis_setup/training-campaign (write-once measured output)
- /workspace/n83_basis_setup/training.stdout.jsonl, training.stderr.log
- /workspace/n83_basis_setup/training.reservation.jsonl
- /workspace/n83_basis_setup/training.exit.txt (only after wrapper/reserve ends)
- /workspace/n83_basis_setup/launch_training.sh operational wrapper only
- fixtures: /workspace/n83_basis_setup/training-fixtures
- preparation: /workspace/n83_basis_setup/preparation
- frozen binaries: /workspace/n83_basis_setup/frozen_binaries

Before postprocessing verify no measured worker/reservation active and check
study/provenance/BINARY_SHA256SUMS and SOURCE_SHA256SUMS. Never overwrite frozen
workers, native_runner, fixtures, protocol, raw records or completed reports.
Frozen checker/generator/selector source under study/verification. Run all
commands with cwd=/workspace/repos/crypto-n83-basis (__FILE__ relative sources).
Use NEW audit directories/receipts, no silent retry/overwrite:

```
/workspace/n83_evaluation/frozen_binaries/isolated_bench busy -- /workspace/n83_basis_setup/frozen_binaries/audit_training /workspace/n83_basis_setup/frozen_binaries/verify_rank /workspace/n83_basis_setup/training-plan/matrix.json /workspace/n83_basis_setup/training-campaign /workspace/n83_basis_setup/training-audit
/workspace/n83_evaluation/frozen_binaries/isolated_bench busy -- /workspace/n83_basis_setup/frozen_binaries/select_training /workspace/n83_basis_setup/training-plan/matrix.json /workspace/n83_basis_setup/training-campaign /workspace/n83_basis_setup/training-audit /workspace/n83_basis_setup/training-fixtures/selection_freeze.json
```

The audit controller separately meters each independent C++ full-rank/log
checker with wait4 wall/CPU/RSS. Audit costs never mix into IC/rho process costs.
A failed measured arm is skipped for full-rank PASS but its partial trace and
failure remain; no full-rank claim for incomplete traces. Completed IC audits
must verify every signed-orbit label/representative, scalar anchor, four-point
witness/intermediate, modular row, rank gain, backsolve, every [log]G, canonical
answer and logs-based target reconstruction. Rho scalars/A/A work replayed by
selector. If a genuine checker problem occurs, preserve failed receipt/source,
write new version and regress it; never edit raw output or silently relax gates.

Individual verifier syntax:
verify_rank BASE TRACE ANSWERS EXPECTED_RANK_SEED NEW_RECEIPT

Training fixture generator original command is documented in preparation and
source. Holdout generator (ONLY after committed selection freeze):
```
/workspace/n83_basis_setup/frozen_binaries/generate_targets /workspace/repos/crypto/experiments/koblitz-single-target-n83-20261006 /workspace/repos/crypto/research/n83_evaluation_20261012/optimization/verification/fixtures /workspace/repos/crypto/research/n83_speed_20261012/verification/fixtures /workspace/repos/crypto/research/n83_matched_20261014/verification/fixtures /workspace/repos/crypto/research/n83_dense_root_20261013/targets.jsonl /workspace/repos/crypto/research/n83_dense_root_20261013/answers.jsonl /workspace/n83_basis_setup/preparation/C4.base.jsonl /workspace/n83_basis_setup/training-fixtures /workspace/n83_basis_setup/holdout-fixtures
```

## Work still required after training

Do not stop at a training status. Complete authorized study through all4goals.
Native holdout runner and final statistics/report implementation are still
pending. After training reservation ends, implement a SEPARATELY VERSIONED
holdout runner (do not edit measured training native_runner). It can reuse the
native monitor source and instrumented frozen IC/rho binaries; it must freeze
all remaining37processes/100answers including2A/A pairs before timing, native
matrix/input/binary identities, write-once directories, phase boundaries that
release reservations for audit. No Python experiment logic, fixture generation,
solver/replay, statistics, bootstrap or selection. Rust/C++ only; Python only
presents recorded values, renders/inspects PDF, packages evidence.

Primary actual cold mean IC/rho on8targetblocks;20k paired percentile bootstrap
seed202610157801. Also selected/C0. Secondary batch8 has3independent seed/prep
clusters,20k bootstrap seed202610157802. Aligned online clocks include identical
generic public-Q validation, solve, generic replay, answer serialization/write/
flush. Timing receipt emission excluded both; input loading excluded online
but INCLUDED actual complete-process clocks. Legacy target_ms is diagnostic.
Report all individual targets/repeats and conventional midpoint median,p90,p95,
max,CPU/RSS,failures,A/A spread, relation attempts/probes/useful relations/rank
curves, exclusive anchor/search/check/row/elimination/backsolve times, phase and
attempt endpoint VmRSS/HWM; wait4 whole-process peak authoritative. Trace footer
records prior formatting/write/flush overhead, not every JSON construction,
memory sample or final footer; retain rank/process residuals rather than claim
full instrumentation calibration. Do not infer zero memory-limit hits.

All ratios IC/rho (<1favorsIC). Local cold win requires all declared holdout
answers/transcripts valid and primary95%upper<1; failures/incomplete arms
withhold unqualified winner. No S/floorcalib,L2,independent-host/global-optimum/
asymptotic/deployed-security claims. New selection support size isn't an
exponent fit. Preserve all failures/slow cases and failed preflight records.

Use report-evidence and research-visuals skills. Native analysis must precede
presentation. Produce stage/source-linked reports, graphsSVG/150dpiPNG/PDF,
inspected standalone addendum and updated comprehensive PDF preserving prior
87pages text+pixels72dpi and independent144dpi. Prior completed report:
/workspace/report_output/unified_curve_basis_report_matched_ic_rho.pdf,
SHAaf881641988a98bec3177e286694233f7b03731e82b6517a51ea78fc5e19793a.
Old70 and prior139launchrho study also unchanged. No timingmix/multipliedgains.
Put deliverables in /workspace/report_output and committed study/report. Copy
PDF/SVG/150dpiPNG to benchmarks/newdatedstudy, chronological index, affected
canonicaldiagnostic panels checked and reasons no ledger promotion. Hashmanifest
and evidenceZIP CRC checks; package existing values without Python recomputation.

Commit/push completed evidence, update new study PR, all current applicable CI
and review gates preserved. No reviewer without owner's name; never bypass
Cursor human review or merge unfinished study. PR1828 parent remains open and
unmerged at59e9c695b4df32bcdb0702b83cff86db1a761e4c; dependent PR base is
codex/n83-rho-speed-20261012 until parent merges, then appropriate retargeting.
Do not duplicate old PR follow-up. That separate task now has global lock guard.
PR follow-up50min/4h cadence, no foreground polling.

After inspected final PDFs/evidence/reviewable PR delivered, disable NEW study
completion watch and notify measured finding+PDF/evidence/PR+real reviewgate.
Never mark planned/running work complete. If study still running leavewatchon.

## Build/validation environment

Rust commands: PATH=/workspace/n83_evaluation/toolchain/rustup/toolchains/stable-x86_64-unknown-linux-gnu/bin:/usr/bin:/bin; RUSTUP_HOME=/workspace/n83_evaluation/toolchain/rustup; CARGO_HOME=/workspace/n83_evaluation/toolchain/cargo; CARGO_TARGET_DIR=/workspace/n83_evaluation/build.
Network-dependent unit tests require exec_command with_additional_permissions
network.enabled=true (local loopback TCP). No tests/guards disabled.
Newworkers19+4pass; runner6pass; independent checker1positive14negative;
selector synthetic tests exclude fast censored arm, select complete C1,
reject existing receipt and missing audit. Full library rerun4089pass/0fail/117ignored; see checks.
Initial missing sparse files, sandbox-loopback failures and JSONmacro recursion
limit failures are preserved. Core library unchanged; no measured timing yet at
initial handoff creation. No subagents unless user/applicable instructions ask.
