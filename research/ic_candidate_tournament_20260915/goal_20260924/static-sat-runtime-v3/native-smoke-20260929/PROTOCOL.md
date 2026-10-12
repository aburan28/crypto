# Complete-source SAT native smoke

This is a new, bounded development control for the v3 controller. It never
redispatches a historical registration. It does not estimate natural yield,
qualify a complete solver, or enter a tournament ranking.

Hypothesis: the retained native inputs, isolated Python controllers and
inherited process-group meter suffice to construct and solve one source-bound
ordinary n17a1 PDP instance and independently account for its verdict. A
verified witness, a source UNSAT consistent with exact group enumeration,
and a retained budget-inconclusive or timeout verdict are admissible outcomes
of this input/source control. A wrong encoding, altered command, missing gate
or invalid witness fails the control; it cannot become a verified solve.

The exact curve is `EC1N17Ckb1hbbe2b5b6b1e6`. The standard dimension-six base
contains 63 geometric points, 62 distinct usable cofactor images and 29 folded
sign/Frobenius columns. `panel.json` freezes one ordinary query, CMS's one
million-conflict cap, 120-second solver watchdog, 60-second exporter watchdog,
one thread/model/random seed, and a 300-second whole-controller watchdog.
The trial-keyed query seed is selected here without a feasibility filter.
The supplied point `[52411,72106]` is disclosed development input, previously
exposed in the closed paired pilot. It is not a fresh qualification target;
no known scalar enters the invocation.

The reference is exhaustive three-sum existence on this exact geometric base,
computed only by the postexecution independent auditor. This is a correctness
reference, not Pollard rho and not a timing baseline. One ordinary query can
give rank at most one, so the anticipated terminal pipeline disposition is
`INCOMPLETE_RELATION_RANK`; there is no target-dependent online interval.
Preserve zero yield, censored outcomes, native errors and all partial output.
Do not retry a registered invocation, select another seed after seeing its
outcome, or extend this control's query cap.

The retained native bundle is `../native-inputs-macos-arm64/`, schema 3,
archive SHA-256
`a1fd5bd49c80076f3b64fd5cb51d891b278afde765b4d3e39853692ac318bd96`,
7,292,102 bytes. Its manifest and source/build receipts bind both native tools
and all retained source archives. This adapter validates the existing physical
macOS ARM64 build only. It does not admit these binaries for Linux execution.

Before dispatch, use `static_sat_registration_v3.py` to freeze the complete
current Python surface, interpreter/stdlib hashes, native bundle, canonical
method, candidate, workload, full invocation and watchdog. Save the generated
registration seal and externally record its execution hash. The final
candidate/run IDs are unresolved until that freeze; no result table may use
the design path as an `IC1` result.

```sh
PYTHONDONTWRITEBYTECODE=1 python3.12 research/ic_candidate_tournament_20260915/static_sat_registration_v3.py \
  --repository "$PWD" \
  --assets research/ic_candidate_tournament_20260915/goal_20260924/static-sat-runtime-v3/native-inputs-macos-arm64 \
  --panel research/ic_candidate_tournament_20260915/goal_20260924/static-sat-runtime-v3/native-smoke-20260929/panel.json \
  --out /absolute/new/registration

PYTHONDONTWRITEBYTECODE=1 python3.12 research/ic_candidate_tournament_20260915/static_sat_pipeline_v3.py \
  --registration /absolute/new/registration \
  --expected-execution-sha256 HASH_FROM_REGISTRATION_SEAL \
  --out /absolute/new/execution

PYTHONDONTWRITEBYTECODE=1 python3.12 research/ic_candidate_tournament_20260915/audit_static_sat_full_v3.py \
  --execution /absolute/new/execution \
  --expected-spec /absolute/new/registration/execution.json \
  --out /absolute/new/independent-audit.json
```

Native child CPU/RSS and whole-controller wall times are retained diagnostics.
They are not isolated/calibrated timings. No speedup, one-target online result,
or natural-yield inference is licensed. Transport the complete registration,
execution, stdout/stderr, each child's pre/post source gates, raw formula/model,
progress and independent audit into a durable hash-verified archive in the
follow-on result PR. A killed controller without its terminal source gate
remains incomplete even if a native verdict appeared before the watchdog.

Execution and independent replay are pending. The distinct follow-on full
development solve must preregister its own caps, source snapshot and seeds.
Fresh paired qualification remains pending on the complete F4/F5 family and a
reviewed exposure, reference, calibration and arm-order adapter. All three
sealed confirmation rounds remain closed.
