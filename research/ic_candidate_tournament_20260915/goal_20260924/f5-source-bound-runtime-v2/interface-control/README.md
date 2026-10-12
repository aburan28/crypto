# Native Job/Config interface control passed

The single parser-only control preregistered at commit
`8bf6181bd0f39ae1296378730267de8d771ec18a`, in
[PR #1060](https://github.com/aburan28/crypto/pull/1060), completed and passed.
The corrected supplied-point job with empty input seeds and its omitted-seed
equivalent parsed through the actual private native `Job` and `Config` types.
Every native field/default matched the declared record. All six adversarial
inputs rejected with retained errors, including the historical `[null]` case.

This diagnostic never called `execute_job` or an IC query. It copied the exact
historical source and added only the frozen parser accessor and new example.
The admitted IC worker SHA-256 `94caf3d67e57dde09488763ec19e791aca54bbb35f1968c73e42d37a210b436a`
and historical worktree were unchanged. The diagnostic overlay has source
manifest SHA-256 `a71625d52addffe5cba74323f09ad17f1940b175f997397ee306c0b7485c9375`.
Build and native processes exited zero without timeout, with their frozen
1800/30-second bounds and owned process-group cleanup. This was a physical
macOS ARM64 correctness control; it establishes no other hardware execution,
natural yield, complete IC family, competitive online cost or speedup.

The Python runner presealed its direct analysis files, not a complete loaded
module execution manifest. The additional `generic_stages.py` DEFAULTS source
is explicitly retained as postexecution independent-auditor context from the
preregistered checkout. This does not upgrade that control to complete Python
execution attestation. The future IC controller separately requires its full
Python/package/interpreter/native pre/post gates. The native parser source,
complete Rust/dependency bytes, exact input and output are bound and replayed.

The v2 registrar requires byte-verified replay of this specific accepted proof
before creating a new measured registration. Its proof hash is in the canonical
candidate's flags. The prospective v2 IC run is still **not dispatched**;
its external invocation receipt must be committed against the accepted adapter
before the one permitted attempt. All v1 and historical registrations remain
closed. This parser control is also closed; replay starts no native binary.

The full [native/evidence.tar.gz](native/evidence.tar.gz) has SHA-256
`05dbb080423e219d7cc9d4fcfb4ae46645a9d09a423ef4f6a76c33b91bc84c58`
and 11,388,134 bytes. [native/receipt.json](native/receipt.json) inventories every
member. It retains source/dependency archives, original and appended example,
binary, compiler/interpreter/environment receipts, exact inputs, raw streams,
commands, results and independent auditor context.
[native/TRANSPORT-REPLAY.json](native/TRANSPORT-REPLAY.json) records the independent
byte/source/input/output replay, which passed without executing a solver.

From this repository:

```sh
python3.12 research/ic_candidate_tournament_20260915/replay_f5_job_schema_control.py \
  --bundle research/ic_candidate_tournament_20260915/goal_20260924/f5-source-bound-runtime-v2/interface-control/native \
  --expected-archive-sha256 05dbb080423e219d7cc9d4fcfb4ae46645a9d09a423ef4f6a76c33b91bc84c58
python3.12 -m unittest discover -s research/ic_candidate_tournament_20260915 \
  -p test_f5_job_schema_control.py -v
```

Controls reject altered sources, inputs, binary, output, independent DEFAULTS
context or external archive hash. Two valid parses and six rejected schemas
are an interface correctness result, with null DLP costs and no promotion.
