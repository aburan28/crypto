# F5 v1: terminal native input-schema failure

The one permitted invocation of the [accepted development protocol](PROTOCOL.md)
is closed. The accepted controller was merged in
[PR #1056](https://github.com/aburan28/crypto/pull/1056) at
`f9e7883f2c6ae85b7b300383381388f6b69ce676`. Its external invocation was committed
before execution at `a596007e8796c0f5bea3dee87bff5a7040daa927`, in
[PR #1058](https://github.com/aburan28/crypto/pull/1058).
[PREEXECUTION.json](PREEXECUTION.json), the registration seal and protocol remain
unchanged; their pending status describes the historical preexecution state.
[TERMINAL.json](TERMINAL.json) records the outcome.

The controller submitted `target_seeds: [null]`. The pinned Rust worker's `Job`
field is `Vec<u64>` with an empty-list default. For an unseeded supplied public
point, its input list must be empty; `execute_job` then constructs `[None]` in
the output fixture's provenance. These two representations were conflated by
the Python registrar. Native JSON deserialization returned exit code 2 with
`invalid type: null, expected u64 at line 1 column 364`. The controller returned
exit code 1. Neither watchdog expired. Parsing failed before `execute_job` and
before any IC query; this is an adapter failure, with no algorithmic yield or
performance result. Build-identity preflight and Python/native source gates
passed. The failed entrypoint is deliberately rejected by complete-IC admission.

| Candidate / workload / run | Actual base / folded columns (registered) | Mathematical result | Online ns | Rho / IC | S / floor |
| --- | --- | --- | --- | --- | --- |
| `IC1N17Ckb1fb62PDP3f5RCsampleLAgaussTDpdpISO0he2ed5ed7bb97` / `1d07932114f4` / R0 | 62 / 29 | Unknown; native input rejected | Unknown | Unknown | Unknown |

The 6,275,834 ns native process interval and 1,863,702,917 ns controller interval
are operational diagnostics. They are not PDP costs or single-target online
times. Rank, relation yield, verified targets, scientific phase costs and
speedup stay `null`. This control does not admit a complete F4/F5 family, qualify
a fresh target or promote a challenger. The exposed n17 public point is
`[52411,72106]`; no target scalar was submitted. All older registrations and the
three sealed confirmation sets remain closed. No retry, resume or budget
extension is permitted for this invocation.

The failure revealed a gap in preflight: the prior controls checked declared
versus observed method records but did not deserialize the newly constructed
stdin through the actual native `Job` interface. Before another measurement,
a new versioned adapter must correct the seed representation and pass a
source-backed native deserialization control. A separate future receipt-writer
defect was also found during replay: raw native resource diagnostics contain
JSON floats, while `write_immutable` permits no floats. The failed-run publisher
does not include those floats and passed. The successful-IC publisher must be
corrected and tested under the new version; this defect did not cause the native
parse failure. Neither fix authorizes redispatching this consumed registration.

The next milestone remains complete source-bound F4/F5 admission, followed by a
separate reviewed exposure/reference/calibration protocol and fresh paired
incumbent/rho comparison. This is an accounting and reproducibility correction,
not a measured speed improvement.

## Durable evidence and replay

[results-20260930/evidence.tar.gz](results-20260930/evidence.tar.gz) retains the
complete registration and execution, exact interpreter/source/native assets,
sealed stdin, all gates and raw streams, and the publisher source. Its SHA-256 is
`9d166ca332e22d90738e62fbdcefd67ef2daa8d16425f162ad8ae3d0ff28962f`;
size is 37,902,366 bytes. The external invocation SHA-256 is
`5969837ddc69d25a3e44954f203077c0e6d682bbce2d1ac8eb0f9c3a10bb826b`.
[receipt.json](results-20260930/receipt.json) inventories every member's bytes,
hash and mode. [FAILURE.json](results-20260930/FAILURE.json) retains the stable
operational assessment; [TRANSPORT-REPLAY.json](results-20260930/TRANSPORT-REPLAY.json)
records fresh extraction and replay. [NATIVE-PARSE-REPLAY.json](results-20260930/NATIVE-PARSE-REPLAY.json)
records the independent native invocation/source/parse audit. This last receipt
contains raw floating resource diagnostics and is not a canonical candidate or
measurement identity record.

From this repository, replay to a new empty destination:

```sh
python3.12 research/ic_candidate_tournament_20260915/publish_sat_runtime_failure_v3.py replay \
  --bundle research/ic_candidate_tournament_20260915/goal_20260924/f5-source-bound-runtime-v1/results-20260930 \
  --out /tmp/ic-f5-v1-failure-replay-NEW \
  --expected-execution-sha256 5969837ddc69d25a3e44954f203077c0e6d682bbce2d1ac8eb0f9c3a10bb826b
python3.12 -m unittest discover -s research/ic_candidate_tournament_20260915 \
  -p test_f5_runtime_v1_terminal_failure.py -v
```

Replay checks byte retention and failure/source/invocation controls. It starts
no native worker, generates no new queries and makes no mathematical admission
from the retained parse error. The controls reject a wrong external seal,
modified stdin/output and attempts to relabel this failure as complete IC.
