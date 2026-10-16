# N83 native known-witness stage control

The single approved full-domain control reached SAT search and returned
`SAT_verified_witness_x`. The pinned macOS arm64 worker emitted all nine
ordered stage markers, reached `solver_entered`, returned a model after **787
conflicts**, and checked its decoded summand coordinates against the fixed
public witness. The independently replayed [launch audit](evidence-01/audit.json)
classifies the result as **`PASS_fixed_witness_encoding_only`**. This resolves
the earlier stage-entry ambiguity for this exact encoding and input.

The worker used the retained K=1,182 small signed-Frobenius base and K=2,048
residual pool in the five-summand ordered, all-affine known-answer control.
Before the fixed units it had 101,636 variables, 714,072 clauses, 4,980 XORs
and 465,738 domain clauses; after pinning eight x-coordinate blocks it had
714,736 clauses. The final [raw result](evidence-01/stdout.raw) reports
`model_verified=true`, `identity_mask=0` and the fixture and input digests.
The [raw markers](evidence-01/stderr.raw) show the ordered construction and
solver-entry sequence, followed by `solver_returned` and `model_checked`.

The [exit receipt](evidence-01/exit.json) records process wall time
**329.129273 seconds**, child user CPU **155.016457 seconds**, child system
CPU **4.258279 seconds**, a zero child exit, no watchdog timeout and verified
process-group cleanup. The launcher imposed an 840-second wall cap and a
5,000-conflict cap. The [authorization](authorization.json) conservatively
charges its full **900-second** allocation, bringing the local N83 charge from
22,200 to **23,100 seconds**. There was exactly one full-domain invocation
under that authorization.

The [prelaunch host receipt](host-prelaunch.txt) recorded load averages above
57 and 0% idle CPU. The [post-run receipt](host-postrun.txt) preserves the
later host state. Wall time is therefore a contended local observation, not
an isolated solver throughput estimate. This control used a planted public
known answer and did not sample natural relations, measure useful rank or
linear algebra, solve a fresh single target, or compare complete cold cost
with matched rho. The selected total-runtime factor base remains unset. The
worker checked its model internally, but the full SAT assignment was not
exported, so the independent replay verifies the raw trace and the fixture's
group/S3 algebra, not every assignment bit.

The launcher copied its authorization and wrote raw stdout, stderr, launch,
exit and audit receipts under `/private/tmp` to avoid the nearly full study
volume during execution. The complete 24-KiB directory was copied without
rewriting into [evidence-01](evidence-01/); replay of the copied files passed:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/known-witness-stage-launch-preflight-20261014/launch.py --check research/koblitz_n83_factor_base_sweep_20261008/verification/known-witness-stage-native-900-20261014/evidence-01
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/known-witness-first-pattern-replay-20261014/replay.py --check
python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008/verification/known-witness-stage-launch-preflight-20261014 -p 'test_*.py' -v
```

The first command recomputes input, source, binary, checker, raw-output and
receipt hashes, stage order, exit classification and process cleanup; it
returned nine markers and `PASS_fixed_witness_encoding_only`. The second
replayed the ordered public witness's group sum and S3 residuals. All four
focused launcher tests passed. All six copied files matched the temporary
source byte-for-byte; [SHA256SUMS](SHA256SUMS) seals those files plus the
authorization and two host receipts, and `shasum -a 256 -c SHA256SUMS` passed.
No new factor-base object was written to S3;
the 57 retained objects remain covered by the separate
[S3 stream integrity audit](../s3-stream-integrity-20261014/RESULT.md).
