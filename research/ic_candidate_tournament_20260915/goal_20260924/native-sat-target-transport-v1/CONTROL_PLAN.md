# Disclosed CryptoMiniSat prestarted-stdin controls

This fixed the three source inputs for the
[prestarted-stdin feasibility probe](../../../../examples/ic_sat_stdin_probe.rs).
The [result](controls-v1/RESULT.md) retains all frozen runs, the failed
query-00 v1 and its separately labeled longer diagnostic. These are
controls from the closed native SAT registration, not fresh targets or new
natural-yield observations. Run only after both 512-query ordinary timing
panels release the shared busy lock. The accepted executable SHA-256 is
`6c509f09622f103d8a3ad90afc151e1c4031275c052d7f8af481465b4f27f2af`.
Use one million conflicts, one solver thread, one model, and a 60,000 ms
watchdog for both file and stdin modes. Preserve failures and original stdout,
stderr and source bytes. The probe source must compile and its hash must be
recorded before the first control. Neither a passing nor a failed control is
a complete IC result or a performance claim.

All paths below are relative to the repository root under
`research/ic_candidate_tournament_20260915/goal_20260924/native-sat-control-registration-v1/result-v1/data/execution/`.
For each `query-NN`, pass its `instance/instance.xor.cnf`, `instance/instance.anf`
and `instance/manifest.json` to the probe with the public point shown. Also
pass `native-sat-control-registration-v1/result-v1/data/preparation.json`;
its canonical hash is pinned in `PreparedState::load`. The probe checks the
manifest's target and all three export digests, independently checks the
point's curve/subgroup membership and enumerates the geometric three-sum
class before either solver invocation. Each SAT model must satisfy ANF, CNF
and a full-point three-sum lift.

| Query | Public point | Prior independently audited class | Manifest SHA-256 | ANF SHA-256 | CNF-XOR SHA-256 | Magma SHA-256 |
| --- | --- | --- | --- | --- | --- | --- |
| `query-00` | `[40991,73355]` | Exact geometric negative; source UNSAT | `d95a7101a3e69759cb1be39363fb305049159c23c5595e0ef82b6511447d7840` | `a7867583d28fe9593fd17dabdf6545d5bf79b0fda08888e5172f051a3d918eda` | `05b0913d27231ae82601b44f042765cf9ba2727aeb66543467e92127b3a2d92f` | `1e7609332d721678d763b1eb48658b505eae39506b0f6e4d0edff8261595eaab` |
| `query-01` | `[73003,104622]` | Exact geometric negative; source UNSAT | `58ffc396a551c424e19aa38af3addce803c6aaa7df3de899732d094fd663ca62` | `eab39a019e5298286e92169efcbe0eb381ce77af2d46749748d0ef202e1478ff` | `03a017d1763f5be128f8069f78158723a456b0803757eaa03ffbb647d9c6d960` | `c5f7baa8ef80be4d749fb3be5ac6c3d957978166d534826090455378ec6e5d7e` |
| `query-02` | `[59775,2910]` | Source-verified point witness; SAT model | `3ecd7c38a81d15b6849e538c89c1673f3e0a979d38e7365f2abb7cd2bfaa9e1e` | `1783698ff0fda61fc056b24b19d12e421f895c03ebbbc63fd399390b82def031` | `2077d013463e6a58270a38f29752ec0395099d1db0abfcfaadb7f6b8970d784c` | `958701bb9e154d2bea03b83e8243aa2b5dfc90278f421ddf3e21ac1d88159a0f` |

For each point, a transport pass requires the exact flushed stdin reader-ready
line before CNF delivery, unchanged accepted binary and copied source input, the
same recognized CMS status in file and stdin modes, and independent ANF/CNF
validation of any SAT model in both modes. A common native error or timeout
never passes. The three expected classes must match the audited controls;
budget-inconclusive and watchdog failure paths need separate disclosed tests
and retained receipts. A prestarted-stdin incompatibility remains an explicit
result and forces a new source-bound solver-server design before online SAT
admission.
The accepted CLI prints the line before stdin stream setup and parser construction.
Passing this transport control does not prove the stricter post-initialization
readiness needed to exclude all target-independent startup from the primary
online interval. That later handshake needs its own source-bound implementation
and disclosed parity tests.

### Explicit post-hoc transport diagnostic after the frozen controls

The frozen 60,000 ms control for `query-00` failed status parity: file mode
timed out after reaching about 375,000 conflicts, while stdin mode proved
UNSAT at 415,862 conflicts. The raw failed result remains in the separately
retained `ic-cms-stdin-control-00-v1` directory and must be published
alongside any follow-up. The other frozen controls passed: `query-01` was
UNSAT in both modes and `query-02` had independently verified source/model
and full-point witnesses in both modes. A 60-second wall watchdog on this
unisolated host cannot distinguish parser incompatibility from scheduling
noise when file mode stops just short of the observed conflict count.

Run one separately labeled `query-00-v2` diagnostic with **the same pinned
accepted executable, CNF/ANF/manifest/preparation bytes, one-million-conflict
budget, single thread and expected UNSAT class**, changing only the wall
watchdog to 120,000 ms, which the already committed probe permits. Use a new
create-only output directory and retain the v1 failure. A v2 parity pass
supports only the narrow claim that the accepted binary can parse and solve
this exact disclosed input from prestarted stdin. It does not retroactively
pass the frozen 60-second panel, show a speed benefit, establish a
post-initialization online boundary, or admit a complete SAT IC target.

Because the frozen accepted-stdin panel passed only two of three points, a
future marked-CMS three-way gate needs a new, explicit control protocol; the
longer post-hoc diagnostic cannot silently replace the failed 60-second row.
Under that new protocol, build a **new** instrumented binary from
the retained CMS source commit `7ae1b4a74259cdce223a584281fb8f090bbd3eed`
plus [the post-buffer marker patch](cms-stdin-ready-postbuffer.patch). The
earlier [one-line sketch](cms-stdin-ready.patch) is superseded because it
leaves `StreamBuffer` allocation after READY. Archive the complete
source, dependencies, build receipt and new binary hash before any target
registration. For each of the same three disclosed inputs, require three-way
agreement among accepted file mode, accepted stdin mode and marked stdin mode
on status and source-model validity. Check that the marked mode emits
`c PREPARED_STDIN_READY_v1` only after stdin stream and parser construction,
before any CNF byte is sent, with a retained PID/argv/stdout receipt. A different
valid SAT model is acceptable if it independently satisfies the same source
ANF and CNF and lifts to a valid three-point decomposition of the disclosed
query; byte-equality of solver stdout is not a mathematical requirement.
Preserve every failed build and control result. This parity gate does not
measure natural yield or solve a fresh public target.

The [prepared exporter mode](../../../../examples/koblitz_pdp_export.rs) has
its own disclosed three-point parity gate. The
[`ic_exporter_prestart_probe`](../../../../examples/ic_exporter_prestart_probe.rs)
executes the accepted exporter and a newly built prepared exporter on each
fixed public point, then checks all three fixed source-file SHA-256 values,
deterministic manifest parity, source-instance identity, readiness before
request delivery and process-group drain. It also sends malformed and
512-byte oversized one-request inputs to fresh prepared children, requiring
nonzero exits and no exported manifest; all raw stdout/stderr and child
receipts remain available. A final never-fed child must cancel with zero
stdin bytes and a confirmed process-group drain. Build it from a recorded new Rust
source/dependency snapshot; it is not the accepted ordinary-exporter binary.
For each row above, prestart a fresh process with static arguments
`17 6 standard 2026100303 1000000 NEW_OUTPUT 1 0 --export-only --prestart-stdin`
and **no** target arguments or target environment values. After the exact
`c EXPORTER_PREPARED_STDIN_READY_v1` marker, send one bounded JSON line with
the row's decimal `target_x`, `target_y` and
`blind_instance_id: native-control-NN`, then close stdin. Require exact ANF,
CNF-XOR and Magma SHA-256 values from the table, the same source-instance
identity and deterministic manifest fields as the accepted export, and an
independent full-point/group check of any SAT witness. Exclude only timing and
the new `transport` metadata from manifest parity; retain both. The exporter
must reject malformed and oversized requests. The target worker must prove by
its original send ledger that it never delivers request bytes before observing
the readiness marker; an OS pipe can otherwise buffer early bytes, which the
exporter cannot distinguish after it begins reading. These controls are
correctness/transport evidence, not natural-yield or online-speed measurements.
