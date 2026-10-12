# Decision: K16 for the next isolated n37 online gate; K8 for cold setup

The [preregistered protocol](PROTOCOL.md) asked whether the preceding
K8 cold whole-solve instruction lead survives when only the previously
unseen public target is charged. It does not. On 16 new, same-point public
workloads, **K8/K16 target-only Callgrind Ir is 2.8860**, with the frozen
20,000-resample target-block 95% interval **[1.7302, 4.6015]**. Its lower
endpoint exceeds the preregistered 1.10 threshold, and the duplicate K16
target-only A/A maximum deviation is exactly zero. The frozen decision is
`prioritize_k16_for_isolated_n37_online_wall_gate`. Both bases remain live
at n41 and n53 until those field sizes are measured.

| Arm on the same 16 public points | Exact method or IC1 candidate | Actual usable base points | Folded columns | Verified profiles | Mean pre-online Ir | Mean target-only Ir | Mean post-online Ir | Mean complete-solve Ir |
|:--|:--|--:|--:|--:|--:|--:|--:|--:|
| Strong signed-Frobenius rho | `ECM1h1d9961ee4601` | — | — | 16/16 | 1,838,265 | 7,164,034 | 24,275 | 9,026,574 |
| IC K8 | `IC1N37Ckb0fb592PDP3mitmfrobeniuscountedRCsampleLAgaussTDpdpISO0hd7ebd078ae0c` | 592 | 8 | 16/16 | 31,013,817 | 3,171,035 | 1,774,255 | 35,959,107 |
| IC K16 | `IC1N37Ckb0fb1184PDP3mitmfrobeniuscountedRCsampleLAgaussTDpdpISO0hbbfdf029e5a2` | 1,184 | 16 | 16/16 | 37,950,369 | 1,098,762 | 3,335,460 | 42,384,590 |
| Identical K16 control | same K16 candidate | 1,184 | 16 | 16/16 | 37,950,332 | 1,098,762 | 3,335,457 | 42,384,550 |

Values are rounded from [DECISION.json](DECISION.json); exact per-point
integers and marker splits are in [CENSUS.json](CENSUS.json). The reported
unit is **Valgrind Callgrind's simulated user-space instructions**, not
retired instructions, group additions, CPU time, or wall time. In target-only
Ir, K8/rho is 0.44263 [0.25105, 0.71138] and K16/rho is **0.15337
[0.09234, 0.24000]**. These are same-point instruction diagnostics; neither
is an admitted IC/rho online speedup. The primary `rho_online_ms /
IC_online_ms` remains **unknown** because all native Mac runs earned L0,
below the required L2 CPU-isolation level.

The complete-solve instruction comparison reverses the base choice:
K8/K16 is **0.84840 [0.81856, 0.88374]**. K8/rho is 3.98369
[3.40432, 4.69251], and K16/rho is 4.69553 [4.02864, 5.53298]. The
K16 complete-solve A/A maximum relative deviation is 0.00282%, below the
frozen 2% limit. The mean pre-online interval is 31.01 million Ir for K8
versus 37.95 million for K16; the larger base also raises post-online
reporting work. Thus target-independent preparation and reporting outweigh
K16's lower target work in this implementation's complete solve. The
previous panel's K8 cold instruction choice is confirmed on new points,
while this panel resolves the previously open target-only instruction
choice in favor of K16.

The [sealed native session](sessions/mac_arm64_l0_01) recovered and checked
384/384 logarithms, including 320 measured executions. The [local full
audit](AUDIT.json) and independent Ubuntu x86-64
[receipt](independent_validation/RECEIPT.json) reproduced all 320 measured
records exactly on different environment classes. The 32 exact
[IC1 one-target claim records](candidate_claims) pair each first-round
K8/rho and K16/rho scalar on the same public point, with every `vs_rho`
correctness check passing. Their wall verdicts remain **descriptive only**
because both arms earned L0. The CI profiler recovered the archived scalar
in all 64 first-round children, and the independent
[Rust analyzer](../../examples/ecbench_n37_online_ir_analyze.rs) re-hashed
every raw file, re-parsed every numbered Callgrind part, checked the
point/method/seed/scalar against the sealed native session, verified the
other-host receipt, and reproduced the decision byte for byte.

The native stage records explain the target result without substituting a
stage prediction for it. Across 80 measured executions per arm, K8 made
6.8625 target PDP attempts on average and K16 made 1.8125; both stopped at
their first verified decomposition or the frozen attempt cap. K8 built 592
usable points, eight folded columns, and 2,344 pair-table entries, reaching
full rank after 40 target-blind rank trials and eight useful rows. K16 built
1,184 usable points, 16 folded columns, and 9,412 pair-table entries,
reaching full rank after 25 trials and 18 verified rows. Their mean
**counted** target group-addition-equivalent charges were 4,214.29 and
1,395.24 respectively. The supplementary cold counted comparison is
K16/K8 = 1.0198 [0.9835, 1.0562] over 80 pairs; K8/rho is 4.7661
[4.3033, 5.2847]. This ledger omits field, hash, allocation, and modular
work, so the quotient of incomplete counted costs does not bound true
IC/rho cost. Full failure and phase ledgers, including base construction,
relation collection, matrix solve, target PDP, and scalar replay, remain in
the [384 execution records](sessions/mac_arm64_l0_01/records.jsonl), the
[table](TABLE.txt), and the [paired comparisons](sessions/mac_arm64_l0_01/comparisons).

The profiles used one release binary under Valgrind 3.22.0 on an Ubuntu
24.04 x86-64 hosted VM reporting AMD EPYC 9V74. The raw archive is
[RAW-CALLGRIND.zip](RAW-CALLGRIND.zip), 6,148,912 bytes, SHA-256
`4adb5b41eef61161c80ef30618c6dcb3abb2dedef702730b659c0daece337c60`.
The frozen first-round [jobs](JOBS.json), native records, compact census,
[provenance](PROVENANCE.txt), independent receipt, and decision have SHA-256
values `a2d63e9063c98955e7fdf49b137d1c9baafc9912b903f8aa79cd69a43035a32f`,
`60838cb661b6df9ba801cca089975ce0d6bca259e1bb18f148461340fbe0f169`,
`d2d8b05ff8c19931bf039deac72db801115b8383b86df4058068c88b4996462f`,
`63141dec02f6016bf80c57705ab8921a613bc5c6934476a29664fc95a593c032`,
`b75740f1a28dd4a5352c15e0e738a82bf9c8d3402e3f6ea797d129c424e05472`,
and `e711979f64d08ea3417844b980e54e0fe095b4802bb803cf4f726c555be0e36a`
respectively. The profiler archived its generated `Cargo.lock`, exact
binary hash, CPU-feature report, child inputs and outputs, Valgrind logs,
and a SHA-256 receipt for every raw part. Extract the ZIP and run
`ecbench_n37_online_ir_analyze ARTIFACT_DIR JOBS.json RECORDS.jsonl OUT.json`
to reproduce the committed decision.

This result changes the **next measurement priority**, not the ECC2K-130
feasibility verdict. It supplies no isolated n37 online wall speedup,
fully charged wall crossover, n41/n53 result, n83 confidence result, or
n131/ECC2K-130 transfer claim. The next primary gate is an auditable
host-isolated, same-point one-target online IC/rho wall comparison for K16,
with K8 retained as a control; only then should the observed instruction
split guide larger-field base selection.
